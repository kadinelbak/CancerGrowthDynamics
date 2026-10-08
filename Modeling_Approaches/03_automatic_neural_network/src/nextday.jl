# Stage 1: next-day prediction benchmark.
# Target: Δlog N from day t to t+1. A row needs measured counts on days t-1, t and t+1 (day 0 is a nominal
# seeding, not a measurement, so the first row is day 2). Co-culture partners are rolled out jointly.
#
# Features are autonomous (user rule: no clock time, the models drive optimal control). Removed on
# 2026-10-07: `day` (clock time), `since_switch` (time since the last dose change = clock time under a
# constant dose) and `cumdose` (∫u dt: under a constant dose it is dose × day, a treatment clock that never
# saturates). The dose history enters through `exposure`, the filtered dose dv/dt = γ(u − v), which is a
# state of the input that saturates and can be restarted from any state.

const FEATURES = [:logN, :dprev, :dose, :exposure,
                  :line_R, :low_resource, :s_frac, :seeding, :coculture, :log_well]
const GAMMA_FEATURE = 0.24

_obs_dict(tr) = Dict(zip(tr.days, log.(max.(tr.mean, 1.0))))

"Feature vector for trajectory `tr` at day t given (possibly predicted) log counts. `t` only selects the dose."
function feature_row(tr::Trajectory, t, logN, logNprev, log_well)
    s = tr.schedule
    return Float64[logN, logN - logNprev, dose_at(s, t), exposure(s, GAMMA_FEATURE, t),
                   tr.line == :R, tr.resource == :low, tr.s_frac, tr.seeding, tr.culture == :co, log_well]
end

"Group trajectories into wells (co-culture partners together)."
wells(trajs) = [filter(t -> t.group == g, trajs) for g in unique(t.group for t in trajs)]

function _well_log(members, obs, t)
    vals = [obs[i][t] for i in eachindex(members) if haskey(obs[i], t)]
    length(vals) == length(members) || return nothing
    return log(sum(exp.(vals)))
end

"One-step rows: (X features x rows, y, meta DataFrame)."
function nextday_rows(trajs)
    X = Vector{Vector{Float64}}(); y = Float64[]
    meta = DataFrame(traj_id = String[], group = String[], split = Symbol[], day = Float64[])
    for members in wells(trajs)
        obs = [_obs_dict(tr) for tr in members]
        for (i, tr) in enumerate(members), t in tr.days
            (t >= 1 && haskey(obs[i], t - 1) && haskey(obs[i], t + 1)) || continue
            lw = _well_log(members, obs, t)
            lw === nothing && continue
            push!(X, feature_row(tr, t, obs[i][t], obs[i][t - 1], lw))
            push!(y, obs[i][t + 1] - obs[i][t])
            push!(meta, (tr.id, tr.group, tr.split, t))
        end
    end
    return reduce(hcat, X), y, meta
end

# ---- models: all are functions X -> Δ ----
struct Ridge
    beta::Vector{Float64}; mu::Vector{Float64}; sd::Vector{Float64}
end
(r::Ridge)(X) = vec(r.beta' * vcat(ones(1, size(X, 2)), (X .- r.mu) ./ r.sd))
function fit_ridge(X, y; lambda = 1.0)
    mu = vec(mean(X; dims = 2)); sd = vec(std(X; dims = 2)); sd[sd .< 1e-8] .= 1.0
    Z = vcat(ones(1, size(X, 2)), (X .- mu) ./ sd)
    P = Diagonal(vcat(0.0, fill(lambda, size(X, 1))))
    return Ridge((Z * Z' + P) \ (Z * y), mu, sd)
end

struct TreeModel
    m::Any
end
_tree_table(X) = (; (f => X[i, :] for (i, f) in enumerate(FEATURES))...)
(t::TreeModel)(X) = EvoTrees.predict(t.m, _tree_table(X))
function fit_trees(X, y; seed = 11)
    cfg = EvoTreeRegressor(; nrounds = 250, eta = 0.04, max_depth = 3, lambda = 1.0, min_weight = 5.0,
                             rowsample = 0.8, colsample = 0.9, seed)
    data = merge(_tree_table(X), (; target = y))
    return TreeModel(EvoTrees.fit(cfg, data; target_name = :target, feature_names = FEATURES, verbosity = 0))
end

struct Fixed
    f::Function
end
(m::Fixed)(X) = m.f(X)

const DPREV = findfirst(==(:dprev), FEATURES)
# (name, fitter(X, y, seed), stochastic?) — stochastic models are refitted with NEXTDAY_SEEDS seeds and
# their scores reported as mean ± SD over seeds; the MLP ensemble averages 5 nets per seed.
const NEXTDAY_SEEDS = parse(Int, get(ENV, "NN_SEEDS", "5"))
const NEXTDAY_MODELS = [
    ("zero (no change)", (X, y, s) -> Fixed(X -> zeros(size(X, 2))), false),
    ("persistence", (X, y, s) -> Fixed(X -> X[DPREV, :]), false),
    ("ridge", (X, y, s) -> fit_ridge(X, y), false),
    ("boosted trees", (X, y, s) -> fit_trees(X, y; seed = 10 + s), true),
    ("MLP ensemble", (X, y, s) -> fit_mlp_ensemble(X, y; n = 5, hidden = 12, lambda = 3e-3, steps = 2500, seed0 = 100 * s), true),
]

"""
Joint recursive rollout of one well from start day t0 for up to `H` days.
Returns Dict((traj index, horizon) => (pred, obs)).
"""
function rollout_well(model, members, t0; H = 7)
    obs = [_obs_dict(tr) for tr in members]
    all(haskey(o, t0) && haskey(o, t0 - 1) for o in obs) || return Dict()
    cur = [o[t0] for o in obs]; prev = [o[t0 - 1] for o in obs]
    out = Dict{Tuple{Int, Int}, Tuple{Float64, Float64}}()
    for h in 1:H
        t = t0 + h - 1
        lw = log(sum(exp.(cur)))
        X = reduce(hcat, [feature_row(tr, t, cur[i], prev[i], lw) for (i, tr) in enumerate(members)])
        nxt = cur .+ model(X)
        for i in eachindex(members)
            haskey(obs[i], t + 1) && (out[(i, h)] = (nxt[i], obs[i][t + 1]))
        end
        prev, cur = cur, nxt
    end
    return out
end

function rollout_errors(model, members; H = 7)
    errs = [Float64[] for _ in 1:H]
    days = sort(unique(vcat((tr.days for tr in members)...)))
    for t0 in days
        t0 >= 1 || continue
        for ((i, h), (p, o)) in rollout_well(model, members, t0; H)
            push!(errs[h], p - o)
        end
    end
    return errs
end

"""
Leave-one-well-out CV on training wells + a final fit on all training rows evaluated on the held-out
validation wells. Stochastic models are refitted for every seed; one-step and rollout RMSEs are given per
seed, and `onestep`/`rollout` hold the mean and SD over seeds. Returns (onestep, rollout, oof, ...).
"""
function run_nextday(trajs; H = 7, seeds = NEXTDAY_SEEDS)
    train = training(trajs)
    X, y, meta = nextday_rows(trajs)
    is_train = meta.split .== :train
    train_wells = wells(train)
    per_seed = DataFrame(model = String[], seed = Int[], split = String[], rmse = Float64[], mae = Float64[], n = Int[])
    roll_seed = DataFrame(model = String[], seed = Int[], split = String[], horizon = Int[], rmse = Float64[], n = Int[])
    oof = copy(meta); oof.y = y
    for (name, fitter, stochastic) in NEXTDAY_MODELS, seed in (stochastic ? (1:seeds) : (1:1))
        pred = fill(NaN, length(y))
        rerr = [Float64[] for _ in 1:H]
        Threads.@threads for w in train_wells
            g = w[1].group
            tr_idx = is_train .& (meta.group .!= g)
            m = fitter(X[:, tr_idx], y[tr_idx], seed)
            te = findall(meta.group .== g)
            isempty(te) || (pred[te] = m(X[:, te]))
            e = rollout_errors(m, w; H)
            lock(_LOCK) do
                for h in 1:H
                    append!(rerr[h], e[h])
                end
            end
        end
        full = fitter(X[:, is_train], y[is_train], seed)
        va = findall(.!is_train)
        pred[va] = full(X[:, va])
        seed == 1 && (oof[!, Symbol(name)] = pred)
        for (split, idx) in [("train (LOWO CV)", findall(is_train)),
                             ("validation_treated_co", findall(meta.split .== :validation_treated_co))]
            r = pred[idx] .- y[idx]
            push!(per_seed, (name, seed, split, sqrt(mean(r .^ 2)), mean(abs.(r)), length(r)))
        end
        for h in 1:H
            push!(roll_seed, (name, seed, "train (LOWO CV)", h, sqrt(mean(rerr[h] .^ 2)), length(rerr[h])))
        end
        ve = [Float64[] for _ in 1:H]
        for w in wells(filter(t -> t.split == :validation_treated_co, trajs))
            e = rollout_errors(full, w; H)
            foreach(h -> append!(ve[h], e[h]), 1:H)
        end
        for h in 1:H
            push!(roll_seed, (name, seed, "validation_treated_co", h, sqrt(mean(ve[h] .^ 2)), length(ve[h])))
        end
        @info "next-day model done" name seed
    end
    sdz(x) = length(x) > 1 ? std(x) : 0.0
    onestep = combine(groupby(per_seed, [:model, :split]), :rmse => mean => :rmse, :rmse => sdz => :rmse_sd_over_seeds,
                      :mae => mean => :mae, :seed => length => :seeds, :n => first => :n)
    rollout = combine(groupby(roll_seed, [:model, :split, :horizon]), :rmse => mean => :rmse, :rmse => sdz => :rmse_sd_over_seeds,
                      :seed => length => :seeds, :n => first => :n)
    return (; onestep, rollout, per_seed, roll_seed, oof, X, y, meta)
end

const _LOCK = ReentrantLock()
