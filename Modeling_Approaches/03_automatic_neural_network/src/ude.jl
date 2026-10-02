# Stage 2: universal differential equation (UDE).
#
#   dS/dt = r_S S (1 - (S + a_SR R)/K_S) - δ_S(t) S - kill_S(v, ρ) S
#   dR/dt = r_R R (1 - (R + a_RS S)/K_R) - δ_R(t) R - kill_R(v, ρ) R
#   dv/dt = γ (u(t) - v)                   (filtered exposure, as in optimal_control_one.jl)
#   δ_L(t) = δ_L,low (1 - exp(-γ_d t))     (delayed death, low-resource wells only)
#   ρ = (S + R) / K_line                   (relative density)
#
# kill term variants
#   :nn    kill = v * softplus(MLP([v, ρ, line]))      (the UDE; no drug -> no kill)
#   :hill  kill = kmax_L * v^n / (EC50_L^n + v^n)       (mechanistic comparator)
#   :none  no drug effect
# Competition a_SR, a_RS is fit separately on untreated co-culture (Step B).
# Everything is fixed-step RK4 in log-space-safe form and ForwardDiff-differentiable.

softplus(x) = log1p(exp(-abs(x))) + max(x, zero(x))

# ---------------- generic tanh MLP on a flat parameter vector ----------------
struct MLPSpec
    sizes::Vector{Int}
end
nparams(s::MLPSpec) = sum(s.sizes[i + 1] * (s.sizes[i] + 1) for i in 1:length(s.sizes) - 1)
function mlp_eval(s::MLPSpec, θ, x)
    k = 0; h = x
    for i in 1:length(s.sizes) - 1
        nin, nout = s.sizes[i], s.sizes[i + 1]
        W = reshape(view(θ, k + 1:k + nin * nout), nout, nin); k += nin * nout
        b = view(θ, k + 1:k + nout); k += nout
        z = W * h .+ b
        h = i < length(s.sizes) - 1 ? tanh.(z) : z
    end
    return h
end
function mlp_init(s::MLPSpec, rng; scale = 1.0)
    θ = Float64[]
    for i in 1:length(s.sizes) - 1
        nin, nout = s.sizes[i], s.sizes[i + 1]
        append!(θ, randn(rng, nout * nin) .* (scale / sqrt(nin)))
        append!(θ, zeros(nout))
    end
    return θ
end

const KILL_NET = MLPSpec([3, 8, 8, 1])
const BLACKBOX_NET = MLPSpec([5, 8, 8, 1])   # pure neural ODE comparator: [ρ, v, line, low, t/14]

"Generic allocation-free 8-8-1 tanh net with NIN inputs."
function net8(θ, off, x::SVector{NIN}) where {NIN}
    W1 = SMatrix{8, NIN}(ntuple(i -> θ[off + i], Val(8 * NIN)))
    b1 = SVector{8}(ntuple(i -> θ[off + 8 * NIN + i], Val(8)))
    h1 = tanh.(W1 * x + b1)
    o2 = off + 8 * NIN + 8
    W2 = SMatrix{8, 8}(ntuple(i -> θ[o2 + i], Val(64)))
    b2 = SVector{8}(ntuple(i -> θ[o2 + 64 + i], Val(8)))
    h2 = tanh.(W2 * h1 + b2)
    o3 = o2 + 72
    w3 = SVector{8}(ntuple(i -> θ[o3 + i], Val(8)))
    return dot(w3, h2) + θ[o3 + 9]
end

"Allocation-free evaluation of the 3-8-8-1 kill network (same layout as mlp_eval)."
function kill_net(θ, off, x1, x2, x3)
    W1 = SMatrix{8, 3}(ntuple(i -> θ[off + i], Val(24)))
    b1 = SVector{8}(ntuple(i -> θ[off + 24 + i], Val(8)))
    h1 = tanh.(W1 * SVector(x1, x2, x3) + b1)
    o2 = off + 32
    W2 = SMatrix{8, 8}(ntuple(i -> θ[o2 + i], Val(64)))
    b2 = SVector{8}(ntuple(i -> θ[o2 + 64 + i], Val(8)))
    h2 = tanh.(W2 * h1 + b2)
    o3 = o2 + 72
    w3 = SVector{8}(ntuple(i -> θ[o3 + i], Val(8)))
    return dot(w3, h2) + θ[o3 + 9]
end

# ---------------- parameter layout ----------------
# backbone (per line L in S,R): log r_L, log K_L; low-resource: log r_L,low, log δ_L,low; shared log γ_d
# drug: log γ (exposure filter); then either NN weights or Hill (log kmax_L, log EC50_L, log n)
const BACKBONE_NAMES = [:logr_S, :logK_S, :logr_R, :logK_R, :logr_S_low, :logδ_S_low,
                        :logr_R_low, :logδ_R_low, :logγd, :logγ, :logfloor, :logfloor_low,
                        :log_r20_over_r30]      # growth-rate multiplier for 20k-seeded wells (density covariate)
const NB = length(BACKBONE_NAMES)
const HILL_NAMES = [:logkmax_S, :logEC50_S, :logkmax_R, :logEC50_R, :logn]

struct UDEModel
    kind::Symbol           # :nn, :hill, :none
    θ::Vector{Float64}     # backbone ++ kill parameters
    a::Vector{Float64}     # competition [a_SR, a_RS] (normal resource) ++ [a_SR, a_RS] (low resource)
end
nkill(kind) = kind == :nn ? nparams(KILL_NET) : kind == :hill ? length(HILL_NAMES) :
              kind == :blackbox ? nparams(BLACKBOX_NET) :
              kind == :symbolic ? length(SYMBOLIC_TERMS[].S) + length(SYMBOLIC_TERMS[].R) : 0

function init_theta(kind; rng = StableRNG(1))
    b = [log(0.9), log(3500.0), log(0.8), log(3500.0), log(0.9), log(0.3), log(0.8), log(0.5), log(0.3), log(0.24), log(5.0), log(10.0), 0.0]
    kind == :blackbox && return vcat(b, mlp_init(BLACKBOX_NET, rng; scale = 0.5))
    k = kind == :nn ? mlp_init(KILL_NET, rng; scale = 0.5) :
        kind == :hill ? [log(0.5), log(1.2), log(0.1), log(3.0), log(3.0)] : Float64[]
    kind == :nn && (k[end] = -2.0)   # start with a small kill
    return vcat(b, k)
end

const KILL_MASK = Ref((1.0, 1.0, 1.0))

"Kill rate (1/day) for line (±1), exposure v and relative density ρ."
function kill_rate(kind, θ, v, ρ, line_sign)
    kind == :none && return zero(v * θ[1])
    if kind == :nn
        m = KILL_MASK[]     # input ablation (Stage 4); (1, 1, 1) = full model
        return v * softplus(kill_net(θ, NB, m[1] * v, m[2] * ρ, m[3] * line_sign))
    end
    kind == :symbolic && return symbolic_kill(θ, v, ρ, line_sign)
    h = view(θ, NB + 1:length(θ))
    kmax, ec50 = line_sign < 0 ? (exp(h[1]), exp(h[2])) : (exp(h[3]), exp(h[4]))
    n = exp(h[5]); vv = max(v, zero(v)) + 1e-9
    return kmax * vv^n / (ec50^n + vv^n)
end

"Right-hand side for one well: x = [S, R, v]."
function well_rhs(kind, θ, a, low::Bool, sched::Schedule, t, x; d20::Bool = false)
    S, R, v = max(x[1], zero(x[1])), max(x[2], zero(x[2])), x[3]
    if kind == :blackbox
        ρ = (S + R) / 3500
        lw = low ? 1.0 : -1.0
        gS = net8(θ, NB, SVector(ρ, v, -1.0, lw, t / 14))
        gR = net8(θ, NB, SVector(ρ, v, 1.0, lw, t / 14))
        return (gS * S, gR * R, exp(θ[10]) * (dose_at(sched, t) - v))
    end
    rS, KS, rR, KR = exp(θ[1]), exp(θ[2]), exp(θ[3]), exp(θ[4])
    γ = exp(θ[10])
    δS = δR = zero(rS)
    if low
        rS, rR = exp(θ[5]), exp(θ[7])
        ramp = 1 - exp(-exp(θ[9]) * t)
        δS, δR = exp(θ[6]) * ramp, exp(θ[8]) * ramp
    end
    if d20                                  # density covariate: 20k wells grow at r × exp(θ[13])
        m20 = exp(θ[13]); rS *= m20; rR *= m20
    end
    kS = kill_rate(kind, θ, v, (S + R) / KS, -1.0)
    kR = kill_rate(kind, θ, v, (S + R) / KR, 1.0)
    aSR, aRS = low ? (a[3], a[4]) : (a[1], a[2])
    dS = rS * S * (1 - (S + aSR * R) / KS) - δS * S - kS * S
    dR = rR * R * (1 - (R + aRS * S) / KR) - δR * R - kR * R
    dv = γ * (dose_at(sched, t) - v)
    return (dS, dR, dv)
end

"""
Integrate a well from t0 with state (S0, R0, v0) and return states at `tobs` (sorted, ≥ t0).
Fixed-step RK4; steps are aligned so dose switches at integer days are hit exactly.
"""
function simulate_well(kind, θ, a, low, sched, t0, S0, R0, v0, tobs; dt = 0.1, d20::Bool = false)
    T = promote_type(eltype(θ), typeof(S0), typeof(R0), typeof(v0), eltype(a), eltype(sched.doses))
    x = (T(S0), T(R0), T(v0))
    t = Float64(t0)
    out = Vector{NTuple{3, T}}(undef, length(tobs))
    j = 1
    while j <= length(tobs) && tobs[j] <= t + 1e-9
        out[j] = x; j += 1
    end
    f(tt, xx) = well_rhs(kind, θ, a, low, sched, tt, xx; d20)
    while j <= length(tobs)
        h = min(dt, tobs[j] - t)
        k1 = f(t, x)
        k2 = f(t + h / 2, x .+ (h / 2) .* k1)
        k3 = f(t + h / 2, x .+ (h / 2) .* k2)
        k4 = f(t + h, x .+ h .* k3)
        x = x .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        x = (max(x[1], T(1e-6)), max(x[2], T(1e-6)), x[3])
        t += h
        while j <= length(tobs) && tobs[j] <= t + 1e-9
            out[j] = x; j += 1
        end
    end
    return out
end

# ---------------- data -> fitting units ----------------
const SIGMA_FLOOR = 0.15   # log-scale noise floor (plate-level artefacts, e.g. shared day-8 dips)

"A fitting unit is one well: its members' observations (log scale) from t0 = first day ≥ 1."
struct Unit
    members::Vector{Trajectory}
    t0::Float64
    low::Bool
    sched::Schedule
    area::Bool
    tobs::Vector{Float64}
    y::Vector{Vector{Float64}}      # per member, log obs on tobs (NaN = missing)
    w::Vector{Vector{Float64}}      # 1/σ²
    logN0::Vector{Float64}          # per member initial log value (count scale, observed)
end

function make_unit(members, apc)
    sched = members[1].schedule
    t0 = minimum(minimum(filter(>=(1.0), tr.days)) for tr in members)
    tobs = sort(unique(vcat((filter(>(t0), tr.days) for tr in members)...)))
    y = Vector{Float64}[]; w = Vector{Float64}[]; l0 = Float64[]
    for tr in members
        d = Dict(zip(tr.days, eachindex(tr.days)))
        yy = [haskey(d, t) ? log(max(tr.mean[d[t]], 1e-1)) : NaN for t in tobs]
        σ = [haskey(d, t) ? sqrt(SIGMA_FLOOR^2 + (tr.sem[d[t]] / max(tr.mean[d[t]], 1e-1))^2) : Inf for t in tobs]
        push!(y, yy); push!(w, 1 ./ σ .^ 2)
        i0 = get(d, t0, 0)
        push!(l0, i0 == 0 ? log(max(tr.mean[1], 1.0)) : log(max(count_scale(tr, apc)[i0], 1.0)))
    end
    return Unit(members, t0, members[1].resource == :low, sched, members[1].units == :area, tobs, y, w, l0)
end

make_units(trajs, apc) = [make_unit(w, apc) for w in wells(trajs)]
is20k(u) = u.members[1].seeding == 20

"Predicted log values per member on `u.tobs` from initial log values `l0` (count scale)."
function predict_unit(kind, θ, a, u::Unit, l0)
    S0 = zero(eltype(l0)); R0 = zero(eltype(l0))
    for (i, tr) in enumerate(u.members)
        tr.line == :S ? (S0 = exp(l0[i])) : (R0 = exp(l0[i]))
    end
    v0 = exposure(u.sched, exp(θ[10]), u.t0)
    xs = simulate_well(kind, θ, a, u.low, u.sched, u.t0, S0, R0, v0, u.tobs; d20 = is20k(u))
    floor = exp(u.low ? θ[12] : θ[11])   # detection floor (image-analysis background), count scale; separate for low resource
    return [[log((tr.line == :S ? x[1] : x[2]) + floor) for x in xs] for tr in u.members]
end

"Weighted SSE for a unit; area-unit wells get a profiled log offset (scale-free)."
function unit_sse(kind, θ, a, u::Unit, l0)
    pred = predict_unit(kind, θ, a, u, l0)
    s = zero(eltype(θ)) + zero(eltype(l0))
    for i in eachindex(u.members)
        r = [isnan(u.y[i][k]) ? zero(s) : u.y[i][k] - pred[i][k] for k in eachindex(u.tobs)]
        wk = [isnan(u.y[i][k]) ? 0.0 : u.w[i][k] for k in eachindex(u.tobs)]
        if u.area
            off = sum(wk .* r) / sum(wk)
            r = r .- off
        end
        s += sum(wk .* r .^ 2)
    end
    return s
end

# ---------------- Step A: fit backbone + kill on :drug units ----------------
const PRIOR_L0_SD = 0.3

function _pack(θ, l0s) ; vcat(θ, reduce(vcat, l0s)) ; end
function _unpack(p, nθ, units)
    θ = view(p, 1:nθ); k = nθ; l0s = []
    for u in units
        n = length(u.members); push!(l0s, view(p, k + 1:k + n)); k += n
    end
    return θ, l0s
end

"""
fit_drug_model(kind, units; ...) jointly fits backbone + kill + per-well initial conditions.
Adam warm-up then L-BFGS, gradients via ForwardDiff. Returns (θ, l0s, loss).
"""
function fit_drug_model(kind, units; θinit = init_theta(kind), l2 = 1e-3, adam_steps = 600, lbfgs_iters = 400,
                        lr = 0.02, verbose = false, weights = ones(length(units)))
    nθ = length(θinit)
    l0prior = [u.logN0 for u in units]
    p0 = _pack(θinit, l0prior)
    a0 = ones(4)
    function loss(p)
        θ, l0s = _unpack(p, nθ, units)
        s = zero(eltype(p))
        for (u, l0, lp, wu) in zip(units, l0s, l0prior, weights)
            s += wu * unit_sse(kind, θ, a0, u, l0)
            u.area || (s += sum(((l0 .- lp) ./ PRIOR_L0_SD) .^ 2))
        end
        kind in (:nn, :blackbox) && (s += l2 * sum(abs2, view(θ, NB + 1:nθ)) * length(units))
        # weak priors keep backbone rates in a biological range
        s += sum(abs2, (θ[[1, 3, 5, 7]] .- log(0.8)) ./ 1.5) + sum(abs2, (θ[[2, 4]] .- log(3500.0)) ./ 1.0)
        s += sum(abs2, (θ[10] - log(0.24)) / 1.0) + sum(abs2, (θ[9] - log(0.3)) / 1.5) + abs2(θ[11] - log(5.0)) + abs2(θ[12] - log(10.0)) + abs2(θ[13])
        return s
    end
    cfg = ForwardDiff.GradientConfig(loss, p0, ForwardDiff.Chunk{min(12, length(p0))}())
    g!(G, p) = ForwardDiff.gradient!(G, loss, p, cfg)
    p = copy(p0); m = zero(p); v = zero(p); G = zero(p)
    for k in 1:adam_steps
        g!(G, p)
        all(isfinite, G) || (G .= 0)
        @. m = 0.9m + 0.1G; @. v = 0.999v + 0.001G^2
        @. p -= lr * (m / (1 - 0.9^k)) / (sqrt(v / (1 - 0.999^k)) + 1e-8)
        verbose && k % 100 == 0 && @info "adam" k loss(p)
    end
    res = Optim.optimize(loss, g!, p, LBFGS(), Optim.Options(iterations = lbfgs_iters, g_tol = 1e-6))
    p = Optim.minimizer(res)
    θ, l0s = _unpack(p, nθ, units)
    return (θ = collect(θ), l0s = [collect(l) for l in l0s], loss = Optim.minimum(res))
end

# ---------------- Step B: competition on untreated co-culture ----------------
"""
Fit competition a (normal and low resource) and re-estimate the low-resource detection floor on
untreated co-culture; drug/backbone parameters stay fixed. Returns the calibrated UDEModel.
"""
function calibrate(kind, θ, units)
    l0prior = [u.logN0 for u in units]
    function loss(p)
        a = exp.(p[1:4])
        θc = vcat(θ[1:11], p[5], θ[13:end])
        s = zero(eltype(p)); k = 5
        for (u, lp) in zip(units, l0prior)
            n = length(u.members); l0 = p[k + 1:k + n]; k += n
            s += unit_sse(kind, θc, a, u, l0) + sum(((l0 .- lp) ./ PRIOR_L0_SD) .^ 2)
        end
        return s + sum(abs2, p[1:4] ./ 2) + abs2(p[5] - θ[12])   # weak priors around a = 1 and the drug-stage floor
    end
    p0 = vcat(zeros(4), θ[12], reduce(vcat, l0prior))
    g!(G, p) = ForwardDiff.gradient!(G, loss, p)
    res = Optim.optimize(loss, g!, p0, LBFGS(), Optim.Options(iterations = 300))
    p = Optim.minimizer(res)
    θc = copy(θ); θc[12] = p[5]
    return UDEModel(kind, θc, exp.(p[1:4]))
end

# ---------------- evaluation (forward prediction, observed initial condition) ----------------
"Per-point log residuals when predicting a unit forward from its observed day-t0 values (no refit)."
function forward_residuals(model::UDEModel, u::Unit)
    pred = predict_unit(model.kind, model.θ, model.a, u, u.logN0)
    rows = DataFrame(traj_id = String[], day = Float64[], obs = Float64[], pred = Float64[])
    for (i, tr) in enumerate(u.members)
        idx = findall(!isnan, u.y[i])
        r = u.y[i][idx] .- pred[i][idx]
        off = u.area ? mean(r) : 0.0     # area units: scale-free comparison
        for k in idx
            push!(rows, (tr.id, u.tobs[k], u.y[i][k], pred[i][k] + off))
        end
    end
    return rows
end

"Rolling one-step (next-day) prediction from each observed day (UDE analogue of Stage 1)."
function nextday_residuals(model::UDEModel, u::Unit)
    out = Float64[]
    γ = exp(model.θ[10])
    for (k, t) in enumerate(u.tobs)
        k == length(u.tobs) && break
        all(!isnan(u.y[i][k]) && !isnan(u.y[i][k + 1]) for i in eachindex(u.members)) || continue
        u.tobs[k + 1] - t == 1 || continue
        S0 = R0 = 0.0
        for (i, tr) in enumerate(u.members)
            val = exp(u.y[i][k])
            tr.line == :S ? (S0 = val) : (R0 = val)
        end
        u.area && (S0 /= _APC[]; R0 /= _APC[])
        x = simulate_well(model.kind, model.θ, model.a, u.low, u.sched, t, S0, R0, exposure(u.sched, γ, t), [t + 1]; d20 = is20k(u))[1]
        for (i, tr) in enumerate(u.members)
            p = log((tr.line == :S ? x[1] : x[2]) + exp(model.θ[u.low ? 12 : 11])) + (u.area ? log(_APC[]) : 0.0)
            push!(out, p - u.y[i][k + 1])
        end
    end
    return out
end
const _APC = Ref(1.0)

# ---------------- ensembles, evaluation helpers, persistence ----------------
"Bayesian-bootstrap ensemble: member 1 uses unit weights 1, others Dirichlet(1) weights × n."
function fit_ensemble(kind, units; n = 8, θwarm = nothing, kw...)
    members = Vector{Any}(undef, n)
    Threads.@threads for i in 1:n
        rng = StableRNG(1000 + i)
        w = i == 1 ? ones(length(units)) : (e = -log.(rand(rng, length(units))); e ./ sum(e) .* length(units))
        θ0 = init_theta(kind; rng)
        θwarm === nothing || (θ0[1:NB] .= θwarm[1:NB]; kind == :hill && (θ0 .= θwarm))
        members[i] = fit_drug_model(kind, units; θinit = θ0, weights = w, kw...)
    end
    return members
end

"Mean log prediction across an ensemble of UDEModels, plus 10/90% band."
function ensemble_forward(models::Vector{UDEModel}, u::Unit)
    rs = [forward_residuals(m, u) for m in models]
    base = rs[1][:, [:traj_id, :day, :obs]]
    P = reduce(hcat, [r.pred for r in rs])
    base.pred = vec(mean(P; dims = 2))
    base.lo = [quantile(P[i, :], 0.1) for i in 1:size(P, 1)]
    base.hi = [quantile(P[i, :], 0.9) for i in 1:size(P, 1)]
    return base
end

rmse(x) = sqrt(mean(abs2, x))

function save_models(path, dict)
    open(path, "w") do io
        JSON3.write(io, Dict(k => [Dict("kind" => String(m.kind), "theta" => m.θ, "a" => m.a) for m in v] for (k, v) in dict))
    end
end
function load_models(path)
    raw = JSON3.read(read(path, String))
    return Dict(String(k) => [UDEModel(Symbol(m.kind), Float64.(m.theta), Float64.(m.a)) for m in v] for (k, v) in raw)
end
