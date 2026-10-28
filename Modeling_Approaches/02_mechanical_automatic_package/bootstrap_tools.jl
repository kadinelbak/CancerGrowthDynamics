# Residual block-bootstrap model selection, shared by run_bootstrap_fits.jl and
# run_lr_mechanism_fits.jl. Include after run_new_data_staged_fits.jl (uses run_joint_fit,
# predict_candidate, dense, ODEProblem/solve/Tsit5).
#
# A group is (id, datasets, u0, cands, ...). A candidate may carry its own `u0` (models with
# hidden states such as resource or damaged cells); otherwise the group's u0 is used.

using Random: Xoshiro
using Statistics: median, quantile

const PER_PARAM_BIC = 3.0   # a simpler model wins if at most this much BIC worse per removed parameter
const BLOCK = 3             # residual block length in days

cu0(c, g) = hasproperty(c, :u0) ? c.u0 : g.u0

function bs_predict(model, u0, p, datasets)
    tmax = maximum(maximum(d.x) for d in datasets)
    sol = solve(ODEProblem(model, u0, (0.0, tmax), p), Tsit5(); reltol=1e-8, abstol=1e-8)
    [[sol(xi)[d.state_index] for xi in d.x] for d in datasets]
end

# circular moving-block resample of the residuals (every fitted point: day 0, the seeding, is not in the
# data); with per-point scales w the standardised residuals e/w are resampled and rescaled by the scale of
# the point they land on
function bs_resample(rng, yhat, e, w=nothing)
    m = length(e)
    m <= 1 && return copy(yhat)
    pool = w === nothing ? e : e ./ w
    out = Float64[]
    while length(out) < m
        s = rand(rng, 1:m)
        for k in 0:BLOCK-1
            push!(out, pool[mod1(s + k, m)])
        end
    end
    w === nothing ? max.(yhat .+ out[1:m], 0.0) : max.(yhat .+ w .* out[1:m], 0.0)
end

# Reproducibility: every optimiser call stops on its iteration cap (BFGS 1500, Nelder–Mead 2500); there is no
# wall-clock limit, so a fit does not depend on machine load, threads or sleep (review, 2026-10-07: the earlier
# 120 s / 300 s limits were hit 83 times). Bootstrap replicates also get a work budget counted in right-hand-side
# evaluations (2026-10-28: one Stage 4 replicate ran for 11+ h): a count, not seconds, so it is reproducible.
# An optimiser run over budget is discarded (its BIC stays Inf for that replicate) and listed in LIMIT_HITS.
const LIMIT_HITS = NamedTuple[]
const LIMIT_LOCK = ReentrantLock()
struct WorkBudgetHit <: Exception end
function counted(model, budget)
    n = Ref(0)
    f!(du, u, p, t) = ((n[] += 1) > budget > 0 && throw(WorkBudgetHit()); model(du, u, p, t))
    f!, n
end
# u0_builder (optional): the starting state as a function of the parameter vector (estimated starting values)
# budget: right-hand-side evaluations allowed per optimiser run (0 = none); evals: per-optimiser counts, filled in
function bs_fit(model, datasets, u0, p_start, bounds; tag="", u0_builder=nothing, budget=0, evals=nothing)
    p0 = [clamp(p, lo + 1e-4*(hi - lo), hi - 1e-4*(hi - lo)) for (p, (lo, hi)) in zip(p_start, bounds)]
    best = (params=fill(NaN, length(p0)), bic=Inf, sse=Inf)
    for opt in (:bfgs, :nelder_mead)
        f!, n = counted(model, budget); t0 = time()
        try
            r = run_joint_fit(f!, gpe_datasets(datasets), u0, isfinite(best.bic) ? best.params : p0; bounds=bounds, optimizer=opt,
                              maxiters = opt == :bfgs ? 1500 : 2500, reltol=1e-7, abstol=1e-7, u0_builder=u0_builder,
                              initial_time=0.0)   # u0 is the state at day 0 (the seeding), not at the first data point
            n[] <= budget || budget == 0 || throw(WorkBudgetHit())
            isfinite(r.bic) && r.bic < best.bic && (best = (params=r.params, bic=r.bic, sse=r.sse))
        catch
        end
        evals === nothing || push!(evals, n[])
        if budget > 0 && n[] > budget
            lock(LIMIT_LOCK) do
                push!(LIMIT_HITS, (fit=tag, optimizer=string(opt), seconds=time() - t0, limit=Float64(budget)))
            end
        end
    end
    best
end

# Starting-value fields of a candidate (zoo_data.jl x0_info): x0nom = the nominal seeding of every observed
# state (in sys_obs order), x0est = the u0 indices whose starting value is estimated, with bounds and starts
function x0_fields(sys, infos)
    obs = sys_obs(sys)
    length(infos) == length(obs) || error("one starting-value record per observed state expected")
    e = [i for i in eachindex(infos) if infos[i].est]
    (x0nom=[x.nom for x in infos], x0est=obs[e], x0bounds=[(infos[i].lo, infos[i].hi) for i in e], x0start=[infos[i].start for i in e])
end

# Fit one candidate from each of its starts and polish the best. Estimated starting values (x0est) are extra
# parameters after the model parameters, counted in the BIC; the fitted starting state is stored as u0 and
# `params` holds the model parameters only, so every later use (bootstrap, analyses, control) is unchanged.
function fit_cand(c, datasets)
    np = length(c.pnames); est = hasproperty(c, :x0est) ? c.x0est : Int[]; nx = length(est)
    builder = nx == 0 ? nothing : (p -> (u = Vector{eltype(p)}(c.u0); u[est] .= p[np+1:end]; u))
    if np + nx == 0
        f = predict_candidate(c.name, c.label, c.model, datasets, c.u0)
        return merge(c, (params=Float64[], bic=f.bic, sse=f.sse, se=Float64[], x0fit=Float64[]))
    end
    bounds = vcat(c.bounds, nx == 0 ? Tuple{Float64,Float64}[] : c.x0bounds)
    starts = np == 0 ? [Float64[]] : c.p0s
    best = (params=fill(NaN, np + nx), bic=Inf, sse=Inf)
    for p0 in starts
        f = bs_fit(c.model, datasets, c.u0, vcat(p0, nx == 0 ? Float64[] : c.x0start), bounds; tag=c.name, u0_builder=builder)
        f.bic < best.bic && (best = f)
    end
    if isfinite(best.bic)
        f = bs_fit(c.model, datasets, c.u0, best.params, bounds; tag=c.name, u0_builder=builder)
        f.bic < best.bic && (best = f)
    end
    isfinite(best.bic) || return merge(c, (params=fill(NaN, np), bic=Inf, sse=Inf, se=fill(NaN, np), x0fit=fill(NaN, nx)))
    u0 = builder === nothing ? c.u0 : Float64.(builder(best.params))
    pm = best.params[1:np]
    se = np == 0 ? Float64[] : (try param_se(c.model, datasets, u0, pm) catch; fill(NaN, np) end)
    merge(c, (params=pm, bic=best.bic, sse=best.sse, se=se, u0=u0, x0fit=best.params[np+1:end]))
end

function choose(bics, ks)
    best = argmin(bics)
    ok = [j for j in eachindex(bics) if isfinite(bics[j]) &&
          (j == best || (ks[j] < ks[best] && bics[j] - bics[best] <= PER_PARAM_BIC * (ks[best] - ks[j])))]
    sort(ok; by=j -> (ks[j], bics[j]))[1]
end

function bootstrap(g; B, seed=20260927)
    cands = g.cands; nc = length(cands)
    ks = [length(c.params) for c in cands]
    # replicates are generated from the most complex candidate (default) or from the carried model
    # (stages.toml [selection] bootstrap_generator = "selected"; candidate 1 is the carried one)
    gen = isdefined(Main, :selection_cfg) ? selection_cfg("bootstrap_generator", "most_complex") : "most_complex"
    gen in ("most_complex", "selected") || error("bootstrap_generator must be most_complex or selected")
    ref = gen == "selected" ? cands[1] : cands[sort(1:nc; by=j -> (-ks[j], cands[j].bic))[1]]
    yhat = bs_predict(ref.model, cu0(ref, g), ref.params, g.datasets)
    resid = [d.y .- yh for (d, yh) in zip(g.datasets, yhat)]
    bic = fill(Inf, B, nc)
    pars = [fill(NaN, B, ks[j]) for j in 1:nc]
    # work budget per candidate: factor × the evaluations of one fit to the original data from the estimates
    # (stages.toml [selection] bootstrap_work_factor; BOOT_WORK_FACTOR overrides; 0 = no budget)
    factor = parse(Float64, get(ENV, "BOOT_WORK_FACTOR", string(isdefined(Main, :selection_cfg) ? selection_cfg("bootstrap_work_factor", 50) : 50)))
    budget = zeros(Int, nc)
    if factor > 0
        Threads.@threads :dynamic for j in 1:nc
            ks[j] == 0 && continue
            ev = Int[]
            bs_fit(cands[j].model, g.datasets, cu0(cands[j], g), cands[j].params, cands[j].bounds; tag="$(g.id) $(cands[j].name) budget", evals=ev)
            budget[j] = max(round(Int, factor * maximum(ev)), 100_000_000)   # floor: a warm-start refit of a small model is cheap
        end
        println("    $(g.id): bootstrap work budget per optimiser run = $(factor)× a refit ($((any(ks .> 0) ? minimum(budget[ks .> 0]) : 0))–$(maximum(budget)) evaluations)"); flush(stdout)
    end
    done = Threads.Atomic{Int}(0)
    Threads.@threads :dynamic for b in 1:B
        rng = Xoshiro(seed + b)
        ds = [hasproperty(d, :w) ? (x=d.x, y=bs_resample(rng, yh, e, d.w), state_index=d.state_index, w=d.w) :
                                   (x=d.x, y=bs_resample(rng, yh, e), state_index=d.state_index) for (d, yh, e) in zip(g.datasets, yhat, resid)]
        for (j, c) in enumerate(cands)
            if ks[j] == 0
                bic[b, j] = predict_candidate(c.name, c.label, c.model, ds, cu0(c, g)).bic
            else
                f = bs_fit(c.model, ds, cu0(c, g), c.params, c.bounds; tag="$(g.id) $(c.name) bootstrap $b", budget=budget[j])
                bic[b, j] = f.bic; pars[j][b, :] = f.params
            end
        end
        n = Threads.atomic_add!(done, 1) + 1
        n % 25 == 0 && (println("    $(g.id): $n / $B replicates"); flush(stdout))
    end
    winner = [argmin(bic[b, :]) for b in 1:B]
    pick = [choose(bic[b, :], ks) for b in 1:B]
    selected = argmax([count(==(j), pick) for j in 1:nc])
    (cands=cands, ks=ks, ref=ref, bic=bic, pars=pars, winner=winner, pick=pick, selected=selected, budget=budget,
     orig_pick=choose([c.bic for c in cands], ks))
end

# 95% band of candidate j's bootstrap fits (default: the selected model), every state it has
function bands(g, res; tmax=XMAX, j=res.selected, dt=0.1)
    c = res.cands[j]; P = res.pars[j]; u0 = cu0(c, g)
    ts = collect(0.0:dt:tmax); nst = length(u0)
    curves = fill(NaN, size(P, 1), length(ts), nst)
    Threads.@threads :dynamic for b in 1:size(P, 1)
        p = size(P, 2) == 0 ? Float64[] : P[b, :]
        any(isnan, p) && continue
        try
            _, sol = dense(c.model, u0, p, tmax; dt=dt)
            length(sol.u) == length(ts) || continue
            for (i, u) in enumerate(sol.u), s in 1:nst
                curves[b, i, s] = u[s]
            end
        catch
        end
    end
    q(a, pr) = (v = filter(!isnan, a); isempty(v) ? NaN : quantile(v, pr))
    lo = [q(curves[:, i, s], 0.025) for i in eachindex(ts), s in 1:nst]
    hi = [q(curves[:, i, s], 0.975) for i in eachindex(ts), s in 1:nst]
    _, sol0 = dense(c.model, u0, c.params, tmax; dt=dt)
    mid = [sol0.u[i][s] for i in eachindex(ts), s in 1:nst]
    (ts=ts, lo=lo, hi=hi, mid=mid)
end

pct(v, j) = 100 * count(==(j), v) / length(v)
qs(v) = (w = filter(isfinite, v); isempty(w) ? (NaN, NaN, NaN) : (quantile(w, 0.025), median(w), quantile(w, 0.975)))
removable(lo, hi, bnd) = bnd[1] == 0.0 && isfinite(lo) && lo <= 0.02 * max(hi, 1e-12)
