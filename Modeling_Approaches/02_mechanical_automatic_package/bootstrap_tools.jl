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

# circular moving-block resample of the residuals for days after day 0; with per-point scales w
# the standardised residuals e/w are resampled and rescaled by the scale of the point they land on
function bs_resample(rng, yhat, e, w=nothing)
    n = length(e)
    n <= 1 && return copy(yhat)
    z = w === nothing ? e : e ./ w
    pool = z[2:end]; m = length(pool)
    out = Float64[]
    while length(out) < m
        s = rand(rng, 1:m)
        for k in 0:BLOCK-1
            push!(out, pool[mod1(s + k, m)])
        end
    end
    w === nothing ? vcat(yhat[1] + e[1], max.(yhat[2:end] .+ out[1:m], 0.0)) :
                    vcat(yhat[1] + e[1], max.(yhat[2:end] .+ w[2:end] .* out[1:m], 0.0))
end

# Reproducibility: every optimiser call stops on its iteration cap; the wall-clock limit is only a
# safety net (without it the bounded BFGS can run for a very long time on a few candidates). Measured
# 2026-10-02 on 20 threads: no primary fit came near it, so the fits are deterministic. A call that
# reaches 95% of its limit may depend on the machine; it is recorded in LIMIT_HITS and run_zoo_fits.jl
# writes the list to outputs/zoo_fits/time_limit_hits.csv.
const LIMIT_HITS = NamedTuple[]
const LIMIT_LOCK = ReentrantLock()
function bs_fit(model, datasets, u0, p_start, bounds; maxtime=120.0, tag="")
    p0 = [clamp(p, lo + 1e-4*(hi - lo), hi - 1e-4*(hi - lo)) for (p, (lo, hi)) in zip(p_start, bounds)]
    best = (params=fill(NaN, length(p0)), bic=Inf, sse=Inf)
    for opt in (:bfgs, :nelder_mead)
        t0 = time()
        try
            r = run_joint_fit(model, gpe_datasets(datasets), u0, isfinite(best.bic) ? best.params : p0; bounds=bounds, optimizer=opt,
                              maxiters = opt == :bfgs ? 1500 : 2500, maxtime=maxtime, reltol=1e-7, abstol=1e-7)
            isfinite(r.bic) && r.bic < best.bic && (best = (params=r.params, bic=r.bic, sse=r.sse))
        catch
        end
        dt = time() - t0
        dt >= 0.95 * maxtime && lock(() -> push!(LIMIT_HITS, (fit=tag, optimizer=opt, seconds=dt, limit=maxtime)), LIMIT_LOCK)
    end
    best
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
    done = Threads.Atomic{Int}(0)
    Threads.@threads :dynamic for b in 1:B
        rng = Xoshiro(seed + b)
        ds = [hasproperty(d, :w) ? (x=d.x, y=bs_resample(rng, yh, e, d.w), state_index=d.state_index, w=d.w) :
                                   (x=d.x, y=bs_resample(rng, yh, e), state_index=d.state_index) for (d, yh, e) in zip(g.datasets, yhat, resid)]
        for (j, c) in enumerate(cands)
            if ks[j] == 0
                bic[b, j] = predict_candidate(c.name, c.label, c.model, ds, cu0(c, g)).bic
            else
                f = bs_fit(c.model, ds, cu0(c, g), c.params, c.bounds; tag="$(g.id) $(c.name) bootstrap $b")
                bic[b, j] = f.bic; pars[j][b, :] = f.params
            end
        end
        n = Threads.atomic_add!(done, 1) + 1
        n % 25 == 0 && (println("    $(g.id): $n / $B replicates"); flush(stdout))
    end
    winner = [argmin(bic[b, :]) for b in 1:B]
    pick = [choose(bic[b, :], ks) for b in 1:B]
    selected = argmax([count(==(j), pick) for j in 1:nc])
    (cands=cands, ks=ks, ref=ref, bic=bic, pars=pars, winner=winner, pick=pick, selected=selected,
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
