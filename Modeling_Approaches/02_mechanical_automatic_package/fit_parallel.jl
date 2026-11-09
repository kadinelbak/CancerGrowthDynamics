# Parallel fitting of a list of candidates (user, 2026-11-09: spread the fitting over every core).
#
# fit_cands(cands, datasets) returns exactly what [fit_cand(c, datasets) for c in cands] returns
# (bootstrap_tools.jl), but every start of every candidate is its own task instead of the starts of one
# candidate running one after another on one thread. Before, a stage with 5–8 candidates kept only 5–8 of
# the cores busy, and the candidate with the slowest starts set the stage's time.
#
# Same numbers as fit_cand, by construction:
#   - the same starts (fixed ones, then the seeded lhs_starts), each fitted with the same bs_fit call and tag;
#   - the seeded starts' work budget needs the evaluations of every fixed start, so a candidate's seeded
#     starts begin once its fixed starts are done (as in the serial loop, where the fixed starts come first);
#   - the best start is chosen in start order with the same strict `<`, so ties go to the same start;
#   - the polish from the best start and the standard errors are the same calls.
# Each start is deterministic (no wall-clock limits, bootstrap_tools.jl), so the order the tasks run in does not
# matter. test/test_fit_parallel.jl compares the two on real stages.
#
# Nothing here changes bootstrap_tools.jl, so the saved fits (fit_cache/, keyed by its hash) stay valid.

function fit_cand_tasks(c, datasets; extra=EXTRA_STARTS)
    np = length(c.pnames); est = hasproperty(c, :x0est) ? c.x0est : Int[]; nx = length(est)
    np + nx == 0 && return fit_cand(c, datasets; extra=extra)          # no starts: a single prediction
    builder = nx == 0 ? nothing : (p -> (u = Vector{eltype(p)}(c.u0); u[est] .= p[np+1:end]; u))
    bounds = vcat(c.bounds, nx == 0 ? Tuple{Float64,Float64}[] : c.x0bounds)
    fixed = np == 0 ? [Float64[]] : c.p0s
    seeded = np == 0 ? Vector{Float64}[] : lhs_starts(c.name, c.bounds, extra, c.p0s)
    starts = vcat(fixed, seeded)
    nf = length(fixed)
    x0s = nx == 0 ? Float64[] : c.x0start
    fits = Vector{Any}(undef, length(starts))
    evs = [Int[] for _ in 1:nf]
    # fixed starts: no budget, their evaluation counts recorded
    @sync for i in 1:nf
        Threads.@spawn fits[i] = bs_fit(c.model, datasets, c.u0, vcat(starts[i], x0s), bounds; tag="$(c.name) start $i",
                                        u0_builder=builder, budget=0, evals=evs[i])
    end
    ev = reduce(vcat, evs; init=Int[])
    # seeded starts: EXTRA_WORK × the slowest fixed-start fit (the serial loop's rule)
    @sync for i in nf+1:length(starts)
        budget = EXTRA_WORK <= 0 || isempty(ev) ? 0 : max(round(Int, EXTRA_WORK * maximum(ev)), 10_000_000)
        Threads.@spawn fits[i] = bs_fit(c.model, datasets, c.u0, vcat(starts[i], x0s), bounds; tag="$(c.name) start $i",
                                        u0_builder=builder, budget=budget, evals=nothing)
    end
    best = (params=fill(NaN, np + nx), bic=Inf, sse=Inf)
    bics = Float64[]
    for f in fits
        push!(bics, f.bic)
        f.bic < best.bic && (best = f)
    end
    nbest = isfinite(best.bic) ? count(b -> b <= best.bic + START_BIC_TOL, bics) : 0
    if isfinite(best.bic)
        f = bs_fit(c.model, datasets, c.u0, best.params, bounds; tag=c.name, u0_builder=builder)
        f.bic < best.bic && (best = f)
    end
    isfinite(best.bic) || return merge(c, (params=fill(NaN, np), bic=Inf, sse=Inf, se=fill(NaN, np), x0fit=fill(NaN, nx),
                                           nstarts=length(starts), nbest=0))
    u0 = builder === nothing ? c.u0 : Float64.(builder(best.params))
    pm = best.params[1:np]
    se = np == 0 ? Float64[] : (try param_se(c.model, datasets, u0, pm) catch; fill(NaN, np) end)
    merge(c, (params=pm, bic=best.bic, sse=best.sse, se=se, u0=u0, x0fit=best.params[np+1:end], nstarts=length(starts), nbest=nbest))
end

# every candidate at once; each spawns its own start tasks, so the cores stay busy across candidates
function fit_cands(cands, datasets)
    tasks = [Threads.@spawn fit_cand_tasks(c, datasets) for c in cands]
    Any[fetch(t) for t in tasks]
end
