# Stage 5: optimal control of the dose schedule on the learned models.
# Same set-up as OptimalControlOne (02_.../src/optimal_control_one.jl) so results are comparable:
#   piecewise-constant daily dose over a 13-day horizon, mean dose fixed (default 1 uM) by the
#   same bisection projection, objective = total cells at the final time.
# Differences: the model is the recovered closed form (primary) or UDE ensemble (check); the
# objective is robust (ensemble mean + λ·sd); optimisation is gradient-based (ForwardDiff +
# L-BFGS, multistart) and cross-checked with BlackBoxOptim.

"""
Same projection as OptimalControlOne._project_mean: shift then clamp so mean(schedule) == target.
The shift is found by bisection on values; its derivative uses the implicit-function result
(∂shift/∂raw_j = -1/|free| for unclamped j) so ForwardDiff gradients are exact.
"""
function project_mean(raw, target, lo, hi)
    rv = ForwardDiff.value.(raw)
    left, right = lo - maximum(rv), hi - minimum(rv)
    for _ in 1:70
        shift = (left + right) / 2
        value = mean(clamp.(rv .+ shift, lo, hi))
        value < target ? (left = shift) : (right = shift)
    end
    s = (left + right) / 2
    free = findall(x -> lo < x + s < hi, rv)
    shift = isempty(free) ? s + zero(eltype(raw)) : s - sum(raw[free] .- rv[free]) / length(free)
    return clamp.(raw .+ shift, lo, hi)
end

daily_schedule(doses) = Schedule(collect(0.0:length(doses) - 1), collect(doses))

Base.@kwdef struct Scenario
    name::String
    S0::Float64
    R0::Float64
    low::Bool = false
    horizon::Float64 = 13.0
end

"Final (S, R) for a model under a daily schedule (no detection floor: true cell numbers)."
function final_state(m::UDEModel, sc::Scenario, doses)
    x = simulate_well(m.kind, m.θ, m.a, sc.low, Schedule(collect(0.0:length(doses) - 1), doses), 0.0,
                      sc.S0, sc.R0, zero(eltype(doses)), [sc.horizon])[1]
    return x[1], x[2]
end

"Trajectory on a fine grid for plotting."
function trajectory(m::UDEModel, sc::Scenario, doses; dt = 0.25)
    ts = collect(0.0:dt:sc.horizon)
    xs = simulate_well(m.kind, m.θ, m.a, sc.low, daily_schedule(doses), 0.0, sc.S0, sc.R0, 0.0, ts)
    return ts, [x[1] for x in xs], [x[2] for x in xs], [x[3] for x in xs]
end

"Robust objective: mean + λ·sd of log total over ensemble members."
function robust_objective(models, sc, doses; λ = 0.5)
    vals = [log(sum(final_state(m, sc, doses))) for m in models]
    μ = sum(vals) / length(vals)
    return length(vals) == 1 ? μ : μ + λ * sqrt(sum(abs2, vals .- μ) / (length(vals) - 1))
end

"""
Gradient-based optimisation over raw daily doses → projected schedule, multistart L-BFGS.
The projection is piecewise linear, so gradients exist almost everywhere.
"""
function optimise_gradient(models, sc; n = 13, target = 1.0, lo = 0.62, hi = 1.47, λ = 0.5, starts = 12, rng = StableRNG(5))
    f(raw) = robust_objective(models, sc, project_mean(raw, target, lo, hi); λ)
    g!(G, raw) = ForwardDiff.gradient!(G, f, raw)
    best = (val = Inf, sched = Float64[])
    inits = vcat([fill(target, n), collect(range(hi, lo; length = n)), collect(range(lo, hi; length = n))],
                 [lo .+ (hi - lo) .* rand(rng, n) for _ in 1:starts - 3])
    for x0 in inits
        res = Optim.optimize(f, g!, x0, LBFGS(), Optim.Options(iterations = 200))
        v = Optim.minimum(res)
        v < best.val && (best = (val = v, sched = project_mean(Optim.minimizer(res), target, lo, hi)))
    end
    return best
end

"BlackBoxOptim cross-check (same settings as OptimalControlOne)."
function optimise_bbo(models, sc; n = 13, target = 1.0, lo = 0.62, hi = 1.47, λ = 0.5, max_time = 20.0, seed = 2025)
    f(raw) = robust_objective(models, sc, project_mean(Float64.(raw), target, lo, hi); λ)
    Random.seed!(seed)
    res = BlackBoxOptim.bboptimize(f; SearchRange = (lo, hi), NumDimensions = n,
        Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = 120, MaxTime = max_time, TraceMode = :silent)
    raw = Vector{Float64}(BlackBoxOptim.best_candidate(res))
    return (val = BlackBoxOptim.best_fitness(res), sched = project_mean(raw, target, lo, hi))
end

"""
Extrapolation flags: doses outside the tested set's range, peak exposure above the training
maximum, and more on/off switches than the training schedules contain (≤ 1).
"""
function extrapolation_flags(doses; vmax_train, γ = 0.24)
    flags = String[]
    nz = filter(>(1e-6), doses)
    any(d -> d > 1.47 + 1e-6, doses) && push!(flags, "dose > 1.47 uM")
    any(d -> 1e-6 < d < 0.62 - 1e-6, doses) && push!(flags, "dose between 0 and 0.62 uM (untested)")
    s = daily_schedule(doses)
    vpeak = maximum(exposure(s, γ, t) for t in 0:0.25:length(doses))
    vpeak > vmax_train + 1e-3 && push!(flags, "peak exposure $(round(vpeak; digits = 2)) > training max $(round(vmax_train; digits = 2))")
    nsw = count(i -> abs(doses[i] - doses[i - 1]) > 0.2, 2:length(doses))
    nsw > 1 && push!(flags, "$(nsw) dose changes (training has ≤ 1)")
    return isempty(flags) ? "within training support" : join(flags, "; ")
end
