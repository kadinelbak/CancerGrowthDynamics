module OptimalControlOne

using ..IOUtils: TREATED_MONOCULTURE_IC_DOSE_MAP
using ..ModelPathTournament
using ..AdaptiveSimulationEngine
using ..A2780AdaptiveAdapter
using BlackBoxOptim
using CSV
using DataFrames
using Plots
using Random
using Statistics

export PaperParameters,
       A2780_DOSE_MAP,
       simulate_two_population,
       simulate_three_population,
       interpolate_parameters,
       run_optimal_control_one!

const PAPER_URL = "https://www.nature.com/articles/s41540-025-00511-3"
# These are legacy processed-folder names. The historical splitter exchanged the
# 0.67 and 1.47 uM source folders; IOUtils owns the canonical correction.
const A2780_DOSE_MAP = Dict(
    "IC25" => TREATED_MONOCULTURE_IC_DOSE_MAP[25],
    "IC50" => TREATED_MONOCULTURE_IC_DOSE_MAP[50],
    "IC75" => TREATED_MONOCULTURE_IC_DOSE_MAP[75],
)
const A2780_SOURCE_ORDER = ("IC75", "IC50", "IC25") # corrected 0.67, 1.0, 1.47 uM order
const DEFAULT_GAMMA = 0.24 # 0.01 / hour converted to 1 / day
const CANDIDATE_MODELS = (
    "dose-specific exponential",
    "shared-K logistic",
    "shared-growth delayed loss",
    "shared-growth delayed logistic",
    "paper two-state shared growth",
    "paper two-state dose-specific",
    "paper three-state shared growth",
)

struct PaperParameters
    rS::Float64
    dS::Float64
    alpha::Float64
    rR::Float64
    dR::Float64
end

function _rk4_step(f, t, x, h)
    k1 = f(t, x)
    k2 = f(t + h / 2, x .+ h .* k1 ./ 2)
    k3 = f(t + h / 2, x .+ h .* k2 ./ 2)
    k4 = f(t + h, x .+ h .* k3)
    return max.(x .+ h .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) ./ 6, 0.0)
end

function _integrate_at(times, x0, rhs; dt = 0.04)
    out = Matrix{Float64}(undef, length(times), length(x0))
    out[1, :] .= x0
    x = Float64.(x0)
    t = Float64(times[1])
    for i in 2:length(times)
        target = Float64(times[i])
        while t < target - 1e-12
            h = min(dt, target - t)
            x = _rk4_step(rhs, t, x, h)
            t += h
        end
        out[i, :] .= x
    end
    return out
end

function simulate_two_population(p::PaperParameters, times; gamma1 = DEFAULT_GAMMA, gamma2 = DEFAULT_GAMMA, dt = 0.04)
    rhs = function (t, x)
        S, R = x
        death_delay = 1 - exp(-gamma1 * t)
        induction_delay = 1 - exp(-gamma2 * t)
        return [
            p.rS * S - p.dS * death_delay * S - p.alpha * induction_delay * S,
            p.rR * R + p.alpha * induction_delay * S - p.dR * death_delay * R,
        ]
    end
    return _integrate_at(times, [1.0, 0.0], rhs; dt = dt)
end

function simulate_three_population(params, times; gamma1 = DEFAULT_GAMMA, gamma2 = DEFAULT_GAMMA, dt = 0.04)
    rS, dS, q, beta, rR, dR = params
    rhs = function (t, x)
        S, Q, R = x
        death_delay = 1 - exp(-gamma1 * t)
        induction_delay = 1 - exp(-gamma2 * t)
        return [
            rS * S - dS * death_delay * S - q * induction_delay * S,
            q * induction_delay * S - beta * Q,
            beta * Q + rR * R - dR * death_delay * R,
        ]
    end
    return _integrate_at(times, [1.0, 0.0, 0.0], rhs; dt = dt)
end

function _decode_two(z)
    rS, dS, alpha, resistant_growth_ratio, resistant_death_ratio = Float64.(z)
    return PaperParameters(rS, dS, alpha, resistant_growth_ratio * rS, resistant_death_ratio * dS)
end

function _decode_three(z)
    rS, dS, q, beta, resistant_growth_ratio, resistant_death_ratio = Float64.(z)
    return [rS, dS, q, beta, resistant_growth_ratio * rS, resistant_death_ratio * dS]
end

function _two_loss(z, times, observed)
    pred = simulate_two_population(_decode_two(z), times)
    total = pred[:, 1] .+ pred[:, 2]
    return all(isfinite, total) ? sum(abs.(total .- observed)) : 1e12
end

function _three_loss(z, times, observed)
    pred = simulate_three_population(_decode_three(z), times)
    total = vec(sum(pred; dims = 2))
    return all(isfinite, total) ? sum(abs.(total .- observed)) : 1e12
end

function _fit_model(loss, bounds; seed, max_time)
    Random.seed!(seed)
    result = bboptimize(
        loss;
        SearchRange = bounds,
        NumDimensions = length(bounds),
        Method = :adaptive_de_rand_1_bin_radiuslimited,
        PopulationSize = 80,
        MaxTime = max_time,
        TraceMode = :silent,
    )
    return Vector{Float64}(result.archive_output.best_candidate), Float64(result.archive_output.best_fitness)
end

function _candidate_definition(name, ndoses)
    if name == "dose-specific exponential"
        return (names = ["g[$i]" for i in 1:ndoses], bounds = fill((-1.5, 1.5), ndoses))
    elseif name == "shared-K logistic"
        return (names = vcat(["r[$i]" for i in 1:ndoses], ["K"]), bounds = vcat(fill((-1.5, 1.5), ndoses), [(1.01, 30.0)]))
    elseif name == "shared-growth delayed loss"
        return (names = vcat(["r"], ["d[$i]" for i in 1:ndoses]), bounds = vcat([(1e-4, 1.5)], fill((1e-5, 2.0), ndoses)))
    elseif name == "shared-growth delayed logistic"
        return (names = vcat(["r", "K"], ["d[$i]" for i in 1:ndoses]), bounds = vcat([(1e-4, 1.5), (1.01, 30.0)], fill((1e-5, 2.0), ndoses)))
    elseif name == "paper two-state shared growth"
        local_names = reduce(vcat, (["dS[$i]", "alpha[$i]", "dR_ratio[$i]"] for i in 1:ndoses))
        local_bounds = reduce(vcat, ([(1e-5, 2.0), (1e-5, 1.5), (1e-4, 1.0)] for _ in 1:ndoses))
        return (names = vcat(["rS", "rR_ratio"], local_names), bounds = vcat([(1e-4, 1.5), (1e-4, 1.0)], local_bounds))
    elseif name == "paper two-state dose-specific"
        names = reduce(vcat, (["rS[$i]", "dS[$i]", "alpha[$i]", "rR_ratio[$i]", "dR_ratio[$i]"] for i in 1:ndoses))
        bounds = reduce(vcat, ([(1e-4, 1.5), (1e-5, 2.0), (1e-5, 1.5), (1e-4, 1.0), (1e-4, 1.0)] for _ in 1:ndoses))
        return (names = names, bounds = bounds)
    elseif name == "paper three-state shared growth"
        local_names = reduce(vcat, (["dS[$i]", "q[$i]", "beta[$i]", "dR_ratio[$i]"] for i in 1:ndoses))
        local_bounds = reduce(vcat, ([(1e-5, 2.0), (1e-5, 1.5), (1e-4, 1.5), (1e-4, 1.0)] for _ in 1:ndoses))
        return (names = vcat(["rS", "rR_ratio"], local_names), bounds = vcat([(1e-4, 1.5), (1e-4, 1.0)], local_bounds))
    end
    error("Unknown candidate model: $name")
end

function _logdose_interpolate(dose, doses, values)
    dose <= first(doses) && return Float64(first(values))
    dose >= last(doses) && return Float64(last(values))
    hi = findfirst(d -> d >= dose, doses)
    lo = hi - 1
    w = (log(dose) - log(doses[lo])) / (log(doses[hi]) - log(doses[lo]))
    return (1 - w) * values[lo] + w * values[hi]
end

function _candidate_rates(name, z, doses, dose)
    n = length(doses)
    interp(values) = _logdose_interpolate(dose, doses, values)
    if name == "dose-specific exponential"
        return (kind = :exponential, g = interp(z[1:n]))
    elseif name == "shared-K logistic"
        return (kind = :logistic, r = interp(z[1:n]), K = z[n + 1])
    elseif name == "shared-growth delayed loss"
        return (kind = :delayed_loss, r = z[1], d = interp(z[2:n + 1]))
    elseif name == "shared-growth delayed logistic"
        return (kind = :delayed_logistic, r = z[1], K = z[2], d = interp(z[3:n + 2]))
    elseif name == "paper two-state shared growth"
        blocks = [z[3i:3i + 2] for i in 1:n]
        dS = interp(first.(blocks)); alpha = interp(getindex.(blocks, 2)); ratio = interp(last.(blocks))
        return (kind = :paper_two, rS = z[1], rR = z[2] * z[1], dS = dS, alpha = alpha, dR = ratio * dS)
    elseif name == "paper two-state dose-specific"
        blocks = [z[5i - 4:5i] for i in 1:n]
        rS = interp(first.(blocks)); dS = interp(getindex.(blocks, 2)); alpha = interp(getindex.(blocks, 3))
        return (kind = :paper_two, rS = rS, rR = interp(getindex.(blocks, 4)) * rS, dS = dS, alpha = alpha, dR = interp(last.(blocks)) * dS)
    elseif name == "paper three-state shared growth"
        blocks = [z[4i - 1:4i + 2] for i in 1:n]
        dS = interp(first.(blocks)); q = interp(getindex.(blocks, 2)); beta = interp(getindex.(blocks, 3)); ratio = interp(last.(blocks))
        return (kind = :paper_three, rS = z[1], rR = z[2] * z[1], dS = dS, q = q, beta = beta, dR = ratio * dS)
    end
    error("Unknown candidate model: $name")
end

function _simulate_candidate(name, z, doses, dose, times; dt = 0.04)
    p = _candidate_rates(name, z, doses, dose)
    if p.kind == :exponential
        total = exp.(p.g .* Float64.(times))
        return (total = total, sensitive = total, latent = zeros(length(total)))
    elseif p.kind == :logistic
        rhs = (t, x) -> [p.r * x[1] * (1 - x[1] / p.K)]
        states = _integrate_at(times, [1.0], rhs; dt = dt)
        return (total = states[:, 1], sensitive = states[:, 1], latent = zeros(length(times)))
    elseif p.kind == :delayed_loss || p.kind == :delayed_logistic
        rhs = function (t, x)
            growth = p.kind == :delayed_logistic ? p.r * x[1] * (1 - x[1] / p.K) : p.r * x[1]
            [growth - p.d * (1 - exp(-DEFAULT_GAMMA * t)) * x[1]]
        end
        states = _integrate_at(times, [1.0], rhs; dt = dt)
        return (total = states[:, 1], sensitive = states[:, 1], latent = zeros(length(times)))
    elseif p.kind == :paper_two
        states = simulate_two_population(PaperParameters(p.rS, p.dS, p.alpha, p.rR, p.dR), times; dt = dt)
        return (total = states[:, 1] .+ states[:, 2], sensitive = states[:, 1], latent = states[:, 2])
    else
        states = simulate_three_population([p.rS, p.dS, p.q, p.beta, p.rR, p.dR], times; dt = dt)
        return (total = vec(sum(states; dims = 2)), sensitive = states[:, 1], latent = states[:, 2] .+ states[:, 3])
    end
end

function _fit_candidate(name, doses, curves; seed, max_time)
    definition = _candidate_definition(name, length(doses))
    objective = function (z)
        sse = 0.0
        for (i, curve) in enumerate(curves)
            prediction = _simulate_candidate(name, z, doses, doses[i], curve.time).total
            all(isfinite, prediction) || return 1e12
            sse += sum((prediction .- curve.normalized) .^ 2)
        end
        return sse
    end
    z, sse = _fit_model(objective, definition.bounds; seed = seed, max_time = max_time)
    nobs = sum(length(curve.time) for curve in curves)
    bic = nobs * log(max(sse, 1e-12) / nobs) + length(z) * log(nobs)
    boundary = [definition.names[i] for i in eachindex(z) if min(z[i] - definition.bounds[i][1], definition.bounds[i][2] - z[i]) <= 0.005 * (definition.bounds[i][2] - definition.bounds[i][1])]
    eligible = isfinite(sse) && sse < 1e11 && all(isfinite, z)
    return (name = name, params = z, doses = Float64.(doses), names = definition.names, bounds = definition.bounds, sse = sse, bic = bic, boundary = boundary, eligible = eligible)
end

function _read_curve(repository_root, lineage, ic_label)
    path = joinpath(
        repository_root,
        "Processed_Datasets",
        "Treated MonoCulture",
        "30k",
        ic_label,
        "Averages",
        "$(lineage)_day_averages.csv",
    )
    df = CSV.read(path, DataFrame)
    means = Float64.(df[!, "Mean Cells"])
    sem = Float64.(df[!, "SEM Cells"])
    return (
        day = Float64.(df.Day),
        time = Float64.(df.Day .- first(df.Day)),
        mean = means,
        normalized = means ./ means[1],
        sem_normalized = sem ./ means[1],
        source = path,
    )
end

function interpolate_parameters(dose, fitted::Dict{Float64, PaperParameters})
    doses = sort(collect(keys(fitted)))
    dose_clamped = clamp(Float64(dose), first(doses), last(doses))
    dose_clamped <= first(doses) && return fitted[first(doses)]
    dose_clamped >= last(doses) && return fitted[last(doses)]
    hi = findfirst(d -> d >= dose_clamped, doses)
    lo = hi - 1
    w = (log10(dose_clamped) - log10(doses[lo])) / (log10(doses[hi]) - log10(doses[lo]))
    a, b = fitted[doses[lo]], fitted[doses[hi]]
    vals = (1 - w) .* [a.rS, a.dS, a.alpha, a.rR, a.dR] .+ w .* [b.rS, b.dS, b.alpha, b.rR, b.dR]
    return PaperParameters(vals...)
end

function _simulate_control(fitted, schedule; horizon = 13.0, dt = 0.025, gamma1 = DEFAULT_GAMMA, gamma2 = DEFAULT_GAMMA)
    n = length(schedule)
    times = collect(0.0:dt:horizon)
    last(times) < horizon && push!(times, horizon)
    rhs = function (t, x)
        index = clamp(floor(Int, t / horizon * n) + 1, 1, n)
        dose = schedule[index]
        p = interpolate_parameters(dose, fitted)
        S, R, vdS, vdR, valpha, z = x
        return [
            p.rS * S - vdS * S - valpha * S,
            p.rR * R + valpha * S - vdR * R,
            gamma1 * (p.dS - vdS),
            gamma1 * (p.dR - vdR),
            gamma2 * (p.alpha - valpha),
            dose,
        ]
    end
    states = _integrate_at(times, zeros(6), rhs; dt = dt)
    states[1, 1] = 1.0
    # Reintegrate after setting the initial sensitive population.
    states = _integrate_at(times, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0], rhs; dt = dt)
    return (time = times, states = states, total = states[:, 1] .+ states[:, 2])
end

function _project_mean(raw, target, lo, hi)
    left, right = lo - maximum(raw), hi - minimum(raw)
    for _ in 1:70
        shift = (left + right) / 2
        value = mean(clamp.(raw .+ shift, lo, hi))
        value < target ? (left = shift) : (right = shift)
    end
    return clamp.(raw .+ (left + right) / 2, lo, hi)
end

function _two_level_schedule(n, lo, hi, target; front = true)
    raw = front ? collect(range(hi, lo; length = n)) : collect(range(lo, hi; length = n))
    return _project_mean(raw, target, lo, hi)
end

function _optimize_control(fitted; horizon = 13.0, intervals = 13, target_dose = 1.0, max_time = 20.0, seed = 2025)
    lo, hi = extrema(collect(keys(fitted)))
    objective = function (raw)
        schedule = _project_mean(Float64.(raw), target_dose, lo, hi)
        sim = _simulate_control(fitted, schedule; horizon = horizon, dt = 0.05)
        return sim.total[end]
    end
    Random.seed!(seed)
    result = bboptimize(
        objective;
        SearchRange = (lo, hi),
        NumDimensions = intervals,
        Method = :adaptive_de_rand_1_bin_radiuslimited,
        PopulationSize = 120,
        MaxTime = max_time,
        TraceMode = :silent,
    )
    raw = Vector{Float64}(result.archive_output.best_candidate)
    return _project_mean(raw, target_dose, lo, hi)
end

function _simulate_selected_control(fit, schedule; horizon = 13.0, dt = 0.025)
    times = collect(0.0:dt:horizon)
    last(times) < horizon && push!(times, horizon)
    dose_at(t) = schedule[clamp(floor(Int, t / horizon * length(schedule)) + 1, 1, length(schedule))]
    initial, rhs = if startswith(fit.name, "paper two-state")
        ([1.0, 0.0, 0.0, 0.0, 0.0], function (t, x)
            p = _candidate_rates(fit.name, fit.params, fit.doses, dose_at(t))
            S, R, vdS, vdR, valpha = x
            [p.rS * S - vdS * S - valpha * S,
             p.rR * R + valpha * S - vdR * R,
             DEFAULT_GAMMA * (p.dS - vdS),
             DEFAULT_GAMMA * (p.dR - vdR),
             DEFAULT_GAMMA * (p.alpha - valpha)]
        end)
    elseif fit.name == "paper three-state shared growth"
        ([1.0, 0.0, 0.0, 0.0, 0.0, 0.0], function (t, x)
            p = _candidate_rates(fit.name, fit.params, fit.doses, dose_at(t))
            S, Q, R, vdS, vdR, vq = x
            [p.rS * S - vdS * S - vq * S,
             vq * S - p.beta * Q,
             p.beta * Q + p.rR * R - vdR * R,
             DEFAULT_GAMMA * (p.dS - vdS),
             DEFAULT_GAMMA * (p.dR - vdR),
             DEFAULT_GAMMA * (p.q - vq)]
        end)
    elseif fit.name == "shared-growth delayed loss" || fit.name == "shared-growth delayed logistic"
        ([1.0, 0.0], function (t, x)
            p = _candidate_rates(fit.name, fit.params, fit.doses, dose_at(t))
            growth = fit.name == "shared-growth delayed logistic" ? p.r * x[1] * (1 - x[1] / p.K) : p.r * x[1]
            [growth - x[2] * x[1], DEFAULT_GAMMA * (p.d - x[2])]
        end)
    elseif fit.name == "shared-K logistic"
        ([1.0], (t, x) -> begin
            p = _candidate_rates(fit.name, fit.params, fit.doses, dose_at(t))
            [p.r * x[1] * (1 - x[1] / p.K)]
        end)
    else
        ([1.0], (t, x) -> begin
            p = _candidate_rates(fit.name, fit.params, fit.doses, dose_at(t))
            [p.g * x[1]]
        end)
    end
    states = _integrate_at(times, initial, rhs; dt = dt)
    if startswith(fit.name, "paper two-state")
        sensitive, latent = states[:, 1], states[:, 2]
    elseif fit.name == "paper three-state shared growth"
        sensitive, latent = states[:, 1], states[:, 2] .+ states[:, 3]
    else
        sensitive, latent = states[:, 1], zeros(length(times))
    end
    return (time = times, states = states, sensitive = sensitive, latent = latent, total = sensitive .+ latent)
end

function _optimize_selected_control(fit; horizon = 13.0, intervals = 13, target_dose = 1.0, max_time = 20.0, seed = 2025)
    lo, hi = extrema(fit.doses)
    objective = raw -> begin
        schedule = _project_mean(Float64.(raw), target_dose, lo, hi)
        _simulate_selected_control(fit, schedule; horizon = horizon, dt = 0.05).total[end]
    end
    Random.seed!(seed)
    result = bboptimize(objective; SearchRange = (lo, hi), NumDimensions = intervals,
        Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = 120,
        MaxTime = max_time, TraceMode = :silent)
    return _project_mean(Vector{Float64}(result.archive_output.best_candidate), target_dose, lo, hi)
end

function _selected_parameter_screen(fit, curves; seed, draws = 4000)
    rng = MersenneTwister(seed)
    optimum = fit.sse
    accepted = Vector{Vector{Float64}}()
    spans = last.(fit.bounds) .- first.(fit.bounds)
    for scale in (0.03, 0.015, 0.0075, 0.00375)
        empty!(accepted)
        for _ in 1:draws
            candidate = clamp.(fit.params .+ scale .* spans .* randn(rng, length(fit.params)), first.(fit.bounds), last.(fit.bounds))
            sse = 0.0
            for (i, curve) in enumerate(curves)
                pred = _simulate_candidate(fit.name, candidate, fit.doses, fit.doses[i], curve.time).total
                sse += sum((pred .- curve.normalized) .^ 2)
            end
            sse <= 1.05 * max(optimum, eps()) && push!(accepted, candidate)
        end
        length(accepted) >= 100 && break
    end
    isempty(accepted) && push!(accepted, fit.params)
    matrix = reduce(hcat, accepted)'
    return DataFrame(
        parameter = fit.names,
        estimate = fit.params,
        lower_95 = [quantile(matrix[:, i], 0.025) for i in axes(matrix, 2)],
        median = [median(matrix[:, i]) for i in axes(matrix, 2)],
        upper_95 = [quantile(matrix[:, i], 0.975) for i in axes(matrix, 2)],
        retained = fill(size(matrix, 1), length(fit.params)),
    )
end

function _schedule_dose_at(t, schedule, horizon)
    return schedule[clamp(floor(Int, t / horizon * length(schedule)) + 1, 1, length(schedule))]
end

function _near_optimal_screen(best, times, observed; seed, draws = 3000)
    rng = MersenneTwister(seed)
    center = [best.rS, best.dS, best.alpha, best.rR / max(best.rS, eps()), best.dR / max(best.dS, eps())]
    optimum = _two_loss(center, times, observed)
    accepted = Vector{Vector{Float64}}()
    costs = Float64[]
    for _ in 1:draws
        candidate = clamp.(center .* exp.(0.18 .* randn(rng, 5)), [1e-4, 1e-4, 1e-5, 1e-4, 1e-4], [1.5, 2.0, 1.5, 1.0, 1.0])
        cost = _two_loss(candidate, times, observed)
        if cost <= 1.05 * max(optimum, eps())
            push!(accepted, candidate)
            push!(costs, cost)
        end
    end
    isempty(accepted) && push!(accepted, center)
    decoded = reduce(hcat, ([p.rS, p.dS, p.alpha, p.rR, p.dR] for p in _decode_two.(accepted)))'
    return (
        count = size(decoded, 1),
        lower = [quantile(decoded[:, j], 0.025) for j in 1:5],
        median = [median(decoded[:, j]) for j in 1:5],
        upper = [quantile(decoded[:, j], 0.975) for j in 1:5],
    )
end

_fmt(x; digits = 4) = isfinite(x) ? string(round(x; digits = digits)) : "NA"

function _table_html(df)
    io = IOBuffer()
    println(io, "<div class=\"table-wrap\"><table><thead><tr>")
    for name in names(df)
        println(io, "<th>", name, "</th>")
    end
    println(io, "</tr></thead><tbody>")
    for row in eachrow(df)
        println(io, "<tr>")
        for value in row
            text = value isa AbstractFloat ? _fmt(value) : string(value)
            println(io, "<td>", text, "</td>")
        end
        println(io, "</tr>")
    end
    println(io, "</tbody></table></div>")
    return String(take!(io))
end

function _write_report(path, fit_table, model_table, validation_table, control_table, near_table, figures)
    html = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Optimal Control One</title>
<script>window.MathJax={tex:{inlineMath:[["\\\\(","\\\\)"]],displayMath:[["\\\\[","\\\\]"]]}};</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
<style>
:root{--ink:#172033;--muted:#596477;--line:#cbd3dd;--soft:#f4f7fa;--accent:#2457d6}*{box-sizing:border-box}body{margin:0;font-family:Arial,sans-serif;color:var(--ink);line-height:1.5}header,main,footer{width:min(1180px,calc(100% - 40px));margin:auto}header{padding:28px 0 18px}a{color:var(--accent)}h1{font-size:34px;margin:18px 0 4px}h2{font-size:25px;margin:0 0 12px}h3{font-size:18px;margin:20px 0 8px}section{padding:25px 0;border-top:1px solid var(--line)}p{max-width:980px}.equation{background:var(--soft);border-left:4px solid var(--accent);padding:8px 16px;overflow:auto}.figure{margin:20px 0}.figure img{display:block;width:100%;height:auto;border:1px solid var(--line)}.caption,.note{color:var(--muted);font-size:14px}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;min-width:720px;font-size:14px}th,td{border:1px solid var(--line);padding:7px 9px;text-align:left}th{background:var(--soft)}.callout{border-left:4px solid #a43b56;padding:8px 16px;background:#fbf5f7}.toc{display:flex;flex-wrap:wrap;gap:8px;margin-top:15px}.toc a{border:1px solid var(--line);padding:6px 9px;text-decoration:none}footer{padding:24px 0 40px;color:var(--muted)}@media(max-width:650px){header,main,footer{width:min(100% - 24px,1180px)}h1{font-size:28px}}
</style></head><body><header><a href="../../../../index.html">&#8592; Back</a><h1>Optimal Control One</h1>
<p>Reproduction and A2780 adaptation of Gevertz et al. (2025), <a href="$PAPER_URL">Understanding therapeutic tolerance through a mathematical model of drug-induced resistance</a>.</p>
<nav class="toc"><a href="#paper">Paper models</a><a href="#adaptation">A2780 adaptation</a><a href="#fits">Fits</a><a href="#validation">Validation</a><a href="#control">Optimal control</a><a href="#interpretation">Interpretation</a></nav></header><main>
<section id="scope"><div class="callout"><strong>Interpretation boundary.</strong> The paper's resistant state is a latent drug-induced phenotype. In this A2780 adaptation, it is fitted to A2780Naive total-count trajectories and is not equated with conversion into the separately established A2780cis cell line. A2780cis is used only as an external resistant-line benchmark.</div></section>
<section id="paper"><h2>Models Used In The Paper</h2><h3>Two-population delayed-induction model</h3><div class="equation">\\[
\\frac{dS}{dt}=r_S S-d_S(1-e^{-\\gamma_1t})S-\\alpha(1-e^{-\\gamma_2t})S
\\]\\[
\\frac{dR}{dt}=r_R R+\\alpha(1-e^{-\\gamma_2t})S-d_R(1-e^{-\\gamma_1t})R
\\]</div><p>The paper fitted five dose-specific rates with \\(r_R\\le r_S\\), \\(d_R\\le d_S\\), fixed delay rates, an L1 objective, and initially \\(S(0)=1,R(0)=0\\).</p>
<h3>Three-population quiescent extension</h3><div class="equation">\\[
\\dot S=r_SS-d_S(1-e^{-\\gamma_1t})S-q(1-e^{-\\gamma_2t})S
\\]\\[
\\dot Q=q(1-e^{-\\gamma_2t})S-\\beta Q
\\]\\[
\\dot R=\\beta Q+r_RR-d_R(1-e^{-\\gamma_1t})R
\\]</div><p>The paper found that adding \\(Q\\) did not visibly improve fit and weakened practical identifiability. Both candidates are fitted below.</p>
<h3>Time-varying control model</h3><div class="equation">\\[
\\dot S=r_S(u)S-v_{d,S}S-v_\\alpha S,\\quad
\\dot R=r_R(u)R+v_\\alpha S-v_{d,R}R
\\]\\[
\\dot v_{d,S}=\\gamma_1[d_S(u)-v_{d,S}],\\quad
\\dot v_{d,R}=\\gamma_1[d_R(u)-v_{d,R}],\\quad
\\dot v_\\alpha=\\gamma_2[\\alpha(u)-v_\\alpha]
\\]</div></section>
<section id="adaptation"><h2>A2780 Data And Adaptation</h2><p>Inputs are the 30k A2780Naive treated-monoculture day averages at 0.67, 1.0, and 1.47 &micro;M cisplatin. Counts are normalized to each condition's first measured day. The delay is fixed at 0.24/day, the unit-converted value of the paper's 0.01/hour. Fits use bounded differential evolution and the paper's L1 objective. The observation horizon is 13 elapsed days.</p><p><strong>Dose-label correction:</strong> the historical processed-data splitter exchanged the endpoint folder names. Legacy folder <code>IC25</code> contains the 1.47 &micro;M source, <code>IC50</code> contains 1.0 &micro;M, and legacy folder <code>IC75</code> contains the 0.67 &micro;M source. This report applies the same canonical correction used by the main fitting workflow.</p><div class="figure"><img src="$(figures["data"])" alt="Normalized A2780Naive and A2780cis observations across cisplatin doses"><p class="caption">Measured A2780Naive trajectories used for fitting; A2780cis trajectories are shown as an external benchmark.</p></div></section>
<section id="fits"><h2>Dose-Specific Fits And Model Comparison</h2><div class="figure"><img src="$(figures["fits"])" alt="Two-population fits at three A2780 cisplatin doses"></div>$( _table_html(fit_table) )
<div class="figure"><img src="$(figures["parameters"])" alt="Fitted paper-model parameters by cisplatin dose"></div><h3>Two versus three populations</h3>$( _table_html(model_table) )<div class="figure"><img src="$(figures["comparison"])" alt="Two and three population model comparison"></div>
<h3>Near-optimal parameter screen</h3><p>Following the paper's 5% criterion, local parameter perturbations with L1 cost no more than 5% above the optimum are retained. Wide intervals indicate practical sensitivity even when an optimum is available.</p>$( _table_html(near_table) )</section>
<section id="validation"><h2>Held-Out Dose Validation</h2><p>The 1.0 &micro;M condition is withheld. Parameters are interpolated linearly in log-dose from fits at 0.67 and 1.47 &micro;M, matching the paper's validation logic.</p>$( _table_html(validation_table) )<div class="figure"><img src="$(figures["validation"])" alt="Held-out one micromolar validation"></div></section>
<section id="control"><h2>Constrained Optimal Control</h2><div class="equation">\\[
\\min_{u(t)} S(t_f)+R(t_f),\\qquad 0.67\\le u(t)\\le1.47,\\qquad \\int_0^{t_f}u(t)dt\\le13
\\]</div><p>The 13-day horizon and 13 &micro;M-day budget equal continuous 1 &micro;M exposure. The dose is piecewise constant over daily intervals. All comparison schedules are projected to the same dose budget.</p>$( _table_html(control_table) )<div class="figure"><img src="$(figures["control"])" alt="Dose schedules and resulting total population"></div><div class="figure"><img src="$(figures["composition"])" alt="Optimal-control sensitive and resistant composition"></div></section>
<section id="interpretation"><h2>Findings And Limits</h2><ul><li>The report tests whether the paper's minimal delayed-induction structure can reproduce the non-monotone A2780Naive dose trajectories.</li><li>The held-out 1 &micro;M prediction distinguishes interpolation from direct fitting.</li><li>The optimized schedule is a model-derived in vitro hypothesis, not a clinical dosing recommendation.</li><li>The latent \\(R\\) trajectory is not directly measured and must not be labeled A2780cis.</li><li>No extrapolation beyond day 14 is used.</li></ul></section>
</main><footer>Generated by Julia from repository data. Paper DOI: 10.1038/s41540-025-00511-3.</footer></body></html>"""
    mkpath(dirname(path))
    write(path, html)
    return path
end

function _selected_equation(name)
    name == "dose-specific exponential" && return "\\dot X=g(u)X"
    name == "shared-K logistic" && return "\\dot X=r(u)X\\left(1-\\frac{X}{K}\\right)"
    name == "shared-growth delayed loss" && return "\\dot X=rX-v_dX,\\qquad \\dot v_d=\\gamma[d(u)-v_d]"
    name == "shared-growth delayed logistic" && return "\\dot X=rX\\left(1-\\frac{X}{K}\\right)-v_dX,\\qquad \\dot v_d=\\gamma[d(u)-v_d]"
    startswith(name, "paper two-state") && return "\\dot S=r_SS-v_{d,S}S-v_\\alpha S,\\quad \\dot R=r_RR+v_\\alpha S-v_{d,R}R"
    return "\\dot S=r_SS-v_{d,S}S-v_qS,\\quad \\dot Q=v_qS-\\beta Q,\\quad \\dot R=\\beta Q+r_RR-v_{d,R}R"
end

function _write_selected_report(path, ranking, selected, parameter_table, validation_table, validation_note, control_table, figures)
    winner = selected.name
    equation = _selected_equation(winner)
    html = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Optimal Control One</title><script>window.MathJax={tex:{inlineMath:[["\\\\(","\\\\)"]],displayMath:[["\\\\[","\\\\]"]]}};</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script><style>
:root{--ink:#172033;--muted:#596477;--line:#cbd3dd;--soft:#f4f7fa;--accent:#2457d6;--good:#176b3a}*{box-sizing:border-box}body{margin:0;font-family:Arial,sans-serif;color:var(--ink);line-height:1.5}header,main,footer{width:min(1180px,calc(100% - 40px));margin:auto}header{padding:28px 0 18px}a{color:var(--accent)}h1{font-size:34px;margin:18px 0 4px}h2{font-size:25px;margin:0 0 12px}h3{font-size:18px;margin:20px 0 8px}section{padding:25px 0;border-top:1px solid var(--line)}p{max-width:980px}.equation{background:var(--soft);border-left:4px solid var(--accent);padding:8px 16px;overflow:auto}.figure{margin:20px 0}.figure img{display:block;width:100%;height:auto;border:1px solid var(--line)}.caption,.note{color:var(--muted);font-size:14px}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;min-width:720px;font-size:14px}th,td{border:1px solid var(--line);padding:7px 9px;text-align:left}th{background:var(--soft)}.callout{border-left:4px solid var(--good);padding:10px 16px;background:#f1f8f3}.warning{border-left-color:#a43b56;background:#fbf5f7}.toc{display:flex;flex-wrap:wrap;gap:8px;margin-top:15px}.toc a{border:1px solid var(--line);padding:6px 9px;text-decoration:none}footer{padding:24px 0 40px;color:var(--muted)}@media(max-width:650px){header,main,footer{width:min(100% - 24px,1180px)}h1{font-size:28px}}
</style></head><body><header><a href="../../../../index.html">&#8592; Back</a><h1>Optimal Control One</h1><p>A2780 model selection, validation, and constrained schedule optimization informed by <a href="$PAPER_URL">Gevertz et al. (2025)</a>.</p><nav class="toc"><a href="#selection">Model selection</a><a href="#fit">Selected fit</a><a href="#validation">Validation</a><a href="#control">Optimal control</a><a href="#limits">Limits</a></nav></header><main>
<section><div class="callout warning"><strong>Interpretation boundary.</strong> Latent drug-induced states are model variables, not measured A2780cis cells. A2780cis remains an external resistant-line benchmark.</div><p><strong>Corrected dose provenance:</strong> legacy folder IC25 contains 1.47 &micro;M, IC50 contains 1.0 &micro;M, and legacy IC75 contains 0.67 &micro;M.</p></section>
<section id="selection"><h2>Model Selection Before Optimal Control</h2><p>Seven candidate families are fitted jointly to the same 42 normalized A2780Naive observations across all three doses. BIC is \\(n\\log(\\mathrm{SSE}/n)+k\\log n\\): lower is better, and every fitted parameter across all doses contributes to \\(k\\). Models within 2 BIC units are treated as similarly supported; within that set, the simplest boundary-free model is carried forward. This prevents a negligible BIC advantage or a bound-limited parameter from deciding the control model.</p><div class="callout"><strong>Selected model:</strong> $winner. This is the simplest boundary-free member of the ΔBIC ≤ 2 supported set and is the only model carried into validation, sensitivity analysis, and optimal control.</div>$(_table_html(ranking))<div class="figure"><img src="$(figures["comparison"])" alt="Joint BIC ranking of all candidate models"></div></section>
<section id="fit"><h2>Selected Model Fit</h2><div class="equation">\\[$equation\\]</div><div class="figure"><img src="$(figures["fits"])" alt="Selected model fitted across all three doses"></div><h3>Selected parameter estimates and 5% near-optimal screen</h3><p>Shared parameters appear once; dose-indexed parameters are labeled in ascending dose order: [1] = 0.67, [2] = 1.0, [3] = 1.47 &micro;M. Wide or one-point intervals indicate practical identifiability limitations.</p>$(_table_html(parameter_table))</section>
<section id="validation"><h2>Held-Out Dose Validation</h2><p>After selecting the family, it is refitted using only 0.67 and 1.47 &micro;M. The 1.0 &micro;M trajectory is then predicted by log-dose interpolation without refitting to that trajectory.</p><div class="callout warning"><strong>Validation status:</strong> $validation_note</div>$(_table_html(validation_table))<div class="figure"><img src="$(figures["validation"])" alt="Held-out 1 micromolar validation for the selected model"></div></section>
<section id="control"><h2>Constrained Optimal Control Using The Selected Model</h2><div class="equation">\\[\\min_{u(t)} X(t_f),\\qquad 0.67\\le u(t)\\le1.47,\\qquad \\int_0^{t_f}u(t)dt=13\\]</div><p>The selected family is refitted to all three doses before schedule optimization. Daily piecewise-constant schedules share the same 13 &micro;M-day exposure budget. Because midpoint validation is weak, these schedules are exploratory model comparisons and should not be interpreted as experimentally validated optima.</p>$(_table_html(control_table))<div class="figure"><img src="$(figures["control"])" alt="Equal-budget control schedules and selected-model responses"></div><div class="figure"><img src="$(figures["composition"])" alt="Selected-model composition under the optimized schedule"></div></section>
<section id="limits"><h2>Interpretation</h2><ul><li>BIC selects the most supported candidate among those tested; it does not prove biological truth.</li><li>Model selection precedes every downstream calculation in this report.</li><li>Held-out error tests interpolation at 1.0 &micro;M; it is not included in that validation refit.</li><li>Near-optimal ranges expose parameter instability rather than presenting one optimizer result as exact.</li><li>The schedule is an in vitro model hypothesis, not a clinical recommendation.</li></ul></section>
</main><footer>Generated by Julia from repository data. Paper DOI: 10.1038/s41540-025-00511-3.</footer></body></html>"""
    mkpath(dirname(path)); write(path, html); return path
end

function _run_model_selected_workflow!(; start, fit_seconds, control_seconds, seed)
    package_root = isdir(joinpath(start, "src")) ? start : normpath(joinpath(start, ".."))
    repository_root = normpath(joinpath(package_root, "..", ".."))
    output_root = joinpath(package_root, "outputs", "csv", "optimal_control_one")
    figure_root = joinpath(output_root, "figures"); mkpath(figure_root)
    report_path = joinpath(package_root, "outputs", "reports", "optimal_control_one.html")
    docs_package = joinpath(repository_root, "docs", "Modeling_Approaches", "02_mechanical_automatic_package")
    docs_output = joinpath(docs_package, "outputs", "csv", "optimal_control_one")
    docs_report = joinpath(docs_package, "outputs", "reports", "optimal_control_one.html")
    stale_outputs = ("dose_specific_fits.csv", "model_comparison.csv", "near_optimal_ranges.csv")
    stale_figures = ("parameters_by_dose.png",)
    for filename in stale_outputs
        rm(joinpath(output_root, filename); force = true)
        rm(joinpath(docs_output, filename); force = true)
    end
    for filename in stale_figures
        rm(joinpath(figure_root, filename); force = true)
        rm(joinpath(docs_output, "figures", filename); force = true)
    end

    all_curves = Dict{Tuple{String, String}, Any}()
    for ic in keys(A2780_DOSE_MAP), lineage in ("A2780Naive", "A2780cis")
        all_curves[(lineage, ic)] = _read_curve(repository_root, lineage, ic)
    end
    doses = [A2780_DOSE_MAP[ic] for ic in A2780_SOURCE_ORDER]
    curves = [all_curves[("A2780Naive", ic)] for ic in A2780_SOURCE_ORDER]

    fits = [_fit_candidate(name, doses, curves; seed = seed + 100i, max_time = fit_seconds) for (i, name) in enumerate(CANDIDATE_MODELS)]
    eligible_bics = [fit.eligible ? fit.bic : Inf for fit in fits]
    best_bic = minimum(eligible_bics)
    supported = [fit for fit in fits if fit.eligible && fit.bic <= best_bic + 2.0]
    boundary_free = [fit for fit in supported if isempty(fit.boundary)]
    selection_pool = isempty(boundary_free) ? supported : boundary_free
    selected = sort(selection_pool; by = fit -> (length(fit.params), fit.bic))[1]
    ranking = sort(DataFrame([
        (model = fit.name, free_parameters = length(fit.params), SSE = fit.sse, RMSE = sqrt(fit.sse / 42), BIC = fit.bic,
         delta_BIC = fit.bic - best_bic, boundary_status = isempty(fit.boundary) ? "none" : join(fit.boundary, ", "),
         eligible = fit.eligible, selected = fit.name == selected.name) for fit in fits
    ]), :BIC)

    screen = _selected_parameter_screen(selected, curves; seed = seed + 900)
    trajectories = NamedTuple[]
    for (i, curve) in enumerate(curves)
        pred = _simulate_candidate(selected.name, selected.params, doses, doses[i], curve.time)
        condition = doses[i] == 0.67 ? "IC25" : doses[i] == 1.0 ? "IC50" : "IC75"
        for j in eachindex(curve.time)
            push!(trajectories, (dose_uM = doses[i], corrected_condition = condition, source_folder = A2780_SOURCE_ORDER[i], day = curve.day[j], observed = curve.normalized[j], sem = curve.sem_normalized[j], sensitive = pred.sensitive[j], latent_resistant = pred.latent[j], total = pred.total[j]))
        end
    end
    trajectory_table = DataFrame(trajectories)

    endpoint_fit = _fit_candidate(selected.name, [0.67, 1.47], [curves[1], curves[3]]; seed = seed + 1000, max_time = 2fit_seconds)
    heldout = curves[2]
    heldout_prediction = _simulate_candidate(endpoint_fit.name, endpoint_fit.params, endpoint_fit.doses, 1.0, heldout.time).total
    validation = DataFrame(day = heldout.day, observed = heldout.normalized, predicted = heldout_prediction, residual = heldout_prediction .- heldout.normalized)
    heldout_rmse = sqrt(mean(validation.residual .^ 2))
    validation_note = heldout_rmse <= 1.0 ? "Passes the prespecified exploratory threshold of RMSE ≤ 1 normalized unit." : "Weak interpolation: RMSE = $(round(heldout_rmse; digits = 3)) normalized units, above the exploratory threshold of 1. Control results are therefore hypothesis-generating."
    validation_summary = DataFrame(metric = ["selected family", "held-out L1", "held-out RMSE", "exploratory threshold"], value = [selected.name, string(sum(abs.(validation.residual))), string(heldout_rmse), "RMSE <= 1"])

    horizon, intervals, target = 13.0, 13, 1.0
    lo, hi = extrema(doses)
    optimized = _optimize_selected_control(selected; horizon = horizon, intervals = intervals, target_dose = target, max_time = control_seconds, seed = seed + 1100)
    schedules = Dict("Optimized" => optimized, "Constant 1 uM" => fill(target, intervals),
        "Front-loaded" => _two_level_schedule(intervals, lo, hi, target; front = true),
        "Back-loaded" => _two_level_schedule(intervals, lo, hi, target; front = false),
        "Pulsed" => _project_mean([isodd(i) ? hi : lo for i in 1:intervals], target, lo, hi))
    simulations = Dict{String, Any}(); control_rows = NamedTuple[]; control_trajectory_rows = NamedTuple[]
    for name in ("Optimized", "Constant 1 uM", "Front-loaded", "Back-loaded", "Pulsed")
        sim = _simulate_selected_control(selected, schedules[name]; horizon = horizon)
        simulations[name] = sim
        push!(control_rows, (schedule = name, selected_model = selected.name, final_total = sim.total[end], final_sensitive = sim.sensitive[end], final_latent_resistant = sim.latent[end], dose_AUC = mean(schedules[name]) * horizon, population_AUC = sum(sim.total) * (sim.time[2] - sim.time[1])))
        for j in eachindex(sim.time)
            push!(control_trajectory_rows, (schedule = name, time_day = sim.time[j], dose_uM = _schedule_dose_at(sim.time[j], schedules[name], horizon), sensitive = sim.sensitive[j], latent_resistant = sim.latent[j], total = sim.total[j]))
        end
    end
    control_table = DataFrame(control_rows); control_trajectories = DataFrame(control_trajectory_rows)

    default(; linewidth = 2, size = (1000, 620), legend = :best, gridalpha = 0.18, foreground_color_legend = nothing, background_color_legend = nothing)
    colors = Dict(0.67 => :steelblue, 1.0 => :darkorange, 1.47 => :firebrick)
    p_data = plot(xlabel = "Elapsed day", ylabel = "Normalized cell count", title = "A2780 30k treated monoculture data")
    for (i, ic) in enumerate(A2780_SOURCE_ORDER)
        n = curves[i]; c = all_curves[("A2780cis", ic)]; dose = doses[i]
        plot!(p_data, n.time, n.normalized; ribbon = n.sem_normalized, color = colors[dose], label = "Naive $(dose) uM")
        plot!(p_data, c.time, c.normalized; linestyle = :dash, color = colors[dose], label = "cis benchmark $(dose) uM")
    end
    data_figure = joinpath(figure_root, "a2780_data.png"); savefig(p_data, data_figure)

    panels = Any[]
    for dose in doses
        sub = filter(:dose_uM => ==(dose), trajectory_table)
        p = scatter(sub.day, sub.observed; yerror = sub.sem, color = :black, label = "observed", title = "$(dose) uM", xlabel = "Day", ylabel = "Normalized count")
        plot!(p, sub.day, sub.total; color = :black, label = "selected model total")
        maximum(sub.latent_resistant) > 1e-8 && plot!(p, sub.day, sub.sensitive; color = :seagreen, label = "S")
        maximum(sub.latent_resistant) > 1e-8 && plot!(p, sub.day, sub.latent_resistant; color = :crimson, linestyle = :dash, label = "latent non-S")
        push!(panels, p)
    end
    fit_figure = joinpath(figure_root, "dose_fits.png"); savefig(plot(panels...; layout = (1, 3), size = (1350, 430)), fit_figure)
    comparison_figure = joinpath(figure_root, "model_comparison.png")
    ranked_plot = sort(ranking, :BIC, rev = true)
    short_labels = replace.(ranked_plot.model,
        "shared-growth delayed logistic" => "Delayed logistic (shared r, K)",
        "shared-growth delayed loss" => "Delayed loss (shared r)",
        "shared-K logistic" => "Logistic (shared K)",
        "paper two-state shared growth" => "Paper 2-state (shared growth)",
        "paper three-state shared growth" => "Paper 3-state (shared growth)",
        "paper two-state dose-specific" => "Paper 2-state (dose-specific)",
        "dose-specific exponential" => "Exponential (dose-specific)")
    positions = collect(1:nrow(ranked_plot))
    p_bic = plot(; xlabel = "Delta BIC", yticks = (positions, short_labels), ylims = (0.5, nrow(ranked_plot) + 0.5),
        legend = false, title = "Joint model selection (lower is better)", size = (1100, 580))
    for i in positions
        color = ranked_plot.selected[i] ? :seagreen : :deepskyblue3
        plot!(p_bic, [0.0, ranked_plot.delta_BIC[i]], [i, i]; linewidth = 7, color = color)
        ranked_plot.selected[i] && annotate!(p_bic, ranked_plot.delta_BIC[i] + 2, i, text("selected", 9, :seagreen))
    end
    scatter!(p_bic, ranked_plot.delta_BIC, positions; markersize = 5, color = :black)
    savefig(p_bic, comparison_figure)
    validation_figure = joinpath(figure_root, "heldout_validation.png")
    p_val = scatter(validation.day, validation.observed; color = :black, label = "held-out observed", xlabel = "Day", ylabel = "Normalized count", title = "Selected-family prediction at 1.0 uM")
    plot!(p_val, validation.day, validation.predicted; color = :purple, label = "endpoint-only prediction"); savefig(p_val, validation_figure)
    control_figure = joinpath(figure_root, "control_comparison.png")
    p_dose = plot(xlabel = "Day", ylabel = "Dose (uM)", title = "Equal-budget schedules"); p_total = plot(xlabel = "Day", ylabel = "Normalized total", title = "Selected-model response")
    palette = Dict("Optimized" => :black, "Constant 1 uM" => :royalblue, "Front-loaded" => :seagreen, "Back-loaded" => :darkorange, "Pulsed" => :crimson)
    for name in keys(schedules)
        sim = simulations[name]; plot!(p_dose, sim.time, [_schedule_dose_at(t, schedules[name], horizon) for t in sim.time]; label = name, color = palette[name]); plot!(p_total, sim.time, sim.total; label = name, color = palette[name])
    end
    savefig(plot(p_dose, p_total; layout = (1, 2), size = (1300, 480)), control_figure)
    composition_figure = joinpath(figure_root, "optimized_composition.png"); sim_opt = simulations["Optimized"]
    p_comp = plot(sim_opt.time, sim_opt.total; color = :black, label = "total", xlabel = "Day", ylabel = "Normalized population", title = "Optimized selected-model state")
    maximum(sim_opt.latent) > 1e-8 && plot!(p_comp, sim_opt.time, sim_opt.sensitive; color = :seagreen, label = "S")
    maximum(sim_opt.latent) > 1e-8 && plot!(p_comp, sim_opt.time, sim_opt.latent; color = :crimson, linestyle = :dash, label = "latent non-S")
    savefig(p_comp, composition_figure)

    CSV.write(joinpath(output_root, "model_selection.csv"), ranking); CSV.write(joinpath(output_root, "selected_parameters.csv"), screen)
    CSV.write(joinpath(output_root, "fitted_trajectories.csv"), trajectory_table); CSV.write(joinpath(output_root, "heldout_validation.csv"), validation)
    CSV.write(joinpath(output_root, "control_metrics.csv"), control_table); CSV.write(joinpath(output_root, "control_trajectories.csv"), control_trajectories)
    figures = Dict("data" => "../csv/optimal_control_one/figures/a2780_data.png", "fits" => "../csv/optimal_control_one/figures/dose_fits.png", "comparison" => "../csv/optimal_control_one/figures/model_comparison.png", "validation" => "../csv/optimal_control_one/figures/heldout_validation.png", "control" => "../csv/optimal_control_one/figures/control_comparison.png", "composition" => "../csv/optimal_control_one/figures/optimized_composition.png")
    _write_selected_report(report_path, ranking, selected, screen, validation_summary, validation_note, control_table, figures)
    mkpath(dirname(docs_report)); mkpath(docs_output); cp(report_path, docs_report; force = true); cp(output_root, docs_output; force = true)
    return (report = report_path, docs_report = docs_report, output = output_root, selected_model = selected.name, ranking = ranking, control_table = control_table)
end

_asbool(value) = lowercase(string(value)) == "true"

function _top_stage_rows(ranking, status; n = 5)
    rows = NamedTuple[]
    group_column = :cell_line in propertynames(ranking) ? :cell_line : nothing
    groups = group_column === nothing ? [ranking] : [DataFrame(group) for group in groupby(ranking, group_column)]
    for group in groups
        ordered = sort(group, :bic)
        best = minimum(Float64.(ordered.bic))
        selected_model = group_column === nothing ? String(status.winning_model[1]) : begin
            key = String(group[1, group_column])
            status_row = filter(row -> String(row.cell_line) == key, status)
            isempty(status_row) ? "" : String(status_row.winning_model[1])
        end
        selected_pooling = group_column === nothing ? String(status.winning_pooling_mode[1]) : begin
            key = String(group[1, group_column]); status_row = filter(row -> String(row.cell_line) == key, status)
            isempty(status_row) ? "" : String(status_row.winning_pooling_mode[1])
        end
        for (candidate_index, row) in enumerate(eachrow(first(ordered, min(n, nrow(ordered)))))
            pooling = :pooling_mode in propertynames(ordered) ? String(row.pooling_mode) : ""
            push!(rows, (group = group_column === nothing ? "all" : String(row[group_column]), candidate = "M$candidate_index", model = String(row.model), pooling = pooling,
                free_parameters = :n_parameters in propertynames(ordered) ? Int(row.n_parameters) : (:parameter_count in propertynames(ordered) ? Int(row.parameter_count) : missing),
                BIC = Float64(row.bic), delta_BIC = Float64(row.bic) - best,
                boundary = :boundary_issue in propertynames(ordered) ? _asbool(row.boundary_issue) : false,
                eligible = :eligible_for_inheritance in propertynames(ordered) ? _asbool(row.eligible_for_inheritance) : true,
                selected = String(row.model) == selected_model && pooling == selected_pooling))
        end
    end
    return DataFrame(rows)
end

function _selected_overlay(overlay, status)
    parts = DataFrame[]
    if :cell_line in propertynames(overlay) && :cell_line in propertynames(status) && nrow(status) > 1
        for row in eachrow(status)
            part = filter(item -> String(item.cell_line) == String(row.cell_line) && String(item.model) == String(row.winning_model) && String(item.pooling_mode) == String(row.winning_pooling_mode), overlay)
            isempty(part) || push!(parts, part)
        end
    else
        model = String(status.winning_model[1]); pooling = String(status.winning_pooling_mode[1])
        part = filter(item -> String(item.model) == model && String(item.pooling_mode) == pooling, overlay)
        isempty(part) || push!(parts, part)
    end
    return isempty(parts) ? DataFrame() : vcat(parts...; cols = :union)
end

function _reproduction_summary(stage, selected)
    isempty(selected) && return DataFrame(stage = [stage], trajectories = [0], points = [0], median_normalized_RMSE = [Inf], worst_normalized_RMSE = [Inf])
    metadata = filter(name -> name in propertynames(selected), [:cell_line, :density, :dose, :mix, :component, :context])
    isempty(metadata) && (metadata = [:model])
    values = Float64[]
    for group in groupby(selected, metadata)
        observed = Float64.(group.observed); predicted = Float64.(group.predicted)
        scale = max(maximum(abs.(observed)), 1.0)
        push!(values, sqrt(mean((predicted .- observed) .^ 2)) / scale)
    end
    return DataFrame(stage = [stage], trajectories = [length(values)], points = [nrow(selected)], median_normalized_RMSE = [median(values)], worst_normalized_RMSE = [maximum(values)])
end

function _attach_candidate_errors(table, overlay)
    median_errors = Float64[]
    worst_errors = Float64[]
    point_counts = Int[]
    for candidate in eachrow(table)
        rows = filter(row -> String(row.model) == candidate.model &&
            String(row.pooling_mode) == candidate.pooling &&
            (candidate.group == "all" || !(:cell_line in propertynames(overlay)) || String(row.cell_line) == candidate.group), overlay)
        if isempty(rows)
            push!(median_errors, Inf); push!(worst_errors, Inf); push!(point_counts, 0)
            continue
        end
        metadata = filter(name -> name in propertynames(rows), [:cell_line, :density, :dose, :mix, :component, :context])
        isempty(metadata) && (metadata = [:model])
        errors = Float64[]
        for trajectory in groupby(rows, metadata)
            observed = Float64.(trajectory.observed); predicted = Float64.(trajectory.predicted)
            push!(errors, sqrt(mean((predicted .- observed) .^ 2)) / max(maximum(abs.(observed)), 1.0))
        end
        push!(median_errors, median(errors)); push!(worst_errors, maximum(errors)); push!(point_counts, nrow(rows))
    end
    result = copy(table)
    result[!, :median_nRMSE] = median_errors
    result[!, :worst_nRMSE] = worst_errors
    result[!, :points] = point_counts
    return result
end

function _stage4_endpoint_audit(ranking, overlays, bootstrap)
    endpoint = filter(row -> String(row.context) == "coculture" && isapprox(Float64(row.time), 14.0; atol = 1e-8), overlays)
    rows = NamedTuple[]
    for model in unique(String.(endpoint.model))
        model_rows = filter(row -> String(row.model) == model, endpoint)
        covered = 0; errors = Float64[]; matched = 0
        for row in eachrow(model_rows)
            lineage = String(row.component) == "sensitive" ? "A2780Naive" : "A2780cis"
            target = filter(item -> String(item.density) == String(row.density) && String(item.mix) == String(row.mix) && String(item.lineage) == lineage, bootstrap)
            isempty(target) && continue
            matched += 1
            prediction = Float64(row.predicted); observed = Float64(target.observed_mean[1])
            Float64(target.ci95_lower[1]) <= prediction <= Float64(target.ci95_upper[1]) && (covered += 1)
            push!(errors, abs(prediction - observed) / max(abs(observed), 1.0))
        end
        rank_row = filter(row -> String(row.model) == model, ranking)
        isempty(rank_row) && continue
        delta_bic = Float64(rank_row.delta_bic[1])
        push!(rows, (model = model, stage_winner = isapprox(delta_bic, 0.0; atol = 1e-8), BIC = Float64(rank_row.bic[1]), delta_BIC = delta_bic,
            free_parameters = Int(rank_row.n_parameters[1]), eligible = _asbool(rank_row.eligible_for_inheritance[1]),
            boundary = _asbool(rank_row.boundary_issue[1]), endpoints_inside_95CI = covered, endpoints_tested = matched,
            coverage = matched == 0 ? 0.0 : covered / matched, mean_relative_endpoint_error = isempty(errors) ? Inf : mean(errors)))
    end
    return sort(DataFrame(rows), [:coverage, :mean_relative_endpoint_error], rev = [true, false])
end

function _save_stage_diagnostic(selected, path, title)
    p = scatter(Float64.(selected.observed), Float64.(selected.predicted); alpha = 0.55, markersize = 4,
        xlabel = "Observed cell count", ylabel = "Predicted cell count", title = title, legend = false, size = (760, 620))
    limit = maximum(vcat(Float64.(selected.observed), Float64.(selected.predicted)))
    plot!(p, [0.0, limit], [0.0, limit]; color = :black, linestyle = :dash, linewidth = 2)
    savefig(p, path)
end

function _save_stage_bic(table, path, title)
    shown = first(sort(table, [:group, :delta_BIC]), min(10, nrow(table)))
    labels = [row.group == "all" ? row.candidate : string(row.group, " ", row.candidate) for row in eachrow(shown)]
    positions = collect(1:nrow(shown))
    p = plot(; xlabel = "Delta BIC within stage/group", yticks = (positions, labels), ylims = (0.5, nrow(shown) + 0.5), legend = false, title = title, size = (1100, 650))
    for i in positions
        color = shown.selected[i] ? :seagreen : :deepskyblue3
        plot!(p, [0.0, shown.delta_BIC[i]], [i, i]; linewidth = 7, color = color)
    end
    scatter!(p, shown.delta_BIC, positions; color = :black, markersize = 4)
    savefig(p, path)
end

function _schedule_protocol(schedule; horizon)
    windows = [TreatmentWindow(index - 1, index, dose, 1.0; label = "day $index") for (index, dose) in enumerate(schedule)]
    return TreatmentProtocol(windows; monitoring_interval = 1.0, allowable_doses = Float64[], horizon = horizon)
end

function _staged_control_simulation(config, schedule; horizon)
    protocol = _schedule_protocol(schedule; horizon = horizon)
    seed_scenario = A2780AdaptiveAdapter.build_a2780_scenario(config, protocol; density = "20k")
    tolerant_fraction = seed_scenario.parameters.initial_tolerant_fraction
    initial_state = [33.5, 33.5 * (1 - tolerant_fraction), 33.5 * tolerant_fraction, 0.0, 0.0, 0.0]
    scenario = A2780AdaptiveAdapter.build_a2780_scenario(config, protocol; density = "20k", initial_state = initial_state)
    return AdaptiveSimulationEngine.simulate_protocol(scenario)
end

function _optional_finite(row, name)
    name in propertynames(row) || return nothing
    value = row[name]
    return ismissing(value) || !isfinite(Float64(value)) ? nothing : Float64(value)
end

function _tournament_control_config(package_root, winner_package)
    config = A2780AdaptiveAdapter.load_a2780_adaptive_config(package_root)
    isdir(winner_package) || return config
    csv_root = joinpath(winner_package, "outputs", "csv")
    baseline_path = joinpath(csv_root, "monoculture_untreated", "untreated_group_baselines.csv")
    stage3_path = joinpath(csv_root, "coculture_untreated", "coculture_untreated_automatic_model_ranking.csv")
    stage4_path = joinpath(csv_root, "coculture_treated", "linked_treatment_model_ranking.csv")
    all(isfile, (baseline_path, stage3_path, stage4_path)) || return config
    baselines = CSV.read(baseline_path, DataFrame)
    baseline_json = [(
        cell_line = String(row.cell_line), density = String(row.density), model = String(row.best_model), pooling = String(row.pooling_mode),
        r = Float64(row.r), K = Float64(row.K),
        shape_parameter = ismissing(row.shape_parameter) ? nothing : String(row.shape_parameter),
        shape_value = _optional_finite(row, :shape_value), theta = _optional_finite(row, :theta),
        lag_time = _optional_finite(row, :lag_time), baranyi_q0 = _optional_finite(row, :baranyi_q0),
        adaptation_rate = _optional_finite(row, :adaptation_rate),
    ) for row in eachrow(baselines)]
    stage3 = A2780AdaptiveAdapter._ranked_candidates(3, stage3_path;
        identifiability_path = joinpath(csv_root, "coculture_untreated", "coculture_untreated_identifiability.csv"),
        evidence_level = "conditional tournament A2780 fit")
    stage4 = A2780AdaptiveAdapter._ranked_candidates(4, stage4_path;
        identifiability_path = joinpath(csv_root, "coculture_treated", "linked_treatment_identifiability.csv"),
        status_allowed = true, evidence_level = "conditional tournament A2780 linked fit")
    stage3_winner, _ = A2780AdaptiveAdapter._winner(stage3, 3)
    stage4_winner, _ = A2780AdaptiveAdapter._winner(stage4, 4)
    winners = merge(config.winners, (
        stage1 = baseline_json,
        stage3 = A2780AdaptiveAdapter._candidate_json(stage3_winner),
        stage4 = A2780AdaptiveAdapter._candidate_json(stage4_winner),
    ))
    return merge(config, (generated_from = "conditional model-tournament winner", winners = winners))
end

function _run_exploratory_staged_control(config, figure_path; max_time, seed)
    horizon = 14.0
    intervals = 14
    target_dose = 1.0
    lower_dose, upper_dose = 0.0, 1.47
    objective = raw -> begin
        schedule = _project_mean(Float64.(raw), target_dose, lower_dose, upper_dose)
        result = _staged_control_simulation(config, schedule; horizon = horizon)
        return Float64(result.observables[end][3])
    end
    Random.seed!(seed)
    optimization = bboptimize(objective; SearchRange = (lower_dose, upper_dose), NumDimensions = intervals,
        Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = 80,
        MaxTime = max_time, TraceMode = :silent)
    optimized = _project_mean(Vector{Float64}(optimization.archive_output.best_candidate), target_dose, lower_dose, upper_dose)
    schedules = Dict(
        "Optimized" => optimized,
        "Constant 1 uM" => fill(target_dose, intervals),
        "Front-loaded" => _two_level_schedule(intervals, lower_dose, upper_dose, target_dose; front = true),
        "Back-loaded" => _two_level_schedule(intervals, lower_dose, upper_dose, target_dose; front = false),
        "Pulsed" => _project_mean([isodd(index) ? upper_dose : lower_dose for index in 1:intervals], target_dose, lower_dose, upper_dose),
    )
    results = Dict(name => _staged_control_simulation(config, schedule; horizon = horizon) for (name, schedule) in schedules)
    rows = NamedTuple[]
    trajectories = NamedTuple[]
    for name in ("Optimized", "Constant 1 uM", "Front-loaded", "Back-loaded", "Pulsed")
        result = results[name]
        naive = Float64[observable[1] for observable in result.observables]
        cis = Float64[observable[2] for observable in result.observables]
        total = Float64[observable[3] for observable in result.observables]
        push!(rows, (schedule = name, final_total = total[end], final_A2780Naive = naive[end], final_total_A2780cis = cis[end],
            total_AUC = AdaptiveSimulationEngine._trapz(result.times, total), cis_AUC = AdaptiveSimulationEngine._trapz(result.times, cis),
            dose_AUC = AdaptiveSimulationEngine._step_auc(result.times, result.commanded_dose),
            daily_dose_uM = join(round.(schedules[name]; digits = 2), ", ")))
        for index in eachindex(result.times)
            push!(trajectories, (schedule = name, time_day = result.times[index], dose_uM = result.commanded_dose[index],
                effective_exposure = result.effective_exposure[index], A2780Naive = naive[index], total_A2780cis = cis[index], total = total[index]))
        end
    end
    palette = Dict("Optimized" => :black, "Constant 1 uM" => :royalblue, "Front-loaded" => :seagreen,
        "Back-loaded" => :darkorange, "Pulsed" => :crimson)
    schedule_order = ("Optimized", "Constant 1 uM", "Front-loaded", "Back-loaded", "Pulsed")
    population_max = maximum(observable[3] for result in values(results) for observable in result.observables)
    panels = Any[]
    for name in schedule_order
        result = results[name]
        naive = [observable[1] for observable in result.observables]
        cis = [observable[2] for observable in result.observables]
        total = [observable[3] for observable in result.observables]
        dose_panel = plot(result.times, result.commanded_dose; color = palette[name], linewidth = 3,
            seriestype = :steppost, fillrange = 0, fillalpha = 0.10, label = false,
            xlabel = "Day", ylabel = "Cisplatin (uM)", title = "$name dose",
            xlims = (0, horizon), ylims = (0, 1.55), yticks = 0:0.5:1.5)
        population_panel = plot(result.times, total; color = :black, linewidth = 3, label = "Total",
            xlabel = "Day", ylabel = "Cell count", title = "$name populations | final total $(round(Int, total[end]))",
            xlims = (0, horizon), ylims = (0, 1.06population_max), legend = :topleft)
        plot!(population_panel, result.times, naive; color = :firebrick, linewidth = 2, label = "A2780Naive")
        plot!(population_panel, result.times, cis; color = :royalblue, linewidth = 2, label = "Total A2780cis")
        scatter!(population_panel, [horizon], [total[end]]; color = :black, markersize = 4, label = false)
        push!(panels, dose_panel, population_panel)
    end
    savefig(plot(panels...; layout = (5, 2), size = (1320, 1780), margin = 5Plots.mm), figure_path)
    return (metrics = DataFrame(rows), trajectories = DataFrame(trajectories))
end

function _write_four_stage_report(path, tournament_section, stage_tables, summaries, endpoint_audit, readiness, control_table, figures)
    stage_names = ["Untreated Monoculture", "Treated Monoculture", "Untreated Coculture", "Treated Coculture"]
    sections = String[]
    for i in 1:4
        validation_text = i < 4 ? "This is a joint-fit reproduction audit across the available environments, not an independent holdout." : "Endpoint bootstrap coverage tests whether fitted day-14 predictions fall inside well-level 95% intervals; it is still an internal validation diagnostic rather than a new experiment."
        push!(sections, "<section id=\"stage$i\"><h2>Stage $i: $(stage_names[i])</h2><p>$validation_text</p>$(_table_html(stage_tables[i]))<div class=\"figure-grid\"><figure><img src=\"$(figures["stage$(i)_bic"])\" alt=\"Stage $i BIC ranking\"></figure><figure><img src=\"$(figures["stage$(i)_fit"])\" alt=\"Stage $i observed versus predicted\"></figure></div>$(_table_html(summaries[i]))</section>")
    end
    html = """<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Optimal Control One</title><style>
:root{--ink:#172033;--muted:#596477;--line:#cbd3dd;--soft:#f4f7fa;--accent:#2457d6;--warn:#a43b56;--good:#176b3a}*{box-sizing:border-box}body{margin:0;font-family:Arial,sans-serif;color:var(--ink);line-height:1.5}header,main,footer{width:min(1240px,calc(100% - 40px));margin:auto}header{padding:28px 0 18px}a{color:var(--accent)}h1{font-size:34px;margin:18px 0 4px}h2{font-size:25px;margin:0 0 12px}h3{font-size:19px;margin:22px 0 8px}section{padding:25px 0;border-top:1px solid var(--line)}p{max-width:1000px}.callout{border-left:4px solid var(--warn);padding:10px 16px;background:#fbf5f7}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;min-width:760px;font-size:14px}th,td{border:1px solid var(--line);padding:7px 9px;text-align:left;vertical-align:top}th{background:var(--soft)}.toc{display:flex;flex-wrap:wrap;gap:8px;margin-top:15px}.toc a{border:1px solid var(--line);padding:6px 9px;text-decoration:none}.figure-grid{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin:18px 0}.figure-grid figure{margin:0}.figure-grid img{display:block;width:100%;height:auto;border:1px solid var(--line)}.control-figure{margin:20px 0}.control-figure img{display:block;width:100%;height:auto;border:1px solid var(--line)}.tournament-summary{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));border:1px solid var(--line);margin:18px 0}.tournament-summary div{padding:12px;border-right:1px solid var(--line)}.tournament-summary div:last-child{border:0}.tournament-summary strong{display:block;font-size:18px}.path-table{min-width:1450px}.tournament details{border-top:1px solid var(--line);padding:12px 0}.tournament small{color:var(--muted)}footer{padding:24px 0 40px;color:var(--muted)}@media(max-width:800px){header,main,footer{width:min(100% - 24px,1240px)}.figure-grid{grid-template-columns:1fr}.tournament-summary{grid-template-columns:1fr 1fr}.tournament-summary div:nth-child(2){border-right:0}h1{font-size:28px}}
</style></head><body><header><a href="../../../../index.html">&#8592; Back</a><h1>Optimal Control One</h1><p>Conditional model selection, four-stage validation, and exploratory optimal control for A2780.</p><nav class="toc"><a href="#tournament">Model tournament</a><a href="#stage1">Stage 1</a><a href="#stage2">Stage 2</a><a href="#stage3">Stage 3</a><a href="#stage4">Stage 4</a><a href="#endpoint">Endpoint audit</a><a href="#control">Optimal control</a></nav></header><main>
$(tournament_section)
<section><div class="callout"><strong>Why the previous validation failed:</strong> the 1.0 &micro;M trajectory is not intermediate between 0.67 and 1.47 &micro;M, so endpoint-only dose interpolation is unsupported. This report now audits the complete inherited four-stage workflow and does not use that failed interpolation as evidence.</div></section>
<section><h2>How To Read The Audit</h2><p><strong>BIC</strong> balances fit against the number of free parameters; lower is better. <strong>Delta BIC</strong> is measured from the lowest BIC in the same stage and lineage group. <strong>Median and worst nRMSE</strong> divide trajectory RMSE by that trajectory's largest observed count, making shape errors comparable across seeding densities. These nRMSE values describe reproduction of fitted data and are not held-out prediction errors. <strong>Boundary</strong> marks a parameter at or near an allowed limit. <strong>Eligible</strong> means the staged workflow permits inheritance; an ineligible low-BIC model is shown for comparison but is not carried forward.</p></section>
$(join(sections, "\n"))
<section id="endpoint"><h2>Stage 4 Candidate Endpoint Audit</h2><p>Every linked treated-coculture candidate is checked against the same 12 day-14 bootstrap intervals. High endpoint coverage cannot rescue an ineligible inheritance model, but it can reveal whether the BIC winner generalizes poorly.</p>$(_table_html(first(endpoint_audit, min(12, nrow(endpoint_audit)))))</section>
<section id="control"><h2>Exploratory Optimal Control</h2><div class="callout"><strong>The optimizer was run despite incomplete validation.</strong> The result is the mathematical optimum found for the selected fitted model under the stated constraints. It is a hypothesis for comparing schedules, not a validated experimental optimum.</div><p>The search minimizes the day-14 total population for the tournament-selected staged model. Dose is piecewise constant in 14 daily intervals, constrained to 0-1.47 &micro;M, and every schedule has the same 14 &micro;M-day commanded-dose budget as constant 1 &micro;M treatment. The initial condition is the measured 20k, 50:50 day-zero composition: 33.5 A2780Naive and 33.5 total A2780cis cells in the image-analysis scale.</p>$(_table_html(control_table))<figure class="control-figure"><img src="$(figures["control"])" alt="Five rows pairing each cisplatin schedule with its corresponding total, A2780Naive, and total A2780cis trajectories"></figure><h3>Why The Parameters Are Not Yet Reliable</h3><ul><li>The selected Stage 4 model covers only 3 of 12 treated-coculture day-14 bootstrap intervals.</li><li>No complete density, mixture, dose, or experiment was prospectively reserved before model selection.</li><li>Several inherited parameters sit at or near fitting limits, including treatment-effect amplitudes, activation rate, onset time, and the Stage 3 loss model.</li><li>The narrow residual-bootstrap ranges show repeated convergence to the same constrained solution; they do not repair the poor endpoint coverage or test structural model error.</li><li>The 1.0 &micro;M response is not intermediate between 0.67 and 1.47 &micro;M, so a smooth dose-response assumption can mis-rank schedules that spend time between measured doses.</li></ul><h3>Validation Status</h3>$(_table_html(readiness))<h3>Alternatives To Test</h3><ul><li>Optimize each well-supported complete model path and report schedule disagreement as model uncertainty.</li><li>Restrict controls to the measured dose set 0, 0.67, 1.0, and 1.47 &micro;M instead of continuous interpolation.</li><li>Use robust optimization against bootstrap and alternate-model ensembles rather than one parameter vector.</li><li>Collect an independent mixture, density, or treatment sequence and compare predictions without refitting.</li></ul></section>
</main><footer>Generated from the staged Julia ranking, inheritance, overlay, and bootstrap artifacts.</footer></body></html>"""
    mkpath(dirname(path)); write(path, html); return path
end

function _run_four_stage_workflow!(; start, control_seconds, seed)
    package_root = isdir(joinpath(start, "src")) ? start : normpath(joinpath(start, ".."))
    repository_root = normpath(joinpath(package_root, "..", "..")); csv_root = joinpath(package_root, "outputs", "csv")
    output_root = joinpath(csv_root, "optimal_control_one"); figure_root = joinpath(output_root, "figures"); mkpath(figure_root)
    report_path = joinpath(package_root, "outputs", "reports", "optimal_control_one.html")
    docs_root = joinpath(repository_root, "docs", "Modeling_Approaches", "02_mechanical_automatic_package")
    docs_output = joinpath(docs_root, "outputs", "csv", "optimal_control_one"); docs_report = joinpath(docs_root, "outputs", "reports", "optimal_control_one.html")

    ranking_paths = [joinpath(csv_root, "monoculture_untreated", "monoculture_untreated_pooling_model_ranking.csv"), joinpath(csv_root, "monoculture_treated", "monoculture_treated_pooling_model_ranking.csv"), joinpath(csv_root, "coculture_untreated", "coculture_untreated_pooling_model_ranking.csv"), joinpath(csv_root, "coculture_treated", "linked_treatment_model_ranking.csv")]
    status_paths = [joinpath(csv_root, "monoculture_untreated", "monoculture_untreated_pooling_status.csv"), joinpath(csv_root, "monoculture_treated", "monoculture_treated_pooling_status.csv"), joinpath(csv_root, "coculture_untreated", "coculture_untreated_pooling_status.csv"), joinpath(csv_root, "coculture_treated", "linked_treatment_status.csv")]
    overlay_paths = [joinpath(csv_root, "monoculture_untreated", "figures", "monoculture_untreated_pooling_overlays.csv"), joinpath(csv_root, "monoculture_treated", "figures", "monoculture_treated_joint_dose_overlays.csv"), joinpath(csv_root, "coculture_untreated", "figures", "coculture_untreated_joint_overlays.csv"), joinpath(csv_root, "coculture_treated", "figures", "linked_treatment_combined_overlays.csv")]
    rankings = CSV.read.(ranking_paths, Ref(DataFrame)); statuses = CSV.read.(status_paths, Ref(DataFrame)); overlays = CSV.read.(overlay_paths, Ref(DataFrame))
    stage_tables = [_attach_candidate_errors(_top_stage_rows(rankings[i], statuses[i]), overlays[i]) for i in 1:4]
    for (stage, table) in enumerate(stage_tables)
        insertcols!(table, 1, :stage => fill(stage, nrow(table)))
    end
    selected = [_selected_overlay(overlays[i], statuses[i]) for i in 1:4]
    summaries = [_reproduction_summary(i, selected[i]) for i in 1:4]
    bootstrap = CSV.read(joinpath(csv_root, "coculture_treated", "linked_treatment_endpoint_bootstrap.csv"), DataFrame)
    endpoint_audit = _stage4_endpoint_audit(rankings[4], overlays[4], bootstrap)
    winner = String(statuses[4].winning_model[1]); winner_audit = only(filter(row -> row.model == winner, endpoint_audit))
    heldout_available = false; coverage_pass = winner_audit.endpoints_inside_95CI >= 10; boundary_pass = !winner_audit.boundary
    inheritance_pass = all(all(_asbool.(status.inheritance_allowed)) for status in statuses)
    inherited_boundary_labels = String[]
    for (stage, table) in enumerate(stage_tables), row in eachrow(table)
        row.selected && row.boundary && push!(inherited_boundary_labels, "Stage $stage: $(row.model)")
    end
    inherited_boundary_pass = isempty(inherited_boundary_labels)
    inherited_boundary_evidence = inherited_boundary_pass ? "no selected stage winner has a boundary flag" : join(inherited_boundary_labels, "; ")
    readiness = DataFrame(criterion = ["All four inherited winners eligible", "All inherited winners boundary-free", "Independent held-out environment", "Stage 4 endpoint coverage >= 10/12", "Stage 4 winner boundary-free", "Overall validation-ready"],
        passed = [inheritance_pass, inherited_boundary_pass, heldout_available, coverage_pass, boundary_pass, false],
        evidence = ["staged status exports", inherited_boundary_evidence, "no predeclared held-out environment was exported", "$(winner_audit.endpoints_inside_95CI)/$(winner_audit.endpoints_tested) endpoints inside bootstrap intervals", winner_audit.boundary ? "winner is boundary-limited" : "no boundary flag", "exploratory optimizer run; model remains unvalidated"])
    CSV.write(joinpath(output_root, "four_stage_model_audit.csv"), vcat(stage_tables...; cols = :union)); CSV.write(joinpath(output_root, "four_stage_reproduction_summary.csv"), vcat(summaries...)); CSV.write(joinpath(output_root, "stage4_endpoint_model_audit.csv"), endpoint_audit); CSV.write(joinpath(output_root, "control_readiness.csv"), readiness)
    figures = Dict{String,String}()
    for i in 1:4
        bic_path = joinpath(figure_root, "stage$(i)_bic.svg"); fit_path = joinpath(figure_root, "stage$(i)_fit.png")
        _save_stage_bic(stage_tables[i], bic_path, "Stage $i model comparison"); _save_stage_diagnostic(selected[i], fit_path, "Stage $i selected-model reproduction")
        figures["stage$(i)_bic"] = "../csv/optimal_control_one/figures/stage$(i)_bic.svg"; figures["stage$(i)_fit"] = "../csv/optimal_control_one/figures/stage$(i)_fit.png"
    end
    tournament_output = joinpath(package_root, "outputs", "model_path_tournament")
    tournament_section = ModelPathTournament.tournament_report_section(tournament_output)
    tournament_ranking = CSV.read(joinpath(tournament_output, "complete_path_ranking.csv"), DataFrame)
    winner_workspace = String(tournament_ranking.workspace[1])
    winner_package = joinpath(winner_workspace, "Modeling_Approaches", "02_mechanical_automatic_package")
    control_config_root = isdir(winner_package) ? winner_package : package_root
    control_config = _tournament_control_config(package_root, control_config_root)
    control_path = joinpath(figure_root, "exploratory_staged_control.png")
    control = _run_exploratory_staged_control(control_config, control_path; max_time = control_seconds, seed = seed)
    CSV.write(joinpath(output_root, "control_metrics.csv"), control.metrics)
    CSV.write(joinpath(output_root, "control_trajectories.csv"), control.trajectories)
    figures["control"] = "../csv/optimal_control_one/figures/exploratory_staged_control.png"
    _write_four_stage_report(report_path, tournament_section, stage_tables, summaries, endpoint_audit, readiness, control.metrics, figures)
    stale = ("fitted_trajectories.csv", "heldout_validation.csv", "model_selection.csv", "selected_parameters.csv")
    stale_figures = ("a2780_data.png", "control_comparison.png", "dose_fits.png", "heldout_validation.png", "model_comparison.png", "optimized_composition.png",
        "stage1_bic.png", "stage2_bic.png", "stage3_bic.png", "stage4_bic.png")
    for name in stale; rm(joinpath(output_root, name); force = true); rm(joinpath(docs_output, name); force = true); end
    for name in stale_figures; rm(joinpath(figure_root, name); force = true); rm(joinpath(docs_output, "figures", name); force = true); end
    mkpath(dirname(docs_report)); mkpath(docs_output); cp(report_path, docs_report; force = true); cp(output_root, docs_output; force = true)
    return (report = report_path, docs_report = docs_report, output = output_root, control_ready = false, control_run = true, readiness = readiness, control_table = control.metrics)
end

function run_optimal_control_one!(; start = @__DIR__, fit_seconds = 8.0, control_seconds = 20.0, seed = 5113)
    return _run_four_stage_workflow!(; start = start, control_seconds = control_seconds, seed = seed)
    package_root = isdir(joinpath(start, "src")) ? start : normpath(joinpath(start, ".."))
    repository_root = normpath(joinpath(package_root, "..", ".."))
    output_root = joinpath(package_root, "outputs", "csv", "optimal_control_one")
    figure_root = joinpath(output_root, "figures")
    report_path = joinpath(package_root, "outputs", "reports", "optimal_control_one.html")
    docs_package = joinpath(repository_root, "docs", "Modeling_Approaches", "02_mechanical_automatic_package")
    docs_output = joinpath(docs_package, "outputs", "csv", "optimal_control_one")
    docs_report = joinpath(docs_package, "outputs", "reports", "optimal_control_one.html")
    mkpath(figure_root)

    curves = Dict{Tuple{String, String}, Any}()
    for ic in keys(A2780_DOSE_MAP), lineage in ("A2780Naive", "A2780cis")
        curves[(lineage, ic)] = _read_curve(repository_root, lineage, ic)
    end

    two_fits = Dict{Float64, PaperParameters}()
    three_fits = Dict{Float64, Vector{Float64}}()
    fit_rows = NamedTuple[]
    comparison_rows = NamedTuple[]
    near_rows = NamedTuple[]
    trajectory_rows = NamedTuple[]
    two_bounds = [(1e-4, 1.5), (1e-4, 2.0), (1e-5, 1.5), (1e-4, 1.0), (1e-4, 1.0)]
    three_bounds = [(1e-4, 1.5), (1e-4, 2.0), (1e-5, 1.5), (1e-4, 1.5), (1e-4, 1.0), (1e-4, 1.0)]

    for (index, ic) in enumerate(A2780_SOURCE_ORDER)
        dose = A2780_DOSE_MAP[ic]
        corrected_condition = dose == 0.67 ? "IC25" : dose == 1.0 ? "IC50" : "IC75"
        curve = curves[("A2780Naive", ic)]
        z2, l1_two = _fit_model(z -> _two_loss(z, curve.time, curve.normalized), two_bounds; seed = seed + index, max_time = fit_seconds)
        p2 = _decode_two(z2)
        two_fits[dose] = p2
        pred2 = simulate_two_population(p2, curve.time)
        total2 = pred2[:, 1] .+ pred2[:, 2]
        rss2 = sum((total2 .- curve.normalized) .^ 2)
        bic2 = length(curve.time) * log(max(rss2, 1e-12) / length(curve.time)) + 5 * log(length(curve.time))
        push!(fit_rows, (dose_uM = dose, corrected_condition = corrected_condition, source_folder = ic, rS_per_day = p2.rS, dS_per_day = p2.dS, alpha_per_day = p2.alpha, rR_per_day = p2.rR, dR_per_day = p2.dR, L1 = l1_two, BIC = bic2))
        for i in eachindex(curve.time)
            push!(trajectory_rows, (dose_uM = dose, corrected_condition = corrected_condition, source_folder = ic, day = curve.day[i], observed = curve.normalized[i], sem = curve.sem_normalized[i], sensitive = pred2[i, 1], latent_resistant = pred2[i, 2], total = total2[i]))
        end

        z3, l1_three = _fit_model(z -> _three_loss(z, curve.time, curve.normalized), three_bounds; seed = seed + 100 + index, max_time = fit_seconds)
        p3 = _decode_three(z3)
        three_fits[dose] = p3
        pred3 = simulate_three_population(p3, curve.time)
        total3 = vec(sum(pred3; dims = 2))
        rss3 = sum((total3 .- curve.normalized) .^ 2)
        bic3 = length(curve.time) * log(max(rss3, 1e-12) / length(curve.time)) + 6 * log(length(curve.time))
        push!(comparison_rows, (dose_uM = dose, corrected_condition = corrected_condition, source_folder = ic, model = "two population", free_parameters = 5, L1 = l1_two, BIC = bic2))
        push!(comparison_rows, (dose_uM = dose, corrected_condition = corrected_condition, source_folder = ic, model = "three population", free_parameters = 6, L1 = l1_three, BIC = bic3))

        screen = _near_optimal_screen(p2, curve.time, curve.normalized; seed = seed + 200 + index)
        for (j, name) in enumerate(("rS", "dS", "alpha", "rR", "dR"))
            push!(near_rows, (dose_uM = dose, corrected_condition = corrected_condition, source_folder = ic, parameter = name, retained = screen.count, lower_95 = screen.lower[j], median = screen.median[j], upper_95 = screen.upper[j]))
        end
    end

    fit_table = DataFrame(fit_rows)
    model_table = DataFrame(comparison_rows)
    near_table = DataFrame(near_rows)
    trajectories = DataFrame(trajectory_rows)

    training = Dict(0.67 => two_fits[0.67], 1.47 => two_fits[1.47])
    predicted_one = interpolate_parameters(1.0, training)
    heldout = curves[("A2780Naive", "IC50")]
    heldout_states = simulate_two_population(predicted_one, heldout.time)
    heldout_total = heldout_states[:, 1] .+ heldout_states[:, 2]
    validation = DataFrame(
        day = heldout.day,
        observed = heldout.normalized,
        predicted = heldout_total,
        residual = heldout_total .- heldout.normalized,
    )
    validation_summary = DataFrame(metric = ["held-out L1", "held-out RMSE"], value = [sum(abs.(validation.residual)), sqrt(mean(validation.residual .^ 2))])

    horizon, intervals, target = 13.0, 13, 1.0
    lo, hi = extrema(collect(keys(two_fits)))
    optimized = _optimize_control(two_fits; horizon = horizon, intervals = intervals, target_dose = target, max_time = control_seconds, seed = seed + 500)
    schedules = Dict(
        "Optimized" => optimized,
        "Constant 1 uM" => fill(target, intervals),
        "Front-loaded" => _two_level_schedule(intervals, lo, hi, target; front = true),
        "Back-loaded" => _two_level_schedule(intervals, lo, hi, target; front = false),
        "Pulsed" => _project_mean([isodd(i) ? hi : lo for i in 1:intervals], target, lo, hi),
    )
    control_rows = NamedTuple[]
    control_trajectory_rows = NamedTuple[]
    simulations = Dict{String, Any}()
    for name in ("Optimized", "Constant 1 uM", "Front-loaded", "Back-loaded", "Pulsed")
        schedule = schedules[name]
        sim = _simulate_control(two_fits, schedule; horizon = horizon)
        simulations[name] = sim
        push!(control_rows, (schedule = name, final_total = sim.total[end], final_sensitive = sim.states[end, 1], final_latent_resistant = sim.states[end, 2], dose_AUC = mean(schedule) * horizon, population_AUC = sum(sim.total) * (sim.time[2] - sim.time[1])))
        for i in eachindex(sim.time)
            push!(control_trajectory_rows, (schedule = name, time_day = sim.time[i], dose_uM = _schedule_dose_at(sim.time[i], schedule, horizon), sensitive = sim.states[i, 1], latent_resistant = sim.states[i, 2], total = sim.total[i]))
        end
    end
    control_table = DataFrame(control_rows)
    control_trajectories = DataFrame(control_trajectory_rows)

    default(; linewidth = 2, size = (1000, 620), legend = :best, gridalpha = 0.18, foreground_color_legend = nothing, background_color_legend = nothing)
    colors = Dict(0.67 => :steelblue, 1.0 => :darkorange, 1.47 => :firebrick)
    p_data = plot(xlabel = "Elapsed day", ylabel = "Normalized cell count", title = "A2780 30k treated monoculture data")
    for ic in A2780_SOURCE_ORDER
        dose = A2780_DOSE_MAP[ic]
        n = curves[("A2780Naive", ic)]
        c = curves[("A2780cis", ic)]
        plot!(p_data, n.time, n.normalized; ribbon = n.sem_normalized, color = colors[dose], label = "Naive $(dose) uM")
        plot!(p_data, c.time, c.normalized; linestyle = :dash, color = colors[dose], label = "cis benchmark $(dose) uM")
    end
    data_figure = joinpath(figure_root, "a2780_data.png"); savefig(p_data, data_figure)

    panels = Any[]
    for ic in A2780_SOURCE_ORDER
        dose = A2780_DOSE_MAP[ic]
        subset = filter(:dose_uM => ==(dose), trajectories)
        p = plot(subset.day, subset.observed; ribbon = subset.sem, seriestype = :scatter, color = :black, label = "observed", title = "$(dose) uM", xlabel = "Day", ylabel = "Normalized count")
        plot!(p, subset.day, subset.total; color = :black, label = "S + R")
        plot!(p, subset.day, subset.sensitive; color = :seagreen, label = "S")
        plot!(p, subset.day, subset.latent_resistant; color = :crimson, linestyle = :dash, label = "latent R")
        push!(panels, p)
    end
    fit_figure = joinpath(figure_root, "dose_fits.png"); savefig(plot(panels...; layout = (1, 3), size = (1350, 430)), fit_figure)

    parameter_figure = joinpath(figure_root, "parameters_by_dose.png")
    p_params = plot(layout = (2, 3), size = (1200, 700))
    for (j, field) in enumerate((:rS_per_day, :dS_per_day, :alpha_per_day, :rR_per_day, :dR_per_day))
        plot!(p_params[j], fit_table.dose_uM, fit_table[!, field]; marker = :circle, label = string(field), xlabel = "Dose (uM)", ylabel = "1/day")
    end
    savefig(p_params, parameter_figure)

    comparison_figure = joinpath(figure_root, "model_comparison.png")
    labels = ["$(row.dose_uM)\n$(row.model == "two population" ? "2-pop" : "3-pop")" for row in eachrow(model_table)]
    savefig(bar(labels, model_table.BIC; xlabel = "Dose and model", ylabel = "BIC", legend = false, title = "Complexity-adjusted comparison", xrotation = 25), comparison_figure)

    validation_figure = joinpath(figure_root, "heldout_validation.png")
    p_val = scatter(validation.day, validation.observed; color = :black, label = "held-out observed", xlabel = "Day", ylabel = "Normalized count", title = "1.0 uM prediction from 0.67 and 1.47 uM fits")
    plot!(p_val, validation.day, validation.predicted; color = :purple, label = "log-dose interpolated prediction")
    savefig(p_val, validation_figure)

    control_figure = joinpath(figure_root, "control_comparison.png")
    p_dose = plot(xlabel = "Day", ylabel = "Dose (uM)", title = "Equal-budget schedules")
    p_total = plot(xlabel = "Day", ylabel = "Normalized total", title = "Model response")
    palette = Dict("Optimized" => :black, "Constant 1 uM" => :royalblue, "Front-loaded" => :seagreen, "Back-loaded" => :darkorange, "Pulsed" => :crimson)
    for name in keys(schedules)
        sim = simulations[name]
        plot!(p_dose, sim.time, [_schedule_dose_at(t, schedules[name], horizon) for t in sim.time]; label = name, color = palette[name])
        plot!(p_total, sim.time, sim.total; label = name, color = palette[name])
    end
    savefig(plot(p_dose, p_total; layout = (1, 2), size = (1300, 480)), control_figure)

    composition_figure = joinpath(figure_root, "optimized_composition.png")
    sim_opt = simulations["Optimized"]
    p_comp = plot(sim_opt.time, sim_opt.total; color = :black, label = "S + R", xlabel = "Day", ylabel = "Normalized population", title = "Optimized latent composition")
    plot!(p_comp, sim_opt.time, sim_opt.states[:, 1]; color = :seagreen, label = "S")
    plot!(p_comp, sim_opt.time, sim_opt.states[:, 2]; color = :crimson, linestyle = :dash, label = "latent R")
    savefig(p_comp, composition_figure)

    CSV.write(joinpath(output_root, "dose_specific_fits.csv"), fit_table)
    CSV.write(joinpath(output_root, "fitted_trajectories.csv"), trajectories)
    CSV.write(joinpath(output_root, "model_comparison.csv"), model_table)
    CSV.write(joinpath(output_root, "near_optimal_ranges.csv"), near_table)
    CSV.write(joinpath(output_root, "heldout_validation.csv"), validation)
    CSV.write(joinpath(output_root, "control_metrics.csv"), control_table)
    CSV.write(joinpath(output_root, "control_trajectories.csv"), control_trajectories)

    figures = Dict(
        "data" => "../csv/optimal_control_one/figures/a2780_data.png",
        "fits" => "../csv/optimal_control_one/figures/dose_fits.png",
        "parameters" => "../csv/optimal_control_one/figures/parameters_by_dose.png",
        "comparison" => "../csv/optimal_control_one/figures/model_comparison.png",
        "validation" => "../csv/optimal_control_one/figures/heldout_validation.png",
        "control" => "../csv/optimal_control_one/figures/control_comparison.png",
        "composition" => "../csv/optimal_control_one/figures/optimized_composition.png",
    )
    _write_report(report_path, fit_table, model_table, validation_summary, control_table, near_table, figures)

    mkpath(dirname(docs_report)); mkpath(docs_output)
    cp(report_path, docs_report; force = true)
    cp(output_root, docs_output; force = true)
    return (report = report_path, docs_report = docs_report, output = output_root, fit_table = fit_table, control_table = control_table)
end

end
