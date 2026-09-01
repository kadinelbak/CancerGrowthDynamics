module OptimalControlOne

using ..IOUtils: TREATED_MONOCULTURE_IC_DOSE_MAP
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

function run_optimal_control_one!(; start = @__DIR__, fit_seconds = 8.0, control_seconds = 20.0, seed = 5113)
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
