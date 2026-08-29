module ReducedStageReport

using CSV
using DataFrames
using Plots

using ..IOUtils
using ..FitWorkflows

export render_reduced_stage_report_html

const REPORT_NAME = "a2780_reduced_stage_comparison.html"
const ROUTE_A = "Stage 1 -> Stage 3 -> Stage 4"
const ROUTE_B = "Stage 1 -> Stage 2 -> Stage 4"
const MODEL_LABELS = Dict(
    "dual_transit_load_scaled" => "Transit damage with load-scaled treatment",
    "dual_time_decay_kill" => "Time-decaying direct loss",
    "dual_constant_kill" => "Constant direct loss",
    "dual_transit_damage" => "Transit damage",
    "dual_transit_competitor_scaled" => "Competitor-scaled transit damage",
    "dual_delayed_ramp_kill" => "Delayed ramp loss",
    "sensitive_tolerant_transition" => "Sensitive/tolerant transition",
    "stage2_direct_context_scaling" => "Inherited treatment with coculture amplitude scaling",
    "stage2_plus_asymmetric_competition" => "Inherited treatment with fitted asymmetric competition",
    "stage2_plus_asymmetric_competition_death" => "Inherited treatment with fitted asymmetric competition and loss",
)
const ROUTE_A_EQUATIONS = raw"""\[\frac{dN}{dt}=G_N\!\left(N,N+\alpha_{NC}C\right)-d_NN-E_N(t,N,C)N\]
\[\frac{dC}{dt}=G_C\!\left(C,C+\alpha_{CN}N\right)-d_CC-E_C(t,N,C)C\]"""
const ROUTE_B_EQUATIONS = raw"""\[\frac{dN}{dt}=G_N\!\left(N,N+\alpha_{NC}C\right)-d_NN-A_N(t)H_N(0.5)N\]
\[\frac{dS}{dt}=G_C\!\left(C,C+\alpha_{CN}N\right)\frac{S}{C}-d_CS-A_C(t)H_{CS}(0.5)S\]
\[\frac{dT}{dt}=G_C\!\left(C,C+\alpha_{CN}N\right)\frac{T}{C}-d_CT-A_C(t)H_{CT}(0.5)T,\qquad C=S+T\]"""
const SIMPLE_ROUTE_A_EQUATION = raw"""\[\frac{dN}{dt}=G_N\!\left(N,N+\alpha_{NC}C\right)-d_NN-k_NN,\qquad \frac{dC}{dt}=G_C\!\left(C,C+\alpha_{CN}N\right)-d_CC-k_CC\]"""
const BIC_EQUATION = raw"""\(n\log(\mathrm{SSE}_{scaled}/n)+k\log n\)"""

_escape(value) = replace(string(value), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
_label(model) = get(MODEL_LABELS, String(model), replace(String(model), "_" => " "))

function _selected_parameters(parameter_text::AbstractString)
    names_match = match(r"names = \[([^\]]+)\]", parameter_text)
    values_match = match(r"values = \[([^\]]+)\]", parameter_text)
    (names_match === nothing || values_match === nothing) && return DataFrame()
    names = strip.(replace.(split(names_match.captures[1], ","), ":" => ""))
    values = tryparse.(Float64, strip.(split(values_match.captures[1], ",")))
    all(value -> value !== nothing, values) || return DataFrame()
    return DataFrame(parameter = names, estimate = round.(Float64[value for value in values]; sigdigits = 5))
end

function _table_html(df::AbstractDataFrame)
    isempty(df) && return "<p>No rows available.</p>"
    headers = join("<th>$(_escape(name))</th>" for name in names(df))
    rows = String[]
    for row in eachrow(df)
        push!(rows, "<tr>" * join("<td>$(_escape(row[column]))</td>" for column in propertynames(row)) * "</tr>")
    end
    return "<div class=\"table-wrap\"><table><thead><tr>$(headers)</tr></thead><tbody>$(join(rows))</tbody></table></div>"
end

function _route_plot(overlay::DataFrame, model::String, pooling::String, path::String, title::String)
    selected = overlay[(String.(overlay.model) .== model) .& (String.(overlay.pooling_mode) .== pooling), :]
    isempty(selected) && error("No overlay rows for $(model), $(pooling)")
    environments = unique(select(selected, :density, :mix))
    sort!(environments, [:density, :mix])
    ymax = 1.05 * maximum(Float64.(vcat(selected.observed, selected.predicted)))
    panels = Any[]
    for environment in eachrow(environments)
        panel = plot(
            title = "$(environment.density), mix $(environment.mix)",
            xlabel = "Time (days)", ylabel = "Measured population",
            legend = false, ylim = (0, ymax), gridalpha = 0.18,
        )
        for (component, color) in (("sensitive", :crimson), ("resistant", :royalblue))
            rows = selected[
                (String.(selected.density) .== String(environment.density)) .&
                (String.(selected.mix) .== String(environment.mix)) .&
                (String.(selected.component) .== component), :,
            ]
            sort!(rows, :time)
            scatter!(panel, rows.time, rows.observed; color = color, markersize = 2.7, label = "")
            plot!(panel, rows.time, rows.predicted; color = color, linewidth = 2.0, label = "")
        end
        push!(panels, panel)
    end
    figure = plot(
        panels...; layout = (2, 3), size = (1180, 690), plot_title = title,
        margin = 4Plots.mm,
    )
    mkpath(dirname(path))
    savefig(figure, path)
    return path
end

function _route_b_plot(overlay::DataFrame, model::String, path::String, title::String)
    adjusted = copy(overlay)
    adjusted.pooling_mode = fill("shared", nrow(adjusted))
    return _route_plot(adjusted, model, "shared", path, title)
end

function _ranking_display(ranking::DataFrame)
    route_a = ranking[String.(ranking.route) .== ROUTE_A, :]
    route_b = ranking[String.(ranking.route) .== ROUTE_B, :]
    a_indices = unique(vcat(collect(1:min(3, nrow(route_a))), findfirst(String.(route_a.model) .== "dual_constant_kill")))
    selected = vcat(route_a[a_indices, :], first(route_b, min(3, nrow(route_b))); cols = :union)
    sort!(selected, :bic)
    return DataFrame(
        Rank = Int.(selected.rank),
        Route = replace.(String.(selected.route), " -> " => " - "),
        Model = [_label(model) for model in selected.model],
        Pooling = String.(selected.pooling_mode),
        Free_parameters = Int.(selected.n_parameters),
        BIC = round.(Float64.(selected.bic); digits = 2),
        Delta_BIC = round.(Float64.(selected.delta_bic); digits = 2),
        Boundary_flag = Bool.(selected.boundary_issue),
    )
end

function _mirror_file(source::String, package_root::String, mirror_directory::String)
    relative = relpath(source, package_root)
    destination = joinpath(mirror_directory, relative)
    mkpath(dirname(destination))
    cp(source, destination; force = true)
    return destination
end

function render_reduced_stage_report_html(
    ; start::AbstractString = pwd(), mirror_directory::Union{Nothing,AbstractString} = nothing,
      refit::Bool = false, maxiters::Int = 900,
)
    root = IOUtils.package_root(start)
    csv_dir = joinpath(root, "outputs", "csv", "reduced_stage_comparison")
    ranking_path = joinpath(csv_dir, "reduced_stage_model_ranking.csv")
    if refit || !isfile(ranking_path)
        FitWorkflows.run_reduced_stage_comparison!(; start = start, maxiters = maxiters)
    end
    ranking = CSV.read(ranking_path, DataFrame)
    route_a = ranking[String.(ranking.route) .== ROUTE_A, :]
    route_b = ranking[String.(ranking.route) .== ROUTE_B, :]
    winner_a = first(sort(route_a, :bic))
    winner_b = first(sort(route_b, :bic))
    simple_a = first(sort(route_a[(String.(route_a.model) .== "dual_constant_kill") .& (String.(route_a.pooling_mode) .== "shared"), :], :bic))

    image_dir = joinpath(root, "outputs", "images", "reduced_stage_comparison")
    route_a_overlay = CSV.read(joinpath(root, "outputs", "csv", "coculture_treated", "figures", "coculture_treated_joint_overlays.csv"), DataFrame)
    route_b_overlay = CSV.read(joinpath(csv_dir, "treatment_first_overlays.csv"), DataFrame)
    winner_a_image = _route_plot(route_a_overlay, String(winner_a.model), String(winner_a.pooling_mode), joinpath(image_dir, "competition_first_winner.png"), "Route A winner: $(_label(winner_a.model))")
    simple_a_image = _route_plot(route_a_overlay, String(simple_a.model), String(simple_a.pooling_mode), joinpath(image_dir, "competition_first_simple.png"), "Route A simplest: $(_label(simple_a.model))")
    winner_b_image = _route_b_plot(route_b_overlay, String(winner_b.model), joinpath(image_dir, "treatment_first_winner.png"), "Route B winner: $(_label(winner_b.model))")

    report_dir = joinpath(root, "outputs", "reports")
    mkpath(report_dir)
    relative_image(path) = replace(relpath(path, report_dir), "\\" => "/")
    display_ranking = _ranking_display(ranking)
    winner_a_params = _selected_parameters(String(winner_a.params))
    simple_a_params = _selected_parameters(String(simple_a.params))
    winner_b_params = _selected_parameters(String(winner_b.params))
    stability_path = joinpath(csv_dir, "treatment_first_parameter_stability.csv")
    stability = CSV.read(stability_path, DataFrame)
    winner_b_stability = stability[String.(stability.model) .== String(winner_b.model), :]
    stability_display = select(winner_b_stability, :parameter, :value_mean, :value_sd, :near_optimization_bound, :stability_class)

    html = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>A2780 reduced-stage comparison</title>
<script>window.MathJax={tex:{inlineMath:[["\\\\(","\\\\)"]],displayMath:[["\\\\[","\\\\]"]]}};</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
<style>
:root{--ink:#172033;--muted:#596579;--line:#cad3df;--paper:#fff;--wash:#f4f7fa;--accent:#315efb;--red:#c81e3a;--blue:#315de8}*{box-sizing:border-box}body{margin:0;background:var(--wash);color:var(--ink);font:16px/1.55 Arial,sans-serif}main{width:min(1160px,100%);margin:auto;background:var(--paper);padding:28px clamp(18px,4vw,48px) 64px}.back{color:var(--accent);font-weight:700;text-decoration:none}h1{font-size:clamp(28px,4vw,44px);line-height:1.1;margin:22px 0 10px;letter-spacing:0}h2{font-size:26px;margin:42px 0 12px;letter-spacing:0}h3{font-size:20px;margin:28px 0 8px;letter-spacing:0}.lead{color:var(--muted);font-size:19px;max-width:850px}.finding{border-left:5px solid var(--accent);padding:14px 18px;background:#eef3ff;margin:24px 0}.route{border-top:1px solid var(--line);padding-top:8px}.equation{overflow-x:auto;background:#f8fafc;border:1px solid var(--line);padding:12px 16px;margin:12px 0}.table-wrap{overflow:auto;border:1px solid var(--line)}table{border-collapse:collapse;width:100%;min-width:760px}th,td{padding:10px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}th{background:#edf1f5}figure{margin:20px 0}img{display:block;width:100%;height:auto;border:1px solid var(--line)}figcaption{color:var(--muted);margin-top:8px}.legend{display:flex;gap:22px;flex-wrap:wrap;margin:8px 0 20px}.dot{display:inline-block;width:12px;height:12px;border-radius:50%;margin-right:6px}.red{background:var(--red)}.blue{background:var(--blue)}details{border:1px solid var(--line);margin:14px 0}summary{cursor:pointer;font-weight:700;padding:12px 14px}details>div{padding:0 14px 14px}.warning{border-left:5px solid #b96b00;background:#fff7e6;padding:12px 16px}code{font-size:.92em}@media(max-width:700px){main{padding:18px 14px 44px}h2{font-size:23px}}
</style></head><body><main>
<a class="back" href="../../../../index.html">&#8592; Back</a>
<h1>Reduced-stage treated coculture comparison</h1>
<p class="lead">Can treated coculture be explained without fitting the treated-monoculture stage first? This report compares two conditional inheritance routes on exactly the same 168 treated-coculture observations.</p>
<div class="finding"><strong>Result.</strong> The competition-first route is favored: BIC $(round(Float64(winner_a.bic);digits=1)) versus $(round(Float64(winner_b.bic);digits=1)) for the best reciprocal route. The simplest two-parameter competition-first model has Delta BIC $(round(Float64(simple_a.delta_bic);digits=1)), so it is useful as a benchmark but not statistically interchangeable with the selected eight-parameter transit model.</div>

<h2>What was compared</h2>
<section class="route"><h3>Route A: Stage 1 &rarr; Stage 3 &rarr; Stage 4</h3><p>Stage 1 growth and the untreated Stage 3 asymmetric competition-with-loss parameters are fixed. Treatment parameters are estimated only from treated coculture. Stage 2 is not used.</p>
<div class="equation">$(ROUTE_A_EQUATIONS)</div></section>
<section class="route"><h3>Route B: Stage 1 &rarr; Stage 2 &rarr; Stage 4</h3><p>Stage 1 growth and Stage 2 timing/Hill treatment parameters are fixed. Stage 3 is not inherited; coculture competition or context terms are estimated directly from treated coculture.</p>
<div class="equation">$(ROUTE_B_EQUATIONS)</div></section>

<h2>Matched conditional BIC</h2>
<p>BIC is calculated as $(BIC_EQUATION), with <em>n</em> = 168 in every row. The scaled SSE divides each residual by that trajectory's observed peak before summing. Only parameters estimated in this conditional treated-coculture fit count toward <em>k</em>; inherited parameters are fixed. Lower BIC is better, and Delta BIC is measured from the best row.</p>
$(_table_html(display_ranking))

<h2>Competition-first fits</h2><div class="legend"><span><i class="dot red"></i>A2780Naive</span><span><i class="dot blue"></i>Total A2780cis</span></div>
<h3>Selected model</h3><p>The selected transit model allows treatment-driven movement into a damaged-but-visible compartment, clearance of that compartment, and load-dependent treatment scaling.</p>
<figure><img src="$(relative_image(winner_a_image))" alt="Competition-first selected model fitted to six treated coculture environments"><figcaption>Model: $(_escape(_label(winner_a.model))); pooling: $(_escape(winner_a.pooling_mode)). Every panel uses the same y-axis scale.</figcaption></figure>
<details><summary>Selected Route A parameters</summary><div>$(_table_html(winner_a_params))</div></details>
<h3>Simplest benchmark</h3><div class="equation">$(SIMPLE_ROUTE_A_EQUATION)</div>
<figure><img src="$(relative_image(simple_a_image))" alt="Simplest competition-first model fitted to six treated coculture environments"><figcaption>Model: constant direct loss with two shared fitted parameters. Delta BIC $(round(Float64(simple_a.delta_bic);digits=1)).</figcaption></figure>
<details><summary>Simplest Route A parameters</summary><div>$(_table_html(simple_a_params))</div></details>

<h2>Treatment-first reciprocal fit</h2>
<figure><img src="$(relative_image(winner_b_image))" alt="Treatment-first reciprocal model fitted to six treated coculture environments"><figcaption>Model: $(_escape(_label(winner_b.model))). Stage 2 response is fixed; four competition/loss parameters are estimated here.</figcaption></figure>
<details><summary>Selected Route B parameters</summary><div>$(_table_html(winner_b_params))</div></details>
<details><summary>Route B multistart stability</summary><div>$(_table_html(stability_display))</div></details>

<h2>Interpretation</h2>
<p>The treated-coculture trajectories retain information that is better captured when untreated coculture competition is established first and treatment is then added. Simply inheriting treated-monoculture response and estimating competition from treated coculture loses substantial support. This does <strong>not</strong> prove the eight-parameter transit mechanism is biologically correct: its winning fit is boundary-limited, including a cis loss estimate near zero and load-scaling terms at bounds.</p>
<div class="warning"><strong>Scientific decision.</strong> Keep the two-parameter constant-loss model as the transparent null model and the transit/load model as the statistical winner. Do not collapse to the null model without accepting a Delta BIC penalty of $(round(Float64(simple_a.delta_bic);digits=1)). Additional washout, viability, or damage-state measurements are needed to identify why the more flexible model fits better.</div>

<h2>Audit trail</h2><p>The reciprocal fits call <code>GrowthParameterEstimation.run_joint_multistart</code>; the competition-first artifacts were produced with <code>profile_joint_fit_bounds</code> and <code>run_joint_fit</code>. Both reject non-finite fits and failure sentinels, retain raw and scaled SSE, use fixed day-zero anchors, and write parameter/stability artifacts.</p>
<details><summary>Complete conditional ranking</summary><div>$(_table_html(select(ranking, :rank, :route, :model, :pooling_mode, :bic, :delta_bic, :n_parameters, :boundary_issue)))</div></details>
</main></body></html>"""

    report_path = joinpath(report_dir, REPORT_NAME)
    write(report_path, html)
    if mirror_directory !== nothing
        for path in (report_path, winner_a_image, simple_a_image, winner_b_image, ranking_path, stability_path, joinpath(csv_dir, "treatment_first_multistart.csv"), joinpath(csv_dir, "inheritance_route_audit.csv"))
            _mirror_file(path, root, String(mirror_directory))
        end
    end
    return report_path
end

end
