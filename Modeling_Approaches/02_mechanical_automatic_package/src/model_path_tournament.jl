module ModelPathTournament

using CSV
using DataFrames
using Dates
using JSON3

using ..IOUtils
using ..HybridRefitService

export TournamentConfig, stage1_beam, run_model_path_tournament!

Base.@kwdef struct TournamentConfig
    candidates_per_lineage::Int = 3
    beam_width::Int = 5
    simple_delta_bic::Float64 = 10.0
    multistarts::Int = 20
    bootstrap_replicates::Int = 200
    blocked_validation::Bool = true
    max_time_per_fit::Float64 = 8.0
    force::Bool = false
end

const LINEAGES = ("A2780Naive", "A2780cis")
const SIMPLE_STAGE1 = ("logistic_growth", "theta_logistic_growth")
const STAGE2_REQUIRED = Dict(
    "A2780Naive" => "joint_ic_effect_hill_ramp_onset",
    "A2780cis" => "joint_ic_effect_two_population",
)
const STAGE3_COMPATIBLE = Set((
    "lv_symmetric_competition",
    "lv_asymmetric_competition",
    "lv_asymmetric_competition_death",
))

_truthy(value::Bool) = value
_truthy(value) = lowercase(string(value)) == "true"
_finite(value) = try isfinite(Float64(value)) catch; false end
_smoke_test() = lowercase(get(ENV, "A2780_TOURNAMENT_SMOKE", "false")) == "true"
_result_field(result, key::Symbol) = result isa AbstractDict ? result[String(key)] : getproperty(result, key)

function _eligible(df)
    mask = _finite.(df.bic)
    :eligible_for_inheritance in propertynames(df) && (mask .&= _truthy.(df.eligible_for_inheritance))
    return sort(df[mask, :], :bic)
end

function _unique_candidates(df, lineage, count, simple_delta)
    eligible = _eligible(df[String.(df.cell_line) .== lineage, :])
    isempty(eligible) && error("No eligible Stage 1 candidates for $(lineage)")
    chosen = DataFrame()
    seen = Set{Tuple{String,String}}()
    for row in eachrow(eligible)
        key = (String(row.model), String(row.pooling_mode))
        key in seen && continue
        push!(seen, key)
        push!(chosen, row; cols = :union)
        nrow(chosen) >= count && break
    end
    best_bic = Float64(first(eligible).bic)
    simple = eligible[
        in.(String.(eligible.model), Ref(SIMPLE_STAGE1)) .&
        (Float64.(eligible.bic) .- best_bic .<= simple_delta),
        :,
    ]
    if !isempty(simple)
        row = first(simple)
        key = (String(row.model), String(row.pooling_mode))
        key in seen || push!(chosen, row; cols = :union)
    end
    return chosen
end

"Select Stage 1 lineage pairs, retaining a simple near-optimal candidate when available."
function stage1_beam(; start = pwd(), config = TournamentConfig())
    root = IOUtils.package_root(start)
    path = joinpath(root, "outputs", "csv", "monoculture_untreated", "monoculture_untreated_pooling_model_ranking.csv")
    ranking = CSV.read(path, DataFrame)
    naive = _unique_candidates(ranking, "A2780Naive", config.candidates_per_lineage, config.simple_delta_bic)
    cis = _unique_candidates(ranking, "A2780cis", config.candidates_per_lineage, config.simple_delta_bic)
    rows = NamedTuple[]
    for n in eachrow(naive), c in eachrow(cis)
        push!(rows, (
            naive_model = String(n.model), naive_pooling = String(n.pooling_mode), naive_bic = Float64(n.bic),
            cis_model = String(c.model), cis_pooling = String(c.pooling_mode), cis_bic = Float64(c.bic),
            stage1_bic = Float64(n.bic) + Float64(c.bic),
            simple_path = String(n.model) in SIMPLE_STAGE1 && String(c.model) in SIMPLE_STAGE1,
        ))
    end
    paths = sort(DataFrame(rows), :stage1_bic)
    keep = collect(1:min(config.beam_width, nrow(paths)))
    best = Float64(first(paths).stage1_bic)
    simple_indices = findall(Bool.(paths.simple_path) .& (Float64.(paths.stage1_bic) .- best .<= config.simple_delta_bic))
    isempty(simple_indices) || push!(keep, first(simple_indices))
    return paths[sort(unique(keep)), :]
end

function _ranking_paths(workspace)
    csv = joinpath(workspace, "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "csv")
    return (
        stage2 = joinpath(csv, "monoculture_treated", "monoculture_treated_pooling_model_ranking.csv"),
        stage3 = joinpath(csv, "coculture_untreated", "coculture_untreated_pooling_model_ranking.csv"),
        stage4 = joinpath(csv, "coculture_treated", "linked_treatment_model_ranking.csv"),
    )
end

function _best_row(df; lineage = nothing, models = nothing)
    eligible = _eligible(df)
    lineage === nothing || (eligible = eligible[String.(eligible.cell_line) .== lineage, :])
    models === nothing || (eligible = eligible[in.(String.(eligible.model), Ref(models)), :])
    isempty(eligible) && error("No eligible candidate satisfies the tournament compatibility rules")
    return first(eligible)
end

_selection(row) = (model = String(row.model), pooling = String(row.pooling_mode))

function _canonical_seed(start)
    root = IOUtils.package_root(start)
    csv = joinpath(root, "outputs", "csv")
    stage2 = CSV.read(joinpath(csv, "monoculture_treated", "monoculture_treated_pooling_model_ranking.csv"), DataFrame)
    stage3 = CSV.read(joinpath(csv, "coculture_untreated", "coculture_untreated_pooling_model_ranking.csv"), DataFrame)
    stage4 = CSV.read(joinpath(csv, "coculture_treated", "linked_treatment_model_ranking.csv"), DataFrame)
    return (
        stage2_A2780Naive = _selection(_best_row(stage2; lineage = "A2780Naive", models = Set([STAGE2_REQUIRED["A2780Naive"]]))),
        stage2_A2780cis = _selection(_best_row(stage2; lineage = "A2780cis", models = Set([STAGE2_REQUIRED["A2780cis"]]))),
        stage3 = _selection(_best_row(stage3; models = STAGE3_COMPATIBLE)),
        stage4 = _selection(_best_row(stage4)),
    )
end

function _request(stage1, downstream, config)
    return (
        selections = (
            stage1_A2780Naive = (model = String(stage1.naive_model), pooling = String(stage1.naive_pooling)),
            stage1_A2780cis = (model = String(stage1.cis_model), pooling = String(stage1.cis_pooling)),
            downstream...,
        ),
        max_time_per_fit = config.max_time_per_fit,
    )
end

function _selected_downstream(workspace)
    paths = _ranking_paths(workspace)
    stage2 = CSV.read(paths.stage2, DataFrame)
    stage3 = CSV.read(paths.stage3, DataFrame)
    stage4 = CSV.read(paths.stage4, DataFrame)
    return (
        stage2_A2780Naive = _selection(_best_row(stage2; lineage = "A2780Naive", models = Set([STAGE2_REQUIRED["A2780Naive"]]))),
        stage2_A2780cis = _selection(_best_row(stage2; lineage = "A2780cis", models = Set([STAGE2_REQUIRED["A2780cis"]]))),
        stage3 = _selection(_best_row(stage3; models = STAGE3_COMPATIBLE)),
        stage4 = _selection(_best_row(stage4)),
    )
end

function _same_selections(a, b)
    return all(getproperty(a, name) == getproperty(b, name) for name in propertynames(a))
end

function _score_path(stage1, workspace)
    paths = _ranking_paths(workspace)
    stage2 = CSV.read(paths.stage2, DataFrame)
    stage3 = CSV.read(paths.stage3, DataFrame)
    stage4 = CSV.read(paths.stage4, DataFrame)
    n2 = _best_row(stage2; lineage = "A2780Naive", models = Set([STAGE2_REQUIRED["A2780Naive"]]))
    c2 = _best_row(stage2; lineage = "A2780cis", models = Set([STAGE2_REQUIRED["A2780cis"]]))
    s3 = _best_row(stage3; models = STAGE3_COMPATIBLE)
    s4 = _best_row(stage4)
    stage2_bic = Float64(n2.bic) + Float64(c2.bic)
    cumulative = Float64(stage1.stage1_bic) + stage2_bic + Float64(s3.bic) + Float64(s4.bic)
    boundary_count = sum(hasproperty(row, :boundary_issue) && _truthy(row.boundary_issue) for row in (n2, c2, s3, s4))
    return (
        naive_model = String(stage1.naive_model), naive_pooling = String(stage1.naive_pooling),
        cis_model = String(stage1.cis_model), cis_pooling = String(stage1.cis_pooling),
        stage1_bic = Float64(stage1.stage1_bic),
        naive_stage2_model = String(n2.model), naive_stage2_pooling = String(n2.pooling_mode),
        cis_stage2_model = String(c2.model), cis_stage2_pooling = String(c2.pooling_mode),
        stage2_bic = stage2_bic,
        stage3_model = String(s3.model), stage3_pooling = String(s3.pooling_mode), stage3_bic = Float64(s3.bic),
        stage4_model = String(s4.model), stage4_pooling = String(s4.pooling_mode), stage4_bic = Float64(s4.bic),
        cumulative_path_score = cumulative,
        boundary_count = boundary_count,
        workspace = workspace,
    )
end

_escape(value) = replace(string(value), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
_fmt(value; digits = 1) = isfinite(Float64(value)) ? string(round(Float64(value); digits = digits)) : "NA"

function _validation_html(output)
    blocked_path = joinpath(output, "linked_treatment_blocked_validation.csv")
    bootstrap_path = joinpath(output, "linked_treatment_whole_fit_bootstrap_parameters.csv")
    blocked = if isfile(blocked_path)
        data = CSV.read(blocked_path, DataFrame)
        body = join([
            "<tr><td>$(_escape(row.holdout_scheme))</td><td>$(_escape(row.block))</td><td>$(_escape(row.status))</td><td>$(_fmt(row.nrmse))</td></tr>"
            for row in eachrow(data)
        ], "")
        "<h3>Blocked validation</h3><div class=\"table\"><table><thead><tr><th>Scheme</th><th>Held-out block</th><th>Status</th><th>Held-out nRMSE</th></tr></thead><tbody>$(body)</tbody></table></div>"
    else
        "<p>Blocked-validation artifacts were not generated for this run.</p>"
    end
    bootstrap = if isfile(bootstrap_path)
        data = CSV.read(bootstrap_path, DataFrame)
        body = join([
            "<tr><td>$(_escape(row.parameter))</td><td>$(_fmt(row.median; digits = 3))</td><td>$(_fmt(row.ci95_lower; digits = 3))</td><td>$(_fmt(row.ci95_upper; digits = 3))</td><td>$(_escape(row.successful_replicates))</td></tr>"
            for row in eachrow(data)
        ], "")
        "<h3>Whole-fit bootstrap</h3><div class=\"table\"><table><thead><tr><th>Parameter</th><th>Median</th><th>2.5%</th><th>97.5%</th><th>Successful refits</th></tr></thead><tbody>$(body)</tbody></table></div>"
    else
        "<p>Whole-fit bootstrap artifacts were not generated for this run.</p>"
    end
    return blocked * bootstrap
end

function _report_html(results, config, generated, output)
    rows = join([
        "<tr><td>$(index)</td><td>$(_escape(row.naive_model))<br><small>$(_escape(row.naive_pooling))</small></td>" *
        "<td>$(_escape(row.cis_model))<br><small>$(_escape(row.cis_pooling))</small></td>" *
        "<td>$(_escape(row.stage3_model))</td><td>$(_escape(row.stage4_model))</td>" *
        "<td>$(_fmt(row.cumulative_path_score))</td><td>$(row.boundary_count)</td></tr>"
        for (index, row) in enumerate(eachrow(results))
    ], "\n")
    return """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>A2780 conditional model-path tournament</title>
<style>:root{--ink:#172033;--muted:#5b6472;--line:#cbd3df;--soft:#f3f6fa;--accent:#3156d3}*{box-sizing:border-box}body{margin:0;font:16px/1.5 "Segoe UI",Arial,sans-serif;color:var(--ink)}main{width:min(1180px,calc(100% - 32px));margin:auto;padding:24px 0 64px}a{color:var(--accent)}h1{font-size:30px;margin:22px 0 8px;letter-spacing:0}h2{font-size:21px;margin:30px 0 8px;letter-spacing:0}p{max-width:920px;color:var(--muted)}.summary{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));border:1px solid var(--line);margin:22px 0}.summary div{padding:14px;border-right:1px solid var(--line)}.summary div:last-child{border:0}.summary strong{display:block;font-size:22px}.table{overflow:auto;border:1px solid var(--line)}table{border-collapse:collapse;width:100%;min-width:900px}th,td{text-align:left;padding:10px 12px;border-bottom:1px solid var(--line);vertical-align:top}th{background:var(--soft);font-size:13px}tr:first-child td{font-weight:650;background:#f7f9ff}small{color:var(--muted)}details{border-top:1px solid var(--line);padding:14px 0}code{font-size:13px}@media(max-width:700px){.summary{grid-template-columns:1fr 1fr}.summary div:nth-child(2){border-right:0}}</style>
</head><body><main><a href="../../../../index.html">&larr; Back</a><h1>Conditional Model-Path Tournament</h1>
<p>Each retained Stage 1 lineage pair was carried into fresh Stage 2 and Stage 3 fits. The compatible winners were then inherited by a fresh linked Stage 4 fit. This is a conditional refit, not a visual substitution of previously exported curves.</p>
<div class="summary"><div><small>Complete paths</small><strong>$(nrow(results))</strong></div><div><small>Starts per fit</small><strong>$(config.multistarts)</strong></div><div><small>Winner bootstrap refits</small><strong>$(config.bootstrap_replicates)</strong></div><div><small>Generated</small><strong>$(generated)</strong></div></div>
<h2>Retained Complete Paths</h2><div class="table"><table><thead><tr><th>Rank</th><th>Stage 1 A2780Naive</th><th>Stage 1 A2780cis</th><th>Stage 3</th><th>Stage 4</th><th>Cumulative path score</th><th>Boundary flags</th></tr></thead><tbody>$(rows)</tbody></table></div>
<h2>How To Read The Score</h2><p>Every stage retains its own BIC. The cumulative path score is the sum used to prune the beam and compare complete inheritance paths. It is not a new global BIC because fitted datasets and inherited information overlap across stages. Smaller values are favored only as a tournament heuristic.</p>
<h2>Fitting Tournament</h2><p>Stages 1-3, linked Stage 4 finalists, and the winning complete chain use $(config.multistarts) deterministic dispersed starts. The broad legacy treated-coculture diagnostic sweep uses three-start successive-halving screening before compatible finalists receive the full budget. Delayed and stiff candidates use bounded Nelder-Mead; smooth candidates can be refined with bounded BFGS. GrowthParameterEstimation records failed starts, boundary profiles, fitted parameters, and prediction artifacts in each path workspace.</p>
<details><summary><strong>Compatibility and screening limits</strong></summary><p>The current population-balance Stage 4 implementation requires delayed Hill-ramp A2780Naive treatment, a sensitive/tolerant A2780cis treatment model, and an LV competition equation. Other Stage 2 and Stage 3 candidates remain in the stage ranking tables but are not silently forced into an incompatible Stage 4 state system. The two seven-parameter interaction-scaled transit diagnostics remain available in the ordinary staged analysis, but are excluded from repeated beam-path refits after exceeding the per-candidate screening budget; the simpler additive transit model remains in the tournament.</p></details>
<h2>Validation</h2><p>Boundary counts are shown above. The winning complete path receives a within-trajectory residual bootstrap with complete model refitting. Density, dose, mixture, and context blocks are also withheld and refitted in turn. These are distinct from endpoint well-resampling.</p>$(_validation_html(output))
</main></body></html>"""
end

function _write_report(results, config, root, output)
    generated = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM")
    report = joinpath(root, "outputs", "reports", "a2780_model_path_tournament.html")
    mkpath(dirname(report))
    write(report, _report_html(results, config, generated, output))
    repo = IOUtils.find_repo_root(root)
    docs = joinpath(repo, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "reports", basename(report))
    mkpath(dirname(docs))
    cp(report, docs; force = true)
    return (report = report, docs_report = docs)
end

"Run the conditional beam tournament and write its CSV, JSON, and HTML artifacts."
function run_model_path_tournament!(; start = pwd(), config = TournamentConfig())
    config.multistarts >= 1 || error("multistarts must be positive")
    root = IOUtils.package_root(start)
    output = joinpath(root, "outputs", "model_path_tournament")
    cache = joinpath(tempdir(), "a2780_model_path_cache")
    mkpath(output)
    beam = stage1_beam(; start = start, config = config)
    CSV.write(joinpath(output, "stage1_beam.csv"), beam)
    seed = _canonical_seed(start)
    scored = NamedTuple[]
    old_multistarts = get(ENV, "A2780_MULTISTARTS", nothing)
    old_refine = get(ENV, "A2780_REFINE_BFGS", nothing)
    old_bootstrap = get(ENV, "A2780_WHOLE_FIT_BOOTSTRAPS", nothing)
    old_blocked = get(ENV, "A2780_BLOCKED_VALIDATION", nothing)
    old_maxtime = get(ENV, "A2780_FIT_MAXTIME", nothing)
    old_conditional = get(ENV, "A2780_CONDITIONAL_TOURNAMENT", nothing)
    ENV["A2780_MULTISTARTS"] = string(config.multistarts)
    ENV["A2780_REFINE_BFGS"] = _smoke_test() ? "false" : "true"
    ENV["A2780_WHOLE_FIT_BOOTSTRAPS"] = "0"
    ENV["A2780_BLOCKED_VALIDATION"] = "false"
    ENV["A2780_FIT_MAXTIME"] = string(config.max_time_per_fit)
    ENV["A2780_CONDITIONAL_TOURNAMENT"] = "true"
    try
        for (index, stage1) in enumerate(eachrow(beam))
            println("Tournament path $(index)/$(nrow(beam)): $(stage1.naive_model) + $(stage1.cis_model)")
            first_result = HybridRefitService.run_conditional_refit!(
                _request(stage1, seed, config); start = start, cache_root = cache, force = config.force,
            )
            workspace = String(_result_field(first_result, :workspace))
            selected = _selected_downstream(workspace)
            if !_same_selections(seed, selected)
                final_result = HybridRefitService.run_conditional_refit!(
                    _request(stage1, selected, config); start = start, cache_root = cache, force = config.force,
                )
                workspace = String(_result_field(final_result, :workspace))
            end
            push!(scored, _score_path(stage1, workspace))
        end
    finally
        old_multistarts === nothing ? delete!(ENV, "A2780_MULTISTARTS") : (ENV["A2780_MULTISTARTS"] = old_multistarts)
        old_refine === nothing ? delete!(ENV, "A2780_REFINE_BFGS") : (ENV["A2780_REFINE_BFGS"] = old_refine)
        old_bootstrap === nothing ? delete!(ENV, "A2780_WHOLE_FIT_BOOTSTRAPS") : (ENV["A2780_WHOLE_FIT_BOOTSTRAPS"] = old_bootstrap)
        old_blocked === nothing ? delete!(ENV, "A2780_BLOCKED_VALIDATION") : (ENV["A2780_BLOCKED_VALIDATION"] = old_blocked)
        old_maxtime === nothing ? delete!(ENV, "A2780_FIT_MAXTIME") : (ENV["A2780_FIT_MAXTIME"] = old_maxtime)
        old_conditional === nothing ? delete!(ENV, "A2780_CONDITIONAL_TOURNAMENT") : (ENV["A2780_CONDITIONAL_TOURNAMENT"] = old_conditional)
    end
    results = sort(DataFrame(scored), :cumulative_path_score)
    results.rank = collect(1:nrow(results))
    results.delta_path_score = results.cumulative_path_score .- minimum(results.cumulative_path_score)
    CSV.write(joinpath(output, "complete_path_ranking.csv"), results)
    if config.bootstrap_replicates > 0 && !isempty(results)
        winner = first(results)
        downstream = (
            stage2_A2780Naive = (model = String(winner.naive_stage2_model), pooling = String(winner.naive_stage2_pooling)),
            stage2_A2780cis = (model = String(winner.cis_stage2_model), pooling = String(winner.cis_stage2_pooling)),
            stage3 = (model = String(winner.stage3_model), pooling = String(winner.stage3_pooling)),
            stage4 = (model = String(winner.stage4_model), pooling = String(winner.stage4_pooling)),
        )
        stage1 = (
            naive_model = String(winner.naive_model), naive_pooling = String(winner.naive_pooling),
            cis_model = String(winner.cis_model), cis_pooling = String(winner.cis_pooling),
        )
        old_bootstrap = get(ENV, "A2780_WHOLE_FIT_BOOTSTRAPS", nothing)
        old_blocked = get(ENV, "A2780_BLOCKED_VALIDATION", nothing)
        old_multistarts = get(ENV, "A2780_MULTISTARTS", nothing)
        old_refine = get(ENV, "A2780_REFINE_BFGS", nothing)
        old_maxtime = get(ENV, "A2780_FIT_MAXTIME", nothing)
        old_conditional = get(ENV, "A2780_CONDITIONAL_TOURNAMENT", nothing)
        ENV["A2780_WHOLE_FIT_BOOTSTRAPS"] = string(config.bootstrap_replicates)
        ENV["A2780_BLOCKED_VALIDATION"] = string(config.blocked_validation)
        ENV["A2780_MULTISTARTS"] = string(config.multistarts)
        ENV["A2780_REFINE_BFGS"] = _smoke_test() ? "false" : "true"
        ENV["A2780_FIT_MAXTIME"] = string(config.max_time_per_fit)
        ENV["A2780_CONDITIONAL_TOURNAMENT"] = "true"
        try
            refitted = HybridRefitService.run_conditional_refit!(
                _request(stage1, downstream, config); start = start, cache_root = cache, force = true,
            )
            source = joinpath(String(_result_field(refitted, :workspace)), "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "csv", "coculture_treated")
            for name in ("linked_treatment_whole_fit_bootstrap_parameters.csv", "linked_treatment_whole_fit_bootstrap_predictions.csv", "linked_treatment_whole_fit_bootstrap_failures.csv", "linked_treatment_blocked_validation.csv")
                isfile(joinpath(source, name)) && cp(joinpath(source, name), joinpath(output, name); force = true)
            end
        finally
            old_bootstrap === nothing ? delete!(ENV, "A2780_WHOLE_FIT_BOOTSTRAPS") : (ENV["A2780_WHOLE_FIT_BOOTSTRAPS"] = old_bootstrap)
            old_blocked === nothing ? delete!(ENV, "A2780_BLOCKED_VALIDATION") : (ENV["A2780_BLOCKED_VALIDATION"] = old_blocked)
            old_multistarts === nothing ? delete!(ENV, "A2780_MULTISTARTS") : (ENV["A2780_MULTISTARTS"] = old_multistarts)
            old_refine === nothing ? delete!(ENV, "A2780_REFINE_BFGS") : (ENV["A2780_REFINE_BFGS"] = old_refine)
            old_maxtime === nothing ? delete!(ENV, "A2780_FIT_MAXTIME") : (ENV["A2780_FIT_MAXTIME"] = old_maxtime)
            old_conditional === nothing ? delete!(ENV, "A2780_CONDITIONAL_TOURNAMENT") : (ENV["A2780_CONDITIONAL_TOURNAMENT"] = old_conditional)
        end
    end
    open(joinpath(output, "tournament_manifest.json"), "w") do io
        JSON3.write(io, (generated_at = string(now()), config = config, paths = nrow(results), score_warning = "additive beam heuristic; not a global BIC"))
    end
    reports = _write_report(results, config, root, output)
    return (ranking = results, output = output, reports...)
end

end
