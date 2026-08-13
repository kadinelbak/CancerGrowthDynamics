using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames
using Dates

include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

start = joinpath(@__DIR__, "notebooks")
max_time = parse(Float64, get(ENV, "A2780_MAX_TIME_PER_FIT", "1.0"))
conditions = ("monoculture_untreated", "monoculture_treated")
results = Dict{String,Any}()

for condition in conditions
    decoded = MechanicalAutomaticModeling.StagedA2780Workflow.decode_a2780_condition(condition; start = start)
    out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs(condition; start = start)
    CSV.write(joinpath(out.csv, "$(condition)_a2780_decoded.csv"), decoded)
    reuse_untreated = condition == "monoculture_untreated" &&
        lowercase(get(ENV, "A2780_REUSE_UNTREATED_POOLING", "true")) in ("1", "true", "yes")
    untreated_ranking_path = joinpath(out.csv, "monoculture_untreated_pooling_model_ranking.csv")
    untreated_status_path = joinpath(out.csv, "monoculture_untreated_pooling_status.csv")
    untreated_overlay_path = joinpath(out.csv, "figures", "monoculture_untreated_pooling_overlays.csv")
    reuse_treated = condition == "monoculture_treated" &&
        lowercase(get(ENV, "A2780_REUSE_TREATED_POOLING", "true")) in ("1", "true", "yes")
    treated_ranking_path = joinpath(out.csv, "monoculture_treated_pooling_model_ranking.csv")
    treated_status_path = joinpath(out.csv, "monoculture_treated_pooling_status.csv")
    treated_overlay_path = joinpath(out.csv, "figures", "monoculture_treated_joint_dose_overlays.csv")
    if reuse_untreated && all(isfile, (untreated_ranking_path, untreated_status_path, untreated_overlay_path))
        ranking = CSV.read(untreated_ranking_path, DataFrame)
        status = CSV.read(untreated_status_path, DataFrame)
        top5 = CSV.read(joinpath(out.csv, "monoculture_untreated_pooling_top5.csv"), DataFrame)
        overlay = CSV.read(untreated_overlay_path, DataFrame)
        MechanicalAutomaticModeling.FitWorkflows._render_untreated_pooling_graph_grid(overlay, top5, out)
        fit = (
            ranking = ranking,
            ranking_path = untreated_ranking_path,
            best_path = joinpath(out.csv, "monoculture_untreated_pooling_top5.csv"),
            inheritance_allowed = all(Bool.(status.inheritance_allowed)),
        )
    elseif reuse_treated && all(isfile, (treated_ranking_path, treated_status_path, treated_overlay_path))
        ranking = CSV.read(treated_ranking_path, DataFrame)
        status = CSV.read(treated_status_path, DataFrame)
        overlay = CSV.read(treated_overlay_path, DataFrame)
        MechanicalAutomaticModeling.FitWorkflows._render_treated_pooling_graph_grid(overlay, status, out)
        fit = (
            ranking = ranking,
            ranking_path = treated_ranking_path,
            best_path = joinpath(out.csv, "monoculture_treated_pooling_top5.csv"),
        )
    else
        fit = MechanicalAutomaticModeling.FitWorkflows.run_condition_fit!(
            decoded,
            condition;
            start = start,
            max_time_per_fit = max_time,
        )
    end
    results[condition] = (decoded = decoded, fit = fit)
end

untreated = results["monoculture_untreated"].fit
hasproperty(untreated, :inheritance_allowed) && untreated.inheritance_allowed ||
    error("Untreated independent diagnostics indicate inadequate pooling")

for condition in conditions
    ranking = results[condition].fit.ranking
    all(isfinite, Float64.(ranking.bic)) || error("$(condition) has non-finite BIC")
    all(isfinite, Float64.(ranking.ssr)) || error("$(condition) has non-finite raw SSR")
    all(Float64.(ranking.ssr) .< 9.99e11) || error("$(condition) contains failure-sentinel SSR")
end

treated_out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs("monoculture_treated"; start = start)
treated_status = CSV.read(joinpath(treated_out.csv, "monoculture_treated_pooling_status.csv"), DataFrame)
all(Bool.(treated_status.inheritance_allowed)) || error("Treated independent diagnostics indicate inadequate pooling")
best_overlays = CSV.read(
    joinpath(treated_out.csv, "figures", "monoculture_treated_best_joint_model_by_environment_overlays.csv"),
    DataFrame,
)
panel_count = nrow(unique(best_overlays[:, [:cell_line, :density, :dose]]))
panel_count == 12 || error("Expected 12 treated winner panels, found $(panel_count)")

report_dir = joinpath(@__DIR__, "outputs", "reports")
mkpath(report_dir)
overview_path = joinpath(report_dir, "a2780_stage_overview.csv")
overview = isfile(overview_path) ? CSV.read(overview_path, DataFrame; stringtype = String) : DataFrame()
for condition in conditions
    result = results[condition]
    ranking = result.fit.ranking
    best = ranking[argmin(ranking.bic), :]
    index = isempty(overview) || !(:condition in propertynames(overview)) ? nothing : findfirst(==(condition), String.(overview.condition))
    index === nothing && continue
    overview[index, :status] = "completed"
    overview[index, :decoded_rows] = nrow(result.decoded)
    overview[index, :fit_rows] = nrow(ranking)
    overview[index, :best_model] = String(best.model)
    overview[index, :best_bic] = Float64(best.bic)
    overview[index, :best_ssr] = Float64(best.ssr)
    overview[index, :message] = "Density-aware shared/partial pooling refit completed."
end
!isempty(overview) && CSV.write(overview_path, overview)

html_path = MechanicalAutomaticModeling.StagedA2780Workflow.render_a2780_report_html(; start = start)
println("untreated_rows=", nrow(results["monoculture_untreated"].fit.ranking))
println("treated_rows=", nrow(results["monoculture_treated"].fit.ranking))
println("treated_panels=", panel_count)
println("html_path=", html_path)
