using Pkg

const PACKAGE_ROOT = @__DIR__
Pkg.activate(PACKAGE_ROOT)

using CSV
using DataFrames
using Dates

include(joinpath(PACKAGE_ROOT, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

start = joinpath(PACKAGE_ROOT, "notebooks")
out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs("monoculture_treated"; start = start)
decoded_path = joinpath(out.csv, "monoculture_treated_a2780_decoded.csv")
decoded = CSV.read(decoded_path, DataFrame)

for (dose, expected_effect, expected_label) in ((0.67, 0.25, "IC25"), (1.0, 0.50, "IC50"), (1.47, 0.75, "IC75"))
    metadata = MechanicalAutomaticModeling.FitWorkflows._treated_dose_metadata(dose)
    metadata == (effect_level = expected_effect, ic_label = expected_label) || error("Incorrect metadata for $(dose) uM: $(metadata)")
end

max_time = parse(Float64, get(ENV, "A2780_MAX_TIME_PER_FIT", "1.0"))
fit = MechanicalAutomaticModeling.FitWorkflows.run_condition_fit!(
    decoded,
    "monoculture_treated";
    start = start,
    max_time_per_fit = max_time,
)

ranking = MechanicalAutomaticModeling.StagedA2780Workflow._append_treated_monoculture_literature_models!(
    fit.ranking,
    decoded;
    start = start,
    max_time_per_fit = max_time,
)
all(isfinite, ranking.bic) || error("Per-dose ranking contains non-finite BIC values")
all(isfinite, ranking.ssr) || error("Per-dose ranking contains non-finite SSR values")
all(ranking.ssr .< 1.0e12) || error("Per-dose ranking contains failure-sentinel SSR values")
CSV.write(fit.ranking_path, ranking)
CSV.write(fit.best_path, first(ranking, min(10, nrow(ranking))))

MechanicalAutomaticModeling.AnalysisWorkflows.run_condition_analysis!(
    decoded,
    merge(fit, (; ranking = ranking)),
    "monoculture_treated";
    start = start,
)

joint_path = joinpath(out.csv, "monoculture_treated_joint_dose_model_ranking.csv")
joint = CSV.read(joint_path, DataFrame; stringtype = String)
all(isfinite, joint.bic) || error("Joint ranking contains non-finite BIC values")
all(isfinite, joint.ssr) || error("Joint ranking contains non-finite SSR values")
all(joint.ssr .< 1.0e12) || error("Joint ranking contains failure-sentinel SSR values")
all(joint.ic_mapping .== "0.67=IC25;1.0=IC50;1.47=IC75") || error("Joint ranking contains stale IC mapping")

report_dir = joinpath(PACKAGE_ROOT, "outputs", "reports")
overview_path = joinpath(report_dir, "a2780_stage_overview.csv")
if isfile(overview_path)
    overview = CSV.read(overview_path, DataFrame; stringtype = String)
    treated_index = findfirst(==("monoculture_treated"), overview.condition)
    if treated_index !== nothing
        best = joint[argmin(joint.bic), :]
        overview[treated_index, :status] = "completed"
        overview[treated_index, :decoded_rows] = nrow(decoded)
        overview[treated_index, :fit_rows] = nrow(ranking) + nrow(joint)
        overview[treated_index, :best_model] = best.model
        overview[treated_index, :best_bic] = best.bic
        overview[treated_index, :best_ssr] = best.ssr
        overview[treated_index, :message] = "Refitted after correcting IC25/IC50/IC75 to 0.67/1.0/1.47 uM."
        CSV.write(overview_path, overview)
    end
end

manifest_path = joinpath(report_dir, "a2780_staged_manifest.csv")
if isfile(manifest_path)
    manifest = CSV.read(manifest_path, DataFrame; stringtype = String)
    push!(manifest, (
        timestamp_utc = now(UTC),
        condition = "monoculture_treated",
        status = "completed",
        message = "Full treated-monoculture refit with corrected dose mapping.",
        output_count = 4,
        outputs = join([fit.ranking_path, fit.best_path, joint_path, joinpath(out.csv, "monoculture_treated_trajectory_behavior.csv")], ";"),
    ))
    CSV.write(manifest_path, manifest)
end

html_path = MechanicalAutomaticModeling.StagedA2780Workflow.render_a2780_report_html(; start = start)
println("per_dose_rows=", nrow(ranking))
println("joint_rows=", nrow(joint))
println("joint_best=", joint[argmin(joint.bic), :model])
println("html_path=", html_path)
