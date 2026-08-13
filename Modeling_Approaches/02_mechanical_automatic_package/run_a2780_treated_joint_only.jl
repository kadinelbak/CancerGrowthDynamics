using Pkg
Pkg.activate(@__DIR__)

using Dates
using CSV
using DataFrames

include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

start = joinpath(@__DIR__, "notebooks")
out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs("monoculture_treated"; start = start)
decoded_path = joinpath(out.csv, "monoculture_treated_a2780_decoded.csv")
isfile(decoded_path) || error("Missing cached A2780 treated monoculture data: $decoded_path")

decoded = CSV.read(decoded_path, DataFrame)
max_time = parse(Float64, get(ENV, "A2780_MAX_TIME_PER_FIT", "1.0"))
joint = MechanicalAutomaticModeling.FitWorkflows._run_joint_treated_monoculture_fitting(
    decoded,
    out;
    start = start,
    max_time_per_fit = max_time,
)

isempty(joint) && error("No finite joint treated-monoculture fits were produced")
all(isfinite, joint.bic) || error("Joint ranking contains non-finite BIC values")
all(isfinite, joint.ssr) || error("Joint ranking contains non-finite SSR values")
all(joint.ssr .< 1.0e12) || error("Joint ranking contains failure-sentinel SSR values")

report_dir = joinpath(@__DIR__, "outputs", "reports")
overview_path = joinpath(report_dir, "a2780_stage_overview.csv")
if isfile(overview_path)
    overview = CSV.read(overview_path, DataFrame; stringtype = String)
    treated_index = findfirst(==("monoculture_treated"), overview.condition)
    if treated_index !== nothing
        best = joint[argmin(joint.bic), :]
        overview[treated_index, :status] = "completed"
        overview[treated_index, :decoded_rows] = nrow(decoded)
        overview[treated_index, :fit_rows] = nrow(joint)
        overview[treated_index, :best_model] = best.model
        overview[treated_index, :best_bic] = best.bic
        overview[treated_index, :best_ssr] = best.ssr
        overview[treated_index, :message] = "Dose-linked treated-monoculture validation completed; coculture stages were not resumed."
        CSV.write(overview_path, overview)
    end
end

manifest_path = joinpath(report_dir, "a2780_staged_manifest.csv")
if isfile(manifest_path)
    manifest = CSV.read(manifest_path, DataFrame; stringtype = String)
    ranking_path = joinpath(out.csv, "monoculture_treated_joint_dose_model_ranking.csv")
    behavior_path = joinpath(out.csv, "monoculture_treated_trajectory_behavior.csv")
    outputs = join([ranking_path, behavior_path], ";")
    push!(manifest, (
        timestamp_utc = now(UTC),
        condition = "monoculture_treated",
        status = "completed",
        message = "Dose-linked joint validation completed with finite metrics.",
        output_count = 2,
        outputs = outputs,
    ))
    CSV.write(manifest_path, manifest)
end
html_path = MechanicalAutomaticModeling.StagedA2780Workflow.render_a2780_report_html(; start = start)

println("joint_rows=", nrow(joint))
println("joint_models=", join(sort(unique(joint.model)), ","))
println("joint_path=", joinpath(out.csv, "monoculture_treated_joint_dose_model_ranking.csv"))
println("html_path=", html_path)



