using Dates
using CSV
using DataFrames

const PKG_ROOT = @__DIR__
const NOTEBOOK_START = joinpath(PKG_ROOT, "notebooks")

include(joinpath(PKG_ROOT, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

conditions = [
    "monoculture_untreated",
    "monoculture_treated",
    "coculture_untreated",
    "coculture_treated",
]

println("Running notebook workflows serially...")

for (i, condition) in enumerate(conditions)
    println("[$i/4] $condition")

    decoded = MechanicalAutomaticModeling.IOUtils.decode_condition_dataframe(condition; start = NOTEBOOK_START)
    out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs(condition; start = NOTEBOOK_START)
    decoded_path = joinpath(out.csv, "$(condition)_automatic_decoded.csv")
    CSV.write(decoded_path, decoded)
    MechanicalAutomaticModeling.IOUtils.write_manifest_row(condition = condition, step = "decode", outputs = [decoded_path], start = NOTEBOOK_START)

    fit_artifacts = MechanicalAutomaticModeling.FitWorkflows.run_condition_fit!(decoded, condition; start = NOTEBOOK_START)
    analysis_artifacts = MechanicalAutomaticModeling.AnalysisWorkflows.run_condition_analysis!(decoded, fit_artifacts, condition; start = NOTEBOOK_START)

    summary = DataFrame(
        condition = [condition],
        decoded_rows = [nrow(decoded)],
        fit_rows = [nrow(fit_artifacts.ranking)],
        sensitivity_rows = [nrow(analysis_artifacts.sensitivity)],
        generated_at = [Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS")],
    )
    summary_path = joinpath(out.metrics, "$(condition)_automatic_summary.csv")
    CSV.write(summary_path, summary)
    MechanicalAutomaticModeling.IOUtils.write_manifest_row(condition = condition, step = "summary", outputs = [summary_path], start = NOTEBOOK_START)

    println("    decoded_rows=$(nrow(decoded)) fit_rows=$(nrow(fit_artifacts.ranking)) sensitivity_rows=$(nrow(analysis_artifacts.sensitivity))")
end

println("Done.")
