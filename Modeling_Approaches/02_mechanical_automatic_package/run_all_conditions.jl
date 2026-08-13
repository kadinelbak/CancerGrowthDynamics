using MechanicalAutomaticModeling
using CSV
using DataFrames
using Dates

conditions = ["monoculture_untreated", "monoculture_treated", "coculture_untreated", "coculture_treated"]
notebook_start = joinpath(@__DIR__, "notebooks")

for cond in conditions
    println("\n" * "="^60)
    println("CONDITION: $cond")
    println("="^60)

    try
        decoded = MechanicalAutomaticModeling.IOUtils.decode_condition_dataframe(cond; start = notebook_start)
        out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs(cond; start = notebook_start)
        decoded_path = joinpath(out.csv, "$(cond)_automatic_decoded.csv")
        CSV.write(decoded_path, decoded)
        MechanicalAutomaticModeling.IOUtils.write_manifest_row(condition = cond, step = "decode", outputs = [decoded_path], start = notebook_start)
        println("  decode_ok  : true  (rows=$(nrow(decoded)))")

        fit_artifacts = MechanicalAutomaticModeling.FitWorkflows.run_condition_fit!(decoded, cond; start = notebook_start)
        ranking = sort(fit_artifacts.ranking, :bic)
        best = first(ranking, 1)
        println("  fit_ok     : true  (fit_rows=$(nrow(ranking)))")
        println("  best_model : $(best.model[1])  BIC=$(round(best.bic[1], digits=2))")

        analysis_artifacts = MechanicalAutomaticModeling.AnalysisWorkflows.run_condition_analysis!(decoded, fit_artifacts, cond; start = notebook_start)
        println("  analysis_ok: true  (sensitivity_rows=$(nrow(analysis_artifacts.sensitivity)))")

        summary = DataFrame(
            condition = [cond],
            decoded_rows = [nrow(decoded)],
            fit_rows = [nrow(ranking)],
            sensitivity_rows = [nrow(analysis_artifacts.sensitivity)],
            generated_at = [Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS")],
        )
        summary_path = joinpath(out.metrics, "$(cond)_automatic_summary.csv")
        CSV.write(summary_path, summary)
        MechanicalAutomaticModeling.IOUtils.write_manifest_row(condition = cond, step = "summary", outputs = [summary_path], start = notebook_start)
    catch e
        println("  condition_ok: false  ($(e))")
    end

    println("  DONE")
end

println("\n" * "="^60)
println("All conditions complete.")
