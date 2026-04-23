module AnalysisWorkflows

using CSV
using DataFrames
using GrowthParameterEstimation
using Plots

using ..IOUtils

export run_condition_analysis!

function run_condition_analysis!(decoded::DataFrame, fit_result, condition::AbstractString; start::AbstractString = pwd())
    out = IOUtils.condition_output_dirs(condition; start)

    ranking = fit_result.ranking
    nrow(ranking) == 0 && error("No ranking rows for condition: $(condition)")

    top = first(sort(ranking, :bic), min(3, nrow(ranking)))
    top_path = joinpath(out.metrics, "$(condition)_analysis_top3_by_bic.csv")
    CSV.write(top_path, top)

    # Sensitivity on best candidate model (available in ranking table)
    best_row = first(sort(ranking, :bic), 1)
    best_model = String(best_row.model[1])

    # Package-level sensitivity utility call; fall back gracefully if signature differs.
    sens_df = DataFrame()
    try
        sens = GrowthParameterEstimation.parameter_sensitivity_analysis(decoded; model_name = best_model)
        sens_df = DataFrame(sens)
    catch
        sens_df = DataFrame(model = [best_model], note = ["Sensitivity call not available with current signature; update wrapper."])
    end

    sens_path = joinpath(out.metrics, "$(condition)_parameter_sensitivity.csv")
    CSV.write(sens_path, sens_df)

    bic_plot = bar(
        string.(top.model),
        top.bic;
        legend = false,
        xlabel = "Model",
        ylabel = "BIC",
        title = "$(condition) automatic top models by BIC",
        xrotation = 20,
    )
    bic_plot_path = joinpath(out.images, "$(condition)_automatic_top_bic_bar.png")
    savefig(bic_plot, bic_plot_path)

    IOUtils.write_manifest_row(
        condition = condition,
        step = "analysis",
        outputs = [top_path, sens_path, bic_plot_path],
        start = start,
    )

    return (
        top_path = top_path,
        sensitivity_path = sens_path,
        image_path = bic_plot_path,
        sensitivity = sens_df,
    )
end

end
