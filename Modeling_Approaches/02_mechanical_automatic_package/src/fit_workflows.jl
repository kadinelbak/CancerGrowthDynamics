module FitWorkflows

using CSV
using DataFrames
using GrowthParameterEstimation

using ..IOUtils
using ..ModelRegistry

export default_model_sets, run_condition_fit!

function default_model_sets()
    return Dict(
        "monoculture_untreated" => ["logistic_growth", "gompertz_growth"],
        "monoculture_treated" => ["theta_logistic_hill_kill", "pkpd_inhibition", "transit_chain_erlang", "adaptive_ic50", "sensitive_resistant"],
        "coculture_untreated" => ["lotka_volterra_competition", "lotka_volterra_hill_competition", "null_coculture"],
        "coculture_treated" => ["lotka_volterra_hill_competition", "theta_logistic_hill_kill", "pkpd_inhibition"],
    )
end

function run_condition_fit!(decoded::DataFrame, condition::AbstractString; start::AbstractString = pwd())
    ModelRegistry.ensure_model_registry!()
    out = IOUtils.condition_output_dirs(condition; start)

    model_map = default_model_sets()
    include_models = get(model_map, condition, ["logistic_growth", "gompertz_growth"])

    cfg = GrowthParameterEstimation.default_config(
        output_dir = out.csv,
    )

    result = GrowthParameterEstimation.run_pipeline(
        decoded;
        config = cfg,
        include_models = include_models,
    )

    ranking_df = DataFrame(result.ranking)
    ranking_path = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
    CSV.write(ranking_path, ranking_df)

    best_df = first(sort(ranking_df, :bic), min(10, nrow(ranking_df)))
    best_path = joinpath(out.csv, "$(condition)_automatic_best_models_top10.csv")
    CSV.write(best_path, best_df)

    IOUtils.write_manifest_row(
        condition = condition,
        step = "fit",
        outputs = [ranking_path, best_path],
        start = start,
    )

    return (result = result, ranking = ranking_df, ranking_path = ranking_path, best_path = best_path)
end

end
