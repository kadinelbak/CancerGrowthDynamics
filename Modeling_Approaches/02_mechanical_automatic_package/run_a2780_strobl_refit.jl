using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using CSV
using DataFrames

include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

const FW = MechanicalAutomaticModeling.FitWorkflows
const SW = MechanicalAutomaticModeling.StagedA2780Workflow
const IO = MechanicalAutomaticModeling.IOUtils

start = joinpath(@__DIR__, "notebooks")
max_time = parse(Float64, get(ENV, "A2780_MAX_TIME_PER_FIT", "2.0"))

function decode_and_write(condition)
    decoded = SW.decode_a2780_condition(condition; start)
    out = IO.condition_output_dirs(condition; start)
    CSV.write(joinpath(out.csv, "$(condition)_a2780_decoded.csv"), decoded)
    return decoded, out
end

untreated, untreated_out = decode_and_write("coculture_untreated")
untreated_fit = FW.run_condition_fit!(
    untreated,
    "coculture_untreated";
    start,
    max_time_per_fit = max_time,
)

decode_and_write("monoculture_treated")
treated, treated_out = decode_and_write("coculture_treated")
recovered, provenance = FW._recover_coculture_design(treated, "coculture_treated")
environments, initial = FW._coculture_environments(recovered)
CSV.write(joinpath(treated_out.csv, "coculture_treated_density_provenance_validation.csv"), provenance)
CSV.write(joinpath(treated_out.csv, "coculture_treated_initial_mix_diagnostics.csv"), initial)

linked_ranking = FW._fit_linked_treatment_joint(
    environments,
    treated_out;
    start,
    max_time_per_fit = max_time,
    resume_from_existing = false,
)
CSV.write(joinpath(treated_out.csv, "coculture_treated_automatic_model_ranking.csv"), linked_ranking)
CSV.write(
    joinpath(treated_out.csv, "coculture_treated_automatic_best_models_top10.csv"),
    first(linked_ranking, min(FW.LINKED_TREATMENT_VISIBLE_LIMIT, nrow(linked_ranking))),
)

for (stage, ranking) in (("untreated", untreated_fit.ranking), ("linked treatment", linked_ranking))
    all(isfinite, Float64.(ranking.bic)) || error("$(stage) ranking contains non-finite BIC")
    all(isfinite, Float64.(ranking.ssr)) || error("$(stage) ranking contains non-finite raw SSE")
    all(Float64.(ranking.ssr) .< 9.99e11) || error("$(stage) ranking contains a failure sentinel")
    strobl_count = count(identity, startswith.(String.(ranking.model), "strobl_"))
    strobl_count == length(FW.STROBL_MODEL_VARIANTS) || error("$(stage) ranking contains $(strobl_count) Strobl rows")
end

summary = SW.refresh_a2780_output_summary!(; start)
html_path = SW.render_a2780_report_html(; start)
println("untreated_rows=", nrow(untreated_fit.ranking))
println("linked_rows=", nrow(linked_ranking))
println("overview_path=", summary.overview_path)
println("html_path=", html_path)
