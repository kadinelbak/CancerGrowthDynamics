using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames

include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

const FW = MechanicalAutomaticModeling.FitWorkflows
const SW = MechanicalAutomaticModeling.StagedA2780Workflow

start = joinpath(@__DIR__, "notebooks")
max_time = parse(Float64, get(ENV, "A2780_MAX_TIME_PER_FIT", "1.5"))
monoculture_decoded = SW.decode_a2780_condition("monoculture_treated"; start = start)
monoculture_out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs("monoculture_treated"; start = start)
CSV.write(joinpath(monoculture_out.csv, "monoculture_treated_a2780_decoded.csv"), monoculture_decoded)
decoded = SW.decode_a2780_condition("coculture_treated"; start = start)
out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs("coculture_treated"; start = start)
CSV.write(joinpath(out.csv, "coculture_treated_a2780_decoded.csv"), decoded)
recovered, provenance = FW._recover_coculture_design(decoded, "coculture_treated")
environments, initial = FW._coculture_environments(recovered)
CSV.write(joinpath(out.csv, "coculture_treated_density_provenance_validation.csv"), provenance)
CSV.write(joinpath(out.csv, "coculture_treated_initial_mix_diagnostics.csv"), initial)

ranking = FW._fit_linked_treatment_joint(
    environments,
    out;
    start = start,
    max_time_per_fit = max_time,
    resume_from_existing = false,
)
all(isfinite, Float64.(ranking.bic)) || error("Linked treatment ranking contains non-finite BIC")
all(isfinite, Float64.(ranking.ssr)) || error("Linked treatment ranking contains non-finite SSR")
CSV.write(joinpath(out.csv, "coculture_treated_automatic_model_ranking.csv"), ranking)
# The filename is legacy, but this report-facing artifact intentionally shows
# the same top-five linked hypotheses as linked_treatment_top5.csv.
CSV.write(joinpath(out.csv, "coculture_treated_automatic_best_models_top10.csv"), first(ranking, min(FW.LINKED_TREATMENT_VISIBLE_LIMIT, nrow(ranking))))

summary = SW.refresh_a2780_output_summary!(; start = start)
html_path = SW.render_a2780_report_html(; start = start)
println("linked_rows=", nrow(ranking))
println("canonical_monoculture_rows=", nrow(monoculture_decoded))
println("canonical_coculture_rows=", nrow(decoded))
println("overview_path=", summary.overview_path)
println("html_path=", html_path)
