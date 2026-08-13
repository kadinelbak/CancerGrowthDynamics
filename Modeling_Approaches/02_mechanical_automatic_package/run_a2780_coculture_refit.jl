using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames

include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

start = joinpath(@__DIR__, "notebooks")
max_time = parse(Float64, get(ENV, "A2780_MAX_TIME_PER_FIT", "1.5"))

for condition in ("coculture_untreated", "coculture_treated")
    decoded = MechanicalAutomaticModeling.StagedA2780Workflow.decode_a2780_condition(condition; start = start)
    out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs(condition; start = start)
    CSV.write(joinpath(out.csv, "$(condition)_a2780_decoded.csv"), decoded)
    fit = MechanicalAutomaticModeling.FitWorkflows.run_condition_fit!(
        decoded,
        condition;
        start = start,
        max_time_per_fit = max_time,
    )
    ranking = fit.ranking
    all(isfinite, Float64.(ranking.bic)) || error("$(condition) has non-finite BIC")
    all(isfinite, Float64.(ranking.ssr)) || error("$(condition) has non-finite raw SSR")
    all(Float64.(ranking.ssr) .< 9.99e11) || error("$(condition) contains failure-sentinel SSR")
    println(condition, "_rows=", nrow(ranking))
end

summary = MechanicalAutomaticModeling.StagedA2780Workflow.refresh_a2780_output_summary!(; start = start)
html_path = MechanicalAutomaticModeling.StagedA2780Workflow.render_a2780_report_html(; start = start)
println("overview_path=", summary.overview_path)
println("html_path=", html_path)
