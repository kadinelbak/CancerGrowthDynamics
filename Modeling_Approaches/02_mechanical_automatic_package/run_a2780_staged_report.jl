using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

max_time = parse(Float64, get(ENV, "A2780_MAX_TIME_PER_FIT", "12.0"))
result = MechanicalAutomaticModeling.StagedA2780Workflow.run_a2780_staged_goal!(
    start = joinpath(@__DIR__, "notebooks"),
    max_time_per_fit = max_time,
)
html_path = MechanicalAutomaticModeling.StagedA2780Workflow.render_a2780_report_html(
    start = joinpath(@__DIR__, "notebooks"),
)

println("overview_path=", result.overview_path)
println("manifest_path=", result.manifest_path)
println("html_path=", html_path)
