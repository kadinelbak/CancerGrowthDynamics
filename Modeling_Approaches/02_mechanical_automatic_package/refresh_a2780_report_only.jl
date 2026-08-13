using Pkg

Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

start = joinpath(@__DIR__, "notebooks")
summary = MechanicalAutomaticModeling.StagedA2780Workflow.refresh_a2780_output_summary!(; start = start)
html_path = MechanicalAutomaticModeling.StagedA2780Workflow.render_a2780_report_html(; start = start)

println("overview_path=", summary.overview_path)
println("html_path=", html_path)
