using Pkg

Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

start = joinpath(@__DIR__, "notebooks")
summary = MechanicalAutomaticModeling.StagedA2780Workflow.refresh_a2780_output_summary!(; start = start)
html_path = MechanicalAutomaticModeling.StagedA2780Workflow.render_a2780_report_html(; start = start)
repository_root = normpath(joinpath(@__DIR__, "..", ".."))
docs_package_root = joinpath(repository_root, "docs", "Modeling_Approaches", "02_mechanical_automatic_package")
docs_report_dir = joinpath(docs_package_root, "outputs", "reports")
docs_image_root = joinpath(docs_package_root, "outputs", "images")
mkpath(docs_report_dir)
mkpath(docs_image_root)
cp(html_path, joinpath(docs_report_dir, basename(html_path)); force = true)
figure_relative_paths = (
    joinpath("monoculture_untreated", "figures", "monoculture_untreated_pooling_model_grid.png"),
    joinpath("monoculture_treated", "figures", "monoculture_treated_best_joint_model_by_environment.png"),
    joinpath("monoculture_treated", "figures", "monoculture_treated_timing_hypothesis_grid.png"),
    joinpath("coculture_untreated", "figures", "coculture_untreated_best_mechanistic_fit_grid.png"),
    joinpath("coculture_treated", "figures", "linked_treatment_coculture_grid.png"),
)
for relative_path in figure_relative_paths
    source = joinpath(@__DIR__, "outputs", "images", relative_path)
    isfile(source) || error("Missing staged report figure: $(source)")
    destination = joinpath(docs_image_root, relative_path)
    mkpath(dirname(destination))
    cp(source, destination; force = true)
end

println("overview_path=", summary.overview_path)
println("html_path=", html_path)
println("docs_html_path=", joinpath(docs_report_dir, basename(html_path)))
