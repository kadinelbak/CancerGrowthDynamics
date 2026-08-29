using MechanicalAutomaticModeling

const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))
const REPOSITORY_ROOT = normpath(joinpath(PACKAGE_ROOT, "..", ".."))
const DOCS_PACKAGE_ROOT = joinpath(REPOSITORY_ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package")

report = MechanicalAutomaticModeling.ReducedStageReport.render_reduced_stage_report_html(
    start = PACKAGE_ROOT,
    mirror_directory = DOCS_PACKAGE_ROOT,
)
println("Wrote reduced-stage report to: ", report)
