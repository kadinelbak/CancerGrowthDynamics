using MechanicalAutomaticModeling

const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))
const REPOSITORY_ROOT = normpath(joinpath(PACKAGE_ROOT, "..", ".."))
const DOCS_PACKAGE_ROOT = joinpath(REPOSITORY_ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package")

report = MechanicalAutomaticModeling.SampleAwareReport.render_a2780_sample_report_html(
    start = joinpath(PACKAGE_ROOT, "notebooks"),
    mirror_directory = DOCS_PACKAGE_ROOT,
)

println("Wrote sample-aware staged report to:")
println("  ", report)
println("  ", joinpath(DOCS_PACKAGE_ROOT, "outputs", "reports", basename(report)))
