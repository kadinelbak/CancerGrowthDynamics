using MechanicalAutomaticModeling

const PACKAGE_ROOT = normpath(joinpath(@__DIR__, ".."))
fit_seconds = parse(Float64, get(ENV, "OPTIMAL_CONTROL_ONE_FIT_SECONDS", "8"))
control_seconds = parse(Float64, get(ENV, "OPTIMAL_CONTROL_ONE_CONTROL_SECONDS", "20"))

result = MechanicalAutomaticModeling.OptimalControlOne.run_optimal_control_one!(
    start = PACKAGE_ROOT,
    fit_seconds = fit_seconds,
    control_seconds = control_seconds,
)

println("report=", result.report)
println("docs_report=", result.docs_report)
println("output=", result.output)
