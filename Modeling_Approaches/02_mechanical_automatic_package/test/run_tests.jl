# The checks run_pipeline.jl does after a run, each in its own process:
#   test_stages.jl             stages.toml rebuilds the candidates of the saved run
#   test_validation_config.jl  validation.toml selects the intended data
#   test_control.jl            plant, linearisation, adjoint and MPC of the 14-day control
#   test_control_fasting.jl    the dose + fasting plant and controller
#   test_regression.jl         the saved run matches test/reference_results.toml (fails without a matching reference)
#   zoo_parity.js              the website's simulator draws the fitted curves (node)
#   julia -t 8 --project=. test/run_tests.jl
const T = @__DIR__
failed = String[]
for f in ("test_stages.jl", "test_validation_config.jl", "test_control.jl", "test_control_fasting.jl", "test_regression.jl")
    isfile(joinpath(T, f)) || (println("\n== $f MISSING"); push!(failed, f); continue)
    println("\n== $f"); flush(stdout)
    ok = success(pipeline(`$(Base.julia_cmd()) --startup-file=no -t $(Threads.nthreads()) --project=$(dirname(T)) $(joinpath(T, f))`; stdout=stdout, stderr=stderr))
    ok || push!(failed, f)
end
# the website draws the same models as the fits (site.js against the curves in bundle.js)
println("\n== zoo_parity.js"); flush(stdout)
success(pipeline(`node $(joinpath(T, "zoo_parity.js"))`; stdout=stdout, stderr=stderr)) || push!(failed, "zoo_parity.js")
isempty(failed) ? println("\nAll test files passed.") : (println("\nFailed: $(join(failed, ", "))"); exit(1))
