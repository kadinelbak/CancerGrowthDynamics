# The checks run_pipeline.jl does after a run, each in its own process:
#   test_stages.jl             stages.toml rebuilds the candidates of the saved run
#   test_validation_config.jl  validation.toml selects the intended data
#   test_control.jl            plant, linearisation, adjoint and MPC of the 14-day control
#   test_control_fasting.jl    the dose + fasting plant and controller
#   test_regression.jl         the saved run matches test/reference_results.toml (fails without a matching reference)
#   zoo_parity.js              the website's simulator draws the fitted curves (node)
#   julia -t 8 --project=. test/run_tests.jl
# The test files are independent, so they run at the same time (user, 2026-11-09; before, one after another took
# ~22 min); each one's output is printed in full, in the order above, once it finishes. TESTS_SERIAL=1 runs them
# one after another as before.
const T = @__DIR__
const FILES = ("test_stages.jl", "test_validation_config.jl", "test_control.jl", "test_control_fasting.jl", "test_regression.jl")
failed = String[]
run_one(f) = begin
    buf = IOBuffer(); t0 = time()
    ok = success(pipeline(`$(Base.julia_cmd()) --startup-file=no -t $(Threads.nthreads()) --project=$(dirname(T)) $(joinpath(T, f))`; stdout=buf, stderr=buf))
    (ok=ok, out=String(take!(buf)), secs=time() - t0)
end
present = [f for f in FILES if isfile(joinpath(T, f))]
for f in FILES
    f in present || (println("\n== $f MISSING"); push!(failed, f))
end
if get(ENV, "TESTS_SERIAL", "0") == "1"
    results = Dict(f => (println("\n== $f (running)"); flush(stdout); run_one(f)) for f in present)
else
    tasks = Dict(f => @async(run_one(f)) for f in present)
    results = Dict(f => fetch(t) for (f, t) in tasks)
end
for f in present
    r = results[f]
    println("\n== $f ($(round(Int, r.secs)) s)"); print(r.out); flush(stdout)
    r.ok || push!(failed, f)
end
# the website draws the same models as the fits (site.js against the curves in bundle.js)
println("\n== zoo_parity.js"); flush(stdout)
success(pipeline(`node $(joinpath(T, "zoo_parity.js"))`; stdout=stdout, stderr=stderr)) || push!(failed, "zoo_parity.js")
isempty(failed) ? println("\nAll test files passed.") : (println("\nFailed: $(join(failed, ", "))"); exit(1))
