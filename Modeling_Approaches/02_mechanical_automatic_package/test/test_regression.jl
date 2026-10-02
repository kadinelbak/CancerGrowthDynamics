# Regression test against test/reference_results.toml (made by test/make_reference.jl).
#   julia -t 8 --project=. test/test_regression.jl          fast: refit Stage 1 and the first bootstrap
#                                                           replicates from scratch (a few minutes) and
#                                                           compare the saved run with the reference
#   REGRESSION=full julia … test/test_regression.jl          also refit Stages 2a and 2b
# Fits stop on iteration limits only (no wall-clock limit), so the numbers must agree whatever the
# machine or thread count: BICs to 1e-6 relative.
const PK = normpath(joinpath(@__DIR__, ".."))
REF = joinpath(@__DIR__, "reference_results.toml")
if !isfile(REF)
    println("No reference yet (test/reference_results.toml); run test/make_reference.jl after an accepted full run.")
    exit(0)
end
ENV["DRYRUN"] = "1"                       # build the candidates exactly as the pipeline does, without fitting
include(joinpath(PK, "run_zoo_fits.jl"))
using Test, TOML
ref = TOML.parsefile(REF)
const RT = 1e-6
close_to(a, b) = isapprox(a, b; rtol=RT, atol=1e-8)

@testset "regression against the reference run" begin
    @testset "saved run = reference" begin
        for g in DRY_PREV.groups
            R = ref["groups"][g.id]
            @test g.cands[1].name == R["carried"]
            @test all(close_to(c.bic, R["bic"][c.name]) for c in g.cands)
            @test all(close_to(v, R["params"][n]) for (n, v) in zip(g.cands[1].pnames, g.cands[1].params))
        end
    end
    full = get(ENV, "REGRESSION", "fast") == "full"
    refit = [g for g in DRY_PREV.groups if startswith(g.id, "s1-") || (full && (startswith(g.id, "s1b-") || startswith(g.id, "s2-")))]
    @testset "refit from scratch reproduces $(g.id)" for g in refit
        rows = fit_zoo(DRY_BUILT[g.id], g.datasets)
        R = ref["groups"][g.id]
        @test all(close_to(r.bic, R["bic"][r.name]) for r in rows)
    end
    @testset "first bootstrap replicates of $(g.id)" for g in filter(g -> startswith(g.id, "s1-"), DRY_PREV.groups)
        R = ref["groups"][g.id]["boot_bic"]; B = length(R)
        res = bootstrap(merge(g, (cands=g.cands,)); B=B)
        @test all(close_to(res.bic[b, j], R[b][c.name]) for b in 1:B for (j, c) in enumerate(res.cands))
    end
end
