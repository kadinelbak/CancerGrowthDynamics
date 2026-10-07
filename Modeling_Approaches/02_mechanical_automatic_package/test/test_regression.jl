# Regression test against test/reference_results.toml (made by test/make_reference.jl).
#   julia -t 8 --project=. test/test_regression.jl          refit Stages 1 and 2b and the first bootstrap
#                                                           replicates from scratch and compare the saved run
#                                                           and the control tables with the reference
#   REGRESSION=quick …   only Stage 1;   REGRESSION=full …   also Stage 2a
# Fits stop on iteration limits only (no wall-clock limit), so the numbers must agree whatever the machine or
# thread count: BICs to 1e-6 relative. The test FAILS when there is no reference, or when the reference was made
# with a different configuration or data (fingerprint): after an intended change, accept the new run with
# test/make_reference.jl and rerun this test, which then checks that a fresh refit reproduces it.
const PK = normpath(joinpath(@__DIR__, ".."))
REF = joinpath(@__DIR__, "reference_results.toml")
isfile(REF) || error("no reference (test/reference_results.toml): run test/make_reference.jl after an accepted full run")
include(joinpath(@__DIR__, "regression_fingerprint.jl"))
ENV["DRYRUN"] = "1"                       # build the candidates exactly as the pipeline does, without fitting
include(joinpath(PK, "run_zoo_fits.jl"))
using Test, TOML
ref = TOML.parsefile(REF)
const RT = 1e-6
close_to(a, b) = isapprox(a, b; rtol=RT, atol=1e-8)

using CSV, DataFrames
@testset "regression against the reference run" begin
    @testset "reference made with this configuration and data" begin
        fp = get(ref, "fingerprint", Dict()); now_fp = config_fingerprint()
        for f in FINGERPRINT_FILES
            @test get(fp, f, "") == now_fp[f]
        end
    end
    @testset "saved run = reference" begin
        for g in DRY_PREV.groups
            R = ref["groups"][g.id]
            @test g.cands[1].name == R["carried"]
            @test all(close_to(c.bic, R["bic"][c.name]) for c in g.cands)
            @test all(close_to(v, R["params"][n]) for (n, v) in zip(g.cands[1].pnames, g.cands[1].params))
        end
    end
    mode = get(ENV, "REGRESSION", "default")
    refit = [g for g in DRY_PREV.groups if startswith(g.id, "s1-") || (mode != "quick" && startswith(g.id, "s2-")) || (mode == "full" && startswith(g.id, "s1b-"))]
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
    # the control outcome tables of the saved run (deterministic: exact search, fixed RK4)
    LONG = round(Int, TOML.parsefile(joinpath(PK, "control_config.toml"))["schedule"]["long_days"])
    @testset "control table $name" for (name, f) in (("control14", "outputs/control/control_outcomes.csv"), ("control_long", "outputs/control_$(LONG)d/fasting_outcomes.csv"))
        haskey(ref, name) || continue
        df = CSV.read(joinpath(PK, f), DataFrame)
        @test nrow(df) == length(ref[name])
        for (row, R) in zip(eachrow(df), ref[name]), (k, v) in R
            x = row[Symbol(k)]
            @test v isa Number ? close_to(Float64(x), v) : string(x) == v
        end
    end
end
