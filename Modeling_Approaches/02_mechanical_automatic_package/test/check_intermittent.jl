# Check: the engine-based intermittent prediction (zoo_analyses.jl) reproduces the previous
# hand-written simulator on the same fits (values stored in outputs/zoo_analyses/analyses_state.jls).
const PK = normpath(joinpath(@__DIR__, ".."))
include(joinpath(PK, "zoo_common.jl")); include(joinpath(PK, "zoo_models.jl")); include(joinpath(PK, "bootstrap_tools.jl"))
using Serialization, Printf, Test
state = deserialize(joinpath(PK, "outputs", "zoo_fits", "zoo_state.jls"))
gid(id) = state.groups[findfirst(g -> g.id == id, state.groups)]
old = deserialize(joinpath(PK, "outputs", "zoo_analyses", "analyses_state.jls"))[:intermittent]
ck = Dict{Symbol,Any}(); save_ck() = nothing
src = read(joinpath(PK, "zoo_analyses.jl"), String)
i = findfirst("# ── intermittent dosing: 1 µM", src)[1]; j = findfirst("# ── HTML", src)[1]
include_string(Main, src[i:j-1])
@testset "engine reproduces the intermittent predictions" begin
    for L in ("A2780", "A2780cis")
        a, b = ck[:intermittent][L], old[L]
        @printf("%s  BIC equal %.6f vs %.6f   k_off %.4g vs %.4g   BIC k_off %.6f vs %.6f\n", L, a.bic_equal, b.bic_equal, a.koff_best, b.koff_best, a.bic_koff, b.bic_koff)
        @test isapprox(a.sse_equal, b.sse_equal; rtol=1e-6)
        @test isapprox(a.sse_koff, b.sse_koff; rtol=1e-6) && a.koff_best == b.koff_best
        for sc in (:onoff, :offon)
            # the scale sits in a flat optimum (golden-section search), so 1e-5; the SSE agrees to ~1e-11
            @test isapprox(a.fits[sc].scale, b.fits[sc].scale; rtol=1e-5) && isapprox(a.fits[sc].pred, b.fits[sc].pred; rtol=1e-5)
            @test isapprox(a.fits[sc].sse, b.fits[sc].sse; rtol=1e-9)
            @test a.xy[sc][1] == b.xy[sc][1] && a.xy[sc][2] == b.xy[sc][2]
        end
    end
end
