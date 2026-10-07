# The tidy table (data/measurements.csv via zoo_data.jl) must reproduce the hard-coded loaders it
# replaced (test/legacy_loaders.jl) exactly: every time, value and exclusion flag of every series.
#   julia --project=. test/test_data.jl
include(joinpath(@__DIR__, "..", "zoo_common.jl"))
include(joinpath(@__DIR__, "legacy_loaders.jl"))
using Test
const O = LegacyData

@testset "tidy data reproduces the legacy loaders" begin
    for L in ("A2780", "A2780cis")
        @test untreated[L].x == O.untreated[L].x && untreated[L].y == O.untreated[L].y
        # A2780cis low resource used day 1 as t = 0 before 2026-10-07; now every series starts at day 0
        sh = L == "A2780cis" && first(O.lrmono[L].x) < 0 ? 1.0 : (L == "A2780cis" && length(O.lrmono[L].x) < length(lrmono[L].x) ? 1.0 : 0.0)
        newd = Dict(zip(lrmono[L].x, lrmono[L].y))
        @test all(get(newd, t + sh, NaN) == v || get(newd, t, NaN) == v for (t, v) in zip(O.lrmono[L].x, O.lrmono[L].y))
        @test treated[L].doses == O.treated[L].doses && treated[L].files == O.treated[L].files
        @test all(a[1] == b[1] && a[2] == b[2] for (a, b) in zip(treated[L].data, O.treated[L].data))
    end
    @test length(co_data) == length(O.co_data) == 6
    for (a, b, la, lb) in zip(co_data, O.co_data, COCULTURE_LAYOUT, O.COCULTURE_LAYOUT)
        @test a.S[1] == b.S[1] && a.S[2] == b.S[2] && a.R[1] == b.R[1] && a.R[2] == b.R[2]
        @test (la.ratio, la.cond, la.file, la.ssheet, la.rsheet, la.dose) == (lb.ratio, lb.cond, lb.file, lb.ssheet, lb.rsheet, lb.dose)
    end
    @test length(lr_cc) == length(O.lr_cc) == 6
    for (a, b) in zip(lr_cc, O.lr_cc)
        @test (a.run, a.ratio, a.label, a.sN, a.sC) == (b.run, b.ratio, b.label, b.sN, b.sC)
        @test a.xN == b.xN && a.yN == b.yN && a.xC == b.xC && a.yC == b.yC
        # QC since 2026-10-07: rules for every source, symmetric, both lines of a well; never fewer exclusions
        @test all(a.bN .| .!b.bN) && all(a.bC .| .!b.bC)
    end
    # Stage 4 now starts at day 0 (the seeding) instead of day 1: the legacy series are the new ones from day 1 on
    for (k, kx) in ((:yN, :xN), (:yC, :xC))
        new = Dict(zip(getfield(combo, kx), getfield(combo, k))); old = Dict(zip(getfield(O.combo, kx) .+ 1, getfield(O.combo, k)))
        @test all(new[t] == v for (t, v) in old)
    end
    @test combo.xN[1] == 0.0 && combo.yN[1] ≈ 49.9
    @test COMBO_DOSE == O.COMBO_DOSE && COCULTURE_DOSE == O.COCULTURE_DOSE
    # intermittent: the four day-average files, days 1–14, dose switched at day 7
    @test sort(collect(keys(intermittent["A2780"]))) == ["off0_on7", "on0_off7"]
    @test all(length(v.x) == 14 for L in values(intermittent) for v in values(L))
end
