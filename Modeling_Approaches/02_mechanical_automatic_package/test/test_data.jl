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
        @test lrmono[L].x == O.lrmono[L].x && lrmono[L].y == O.lrmono[L].y && lrmono[L].bad == O.lrmono[L].bad && lrmono[L].why == O.lrmono[L].why
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
        @test a.xN == b.xN && a.yN == b.yN && a.bN == b.bN && a.xC == b.xC && a.yC == b.yC && a.bC == b.bC
    end
    for k in (:xN, :yN, :bN, :xC, :yC, :bC)
        @test getfield(combo, k) == getfield(O.combo, k)
    end
    @test COMBO_DOSE == O.COMBO_DOSE && COCULTURE_DOSE == O.COCULTURE_DOSE
    # intermittent: the four day-average files, days 1–14, dose switched at day 7
    @test sort(collect(keys(intermittent["A2780"]))) == ["off0_on7", "on0_off7"]
    @test all(length(v.x) == 14 for L in values(intermittent) for v in values(L))
end
