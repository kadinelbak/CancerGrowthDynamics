# stages.toml must rebuild every candidate of the last zoo run exactly: names, labels, free
# parameters, bounds, starts, initial states, fixed values, crowding and inputs. DRYRUN builds each
# stage without fitting, inheriting the fitted rows of outputs/zoo_fits/zoo_state.jls.
#   julia -t 4 --project=. test/test_stages.jl
ENV["DRYRUN"] = "1"
include(joinpath(@__DIR__, "..", "run_zoo_fits.jl"))
using Test

same_line(a, b) = all(isequal(getfield(a, k), getfield(b, k)) for k in keys(b) if k != :s) && isequal(a.s, b.s) && keys(a) == keys(b)
function same_cand(c, p)
    ok = c.name == p.name && c.label == p.label && c.family == p.family && c.pnames == p.pnames &&
         isequal(c.bounds, p.bounds) && isequal(c.p0s, p.p0s) && isequal(c.u0, p.u0) &&
         length(c.sys.lines) == length(p.sys.lines) && all(same_line(a, b) for (a, b) in zip(c.sys.lines, p.sys.lines)) &&
         isequal(c.sys.alpha, p.sys.alpha) && isequal(c.sys.envs, p.sys.envs) && c.sys.xi == p.sys.xi && c.sys.hb == p.sys.hb
    for k in (:slots, :fixedcc, :noinput, :crowd, :reskill, :βC2)
        hasproperty(p, k) && (ok &= hasproperty(c, k) && isequal(getfield(c, k), getfield(p, k)))
    end
    hasproperty(p, :mech) && (ok &= [m.name for m in c.mech] == [m.name for m in p.mech])
    ok || @warn "candidate differs" c.name
    ok
end

@testset "stages.toml rebuilds the last run's candidates" begin
    @test Set(keys(DRY_BUILT)) == Set(g.id for g in DRY_PREV.groups)
    for g in DRY_PREV.groups
        built = Dict(c.name => c for c in DRY_BUILT[g.id])
        @test Set(keys(built)) == Set(c.name for c in g.cands)
        @test all(same_cand(built[p.name], p) for p in g.cands)
    end
    @test [g.id for g in groups] == [g.id for g in DRY_PREV.groups]
    @test all(a.stage == b.stage && a.title == b.title && a.line == b.line && a.kind == b.kind && isequal(a.datasets, b.datasets)
              for (a, b) in zip(groups, DRY_PREV.groups))
    @test JOINT_RULES == [("s1b", "ab"), ("s2", "shared")]
end
