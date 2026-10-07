# stages.toml must rebuild every candidate of the last zoo run exactly: names, labels, free
# parameters, bounds, starts, initial states, fixed values, crowding and inputs. DRYRUN builds each
# stage without fitting, inheriting the fitted rows of outputs/zoo_fits/zoo_state.jls.
#   julia -t 4 --project=. test/test_stages.jl
ENV["DRYRUN"] = "1"
include(joinpath(@__DIR__, "..", "run_zoo_fits.jl"))
using Test

# slots are compared on the keys the saved line has: lines saved before 2026-10-02 have no EC50 slot
# (they used EC50 = 1 µM for both lines); the new EC50 slot is checked against stages.toml below
same_slots(a, b) = all(isequal(getfield(a, k), getfield(b, k)) for k in keys(b))
same_line(a, b) = all(k in (:s, :extra) || isequal(getfield(a, k), getfield(b, k)) for k in keys(b)) && same_slots(a.s, b.s) && keys(a) == keys(b) &&
                  (!hasproperty(b, :extra) || all(same_line(x, y) for (x, y) in zip(a.extra, b.extra)))
# starting states agree except where the saved run estimated a starting value (x0est), which the built candidate
# holds at the nominal seeding; the estimated indices and nominal values must agree too
function same_u0(c, p)
    est = hasproperty(p, :x0est) ? p.x0est : Int[]
    keepi = setdiff(eachindex(p.u0), est)
    length(c.u0) == length(p.u0) && isequal(c.u0[keepi], p.u0[keepi]) &&
        (!hasproperty(p, :x0est) || (c.x0est == p.x0est && isequal(c.x0nom, p.x0nom)))
end
function same_cand(c, p)
    ok = c.name == p.name && c.label == p.label && c.family == p.family && c.pnames == p.pnames &&
         isequal(c.bounds, p.bounds) && isequal(c.p0s, p.p0s) && same_u0(c, p) &&
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
    # every Hill kill line carries its line's measured IC50 (stages.toml [lines])
    for g in DRY_PREV.groups, c in DRY_BUILT[g.id]
        g.kind == :s4 && continue                       # Stage 4 lines are frozen copies of the saved Stage 3a lines
        lines = g.kind in (:s3, :s3b) ? ["A2780", "A2780cis"] : [g.line]
        for (L, name) in zip(c.sys.lines, lines)
            @test L.s.EC50 == Slot(0, LINE_CFG[name]["ec50"])
        end
    end
end
