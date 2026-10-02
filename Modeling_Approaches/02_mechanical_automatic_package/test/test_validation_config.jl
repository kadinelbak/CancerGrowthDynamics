# validation.toml must select exactly the data the hard-coded validation used: the Ce0 / Ce1
# datasets of Stage 3a and the held-out label of every dataset of each leave-one-out block.
#   julia --project=. test/test_validation_config.jl
include(joinpath(@__DIR__, "..", "zoo_common.jl")); include(joinpath(@__DIR__, "..", "zoo_models.jl"))
using Serialization, Test
state = deserialize(joinpath(@__DIR__, "..", "outputs", "zoo_fits", "zoo_state.jls"))
gid(id) = state.groups[findfirst(g -> g.id == id, state.groups)]
VAL = TOML.parsefile(joinpath(@__DIR__, "..", "validation.toml"))
tag_match(s, sel) = all(string(get(s, k, "")) == string(v) for (k, v) in sel)
tag_label(s, tag) = tag == "schedule" ? "$(constant_dose(s["schedule"])) µM" : tag == "run" ? "run $(s["run"])" : string(s[tag])
@testset "validation.toml reproduces the hard-coded selections" begin
    g3 = gid("s3"); T = dataset_tags("s3")
    @test length(T) == length(g3.datasets)
    @test [k for (k, s) in enumerate(T) if tag_match(s, VAL["heldout"]["train"])] == [k for k in eachindex(g3.datasets) if COCULTURE_LAYOUT[(k + 1) ÷ 2].cond == "Ce0"]
    @test [k for (k, s) in enumerate(T) if tag_match(s, VAL["heldout"]["predict"])] == [k for k in eachindex(g3.datasets) if COCULTURE_LAYOUT[(k + 1) ÷ 2].cond == "Ce1"]
    old = Dict(("s2-a2780", "leave one dose out") => ["$(d) µM" for d in treated["A2780"].doses],
               ("s2-a2780cis", "leave one dose out") => ["$(d) µM" for d in treated["A2780cis"].doses],
               ("s3", "leave one ratio out") => [COCULTURE_LAYOUT[(k + 1) ÷ 2].ratio for k in eachindex(g3.datasets)],
               ("s3b", "leave one run out") => ["run $(lr_cc[(k + 1) ÷ 2].run)" for k in eachindex(gid("s3b").datasets)],
               ("s3b", "leave one ratio out") => [lr_cc[(k + 1) ÷ 2].label for k in eachindex(gid("s3b").datasets)])
    @test length(VAL["leave_out"]) == length(old)
    for lo in VAL["leave_out"]
        @test [tag_label(s, lo["tag"]) for s in dataset_tags(lo["group"])] == old[(lo["group"], lo["scheme"])]
    end
    # every dataset's tags give back its own series (same values, same order)
    for g in state.groups
        g.kind == :s4 && continue
        for (d, s) in zip(g.datasets, dataset_tags(g.id))
            S = series(s["id"]); @test d.x == S.x[.!S.bad] && d.y == S.y[.!S.bad]
        end
    end
end
