# The held-out treated-coculture test (validation.toml [heldout], the same code as the Validation
# page) on any zoo run, e.g. a sensitivity variant fitted with ZOO_OUT:
#   ZOO_OUT=<dir> [STAGES_TOML=… SOURCES_TOML=… DATA_TABLE=…] julia -t 8 --project=. run_heldout.jl
# Reads <dir>/zoo_state.jls (default outputs/zoo_fits) and writes <dir>/heldout_ce1.csv (the same file
# zoo_analyses.jl writes for the main run): the prediction of Ce1 from crowding fitted (weighted) to Ce0
# with the Stage 2b kills and nominal Ce1 starting values, the in-sample fit, the untreated reference and
# a paired bootstrap of the RMSE difference over the Ce1 trajectories.
include(joinpath(@__DIR__, "zoo_common.jl")); include(joinpath(@__DIR__, "zoo_models.jl")); include(joinpath(@__DIR__, "bootstrap_tools.jl"))
using Serialization, Statistics
const GPE = GrowthParameterEstimation
const RUN_DIR = get(ENV, "ZOO_OUT", joinpath(@__DIR__, "outputs", "zoo_fits"))
state = deserialize(joinpath(RUN_DIR, "zoo_state.jls"))
const groups = state.groups
gid(id) = groups[findfirst(g -> g.id == id, groups)]
include(joinpath(@__DIR__, "zoo_heldout.jl"))
r = heldout_ce1()
write_heldout_csv(r, RUN_DIR)
m1 = mean(reduce(vcat, [d.y for d in g3.datasets[CE1]]))
p = r.ce1[1]; b = r.ce1_boot
println("Stage 3a model $(c3.name): held-out Ce1 RMSE $(round(p.pred.rmse; digits=1)) cells (nRMSE $(round(p.pred.rmse / m1; digits=3)), ",
        "weighted RMSE $(round(p.pred.wrmse; digits=3))); in-sample $(round(r.ce1_full.rmse; digits=1)); untreated reference $(round(r.ce1_baseline_rmse; digits=1))")
println("RMSE difference prediction − reference $(round(b.diff; digits=1)) cells (95% paired bootstrap $(round(b.lo; digits=1)) – $(round(b.hi; digits=1))); Ce1 starts: $(r.ce1_u0_note)")
