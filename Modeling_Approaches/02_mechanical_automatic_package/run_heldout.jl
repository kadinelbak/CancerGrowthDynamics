# The held-out treated-coculture test (validation.toml [heldout], the same code as the Validation
# page) on any zoo run, e.g. a sensitivity variant fitted with ZOO_OUT:
#   ZOO_OUT=<dir> [STAGES_TOML=… SOURCES_TOML=… DATA_TABLE=…] julia -t 8 --project=. run_heldout.jl
# Reads <dir>/zoo_state.jls (default outputs/zoo_fits) and writes <dir>/heldout_ce1.csv: the prediction
# of Ce1 from crowding fitted to Ce0 with the Stage 2b kills, the in-sample fit and the untreated reference.
include(joinpath(@__DIR__, "zoo_common.jl")); include(joinpath(@__DIR__, "zoo_models.jl")); include(joinpath(@__DIR__, "bootstrap_tools.jl"))
using Serialization, Statistics
const GPE = GrowthParameterEstimation
const RUN_DIR = get(ENV, "ZOO_OUT", joinpath(@__DIR__, "outputs", "zoo_fits"))
state = deserialize(joinpath(RUN_DIR, "zoo_state.jls"))
const groups = state.groups
gid(id) = groups[findfirst(g -> g.id == id, groups)]
include(joinpath(@__DIR__, "zoo_heldout.jl"))
r = heldout_ce1()
m1 = mean(reduce(vcat, [d.y for d in g3.datasets[CE1]]))
p = r.ce1[1]
rows = [(row="prediction (crowding from Ce0, kills from Stage 2b)", rmse=p.pred.rmse, nrmse=p.pred.rmse / m1, n=nobs(g3.datasets[CE1]),
         params=join(["$(n)=$(v)" for (n, v) in zip(c3.pnames[HO_RI], p.params)], "; ")),
        (row="in-sample (Stage 3a fitted with Ce1)", rmse=r.ce1_full.rmse, nrmse=r.ce1_full.rmse / m1, n=nobs(g3.datasets[CE1]), params=""),
        (row="untreated reference (same crowding, no drug)", rmse=r.ce1_baseline_rmse, nrmse=r.ce1_baseline_rmse / m1, n=nobs(g3.datasets[CE1]), params="")]
CSV.write(joinpath(RUN_DIR, "heldout_ce1.csv"), DataFrame(rows))
println("Stage 3a model $(c3.name): held-out Ce1 RMSE $(round(p.pred.rmse; digits=1)) (nRMSE $(round(p.pred.rmse / m1; digits=3))); ",
        "in-sample $(round(r.ce1_full.rmse; digits=1)); untreated reference $(round(r.ce1_baseline_rmse; digits=1))")
println("Wrote $(joinpath(RUN_DIR, "heldout_ce1.csv"))")
