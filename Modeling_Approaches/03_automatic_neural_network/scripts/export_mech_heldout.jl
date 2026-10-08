# Per-point export of the mechanistic zoo's held-out Ce1 prediction, for the three-way comparison
# (run_07_comparison.jl). It runs the zoo's own held-out code (02_mechanical_automatic_package/zoo_heldout.jl,
# the same as run_heldout.jl and the Validation page) on the zoo's latest fitted state and writes, into this
# track's outputs/comparison/:
#   mech_heldout_points.csv   one row per held-out Ce1 point: ratio, line, day, observed, predicted
#   mech_heldout_params.csv   the free parameters along the zoo's path to the prediction
# It only READS the zoo; run it in the zoo's environment from this folder:
#   julia -t 4 --project=../02_mechanical_automatic_package scripts/export_mech_heldout.jl
const PKG = normpath(joinpath(@__DIR__, "..", "..", "02_mechanical_automatic_package"))
include(joinpath(PKG, "zoo_common.jl")); include(joinpath(PKG, "zoo_models.jl")); include(joinpath(PKG, "bootstrap_tools.jl"))
using Serialization, Statistics, Dates
const GPE = GrowthParameterEstimation
const RUN_DIR = get(ENV, "ZOO_OUT", joinpath(PKG, "outputs", "zoo_fits"))
state = deserialize(joinpath(RUN_DIR, "zoo_state.jls"))
const groups = state.groups
gid(id) = groups[findfirst(g -> g.id == id, groups)]
include(joinpath(PKG, "zoo_heldout.jl"))

r = heldout_ce1()
p = r.ce1[1]
D = g3.datasets[CE1]
tags = HO_TAGS[CE1]
rows = DataFrame(source_id = String[], ratio = String[], line = String[], day = Float64[], observed = Float64[], predicted = Float64[])
for (k, d) in enumerate(D)
    for (x, y, yhat) in zip(d.x, d.y, p.pred.predictions[k])
        push!(rows, (tags[k]["id"], string(tags[k]["ratio"]), tags[k]["line"], Float64(x), Float64(y), Float64(yhat)))
    end
end
out = mkpath(joinpath(@__DIR__, "..", get(ENV, "NN_OUT_REL", "outputs"), "comparison"))
CSV.write(joinpath(out, "mech_heldout_points.csv"), rows)

# free parameters along the path: Stage 1 growth and Stage 2b drug response of each line (carried models),
# plus the crowding parameters refitted to Ce0 for the prediction; estimated starting values listed apart
prow = DataFrame(group = String[], model = String[], params = Int[], estimated_starts = Int[], names = String[])
for id in ("s1-a2780", "s1-a2780cis", "s2-a2780", "s2-a2780cis")
    any(g -> g.id == id, groups) || continue
    c = gid(id).cands[1]
    push!(prow, (id, String(c.name), length(c.params), hasproperty(c, :x0est) ? length(c.x0est) : 0, join(string.(c.pnames), ", ")))
end
push!(prow, ("s3 held-out refit (Ce0)", String(c3.name), length(HO_NAMES), 0, join(string.(HO_NAMES), ", ")))
CSV.write(joinpath(out, "mech_heldout_params.csv"), prow)
m1 = mean(reduce(vcat, [d.y for d in D]))
open(joinpath(out, "mech_heldout_source.txt"), "w") do io
    println(io, "exported $(now()) from $(joinpath(RUN_DIR, "zoo_state.jls")) (modified $(Dates.unix2datetime(mtime(joinpath(RUN_DIR, "zoo_state.jls")))) UTC)")
    println(io, "Stage 3a carried model: $(c3.name); held-out RMSE over all $(nrow(rows)) exported points = $(p.pred.rmse) (nRMSE $(p.pred.rmse / m1))")
end
println("held-out Ce1: $(nrow(rows)) points, RMSE $(round(p.pred.rmse; digits = 1)); wrote $(out)")
