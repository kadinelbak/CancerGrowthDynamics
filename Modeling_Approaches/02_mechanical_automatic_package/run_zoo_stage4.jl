# Adds Stage 4 (low resource + cisplatin coculture, zoo_stage4.jl) to an existing zoo run without
# refitting Stages 1–3b: loads outputs/zoo_fits/zoo_state.jls, fits and bootstraps Stage 4 from the
# carried Stage 3 / 3b models, appends it (replacing an older Stage 4 if present) and saves the state
# (the previous file is kept as zoo_state_before_s4.jls). Then run zoo_rebuild_pages.jl.
#   julia -t 12 --project=. run_zoo_stage4.jl
# run_zoo_fits.jl runs the same Stage 4 code in a full run.

include(joinpath(@__DIR__, "zoo_common.jl"))
include(joinpath(@__DIR__, "zoo_models.jl"))
include(joinpath(@__DIR__, "bootstrap_tools.jl"))
using Serialization
const ZOO_DIR = joinpath(@__DIR__, "outputs", "zoo_fits")
state = deserialize(joinpath(ZOO_DIR, "zoo_state.jls"))
const BOOT_B = state.B
include(joinpath(@__DIR__, "zoo_selection.jl"))
include(joinpath(@__DIR__, "zoo_stage4.jl"))

gid(id) = state.groups[findfirst(g -> g.id == id, state.groups)]
c3, c3b = gid("s3").cands[1], gid("s3b").cands[1]
βC2 = gid("s2-a2780cis").cands[1].params[1]
println("=== Stage 4 add-on: $(Threads.nthreads()) threads, B = $BOOT_B ==="); flush(stdout)
t0 = time()
g4 = stage4_group(c3, c3b, βC2)
for r in g4.cands
    @printf("  %-22s k=%d BIC=%8.2f  %s\n", r.name, length(r.params), r.bic, join(["$n = $(fmt(v))" for (n, v) in zip(r.pnames, r.params)], ", "))
end
flush(stdout)
r4 = stage4_results(g4; B=BOOT_B)
@printf("  s4 optimal %s, chosen in %.1f%% of replicates (lowest BIC in %.1f%%)  [%.0fs]\n", g4.cands[1].name, pct(r4.res.pick, 1), pct(r4.res.winner, 1), time() - t0)

groups = vcat([g for g in state.groups if g.id != "s4"], [g4])
results = copy(state.results); results["s4"] = r4
cp(joinpath(ZOO_DIR, "zoo_state.jls"), joinpath(ZOO_DIR, "zoo_state_before_s4.jls"); force=true)
serialize(joinpath(ZOO_DIR, "zoo_state.jls"), merge(state, (groups=groups, results=results)))

# CSV rows for Stage 4, appended to the zoo selection / parameter tables
sel = CSV.read(joinpath(ZOO_DIR, "zoo_model_selection.csv"), DataFrame); par = CSV.read(joinpath(ZOO_DIR, "zoo_parameters.csv"), DataFrame)
sel = sel[string.(sel.group) .!= "s4", :]; par = par[string.(par.group) .!= "s4", :]
for (j, c) in enumerate(g4.cands)
    push!(sel, (group="s4", stage=g4.stage, candidate=c.name, label=c.label, family=c.family, k=length(c.params), bic=c.bic,
                delta_bic=c.bic - g4.cands[1].bic, pct_lowest_bic=pct(r4.res.winner, j), pct_chosen=pct(r4.res.pick, j),
                optimal=(j == 1), params=join(["$n=$(v)" for (n, v) in zip(c.pnames, c.params)], "; ")); cols=:subset)
end
c = g4.cands[1]
for (i, n) in enumerate(c.pnames)
    l, m, h = qs(r4.res.pars[1][:, i])
    push!(par, (group="s4", candidate=c.name, param=n, estimate=c.params[i], se=c.se[i], boot_median=m, boot_lo=l, boot_hi=h); cols=:subset)
end
CSV.write(joinpath(ZOO_DIR, "zoo_model_selection.csv"), sel); CSV.write(joinpath(ZOO_DIR, "zoo_parameters.csv"), par)
println("Saved Stage 4 into zoo_state.jls. Next: julia --project=. zoo_rebuild_pages.jl")
