# Stage 2b: hold-out folds of Stage 2 (control-relevant generalisation). Each fold refits no drug, Hill, the
# UDE and the pure neural ODE (NN_SEEDS seeds each) without the held-out drug wells and predicts them from
# their nominal seeding. Run after Stage 7 (run_all.ps1 does), so the comparison does not wait for it.
# The learned-resource UDE is not refitted here (runtime: its fits cost about twice a UDE fit); it differs
# from the UDE only in low-resource wells, which only the 20k fold holds out (2 of its 12 wells).
# Fits are checkpointed in outputs/ude/checkpoints (shared with Stage 2) and always reused.
#   julia -t 6 --project=. scripts/run_03b_folds.jl
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Statistics, Serialization
const N = NNTrack
const FAST = get(ENV, "NN_FAST", "0") == "1"
const ITERS = FAST ? (adam_steps = 30, lbfgs_iters = 30) : (adam_steps = 300, lbfgs_iters = 300)
const NSEEDS = FAST ? 2 : parse(Int, get(ENV, "NN_SEEDS", "5"))
t_start = time()
progress(msg) = (println("[$(round((time() - t_start) / 60; digits = 1)) min] ", msg); flush(stdout))

trajs = N.load_all()
out = N.outdir("ude")
ckdir = mkpath(joinpath(out, "checkpoints"))
function cached(f, name)
    path = joinpath(ckdir, name * ".jls")
    isfile(path) && return deserialize(path)
    r = f(); serialize(path, r)
    return r
end
drug = N.make_units(N.training(trajs; role = :drug))
comp = N.make_units(N.training(trajs; role = :competition))
folds = [
    "hold out on→off intermittent" => u -> u.members[1].source == :intermittent && u.sched.doses[1] > 0,
    "hold out off→on intermittent" => u -> u.members[1].source == :intermittent && u.sched.doses[1] == 0,
    "hold out 1 uM monoculture" => u -> u.members[1].source == :mono && u.sched.doses[1] == 1.0,
    "hold out all 20k wells" => u -> N.is20k(u),
]
const KINDS = Dict("UDE (NN kill)" => :nn, "pure neural ODE" => :blackbox)

p1 = [(fi, kind) for fi in eachindex(folds) for kind in (:none, :hill)]
progress("fold phase 1: $(length(p1)) mechanistic fits")
r1 = N.pmap_dynamic(p1) do (fi, kind)
    cached("p1_$(fi)_$(kind)_0") do
        N.fit_drug_model(kind, filter(!folds[fi].second, drug); ITERS...)
    end
end
get1(fi, kind) = r1[findfirst(==((fi, kind)), p1)]

# neural: pure neural ODE first (longest), then the UDE
p2 = [(fi, m, s) for m in ("pure neural ODE", "UDE (NN kill)") for fi in eachindex(folds) for s in 1:NSEEDS]
progress("fold phase 2: $(length(p2)) neural fits ($(count(t -> isfile(joinpath(ckdir, "p2_$(t[1])_$(KINDS[t[2]])_$(t[3]).jls")), p2)) already done)")
r2 = N.pmap_dynamic(p2) do (fi, mname, seed)
    kind = KINDS[mname]; held = folds[fi].second
    tr = filter(!held, drug)
    tr = kind == :blackbox ? vcat(tr, filter(!held, comp)) : tr
    f = cached("p2_$(fi)_$(kind)_$(seed)") do
        θ0 = N.init_theta(kind; rng = N.StableRNG(1000 + seed)); θ0[1:N.NB] .= get1(fi, :hill).θ[1:N.NB]
        N.fit_drug_model(kind, tr; θinit = θ0, ITERS...)
    end
    progress("done fold $(fi) $(mname) seed $(seed)")
    f
end

fold_rows = DataFrame(fold = String[], model = String[], seed = Int[], rmse = Float64[], n = Int[])
for (fi, (name, held)) in enumerate(folds)
    te = filter(held, drug)
    score(kind, f) = (r = reduce(vcat, [N.forward_residuals(N.UDEModel(kind, f.θ, ones(4)), u; nominal = true) for u in te]); (N.rmse(r.obs .- r.pred), nrow(r)))
    push!(fold_rows, (name, "no drug", 1, score(:none, get1(fi, :none))...))
    push!(fold_rows, (name, "Hill (mechanistic)", 1, score(:hill, get1(fi, :hill))...))
    for (i, (fj, mname, s)) in enumerate(p2)
        fj == fi && push!(fold_rows, (name, mname, s, score(KINDS[mname], r2[i])...))
    end
end
sdz(x) = length(x) > 1 ? std(x) : 0.0
fold_tbl = combine(groupby(fold_rows, [:fold, :model]), :rmse => mean => :rmse, :rmse => sdz => :rmse_sd_over_seeds,
                   :seed => length => :seeds, :n => first => :n)
CSV.write(joinpath(out, "holdout_folds_per_seed.csv"), fold_rows)
CSV.write(joinpath(out, "holdout_folds.csv"), fold_tbl)
println(fold_tbl)

fmtpm(m, s) = s > 0 ? "$(round(m; digits = 3)) ± $(round(s; digits = 3))" : "$(round(m; digits = 3))"
fw = unstack(transform(fold_tbl, [:rmse, :rmse_sd_over_seeds] => ByRow(fmtpm) => :v)[:, [:model, :fold, :v]], :model, :fold, :v)
runtime = round((time() - t_start) / 60; digits = 1)
N.write_report(joinpath(out, "folds.html"), "Stage 2b — hold-out folds",
    [("What this is", N.para("Each fold refits the models without some training drug wells and predicts those wells from their nominal seeding (no refit). The intermittent folds test an unseen dose schedule — what optimal control relies on; the 1 µM fold an unseen dose between tested doses; the 20k fold an unseen seeding density. The UDE and the pure neural ODE are refitted from $(NSEEDS) seeds each (mean ± SD over seeds); Hill and no drug are deterministic single fits. The learned-resource UDE is not refitted here (runtime). Runtime: $(runtime) min (fits reused from checkpoints where present). <a href=\"report.html\">Back to Stage 2</a>.")),
     ("Held-out error by fold", N.html_table(fw; tip = ("Log-scale RMSE of the held-out wells of each fold (0.2 ≈ 22% typical error); neural models mean ± SD over seeds.", :lower)))])
println("wrote ", joinpath(out, "folds.html"), " in $(runtime) min")
