# Stage 4: feature importance.
#  A. Next-day models (Stage 1): permutation importance (train in-sample and held-out validation)
#     and TreeSHAP for boosted trees.
#  B. UDE kill term: gradient sensitivity at visited states, partial dependence, and
#     input-ablation retraining (the most honest measure with this little data).
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Plots, Statistics, StableRNGs, EvoTrees, ForwardDiff, Random
const N = NNTrack
const FAST = get(ENV, "NN_FAST", "0") == "1"
const ITERS = FAST ? (adam_steps = 40, lbfgs_iters = 40) : (adam_steps = 300, lbfgs_iters = 300)
# ablation retrains: 3 seeds per variant (runtime, 2026-10-07: the full rerun shared the machine with the
# mechanistic refit); the full model's 5-seed fits are in Stage 2
const NSEEDS = FAST ? 2 : parse(Int, get(ENV, "NN_ABL_SEEDS", "3"))
t_start = time()

trajs = N.load_all()
out = N.outdir("importance")
gr()

# ================= A. next-day models =================
X, y, meta = N.nextday_rows(trajs)
tr = meta.split .== :train; va = .!tr
trees = N.fit_trees(X[:, tr], y[tr])
mlp = N.fit_mlp_ensemble(X[:, tr], y[tr]; n = 5, hidden = 12, lambda = 3e-3, steps = 2500)

"Grouped-by-row permutation importance: mean ΔRMSE over `reps` shuffles."
function perm_importance(model, X, y; reps = 20, rng = StableRNG(7))
    base = N.rmse(model(X) .- y)
    out = zeros(size(X, 1))
    for j in 1:size(X, 1), _ in 1:reps
        Xp = copy(X); Xp[j, :] = Xp[j, randperm(rng, size(X, 2))]
        out[j] += (N.rmse(model(Xp) .- y) - base) / reps
    end
    return out
end
imp = DataFrame(feature = String.(N.FEATURES))
imp.trees_train = perm_importance(trees, X[:, tr], y[tr])
imp.trees_validation = perm_importance(trees, X[:, va], y[va])
imp.mlp_train = perm_importance(mlp, X[:, tr], y[tr])
imp.mlp_validation = perm_importance(mlp, X[:, va], y[va])
S = EvoTrees.Shap.shap(trees.m, N._tree_table(X[:, tr]))
imp.trees_mean_abs_shap = vec(mean(abs.(S); dims = 1))
sort!(imp, :trees_mean_abs_shap; rev = true)
CSV.write(joinpath(out, "nextday_importance.csv"), imp)
println(imp)

p1 = bar(imp.feature, imp.trees_mean_abs_shap; orientation = :h, yflip = true, legend = false,
         title = "Boosted trees: mean |SHAP| (Δlog N per day)", xlabel = "mean |SHAP|", titlefontsize = 9)
sp = []
for f in (:logN, :exposure, :dose, :dprev)
    j = findfirst(==(f), N.FEATURES)
    push!(sp, scatter(X[j, tr], S[:, j]; ms = 2, alpha = 0.5, marker_z = X[findfirst(==(:line_R), N.FEATURES), tr],
                      color = cgrad([:crimson, :royalblue]), colorbar = false, legend = false,
                      xlabel = string(f), ylabel = "SHAP", title = "SHAP dependence: $(f) (red S, blue R)", titlefontsize = 8))
end
savefig(plot(p1, sp...; layout = @layout([a{0.35w} grid(2, 2)]), size = (1400, 600), margin = 4Plots.mm), joinpath(out, "nextday_shap.png"))

# ================= B. UDE kill term =================
models = N.load_models(joinpath(N.OUT_ROOT, "ude", "models.json"))
ude = models["UDE (NN kill)"]
drug = N.make_units(N.training(trajs; role = :drug))
comp = N.make_units(N.training(trajs; role = :competition))
val = N.make_units(N.validation(trajs))

# gradient sensitivity at visited states: E|∂kill/∂x| · sd(x), and the line effect
sens = DataFrame(member = Int[], input = String[], sensitivity = Float64[])
for (i, m) in enumerate(ude)
    pts = N.visited_states(m, drug)
    allp = vcat([(v, ρ, l) for l in (-1.0, 1.0) for (v, ρ) in pts[l]]...)
    V = [p[1] for p in allp]; P = [p[2] for p in allp]
    gv = [ForwardDiff.derivative(x -> N.kill_rate(:nn, m.θ, x, p[2], p[3]), p[1]) for p in allp]
    gρ = [ForwardDiff.derivative(x -> N.kill_rate(:nn, m.θ, p[1], x, p[3]), p[2]) for p in allp]
    line_eff = [abs(N.kill_rate(:nn, m.θ, p[1], p[2], 1.0) - N.kill_rate(:nn, m.θ, p[1], p[2], -1.0)) for p in allp]
    push!(sens, (i, "exposure v", mean(abs.(gv)) * std(V)))
    push!(sens, (i, "density ρ", mean(abs.(gρ)) * std(P)))
    push!(sens, (i, "line (S vs R)", mean(line_eff)))
end
sens_summary = combine(groupby(sens, :input), :sensitivity => mean => :mean, :sensitivity => std => :sd)
CSV.write(joinpath(out, "ude_sensitivity.csv"), sens_summary)

# partial dependence
pd = plot(; layout = (1, 2), size = (1200, 420), margin = 5Plots.mm)
vg = range(0, 1.5; length = 50)
for (line, c, lab) in [(-1.0, :crimson, "A2780 (S)"), (1.0, :royalblue, "A2780cis (R)")]
    ρs = [p[2] for p in N.visited_states(ude[1], drug)[line]]
    vs = [p[1] for p in N.visited_states(ude[1], drug)[line]]
    pdv = [mean(mean(N.kill_rate(:nn, m.θ, v, ρ, line) for ρ in ρs) for m in ude) for v in vg]
    plot!(pd[1], vg, pdv; color = c, lw = 2, label = lab, xlabel = "exposure v (uM)", ylabel = "kill (1/day)", title = "Partial dependence on v")
    ρg = range(extrema(ρs)...; length = 40)
    pdr = [mean(mean(N.kill_rate(:nn, m.θ, v, ρ, line) for v in vs) for m in ude) for ρ in ρg]
    plot!(pd[2], ρg, pdr; color = c, lw = 2, label = lab, xlabel = "relative density ρ", ylabel = "kill (1/day)", title = "Partial dependence on ρ")
end
savefig(pd, joinpath(out, "ude_partial_dependence.png"))

# input-ablation retrain, NSEEDS seeds per variant
variants = ["full" => (1.0, 1.0, 1.0), "without density ρ" => (1.0, 0.0, 1.0), "without line" => (1.0, 1.0, 0.0)]
abl_seed = DataFrame(variant = String[], seed = Int[], weighted_sse = Float64[], val_treated_co_rmse = Float64[])
# KILL_MASK is global, so variants run one after another; the seeds of a variant run in parallel
for (name, mask) in variants
    N.KILL_MASK[] = mask
    fs = N.fit_seeds(:nn, drug; n = NSEEDS, θwarm = ude[1].θ, ITERS...)
    for (i, f) in enumerate(fs)
        m = N.calibrate(:nn, f.θ, comp)
        sse = sum(N.unit_sse(:nn, f.θ, ones(4), u, l0) for (u, l0) in zip(drug, f.l0s))
        r2 = reduce(vcat, [N.forward_residuals(m, u; nominal = true) for u in val])
        push!(abl_seed, (name, i, sse, N.rmse(r2.obs .- r2.pred)))
    end
    @info "ablation" name
end
N.KILL_MASK[] = (1.0, 1.0, 1.0)
sdz(x) = length(x) > 1 ? std(x) : 0.0
abl = combine(groupby(abl_seed, :variant; sort = false), :weighted_sse => mean => :weighted_sse, :weighted_sse => sdz => :weighted_sse_sd,
              :weighted_sse => minimum => :weighted_sse_best_seed, :val_treated_co_rmse => mean => :val_treated_co_rmse,
              :val_treated_co_rmse => sdz => :val_treated_co_rmse_sd, :seed => length => :seeds)
abl.delta_sse_best_seed = abl.weighted_sse_best_seed .- abl.weighted_sse_best_seed[1]
CSV.write(joinpath(out, "ude_ablation_per_seed.csv"), abl_seed)
CSV.write(joinpath(out, "ude_ablation.csv"), abl)
println(sens_summary); println(abl)

runtime = round((time() - t_start) / 60; digits = 1)
N.write_report(joinpath(out, "report.html"), "Stage 4 — feature importance",
    [("Next-day models", N.para("Permutation importance is the rise in one-step RMSE (log scale) when a feature is shuffled: on training rows (in-sample) and on the held-out validation rows. TreeSHAP is exact for the boosted trees. No clock-time features are used (see Stage 1). Runtime of this stage: $(runtime) min.") *
        N.html_table(imp; tip = ("Rise in RMSE when a feature is shuffled (permutation) and mean |SHAP| of the boosted trees. A larger value = the model relies more on that feature; near 0 or negative = not used.", "Neither: larger means more important, not better.")) *
        N.figure("nextday_shap.png", "TreeSHAP ranking and dependence plots"; tip = ("Mean |SHAP| per feature and how each feature's value moves the prediction.", :none))),
     ("UDE kill term: gradient sensitivity", N.para("Mean |∂kill/∂input| × sd(input) over the states the fitted model visits; for line, the mean |kill<sub>R</sub> − kill<sub>S</sub>| at the same state. Mean and sd over the $(length(ude)) seeds.") *
        N.html_table(sens_summary; tip = ("Sensitivity of the learned kill rate to each input, mean and SD over seeds.", "Neither: larger means the kill depends more on that input.")) ),
     ("UDE kill term: partial dependence", N.figure("ude_partial_dependence.png", "Kill rate averaged over the other input (mean over seeds)"; tip = ("Kill rate as one input varies, averaged over the visited values of the other.", :none))),
     ("UDE kill term: input ablation (retrained)", N.para("Each variant removes one NN input and is retrained from $(NSEEDS) seeds. Δ weighted SSE (best seed) on training (a χ², so ~2 per unit of log-likelihood) is the in-sample cost; validation RMSE (mean ± SD over seeds, from the nominal seeding) is the out-of-sample cost.") *
        N.html_table(abl; tip = ("Training fit and held-out error when one input is removed from the kill network. If removing an input does not raise the validation RMSE, the data do not show it is needed.", :lower)))])
println("wrote ", joinpath(out, "report.html"))
