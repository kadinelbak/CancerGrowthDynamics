# Stage 2: UDE (logistic backbone + NN drug-kill), with Hill / no-drug / pure-neural-ODE comparators.
# Run with threads, e.g. `julia -t 12 --project=. scripts/run_03_ude.jl`.
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Plots, Statistics
const N = NNTrack
const FAST = get(ENV, "NN_FAST", "0") == "1"      # quick smoke run
const ITERS = FAST ? (adam_steps = 60, lbfgs_iters = 60) : (adam_steps = 300, lbfgs_iters = 300)
const NENS = FAST ? 2 : 8

trajs = N.load_all()
apc = N.area_per_count(trajs); N._APC[] = apc
out = N.outdir("ude")
drug = N.make_units(N.training(trajs; role = :drug), apc)
comp = N.make_units(N.training(trajs; role = :competition), apc)
val = N.make_units(N.validation(trajs), apc)
unit_id(u) = u.members[1].group

# ---------- Step A: fit on drug-role units ----------
@info "fitting no-drug and Hill ensembles"
none_fit = N.fit_drug_model(:none, drug; ITERS...)
hill_ens = N.fit_ensemble(:hill, drug; n = NENS, θwarm = N.init_theta(:hill), ITERS...)
@info "fitting NN (UDE) ensemble"
nn_ens = N.fit_ensemble(:nn, drug; n = NENS, θwarm = hill_ens[1].θ, ITERS...)
@info "fitting pure neural ODE"
bb_fit = N.fit_drug_model(:blackbox, drug; θinit = N.init_theta(:blackbox), ITERS...)

# ---------- Step B: competition on untreated co-culture ----------
mk(kind, f) = N.calibrate(kind, f.θ, comp)
models = Dict(
    "no drug" => [mk(:none, none_fit)],
    "Hill (mechanistic)" => [mk(:hill, f) for f in hill_ens],
    "UDE (NN kill)" => [mk(:nn, f) for f in nn_ens],
    "pure neural ODE" => [N.UDEModel(:blackbox, bb_fit.θ, ones(4))],
)
N.save_models(joinpath(out, "models.json"), models)
@info "competition" Dict(k => v[1].a for (k, v) in models)

# ---------- hold-out folds (control-relevant generalisation) ----------
folds = [
    "hold out on→off intermittent" => u -> u.members[1].source == :intermittent && u.sched.doses[1] > 0,
    "hold out off→on intermittent" => u -> u.members[1].source == :intermittent && u.sched.doses[1] == 0,
    "hold out 1 uM monoculture" => u -> u.members[1].source == :mono && u.sched.doses[1] == 1.0,
    "hold out all 20k wells" => u -> N.is20k(u),
]
fold_rows = Vector{DataFrame}(undef, length(folds))
Threads.@threads for fi in eachindex(folds)
    name, held = folds[fi]
    tr = filter(!held, drug); te = filter(held, drug)
    fn = N.fit_drug_model(:none, tr; ITERS...)
    fh = N.fit_drug_model(:hill, tr; ITERS...)
    fnn = N.fit_drug_model(:nn, tr; θinit = vcat(fh.θ[1:N.NB], N.init_theta(:nn)[N.NB + 1:end]), ITERS...)
    fbb = N.fit_drug_model(:blackbox, tr; ITERS...)
    rows = DataFrame(fold = String[], model = String[], rmse = Float64[], n = Int[])
    for (mname, kind, f) in [("no drug", :none, fn), ("Hill (mechanistic)", :hill, fh), ("UDE (NN kill)", :nn, fnn), ("pure neural ODE", :blackbox, fbb)]
        m = N.UDEModel(kind, f.θ, ones(4))
        r = reduce(vcat, [N.forward_residuals(m, u) for u in te])
        push!(rows, (name, mname, N.rmse(r.obs .- r.pred), nrow(r)))
    end
    fold_rows[fi] = rows
end
fold_tbl = reduce(vcat, fold_rows)
CSV.write(joinpath(out, "holdout_folds.csv"), fold_tbl)

# ---------- metrics ----------
sets = ["train: drug units" => drug, "train: untreated co-culture" => comp,
        "VALIDATION: ce1 treated co-culture" => filter(u -> u.members[1].split == :validation_treated_co, val)]
metrics = DataFrame(model = String[], set = String[], forward_rmse = Float64[], nextday_rmse = Float64[], n = Int[])
preds = DataFrame()
for (mname, ms) in models, (sname, us) in sets
    r = reduce(vcat, [N.ensemble_forward(ms, u) for u in us])
    nd = reduce(vcat, [N.nextday_residuals(ms[1], u) for u in us])
    push!(metrics, (mname, sname, N.rmse(r.obs .- r.pred), N.rmse(nd), nrow(r)))
    r.model .= mname; r.set .= sname
    append!(preds, r)
end
CSV.write(joinpath(out, "metrics.csv"), metrics)
CSV.write(joinpath(out, "predictions.csv"), preds)

# backbone parameters with bootstrap intervals (natural scale), incl. the 20k/30k growth-rate ratio
bb = DataFrame(model = String[], parameter = String[], median = Float64[], lo_5 = Float64[], hi_95 = Float64[])
for mname in ("Hill (mechanistic)", "UDE (NN kill)")
    Θ = reduce(hcat, [m.θ[1:N.NB] for m in models[mname]])
    for (j, nm) in enumerate(N.BACKBONE_NAMES)
        v = exp.(Θ[j, :])
        push!(bb, (mname, replace(String(nm), "log_" => "", "log" => ""), median(v), quantile(v, 0.05), quantile(v, 0.95)))
    end
end
CSV.write(joinpath(out, "backbone_parameters.csv"), bb)

# BIC on the drug units (weighted SSE is a chi-square, so BIC = SSE + k log n)
npts = sum(sum(!isnan, y) for u in drug for y in u.y)
bic = DataFrame(model = String[], kill_params = Int[], weighted_sse = Float64[], bic = Float64[])
for (mname, kind, f) in [("no drug", :none, none_fit), ("Hill (mechanistic)", :hill, hill_ens[1]),
                         ("UDE (NN kill)", :nn, nn_ens[1]), ("pure neural ODE", :blackbox, bb_fit)]
    sse = sum(N.unit_sse(kind, f.θ, ones(4), u, l0) for (u, l0) in zip(drug, f.l0s))
    k = N.NB + N.nkill(kind)
    push!(bic, (mname, N.nkill(kind), sse, sse + k * log(npts)))
end
CSV.write(joinpath(out, "bic.csv"), bic)
println(metrics); println(fold_tbl); println(bic)

# ---------- figures ----------
gr()
colors = Dict("no drug" => :gray, "Hill (mechanistic)" => :darkorange, "UDE (NN kill)" => :royalblue, "pure neural ODE" => :seagreen)
"""
Grid of per-trajectory panels in cell counts (log axis) with the SAME y-range on every panel
and one large shared key across the top.
"""
function panel_plot(set_name, fname; ncol = 4)
    d = preds[preds.set .== set_name, :]
    ids = unique(d.traj_id)
    order = ["UDE (NN kill)", "Hill (mechanistic)", "pure neural ODE", "no drug"]
    allv = exp.(vcat(d.obs, d.lo, d.hi))
    ylims_ = (10.0^floor(log10(minimum(allv))), 10.0^ceil(log10(maximum(allv))))
    ps = []
    for id in ids
        p = plot(; title = id, titlefontsize = 9, legend = false, tickfontsize = 7, yscale = :log10,
                 ylims = ylims_, xlims = (0, 14.5), xlabel = "day", ylabel = "cells", guidefontsize = 7)
        for mname in order
            s = d[(d.traj_id .== id) .& (d.model .== mname), :]
            plot!(p, s.day, exp.(s.pred); ribbon = (exp.(s.pred) .- exp.(s.lo), exp.(s.hi) .- exp.(s.pred)),
                  color = colors[mname], fillalpha = 0.2, lw = 2)
        end
        s = d[(d.traj_id .== id) .& (d.model .== "UDE (NN kill)"), :]
        scatter!(p, s.day, exp.(s.obs); color = :black, ms = 3.5)
        push!(ps, p)
    end
    key = plot(; framestyle = :none, legend = :top, legend_column = length(order) + 1, legendfontsize = 12,
               foreground_color_legend = nothing, xlims = (0, 1), ylims = (0, 1))
    for mname in order
        plot!(key, [-1, -1], [-1, -1]; color = colors[mname], lw = 5, label = mname)
    end
    scatter!(key, [-1], [-1]; color = :black, ms = 7, label = "observed")
    nrow_ = cld(length(ps), ncol)
    grid_ = plot(ps...; layout = (nrow_, ncol))
    fig = plot(key, grid_; layout = @layout([a{0.06h}; b]), size = (1400, 90 + 300 * nrow_), left_margin = 5Plots.mm, bottom_margin = 3Plots.mm)
    savefig(fig, joinpath(out, fname))
    return fname
end
f_train = panel_plot("train: drug units", "fit_drug_units.png")
f_val2 = panel_plot("VALIDATION: ce1 treated co-culture", "validation_treated_co.png")
f_comp = panel_plot("train: untreated co-culture", "fit_untreated_coculture.png"; ncol = 4)

vgrid = range(0, 1.6; length = 60)
pk = plot(; xlabel = "exposure v (uM)", ylabel = "kill rate (1/day)", title = "Learned drug-kill term (ρ = 0.3)", legend = :topleft)
for (line, ls, lab) in [(-1.0, :solid, "A2780 (S)"), (1.0, :dash, "A2780cis (R)")]
    for (i, m) in enumerate(models["UDE (NN kill)"])
        plot!(pk, vgrid, [N.kill_rate(:nn, m.θ, v, 0.3, line) for v in vgrid]; color = :royalblue, alpha = 0.35, ls, label = i == 1 ? "UDE members, $lab" : "")
    end
    m = models["Hill (mechanistic)"][1]
    plot!(pk, vgrid, [N.kill_rate(:hill, m.θ, v, 0.3, line) for v in vgrid]; color = :darkorange, lw = 2.5, ls, label = "Hill, $lab")
end
vline!(pk, [0.62, 1.0, 1.47]; color = :gray, alpha = 0.4, label = "doses in training")
savefig(pk, joinpath(out, "kill_curves.png"))

wide = unstack(metrics, :model, :set, :forward_rmse)
N.write_report(joinpath(out, "report.html"), "Stage 2 — universal differential equation",
    [("Model", N.para("dS/dt = r<sub>S</sub>S(1 − (S + a<sub>SR</sub>R)/K<sub>S</sub>) − δ<sub>S</sub>(t)S − kill<sub>S</sub>(v, ρ)S, same for R; dv/dt = γ(u(t) − v). δ(t) is a delayed death active only in low-resource wells. UDE: kill = v·softplus(NN(v, ρ, line)), a 3-8-8-1 tanh net, so no drug ⇒ no kill. Comparators: Hill kill (mechanistic), no drug, and a pure neural ODE (growth rate = NN(ρ, v, line, low, t)). Competition a is fit afterwards on untreated co-culture only. Ensembles are $(NENS)-member Bayesian bootstraps over wells.")),
     ("Headline: forward prediction RMSE (log scale)", N.para("Every trajectory is predicted forward from its observed day-1 value with no refitting. The validation rows are data never used in fitting.") * N.html_table(wide)),
     ("Next-day RMSE (compare Stage 1)", N.html_table(unstack(metrics, :model, :set, :nextday_rmse))),
     ("Hold-out folds", N.para("Refit without the held-out trajectories and predict them; the intermittent folds test an unseen dose schedule, which is what optimal control relies on.") * N.html_table(unstack(fold_tbl, :model, :fold, :rmse))),
     ("BIC on training drug units", N.html_table(bic)),
     ("Backbone parameters (bootstrap median, 5–95%)", N.para("<code>r20_over_r30</code> is the growth-rate multiplier for 20k-seeded wells: an interval including 1 means density acts only through the starting cell number; the <code>hold out all 20k wells</code> fold above tests whether 30k-only fits predict 20k wells.") * N.html_table(bb)),
     ("Learned kill term", N.figure("kill_curves.png", "NN kill rate vs exposure, all ensemble members, against the Hill fit")),
     ("Validation predictions", N.figure(f_val2, "ce1 treated co-culture (held out)")),
     ("Training fits", N.figure(f_train, "Drug-role units") * N.figure(f_comp, "Untreated co-culture (competition fit)"))])
println("wrote ", joinpath(out, "report.html"))
