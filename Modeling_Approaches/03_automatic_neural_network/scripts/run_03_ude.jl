# Stage 2: UDE (logistic backbone + NN drug-kill), with Hill / no-drug / pure-neural-ODE comparators.
# Run with threads, e.g. `julia -t 6 --project=. scripts/run_03_ude.jl`.
# Every neural model is fitted from NN_SEEDS (default 5) random initialisations; scores are reported as
# mean ± SD over seeds and for the seed-ensemble mean. Held-out wells start from the nominal seeding.
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Plots, Statistics, Serialization
const N = NNTrack
const FAST = get(ENV, "NN_FAST", "0") == "1"      # quick smoke run
const ITERS = FAST ? (adam_steps = 30, lbfgs_iters = 30) : (adam_steps = 300, lbfgs_iters = 300)
const NSEEDS = FAST ? 2 : parse(Int, get(ENV, "NN_SEEDS", "5"))
const NBOOT = FAST ? 2 : 8
const RESUME = get(ENV, "NN_RESUME", "0") == "1"  # reuse finished fits from outputs/ude/checkpoints
t_start = time()
progress(msg) = (println("[$(round((time() - t_start) / 60; digits = 1)) min] ", msg); flush(stdout))

trajs = N.load_all()
out = N.outdir("ude")
ckdir = mkpath(joinpath(out, "checkpoints"))
RESUME || foreach(f -> rm(joinpath(ckdir, f)), readdir(ckdir))
"Run `f()` or, when resuming, load its result from the checkpoint `name`."
function cached(f, name)
    path = joinpath(ckdir, name * ".jls")
    RESUME && isfile(path) && return deserialize(path)
    r = f(); serialize(path, r)
    return r
end
drug = N.make_units(N.training(trajs; role = :drug))
comp = N.make_units(N.training(trajs; role = :competition))
val = N.make_units(N.validation(trajs))
const NEURAL_KINDS = Dict("UDE (NN kill)" => :nn, "UDE (NN kill + NN resource)" => :nn_res, "pure neural ODE" => :blackbox)

# ---------- phase 1 (mechanistic): no drug + Hill (Bayesian bootstrap; member 1 = unit weights) ----------
# The hold-out folds are fitted by run_03b_folds.jl (run after Stage 7, so Stages 3-7 do not wait for them).
p1 = vcat([(:main, :none, 0)], [(:main, :hill, b) for b in 1:NBOOT])
progress("phase 1: $(length(p1)) mechanistic fits")
r1 = N.pmap_dynamic(p1) do (where_, kind, b)
    cached("p1_$(where_)_$(kind)_$(b)") do
        if kind == :hill
            w = b == 1 ? ones(length(drug)) : (rng = N.StableRNG(1000 + b); e = -log.(rand(rng, length(drug))); e ./ sum(e) .* length(drug))
            N.fit_drug_model(:hill, drug; weights = w, ITERS...)
        else
            N.fit_drug_model(kind, drug; ITERS...)
        end
    end
end
none_fit = r1[1]; hill_ens = r1[2:end]

# ---------- phase 2 (neural): every seed of every neural model, one dynamic queue ----------
p2 = [(m, s) for m in ("UDE (NN kill + NN resource)", "pure neural ODE", "UDE (NN kill)") for s in 1:NSEEDS]   # longest first
progress("phase 2: $(length(p2)) neural fits")
r2 = N.pmap_dynamic(p2) do (mname, seed)
    kind = NEURAL_KINDS[mname]
    tr = kind == :blackbox ? vcat(drug, comp) : drug      # the pure neural ODE also sees the untreated co-culture
    f = cached("p2_main_$(kind)_$(seed)") do
        θ0 = N.init_theta(kind; rng = N.StableRNG(1000 + seed)); θ0[1:N.NB] .= hill_ens[1].θ[1:N.NB]
        N.fit_drug_model(kind, tr; θinit = θ0, ITERS...)
    end
    progress("done main $(mname) seed $(seed)")
    f
end
get2(mname) = [r2[i] for i in eachindex(p2) if p2[i][1] == mname]
nn_fits = get2("UDE (NN kill)"); nnres_fits = get2("UDE (NN kill + NN resource)"); bb_fits = get2("pure neural ODE")
progress("all fits done")

# ---------- Step B: competition on untreated co-culture ----------
mk(kind, f) = N.calibrate(kind, f.θ, comp)
ctasks = vcat([(:none, none_fit)], [(:hill, f) for f in hill_ens], [(:nn, f) for f in nn_fits], [(:nn_res, f) for f in nnres_fits])
calib = N.UDEModel[m for m in N.pmap_dynamic(t -> mk(t...), ctasks)]
models = Dict(
    "no drug" => [calib[1]],
    "Hill (mechanistic)" => calib[2:1 + NBOOT],
    "UDE (NN kill)" => calib[2 + NBOOT:1 + NBOOT + NSEEDS],
    "UDE (NN kill + NN resource)" => calib[2 + NBOOT + NSEEDS:1 + NBOOT + 2 * NSEEDS],
    "pure neural ODE" => [N.UDEModel(:blackbox, f.θ, ones(4)) for f in bb_fits],
)
N.save_models(joinpath(out, "models.json"), models)
progress("competition calibrated")
# the models scored as single fits: every seed of a neural model; the unit-weight fit of Hill / no drug
scored(mname) = mname == "Hill (mechanistic)" ? models[mname][1:1] : models[mname]

sdz(x) = length(x) > 1 ? std(x) : 0.0

# ---------- metrics ----------
sets = ["train: drug units" => (drug, false), "train: untreated co-culture" => (comp, false),
        "VALIDATION: ce1 treated co-culture" => (val, true)]
metrics = DataFrame(model = String[], set = String[], forward_rmse = Float64[], forward_rmse_sd_over_seeds = Float64[],
                    forward_rmse_seed_ensemble = Float64[], nextday_rmse = Float64[], seeds = Int[], n = Int[])
preds = DataFrame()
for (mname, ms) in models, (sname, (us, nominal)) in sets
    sm = scored(mname)
    per = [(r = reduce(vcat, [N.forward_residuals(m, u; nominal) for u in us]); N.rmse(r.obs .- r.pred)) for m in sm]
    r = reduce(vcat, [N.ensemble_forward(sm, u; nominal) for u in us])
    nd = [N.rmse(reduce(vcat, [N.nextday_residuals(m, u; nominal) for u in us])) for m in sm]
    push!(metrics, (mname, sname, mean(per), sdz(per), N.rmse(r.obs .- r.pred), mean(nd), length(sm), nrow(r)))
    r.model .= mname; r.set .= sname
    append!(preds, r)
end
CSV.write(joinpath(out, "metrics.csv"), metrics)
CSV.write(joinpath(out, "predictions.csv"), preds)

# backbone parameters (natural scale): Hill over bootstrap members, neural models over seeds
bb = DataFrame(model = String[], over = String[], parameter = String[], median = Float64[], lo_5 = Float64[], hi_95 = Float64[])
for (mname, over) in (("Hill (mechanistic)", "bootstrap"), ("UDE (NN kill)", "seeds"), ("UDE (NN kill + NN resource)", "seeds"))
    Θ = reduce(hcat, [m.θ[1:N.NB] for m in models[mname]])
    for (j, nm) in enumerate(N.BACKBONE_NAMES)
        v = exp.(Θ[j, :])
        push!(bb, (mname, over, replace(String(nm), "log_" => "", "log" => ""), median(v), quantile(v, 0.05), quantile(v, 0.95)))
    end
end
CSV.write(joinpath(out, "backbone_parameters.csv"), bb)

# BIC on the drug units (weighted SSE is a chi-square, so BIC = SSE + k log n); k counts the estimated
# starting values too. Neural models: every seed; `bic` = the best seed (lowest training loss).
npts = sum(sum(!isnan, y) for u in drug for y in u.y)
nstart = sum(N.nfree(u) for u in drug)
bic_rows = DataFrame(model = String[], seed = Int[], kill_params = Int[], weighted_sse = Float64[], bic = Float64[])
for (mname, kind, fs) in [("no drug", :none, [none_fit]), ("Hill (mechanistic)", :hill, [hill_ens[1]]),
                          ("UDE (NN kill)", :nn, nn_fits), ("UDE (NN kill + NN resource)", :nn_res, nnres_fits),
                          ("pure neural ODE", :blackbox, bb_fits)]
    for (i, f) in enumerate(fs)
        sse = sum(N.unit_sse(kind, f.θ, ones(4), u, l0) for (u, l0) in zip(drug, f.l0s))
        k = N.NB + N.nkill(kind) + nstart
        push!(bic_rows, (mname, i, N.nkill(kind), sse, sse + k * log(npts)))
    end
end
bic = combine(groupby(bic_rows, :model), :kill_params => first => :kill_params, :weighted_sse => minimum => :weighted_sse,
              :bic => minimum => :bic, :bic => mean => :bic_mean_over_seeds, :bic => sdz => :bic_sd_over_seeds, :seed => length => :seeds)
CSV.write(joinpath(out, "bic_per_seed.csv"), bic_rows)
CSV.write(joinpath(out, "bic.csv"), bic)
println(metrics); println(bic)

# ---------- figures ----------
gr()
colors = Dict("no drug" => :gray, "Hill (mechanistic)" => :darkorange, "UDE (NN kill)" => :royalblue,
              "UDE (NN kill + NN resource)" => :purple, "pure neural ODE" => :seagreen)
"""
Grid of per-trajectory panels in cell counts (log axis) with the SAME y-range on every panel
and one large shared key across the top.
"""
function panel_plot(set_name, fname; ncol = 4)
    d = preds[preds.set .== set_name, :]
    ids = unique(d.traj_id)
    order = ["UDE (NN kill)", "UDE (NN kill + NN resource)", "Hill (mechanistic)", "pure neural ODE", "no drug"]
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
        plot!(pk, vgrid, [N.kill_rate(:nn, m.θ, v, 0.3, line) for v in vgrid]; color = :royalblue, alpha = 0.35, ls, label = i == 1 ? "UDE seeds, $lab" : "")
    end
    m = models["Hill (mechanistic)"][1]
    plot!(pk, vgrid, [N.kill_rate(:hill, m.θ, v, 0.3, line) for v in vgrid]; color = :darkorange, lw = 2.5, ls, label = "Hill, $lab")
end
vline!(pk, [0.62, 1.0, 1.47]; color = :gray, alpha = 0.4, label = "doses in training")
savefig(pk, joinpath(out, "kill_curves.png"))

T_RMSE = "Root-mean-square error of log cell counts (0.2 ≈ 22% typical error). Mean over seeds of each seed's single fit; ± SD over seeds; 'seed ensemble' = the error of the mean prediction of all seeds."
fmtpm(m, s) = s > 0 ? "$(round(m; digits = 3)) ± $(round(s; digits = 3))" : "$(round(m; digits = 3))"
wide = unstack(transform(metrics, [:forward_rmse, :forward_rmse_sd_over_seeds] => ByRow(fmtpm) => :v)[:, [:model, :set, :v]], :model, :set, :v)
wide_e = unstack(metrics, :model, :set, :forward_rmse_seed_ensemble)
runtime = round((time() - t_start) / 60; digits = 1)
N.write_report(joinpath(out, "report.html"), "Stage 2 — universal differential equation",
    [("Model", N.para("<span class='eqg-text' data-vars='S,r_S,a_SR,R,K_S,d_S,Q,kill_S,v,ρ'>dS/dt = r<sub>S</sub>S(1 − (S + a<sub>SR</sub>R)/K<sub>S</sub>) − d<sub>S</sub>(Q)S − kill<sub>S</sub>(v, ρ)S</span>, same for R; <span class='eqg-text' data-vars='v,γ,u'>dv/dt = γ(u − v)</span> (drug exposure driven by the dose u); <span class='eqg-text' data-vars='Q,c,S,R'>dQ/dt = −cQ(S + R)/3500</span> (resource left, low-resource wells only); <span class='eqg-text' data-vars='d_L,Q,δ_L'>d<sub>L</sub>(Q) = δ<sub>L</sub>(1 − Q)</span>. <b>Every model is autonomous</b> — no clock-time term or input anywhere; delays come from v and Q. UDE: <span class='eqg-text' data-vars='kill,v,NN,ρ,line'>kill = v·softplus(NN(v, ρ, line))</span>, a 3-8-8-1 tanh net, so no drug ⇒ no kill. Comparators: Hill kill (mechanistic), no drug, and a pure neural ODE (per-capita growth = NN(ρ, v, line, low, Q)), which has no competition term and is therefore fitted to the drug and untreated co-culture units together. The UDE and Hill fit the competition a afterwards on untreated co-culture only.") *
        N.para("<b>Starting values</b>: every well starts at t = 0 from its nominal seeding (100 cells per 30k well, 66.7 per 20k well, times the line's share); day 0 is not a measurement. A training series whose first count is more than 3<sup>d</sup>-fold off nominal has its starting value estimated ($(nstart) of the training drug series). Held-out wells always start from the nominal seeding.") *
        N.para("<b>Seeds</b>: each neural model is fitted from $(NSEEDS) random initialisations (seeds 1–$(NSEEDS)); Hill and no-drug are deterministic single fits (an $(NBOOT)-member Bayesian bootstrap of Hill gives parameter intervals). Runtime of this run: $(runtime) min$(RESUME ? " (fits loaded from checkpoints are not timed; the full set of fits took about 2.5 h on 6 threads alongside the mechanistic refit, 2026-10-07)" : "").")),
     ("Headline: forward prediction RMSE (log scale)", N.para("Every trajectory is predicted forward from t = 0 with no refitting. The validation rows are data never used in fitting.") *
        N.html_table(wide; tip = ("Forward-prediction error on each set. " * T_RMSE * " Training rows show fit, not generalisation; the VALIDATION column is the out-of-sample test.", :lower)) *
        N.para("Seed-ensemble mean prediction:") *
        N.html_table(wide_e; tip = ("Error of the mean (log-scale) prediction over all seeds of a neural model; for Hill and no drug the bootstrap-ensemble mean.", :lower))),
     ("Next-day RMSE (compare Stage 1)", N.html_table(unstack(metrics, :model, :set, :nextday_rmse); tip = ("One-day-ahead error: restart from each observed day and predict the next day (log scale), mean over seeds.", :lower))),
     ("Hold-out folds", N.para("Refits without held-out trajectories (unseen schedules, an unseen dose, an unseen seeding density) are on the <a href=\"folds.html\">hold-out folds page</a> (run_03b_folds.jl).")),
     ("BIC on training drug units", N.para("The pure neural ODE was also fitted to the untreated co-culture, so its drug-unit BIC is not strictly comparable.") *
        N.html_table(bic; tip = ("BIC = weighted SSE + k·log(n), with k = all free parameters including the estimated starting values. 'bic' is the best seed; mean and SD over seeds are shown too. Differences under 3 are not meaningful.", :lower))),
     ("Backbone parameters (median, 5–95%)", N.para("<code>r20_over_r30</code> is the growth-rate multiplier for 20k-seeded wells: an interval including 1 means density acts only through the starting cell number.") *
        N.html_table(bb; tip = ("Fitted backbone parameters on the natural scale; Hill over bootstrap members, neural models over seeds. A narrow interval = well determined.", :none))),
     ("Learned kill term", N.figure("kill_curves.png", "NN kill rate vs exposure, every seed, against the Hill fit"; tip = ("Kill rate per day vs drug exposure at relative density 0.3. Agreement between seeds shows the kill term is determined by the data; spread shows it is not.", "Closer agreement between seeds is better (better determined)."))),
     ("Validation predictions", N.figure(f_val2, "ce1 treated co-culture (held out), from the nominal seeding"; tip = ("Points: measured held-out counts. Lines: predictions (seed-ensemble mean; band = 10–90% over seeds or bootstrap members).", "Better: lines run through the points."))),
     ("Training fits", N.figure(f_train, "Drug-role units"; tip = ("In-sample fits to the training drug units.", "Better: lines run through the points (in-sample, so it shows flexibility, not generalisation).")) *
        N.figure(f_comp, "Untreated co-culture (competition fit)"; tip = ("In-sample fits to the untreated co-culture (competition step).", "Better: lines run through the points.")))])
println("wrote ", joinpath(out, "report.html"), " in $(runtime) min")
