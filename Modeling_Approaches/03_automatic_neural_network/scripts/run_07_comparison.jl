# Stage 6: fair three-way comparison on the SAME held-out points.
#   held out : the 30k treated cocultures Ce1 (3 ratios × 2 lines = 6 trajectories, days 1-14, 84 points) of
#              the mechanistic zoo's tidy table (02_mechanical_automatic_package/data/measurements.csv)
#   training : the zoo's own training data for that prediction — the 30k monocultures (untreated, 0.62, 1,
#              1.47 µM) for growth and drug response, the 30k untreated cocultures Ce0 for competition —
#              read from the same tidy table (same points, same QC flags)
#   starts   : every held-out series starts from its nominal seeding (as validation.toml); training series
#              whose first count is > 3^d-fold off nominal have an estimated start (as zoo_data.jl x0_info)
# Models: (a) the zoo's carried held-out prediction (per point, exported by export_mech_heldout.jl),
#         (b) this track's Hill kill and recovered closed form, (c) the UDE, (d) the pure neural ODE, plus a
#         no-drug reference; neural models from NN_SEEDS seeds. The same models fitted to the full NN-track
#         training set (Stage 2/3) are scored on the same points as a secondary comparison.
# Needs outputs/symbolic/symbolic_models.json and outputs/ude/models.json (Stages 2-3) and, for (a),
# outputs/comparison/mech_heldout_points.csv. Run: julia -t 6 --project=. scripts/run_07_comparison.jl
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Plots, Statistics, JSON3, Random, Printf
const N = NNTrack
const FAST = get(ENV, "NN_FAST", "0") == "1"
const ITERS = FAST ? (adam_steps = 30, lbfgs_iters = 30) : (adam_steps = 300, lbfgs_iters = 300)
const NSEEDS = FAST ? 2 : parse(Int, get(ENV, "NN_SEEDS", "5"))
const NBOOT = FAST ? 500 : 10_000
t_start = time()
out = N.outdir("comparison")

tidy = N.load_tidy()
drug = N.make_units(filter(t -> t.role == :drug, tidy))
comp = N.make_units(filter(t -> t.role == :competition, tidy))
ho = N.make_units(filter(t -> t.source == :ce1, tidy))

sym = JSON3.read(read(joinpath(N.OUT_ROOT, "symbolic", "symbolic_models.json"), String))["optimal"]
N.SYMBOLIC_TERMS[] = (S = Int.(sym.termsS), R = Int.(sym.termsR))
closed_formula = "kill_S = $(sym.formulaS == "" ? "0" : sym.formulaS); kill_R = $(sym.formulaR == "" ? "0" : sym.formulaR)"
coef0 = Float64.(sym.members[1].theta)[N.NB + 1:end]

# ---------------- (b)-(d): fits on the matched training data ----------------
@info "matched fits: no drug, Hill, closed form"
f_none = N.fit_drug_model(:none, drug; ITERS...)
f_hill = N.fit_drug_model(:hill, drug; ITERS...)
f_sym = N.fit_drug_model(:symbolic, drug; θinit = vcat(f_hill.θ[1:N.NB], coef0), ITERS...)
tasks = [(kind, seed) for kind in (:nn, :blackbox) for seed in 1:NSEEDS]
nnf = Vector{Any}(undef, length(tasks))
@info "matched fits: UDE and pure neural ODE, $(NSEEDS) seeds each"
Threads.@threads for i in eachindex(tasks)
    kind, seed = tasks[i]
    θ0 = N.init_theta(kind; rng = N.StableRNG(1000 + seed)); θ0[1:N.NB] .= f_hill.θ[1:N.NB]
    nnf[i] = N.fit_drug_model(kind, kind == :blackbox ? vcat(drug, comp) : drug; θinit = θ0, ITERS...)
end
matched = Dict(
    "no drug (reference)" => [N.calibrate(:none, f_none.θ, comp)],
    "Hill kill (NN track, mechanistic)" => [N.calibrate(:hill, f_hill.θ, comp)],
    "closed form (NN track, mechanistic)" => [N.calibrate(:symbolic, f_sym.θ, comp)],
    "UDE (NN kill)" => [N.calibrate(:nn, nnf[i].θ, comp) for i in eachindex(tasks) if tasks[i][1] == :nn],
    "pure neural ODE (autonomous)" => [N.UDEModel(:blackbox, nnf[i].θ, ones(4)) for i in eachindex(tasks) if tasks[i][1] == :blackbox],
)
N.save_models(joinpath(out, "matched_models.json"), matched)
@info "matched fits done" minutes = round((time() - t_start) / 60; digits = 1)

# the same model families fitted to the full NN-track training set (Stages 2-3), no refit
full = N.load_models(joinpath(N.OUT_ROOT, "ude", "models.json"))
fullset = Dict(
    "Hill kill (NN track, mechanistic)" => full["Hill (mechanistic)"][1:1],
    "closed form (NN track, mechanistic)" => [N.UDEModel(:symbolic, Float64.(sym.members[1].theta), Float64.(sym.members[1].a))],
    "UDE (NN kill)" => full["UDE (NN kill)"],
    "pure neural ODE (autonomous)" => full["pure neural ODE"],
)

# ---------------- predictions on the held-out points ----------------
# point table: one row per held-out point (tidy table, days ≥ 1)
pts = reduce(vcat, [DataFrame(traj_id = tr.id, day = tr.days, obs = tr.mean) for u in ho for tr in u.members])
key(id, d) = (id, d)
idx = Dict(key(r.traj_id, r.day) => i for (i, r) in enumerate(eachrow(pts)))
"Cell-scale predictions (points × seeds) of a model set on the held-out points, from the nominal seeding."
function predict_cells(ms)
    P = fill(NaN, nrow(pts), length(ms))
    for (j, m) in enumerate(ms), u in ho
        r = N.forward_residuals(m, u; nominal = true)
        for row in eachrow(r)
            P[idx[key(row.traj_id, row.day)], j] = exp(row.pred)
        end
    end
    @assert !any(isnan, P)
    return P
end
geomean_rows(P) = exp.(vec(mean(log.(P); dims = 2)))

const ORDER = ["mechanistic zoo (carried model)", "Hill kill (NN track, mechanistic)", "closed form (NN track, mechanistic)",
               "UDE (NN kill)", "pure neural ODE (autonomous)", "no drug (reference)"]
preds = Dict{String, Vector{Float64}}()       # point predictions (seed ensemble = geometric mean)
seedpreds = Dict{String, Matrix{Float64}}()
for (name, ms) in matched
    P = predict_cells(ms); seedpreds[name] = P; preds[name] = geomean_rows(P)
end
fullpreds = Dict(name => geomean_rows(predict_cells(ms)) for (name, ms) in fullset)
fullseed = Dict(name => predict_cells(ms) for (name, ms) in fullset)

# (a) the mechanistic zoo, per point
mech_file = joinpath(out, "mech_heldout_points.csv")
mech_src = isfile(joinpath(out, "mech_heldout_source.txt")) ? strip(read(joinpath(out, "mech_heldout_source.txt"), String)) : ""
have_mech = isfile(mech_file)
if have_mech
    mp = CSV.read(mech_file, DataFrame; stringtype = String, types = Dict(:ratio => String))
    mp = mp[mp.day .>= 1, :]
    mid(r) = "ce1_$(r.ratio)_30k_$(r.line == "A2780cis" ? "R" : "S")"
    v = fill(NaN, nrow(pts))
    for r in eachrow(mp)
        i = idx[key(mid(r), r.day)]
        @assert abs(r.observed - pts.obs[i]) <= 1e-6 * max(1.0, pts.obs[i]) "zoo point $(mid(r)) day $(r.day) does not match the tidy table"
        v[i] = r.predicted
    end
    @assert !any(isnan, v) "the zoo export does not cover every held-out point"
    preds[ORDER[1]] = v
else
    @warn "no per-point zoo export; run scripts/export_mech_heldout.jl. Only the zoo's summary RMSE is reported."
end
zoo_summary = let f = joinpath(N.REPO_ROOT, "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "zoo_fits", "heldout_ce1.csv")
    isfile(f) ? CSV.read(f, DataFrame) : nothing
end

# ---------------- metrics, bootstrap CIs, parameter counts ----------------
logc(x) = log(max(x, 1.0))
function metrics(o, p)
    r = sqrt(mean(abs2, p .- o))
    return (rmse = r, nrmse = r / mean(o), logrmse = sqrt(mean(abs2, logc.(p) .- logc.(o))))
end
o = pts.obs
names_present = [n for n in ORDER if haskey(preds, n)]
M = Dict(n => metrics(o, preds[n]) for n in names_present)
seedM = Dict(n => [metrics(o, seedpreds[n][:, j]) for j in 1:size(seedpreds[n], 2)] for n in keys(seedpreds))

# paired bootstrap over the 6 held-out trajectories: resample trajectories with replacement, recompute every
# model's metric on the pooled points, take the difference against the model that is best on the full data
trajs_ho = unique(pts.traj_id)
rows_of = Dict(t => findall(==(t), pts.traj_id) for t in trajs_ho)
rng = N.StableRNG(2026)
cands = [n for n in names_present if n != "no drug (reference)"]
best = Dict(m => cands[argmin([getfield(M[n], m) for n in cands])] for m in (:rmse, :nrmse, :logrmse))
D = Dict((n, m) => Float64[] for n in names_present for m in (:rmse, :logrmse))
for b in 1:NBOOT
    ii = reduce(vcat, [rows_of[t] for t in rand(rng, trajs_ho, length(trajs_ho))])
    mb = Dict(n => metrics(o[ii], preds[n][ii]) for n in names_present)
    for n in names_present, m in (:rmse, :logrmse)
        push!(D[(n, m)], getfield(mb[n], m) - getfield(mb[best[m]], m))
    end
end
ci(x) = (quantile(x, 0.025), quantile(x, 0.975))

"Free parameters that act on these predictions (normal medium, 30k): growth r, K per line, exposure γ, detection floor; the drug term; competition a_SR, a_RS."
function nparams_used(name)
    name == "pure neural ODE (autonomous)" && return 2 + N.nparams(N.BLACKBOX_NET)          # γ, floor + network
    kill = name == "Hill kill (NN track, mechanistic)" ? length(N.HILL_NAMES) :
           name == "closed form (NN track, mechanistic)" ? length(N.SYMBOLIC_TERMS[].S) + length(N.SYMBOLIC_TERMS[].R) :
           name == "UDE (NN kill)" ? N.nparams(N.KILL_NET) : 0
    return 4 + (name == "no drug (reference)" ? 0 : 1) + 1 + kill + 2                        # r,K ×2, γ (if drug), floor, kill, a ×2
end
mech_params = isfile(joinpath(out, "mech_heldout_params.csv")) ? CSV.read(joinpath(out, "mech_heldout_params.csv"), DataFrame) : nothing
nstart_train = sum(N.nfree(u) for u in drug) + sum(N.nfree(u) for u in comp)
fmt3(x) = @sprintf("%.3f", x)
fmtci(c) = "[$(@sprintf("%+.0f", c[1])), $(@sprintf("%+.0f", c[2]))]"
fmtcil(c) = "[$(@sprintf("%+.3f", c[1])), $(@sprintf("%+.3f", c[2]))]"
pm(v) = length(v) > 1 ? "$(@sprintf("%.0f", mean(v))) ± $(@sprintf("%.0f", std(v)))" : @sprintf("%.0f", v[1])
pml(v) = length(v) > 1 ? "$(@sprintf("%.3f", mean(v))) ± $(@sprintf("%.3f", std(v)))" : @sprintf("%.3f", v[1])

tbl = DataFrame(model = String[], kind = String[], parameters = String[], seeds = Int[], rmse_cells = Float64[],
                rmse_cells_per_seed = String[], nrmse = Float64[], log_rmse = Float64[], log_rmse_per_seed = String[],
                d_rmse_vs_best = Float64[], d_rmse_ci95 = String[], d_logrmse_vs_best = Float64[], d_logrmse_ci95 = String[])
for n in names_present
    m = M[n]
    kind = n == ORDER[1] ? "mechanistic (zoo)" : occursin("mechanistic", n) ? "mechanistic (NN track)" : n == "no drug (reference)" ? "reference" : "neural"
    np_ = n == ORDER[1] ? (mech_params === nothing ? "?" : string(sum(mech_params.params))) : string(nparams_used(n))
    sm = get(seedM, n, nothing)
    push!(tbl, (n, kind, np_, sm === nothing ? 1 : length(sm), m.rmse,
                sm === nothing ? "—" : pm([s.rmse for s in sm]), m.nrmse, m.logrmse, sm === nothing ? "—" : pml([s.logrmse for s in sm]),
                m.rmse - M[best[:rmse]].rmse, n == best[:rmse] ? "best" : fmtci(ci(D[(n, :rmse)])),
                m.logrmse - M[best[:logrmse]].logrmse, n == best[:logrmse] ? "best" : fmtcil(ci(D[(n, :logrmse)]))))
end
CSV.write(joinpath(out, "comparison.csv"), tbl)

# per-point predictions (for re-analysis) and per-trajectory RMSE
pp = copy(pts)
for n in names_present
    pp[!, n] = preds[n]
end
CSV.write(joinpath(out, "heldout_point_predictions.csv"), pp)
pertraj = DataFrame(traj_id = trajs_ho, mean_obs = [mean(o[rows_of[t]]) for t in trajs_ho])
for n in names_present
    pertraj[!, n] = [sqrt(mean(abs2, preds[n][rows_of[t]] .- o[rows_of[t]])) for t in trajs_ho]
end
CSV.write(joinpath(out, "per_trajectory_rmse.csv"), pertraj)

# secondary: the same families fitted to the full NN-track training set
ftbl = DataFrame(model = String[], training = String[], rmse_cells = Float64[], rmse_cells_per_seed = String[], nrmse = Float64[], log_rmse = Float64[])
for n in ORDER
    haskey(fullpreds, n) || continue
    m = metrics(o, fullpreds[n]); sm = [metrics(o, fullseed[n][:, j]) for j in 1:size(fullseed[n], 2)]
    push!(ftbl, (n, "full NN-track set (20k + 30k, intermittent, low resource)", m.rmse, pm([s.rmse for s in sm]), m.nrmse, m.logrmse))
    haskey(M, n) && push!(ftbl, (n, "matched (zoo's 30k monocultures + Ce0)", M[n].rmse, haskey(seedM, n) ? pm([s.rmse for s in seedM[n]]) : "—", M[n].nrmse, M[n].logrmse))
end
CSV.write(joinpath(out, "training_set_effect.csv"), ftbl)
println(tbl); println(pertraj); println(ftbl)

# ---------------- figure ----------------
gr()
cols = Dict(ORDER[1] => :black, ORDER[2] => :darkorange, ORDER[3] => :crimson, ORDER[4] => :royalblue, ORDER[5] => :seagreen, ORDER[6] => :gray)
ps = []
for t in trajs_ho
    ii = rows_of[t]
    p = plot(; title = t, titlefontsize = 9, legend = false, xlabel = "day", ylabel = "cells", xlims = (0, 14.5), tickfontsize = 7)
    for n in names_present
        plot!(p, pts.day[ii], preds[n][ii]; color = cols[n], lw = n == ORDER[1] ? 2.5 : 1.8, ls = n == ORDER[6] ? :dash : :solid)
    end
    scatter!(p, pts.day[ii], o[ii]; color = :black, ms = 3.5, mc = :white, msc = :black)
    push!(ps, p)
end
key_ = plot(; framestyle = :none, legend = :top, legend_column = 4, legendfontsize = 10, foreground_color_legend = nothing, xlims = (0, 1), ylims = (0, 1))
for n in names_present
    plot!(key_, [-1, -1], [-1, -1]; color = cols[n], lw = 4, ls = n == ORDER[6] ? :dash : :solid, label = n)
end
scatter!(key_, [-1], [-1]; mc = :white, msc = :black, label = "observed (held out)")
savefig(plot(key_, plot(ps...; layout = (2, 3)); layout = @layout([a{0.1h}; b]), size = (1400, 760), margin = 4Plots.mm),
        joinpath(out, "heldout_predictions.png"))

# ---------------- report ----------------
runtime = round((time() - t_start) / 60; digits = 1)
nobs = nrow(pts)
mech_note = have_mech ? "per-point predictions exported from the zoo's own held-out code (<code>scripts/export_mech_heldout.jl</code>; $(N._esc(mech_src)))" :
            "<b>per-point predictions were not available</b>; only the zoo's summary RMSE (below) can be reported, so no bootstrap CI for it"
zoo_line = zoo_summary === nothing ? "" : N.para("For reference, the zoo's own summary file <code>heldout_ce1.csv</code> reports held-out RMSE $(round(zoo_summary.rmse[1]; digits = 1)) cells (nRMSE $(round(zoo_summary.nrmse[1]; digits = 3))) on $(zoo_summary.n[1]) points; it counts the day-0 rows too, which are the starting values (zero error by construction), so it is slightly lower than the days ≥ 1 value in the table.")
T_TBL = "RMSE in cells and nRMSE (RMSE ÷ mean of the $(nobs) held-out counts, $(round(mean(o); digits = 0)) cells) weight the large late counts most; log-RMSE weights every point's relative error equally (0.3 ≈ 35% typical error). Neural models: the seed-ensemble (geometric-mean) prediction, with mean ± SD of the single-seed scores. Δ vs best = this model minus the best model on that metric, with a paired-bootstrap 95% CI over the 6 trajectories; a CI that includes 0 means the data cannot separate the two."
const W_ZOO = "The zoo's Stage 3a structure, weighting scheme and drug model were chosen with Ce1 in the fitted data, so its score on Ce1 is somewhat optimistic; only its crowding parameters are refitted without Ce1 here (see Remaining differences between the arms)."
N.write_report(joinpath(out, "report.html"), "Three-way comparison on the held-out treated cocultures",
    [("What is compared", N.para("All models predict the <b>same $(nobs) held-out points</b>: the 30k treated (1 µM) cocultures Ce1 — 3 seeding ratios × 2 lines, days 1–14 — from the mechanistic zoo's tidy table (<code>02_mechanical_automatic_package/data/measurements.csv</code>). None of them saw Ce1 while fitting. Every held-out series starts at day 0 from its nominal seeding (100 cells per 30k well × the line's share), as in the zoo's validation; day 0 is not scored.") *
        N.para("<b>Training data (matched)</b>: the zoo's own training data for this prediction, read from the same table — the 30k monocultures (untreated, 0.62, 1 and 1.47 µM) for growth and drug response, and the 30k untreated cocultures Ce0 for competition. Training series whose first count is more than 3<sup>d</sup>-fold off nominal get an estimated starting value ($(nstart_train) series), the zoo's rule.") *
        N.para("<b>Models</b>: (a) the zoo's carried held-out prediction — crowding refitted to Ce0, kills from the monoculture drug stage — $(mech_note); (b) this track's Hill kill and the closed form recovered from the UDE in Stage 3 ($(N._esc(closed_formula))), refitted here on the matched data; (c) the UDE (logistic backbone + 3-8-8-1 NN kill), $(NSEEDS) seeds; (d) the pure neural ODE (per-capita growth = NN(ρ, v, line, low, Q), no time input, no competition term, so it is fitted to the monocultures and Ce0 together), $(NSEEDS) seeds; and a no-drug reference. Bootstrap: $(NBOOT) resamples. Runtime: $(runtime) min.")),
     ("Result", N.html_table(tbl; tip = (T_TBL, "Lower RMSE, nRMSE and log-RMSE are better; Δ vs best closer to 0 is better.", W_ZOO)) * zoo_line *
        N.figure("heldout_predictions.png", "Held-out Ce1 predictions (lines) against the data (open points), cells, linear axis"; tip = ("Each panel is one held-out trajectory (S = A2780, R = A2780cis). Lines are predictions made without Ce1.", "Better: the line runs through the points."))),
     ("Per-trajectory RMSE (cells)", N.html_table(pertraj; tip = ("RMSE in cells of each model on each held-out trajectory; mean_obs is that trajectory's mean count. Shows which trajectories drive the overall ranking.", :lower))),
     ("Does more NN-track training data change the ranking?", N.para("The same model families fitted to the full NN-track training set of Stages 2–3 (20k and 30k wells, intermittent schedules, low-resource runs; raw-image rebuild with the NN track's QC) and scored on the same points, next to the matched fits. See also the learning curve (Stage 7).") *
        N.html_table(ftbl; tip = ("Held-out error of each model family under the two training sets.", :lower))),
     ("Remaining differences between the arms", N.para("""
        <b>Loss</b>: the zoo fits squared errors in cells (with the run's residual weighting); the NN track fits squared errors of log counts (σ = √(0.15² + (SEM/mean)²)). Each arm is therefore favoured on the metric matching its loss: cells-RMSE for the zoo, log-RMSE for the NN track. Both metrics are reported for that reason.<br>
        <b>Model choice</b>: the zoo's Stage 3a structure, its weighting scheme and drug model were chosen with Ce1 in the fitted data (stated in the zoo's validation.toml); only its crowding parameters are refitted without Ce1 here. The NN-track structures (and the closed form's terms, chosen in Stage 3 on the full NN-track training set) were fixed without Ce1 by any automated step.<br>
        <b>What is refitted on Ce0</b>: the zoo refits a shared K, both crowding coefficients and a sensitive death rate; the NN-track UDE / Hill / closed form refit only the two competition coefficients (K stays per line from the monocultures); the pure neural ODE fits everything jointly on monocultures + Ce0.<br>
        <b>Observation model</b>: the NN-track models predict counts + a fitted detection floor (≈5 cells); the zoo predicts counts.<br>
        <b>Parameter counts</b>: zoo = free parameters along its path (Stage 1 growth + Stage 2b drug of each line + the Ce0 refit); NN track = parameters that act on these normal-medium 30k predictions (r, K per line, γ, floor, drug term, 2 competition coefficients); estimated starting values are excluded from both."""))])
println("wrote ", joinpath(out, "report.html"), " in $(runtime) min")
