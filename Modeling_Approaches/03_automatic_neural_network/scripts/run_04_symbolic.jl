# Stage 3: recover a closed-form kill term from the UDE, refit candidates as ODEs, keep the lowest BIC.
# Needs outputs/ude/models.json from run_03_ude.jl. Run with threads.
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Plots, Statistics, StableRNGs, JSON3
const N = NNTrack
const FAST = get(ENV, "NN_FAST", "0") == "1"
const ITERS = FAST ? (adam_steps = 40, lbfgs_iters = 60) : (adam_steps = 150, lbfgs_iters = 300)
const NBOOT = FAST ? 3 : 8
const SCREEN_ITERS = FAST ? (adam_steps = 20, lbfgs_iters = 40) : (adam_steps = 80, lbfgs_iters = 200)

t_start = time()
trajs = N.load_all()
out = N.outdir("symbolic")
drug = N.make_units(N.training(trajs; role = :drug))
comp = N.make_units(N.training(trajs; role = :competition))
val = N.make_units(N.validation(trajs))
models = N.load_models(joinpath(N.OUT_ROOT, "ude", "models.json"))
ude = models["UDE (NN kill)"]

# ---------- 1-2: sample NN kill, sparse regression per line ----------
samples = N.kill_samples(ude, drug)
fronts = Dict{Float64, Any}()
tables = String[]
pf = plot(; xlabel = "number of terms", ylabel = "RMSE vs NN kill (1/day)", title = "Best fit to the NN kill per number of terms", yscale = :log10)
for (line, lname) in [(-1.0, "A2780 (S)"), (1.0, "A2780cis (R)")]
    s = samples[line]
    subs = N.best_subsets(s; maxterms = 3)
    front = combine(groupby(subs, :nterms), first)
    sweep = N.stlsq_sweep(s)
    fronts[line] = (front = front, vmax = s.vmax, ρrange = s.ρrange)
    CSV.write(joinpath(out, "subsets_$(line < 0 ? "S" : "R").csv"), select(subs, Not([:terms, :coefs])))
    plot!(pf, front.nterms, max.(front.rmse, 1e-5); marker = :circle, label = lname)
    push!(tables, "<h3>$lname</h3>" *
        N.para("Sampled region: v ∈ [0, $(round(s.vmax; digits = 2))] uM, ρ ∈ [$(round(s.ρrange[1]; digits = 2)), $(round(s.ρrange[2]; digits = 2))]. Ensemble spread (mean sd of NN kill): $(round(mean(s.ksd); sigdigits = 3)) /day.") *
        N.html_table(select(front, :nterms, :formula, :rmse, :r2); tip = ("Best formula per number of terms, fitted to the sampled NN kill (RMSE in 1/day; r2 = share of variance explained).", "Lower RMSE / higher r2 is better.")) *
        "<p class=\"note\">STLSQ (SINDy) sweep, for comparison:</p>" * N.html_table(sweep; tip = ("Sparse regression path over the threshold λ.", "Lower RMSE with fewer terms is better.")))
end
savefig(pf, joinpath(out, "pareto.png"))

# ---------- 3: refit every (S form, R form) candidate as an ODE on the data; rank by BIC ----------
npts = sum(sum(!isnan, y) for u in drug for y in u.y)
nstart = sum(N.nfree(u) for u in drug)      # estimated starting values, counted in k (same for every model)
function fit_one(tS, tR, c0S, c0R, θback; iters = ITERS)
    N.SYMBOLIC_TERMS[] = (S = tS, R = tR)
    return N.fit_drug_model(:symbolic, drug; θinit = vcat(θback[1:N.NB], c0S, c0R), iters...)
end
bic_of(f, k) = sum(N.unit_sse(:symbolic, f.θ, ones(4), u, l0) for (u, l0) in zip(drug, f.l0s)) + (N.NB + k + nstart) * log(npts)
screen = DataFrame(kill_S = String[], kill_R = String[], kill_params = Int[], bic_drug = Float64[])
cands = []
for rs in eachrow(fronts[-1.0].front), rr in eachrow(fronts[1.0].front)
    rs.nterms == 0 && continue
    f = fit_one(rs.terms, rr.terms, rs.coefs, rr.coefs, ude[1].θ; iters = SCREEN_ITERS)
    k = length(rs.terms) + length(rr.terms)
    push!(screen, (rs.formula, rr.formula, k, bic_of(f, k)))
    push!(cands, (tS = rs.terms, tR = rr.terms, f = f))
    @info "screened" rs.formula rr.formula screen.bic_drug[end]
end
order = sortperm(screen.bic_drug); screen = screen[order, :]; cands = cands[order]
# selection rule (user, 2026-09-29): fewest free parameters among candidates within 3 BIC of the lowest;
# a lower BIC breaks ties
const RULE_BIC = 3.0
screen.delta_bic = screen.bic_drug .- minimum(screen.bic_drug)
eligible = findall(screen.delta_bic .<= RULE_BIC)
pick = eligible[argmin([(screen.kill_params[i], screen.bic_drug[i]) for i in eligible])]
screen.chosen = [i == pick ? "chosen" : "" for i in 1:nrow(screen)]
CSV.write(joinpath(out, "candidate_bic.csv"), screen)
win = cands[pick]

# bootstrap the winner
N.SYMBOLIC_TERMS[] = (S = win.tS, R = win.tR)
fits = Vector{Any}(undef, NBOOT)
Threads.@threads for i in 1:NBOOT
    rng = StableRNG(2000 + i)
    w = i == 1 ? ones(length(drug)) : (e = -log.(rand(rng, length(drug))); e ./ sum(e) .* length(drug))
    fits[i] = N.fit_drug_model(:symbolic, drug; θinit = copy(win.f.θ), weights = w, ITERS...)
end
C = reduce(hcat, [f.θ[N.NB + 1:end] for f in fits])
lo = [quantile(C[j, :], 0.05) for j in 1:size(C, 1)]; hi = [quantile(C[j, :], 0.95) for j in 1:size(C, 1)]
coef_rows = DataFrame(line = String[], term = String[], estimate = Float64[], lo = Float64[], hi = Float64[], includes_zero = Bool[])
for (j, (l, nm)) in enumerate(vcat([("S", N.LIBRARY[t].name) for t in win.tS], [("R", N.LIBRARY[t].name) for t in win.tR]))
    push!(coef_rows, (l, nm, fits[1].θ[N.NB + j], lo[j], hi[j], lo[j] <= 0 <= hi[j]))
end
CSV.write(joinpath(out, "coefficients.csv"), coef_rows)
opt = (terms = (S = win.tS, R = win.tR), models = [N.calibrate(:symbolic, f.θ, comp) for f in fits], fit = fits[1])

# ---------- 4: compare with UDE / Hill ----------
# held-out wells start from the nominal seeding; the closed form is a deterministic single fit (bootstrap
# member 1, unit weights); the other models' scores come from Stage 2 (mean over seeds for neural models)
sets = ["train: drug units" => (drug, false), "train: untreated co-culture" => (comp, false),
        "VALIDATION: ce1 treated co-culture" => (val, true)]
cmp = DataFrame(model = String[], kill_params = Int[], bic_drug = Float64[])
metric = DataFrame(model = String[], set = String[], forward_rmse = Float64[])
N.SYMBOLIC_TERMS[] = opt.terms
k = length(opt.terms.S) + length(opt.terms.R)
push!(cmp, ("closed form (optimal)", k, bic_of(opt.fit, k)))
for (sname, (us, nominal)) in sets
    r = reduce(vcat, [N.forward_residuals(opt.models[1], u; nominal) for u in us])
    push!(metric, ("closed form (optimal)", sname, N.rmse(r.obs .- r.pred)))
end
for r in eachrow(CSV.read(joinpath(N.OUT_ROOT, "ude", "metrics.csv"), DataFrame))
    push!(metric, (r.model, r.set, r.forward_rmse))
end
for r in eachrow(CSV.read(joinpath(N.OUT_ROOT, "ude", "bic.csv"), DataFrame))
    push!(cmp, (r.model, r.kill_params, r.bic))
end
sort!(cmp, :bic_drug)
CSV.write(joinpath(out, "comparison_bic.csv"), cmp)
CSV.write(joinpath(out, "comparison_forward.csv"), metric)
println(screen); println(cmp); println(unstack(metric, :model, :set, :forward_rmse)); println(coef_rows)

open(joinpath(out, "symbolic_models.json"), "w") do io
    JSON3.write(io, Dict("optimal" => Dict("termsS" => opt.terms.S, "termsR" => opt.terms.R,
        "formulaS" => join([N.LIBRARY[j].name for j in opt.terms.S], " + "),
        "formulaR" => join([N.LIBRARY[j].name for j in opt.terms.R], " + "),
        "members" => [Dict("theta" => m.θ, "a" => m.a) for m in opt.models])))
end

# ---------- figures ----------
gr()
ps = []
for (line, lname) in [(-1.0, "A2780 (S)"), (1.0, "A2780cis (R)")]
    ρlo, ρhi = fronts[line].ρrange
    vg = range(0, max(1.6, fronts[line].vmax); length = 80)
    for ρ in (ρlo, (ρlo + ρhi) / 2, ρhi)
        p = plot(; title = "$lname, ρ = $(round(ρ; digits = 2))", titlefontsize = 8, xlabel = "v (uM)", ylabel = "kill (1/day)", legend = :topleft, legendfontsize = 6)
        K = reduce(hcat, [[N.kill_rate(:nn, m.θ, v, ρ, line) for v in vg] for m in ude])
        plot!(p, vg, vec(mean(K; dims = 2)); ribbon = vec(std(K; dims = 2)), color = :royalblue, fillalpha = 0.2, label = "UDE (mean ± sd over seeds)")
        Kc = reduce(hcat, [[N.kill_rate(:symbolic, m.θ, v, ρ, line) for v in vg] for m in opt.models])
        plot!(p, vg, vec(mean(Kc; dims = 2)); ribbon = vec(std(Kc; dims = 2)), color = :crimson, fillalpha = 0.2, lw = 2, label = "closed form (optimal, refit)")
        vspan!(p, [fronts[line].vmax, maximum(vg)]; color = :gray, alpha = 0.15, label = "outside data")
        push!(ps, p)
    end
end
savefig(plot(ps...; layout = (2, 3), size = (1300, 700), margin = 4Plots.mm), joinpath(out, "kill_closed_forms.png"))

eq(sm) = "kill<sub>S</sub> = $(isempty(sm.terms.S) ? "0" : join(["$(round(sm.fit.θ[N.NB + i]; sigdigits = 3))·$(N.LIBRARY[j].name)" for (i, j) in enumerate(sm.terms.S)], " + "))<br>" *
         "kill<sub>R</sub> = $(isempty(sm.terms.R) ? "0" : join(["$(round(sm.fit.θ[N.NB + length(sm.terms.S) + i]; sigdigits = 3))·$(N.LIBRARY[j].name)" for (i, j) in enumerate(sm.terms.R)], " + "))"
runtime = round((time() - t_start) / 60; digits = 1)
N.write_report(joinpath(out, "report.html"), "Stage 3 — analytical form of the learned kill term",
    [("Optimal closed form",
      N.para("<span class='eqg-text'>" * eq(opt) * "</span>") *
      N.para("v = filtered exposure (uM), <span class='eqg-text' data-vars='v,γ,u'>dv/dt = γ(u − v)</span>; ρ = relative density; <span class='eqg-text' data-vars='H,v,c,n'>H(v; c, n) = vⁿ/(cⁿ + vⁿ)</span>. The NN kill term was sampled from the mean of the $(length(ude)) UDE seeds. Every combination of the best 1–3-term A2780 form and 0–3-term A2780cis form was refit as an ODE on the data and ranked by BIC. The optimal form is the one with the fewest free parameters among candidates within 3 BIC of the lowest (a lower BIC breaks ties). Its coefficients carry Bayesian-bootstrap 5–95% intervals; terms whose interval includes 0 are flagged, not removed. Runtime: $(runtime) min.") *
      N.html_table(coef_rows; tip = ("Coefficient of each kill term with its bootstrap 5–95% interval; includes_zero = the data do not show the term is needed.", :none))),
     ("Candidate ranking (BIC on training drug units; lower is better)", N.html_table(screen; tip = ("Each (A2780, A2780cis) closed-form pair refit as an ODE. delta_bic = BIC minus the lowest; the chosen form has the fewest parameters within 3 of the lowest.", :lower))),
     ("Forward RMSE against the NN and Hill (log scale)", N.html_table(unstack(metric, :model, :set, :forward_rmse); tip = ("Forward-prediction RMSE of log counts (neural models: mean over seeds). Training columns show fit; the VALIDATION column is out of sample, from the nominal seeding.", :lower))),
     ("BIC against the other models", N.html_table(cmp; tip = ("BIC on the training drug units (best seed for neural models). The pure neural ODE also saw the untreated co-culture.", :lower))),
     ("Closed form against the NN", N.figure("kill_closed_forms.png", "Kill rate vs exposure at low/mid/high density; bands = sd over seeds (UDE) or bootstrap members (closed form); shaded = exposures never reached in the data"; tip = ("Kill rate per day of the learned NN term and of the refitted closed form.", "Closer agreement between the two curves inside the unshaded region is better."))),
     ("Sparse regression detail", N.figure("pareto.png", "Best RMSE per number of terms"; tip = ("How closely the best k-term formula reproduces the NN kill surface.", :lower)) * join(tables))])
println("wrote ", joinpath(out, "report.html"))
