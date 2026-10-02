# Verdict: the neural models are a competing method to the mechanistic ones. For every test in the
# pipeline, compare the best neural model with the best mechanistic model and say where each should be
# used. Reads only the CSVs written by run_02 … run_06, so it stays in step with the latest run.
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Printf
const N = NNTrack

rd(stage, f) = CSV.read(joinpath(N.OUT_ROOT, stage, f), DataFrame)
out = N.outdir("verdict")

const NEURAL = Set(["UDE (NN kill)", "pure neural ODE", "MLP ensemble"])
const MECH = Set(["Hill (mechanistic)", "closed form (optimal)", "no drug"])
fold(x) = @sprintf("%.0f%%", 100 * (exp(x) - 1))    # log-scale RMSE -> typical fold error

rows = DataFrame(test = String[], why_it_matters = String[], best_neural = String[], neural = Float64[],
                 best_mechanistic = String[], mechanistic = Float64[], winner = String[], margin = String[])
function add!(test, why, df, key, val)
    nn = df[in.(df[!, key], Ref(NEURAL)), :]; mm = df[in.(df[!, key], Ref(MECH)), :]
    (nrow(nn) == 0 || nrow(mm) == 0) && return
    bn = nn[argmin(nn[!, val]), :]; bm = mm[argmin(mm[!, val]), :]
    d = bm[val] - bn[val]
    winner = abs(d) < 0.01 ? "tie" : d > 0 ? "neural" : "mechanistic"
    push!(rows, (test, why, bn[key], bn[val], bm[key], bm[val], winner, @sprintf("%+.3f", d)))
end

# Stage 1 (one-step) and the ODE next-day errors on the held-out co-culture
nd = rd("nextday", "onestep_metrics.csv")
ude = rd("ude", "metrics.csv")
ndv = vcat(select(nd[nd.split .== "validation_treated_co", :], :model, :rmse),
           rename(select(ude[occursin.("VALIDATION", ude.set), :], :model, :nextday_rmse), :nextday_rmse => :rmse))
add!("Next-day prediction, held-out treated co-culture", "short-term forecasting", ndv, :model, :rmse)

# Stage 2/3 forward prediction, all models incl. the recovered closed form
fw = rd("symbolic", "comparison_forward.csv")
bb = ude[ude.model .== "pure neural ODE", :]
fw = vcat(fw, DataFrame(model = bb.model, set = bb.set, forward_rmse = bb.forward_rmse))
for (set, why) in [("VALIDATION: ce1 treated co-culture", "the headline out-of-sample test"),
                   ("train: untreated co-culture", "competition dynamics (fit, not held out)"),
                   ("train: drug units", "in-sample fit (flexibility, not generalisation)")]
    add!("Forward prediction — $(replace(set, "VALIDATION: " => "held out: "))", why, fw[fw.set .== set, :], :model, :forward_rmse)
end

# Hold-out folds: the control-relevant generalisation tests
hf = rd("ude", "holdout_folds.csv")
why_fold = Dict("hold out on→off intermittent" => "unseen schedule: drug first, then removed",
                "hold out off→on intermittent" => "unseen schedule: drug added late (delayed-onset dynamics)",
                "hold out 1 uM monoculture" => "unseen dose (interpolation) — control explores doses",
                "hold out all 20k wells" => "unseen seeding density")
for f in unique(hf.fold)
    add!("Hold-out: $f", get(why_fold, f, ""), rename(hf[hf.fold .== f, :], :rmse => :val), :model, :val)
end

# Parameter efficiency
bic = rd("symbolic", "comparison_bic.csv")
add!("BIC on training drug units (lower = better fit per parameter)", "evidence per parameter; 12–28 wells", rename(bic, :bic_drug => :val), :model, :val)

CSV.write(joinpath(out, "verdict.csv"), rows)
println(rows)

rows.evidence = [occursin("train:", t) ? "in-sample fit" : occursin("BIC", t) ? "fit per parameter" : "held out" for t in rows.test]
ho = rows[rows.evidence .== "held out", :]
nw = count(==("neural"), ho.winner); mw = count(==("mechanistic"), ho.winner); tw = count(==("tie"), ho.winner)
tag(r) = r.evidence == "held out" ? r.test : "$(r.test) ($(r.evidence))"
neural_wins = [tag(r) for r in eachrow(rows) if r.winner == "neural"]
mech_wins = [tag(r) for r in eachrow(rows) if r.winner == "mechanistic"]
li(v) = isempty(v) ? "<li>none</li>" : join(("<li>$(N._esc(x))</li>" for x in v))

tbl = copy(rows)
tbl.neural_fold_error = fold.(tbl.neural); tbl.mechanistic_fold_error = fold.(tbl.mechanistic)
tbl[occursin.("BIC", tbl.test), :neural_fold_error] .= "—"; tbl[occursin.("BIC", tbl.test), :mechanistic_fold_error] .= "—"

N.write_report(joinpath(out, "report.html"), "Neural vs mechanistic — where each should be used",
    [("How to read this", N.para("The neural models (MLP, UDE with a neural kill term, pure neural ODE) are a <b>competing method</b>, not a step towards the mechanistic model. Each row takes the best neural and the best mechanistic model on one test. Errors are forward-prediction RMSE on the log scale (0.2 ≈ 22% typical error); <code>margin</code> is how much the neural model is better (positive) or worse (negative). Differences under 0.01 count as a tie. The closed form counts as mechanistic even though its terms were suggested by the neural kill term.")),
     ("Scorecard", N.para("<b>On held-out tests</b> (data the model never saw): neural wins $nw, mechanistic wins $mw, ties $tw. Rows marked <i>in-sample fit</i> only show flexibility — a more flexible model always fits its own training data better — and BIC weighs fit against the number of parameters.") * N.html_table(tbl)),
     ("Where the neural model is the better choice", "<ul>" * li(neural_wins) * "</ul>" * N.para("""
        <b>Interpretation.</b> The neural kill term (UDE) is most useful where the drug response has a shape the Hill curve or a
        single exposure filter cannot express — for example a delayed onset when the drug is added late (the off→on hold-out) and
        the density dependence of the kill in co-culture. In those settings it predicts held-out wells as well as or slightly better
        than the mechanistic models. It was also useful as a <b>discovery tool</b>: its learned kill surface suggested the
        density-dependent terms (v·ρ, v²·ρ) of the closed form, which then beat the Hill model on BIC.""")),
     ("Where the mechanistic model is the better choice", "<ul>" * li(mech_wins) * "</ul>" * N.para("""
        <b>Interpretation.</b> Mechanistic models win when the model must <b>extrapolate</b>: an unseen dose, an unseen seeding
        density, or an unseen schedule where the drug is given first. With 12–28 training wells, 4 dose levels and 2 switching
        schedules, the ~110 neural weights are weakly constrained outside the training conditions (the pure neural ODE fails on the
        unseen 1 uM dose), and BIC penalises them heavily. Mechanistic parameters (growth rate, carrying capacity, EC50-like
        terms) also have biological meaning and bootstrap intervals, which the neural weights do not.""")),
     ("Recommendation by task", N.para("""
        <b>Short-term (next-day) forecasting:</b> any of the ODE models or boosted trees; neural adds nothing measurable.<br>
        <b>Predicting held-out wells under tested conditions (same doses/schedules):</b> the UDE is competitive and sometimes best —
        a reasonable choice, ideally as an ensemble alongside the closed form.<br>
        <b>Optimal control (new schedules and doses):</b> use the mechanistic closed form as the primary model, and the UDE as a
        cross-check — where they disagree (e.g. A2780 monoculture schedules), treat the prediction as uncertain.<br>
        <b>Understanding mechanism / reporting parameters:</b> mechanistic.<br>
        <b>When the neural model would start to pay off:</b> many more conditions (doses, schedules, densities, time-varying
        dosing), or inputs with no obvious mechanistic form (e.g. image features).""")),
     ("Sources", N.para("Built from <code>nextday/onestep_metrics.csv</code>, <code>ude/metrics.csv</code>, <code>ude/holdout_folds.csv</code>, <code>symbolic/comparison_forward.csv</code> and <code>symbolic/comparison_bic.csv</code> of the latest run."))])
println("wrote ", joinpath(out, "report.html"))
