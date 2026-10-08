# Verdict: the neural models are a competing method to the mechanistic ones. For every HELD-OUT test in the
# pipeline, compare the best neural model with the best mechanistic model. Reads only the CSVs written by
# run_02 … run_08, so it stays in step with the latest run.
# Counting rules (2026-10-07): only held-out tests are scored; in-sample fits and BIC are shown separately and
# never counted as wins; a difference smaller than the noise is a tie (see `tie_rule`).
include(joinpath(@__DIR__, "..", "src", "NNTrack.jl"))
using .NNTrack, DataFrames, CSV, Printf, Statistics
const N = NNTrack

rd(stage, f) = (p = joinpath(N.OUT_ROOT, stage, f); isfile(p) ? CSV.read(p, DataFrame; stringtype = String) : nothing)
out = N.outdir("verdict")

const NEURAL = Set(["UDE (NN kill)", "UDE (NN kill + NN resource)", "pure neural ODE", "pure neural ODE (autonomous)", "MLP ensemble", "boosted trees"])
# non-neural side: mechanistic ODEs including the no-drug model (the simplest mechanistic baseline), and for
# next-day prediction the persistence and ridge baselines
const MECH = Set(["Hill (mechanistic)", "closed form (optimal)", "no drug", "Hill kill (NN track, mechanistic)", "closed form (NN track, mechanistic)",
                  "mechanistic zoo (carried model)", "no drug (reference)", "persistence", "ridge"])
const TIE_ABS = 0.02     # log units (≈ 2%): differences below this are never called a win
tie_rule = "a difference is a tie when it is below $(TIE_ABS) log units or below the seed-to-seed SD of the neural score (or, for Stage 6, when the paired-bootstrap 95% CI includes 0)"
fold(x) = @sprintf("%.0f%%", 100 * (exp(x) - 1))    # log-scale RMSE -> typical fold error

rows = DataFrame(test = String[], why_it_matters = String[], metric = String[], best_neural = String[], neural = Float64[], neural_sd = Float64[],
                 best_mechanistic = String[], mechanistic = Float64[], difference = Float64[], winner = String[])
"Add one held-out test: df has columns key, val and (optionally) sd = SD over seeds."
function add!(test, why, metric, df, key, val; sd = nothing, ci_tie = nothing)
    nn = df[in.(df[!, key], Ref(NEURAL)), :]; mm = df[in.(df[!, key], Ref(MECH)), :]
    (nrow(nn) == 0 || nrow(mm) == 0) && return
    bn = nn[argmin(nn[!, val]), :]; bm = mm[argmin(mm[!, val]), :]
    s = sd === nothing ? 0.0 : coalesce(bn[sd], 0.0)
    d = bm[val] - bn[val]                     # positive = neural better
    tie = ci_tie === nothing ? abs(d) < max(TIE_ABS, s) : ci_tie(bn[key], bm[key])
    winner = tie ? "tie" : d > 0 ? "neural" : "mechanistic"
    push!(rows, (test, why, metric, bn[key], bn[val], s, bm[key], bm[val], d, winner))
end

# Stage 1 next-day (held-out treated co-culture): regression models and the ODE models' next-day error
nd = rd("nextday", "onestep_metrics.csv"); ude = rd("ude", "metrics.csv")
if nd !== nothing && ude !== nothing
    ndv = vcat(select(nd[nd.split .== "validation_treated_co", :], :model, :rmse, :rmse_sd_over_seeds => :sd),
               transform(select(ude[occursin.("VALIDATION", ude.set), :], :model, :nextday_rmse => :rmse), :rmse => (x -> zeros(length(x))) => :sd))
    add!("Next-day prediction, held-out treated co-culture", "short-term forecasting", "log RMSE", ndv, :model, :rmse; sd = :sd)
end

# Stage 2/3 forward prediction of the held-out ce1 co-culture (NN-track data), all models incl. the closed form
fw = rd("symbolic", "comparison_forward.csv")
if fw !== nothing && ude !== nothing
    fv = fw[occursin.("VALIDATION", fw.set), :]
    fv.sd = [let k = findfirst((ude.model .== m) .& occursin.("VALIDATION", ude.set)); k === nothing ? 0.0 : ude.forward_rmse_sd_over_seeds[k] end for m in fv.model]
    add!("Forward prediction, held-out ce1 co-culture (NN-track data, 20k + 30k)", "out-of-sample test on the NN track's own data", "log RMSE", fv, :model, :forward_rmse; sd = :sd)
end

# Hold-out folds: the control-relevant generalisation tests
hf = rd("ude", "holdout_folds.csv")
why_fold = Dict("hold out on→off intermittent" => "unseen schedule: drug first, then removed",
                "hold out off→on intermittent" => "unseen schedule: drug added late",
                "hold out 1 uM monoculture" => "unseen dose (interpolation)",
                "hold out all 20k wells" => "unseen seeding density")
if hf !== nothing
    for f in unique(hf.fold)
        add!("Hold-out: $f", get(why_fold, f, ""), "log RMSE", hf[hf.fold .== f, :], :model, :rmse; sd = :rmse_sd_over_seeds)
    end
end

# Stage 6: the fair three-way comparison on the zoo's held-out points (cells and log), CI-based ties
cmpt = rd("comparison", "comparison.csv")
if cmpt !== nothing
    ciinc(s) = s == "best" ? true : (m = match(r"\[([-+0-9.]+), ([-+0-9.]+)\]", s); m === nothing ? false : parse(Float64, m[1]) <= 0 <= parse(Float64, m[2]))
    for (metric, val, cicol) in (("RMSE (cells)", :rmse_cells, :d_rmse_ci95), ("log RMSE", :log_rmse, :d_logrmse_ci95))
        c = cmpt[cmpt.model .!= "no drug (reference)", :]
        lookup(n) = c[findfirst(==(n), c.model), cicol]
        # tie when the worse of the two has a CI (vs the best model) that includes 0
        # (when neither is the best model, both CIs against the best must include 0: an approximation)
        tie(a, b) = begin
            ca, cb = lookup(a), lookup(b)
            ca == "best" ? ciinc(cb) : cb == "best" ? ciinc(ca) : (ciinc(ca) && ciinc(cb))
        end
        add!("Stage 6 three-way: held-out 30k Ce1, matched training data", "same points, same training data, same starting values", metric,
             c, :model, val; ci_tie = tie)
    end
end

CSV.write(joinpath(out, "verdict.csv"), rows)
println(rows)

# in-sample diagnostics: shown, never counted
diag = DataFrame(item = String[], model = String[], value = Float64[], note = String[])
if fw !== nothing
    for r in eachrow(fw[.!occursin.("VALIDATION", fw.set), :])
        push!(diag, ("forward fit, $(r.set)", r.model, r.forward_rmse, "in-sample: a more flexible model always fits its own training data better"))
    end
end
bic = rd("symbolic", "comparison_bic.csv")
bic === nothing || for r in eachrow(bic)
    push!(diag, ("BIC on training drug units", r.model, r.bic_drug, "fit per parameter, in-sample; lower is better"))
end
CSV.write(joinpath(out, "in_sample_diagnostics.csv"), diag)

# learning curve (Stage 7)
lc = rd("learning_curve", "learning_curve_paired.csv")
lc_txt = lc === nothing ? "not run" : join([begin
        d = sort(lc[lc.vs .== vs, :], :frac)
        "UDE − $(vs): $(join(["$(Int(100 * r.frac))%: $(@sprintf("%+.3f", r.d_log_rmse)) ± $(@sprintf("%.3f", r.d_log_rmse_sd))" for r in eachrow(d)], ", "))"
    end for vs in unique(lc.vs)], "<br>")

nw = count(==("neural"), rows.winner); mw = count(==("mechanistic"), rows.winner); tw = count(==("tie"), rows.winner)
li(v) = isempty(v) ? "<li>none</li>" : join(("<li>$(N._esc(x))</li>" for x in v))
neural_wins = [r.test * " — " * r.metric for r in eachrow(rows) if r.winner == "neural"]
mech_wins = [r.test * " — " * r.metric for r in eachrow(rows) if r.winner == "mechanistic"]
ties = [r.test * " — " * r.metric for r in eachrow(rows) if r.winner == "tie"]
tbl = copy(rows)
tbl.neural_fold_error = [r.metric == "log RMSE" ? fold(r.neural) : "—" for r in eachrow(tbl)]
tbl.mechanistic_fold_error = [r.metric == "log RMSE" ? fold(r.mechanistic) : "—" for r in eachrow(tbl)]

N.write_report(joinpath(out, "report.html"), "Neural vs mechanistic — held-out scorecard",
    [("How to read this", N.para("The neural models (boosted trees, MLP, a UDE with a neural kill term, a pure neural ODE) are a <b>competing method</b>. Each row takes the best neural and the best non-neural model on one <b>held-out</b> test. The non-neural side holds the mechanistic ODEs, including the no-drug model, and for next-day prediction also the persistence and ridge baselines. <code>difference</code> = mechanistic minus neural error (positive = neural better). Rule for ties: $(tie_rule). Every model is autonomous (no clock-time inputs or terms); the next-day regression features carry no clock time either.")),
     ("Scorecard (held-out tests only)", N.para("<b>Neural wins $nw, mechanistic wins $mw, ties $tw.</b>") *
        N.html_table(tbl; tip = ("Best neural vs best mechanistic error on each held-out test (log RMSE: 0.2 ≈ 22% typical error; cells: RMSE in cells). neural_sd = SD over seeds. Ties follow the rule above.", "Lower error is better; a positive difference means the neural model was better."))),
     ("Where each was better", "<p><b>Neural better (beyond the noise):</b></p><ul>" * li(neural_wins) * "</ul><p><b>Mechanistic better:</b></p><ul>" * li(mech_wins) *
        "</ul><p><b>Ties (difference within the noise):</b></p><ul>" * li(ties) * "</ul>"),
     ("Does the UDE gain with more data? (Stage 7)", N.para("Paired held-out log-RMSE difference, UDE minus mechanistic, on the same training subsets (negative = UDE better), mean ± SD over subsets:<br>$(lc_txt)") *
        N.para("See <a href=\"../learning_curve/report.html\">the learning curve</a>. The largest training set here is small (≤ 28 drug wells, 4 doses, 2 switching schedules), so the curve tests the trend within these data; it cannot show what happens with much larger data sets.")),
     ("In-sample diagnostics (not counted)", N.para("Training-set fits and BIC show flexibility and fit per parameter, not generalisation; a more flexible model always fits its own training data better. They are listed for completeness and are never counted as wins.") *
        N.html_table(diag; tip = ("Training-set forward RMSE (log) and BIC on the training drug units.", "Lower is better, but only as a description of the fit, not of predictive skill."))),
     ("What these data can and cannot say", N.para("""
        These tests cover 4 constant doses (0–1.47 µM), two on/off schedules at 1 µM, two seeding densities and one held-out treated co-culture condition.
        They do <b>not</b> test schedules with more than one switch, doses outside 0–1.47 µM, other cell lines, or data sets much larger than these;
        no statement here should be read as covering those. The three-way comparison (<a href="../comparison/report.html">Stage 6</a>) is the
        fairest single test: same held-out points, same training data, same starting-value rule.""")),
     ("Sources", N.para("Built from <code>nextday/onestep_metrics.csv</code>, <code>ude/metrics.csv</code>, <code>ude/holdout_folds.csv</code>, <code>symbolic/comparison_forward.csv</code>, <code>symbolic/comparison_bic.csv</code>, <code>comparison/comparison.csv</code> and <code>learning_curve/learning_curve_paired.csv</code> of the latest run."))])
println("wrote ", joinpath(out, "report.html"))
