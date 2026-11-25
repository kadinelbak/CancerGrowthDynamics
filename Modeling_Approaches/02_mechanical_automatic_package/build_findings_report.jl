# Key findings, assumptions register and limitations — one page generated from the outputs of the latest run
# (Final Report chapter 8). Every number is read from the run's files, so the page never goes stale.
#   julia --project=. build_findings_report.jl
# Reads outputs/zoo_fits/{zoo_model_selection,zoo_parameters,starting_values,heldout_ce1}.csv, zoo_state.jls,
# outputs/zoo_analyses/*.csv, outputs/control/control_outcomes.csv, outputs/control_<long>d/fasting_outcomes.csv,
# outputs/sensitivity/hill/** (drug-response shape) and the neural-network comparison, when present.
include(joinpath(@__DIR__, "zoo_common.jl")); include(joinpath(@__DIR__, "zoo_models.jl")); include(joinpath(@__DIR__, "bootstrap_tools.jl"))
include(joinpath(@__DIR__, "info_tips.jl"))
include(joinpath(@__DIR__, "fit_set.jl"))              # which saved fit set the outputs use (user, 2026-11-23)
using Serialization, Printf, Statistics

const OUTD = joinpath(@__DIR__, "outputs")
const DOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "findings")
const CCFG = TOML.parsefile(joinpath(@__DIR__, "control_config.toml"))
const LONG = round(Int, CCFG["schedule"]["long_days"])
const HEADLINE = String(get(get(CCFG, "scenarios", Dict{String,Any}()), "headline", "A"))
readcsv(p) = isfile(p) ? CSV.read(p, DataFrame; stringtype=String) : nothing
h(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
f1(x) = x isa Number && isfinite(x) ? @sprintf("%.1f", x) : "—"
f0(x) = x isa Number && isfinite(x) ? @sprintf("%.0f", x) : "—"
g3(x) = x isa Number && isfinite(x) ? @sprintf("%.3g", x) : "—"
th(label, text, better) = "<th>$(h(label)) $(tip(text, better))</th>"
table(head, rows) = "<div class='tscroll'><table><tr>" * join(head) * "</tr>" * join(rows) * "</table></div>"
tr(cells) = "<tr>" * join("<td>$c</td>" for c in cells) * "</tr>"

sel = readcsv(joinpath(OUTD, "zoo_fits", "zoo_model_selection.csv"))
par = readcsv(joinpath(OUTD, "zoo_fits", "zoo_parameters.csv"))
x0t = readcsv(joinpath(OUTD, "zoo_fits", "starting_values.csv"))
hel = readcsv(joinpath(OUTD, "zoo_fits", "heldout_ce1.csv"))
st = deserialize(joinpath(OUTD, "zoo_fits", "zoo_state.jls"))
grp(id) = st.groups[findfirst(g -> g.id == id, st.groups)]
pval(g, n) = (r = par[(par.group .== g) .& (par.param .== n), :]; nrow(r) == 0 ? NaN : r.estimate[1])

# ── 1. carried models ─────────────────────────────────────────────────────────
function glance()
    rows = String[]
    for r in eachrow(sel[sel.optimal .== true, :])
        line = endswith(r.group, "-a2780") ? " (A2780)" : endswith(r.group, "-a2780cis") ? " (A2780cis)" : ""
        push!(rows, tr([h(stage_label(r.stage)) * line, h(r.label), string(r.k), f1(r.bic), f0(r.pct_chosen) * "%"]))
    end
    "<h2>Carried model at each stage $(tip("The model each stage carries forward: the fewest free parameters among the candidates within 2 BIC of the lowest. 'Chosen in' is how often the same rule picks it when the data are resampled (bootstrap generated from the carried model).", :none))</h2>" *
    table([th("Stage", "Fitting stage (thesis numbering).", :none), th("Carried model", "Selected model.", :none),
           th("k", "Free model parameters (estimated starting values are counted in the BIC but not here).", "Fewer is better at equal fit."),
           th("BIC", "Weighted BIC of the carried model; compare only within a stage.", :lower),
           th("Chosen in", "Share of bootstrap replicates in which the rule picks the same model.", :higher)], rows)
end

# ── 2. transportability ───────────────────────────────────────────────────────
# BIC of a monoculture drug model with its kill set to zero (k = 0), on the fit's own weighted scale
function nokill(gid)
    g = grp(gid); c = g.cands[1]
    yh = bs_predict(c.model, c.u0, begin p = copy(c.params); i = findfirst(n -> startswith(n, "β"), c.pnames); i === nothing || (p[i] = 0.0); p end, g.datasets)
    s = sum(sum(((d.y .- y) ./ collect(dweights(d))) .^ 2) for (d, y) in zip(g.datasets, yh)); n = sum(length(d.x) for d in g.datasets)
    (fit=c.bic, none=n * log(s / n) + (hasproperty(c, :x0est) ? length(c.x0est) : 0) * log(n))
end
function transport()
    items = String[]
    s3 = sel[sel.group .== "s3", :]
    carried = s3[s3.optimal .== true, :candidate][1]
    nc = s3[endswith.(s3.candidate, "_noC"), :]; wc = s3[.!endswith.(s3.candidate, "_noC"), :]
    if nrow(nc) > 0 && nrow(wc) > 0
        bnc = minimum(nc.bic); bwc = minimum(wc.bic)
        push!(items, "<li><b>Resistant kill in coculture.</b> The best coculture model without any cisplatin term for A2780cis has BIC $(f1(bnc)); the best with one $(f1(bwc)) (Δ $(@sprintf("%+.1f", bnc - bwc))). The carried model is <code>$(h(carried))</code>$(endswith(carried, "_noC") ? ": in the treated cocultures cisplatin has no detectable effect on A2780cis at 1 µM" : ""). In monoculture the A2780cis kill is needed: without it the Stage 2b fit is $(@sprintf("%+.1f", (k = nokill("s2-a2780cis"); k.none - k.fit))) BIC worse. The monoculture drug response therefore does not transfer to coculture.</li>")
    end
    s4 = sel[sel.group .== "s4", :]
    if nrow(s4) > 0
        inh = s4[s4.candidate .== "s4_inherited", :]; c4 = s4[s4.optimal .== true, :]
        nrow(inh) > 0 && push!(items, "<li><b>Low resource + cisplatin (Stage 4).</b> The model assembled from the earlier stages, with nothing refitted, is $(@sprintf("%+.1f", inh.bic[1] - c4.bic[1])) BIC worse than the carried <code>$(h(c4.candidate[1]))</code> ($(h(join([(kv = split(x, "="); length(kv) == 2 ? kv[1] * " = " * g3(parse(Float64, kv[2])) : x) for x in split(c4.params[1], "; ")], ", ")))). The drug and low-resource mechanisms fitted separately do not combine additively.</li>")
    end
    if x0t !== nothing && nrow(x0t) > 0
        lr = x0t[occursin.("lowres", x0t.series) .| occursin.("combination", x0t.series), :]
        lo = lr[lr.fold .> 0.11, :]; atb = lr[lr.fold .<= 0.11, :]
        nrow(lr) > 0 && push!(items, "<li><b>Seeding scale of the low-resource experiments.</b> $(nrow(lo)) of the $(nrow(lr)) low-resource and combination series start $(@sprintf("%.1f", minimum(lo.fold)))–$(@sprintf("%.1f", maximum(lo.fold)))-fold above their nominal seeding (median $(@sprintf("%.1f", median(lo.fold)))): their first counts are far beyond one day of growth, so those plates were seeded or imaged at a different effective density. Their starting values are estimated; their capacity and low-resource effects should be compared with the normal-medium experiments with this in mind.$(nrow(atb) > 0 ? " The start of " * join(h.(atb.series), ", ") * " sits at its lower bound (0.1 × nominal): the model does not fit that series' early counts." : "")</li>")
    end
    "<h2>What did not carry over between conditions $(tip("Places where a mechanism fitted in a simpler condition fails in a more complex one. These are results, not fitting problems.", :none))</h2><div class='intro'><ul>" * join(items) * "</ul></div>"
end

# ── 3. starting values ────────────────────────────────────────────────────────
function starts()
    x0t === nothing && return ""
    e = x0t[x0t.estimated .== true, :]
    rows = [tr([h(r.group), h(r.series), f1(r.nominal), f1(r.start), @sprintf("%.2f", r.fold)]) for r in eachrow(e)]
    "<h2>Estimated starting values $(tip("Day 0 is the nominal seeding (100 cells at 30k, split by the seeding ratio). A series whose first count is more than 3-fold off it gets its starting value estimated; all others start at the nominal value.", :none))</h2>" *
    "<p>$(nrow(e)) of $(nrow(x0t)) series have an estimated starting value; the others start at the nominal seeding.</p>" *
    table([th("Stage group", "Fit group.", :none), th("Series", "Data series (source id).", :none), th("Nominal", "Nominal seeding of the line in that well.", :none),
           th("Estimated", "Fitted starting value.", :none), th("Fold", "Estimated ÷ nominal.", "Closer to 1 is better (consistent with the seeding).")], rows)
end

# ── 4. validation ─────────────────────────────────────────────────────────────
function validation()
    parts = String[]
    if hel !== nothing
        rows = [tr([h(r.row), f0(r.rmse), @sprintf("%.2f", r.nrmse)]) for r in eachrow(hel)]
        pr = hel[startswith.(hel.row, "prediction ("), :]; rf = hel[startswith.(hel.row, "untreated reference"), :]; df = hel[startswith.(hel.row, "prediction minus"), :]
        verdict = (nrow(pr) == 1 && nrow(rf) == 1) ? (pr.rmse[1] > rf.rmse[1] ?
            "<p><b>Transferring the monoculture drug response makes the prediction worse than ignoring the drug</b> ($(f0(pr.rmse[1])) vs $(f0(rf.rmse[1])) cells$(nrow(df) == 1 && isfinite(df.diff_lo[1]) ? "; difference 95% CI $(f0(df.diff_lo[1])) to $(f0(df.diff_hi[1])) over the 6 trajectories" : "")): the monoculture kills, which include a kill of A2780cis, do not describe the treated cocultures, where the fitted model has none.</p>" :
            "<p>The transferred monoculture drug response improves on ignoring the drug ($(f0(pr.rmse[1])) vs $(f0(rf.rmse[1])) cells).</p>") : ""
        push!(parts, verdict)
        push!(parts, "<h3>Held-out treated cocultures (Ce1) $(tip("Ce1 predicted from crowding fitted to the untreated cocultures and the monoculture kills, with no refit.", :none; warn=:ce1_selection))</h3>" *
              table([th("Prediction", "What is compared.", :none), th("RMSE (cells)", "Root-mean-square error.", :lower), th("nRMSE", "RMSE ÷ mean of the held-out data.", :lower)], rows))
    end
    an = joinpath(OUTD, "zoo_analyses", "analyses_state.jls")
    I = isfile(an) ? get(deserialize(an), :intermittent, nothing) : nothing
    if I !== nothing
        rows = [tr([h(L), h(string(sc)), f0(I[L].sc_rmse[sc]), @sprintf("%.2f", I[L].sc_nrmse[sc])]) for L in sort(collect(keys(I))) for sc in sort(collect(keys(I[L].sc_rmse)))]
        push!(parts, "<h3>Intermittent dosing (untouched validation) $(tip("Monocultures with the drug switched on or off at day 7, predicted from the Stage 1 + 2b model with no refit, the scale fixed at 1 and the nominal seeding. These data informed no modelling choice.", :none))</h3>" *
              table([th("Line", "Cell line.", :none), th("Schedule", "on0_off7: 1 µM for days 0–7; off0_on7: 1 µM from day 7.", :none),
                     th("RMSE (cells)", "Error of the no-refit prediction.", :lower), th("nRMSE", "RMSE ÷ mean of the data.", :lower)], rows) *
              "<p>Washout test $(warn(:hypothesis_test)): equal build-up and fade rates are kept for $(join(["$L (ΔBIC $(@sprintf("%+.1f", I[L].bic_koff - I[L].bic_equal)) for a separate fade rate)" for L in sort(collect(keys(I)))], ", ")).</p>")
    end
    "<h2>Validation $(tip("How well the models predict data they were not fitted to.", :none))</h2>" * join(parts)
end

# ── 5. control headline ───────────────────────────────────────────────────────
function control()
    c14 = readcsv(joinpath(OUTD, "control", "control_outcomes.csv"))
    cl = readcsv(joinpath(OUTD, "control_$(LONG)d", "fasting_outcomes.csv"))
    out = String[]
    function block(df, title, mpcname)
        rows = String[]
        for g in groupby(df[df.ratio .== "50:50", :], intersect(["environment", "scenario"], names(df)))
            m = g[g.policy .== mpcname, :]; nrow(m) == 0 && continue
            fixed = g[.!startswith.(g.policy, "MPC"), :]
            bt = fixed[argmin(fixed.total_final), :]; bc = fixed[argmin(fixed.C_final), :]
            lab = join([string(g[1, k]) for k in intersect(["environment", "scenario"], names(df))], " · ")
            sc = string(g[1, "scenario"]); lab *= sc == HEADLINE ? " (main result)" : sc == "B" ? " (what-if: monoculture drug response)" : " (sensitivity case)"
            ttp(x) = x isa Number && isfinite(x) ? f1(x) : "not progressed"
            push!(rows, tr([h(lab), ttp(m.time_to_progression_day[1]), f0(m.total_final[1]) * " (" * f0(m.C_final[1]) * ")", f1(m.cumulative_dose_uM_days[1]),
                            h(bt.policy) * ": " * f0(bt.total_final), h(bc.policy) * ": " * f0(bc.C_final)]))
        end
        "<h3>$title</h3>" * table([th("Run", "Environment / scenario, 50:50 start.", :none), th("MPC time to progression (day)", "When the controller's culture progresses; 'not progressed' within the run is best.", :higher),
                                    th("MPC: total (resistant)", "Cells at the end under the controller.", :lower),
                                    th("MPC dose (µM·d)", "Cumulative cisplatin.", :lower), th("Best fixed policy, total", "The fixed schedule with the fewest total cells.", :lower),
                                    th("Best fixed policy, resistant", "The fixed schedule with the fewest resistant cells.", :lower)], rows)
    end
    c14 !== nothing && push!(out, block(c14, "14 days, dose only", "MPC, measured doses only"))
    cl !== nothing && push!(out, block(cl, "$LONG days, dose + fasting", "MPC, dose + fasting (measured levels)"))
    "<h2>Optimal control against the best fixed policies $(tip("The controller maximises the time to progression (resistant cells at K/2 or sensitive cells at 85% of their untreated plateau) first; total cells, resistant cells, dose and host costs form a second tier, so a fixed schedule can beat it on a single end-of-run metric. Scenario A (the coculture fit) is the main result; B applies the monoculture drug response in coculture and is contradicted by the treated cocultures. Days after 14 are model extrapolation.", :none))</h2>" * join(out)
end

# ── 6. drug-response shape ────────────────────────────────────────────────────
function shape()
    hc = readcsv(joinpath(OUTD, "sensitivity", "hill", "control", "control_outcomes.csv"))
    hl = readcsv(joinpath(OUTD, "sensitivity", "hill", "control_$(LONG)d", "fasting_outcomes.csv"))
    hc === nothing && hl === nothing && return "<h2>Drug-response shape</h2><p class='muted'>The Hill-kill variant has not been run (sensitivity/run_hill_variant.sh).</p>"
    m14 = readcsv(joinpath(OUTD, "control", "control_outcomes.csv")); ml = readcsv(joinpath(OUTD, "control_$(LONG)d", "fasting_outcomes.csv"))
    rows = String[]
    for (lab, main, alt, pol) in (("14 days", m14, hc, "MPC, measured doses only"), ("$LONG days", ml, hl, "MPC, dose + fasting (measured levels)"))
        (main === nothing || alt === nothing) && continue
        for scen in ("A", "B")
            pick(df) = (r = df[(df.ratio .== "50:50") .& (string.(df.scenario) .== scen) .& (df.policy .== pol) .& ((:environment in propertynames(df)) ? (df.environment .== "normal") : trues(nrow(df))), :]; nrow(r) == 0 ? nothing : r[1, :])
            a, b = pick(main), pick(alt)
            (a === nothing || b === nothing) && continue
            push!(rows, tr([lab, scen, f0(a.total_final) * " (" * f0(a.C_final) * ")", f0(b.total_final) * " (" * f0(b.C_final) * ")"]))
        end
    end
    "<h2>Drug-response shape: linear (Norton–Simon) vs Hill $(tip("The bootstrap cannot separate a linear and a saturating (Hill) kill in Stage 2b, so the controller is run on both. Cells at the end under MPC, 50:50 start, normal medium.", :none))</h2>" *
    table([th("Horizon", "Run length.", :none), th("Scenario", "A: coculture fit; B: monoculture kills.", :none),
           th("Linear kill: total (resistant)", "Main model.", :lower), th("Hill kill: total (resistant)", "Variant with Stage 2b forced to the Hill kill.", :lower)], rows)
end

# ── 7. neural-network comparison ──────────────────────────────────────────────
function nn()
    nnd = joinpath(ROOT, "Modeling_Approaches", "03_automatic_neural_network", "outputs", "comparison")
    f = isdir(nnd) ? filter(x -> endswith(x, ".csv"), readdir(nnd)) : String[]
    link = "<p><a href='../../../../outputs/nn_track/comparison/report.html'>Full mechanistic vs neural-network comparison</a></p>"
    isempty(f) && return "<h2>Mechanistic vs neural-network models</h2><p class='muted'>The matched comparison has not been produced yet.</p>" * link
    df = readcsv(joinpath(nnd, first(sort(f))))
    "<h2>Mechanistic vs neural-network models $(tip("The same held-out treated-coculture points for every method (see the linked page for training-data differences and confidence intervals).", :lower))</h2>" *
    table([th(string(n), string(n), :none) for n in names(df)], [tr([x isa Number ? g3(x) : h(x) for x in r]) for r in eachrow(df)]) * link
end

# ── 8. assumptions register ───────────────────────────────────────────────────
function register()
    R = [("No explicit time terms; ≤ 3 compartments per mechanism", "Design choice (optimal control needs autonomous models)", "Not compared with time-based models in the results; the earlier time-based zoo is archived."),
         ("Fewest parameters within 2 BIC of the lowest", "Partly", "Bootstrap shows how often the rule re-picks each model (table above)."),
         ("Residual variance ∝ max(y, 50)^1.1 (pooled power law)", "Partly", "Earlier SD-weighted and unweighted fits were compared on the held-out cocultures; residuals are correlated day to day, so the effective number of points is below n."),
         ("Day 0 = nominal seeding; estimated when > 3-fold off", "Yes, per series", "See 'Estimated starting values'."),
         ("Shared cisplatin structure for both lines", "Imposed for the drug; tested for low resource", "Stage 2b chooses one structure on the joint BIC."),
         ("Drug kill linear and scaled by growth (Norton–Simon)", "Compared with Hill and with no coupling", "The Hill variant is carried through the control (above)."),
         ("Drug effect wears off at the rate it builds up", "Yes (intermittent data)", "A separate washout rate is not supported, but the test has little power."),
         ("Low resource not reset at media changes", "Yes", "Reset models not supported (14-day control report)."),
         ("Inherited parameters treated as known in later stages", "No", "Uncertainty is not propagated between stages."),
         ("Drug concentration constant between media changes (no breakdown)", "No", "—"),
         ("Toxicity, weight and fasting protection from the literature", "Swept", "The toxicity limit cannot bind at the default host recovery rate."),
         ("Results beyond day 14", "No data", "Days 15–$LONG are extrapolation; the model has no passaging.")]
    "<h2>Assumptions register $(tip("Each major assumption and whether the analysis tested it. Reviewers look for exactly this.", :none))</h2>" *
    table([th("Assumption", "What the models assume.", :none), th("Tested?", "Whether an alternative was fitted or swept.", "'Yes' is better."), th("Result / note", "Outcome of the test, or why it is untested.", :none)],
          [tr([h(a), h(b), h(c)]) for (a, b, c) in R])
end

function limits()
    """<h2>Limitations and the experiments that would resolve them $(tip("Every limitation in full; the ⚠ buttons in the chapters show the one that applies there.", :none))</h2><div class='intro'><ul>
<li>One biological run per condition (two for low resource); residuals are correlated day to day (shared imaging artefacts), so intervals are conditional and probably too narrow.</li>
<li>The treated cocultures informed the choice of weighting and drug model; the intermittent-dosing data are the only untouched validation set.</li>
<li>The cocultures were treated at one dose only; the coculture dose response and the resistant kill in coculture rest on that single dose. <b>Experiment:</b> a coculture dose series (e.g. 0.25, 0.5, 2, 4 µM).</li>
<li>Recovery from starvation stress after refeeding was never measured and drives the fasting result. <b>Experiment:</b> fast-then-refeed cultures with daily counts.</li>
<li>No data beyond day 14; long control runs are extrapolations without passaging. <b>Experiment:</b> 28-day cultures with media changes.</li>
<li>The low-resource plates start far above their nominal seeding (see above). <b>Check:</b> the seeding density and imaging field of those plates.</li>
<li>The 30k Ce0 75:25 resistant series reads about 0.8× of the raw-image rebuild; it is used as supplied.</li></ul></div>"""
end

const CSS = "body{font:14px/1.5 system-ui,sans-serif;margin:0;background:#f4f7f8;color:#17212b}main{max-width:1200px;margin:auto;padding:24px}h1{color:#124b6e}h2{color:#124b6e;margin-top:28px}" *
            ".intro{background:#fff;border:1px solid #d9e2e6;border-left:6px solid #176b87;border-radius:12px;padding:14px 18px;margin:12px 0}.muted{color:#5b6b73}" *
            "table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0;background:#fff}td,th{border:1px solid #d3dde2;padding:6px 8px;text-align:left;vertical-align:top}th{background:#e9f3f8}.tscroll{overflow-x:auto}code{background:#eef5f6;padding:1px 4px;border-radius:4px}"
html = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Key findings</title>
<style>$CSS</style>$INFO_TIPS_HEAD</head><body><main>
<p><a href='../../../../index.html'>Reports home</a> · <a href='../../../../outputs/final_report/report.html'>Final report</a></p>
<h1>Key findings, assumptions and limitations</h1>
<div class='intro'><p><b>$(fits_used_text())</b></p><p>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")) from the outputs of the latest pipeline run (zoo fits of $(Libc.strftime("%Y-%m-%d %H:%M", mtime(joinpath(OUTD, "zoo_fits", "zoo_state.jls"))))). Every number below is read from that run's files.</p></div>
$(glance())
$(transport())
$(starts())
$(validation())
$(control())
$(shape())
$(nn())
$(register())
$(limits())
</main></body></html>"""
mkpath(joinpath(OUTD, "findings")); mkpath(DOCS)
write(joinpath(OUTD, "findings", "report.html"), stage_names(html)); cp(joinpath(OUTD, "findings", "report.html"), joinpath(DOCS, "report.html"); force=true)
println("Wrote $(joinpath(DOCS, "report.html"))")
