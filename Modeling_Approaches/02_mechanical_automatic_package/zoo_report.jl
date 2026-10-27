# HTML report for run_zoo_fits.jl: the model catalogue, then for every fit group the ranked
# candidates, the optimal model's equations and parameters (with bootstrap intervals) and its
# fit with a 95% bootstrap band. Written to outputs/zoo_fits/report.html and copied into docs/.

hesc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "'" => "&#39;")
include(joinpath(@__DIR__, "info_tips.jl"))   # tip(), INFO_TIPS_HEAD: the "i" buttons
thtip(label, text, better; warn=nothing) = "<th>$label$(tip(text, better; warn=warn))</th>"
const ZT_BIC = "Bayesian information criterion, n·ln(SSE/n) + k·ln(n): the fit error penalised for the number of parameters. Lower is better; differences under about 2 are negligible, above 10 decisive."
const ZT_BAND = "Points: measured mean cells (grey × = excluded outlier; in cocultures solid/filled = untreated Ce0, dashed/open = treated Ce1). Line: the carried model fitted to the data. Shaded band: where 95% of the curves refitted to the bootstrap data sets lie, i.e. the uncertainty of the fitted curve."

function zoo_band_svg(title, series, curve_ts, curve, band, tmax, ymax; W=760, H=330)
    L, R, T, B = 64, 20, 30, 62
    sx(v) = L + v / tmax * (W - L - R); sy(v) = H - B - clamp(v, 0, ymax) / ymax * (H - T - B)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(hesc(title))'><rect width='100%' height='100%' fill='white'/>",
         "<text x='$L' y='18' font-size='14' font-weight='bold'>$(hesc(title))</text>"]
    for v in range(0, ymax; length=5)
        push!(p, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#e5e9ec'/><text x='%d' y='%.1f' font-size='10' text-anchor='end'>%.0f</text>", L, W - R, sy(v), sy(v), L - 6, sy(v) + 3, v))
    end
    for v in 0:2:floor(Int, tmax)
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%d</text>", sx(v), H - B + 14, v))
    end
    push!(p, "<line x1='$L' y1='$(H-B)' x2='$(W-R)' y2='$(H-B)' stroke='#555'/><line x1='$L' y1='$T' x2='$L' y2='$(H-B)' stroke='#555'/>")
    push!(p, "<text x='$((L+W-R)÷2)' y='$(H-B+30)' font-size='11' text-anchor='middle'>Day</text><text transform='translate(14 $((T+H-B)÷2)) rotate(-90)' font-size='11' text-anchor='middle'>Mean Cells</text>")
    for s in series
        if band !== nothing
            up = [@sprintf("%.1f,%.1f", sx(t), sy(band.hi[i, s.state])) for (i, t) in enumerate(band.ts) if !isnan(band.hi[i, s.state])]
            dn = [@sprintf("%.1f,%.1f", sx(t), sy(band.lo[i, s.state])) for (i, t) in reverse(collect(enumerate(band.ts))) if !isnan(band.lo[i, s.state])]
            isempty(up) || push!(p, "<polygon points='$(join(vcat(up, dn), " "))' fill='$(s.color)' opacity='0.18'/>")
        end
        pts = join([@sprintf("%.1f,%.1f", sx(t), sy(curve[i][s.state])) for (i, t) in enumerate(curve_ts)], " ")
        push!(p, "<polyline points='$pts' fill='none' stroke='$(s.color)' stroke-width='2'$(s.dashed ? " stroke-dasharray='7,4'" : "")/>")
        fill = s.filled ? s.color : "white"
        push!(p, join([@sprintf("<circle cx='%.1f' cy='%.1f' r='3.2' fill='%s' stroke='%s' stroke-width='1.4'/>", sx(x), sy(y), fill, s.color) for (x, y) in zip(s.x, s.y)]))
        for (x, y) in zip(s.xbad, s.ybad)
            push!(p, @sprintf("<path d='M%.1f %.1fL%.1f %.1fM%.1f %.1fL%.1f %.1f' stroke='#9ca3af' stroke-width='1.8'/>", sx(x)-4, sy(y)-4, sx(x)+4, sy(y)+4, sx(x)-4, sy(y)+4, sx(x)+4, sy(y)-4))
        end
    end
    for (i, s) in enumerate(unique(s -> s.name, series))
        lx = L + (i - 1) * 150; fill = s.filled ? s.color : "white"
        push!(p, "<circle cx='$(lx+4)' cy='$(H-16)' r='4' fill='$fill' stroke='$(s.color)' stroke-width='1.4'/><text x='$(lx+12)' y='$(H-12)' font-size='11' fill='$(s.color)'>$(hesc(s.name))</text>")
    end
    push!(p, "</svg>"); join(p)
end

const ZOO_CATALOGUE = """
<table><tr><th>Delay / accumulation (how the input w = D or u<sub>R</sub> becomes the effect E)</th><th>Equations</th></tr>
<tr><td>Immediate</td><td>\\(E = w\\)</td></tr>
<tr><td>Effect compartment</td><td>\\(\\frac{dA}{dt}=k\\,(w-A),\\quad E=A\\)</td></tr>
<tr><td>2 or 3 transit compartments (gamma-distributed delay)</td><td>\\(\\frac{dA_1}{dt}=k(w-A_1),\\ \\frac{dA_2}{dt}=k(A_1-A_2),\\ \\frac{dA_3}{dt}=k(A_2-A_3),\\quad E=A_3\\)</td></tr>
<tr><td>Second order (double derivative)</td><td>\\(\\frac{d^2A}{dt^2}=\\omega^2(w-A)-2\\zeta\\omega\\frac{dA}{dt},\\quad E=\\max(A,0)\\)</td></tr>
<tr><td>Accumulation (never cleared)</td><td>\\(\\frac{dA}{dt}=w,\\quad E=A\\)</td></tr>
<tr><th>What the effect does to the cells (X observed, V viable)</th><th></th></tr>
<tr><td>Linear kill / loss</td><td>\\(-\\beta\\,E\\,V\\) (optionally \\(-(d^{LR}w+d^{acc}E)V\\))</td></tr>
<tr><td>Hill kill</td><td>\\(-\\beta\\,\\frac{E^h}{EC_{50}^h+E^h}\\,V\\)</td></tr>
<tr><td>Damaged cells</td><td>\\(\\frac{dX_D}{dt}=\\beta E V-\\delta X_D,\\ V=X-X_D\\), deaths \\(\\delta X_D\\) (damaged cells are counted and crowd, but do not divide)</td></tr>
<tr><td>Two-stage damage</td><td>\\(X_{D1}\\to X_{D2}\\to\\) death at rate \\(\\tau\\)</td></tr>
<tr><td>Lower capacity / slower growth</td><td>\\(K_{eff}=K/(1+cE)\\) or \\(r_{eff}=r/(1+cE)\\)</td></tr>
<tr><td>Resource depletion / toxic by-product (low resource)</td><td>\\(\\frac{dR}{dt}=-cER\\frac{X}{1000}\\), growth \\(\\times R\\); \\(\\frac{dW}{dt}=E\\frac{X}{1000}-\\kappa W\\), loss \\(\\gamma W\\)</td></tr>
</table>"""

# the paired-line rules (Stage 2 shared structure, Stage 1b A/B test), shown in both lines' sections
function joint_html(g)
    (isdefined(Main, :JOINT_INFO) && g.kind in (:s1b, :s2)) || return ""
    info = JOINT_INFO[g.kind == :s2 ? "s2" : "s1b"]
    top = join(["<tr$(x.name == (g.kind == :s2 ? info.pick.name : info.B.name) ? " class='win'" : "")><td>$(hesc(x.label))</td>" *
                @sprintf("<td class='num'>%.2f</td><td class='num'>%.2f</td><td class='num'>%.2f</td><td class='num'>%d</td></tr>", x.bicN, x.bicC, x.bic, x.k)
                for x in info.ranking[1:min(6, end)]])
    tbl = "<div class='tscroll'><table><tr>$(thtip("Shared structure", "A mechanism structure used by both lines, each with its own rates.", :none))$(thtip("BIC A2780", "BIC of this structure fitted to the A2780 data.", :lower))$(thtip("BIC A2780cis", "BIC of this structure fitted to the A2780cis data.", :lower))$(thtip("Joint BIC", "Sum of the two lines' BICs. The selection rule (fewest parameters within 2 BIC of the lowest) is applied to this.", :lower))$(thtip("Joint k", "Free parameters of both lines together.", "Lower is better (fewer parameters), as long as the fit stays within 2 BIC."))</tr>$top</table></div>"
    if g.kind == :s2
        ex = [e for e in excluded_sources() if e.stage == "s2"]
        exn = isempty(ex) ? "" : "<div class='note'><b>Left out of this fit</b> (user, 2026-10-27; stages.toml <code>exclude</code>): " *
              join(["$(hesc(e.line)) $(hesc(e.file))" for e in ex], ", ") * ". The A2780cis 1 µM monoculture shows no kill by day 14 (it ends above untreated), so A2780cis is fitted on the other doses. The fit with every dose is kept as a sensitivity variant (<code>ZOO_NO_EXCLUDE=1</code>).</div>"
        return exn * "<div class='note'><b>Shared mechanism (both lines).</b> A2780cis is derived from A2780, and cisplatin kills both through DNA damage followed by delayed death; resistance changes the rates (less uptake, glutathione, faster repair, more damage tolerance), not the mechanism. So both lines use one structure with their own rates: the rule is applied to the joint BIC (sum of the two lines' BICs, each fitted with its own error variance) and joint k. Chosen: <b>$(hesc(info.pick.label))</b> (joint BIC $(fmt(info.pick.bic)), k = $(info.pick.k)).</div>$tbl"
    end
    A, B = info.A, info.B
    "<div class='note'><b>Shared low-resource mechanism: tested, not assumed.</b> A = each line's own pick ($(hesc(A.labelN)) for A2780, $(hesc(A.labelC)) for A2780cis): summed BIC $(fmt(A.bic)), k = $(A.k). B = best shared structure ($(hesc(B.label))): joint BIC $(fmt(B.bic)), k = $(B.k), ΔBIC(B − A) = $(@sprintf("%+.2f", info.dbic)). B is used only if within $(Int(RULE_BIC)) BIC of A with no more parameters, so <b>$(info.useB ? "the shared structure (B)" : "each line keeps its own model (A)")</b>.</div>$tbl"
end

function write_zoo_report(groups, results)
    home_svg = "<svg width='22' height='22' viewBox='0 0 24 24' aria-hidden='true'><path fill='currentColor' d='M12 3.1 2.8 10.6v10.3h6.1v-6.2h6.2v6.2h6.1V10.6L12 3.1Zm7.1 15.8h-2v-6.2H6.9v6.2h-2V11.6L12 5.8l7.1 5.8v7.3Z'/></svg>"
    nl(label, href, on) = "<a href='$href' class='nav$(on ? " nav-active" : "")'>$label</a>"
    nav = "<nav class='topnav' aria-label='Report navigation'><a href='../../../../index.html' class='nav nav-home' aria-label='Reports home' title='Reports home'>$home_svg</a>" *
          nl("Summary", "../../../../outputs/results/summary.html", false) * nl("New data", "../../../../outputs/new_data_report/report.html", false) *
          nl("Model zoo fits", "report.html", true) * nl("Archive: time-based models", "../../../../outputs/archive_time/index.html", false) * "</nav>"
    summary = String[]; sections = String[]
    for g in groups
        r = results[g.id].res; band = results[g.id].band; c = g.cands[1]; hasb = size(r.bic, 1) > 0
        push!(summary, "<tr><td><a href='#$(g.id)'>$(hesc(g.title))</a></td><td><b>$(hesc(c.label))</b>" *
              join(["<div class='eq glance'>\\[$t\\]</div>" for t in zoo_tex(g, c)]) * "</td><td>$(hesc(c.family))</td><td class='num'>$(length(c.params))</td>" *
              @sprintf("<td class='num'>%.1f</td><td class='num'>%s</td></tr>", c.bic, hasb ? @sprintf("%.0f%%", pct(r.pick, 1)) : "—"))
        trs = join([begin
                        low = hasb ? @sprintf("%.0f%% / %.0f%%", pct(r.pick, j), pct(r.winner, j)) : "—"
                        "<tr$(j == 1 ? " class='win'" : "")><td><b>$(hesc(x.label))</b>" * join(["<div class='eq glance'>\\[$t\\]</div>" for t in zoo_tex(g, x)]) *
                        "</td><td>$(hesc(x.family))</td><td class='num'>$(length(x.params))</td>" *
                        @sprintf("<td class='num'>%.2f</td><td class='num'>%+.2f</td><td class='num'>%s</td></tr>", x.bic, x.bic - c.bic, low)
                    end for (j, x) in enumerate(g.cands)])
        ptrs = join([begin
                         l, m, h = hasb ? qs(r.pars[1][:, i]) : (NaN, NaN, NaN)
                         "<tr><td>\\($(ptex(n))\\)</td><td class='num'>$(fmt(c.params[i]))</td><td class='num'>$(fmt(1.96c.se[i]))</td><td class='num'>$(hasb ? "$(fmt(l)) to $(fmt(h))" : "—")</td></tr>"
                     end for (i, n) in enumerate(c.pnames)])
        inh = zoo_inherited(g, c)
        itrs = join(["<tr class='inh'><td>\\($(x.sym)\\)</td><td class='num'>$(fmt(x.value))</td><td colspan='2' class='muted'>$(hesc(x.note))</td></tr>" for x in inh])
        tmax = g.kind == :s1b ? 20.0 : XMAX
        ts, sol = dense(c.model, c.u0, c.params, tmax; dt=0.1)
        series = zoo_series(g); panels = zoo_panels(g)
        ymax = g.kind == :s1b ? ceil(maximum(maximum(s.y) for s in series) * 1.2; sigdigits=2) : YMAX
        plots = join(["<div class='tipwrap'>$(tip(ZT_BAND, "Better: the line runs through the points and the shaded band is narrow."))" * zoo_band_svg(isempty(pt) ? "$(g.title): optimal fit" * (hasb ? " with 95% bootstrap band" : "") : pt,
                                   [s for s in series if s.panel == pn], ts, sol.u, band, tmax, ymax) * "</div>" for (pn, pt) in enumerate(panels)])
        extra = g.kind == :s1b ? "<p class='muted'>Long-term check of the optimal model$(tip("The carried model simulated to day 30 and day 200, to check it behaves sensibly, e.g. whether low resource alone would eventually kill the culture.", :none; warn="Far beyond the 14 days of data: a qualitative check of the equations, not a prediction.")): $(fmt(c.x30)) cells at day 30, $(fmt(c.x200)) at day 200 ($(c.extinct ? "dies out" : "does not die out")).</p>" : ""
        push!(sections, "<section class='panel' id='$(g.id)'><h3>$(hesc(g.title))</h3><p class='muted'>$(hesc(ZOO_GROUP_TEXT[g.kind]))</p>" *
              "<div class='tscroll'><table><tr>" * thtip("Candidate", "Every model fitted at this stage, with its equations. The highlighted row is the model carried forward.", :none) *
              thtip("Family", "The delay mechanism the model uses: immediate, effect compartment, transit compartments, second order, accumulation, or damaged cells.", :none) *
              thtip("k", "Free parameters fitted at this stage; values carried over from earlier stages are not counted.", "Lower is better (fewer parameters), as long as the fit stays within 2 BIC.") * thtip("BIC", ZT_BIC, :lower) *
              thtip("ΔBIC vs optimal", "BIC minus the BIC of the carried model. Negative = fits slightly better but needs more parameters (within 2 BIC the rule prefers fewer); positive = fits worse.", "Lower means a better fit, but within 2 of 0 the carried model wins by having fewer parameters; above +2 is clearly worse.") *
              thtip("Bootstrap: chosen / lowest BIC", "First number: share of the bootstrap data sets in which the selection rule picks this model. Second: share in which it has the lowest BIC outright.", :higher) * "</tr>$trs</table></div>" *
              joint_html(g) * "<h4>Optimal model: $(hesc(c.label))$(tip("The carried model's equations. N, C: sensitive (A2780) and resistant (A2780cis) counts; V: viable, non-damaged cells; A or S: delay compartments driven by the dose D or the low-resource level u_R. Parameter values are in the table below.", :none))</h4>" * join(["<div class='eq'>\\[$t\\]</div>" for t in zoo_tex(g, c)]) *
              "<div class='tscroll'><table class='params'><tr>" * thtip("Parameter", "Parameters fitted at this stage. Rows under 'Carried over / fixed' come from earlier stages and are held at their values.", :none) *
              thtip("Estimate", "Value fitted to the data.", :none) *
              thtip("±95% (asymptotic)", "1.96 × the standard error from the curvature of the fit.", "Narrower is better: the parameter is more precisely determined."; warn="It assumes independent errors, so it is usually too narrow; prefer the bootstrap interval.") *
              thtip("95% bootstrap interval", "95% of the values refitted to the block-bootstrap data sets fall in this range. Resampling 3-day blocks keeps day-to-day correlation.", "Narrower is better: the parameter is more precisely determined.") * "</tr>$ptrs" *
              (isempty(itrs) ? "" : "<tr><th colspan='4'>Carried over / fixed</th></tr>$itrs") * "</table></div>$extra$plots</section>")
    end
    html = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Model zoo fits</title>
<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}};</script><script defer src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>
<style>body{font:14px/1.5 system-ui,sans-serif;margin:0;background:#f4f7f8;color:#17212b}main{max-width:1250px;margin:auto;padding:24px}h1{color:#124b6e;margin:6px 0}h2{color:#124b6e;margin-top:30px}
.topnav{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:0 0 16px}.nav{display:inline-flex;align-items:center;min-height:34px;padding:0 11px;border-radius:8px;color:#176b87;background:#eef5f6;text-decoration:none;font:650 14px system-ui,sans-serif;white-space:nowrap}.nav-home{padding:6px 9px}.nav-active{color:#fff;background:#176b87}
.intro,.panel{background:#fff;border:1px solid #d9e2e6;border-radius:12px;padding:16px 18px;margin:14px 0;min-width:0;overflow-wrap:anywhere}.intro{border-left:6px solid #176b87}.note{background:#fff6d9;border-radius:8px;padding:10px 14px}
.eq{overflow-x:auto;overflow-y:hidden;background:#fbfdfd;border:1px solid #e4ecef;border-radius:8px;padding:4px 12px;margin:8px 0}.eq.glance{font-size:13px;margin:4px 0;padding:0 8px}svg{display:block;width:100%;height:auto;margin:10px 0}.tscroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}td,th{border:1px solid #d3dde2;padding:6px 8px;text-align:left;vertical-align:top}th{background:#e9f3f8}td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}tr.win td{background:#eef8f1;font-weight:600}tr.inh td{background:#f7f9fa}.muted{color:#5b6b73}h3{margin:0 0 6px}h4{margin:16px 0 4px}
mjx-container:not([display]){display:inline-block;max-width:100%;overflow-x:auto;vertical-align:middle}.tipwrap{position:relative}.tipwrap>.ii{position:absolute;top:8px;right:8px;z-index:2}</style>$INFO_TIPS_HEAD</head><body><main>
$nav
<h1>Model zoo fits — no time-based terms</h1>
<div class='intro'><p>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")). <b>No model here uses time t directly.</b> Every right-hand side depends only on the cell states and the two inputs, the cisplatin dose \\(D\\) and the low-resource level \\(u_R\\). Delays and accumulation come from extra states driven by those inputs, so any of these models can be driven by a changing dose schedule \\(D(t)\\) for optimal control. For a constant dose switched on at day 0, an effect compartment gives \\(D(1-e^{-kt})\\), the same curve as the archived time-delay term, but it stays correct when the dose changes.</p>
<ul><li><b>Optimal model</b> at each stage: among the candidates within $(Int(RULE_BIC)) BIC of the lowest, the one with the fewest free parameters (lower BIC breaks ties). It is carried into the later stages, so a negative ΔBIC in a table marks a model that fits slightly better but needs more terms.</li>
<li><b>Stage 2</b>: both lines share one cisplatin mechanism structure with their own rates (the rule is applied to the joint BIC of the two lines). <b>Stage 1b</b>: a shared low-resource structure is tested against each line's own pick and used only if it is within $(Int(RULE_BIC)) BIC with no more parameters.</li>
<li><b>At most $MAX_COMPARTMENTS compartments per cell line</b> (delay states plus damaged-cell states), so every model stays small enough for optimal control.</li>
<li><b>Bootstrap</b>: every candidate is refitted to $BOOT_B residual block-bootstrap replicates (blocks of $BLOCK days, generated from the carried model; starting values held at their fitted values) for 95% parameter intervals, the band around the fit, and how often the same rule picks each model (and how often it has the lowest BIC) when the data are resampled.</li>
<li><b>Starting values</b>: day 0 is the nominal seeding (100 cells at 30k, split by the seeding ratio in cocultures), used as the initial condition, not as a data point. Where the first measured count is more than 3-fold above or below it, that series' starting value is estimated (one extra parameter per series, counted in the BIC).</li>
<li>Values carried over from earlier stages are held at their point estimates. The older time-based models are kept in the <a href='../../../../outputs/archive_time/index.html'>archive</a>.</li></ul>
<h3>Model zoo $(tip("The building blocks of every candidate: a delay (how the dose or low-resource level becomes an effect) combined with an action (what the effect does to the cells). Every combination with at most 3 compartments per cell line is fitted.", :none))</h3>$ZOO_CATALOGUE</div>
<h2>Optimal model at each stage</h2>
<div class='tscroll'><table><tr>$(thtip("Fit", "The stage and data set.", :none))$(thtip("Optimal model and its equations", "The carried model: fewest free parameters within 2 BIC of the lowest. N, C: sensitive and resistant counts; D: cisplatin dose (µM); u_R: low-resource level; A, S: delay compartments; K: carrying capacity; α: crowding coefficients.", :none))$(thtip("Family", "The delay mechanism of the model.", :none))$(thtip("k", "Free parameters fitted at this stage.", "Lower is better (fewer parameters), as long as the fit stays within 2 BIC."))$(thtip("BIC", ZT_BIC, :lower))$(thtip("Chosen in bootstrap", "Share of the bootstrap data sets in which the selection rule picks this same model. Above about 70% = a robust choice.", :higher))</tr>$(join(summary))</table></div>
<h2>Stages</h2>
$(join(sections, "\n"))
<p class='muted'>Outputs: <code>zoo_model_selection.csv</code>, <code>zoo_parameters.csv</code> in <code>Modeling_Approaches/02_mechanical_automatic_package/outputs/zoo_fits/</code>. Generated by <code>run_zoo_fits.jl</code> (models in <code>zoo_models.jl</code>).</p>
</main></body></html>"""
    write(joinpath(ZOO_DIR, "report.html"), stage_names(html))
    if !isdefined(Main, :ZOO_PUBLISH) || ZOO_PUBLISH          # sensitivity runs (ZOO_OUT) are not published
        mkpath(ZOO_DOCS)
        cp(joinpath(ZOO_DIR, "report.html"), joinpath(ZOO_DOCS, "report.html"); force=true)
        println("Wrote $(joinpath(ZOO_DIR, "report.html"))\n  and $(joinpath(ZOO_DOCS, "report.html"))")
    else
        println("Wrote $(joinpath(ZOO_DIR, "report.html"))")
    end
end
