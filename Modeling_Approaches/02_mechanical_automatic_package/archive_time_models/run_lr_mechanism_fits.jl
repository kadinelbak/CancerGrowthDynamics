# Low-resource mechanisms: fit 13 candidate descriptions of the 0.5% FBS monoculture effect,
# check each fit's long-term fate (does the population reach zero?), and pick the
# bootstrap-optimal model with the parsimony rule from bootstrap_tools.jl.
#
# Run with threads:  julia -t auto run_lr_mechanism_fits.jl      (B from env BOOT_B, default 200)

include(joinpath(@__DIR__, "run_new_data_staged_fits.jl"))
include(joinpath(@__DIR__, "bootstrap_tools.jl"))
include(joinpath(@__DIR__, "lr_mechanisms.jl"))

const LM_B = parse(Int, get(ENV, "BOOT_B", "200"))
const LM_DIR = joinpath(@__DIR__, "outputs", "lr_mechanisms")
const LM_DOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "lr_mechanisms")
const LM_HORIZON = 30.0
mkpath(LM_DIR)
println("\n=== LR mechanisms: B = $LM_B, $(Threads.nthreads()) threads ==="); flush(stdout)

lm = Dict{String,Any}()
for line in ("A2780cis", "A2780")
    t0 = time()
    rows = fit_lr_mechanisms(line)
    s = stage1b[line]
    g = (id="lr-$(lowercase(line))", datasets=[(x=s.x[.!s.bad], y=s.y[.!s.bad], state_index=1)], u0=[s.y[1]], cands=rows)
    res = bootstrap(g; B=LM_B)
    curves = Dict(c.name => dense(c.model, c.u0, c.params, LM_HORIZON; dt=0.1) for c in rows)
    lm[line] = (g=g, res=res, band=bands(g, res; tmax=LM_HORIZON), curves=curves)
    sel = res.cands[res.selected]
    @printf("  %-9s BIC-best %-24s | bootstrap %-24s k=%d picked %5.1f%% | reaches zero: %s  [%.0fs]\n", line, rows[1].name,
            sel.name, res.ks[res.selected], pct(res.pick, res.selected), sel.extinct ? "yes" : "no", time() - t0)
    for c in rows
        @printf("      %-24s k=%d BIC=%8.2f  x30=%8.1f x200=%8.2f  %s\n", c.name, length(c.params), c.bic, c.x30, c.x200, param_str(c))
    end
    flush(stdout)
end

# ── CSV ──────────────────────────────────────────────────────────────────────
csv_rows = NamedTuple[]
for (line, d) in lm
    r = d.res; best = minimum(c.bic for c in r.cands)
    for (j, c) in enumerate(r.cands)
        push!(csv_rows, (line=line, candidate=c.name, family=c.family, k=r.ks[j], bic=c.bic, delta_bic=c.bic - best,
            params=param_str(c), std_errors=join(["$(n) ± $(fmt(s))" for (n, s) in zip(c.pnames, c.se)], ", "),
            day30=c.x30, day60=c.x60, day200=c.x200, reaches_zero=c.extinct,
            pct_parsimony_pick=pct(r.pick, j), pct_lowest_bic=pct(r.winner, j), bootstrap_optimal=(j == r.selected)))
    end
end
CSV.write(joinpath(LM_DIR, "lr_mechanism_fits.csv"), DataFrame(csv_rows))

# ── report ───────────────────────────────────────────────────────────────────
const LM_COLORS = ["#0b6e99", "#d95f02", "#3b8b4f", "#8b5cf6", "#b8860b", "#c0392b"]
const LM_DASH = ["", "7,4", "2,3", "10,3,2,3", "4,4", "1,3"]

# obs with excluded points; curves: (label, color, dash, ts, ys); optional band (ts, lo, hi, color)
function lm_plot(title, x, y, bad, curves; band=nothing, xmax=LM_HORIZON, ymax=100.0, W=760, H=360)
    L, R, T, B = 64, 20, 30, 40 + 18 * cld(length(curves) + 1, 3)
    sx(v) = L + v / xmax * (W - L - R); sy(v) = H - B - clamp(v, 0, ymax) / ymax * (H - T - B)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(esc(title))'><rect width='100%' height='100%' fill='white'/>",
         "<text x='$L' y='18' font-size='14' font-weight='bold'>$(esc(title))</text>"]
    for v in range(0, ymax; length=5)
        push!(p, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#e5e9ec'/><text x='%d' y='%.1f' font-size='10' text-anchor='end'>%.0f</text>", L, W - R, sy(v), sy(v), L - 6, sy(v) + 3, v))
    end
    for v in 0:5:floor(Int, xmax)
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%d</text>", sx(v), H - B + 14, v))
    end
    xd = maximum(x)
    push!(p, @sprintf("<rect x='%.1f' y='%d' width='%.1f' height='%d' fill='#f3f5f7'/><text x='%.1f' y='%d' font-size='10' fill='#5b6b73'>beyond the data: extrapolation</text>", sx(xd), T, sx(xmax) - sx(xd), H - B - T, sx(xd) + 6, T + 12))
    push!(p, "<line x1='$L' y1='$(H-B)' x2='$(W-R)' y2='$(H-B)' stroke='#555'/><line x1='$L' y1='$T' x2='$L' y2='$(H-B)' stroke='#555'/>")
    push!(p, "<text x='$((L+W-R)÷2)' y='$(H-B+28)' font-size='11' text-anchor='middle'>Day</text><text transform='translate(14 $((T+H-B)÷2)) rotate(-90)' font-size='11' text-anchor='middle'>Mean Cells</text>")
    if band !== nothing
        up = [@sprintf("%.1f,%.1f", sx(t), sy(band.hi[i])) for (i, t) in enumerate(band.ts) if !isnan(band.hi[i])]
        dn = [@sprintf("%.1f,%.1f", sx(t), sy(band.lo[i])) for (i, t) in reverse(collect(enumerate(band.ts))) if !isnan(band.lo[i])]
        isempty(up) || push!(p, "<polygon points='$(join(vcat(up, dn), " "))' fill='$(band.color)' opacity='0.2'/>")
    end
    for c in curves
        pts = join([@sprintf("%.1f,%.1f", sx(t), sy(v)) for (t, v) in zip(c.ts, c.ys)], " ")
        push!(p, "<polyline points='$pts' fill='none' stroke='$(c.color)' stroke-width='2.2'$(isempty(c.dash) ? "" : " stroke-dasharray='$(c.dash)'")/>")
    end
    push!(p, join([@sprintf("<circle cx='%.1f' cy='%.1f' r='3.4' fill='#374151'/>", sx(a), sy(b)) for (a, b, z) in zip(x, y, bad) if !z]))
    push!(p, join([@sprintf("<path d='M%.1f %.1fL%.1f %.1fM%.1f %.1fL%.1f %.1f' stroke='#9ca3af' stroke-width='1.8'/>", sx(a)-4, sy(b)-4, sx(a)+4, sy(b)+4, sx(a)-4, sy(b)+4, sx(a)+4, sy(b)-4) for (a, b, z) in zip(x, y, bad) if z]))
    items = vcat([("data", "#374151", "", :dot)], [(c.label, c.color, c.dash, :line) for c in curves])
    for (i, (lab, col, dash, kind)) in enumerate(items)
        lx = L + ((i - 1) % 3) * 230; ly = H - B + 44 + ((i - 1) ÷ 3) * 18
        push!(p, kind == :dot ? "<circle cx='$(lx+12)' cy='$(ly-4)' r='4' fill='$col'/>" :
            "<line x1='$lx' x2='$(lx+24)' y1='$(ly-4)' y2='$(ly-4)' stroke='$col' stroke-width='2.2'$(isempty(dash) ? "" : " stroke-dasharray='$dash'")/>")
        push!(p, "<text x='$(lx+30)' y='$ly' font-size='11'>$(esc(lab))</text>")
    end
    push!(p, "</svg>"); join(p)
end

yesno(b) = b ? "<b>yes</b>" : "no"
sections = String[]; lm_summary = String[]
for line in ("A2780cis", "A2780")
    d = lm[line]; r = d.res; s = stage1b[line]; best = r.cands[1]; sel = r.cands[r.selected]
    ymax = ceil(maximum(s.y) * 1.2; sigdigits=2)
    push!(lm_summary, "<tr><td>$line</td><td>$(esc(best.label)) (k=$(length(best.params)))</td><td>$(yesno(best.extinct))</td>" *
        "<td><b>$(esc(sel.label))</b> (k=$(length(sel.params)))</td>" * @sprintf("<td class='num'>%.0f%%</td><td class='num'>%+.1f</td>", pct(r.pick, r.selected), sel.bic - best.bic) *
        "<td>$(yesno(sel.extinct))</td></tr>")
    trs = String[]
    for (j, c) in enumerate(r.cands)
        push!(trs, "<tr$(j == r.selected ? " class='win'" : "")><td>$(esc(c.label))</td><td>$(esc(c.family))</td><td class='num'>$(r.ks[j])</td>" *
            @sprintf("<td class='num'>%.2f</td><td class='num'>%.2f</td><td class='num'>%.0f%%</td><td class='num'>%.0f%%</td><td class='num'>%.0f</td><td class='num'>%.1f</td>",
                     c.bic, c.bic - best.bic, pct(r.pick, j), pct(r.winner, j), c.x30, c.x200) * "<td>$(yesno(c.extinct))</td></tr>")
    end
    top = r.cands[1:min(5, end)]
    r.selected > 5 && (top = vcat(top, [sel]))
    curves = [(label=c.label, color=LM_COLORS[mod1(i, 6)], dash=LM_DASH[mod1(i, 6)], ts=d.curves[c.name][1], ys=[u[1] for u in d.curves[c.name][2].u]) for (i, c) in enumerate(top)]
    selcurve = [(label="$(sel.label) (bootstrap-optimal)", color="#176b87", dash="", ts=d.band.ts, ys=d.band.mid[:, 1])]
    ptrs = String[]
    for (i, n) in enumerate(sel.pnames)
        l, m, h = qs(r.pars[r.selected][:, i])
        push!(ptrs, "<tr><td>$(esc(n))</td><td class='num'>$(shown(sel.params[i]))</td><td class='num'>$(isnan(sel.se[i]) ? "—" : "± " * fmt(1.96 * sel.se[i]))</td><td class='num'>$(fmt(l)) to $(fmt(h))</td><td>$(removable(l, h, sel.bounds[i]) ? "<b>interval reaches 0</b>" : "")</td></tr>")
    end
    push!(sections, "<section class='panel' id='$(lowercase(line))'><h2>$line at 0.5% FBS</h2>" *
        "<p class='muted'>Sheet <code>$(esc(s.sheet))</code>; r and K from the Stage 1 logistic fit. $(count(s.bad)) point(s) excluded as transient drops (grey ×).</p>" *
        "<div class='tscroll'><table><tr><th>Candidate</th><th>Family</th><th>k</th><th>BIC</th><th>ΔBIC</th><th>Bootstrap pick</th><th>Lowest BIC</th><th>Cells at day 30</th><th>Cells at day 200</th><th>Reaches zero?</th></tr>$(join(trs))</table></div>" *
        bic_svg(r.cands) *
        lm_plot("$line — five lowest-BIC models, extended to day $(Int(LM_HORIZON))", s.x, s.y, s.bad, curves; ymax=ymax) *
        "<h3>Bootstrap-optimal: $(esc(sel.label))</h3><div class='eq'>\\[$(sel.tex)\\]</div>" *
        "<p class='muted'>Picked in $(round(pct(r.pick, r.selected); digits=1))% of $LM_B replicates (lowest BIC in $(round(pct(r.winner, r.selected); digits=1))%). Cells at day 30: $(round(sel.x30; digits=1)); day 200: $(round(sel.x200; digits=2)).</p>" *
        "<div class='tscroll'><table class='params'><tr><th>Parameter</th><th>Estimate</th><th>95% (asymptotic)</th><th>95% bootstrap interval</th><th></th></tr>$(join(ptrs))</table></div>" *
        lm_plot("$line — bootstrap-optimal model with 95% bootstrap band", s.x, s.y, s.bad, selcurve;
                band=(ts=d.band.ts, lo=d.band.lo[:, 1], hi=d.band.hi[:, 1], color="#176b87"), ymax=ymax) *
        "</section>")
end

specs = lr_mechanism_specs("X", 1.0, 1.0, 1.0)
nmech = length(lm["A2780cis"].res.cands)
extra_eqs = [("Accumulated damage", "Accumulated low-resource damage, 4 steps", 2, "\\frac{dX}{dt}=r_X X\\left(1-\\frac{X}{K_X}\\right)-d^{LR}\\,u_R\\,E_4(\\lambda t)\\,X"),
             ("Accumulated damage", "Immediate loss + accumulated damage, 4 steps", 3, "\\frac{dX}{dt}=r_X X\\left(1-\\frac{X}{K_X}\\right)-\\left(d^{LR}+d^{acc}E_4(\\lambda t)\\right)u_R X"),
             ("Accumulated damage", "Immediate loss + accumulated damage, 8 steps", 3, "\\frac{dX}{dt}=r_X X\\left(1-\\frac{X}{K_X}\\right)-\\left(d^{LR}+d^{acc}E_8(\\lambda t)\\right)u_R X")]
eq_rows = join(["<tr><td>$(esc(c.family))</td><td>$(esc(c.label))</td><td class='num'>$(length(c.pnames))</td><td><div class='eq'>\\[$(c.tex)\\]</div></td></tr>" for c in specs]) *
    join(["<tr><td>$(esc(f))</td><td>$(esc(l))</td><td class='num'>$k</td><td><div class='eq'>\\[$t\\]</div></td></tr>" for (f, l, k, t) in extra_eqs])

lnav = "<nav aria-label='New data report navigation' class='topnav'><a href='../../../../index.html' class='nav nav-home' aria-label='Reports home' title='Reports home'>$home_svg</a>" *
    nav_link("Summary", "../../../../outputs/results/summary.html", false) *
    nav_link("New data", "../../../../outputs/new_data_report/report.html", false) *
    nav_link("Candidate models", "../reports/a2780_low_resource_modeling_plan.html", false) *
    nav_link("New data fitting", "../new_data_fits/report.html", false) *
    nav_link("Bootstrapping", "../bootstrap_fits/report.html", false) *
    nav_link("LR mechanisms", "report.html", true) *
    nav_link("Low-resource fitting", "../low_resource_staged_fits/report.html", false) * "</nav>"

lhtml = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>LR mechanisms</title>
<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}};</script><script defer src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>
<style>body{font:14px/1.5 system-ui,sans-serif;margin:0;background:#f4f7f8;color:#17212b}main{max-width:1250px;margin:auto;padding:24px}h1{color:#124b6e;margin:6px 0}h2{color:#124b6e;margin:0 0 6px}
.topnav{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:0 0 16px}.nav{display:inline-flex;align-items:center;min-height:34px;padding:0 11px;border-radius:8px;color:#176b87;background:#eef5f6;text-decoration:none;font:650 14px system-ui,sans-serif;white-space:nowrap}.nav-home{padding:6px 9px}.nav-active{color:#fff;background:#176b87}
.intro,.panel{background:#fff;border:1px solid #d9e2e6;border-radius:12px;padding:16px 18px;margin:14px 0;min-width:0;overflow-wrap:anywhere}.intro{border-left:6px solid #176b87}.note{background:#fff6d9;border-radius:8px;padding:10px 14px}
.eq{overflow-x:auto;overflow-y:hidden;background:#fbfdfd;border:1px solid #e4ecef;border-radius:8px;padding:4px 12px;margin:6px 0}svg{display:block;width:100%;height:auto;margin:10px 0}svg.bic{max-width:900px}.tscroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}td,th{border:1px solid #d3dde2;padding:6px 8px;text-align:left;vertical-align:top}th{background:#e9f3f8}td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}tr.win td{background:#eef8f1;font-weight:600}.muted{color:#5b6b73}h3{margin:16px 0 4px}
code{overflow-wrap:anywhere}mjx-container:not([display]){display:inline-block;max-width:100%;overflow-x:auto;vertical-align:middle}</style></head><body><main>
$lnav
<h1>Low-resource mechanisms</h1>
<div class='intro'><p>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")). At 0.5% FBS, A2780cis rises quickly, holds for a few days and then dies off to near zero. The earlier low-resource candidates all settle at a positive level, so they cannot reach zero. This report fits <b>$nmech candidate mechanisms</b> to each low-resource monoculture (30k sheet) and asks which describe the data <i>and</i> send the population to zero.</p>
<ul><li>As in the New data fitting report, the A2780cis day-0 value (the nominal seeding count, 5× below day 1) is dropped and day 1 is treated as day 0.</li>
<li>r and K are always the Stage 1 logistic values; only the low-resource terms are fitted. Hidden states (resource R, by-product W, quiescent Q, damaged D) start at the values shown in the equations.</li>
<li>Each mechanism is fitted from several starting points. The <b>bootstrap-optimal</b> model uses the same rule as the <a href='../bootstrap_fits/report.html'>Bootstrapping</a> report: $LM_B residual block-bootstrap replicates, and a simpler model is preferred when it is no more than $(Int(PER_PARAM_BIC)) BIC worse per removed parameter.</li>
<li><b>Reaches zero?</b> means the fitted model has fewer than 1 cell by day 200. Plots run to day $(Int(LM_HORIZON)); the shaded region is extrapolation beyond the data.</li></ul></div>
<h2>Summary</h2>
<div class='tscroll'><table><tr><th>Line</th><th>Lowest BIC</th><th>Reaches zero?</th><th>Bootstrap-optimal (parsimony)</th><th>Picked in</th><th>BIC cost</th><th>Reaches zero?</th></tr>$(join(lm_summary))</table></div>
$(join(sections, "\n"))
<section class='panel'><h2>The $nmech mechanisms</h2><p class='muted'>X is the observed population (N for A2780, C for A2780cis). ★ = the model carried forward on the Candidate models page. \\(E_n(z)=1-e^{-z}\\sum_{k=0}^{n-1}z^k/k!\\) is the Erlang onset: damage passes through n steps before it kills. The accumulated-damage models come from the New data fitting report.</p>
<div class='tscroll'><table><tr><th>Family</th><th>Model</th><th>k</th><th>Equation</th></tr>$eq_rows</table></div></section>
<p class='muted'>Outputs: <code>lr_mechanism_fits.csv</code> in <code>Modeling_Approaches/02_mechanical_automatic_package/outputs/lr_mechanisms/</code>. Generated by <code>run_lr_mechanism_fits.jl</code>.</p>
</main></body></html>"""

write(joinpath(LM_DIR, "report.html"), lhtml)
mkpath(LM_DOCS)
cp(joinpath(LM_DIR, "report.html"), joinpath(LM_DOCS, "report.html"); force=true)
println("\nWrote $(joinpath(LM_DIR, "report.html"))\n  and $(joinpath(LM_DOCS, "report.html"))")
