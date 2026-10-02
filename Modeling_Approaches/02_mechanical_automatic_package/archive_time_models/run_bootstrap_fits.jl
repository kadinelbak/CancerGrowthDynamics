# Bootstrap model selection for every stage of run_new_data_staged_fits.jl.
#
# Each dataset is regenerated B times as (fit of the most complex candidate) + (its residuals,
# resampled in circular blocks of 3 days within each trajectory; day 0 kept as observed). Every
# candidate is refitted to every replicate. Within a replicate the parsimony rule picks the
# candidate with the fewest free parameters whose BIC is no more than PER_PARAM_BIC per removed
# parameter above the lowest BIC. The bootstrap-optimal model is the one picked most often.
# Inherited parameters are held at their point estimates (upstream uncertainty is not propagated).
#
# Run with threads:  julia -t auto run_bootstrap_fits.jl      (B from env BOOT_B, default 200)

include(joinpath(@__DIR__, "run_new_data_staged_fits.jl"))

include(joinpath(@__DIR__, "bootstrap_tools.jl"))
include(joinpath(@__DIR__, "lr_mechanisms.jl"))

const BOOT_B = parse(Int, get(ENV, "BOOT_B", "200"))
const BOOT_DIR = joinpath(@__DIR__, "outputs", "bootstrap_fits")
const BOOT_DOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "bootstrap_fits")
mkpath(BOOT_DIR)
println("\n=== Bootstrap: B = $BOOT_B, $(Threads.nthreads()) threads ==="); flush(stdout)

# ── extra, pruned candidates (fitted to the original data here, not in the main report) ──

function make_treated_immediate(kind, r, K, doses)
    kind == :linear && return (du, u, p, t) -> begin
        for i in eachindex(doses)
            X = max(u[i], zero(p[1])); du[i] = r*X*(1 - X/K) - p[1]*doses[i]*X
        end
    end
    (du, u, p, t) -> begin
        for i in eachindex(doses)
            X = max(u[i], zero(p[1])); du[i] = r*X*(1 - X/K) - p[1]*hill(doses[i], p[2])*X
        end
    end
end

# ── fit groups ────────────────────────────────────────────────────────────────

groups = NamedTuple[]
for (line, _, _) in stage1_lines
    s = stage1[line]
    push!(groups, (id="s1-$(lowercase(line))", stage="1 Untreated monoculture", title="Stage 1 · $line untreated",
        datasets=[(x=s.x, y=s.y, state_index=1)], u0=[s.y[1]], cands=copy(s.rows), line=line, kind=:s1))
end
for (line, _, _) in stage1_lines
    s = stage1b[line]
    push!(groups, (id="s1b-$(lowercase(line))", stage="1b Low-resource monoculture", title="Stage 1b · $line low resource",
        datasets=[(x=s.x[.!s.bad], y=s.y[.!s.bad], state_index=1)], u0=[s.y[1]], cands=fit_lr_mechanisms(line), line=line, kind=:s1b))
end
for (line, _, _) in stage1_lines
    s = stage2[line]; s1 = stage1[line]; X = s.sym
    datasets = [(x=x, y=y, state_index=i) for (i, (x, y)) in enumerate(s.data)]
    extra = fit_candidates([
        (name="immediate_hill_kill", label="Immediate Hill kill (no delay)", model=make_treated_immediate(:hill, s1.r, s1.K, s.doses),
         pnames=["β_$X", "h_$X"], p0=[0.3, 2.0], bounds=[(0.0, 10.0), (0.1, 10.0)]),
        (name="immediate_linear_kill", label="Immediate kill, no delay, no Hill", model=make_treated_immediate(:linear, s1.r, s1.K, s.doses),
         pnames=["β_$X"], p0=[0.1], bounds=[(0.0, 10.0)]),
    ], datasets, s.u0)
    push!(groups, (id="s2-$(lowercase(line))", stage="2 Treated monoculture", title="Stage 2 · $line treated (3 doses jointly)",
        datasets=datasets, u0=s.u0, cands=sort(vcat(copy(s.rows), extra); by=r -> r.bic), line=line, kind=:s2))
end
push!(groups, (id="s3", stage="3 Ce0 + Ce1 coculture", title="Stage 3 · Ce0 + Ce1 coculture (shared K)",
    datasets=co_datasets, u0=u0_s3, cands=copy(s3_rows), line="coculture", kind=:s3))
extra3b = fit_candidates([
    (name="lr_cc_crowd_no_loss", label="LR coculture: K + α_NC + α_CN, no low-resource loss",
     model=make_lr_cc(p -> (p[1], p[2], p[3], 0.0, 0.0), rN, rC), pnames=["K", "α_NC", "α_CN"], p0=[K3, 0.5, 2.0],
     bounds=[(0.5lrmax, 20lrmax), (0.0, 5.0), (0.0, 5.0)]),
    (name="lr_cc_alphaCN_only", label="LR coculture: K + α_CN only (α_NC = 0, no low-resource loss)",
     model=make_lr_cc(p -> (p[1], 0.0, p[2], 0.0, 0.0), rN, rC), pnames=["K", "α_CN"], p0=[K3, 2.0],
     bounds=[(0.5lrmax, 20lrmax), (0.0, 10.0)]),
], lr_cc_datasets, u0_3b)
push!(groups, (id="s3b", stage="3b Low-resource coculture", title="Stage 3b · low-resource coculture",
    datasets=lr_cc_datasets, u0=u0_3b, cands=sort(vcat(copy(s3b_rows), extra3b); by=r -> r.bic), line="LR coculture", kind=:s3b))

results = Dict{String,Any}()
for g in groups
    t0 = time()
    res = bootstrap(g; B=BOOT_B)
    results[g.id] = (res=res, band=bands(g, res))
    sel = res.cands[res.selected]
    @printf("  %-8s selected %-30s k=%d  picked in %5.1f%%  (lowest BIC in %5.1f%%)  [%.0fs]\n", g.id, sel.name, res.ks[res.selected],
            100 * count(==(res.selected), res.pick) / BOOT_B, 100 * count(==(res.selected), res.winner) / BOOT_B, time() - t0)
    flush(stdout)
end

# ── CSV outputs ───────────────────────────────────────────────────────────────

sel_rows = NamedTuple[]; par_rows = NamedTuple[]
for g in groups
    r = results[g.id].res; B = size(r.bic, 1)
    dbic = r.bic .- minimum(r.bic; dims=2)
    for (j, c) in enumerate(r.cands)
        l, m, h = qs(dbic[:, j])
        push!(sel_rows, (group=g.id, stage=g.stage, candidate=c.name, k=r.ks[j], bic_original=c.bic,
            dbic_original=c.bic - minimum(x.bic for x in r.cands), dbic_boot_median=m, dbic_boot_lo=l, dbic_boot_hi=h,
            pct_lowest_bic=pct(r.winner, j), pct_parsimony_pick=pct(r.pick, j),
            selected=(j == r.selected), bic_best_original=(j == 1), parsimony_pick_original=(j == r.orig_pick)))
    end
    c = r.cands[r.selected]
    for (i, n) in enumerate(c.pnames)
        l, m, h = qs(r.pars[r.selected][:, i])
        push!(par_rows, (group=g.id, candidate=c.name, param=n, estimate=c.params[i], boot_median=m, boot_lo=l, boot_hi=h,
            removable=removable(l, h, c.bounds[i])))
    end
end
CSV.write(joinpath(BOOT_DIR, "bootstrap_model_selection.csv"), DataFrame(sel_rows))
CSV.write(joinpath(BOOT_DIR, "bootstrap_parameters.csv"), DataFrame(par_rows))

# ── report ────────────────────────────────────────────────────────────────────

function bs_tex(g, c)
    p = c.params
    if g.kind == :s1
        return [tex_stage1(c.name, stage1[g.line].sym, p)]
    elseif g.kind == :s1b
        s1 = stage1[g.line]
        c.name in ("lr_subtractive_loss", "lr_rho_modifiers", "lr_separate_r") || return [c.tex]
        return [tex_stage1b(c.name, s1.sym, p, s1.r, s1.K)]
    elseif g.kind == :s2
        s1 = stage1[g.line]; X = s1.sym
        grow = "$(tn(s1.r))\\,$X\\left(1-\\frac{$X}{$(tn(s1.K))}\\right)"
        c.name == "immediate_linear_kill" && return ["\\frac{d$X}{dt}=$grow-$(tn(p[1]))\\,D\\,$X"]
        c.name == "immediate_hill_kill" && return ["\\frac{d$X}{dt}=$grow-$(tn(p[1]))\\,\\frac{D^{$(tn(p[2]))}}{1+D^{$(tn(p[2]))}}\\,$X"]
        return [tex_stage2(c.name, X, p, s1.r, s1.K)]
    elseif g.kind == :s3
        return tex_stage3(c)
    else
        if c.name == "lr_cc_crowd_no_loss" || c.name == "lr_cc_alphaCN_only"
            K = p[1]; aNC = c.name == "lr_cc_crowd_no_loss" ? p[2] : 0.0; aCN = p[end]
            return ["\\frac{dN}{dt}=$(tn(rN))\\,N\\left(1-\\frac{N+$(tn(aNC))\\,C}{$(tn(K))}\\right)",
                    "\\frac{dC}{dt}=$(tn(rC))\\,C\\left(1-\\frac{C+$(tn(aCN))\\,N}{$(tn(K))}\\right)"]
        end
        return tex_stage3b(c)
    end
end

function freq_svg(r; W=760)
    nc = length(r.cands); H = 40 + 46nc
    p = ["<svg viewBox='0 0 $W $H' class='bic' role='img' aria-label='bootstrap selection frequencies'>",
         "<text x='8' y='16' font-size='13' font-weight='bold'>How often each candidate is chosen across $(size(r.bic, 1)) bootstrap replicates</text>",
         "<rect x='470' y='24' width='12' height='10' fill='#9fb7c3'/><text x='488' y='33' font-size='11'>lowest BIC</text>",
         "<rect x='580' y='24' width='12' height='10' fill='#176b87'/><text x='598' y='33' font-size='11'>parsimony rule</text>"]
    for (j, c) in enumerate(r.cands)
        y = 44 + (j - 1) * 46
        a, b = pct(r.winner, j), pct(r.pick, j)
        push!(p, @sprintf("<text x='8' y='%d' font-size='12'%s>%s (k=%d)</text>", y + 10, j == r.selected ? " font-weight='bold'" : "", esc(c.label), r.ks[j]))
        push!(p, @sprintf("<rect x='8' y='%d' width='%.1f' height='10' fill='#9fb7c3'/><text x='%.1f' y='%d' font-size='11'>%.0f%%</text>", y + 15, 4a, 14 + 4a, y + 24, a))
        push!(p, @sprintf("<rect x='8' y='%d' width='%.1f' height='10' fill='#176b87'/><text x='%.1f' y='%d' font-size='11'>%.0f%%</text>", y + 27, 4b, 14 + 4b, y + 36, b))
    end
    push!(p, "</svg>"); join(p)
end

function band_svg(title, series, band; W=760, H=330)
    L, R, T, B = 64, 20, 30, 62
    sx(v) = L + v / XMAX * (W - L - R); sy(v) = H - B - clamp(v, 0, YMAX) / YMAX * (H - T - B)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(esc(title))'><rect width='100%' height='100%' fill='white'/>",
         "<text x='$L' y='18' font-size='14' font-weight='bold'>$(esc(title))</text>"]
    for v in range(0, YMAX; length=5)
        push!(p, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#e5e9ec'/><text x='%d' y='%.1f' font-size='10' text-anchor='end'>%.0f</text>", L, W - R, sy(v), sy(v), L - 6, sy(v) + 3, v))
    end
    for v in 0:2:floor(Int, XMAX)
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%d</text>", sx(v), H - B + 14, v))
    end
    push!(p, "<line x1='$L' y1='$(H-B)' x2='$(W-R)' y2='$(H-B)' stroke='#555'/><line x1='$L' y1='$T' x2='$L' y2='$(H-B)' stroke='#555'/>")
    push!(p, "<text x='$((L+W-R)÷2)' y='$(H-B+30)' font-size='11' text-anchor='middle'>Day</text><text transform='translate(14 $((T+H-B)÷2)) rotate(-90)' font-size='11' text-anchor='middle'>Mean Cells</text>")
    for s in series
        up = [@sprintf("%.1f,%.1f", sx(t), sy(band.hi[i, s.state])) for (i, t) in enumerate(band.ts) if !isnan(band.hi[i, s.state])]
        dn = [@sprintf("%.1f,%.1f", sx(t), sy(band.lo[i, s.state])) for (i, t) in reverse(collect(enumerate(band.ts))) if !isnan(band.lo[i, s.state])]
        isempty(up) || push!(p, "<polygon points='$(join(vcat(up, dn), " "))' fill='$(s.color)' opacity='0.18'/>")
        pts = join([@sprintf("%.1f,%.1f", sx(t), sy(band.mid[i, s.state])) for (i, t) in enumerate(band.ts)], " ")
        push!(p, "<polyline points='$pts' fill='none' stroke='$(s.color)' stroke-width='2'$(s.dashed ? " stroke-dasharray='7,4'" : "")/>")
        fill = s.filled ? s.color : "white"
        push!(p, join([@sprintf("<circle cx='%.1f' cy='%.1f' r='3.2' fill='%s' stroke='%s' stroke-width='1.4'/>", sx(x), sy(y), fill, s.color) for (x, y) in zip(s.x, s.y)]))
    end
    for (i, s) in enumerate(series)
        lx = L + (i - 1) * 170; fill = s.filled ? s.color : "white"
        push!(p, "<circle cx='$(lx+4)' cy='$(H-16)' r='4' fill='$fill' stroke='$(s.color)' stroke-width='1.4'/><text x='$(lx+12)' y='$(H-12)' font-size='11' fill='$(s.color)'>$(esc(s.name))</text>")
    end
    push!(p, "</svg>"); join(p)
end

function group_plots(g, band)
    if g.kind in (:s1, :s1b)
        d = g.datasets[1]
        return [band_svg("$(g.title) — bootstrap-optimal fit with 95% band", [(name="observed", x=d.x, y=d.y, state=1, color=BLUE, filled=true, dashed=false)], band)]
    elseif g.kind == :s2
        cols = [PURPLE, GREEN, ORANGE]
        return [band_svg("$(g.title) — bootstrap-optimal fit with 95% band",
            [(name="$(dose) µM", x=d.x, y=d.y, state=i, color=cols[i], filled=true, dashed=false) for (i, (dose, d)) in enumerate(zip(stage2[g.line].doses, g.datasets))], band)]
    elseif g.kind == :s3
        out = String[]
        for ratio in ["25S75R", "50S50R", "75S25R"]
            i0 = findfirst(e -> e.ratio == ratio && e.cond == "Ce0", coculture_layout); i1 = findfirst(e -> e.ratio == ratio && e.cond == "Ce1", coculture_layout)
            d0, d1 = co_data[i0], co_data[i1]
            push!(out, band_svg("$ratio — Ce0 solid / filled, Ce1 dashed / open", [
                (name="Ce0 sensitive", x=d0.S[1], y=d0.S[2], state=2i0-1, color=BLUE, filled=true, dashed=false),
                (name="Ce0 resistant", x=d0.R[1], y=d0.R[2], state=2i0, color=ORANGE, filled=true, dashed=false),
                (name="Ce1 sensitive", x=d1.S[1], y=d1.S[2], state=2i1-1, color=BLUE, filled=false, dashed=true),
                (name="Ce1 resistant", x=d1.R[1], y=d1.R[2], state=2i1, color=ORANGE, filled=false, dashed=true)], band))
        end
        return out
    else
        return [band_svg("Run $(c.run), ratio $(c.ratio) S:R", [
            (name="sensitive", x=c.xN[.!c.bN], y=c.yN[.!c.bN], state=2i-1, color=BLUE, filled=true, dashed=false),
            (name="resistant", x=c.xC[.!c.bC], y=c.yC[.!c.bC], state=2i, color=ORANGE, filled=true, dashed=false)], band) for (i, c) in enumerate(lr_cc)]
    end
end

inherited_note = Dict(
    :s1 => "No inherited parameters.",
    :s1b => "r and K held at their Stage 1 logistic estimates; u_R = 0.955 (0.5% FBS). Candidates are the 13 low-resource mechanisms from the <a href='../lr_mechanisms/report.html'>LR mechanisms</a> report.",
    :s2 => "r and K held at their Stage 1 logistic estimates; EC50 = 1 µM.",
    :s3 => "r from Stage 1 and λ, H (and β for the inherited-kill models) from Stage 2, held at their point estimates.",
    :s3b => "r from Stage 1; K, α from Stage 3 and d^LR from Stage 1b where a model inherits them; u_R = 0.955.")

summary_rows = String[]; sections = String[]
for g in groups
    r = results[g.id].res; band = results[g.id].band; B = size(r.bic, 1)
    best = r.cands[1]; sel = r.cands[r.selected]
    removed = setdiff(best.pnames, sel.pnames); added = setdiff(sel.pnames, best.pnames)
    push!(summary_rows, "<tr><td>$(esc(g.title))</td><td>$(esc(best.label)) (k=$(length(best.params)))</td><td><b>$(esc(sel.label))</b> (k=$(length(sel.params)))</td>" *
        @sprintf("<td class='num'>%.0f%%</td><td class='num'>%+.1f</td>", pct(r.pick, r.selected), sel.bic - best.bic) *
        "<td>$(r.selected == 1 ? "same as BIC-optimal" : isempty(added) ? "removes " * esc(join(removed, ", ")) : "replaces " * esc(join(removed, ", ")) * " with " * esc(join(added, ", ")))</td></tr>")
    dbic = r.bic .- minimum(r.bic; dims=2)
    trs = String[]
    for (j, c) in enumerate(r.cands)
        l, m, h = qs(dbic[:, j])
        push!(trs, "<tr$(j == r.selected ? " class='win'" : "")><td>$(esc(c.label))</td><td class='num'>$(r.ks[j])</td>" *
            @sprintf("<td class='num'>%.2f</td><td class='num'>%.2f</td><td class='num'>%.1f (%.1f to %.1f)</td><td class='num'>%.0f%%</td><td class='num'>%.0f%%</td></tr>",
                     c.bic, c.bic - best.bic, m, l, h, pct(r.winner, j), pct(r.pick, j)))
    end
    ptrs = String[]
    for (i, n) in enumerate(sel.pnames)
        l, m, h = qs(r.pars[r.selected][:, i])
        rm_ = removable(l, h, sel.bounds[i])
        push!(ptrs, "<tr><td>$(esc(n))</td><td class='num'>$(shown(sel.params[i]))</td><td class='num'>$(fmt(m))</td><td class='num'>$(fmt(l)) to $(fmt(h))</td><td>$(rm_ ? "<b>interval reaches 0 — candidate for removal</b>" : "")</td></tr>")
    end
    ptable = isempty(sel.pnames) ? "<p class='muted'>No free parameters: every value is inherited.</p>" :
        "<div class='tscroll'><table class='params'><tr><th>Parameter</th><th>Estimate</th><th>Bootstrap median</th><th>95% bootstrap interval</th><th></th></tr>$(join(ptrs))</table></div>"
    push!(sections, "<section class='panel' id='$(g.id)'><h3>$(esc(g.title))</h3>" *
        "<p class='muted'>Replicates generated from: $(esc(r.ref.label)) (the most complex candidate). $(inherited_note[g.kind])</p>" *
        freq_svg(r) *
        "<div class='tscroll'><table><tr><th>Candidate</th><th>k</th><th>BIC</th><th>ΔBIC</th><th>Bootstrap ΔBIC, median (95%)</th><th>Lowest BIC</th><th>Parsimony pick</th></tr>$(join(trs))</table></div>" *
        "<h4>Bootstrap-optimal model: $(esc(sel.label))</h4>" * join(["<div class='eq'>\\[$t\\]</div>" for t in bs_tex(g, sel)]) * ptable *
        join(group_plots(g, band)) * "</section>")
end

bnav = "<nav aria-label='New data report navigation' class='topnav'><a href='../../../../index.html' class='nav nav-home' aria-label='Reports home' title='Reports home'>$home_svg</a>" *
    nav_link("Summary", "../../../../outputs/results/summary.html", false) *
    nav_link("New data", "../../../../outputs/new_data_report/report.html", false) *
    nav_link("Candidate models", "../reports/a2780_low_resource_modeling_plan.html", false) *
    nav_link("New data fitting", "../new_data_fits/report.html", false) *
    nav_link("Bootstrapping", "report.html", true) *
    nav_link("LR mechanisms", "../lr_mechanisms/report.html", false) *
    nav_link("Low-resource fitting", "../low_resource_staged_fits/report.html", false) * "</nav>"

bhtml = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Bootstrapping</title>
<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}};</script><script defer src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>
<style>body{font:14px/1.5 system-ui,sans-serif;margin:0;background:#f4f7f8;color:#17212b}main{max-width:1250px;margin:auto;padding:24px}h1{color:#124b6e;margin:6px 0}h2{color:#124b6e;margin-top:30px}
.topnav{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:0 0 16px}.nav{display:inline-flex;align-items:center;min-height:34px;padding:0 11px;border-radius:8px;color:#176b87;background:#eef5f6;text-decoration:none;font:650 14px system-ui,sans-serif;white-space:nowrap}.nav-home{padding:6px 9px}.nav-active{color:#fff;background:#176b87}
.intro,.panel{background:#fff;border:1px solid #d9e2e6;border-radius:12px;padding:16px 18px;margin:14px 0;min-width:0;overflow-wrap:anywhere}.intro{border-left:6px solid #176b87}.note{background:#fff6d9;border-radius:8px;padding:10px 14px}
.eq{overflow-x:auto;overflow-y:hidden;background:#fbfdfd;border:1px solid #e4ecef;border-radius:8px;padding:4px 12px;margin:8px 0}svg{display:block;width:100%;height:auto;margin:10px 0}svg.bic{max-width:900px}.tscroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}td,th{border:1px solid #d3dde2;padding:6px 8px;text-align:left;vertical-align:top}th{background:#e9f3f8}td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}tr.win td{background:#eef8f1;font-weight:600}.muted{color:#5b6b73}h3{margin:0 0 6px}h4{margin:16px 0 4px}
code{overflow-wrap:anywhere}mjx-container:not([display]){display:inline-block;max-width:100%;overflow-x:auto;vertical-align:middle}</style></head><body><main>
$bnav
<h1>Bootstrapping — model selection that favours fewer terms</h1>
<div class='intro'><p>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")). Every candidate from the <a href='../new_data_fits/report.html'>New data fitting</a> report is refitted to <b>$BOOT_B bootstrap replicates</b> of its data, plus pruned candidates that remove terms (Stage 2: no delay; Stage 3b: no low-resource loss).</p>
<ul><li><b>Replicates:</b> the most complex candidate's fit plus its residuals, resampled in circular blocks of $BLOCK days within each trajectory so day-to-day correlation is kept; day 0 stays as observed.</li>
<li><b>Parsimony rule (per replicate):</b> take the lowest-BIC candidate, then prefer any simpler candidate whose BIC is no more than <b>$(Int(PER_PARAM_BIC)) per removed parameter</b> higher (removing two terms may cost up to $(Int(2PER_PARAM_BIC)) BIC). Among eligible models the one with the fewest parameters wins.</li>
<li><b>Bootstrap-optimal model:</b> the candidate the rule picks most often. Its parameters get 95% bootstrap percentile intervals; a parameter whose interval reaches zero is flagged as removable.</li></ul>
<p class='note'>Inherited parameters are held at their point estimates, so upstream uncertainty is not carried forward; the bootstrap intervals describe each stage on its own. Low-resource fits use the same excluded points as the fitting report.</p></div>
<h2>Summary</h2>
<div class='tscroll'><table><tr><th>Fit</th><th>BIC-optimal (original data)</th><th>Bootstrap-optimal (parsimony)</th><th>Picked in</th><th>BIC cost</th><th>Change</th></tr>$(join(summary_rows))</table></div>
<h2>Stages</h2>
$(join(sections, "\n"))
<p class='muted'>Outputs: <code>bootstrap_model_selection.csv</code>, <code>bootstrap_parameters.csv</code> in <code>Modeling_Approaches/02_mechanical_automatic_package/outputs/bootstrap_fits/</code>. Generated by <code>run_bootstrap_fits.jl</code>.</p>
</main></body></html>"""

write(joinpath(BOOT_DIR, "report.html"), bhtml)
mkpath(BOOT_DOCS)
cp(joinpath(BOOT_DIR, "report.html"), joinpath(BOOT_DOCS, "report.html"); force=true)
println("\nWrote $(joinpath(BOOT_DIR, "report.html"))\n  and $(joinpath(BOOT_DOCS, "report.html"))")

# refresh the Summary / Presentation / Explorer data so the site shows this run
include(joinpath(@__DIR__, "export_results_bundle.jl"))
export_results_bundle(groups, results)
# the Final Report embeds the chapter pages at build time, so rebuild it too
try
    run(`$(Base.julia_cmd()) --startup-file=no $(joinpath(ROOT, "build_final_report.jl"))`)
catch err
    @warn "Final report rebuild failed" err
end
