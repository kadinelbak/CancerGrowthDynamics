# HTML report for the optimal-control experiments (run_control.jl checkpoint).
#   julia --project=. control_report.jl
# Writes outputs/control/report.html and docs/Modeling_Approaches/02_mechanical_automatic_package/outputs/control/report.html

include(joinpath(@__DIR__, "control_fasting.jl"))     # includes control_model.jl and control_mpc.jl; lr_plant lives there
const ZOO_DIR = joinpath(@__DIR__, "outputs", "zoo_fits")
const ZOO_DOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "zoo_fits")
const BOOT_B = ZOO_STATE.B
include(joinpath(@__DIR__, "zoo_report.jl"))           # hesc, tip(), INFO_TIPS_HEAD

# Run length: CONTROL_DAYS (default closed_loop_days = 14). CONTROL_DAYS=30 renders the long report from
# outputs/control_30d/ (run_control.jl with the same setting) into docs/.../outputs/control_30d/.
const DAYS = parse(Float64, get(ENV, "CONTROL_DAYS", string(CFG["schedule"]["closed_loop_days"])))
const LONG = DAYS != Float64(CFG["schedule"]["closed_loop_days"])
const DATA_DAYS = Float64(CFG["schedule"]["data_days"])
const SUB = LONG ? "control_$(round(Int, DAYS))d" : "control"
const DL = "day $(round(Int, DAYS))"                          # table headers and tips: "N day 14" / "N day 30"
const CTRL14 = joinpath(@__DIR__, "outputs", "control")       # tests and assumption checks live with the 14-day run
const COUT = joinpath(@__DIR__, "outputs", SUB)
const CDOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", SUB)
C = deserialize(joinpath(COUT, "control_state.jls"))
AN = AN_STATE
const Δ0 = Float64(CFG["schedule"]["decision_interval_days"])
const NRUN = round(Int, DAYS / Δ0)
const NWORD = get(Dict(7 => "seven"), NRUN, string(NRUN))
# long runs: a dashed line where the data end on every time plot
data_end(sx, top, bottom) = LONG ? @sprintf("<line x1='%.1f' x2='%.1f' y1='%d' y2='%d' stroke='#7c3aed' stroke-dasharray='4,4'/><text x='%.1f' y='%d' font-size='10' fill='#7c3aed'>end of data (day %d)</text>",
                                            sx(DATA_DAYS), sx(DATA_DAYS), top, bottom, sx(DATA_DAYS) + 4, top + 10, round(Int, DATA_DAYS)) : ""
const POL_COL = Dict("No treatment" => "#6b7280", "Constant 1 µM" => "#d95f02", "Maximum measured dose (1.47 µM)" => "#b91c1c",
                     "Periodic 1 µM, 2 days on / 2 off" => "#8b5cf6", "Threshold adaptive (Zhang-style)" => "#3b8b4f",
                     "MPC, measured doses only" => "#0b6e99")
polcol(n) = get(POL_COL, n, "#0f766e")
poldash(n) = startswith(n, "MPC demonstration") ? "6,4" : ""
f2(x) = isfinite(x) ? fmt(x) : "—"
cells(x) = isfinite(x) ? @sprintf("%.0f", x) : "—"
yn(b) = b ? "<span class='flag'>yes</span>" : "no"

# ── plotting ─────────────────────────────────────────────────────────────────
function ts_svg(title, series; xmax=DAYS, ymax=nothing, hline=nothing, W=760, H=320, ylabel="Mean Cells")
    ym = ymax === nothing ? ceil(maximum(maximum(s.y) for s in series) * 1.12; sigdigits=2) : ymax
    ym <= 0 && (ym = 1.0)
    L, R, T, B = 64, 16, 28, 72
    sx(v) = L + v / xmax * (W - L - R); sy(v) = H - B - clamp(v, 0, ym) / ym * (H - T - B)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(hesc(title))'><rect width='100%' height='100%' fill='white'/><text x='$L' y='17' font-size='13' font-weight='bold'>$(hesc(title))</text>"]
    for v in range(0, ym; length=5)
        push!(p, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#e5e9ec'/><text x='%d' y='%.1f' font-size='10' text-anchor='end'>%s</text>", L, W - R, sy(v), sy(v), L - 6, sy(v) + 3, ym < 5 ? @sprintf("%.2f", v) : @sprintf("%.0f", v)))
    end
    for v in 0:(xmax > 16 ? 5 : 2):floor(Int, xmax)
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%d</text>", sx(v), H - B + 14, v))
    end
    push!(p, data_end(sx, T, H - B))
    push!(p, "<line x1='$L' y1='$(H-B)' x2='$(W-R)' y2='$(H-B)' stroke='#555'/><text x='$((L+W-R)÷2)' y='$(H-B+28)' font-size='11' text-anchor='middle'>Day</text><text transform='translate(13 $((T+H-B)÷2)) rotate(-90)' font-size='11' text-anchor='middle'>$(hesc(ylabel))</text>")
    if hline !== nothing
        push!(p, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#b45309' stroke-dasharray='5,4'/><text x='%d' y='%.1f' font-size='10' fill='#b45309' text-anchor='end'>%s</text>", L, W - R, sy(hline[1]), sy(hline[1]), W - R - 4, sy(hline[1]) - 4, hesc(hline[2])))
    end
    for s in series
        hasproperty(s, :lo) && push!(p, "<polygon points='" * join(vcat([@sprintf("%.1f,%.1f", sx(t), sy(v)) for (t, v) in zip(s.t, s.hi)],
                                                                            [@sprintf("%.1f,%.1f", sx(t), sy(v)) for (t, v) in zip(reverse(s.t), reverse(s.lo))]), " ") * "' fill='$(s.color)' opacity='0.16'/>")
        push!(p, "<polyline points='" * join([@sprintf("%.1f,%.1f", sx(t), sy(v)) for (t, v) in zip(s.t, s.y)], " ") * "' fill='none' stroke='$(s.color)' stroke-width='2'$(isempty(s.dash) ? "" : " stroke-dasharray='$(s.dash)'")/>")
    end
    for (i, s) in enumerate(unique(s -> s.label, series))
        x = L + ((i - 1) % 3) * 230; y = H - 24 + ((i - 1) ÷ 3) * 13
        push!(p, "<line x1='$x' x2='$(x+18)' y1='$(y-4)' y2='$(y-4)' stroke='$(s.color)' stroke-width='2'$(isempty(s.dash) ? "" : " stroke-dasharray='$(s.dash)'")/><text x='$(x+22)' y='$y' font-size='10'>$(hesc(first(s.label, 38)))</text>")
    end
    push!(p, "</svg>"); join(p)
end
# dose schedules as step plots, one row per policy
function dose_svg(title, pols; W=760, umax=3.0)
    rowH = 26; H = 30 + rowH * length(pols) + 20; L = 210; R = 16
    sx(t) = L + t / DAYS * (W - L - R)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(hesc(title))'><rect width='100%' height='100%' fill='white'/><text x='8' y='17' font-size='13' font-weight='bold'>$(hesc(title))</text>"]
    for (i, (name, U)) in enumerate(pols)
        y0 = 30 + rowH * (i - 1)
        push!(p, "<text x='8' y='$(y0+15)' font-size='10'>$(hesc(first(name, 36)))</text>")
        for (k, u) in enumerate(U)
            u > 0 || continue
            h = (rowH - 6) * min(u / umax, 1.0)
            push!(p, @sprintf("<rect x='%.1f' y='%.1f' width='%.1f' height='%.1f' fill='%s'><title>days %g–%g: %.2f µM</title></rect>", sx((k - 1) * Δ0), y0 + rowH - 4 - h, sx(k * Δ0) - sx((k - 1) * Δ0) - 1, h, polcol(name), (k - 1) * Δ0, k * Δ0, u))
        end
        push!(p, "<line x1='$L' x2='$(W-R)' y1='$(y0+rowH-4)' y2='$(y0+rowH-4)' stroke='#d3dde2'/>")
    end
    for v in 0:(DAYS > 16 ? 4 : 2):round(Int, DAYS)
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%d</text>", sx(v), H - 4, v))
    end
    push!(p, data_end(sx, 24, H - 16))
    push!(p, "</svg>"); join(p)
end
function flow_svg()
    boxes = [("Stage 1", "untreated monocultures", "r, K, θ"), ("Stage 1b", "low-resource monocultures", "low-resource mechanism per line"),
             ("Stage 2", "treated monocultures", DRUG_HILL ? "drug chain k, Hill h, kill β" : "drug chain k, linear kill β" * (DRUG_NS ? " (Norton–Simon)" : "")), ("Stage 3", "Ce0 + Ce1 cocultures", "shared K, " * join(filter(n -> n in zgroup("s3").cands[1].pnames, ["α_NC", "α_CN", "d_N", "β_N", "β_C"]), ", ")),
             ("Stage 3b", "low-resource cocultures", "shared K (low resource)"), ("Stage 4", "low resource + cisplatin", "ε_N, ε_C (low-resource strength)"),
             ("Control", "MPC plant", "+ dose u, toxicity B, T")]
    pos = [(20, 40), (20, 150), (230, 40), (440, 40), (440, 150), (650, 150), (650, 40)]
    arrows = [(1, 2), (1, 3), (3, 4), (1, 4), (2, 5), (4, 5), (4, 6), (5, 6), (4, 7), (6, 7)]
    W, H, bw, bh = 860, 240, 190, 76
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='Parameter flow between stages'><rect width='100%' height='100%' fill='white'/><defs><marker id='ah' viewBox='0 0 10 10' refX='9' refY='5' markerWidth='7' markerHeight='7' orient='auto'><path d='M0,0L10,5L0,10z' fill='#52636b'/></marker></defs>"]
    cx(i) = pos[i][1] + bw / 2; cy(i) = pos[i][2] + bh / 2
    for (a, b) in arrows
        x1, y1, x2, y2 = cx(a), cy(a), cx(b), cy(b)
        dx, dy = x2 - x1, y2 - y1; t1 = min(abs(bw / 2 / (dx == 0 ? 1e-9 : dx)), abs(bh / 2 / (dy == 0 ? 1e-9 : dy)))
        push!(p, @sprintf("<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' stroke='#52636b' stroke-width='1.4' marker-end='url(#ah)'/>", x1 + dx * t1, y1 + dy * t1, x2 - dx * t1, y2 - dy * t1))
    end
    for (i, (h, s, q)) in enumerate(boxes)
        x, y = pos[i]
        push!(p, "<rect x='$x' y='$y' width='$bw' height='$bh' rx='9' fill='$(i == 7 ? "#e6f1f5" : "#f7fafb")' stroke='#176b87'/><text x='$(x+10)' y='$(y+20)' font-size='13' font-weight='bold' fill='#124b6e'>$h</text><text x='$(x+10)' y='$(y+37)' font-size='11'>$(hesc(s))</text><text x='$(x+10)' y='$(y+56)' font-size='10.5' fill='#52636b'>→ $(hesc(q))</text>")
    end
    push!(p, "</svg>"); join(p)
end

# ── pieces ──────────────────────────────────────────────────────────────────
main_of(env, scen, ratio) = only(filter(x -> x.env == env && x.scen == scen && x.ratio == ratio, C[:main]))
const POL_ORDER = ["No treatment", "Constant 1 µM", "Maximum measured dose (1.47 µM)", "Periodic 1 µM, 2 days on / 2 off", "Threshold adaptive (Zhang-style)", "MPC, measured doses only"]
pol_order(m) = vcat([p for p in POL_ORDER if haskey(m.policies, p)], sort([p for p in keys(m.policies) if !(p in POL_ORDER)]))

const T_TOTAL = "Total cells N + C under each dosing policy (model simulation; no data were collected for these schedules). The dashed orange line is the cancer-burden bound B_cancer (an assumption: 80% of the fitted capacity K). Lower is less cancer burden."
const T_NC = "Sensitive (A2780, solid) and resistant (A2780cis, dashed) cells under the data-supported MPC schedule and under constant 1 µM. Whether MPC keeps sensitive cells to hold resistant cells back shows here: the resistant curve under MPC vs constant dosing."
const T_TOX = "Toxicity burden T = 1 − H (fraction of healthy host tissue lost; literature model of de Pillis & Radunskaya 2001). The dashed line is the bound B_tox (0.25: healthy tissue kept above 75%). Illustrative: no toxicity was measured in these experiments."
const T_DOSE = "The cisplatin concentration (µM) each policy puts in the medium for every 2-day interval between media changes. Bar height is the dose (hover a bar for its value)."
const T_OUT = "Outcomes on $DL of the model simulation. Resistant fraction = C/(N+C). Max toxicity is the highest T reached; cumulative dose is ∫u dt in µM·days. A violation means the trajectory touched that bound at any time."

function outcome_table(m)
    rows = String[]
    for pn in pol_order(m)
        v = m.policies[pn]; o = v.o
        push!(rows, trow([td("<span style='color:$(polcol(pn))'>■</span> " * hesc(pn)), tn(cells(o.N_final)), tn(cells(o.C_final)), tn(cells(o.total_final)),
                          tn(@sprintf("%.2f", o.resistant_fraction)), tn(@sprintf("%.3f", o.max_tox)), tn(@sprintf("%.1f", o.cum_dose)), td(yn(o.cancer_violation)), td(yn(o.tox_violation))];
                     cls=startswith(pn, "MPC, measured") ? "win" : ""))
    end
    table([("Policy", "The dosing rule. MPC = model-predictive control; 'measured doses only' restricts it to 0, 0.62, 1 and 1.47 µM, the doses that were measured."),
           ("N $DL", "Sensitive A2780 cells at $DL."), ("C $DL", "Resistant A2780cis cells at $DL."), ("N + C", "Total cells at $DL."),
           ("Resistant fraction", "C / (N + C) at $DL. Lower keeps the population drug-sensitive."), ("Max toxicity", T_TOX),
           ("Cumulative dose (µM·d)", "∫u dt over $(round(Int, DAYS)) days."), ("Cancer bound hit", "Whether N + C reached B_cancer at any time."), ("Toxicity bound hit", "Whether T reached B_tox at any time.")], rows) |> x -> withtip(x, T_OUT, "Lower is better for cell counts, resistant fraction, toxicity and dose; No is better for the bound columns.")
end

function env_section(env)
    envname = env == :normal ? "normal medium (10% FBS), Stage 3 model" : "low resource (0.5% FBS), Stage 3 model + low-resource modules, ε(u)"
    secs = String[]
    for scen in (:A, :B)
        m = main_of(env, scen, "50:50")
        ts = m.policies["No treatment"].o.ts
        tot = [(label=pn, color=polcol(pn), dash=poldash(pn), t=ts, y=m.policies[pn].o.N .+ m.policies[pn].o.C) for pn in pol_order(m)]
        mp = m.policies["MPC, measured doses only"].o; cst = m.policies["Constant 1 µM"].o
        nc = [(label="MPC: sensitive N", color="#0b6e99", dash="", t=ts, y=mp.N), (label="MPC: resistant C", color="#0b6e99", dash="6,4", t=ts, y=mp.C),
              (label="Constant 1 µM: sensitive N", color="#d95f02", dash="", t=ts, y=cst.N), (label="Constant 1 µM: resistant C", color="#d95f02", dash="6,4", t=ts, y=cst.C)]
        tox = [(label=pn, color=polcol(pn), dash=poldash(pn), t=ts, y=m.policies[pn].o.T) for pn in pol_order(m)]
        W = default_weights(plant_for_report(env, scen))
        ratios = join(["<h4>Seeding ratio $(r)</h4>" * outcome_table(main_of(env, scen, r)) for r in ("25:75", "50:50", "75:25")])
        push!(secs, """<section class='panel' id='$(env)-$(scen)'><h3>Scenario $(hesc(m.name))</h3>
<div class='grid'><div>$(withtip(ts_svg("Total N + C, 50:50 start", tot; hline=(m.Bc, "B_cancer"), ymax=ceil(max(m.Bc * 1.15, maximum(maximum(s.y) for s in tot) * 1.05); sigdigits=2)), T_TOTAL, "Lower is better, and the curve must stay under the dashed line."))</div>
<div>$(withtip(ts_svg("Sensitive and resistant cells", nc), T_NC, "Lower is better, above all for the resistant (dashed) curve."))</div>
<div>$(withtip(dose_svg("Dose schedules (µM per 2-day interval)", [(pn, m.policies[pn].U) for pn in pol_order(m)]), T_DOSE, "Neither on its own: less drug is better when the cell counts come out similar."))</div>
<div>$(withtip(ts_svg("Toxicity burden T", tox; hline=(W.Btox, "B_tox"), ymax=max(W.Btox * 1.3, maximum(maximum(s.y) for s in tox) * 1.1), ylabel="T (fraction)"), T_TOX, "Lower is better, and it must stay under the dashed line."))</div></div>
$ratios</section>""")
    end
    "<h2>Results: $(envname)</h2>" * join(secs)
end
plant_for_report(env, scen) = env == :normal ? scenario_plant(scen) : lr_plant(scen)

function robust_section()
    secs = String[]; rrows = String[]
    for env in (:normal, :lowres), scen in (:A, :B)
        R = C[:robust][(env, scen)]; ts = R.ts
        m = main_of(env, scen, "50:50")
        q(v, pr) = isempty(v) ? NaN : quantile(v, pr)
        band(pn, f) = begin
            M = hcat([f(d[pn]) for d in R.draws if d[pn] !== nothing]...)
            (t=ts, y=[median(M[i, :]) for i in axes(M, 1)], lo=[q(M[i, :], 0.025) for i in axes(M, 1)], hi=[q(M[i, :], 0.975) for i in axes(M, 1)])
        end
        ser = [merge(band(pn, o -> o.N .+ o.C), (label=pn * " (median, 95%)", color=polcol(pn), dash="")) for pn in ("Constant 1 µM", "MPC, measured doses only")]
        serC = [merge(band(pn, o -> o.C), (label=pn * ": resistant C", color=polcol(pn), dash="6,4")) for pn in ("Constant 1 µM", "MPC, measured doses only")]
        push!(secs, "<div>" * withtip(ts_svg("$(env == :normal ? "Normal medium" : "Low resource") · scenario $scen: N + C over bootstrap parameter sets", vcat(ser, serC); hline=(m.Bc, "B_cancer")),
            "Each schedule (fixed, planned with the best-fit model) re-simulated on $(length(R.draws)) bootstrap parameter sets. Line: median; band: 2.5–97.5 percentiles. A wide band means the fitted parameters leave the outcome uncertain; a band crossing the dashed line means the cancer bound may be violated under plausible parameters.", "Lower and narrower is better; the band should stay under the dashed line.") * "</div>")
        for pn in ("No treatment", "Constant 1 µM", "Periodic 1 µM, 2 days on / 2 off", "Threshold adaptive (Zhang-style)", "MPC, measured doses only")
            ds = [d[pn] for d in R.draws if d[pn] !== nothing]; isempty(ds) && continue
            tot = [d.tot for d in ds]; cf = [d.Cf for d in ds]; nf = [d.Nf for d in ds]
            push!(rrows, trow([td("$(env == :normal ? "normal" : "low resource")"), td(string(scen)), td(hesc(pn)), tn(length(ds)),
                               tn(@sprintf("%.0f%%", 100 * mean(d.cv for d in ds))), tn(@sprintf("%.0f%%", 100 * mean(d.tv for d in ds))),
                               tn("$(cells(median(nf))) ($(cells(q(nf, .025)))–$(cells(q(nf, .975))))"), tn("$(cells(median(cf))) ($(cells(q(cf, .025)))–$(cells(q(cf, .975))))"),
                               tn("$(cells(median(tot))) ($(cells(q(tot, .025)))–$(cells(q(tot, .975))))"), tn(@sprintf("%.1f", median(d.cum for d in ds)))]))
        end
        if !isempty(R.closed)
            cl = R.closed; tot = [c.tot for c in cl]
            push!(rrows, trow([td("$(env == :normal ? "normal" : "low resource")"), td(string(scen)), td("<b>Closed-loop MPC</b> (re-plans every 2 days on the true plant)"), tn(length(cl)),
                               tn(@sprintf("%.0f%%", 100 * mean(c.cv for c in cl))), tn(@sprintf("%.0f%%", 100 * mean(c.tv for c in cl))), tn("—"), tn("—"),
                               tn("$(cells(median(tot))) ($(cells(q(tot, .025)))–$(cells(q(tot, .975))))"), tn(@sprintf("%.1f", median(c.cum for c in cl)))]; cls="win"))
        end
    end
    """<h2>Robustness over bootstrap parameter sets</h2>
<div class='intro'><p>The schedules above are planned with the best-fit parameters. Here each one is replayed on bootstrap parameter sets: the growth, drug and coculture parameters of the same bootstrap replicate of their own stages (an approximation: each stage was bootstrapped with earlier-stage values held at their estimates). Threshold adaptive therapy re-decides on each plant. The closed-loop row lets MPC plan with the best-fit model but observe the bootstrap plant every 2 days, as it would in an experiment.</p></div>
<div class='grid'>$(join(secs))</div>
$(table([("Environment", "Normal medium (Stage 3 model) or low resource (Stage 3 model with the fitted low-resource modules; their strength ε(u) depends on the dose and equals the Stage 4 fit at 1 µM)."), ("Scenario", "A: the coculture fit as selected (Stage 3); B: monoculture drug response for both lines with crowding from the untreated cocultures."),
         ("Policy", "The dosing rule; its schedule is planned with the best-fit model."), ("Draws", "Bootstrap parameter sets simulated."),
         ("Cancer bound hit", "Share of parameter sets in which N + C reaches B_cancer at any time."), ("Toxicity bound hit", "Share of parameter sets in which T reaches B_tox."),
         ("N $DL", "Median sensitive cells (95% range) across parameter sets."), ("C $DL", "Median resistant cells (95% range)."),
         ("N + C $DL", "Median total (95% range)."), ("Dose (µM·d)", "Median cumulative dose.")], rrows))"""
end

function sweep_section()
    rows = String[]
    for s in C[:sweeps]
        o = s.res.o
        push!(rows, trow([td(hesc(s.kind)), td(string(s.scen)), td(hesc(s.label)), td(join([@sprintf("%.3g", u) for u in s.res.U], " ")),
                          tn(cells(o.N_final)), tn(cells(o.C_final)), tn(cells(o.total_final)), tn(@sprintf("%.2f", o.resistant_fraction)), tn(@sprintf("%.3f", o.max_tox)), tn(@sprintf("%.1f", o.cum_dose))]))
    end
    srows = [trow([td(string(s.scen)), tn(string(s.h)), tn(cells(s.N)), tn(cells(s.C)), tn(@sprintf("%.5f", s.T))]) for s in C[:steps]]
    """<h2>Sensitivity of the controller</h2>
<div class='intro'><p>Data-supported MPC in normal medium from the 50:50 start, varying one setting at a time: the cost weights, the toxicity model (host recovery rate r_H, dose mapping s_B, bound B_tox), the planning horizon, the decision interval and the simulation length. If the schedule barely changes across plausible weights, the recommendation does not hinge on an arbitrary weight choice.</p></div>
$(table([("Setting", "Which assumption is varied."), ("Scen.", "Drug-response scenario (A or B)."), ("Value", "The value used; everything else is at its default."),
         ("Schedule (µM)", "The dose applied in each 2-day interval (1-day intervals for the 'every 1 day' row" * (LONG ? ")." : "; 14 intervals for the 28-day row).")),
         ("N $DL", "Sensitive cells at the end."), ("C $DL", "Resistant cells at the end."), ("N + C", "Total at the end."),
         ("Resistant fraction", "C/(N + C) at the end."), ("Max toxicity", "Highest T."), ("Dose (µM·d)", "Cumulative dose.")], rows))
<h3>Integration step $(tip("The data-supported MPC schedule re-simulated with smaller fixed RK4 steps. Matching digits mean the step used by the controller (0.05 day) is accurate enough; the step is also checked against an adaptive solver in the tests.", "Better: the rows agree to several digits."))</h3>
$(table([("Scen.", "Scenario."), ("Step (days)", "Fixed RK4 step."), ("N $DL", "Sensitive cells at the end.", AGREE), ("C $DL", "Resistant cells at the end.", AGREE), ("T $DL", "Toxicity at the end.", AGREE)], srows))"""
end

# the "i" on the toxicity equations: what they are, the source, and why the proxy is appropriate
# highest toxicity over every main run that stays within the measured doses (the demonstration MPC is excluded)
const TMAX_DATA = maximum(v.o.max_tox for m in C[:main] for (pn, v) in m.policies if !startswith(pn, "MPC demonstration"))
# the carried drug model (Stage 2b), for the text below: Hill or linear kill, Norton–Simon or not
const DRUG_HILL = zgroup("s2-a2780").cands[1].sys.lines[1].action == :hill
const DRUG_NS = ns_on(zgroup("s2-a2780").cands[1].sys.lines[1])
const DRUG_CURVE = DRUG_HILL ? "the Hill curve" : "the linear kill"
const TOX_TIP = "Toxicity proxy: the damage the same drug exposure would do to healthy host tissue. No toxicity was measured in these in-vitro experiments, so it comes from the literature. " *
    "B = systemic drug burden (built up by the dose, cleared at d_B); T = 1 − H = fraction of healthy tissue lost (0 = none). " *
    "SOURCE: the normal-cell compartment of de Pillis LG & Radunskaya AE (2001), 'A mathematical tumor model with immune resistance and drug therapy: an optimal control approach', " *
    "J Theor Med (now Comput Math Methods Med) 3:79–100. Their drug equation dM/dt = −d₂M + v(t) is dB/dt here, and their normal-cell equation " *
    "(logistic regrowth r₂, fractional kill a₃(1 − e^(−M))) is rewritten for T = 1 − H. Parameter values r_H = 1/day, a_H = 0.1, d_B = 1/day and the ≥ 75% healthy-tissue limit " *
    "are those of de Pillis & Radunskaya as tabulated in DTR-Bench (Luo et al. 2024, arXiv:2405.18610, Table IX). " *
    "WHY APPROPRIATE: it is the standard benchmark for optimal control of chemotherapy. Its saturating kill (1 − e^(−B)) and logistic recovery make toxicity rise with dose " *
    "and recover during drug-free intervals, which is exactly the trade-off a dosing schedule has to balance, and it adds only two states. " *
    "LIMITS: it is a proxy, not a fitted quantity: it constrains the schedule but nothing here measured it. It is a generic host model, not calibrated to cisplatin or to these cells; " *
    "s_B (1 µM in the medium = 1 unit of burden) is an assumption. The report varies s_B, the limit, and r_H down to 0.2/day (about a 5-day recovery, the order of the " *
    "maturation time in Friberg et al. 2002, J Clin Oncol 20:4713) to show how the schedules depend on them. " *
    (TMAX_DATA < CFG["toxicity"]["B_tox"] ? "In these runs T peaks at $(@sprintf("%.3f", TMAX_DATA)) with measured doses, so the $(CFG["toxicity"]["B_tox"]) limit never binds." :
                                            "In these runs T reaches $(@sprintf("%.3f", TMAX_DATA)) with measured doses, so the $(CFG["toxicity"]["B_tox"]) limit does bind.")

function model_section()
    P = scenario_plant(:A); c3 = zgroup("s3").cands[1]
    # normal-medium cell equations generated from the carried Stage 3 model itself (D = u), as on the zoo pages,
    # so they follow whatever drug model was selected (Hill or linear kill, Norton–Simon or not)
    ctex(a, X, Y, nm) = (a.i == 0 && a.v == 0) ? X : (a.i == 0 && a.v == 1) ? "$X+$Y" : "$X+\\alpha_{$nm}$Y"
    S3 = c3.sys
    eqs = vcat(line_tex(S3.lines[1], "N", ctex(S3.alpha[1][2], "N", "C", "NC"), "K"), line_tex(S3.lines[2], "C", ctex(S3.alpha[2][1], "C", "N", "CN"), "K"),
               ["\\frac{dB}{dt}=s_B\\,u-d_B\\,B,\\qquad \\frac{dT}{dt}=a_H\\left(1-e^{-B}\\right)(1-T)-r_H\\,T(1-T)"])
    # low-resource (Stage 4) equations generated from the fitted model itself, as on the zoo pages (D = u)
    c4 = zgroup("s4").cands[1]; S4 = c4.sys
    eqlr = vcat(line_tex(S4.lines[1], "N", ctex(S4.alpha[1][2], "N", "C", "NC"), "K"), line_tex(S4.lines[2], "C", ctex(S4.alpha[2][1], "C", "N", "CN"), "K"))
    v = Dict(zip(c3.pnames, c3.params))
    Ln, Lc = P.sys.lines
    hill = Lc.action == :hill; ecC = ec50_of(Lc.s, P.p, Float64)
    kill1(β) = hill ? β / (ecC^Lc.s.h.v + 1) : β                       # resistant kill at a steady 1 µM effect (before any Norton–Simon factor)
    prow(sym, val, unit, src, note) = trow([td("\\($sym\\)"), tn(val isa Number ? fmt(val) : val), td(unit), td(src), td(hesc(note))])
    t = default_tox()
    prows = [prow("r_N", Ln.s.r.v, "/day", "DATA · Stage 1", "sensitive growth rate"), prow("r_C", Lc.s.r.v, "/day", "DATA · Stage 1", "resistant growth rate"),
             prow("\\theta_C", Lc.s.θ.v, "—", "DATA · Stage 1", "resistant crowding shape"), prow("K", v["K"], "cells", "DATA · Stage 3", "shared carrying capacity"),
             prow("\\alpha_{NC}", v["α_NC"], "—", "DATA · Stage 3", "crowding of N by C"), prow("\\alpha_{CN}", v["α_CN"], "—", "DATA · Stage 3", "crowding of C by N"),
             (haskey(v, "d_N") ? prow("d_N", v["d_N"], "/day", "DATA · Stage 3", "sensitive death in coculture") : ""), prow("\\beta_N", v["β_N"], "/day", "DATA · Stage 3", "max sensitive kill (scenario A; B uses Stage 2)"),
             (hill ? prow("h_N", Ln.s.h.v, "—", "DATA · Stage 2", "Hill slope, sensitive") : ""), prow("k_N", Ln.s.k.v, "/day", "DATA · Stage 2", "drug-effect transit rate, sensitive"),
             prow("\\beta_C", stage3_betaC().value, "/day", "DATA · $(stage3_betaC().source) (A) / Stage 2 (B)", "max resistant kill: $(fmt(stage3_betaC().value)) in A, $(fmt(zgroup("s2-a2780cis").cands[1].params[1])) in B. $(hill ? "Not identifiable on its own: every measured dose is far below EC50,C, so the data fix only the kill at those doses, β_C·E^h_C/(EC50,C^h_C + E^h_C)" : "The kill is linear in the effect, β_C·E (per µM-equivalent per day)") = $(fmt(kill1(stage3_betaC().value))) /day at a steady 1 µM in A and $(fmt(kill1(zgroup("s2-a2780cis").cands[1].params[1]))) /day in B$(ns_on(Lc) ? ", before the Norton–Simon growth factor, which lowers it as the culture fills" : "")"),
             (hill ? prow("h_C", Lc.s.h.v, "—", "DATA · Stage 2", "Hill slope, resistant") : ""), prow("k_C", Lc.s.k.v, "/day", "DATA · Stage 2", "drug-effect transit rate, resistant"),
             (hill ? prow("EC_{50,N}", ec50_of(Ln.s, P.p, Float64), "µM", "DATA (fixed)", "measured IC50 of A2780") : ""),
             (hill ? prow("EC_{50,C}", ec50_of(Lc.s, P.p, Float64), "µM", "DATA (fixed)", "measured IC50 of A2780cis") : ""),
             prow("s_B", t.s_B, "per µM", "ASSUMPTION", "medium µM → systemic burden"), prow("d_B", t.d_B, "/day", "LITERATURE", "drug decay (de Pillis & Radunskaya 2001)"),
             prow("a_H", t.a_H, "—", "LITERATURE", "fractional kill of normal tissue (de Pillis & Radunskaya 2001)"), prow("r_H", t.r_H, "/day", "LITERATURE", "normal-tissue regrowth (de Pillis & Radunskaya 2001)"),
             prow("B_{tox}", CFG["toxicity"]["B_tox"], "—", "LITERATURE", "normal tissue ≥ 75% (de Pillis & Radunskaya 2001)"),
             prow("B_{cancer}", "$(CFG["cancer"]["B_cancer_fraction_of_K"]) K", "cells", "ASSUMPTION", "maximum tolerated N + C"),
             prow("u_{max}", "$(CFG["dose"]["u_max_data"]) (data) / $(CFG["dose"]["u_max_demo"]) (demo)", "µM", "DATA / ASSUMPTION", "highest measured dose / demonstration only")]
    srows = [trow([td("\\($a\\)"), td(b), td(c)]) for (a, b, c) in [("N, C", "cells (Mean Cells units)", "sensitive A2780 and resistant A2780cis counts"),
             ("A_{X,j}", "µM-equivalent", "drug-effect compartments of line X; E_X = A_{X,3} drives the kill"), ("S_{C,j}", "fraction of u_R", "low-resource stress compartments of A2780cis (low resource only)"),
             ("B", "dose units", "systemic drug burden (toxicity side)"), ("T", "fraction 0–1", "toxicity burden = healthy tissue lost"), ("u", "µM", "cisplatin in the medium (control)")]]
    """<h2>The model the controller uses</h2>
<div class='intro'><p>The cell equations are exactly the carried compartment models: the control plant calls the same right-hand side as the fits, with the same parameter values (a test checks that it reproduces the fitted Stage 3 curves at 1 µM to 10⁻⁵). The controller only adds the dose input \\(u\\) and a toxicity proxy. Drug → effect → kill follows \\(u \\to A_{X,1} \\to A_{X,2} \\to A_{X,3}=E_X \\to\\) kill, and separately \\(u \\to B \\to T\\).</p>
$(join(["<div class='eq'>$(i == length(eqs) ? tip(TOX_TIP, "Lower T is better, and it must stay under the 0.25 bound (at least 75% of healthy tissue kept).") : "")\\[$e\\]</div>" for (i, e) in enumerate(eqs)]))
<p><b>Low resource (D = u, u<sub>R</sub> = $(U_R)).</b> $(tip("The low-resource equations as the zoo writes them for the carried Stage 4 model: the cisplatin module and the low-resource module of each line. The drug-effect compartments and the toxicity states B, T are as above.", :none)) The low-resource strengths depend on the dose, \\(\\varepsilon_X(u)=1+(\\varepsilon_X^{S4}-1)\\min(u/1\\,\\mu\\mathrm{M},1)\\): with no drug the model is the low-resource fit alone (\\(\\varepsilon=1\\)), and at 1 µM it is exactly the Stage 4 fit. This is the same plant as the 30-day fasting controller held in low resource throughout, so the two reports agree.</p>
$(join(["<div class='eq'>\\[$e\\]</div>" for e in eqlr]))
<p>Discrete transition over one decision interval (Δ = $(Δ0) days, u constant): \\(x_{k+1}=F(x_k,u_k)\\), fixed-step RK4 with step 0.05 day. Local linearisation \\(A_k=\\partial F/\\partial x\\), \\(B_k=\\partial F/\\partial u\\) by forward-mode automatic differentiation through the RK4 steps (one Jacobian of F over [x; u]: nx + 1 dual directions, about the cost of reverse mode at this size, with no tape).</p></div>
<div class='grid'><div>$(table([("State", "Symbol."), ("Units", "Units."), ("Meaning", "What it represents.")], srows))</div>
<div>$(withtip(flow_svg(), "How parameters flow between stages: each stage fits only its new parameters and inherits the rest. The control plant inherits everything and adds the dose and the toxicity proxy.", :none))</div></div>
$(table([("Parameter", "Symbol in the equations."), ("Value", "Value used by the controller (scenario A, normal medium)."), ("Units", "Units."),
         ("Source", "DATA = fitted to these experiments (stage named); LITERATURE = taken from a cited model; ASSUMPTION = exploratory, swept below, to be replaced by measurements."),
         ("Meaning", "What the parameter does.")], prows))"""
end

function objective_section()
    """<h2>Objective, constraints and barriers</h2>
<div class='intro'>
<div class='eq'>\\[\\min_{u_1,\\dots,u_H}\\ \\sum_{k=1}^{H}\\ell(x_{k+1},u_k)+\\phi(x_{H+1}),\\qquad 0\\le u_k\\le u_{max}\\]</div>
<div class='eq'>\\[\\ell=q_B\\,\\tfrac{N+C}{K}+q_R\\,\\tfrac{C}{K}+q_T\\Big(\\tfrac{T}{B_{tox}}\\Big)^2+r_u\\Big(\\tfrac{u}{u_{max}}\\Big)^2+b_{\\mu_C}\\!\\Big(1-\\tfrac{N+C}{B_{cancer}}\\Big)+b_{\\mu_T}\\!\\Big(1-\\tfrac{T}{B_{tox}}\\Big)\\]</div>
<div class='eq'>\\[b_\\mu(z)=\\begin{cases}-\\mu\\log z, & z\\ge\\delta\\\\ \\mu\\left[\\tfrac12\\left(\\tfrac{z-2\\delta}{\\delta}\\right)^2-\\tfrac12-\\log\\delta\\right], & z<\\delta\\end{cases}\\qquad \\phi=$(CFG["weights"]["terminal_factor"])\\left[q_B\\,\\tfrac{N+C}{K}+q_R\\,\\tfrac{C}{K}+q_T\\Big(\\tfrac{T}{B_{tox}}\\Big)^2\\right]\\]</div>
<p><b>Goal.</b> Keep the total tumour burden low <i>and</i> hold resistant cells back. The burden term \\(q_B(N+C)/K\\) is linear: every cell costs the same however large the tumour is, so the controller still treats a small culture when the drug removes cells. The resistance premium \\(q_R\\,C/K\\) makes a resistant cell cost more than a sensitive one (with \\(q_R=q_B\\) it counts as two cells), because further cisplatin cannot remove it. Sensitive cells are never rewarded: a schedule that leaves more cells in total is worse unless it buys a larger drop in resistant cells. The hard limit on the burden is the cancer barrier at B_cancer. Precedents: tumour burden as a linear cost (Schättler &amp; Ledzewicz 2015) and in-vitro chemotherapy optimised to avoid resistant growth (Carrère 2017).</p>
<p>States are normalised (by K and B_tox) so that numbers of cells do not swamp the toxicity terms. The cancer barrier grows as N + C nears B_cancer and favours treating; the toxicity barrier grows as T nears B_tox and favours stopping; the schedule is the balance of the two, not a hard-coded rule. The relaxed barrier (Feller &amp; Ebenbauer 2017) is the exact log barrier inside the feasible set and a steep quadratic below z = δ, so a start where the bound cannot be met still gives a finite, informative cost; violations are counted in the tables, never hidden. Defaults: q_B = q_R = q_T = 1, r_u = 0.1, μ = 0.05, δ = 0.02; the weights (including q_R = 0, burden only, and q_R = 3) are swept below.</p>
<p><b>Controller: sequential (SQP) MPC.</b> At each decision time: (1) simulate the nonlinear model over the horizon with the current plan; (2) linearise F along it (A_k, B_k by automatic differentiation); (3) gradient by the adjoint recursion \\(\\mu_k=\\nabla\\ell_k+A_k^\\top\\mu_{k+1}\\) (vector–Jacobian products); (4) Gauss–Newton Hessian from the forward sensitivities \\(S_{k+1}=A_kS_k+B_ke_k^\\top\\); (5) solve the box-constrained quadratic subproblem; (6) line search on the true nonlinear cost; repeat to convergence; apply the first dose; re-plan at the next media change. The local subproblem has the LQR structure \\(\\delta x_{k+1}=A_k\\delta x_k+B_k\\delta u_k\\), but it is re-linearised every iteration, bounded, and has non-quadratic barriers, so this is successive-linearisation MPC, not a fixed LQR.</p>
<p><b>Two modes.</b> <i>Data-supported</i>: doses restricted to 0, 0.62, 1.0 and 1.47 µM (the measured doses); the finite-horizon problem is then solved exactly by searching every dose sequence (prefix sharing and pruning), $(LONG ? "4⁷ sequences over each 14-day look-ahead, re-planned at each of the $NRUN media changes (the last plans shrink to end at day $(round(Int, DAYS)))." : "4⁷ sequences over 14 days.") <i>Demonstration</i>: continuous doses up to $(CFG["dose"]["u_max_demo"]) µM with SQP; above 1.47 µM $DRUG_CURVE is extrapolated and the result is not supported by data.</p></div>"""
end

function tests_section()
    f = joinpath(CTRL14, "tests.txt")
    body = isfile(f) ? "<pre style='white-space:pre-wrap;background:#fbfdfd;border:1px solid #e4ecef;border-radius:8px;padding:10px'>$(hesc(read(f, String)))</pre>" : "<p class='muted'>Run test/test_control.jl to fill this in.</p>"
    """<h2>Tests $(tip("From test/test_control.jl. Each check guards a step the controller relies on: that it simulates the fitted models, that the integrator is accurate, that the automatic derivatives are right, that the linearisation is locally exact, that the adjoint gradient is right, and that the exact search really finds the optimum.", "Better: every check passes."))</h2>
<div class='intro'><ul><li>The control plant at 1 µM reproduces the fitted Stage 3 coculture model (relative difference &lt; 10⁻⁵).</li>
<li>Fixed-step RK4 F(x, u) matches an adaptive Tsit5 solve (&lt; 10⁻⁶) at random biologically valid states.</li>
<li>Automatic-differentiation Jacobians A, B match central finite differences (&lt; 10⁻⁵ relative); finite differences are used only here, never in the optimiser.</li>
<li>Linearisation error falls 4× when the perturbation halves (second order).</li>
<li>The adjoint (vector–Jacobian) gradient of the cost matches finite differences.</li>
<li>The exact discrete search equals brute force; SQP never breaks the dose bounds and never increases the cost.</li></ul>$body</div>"""
end

function assumptions_section()
    I = get(AN, :intermittent, nothing)
    wash = I === nothing ? "" : join(["$(l): ΔBIC $(@sprintf("%+.2f", I[l].bic_koff - I[l].bic_equal))" for l in ("A2780", "A2780cis")], "; ")
    rr = joinpath(CTRL14, "resource_reset_test.csv")
    rrows = isfile(rr) ? [trow([td(r.line), td(hesc(r.model)), tn(r.k), tn(@sprintf("%.2f", r.bic)), tn(@sprintf("%+.2f", r.delta_bic)), td(hesc(r.params))]) for r in eachrow(CSV.read(rr, DataFrame))] : String[]
    """<h2>Tests of two model assumptions</h2>
<div class='intro'><p><b>Drug washout.</b> The drug-effect compartments assume the effect decays after the drug is removed at the same rate it built up. The intermittent experiments (1 µM on then off, or off then on, at day 7) test this: a separate washout rate is <b>not</b> supported ($wash; a negative value beyond −3 would be needed). The equal-rate model is kept, and its no-refit predictions of those experiments are on the <a href='../zoo_analyses/validation/report.html'>Validation</a> page.</p>
<p><b>Low resource with media changes.</b> Media are changed every 2 days, so a resource state R that is depleted between changes and reset at each change was compared with the carried constant effect (both a simple exponential depletion and a cell-driven depletion \\(dR/dt=-\\lambda X R/1000\\)). Neither passes the 3-BIC rule, so the constant low-resource effect is kept.</p>
$(isempty(rrows) ? "" : table([("Line", "Cell line (low-resource monoculture)."), ("Model", "Resource-reset candidate."), ("k", "Free parameters."), ("BIC", "BIC of the fit."),
                                ("ΔBIC vs carried", "BIC minus the BIC of the carried Stage 1b model; the rule keeps the carried model unless this is below −3, or within 3 with fewer parameters."), ("Fitted values", "λ (and d where present).")], rrows))</div>"""
end

const BEYOND_LONG = "<li><b>Days 15–$(round(Int, DAYS)) are extrapolation.</b> Every stage was fitted on 14 days of data or fewer, so everything after day 14 follows the fitted equations " *
    "without any data to check them, like the continuous-dose demonstration mode. Most normal-medium runs are near the carrying capacity K by about day 10, so the later days are " *
    "mostly a plateau plus competitive release; resistant takeover under constant dosing is much larger by day $(round(Int, DAYS)) than by day 14. The model has no passaging or " *
    "confluence handling, the low-resource level is held constant, and the bootstrap bands widen after day 14.</li>"
function warnings_section()
    P = scenario_plant(:A); Ln, Lc = P.sys.lines
    ecN, ecC = ec50_of(Ln.s, P.p, Float64), ec50_of(Lc.s, P.p, Float64)
    """<h2>What the data do not identify</h2>
<div class='intro'><ul>
$(LONG ? BEYOND_LONG : "")<li><b>The toxicity limit is inactive at the default settings.</b> Since \\(dT/dt=(1-T)\\big[a_H(1-e^{-B})-r_H\\,T\\big]\\), T can never exceed \\(a_H/r_H\\) = $(fmt(default_tox().a_H / default_tox().r_H)), below the bound B_tox = $(CFG["toxicity"]["B_tox"]), whatever the dose. The toxicity term in the cost still discourages dosing, but the hard limit only matters in the sweep with slower host recovery (r_H = $(minimum(CFG["toxicity"]["r_H_sweep"])) /day, bound \\(a_H/r_H\\) = $(fmt(default_tox().a_H / minimum(CFG["toxicity"]["r_H_sweep"])))). The schedules are therefore not toxicity-constrained in the main runs.</li>
<li><b>How each scenario is validated.</b> Scenario A (the Stage 3a coculture fit) passes leave-one-ratio-out validation, but it was fitted to the treated cocultures, so there is no held-out test of its drug response in coculture. Scenario B is the construction that was tested on held-out data: crowding fitted to the untreated cocultures only, then the treated cocultures predicted with no refit (<a href='../zoo_analyses/validation/report.html'>Validation</a>).</li>
<li><b>Dose response outside 0.62–1.47 µM.</b> The monocultures were treated at three doses and the cocultures at 1 µM only. $(DRUG_HILL ? "EC50 is fixed at each line's measured IC50 ($(fmt(ecN)) µM for A2780, $(fmt(ecC)) µM for A2780cis), not estimated." : "The carried kill is linear in the drug effect, so it never saturates: it cannot show a maximum kill at high doses.") Continuous-dose control above 1.47 µM (demonstration mode) extrapolates $DRUG_CURVE.</li>
<li><b>Resistant kill in coculture.</b> The cocultures were treated at one dose only, so how strongly cisplatin kills A2780cis in coculture rests on the monoculture fit: scenario A uses $(stage3_betaC().free ? "the coculture refit (bounded below by half the monoculture value, so it is never switched off)" : "the monoculture value inherited by the coculture fit") and scenario B the monoculture drug response for both lines; every result is given for both.</li>
<li><b>Toxicity.</b> Nothing in these in-vitro data measures host toxicity. Its structure and values come from de Pillis &amp; Radunskaya (2001); the mapping from medium concentration to systemic burden (s_B) and B_cancer are assumptions. Treat toxicity-limited results as illustrative.</li>
$(LONG ? "" : "<li><b>Beyond 14 days.</b> The 28-day runs extrapolate every stage beyond the data.</li>")
<li><b>Low resource + drug.</b> Stage 4 is fitted to a single 1:1 culture at one dose; its low-resource strengths ε_N, ε_C are identified only for that condition. Its first observation is one day after exposure started; the hidden drug and stress states at that point are simulated once from the candidates' starting values and held fixed during fitting (a small approximation).</li>
<li><b>Parameter uncertainty propagation.</b> The bootstrap combines replicates across stages, but each stage was bootstrapped with earlier-stage values held fixed, so the bands are approximate.</li></ul></div>
<h2>Experiments that would most improve control</h2>
<div class='intro'><ol>
<li><b>Coculture dose series</b> (e.g. 0.25, 0.5, 2 and 4 µM at 50:50, 14 days): identifies the coculture dose response and whether A2780cis is killed in coculture (scenario A vs B), the largest structural uncertainty for control, and lets the shape of the dose response (linear or saturating) be estimated in the condition that is being controlled.</li>
<li><b>Coculture washout / intermittent schedules</b> (1 µM on for 2 or 4 days then off, repeated): tests the drug-effect decay and the regrowth after removal in the controlled system; the monoculture intermittent data only test it for single lines.</li>
<li><b>Replicate-level counts over time.</b> The sensitivity analysis ranks sensitive growth r_N, coculture death d_N and the crowding coefficients as the main drivers of the controlled outcome; more replicates of the untreated and treated cocultures narrow exactly those.</li>
<li><b>Low resource + cisplatin at a second dose.</b> Stage 4 has one condition; a second dose separates the low-resource strengths from a drug × low-resource interaction.</li>
<li><b>A toxicity readout</b> (e.g. viability of a non-cancer line such as fibroblasts under the same schedules) would replace the literature toxicity proxy.</li></ol></div>
<h2>References</h2>
<div class='intro'><ul class='muted'>
<li>de Pillis LG, Radunskaya AE (2001). A mathematical tumor model with immune resistance and drug therapy: an optimal control approach. Comput Math Methods Med 3:79–100. Parameter values as tabulated in Luo et al. (2024) DTR-Bench, arXiv:2405.18610 (Table IX).</li>
<li>Friberg LE et al. (2002). Model of chemotherapy-induced myelosuppression with parameter consistency across drugs. J Clin Oncol 20:4713–4721.</li>
<li>Feller C, Ebenbauer C (2017). Relaxed logarithmic barrier function based model predictive control of linear systems. IEEE Trans Autom Control 62:1223–1238.</li>
<li>Carrère C (2017). Optimization of an in vitro chemotherapy to avoid resistant tumours. J Theor Biol 413:24–33.</li>
<li>Zhang J, Cunningham JJ, Brown JS, Gatenby RA (2017). Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. Nat Commun 8:1816.</li>
<li>Schättler H, Ledzewicz U (2015). Optimal Control for Mathematical Models of Cancer Therapies. Springer, New York.</li>
<li>Carrère C (2017). Optimization of an in vitro chemotherapy to avoid resistant tumours. J Theor Biol 413:24–33.</li>
<li>Strobl MAR et al. (2021). Turnover modulates the need for a cost of resistance in adaptive therapy. Cancer Res 81:1135–1147.</li></ul></div>"""
end

function key_findings()
    out = String[]
    for env in (:normal, :lowres), scen in (:A, :B)
        m = main_of(env, scen, "50:50")
        mp = m.policies["MPC, measured doses only"]; c1 = m.policies["Constant 1 µM"]; mx = m.policies["Maximum measured dose (1.47 µM)"]
        sched = join([@sprintf("%.3g", u) for u in mp.U], ", ")
        push!(out, "<li><b>$(env == :normal ? "Normal medium" : "Low resource") · scenario $scen.</b> MPC (measured doses) chooses [$sched] µM for the $NWORD 2-day intervals: $(LONG ? DL : "day-14") total $(cells(mp.o.total_final)) vs $(cells(c1.o.total_final)) under constant 1 µM and $(cells(mx.o.total_final)) under 1.47 µM, with a cumulative dose of $(@sprintf("%.1f", mp.o.cum_dose)) vs $(@sprintf("%.1f", c1.o.cum_dose)) µM·d; resistant fraction $(@sprintf("%.2f", mp.o.resistant_fraction)) vs $(@sprintf("%.2f", c1.o.resistant_fraction)).</li>")
    end
    "<ul>" * join(out) * "</ul>"
end

const CSS = """body{font:14px/1.5 system-ui,sans-serif;margin:0;background:#f4f7f8;color:#17212b}main{max-width:1250px;margin:auto;padding:24px}h1{color:#124b6e;margin:6px 0}h2{color:#124b6e;margin-top:30px}h3{margin:14px 0 6px}h4{margin:14px 0 4px}
.topnav{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:0 0 16px}.nav{display:inline-flex;align-items:center;min-height:34px;padding:0 11px;border-radius:8px;color:#176b87;background:#eef5f6;text-decoration:none;font:650 14px system-ui,sans-serif;white-space:nowrap}.nav-home{padding:6px 9px}.nav-active{color:#fff;background:#176b87}
.intro,.panel{background:#fff;border:1px solid #d9e2e6;border-radius:12px;padding:16px 18px;margin:14px 0;min-width:0;overflow-wrap:anywhere}.intro{border-left:6px solid #176b87}.note{background:#fff6d9;border-radius:8px;padding:10px 14px}
.eq{overflow-x:auto;overflow-y:hidden;background:#fbfdfd;border:1px solid #e4ecef;border-radius:8px;padding:4px 12px;margin:8px 0}svg{display:block;width:100%;height:auto;margin:10px 0}.tscroll{overflow-x:auto}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:12px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}td,th{border:1px solid #d3dde2;padding:6px 8px;text-align:left;vertical-align:top}th{background:#e9f3f8}td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}tr.win td{background:#eef8f1;font-weight:600}.muted{color:#5b6b73}.flag{color:#b45309;font-weight:600}.ok{color:#2f7d4a;font-weight:600}
code{overflow-wrap:anywhere;word-break:break-word}.tipwrap{position:relative}.tipwrap>.ii{position:absolute;top:8px;right:8px;z-index:2}mjx-container:not([display]){display:inline-block;max-width:100%;overflow-x:auto;vertical-align:middle}"""
# what is better for each table column (info_tips.jl: header_better); unknown columns stop the build
const BETTER = Dict{String,Any}(
    "Policy" => :none,
    "State" => :none,
    "Environment" => :none,
    "Scenario" => :none,
    "Setting" => :none,
    "Scen." => :none,
    "Value" => :none,
    "Schedule (µM)" => :none,
    "Step (days)" => :none,
    "Parameter" => :none,
    "Units" => :none,
    "Source" => :none,
    "Meaning" => :none,
    "Line" => :none,
    "Model" => :none,
    "Fitted values" => :none,
    "N $DL" => :lower,
    "C $DL" => :lower,
    "N + C" => :lower,
    "N + C $DL" => :lower,
    "Resistant fraction" => :lower,
    "Max toxicity" => :lower,
    "Cumulative dose (µM·d)" => "Lower is better when the cell counts are similar.",
    "Dose (µM·d)" => "Lower is better when the cell counts are similar.",
    "Cancer bound hit" => "Better: No / 0% (the bound is never reached).",
    "Toxicity bound hit" => "Better: No / 0% (the bound is never reached).",
    "Draws" => "Higher gives more reliable percentages.",
    "T $DL" => :lower,
    "k" => "Lower is better (fewer parameters), as long as the fit stays within 3 BIC.",
    "BIC" => :lower,
    "ΔBIC vs carried" => "Below −3 would favour the resource-reset model; otherwise the carried model stays.")
const AGREE = "Better: the rows agree; the smallest step is the reference."
hcell(h::AbstractString) = "<th>$h</th>"; hcell(h::Tuple) = "<th>$(h[1])$(tip(h[2], header_better(BETTER, h)))</th>"
table(head, rows) = "<div class='tscroll'><table><tr>" * join(hcell.(head)) * "</tr>" * join(rows) * "</table></div>"
trow(cells; cls="") = "<tr$(isempty(cls) ? "" : " class='$cls'")>" * join(cells) * "</tr>"
td(x) = "<td>$x</td>"; tn(x) = "<td class='num'>$x</td>"
withtip(html, text, better) = "<div class='tipwrap'>$(tip(text, better))$html</div>"

const TITLE = LONG ? "Optimal control of cisplatin dosing over $(round(Int, DAYS)) days" : "Optimal control of cisplatin dosing"
# the long report says it runs past the data; the 14-day report links to the long one when it exists
function intro_extra()
    if LONG
        return "</p><p class='note'><b>Runs to day $(round(Int, DAYS)); the data end at day $(round(Int, DATA_DAYS)).</b> The controller decides at each of the $NRUN media changes, " *
               "looking 14 days ahead each time. Days $(round(Int, DATA_DAYS) + 1)–$(round(Int, DAYS)) are model extrapolation (marked on every time plot); see " *
               "'What the data do not identify' below. The <a href='../control/report.html'>14-day report</a> stays within the data."
    end
    long = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "control_30d", "report.html")
    isfile(long) ? " A <a href='../control_30d/report.html'>30-day version</a> runs the same controller past the end of the data." : ""
end
const INTRO_EXTRA = intro_extra()
function write_control_report()
    home = "<svg width='22' height='22' viewBox='0 0 24 24' aria-hidden='true'><path fill='currentColor' d='M12 3.1 2.8 10.6v10.3h6.1v-6.2h6.2v6.2h6.1V10.6L12 3.1Zm7.1 15.8h-2v-6.2H6.9v6.2h-2V11.6L12 5.8l7.1 5.8v7.3Z'/></svg>"
    D = "../../../../"
    nl(l, h, a) = "<a href='$h' class='nav$(a ? " nav-active" : "")'>$l</a>"
    nav = "<nav class='topnav' aria-label='Report navigation'><a href='$(D)index.html' class='nav nav-home' aria-label='Reports home' title='Reports home'>$home</a>" *
          nl("Summary", "$(D)outputs/results/summary.html", false) * nl("Model zoo fits", "../zoo_fits/report.html", false) *
          nl("Validation", "../zoo_analyses/validation/report.html", false) * (LONG ? nl("Optimal control", "../control/report.html", false) * nl("Optimal control, $(round(Int, DAYS)) days", "report.html", true) :
                    nl("Optimal control", "report.html", true)) * nl("Final report", "$(D)outputs/final_report/report.html", false) * "</nav>"
    html = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>$(TITLE)</title>
<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}};</script><script defer src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>
<style>$CSS</style>$INFO_TIPS_HEAD</head><body><main>$nav
<h1>$(TITLE)</h1>
<div class='intro'><p>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")) by <code>run_control.jl</code> and <code>control_report.jl</code>. Sequential model-predictive control (MPC) chooses the cisplatin concentration for every 2-day interval between media changes, using the carried compartment models: normal medium (Stage 3) and low resource (Stage 4), starting from the 25:75, 50:50 and 75:25 seedings. Two drug-response scenarios are run throughout: <b>A</b> the coculture fit as selected (Stage 3) and <b>B</b> the monoculture drug response for both lines. The control goal is a low total tumour burden with resistant cells held back (a resistant cell counts as two), within a literature-based toxicity limit and a hard cap on the burden; everything it adds to the fitted biology is listed with its source below. Every result is a model simulation: none of these schedules has been tested experimentally yet.$(INTRO_EXTRA)</p>
<h3>Key findings (50:50 start) $(tip("Auto-generated from the runs below. The schedule lists the dose (µM) chosen for each of the $NWORD 2-day intervals.", :none))</h3>$(key_findings())<p>$(hesc(adaptive_finding()))</p></div>
$(model_section())
$(objective_section())
$(env_section(:normal))
$(env_section(:lowres))
$(robust_section())
$(sweep_section())
$(tests_section())
$(assumptions_section())
$(warnings_section())
<p class='muted'>Outputs: <code>control_outcomes.csv</code>, <code>resource_reset_test.csv</code> in <code>Modeling_Approaches/02_mechanical_automatic_package/outputs/$SUB/</code>. Code: control_config.toml, control_model.jl, control_mpc.jl, run_control.jl, control_report.jl, test/test_control.jl.</p>
</main></body></html>"""
    write(joinpath(COUT, "report.html"), stage_names(html)); mkpath(CDOCS); cp(joinpath(COUT, "report.html"), joinpath(CDOCS, "report.html"); force=true)
    println("Wrote $(joinpath(CDOCS, "report.html"))")
end

# pattern of the data-supported MPC schedules, stated from the runs (not assumed)
function adaptive_finding()
    parts = String[]
    for env in (:normal, :lowres), scen in (:A, :B)
        m = main_of(env, scen, "50:50"); mp = m.policies["MPC, measured doses only"]; c1 = m.policies["Constant 1 µM"]
        U = mp.U
        blocks = Tuple{Int,Int}[]                       # runs of treated intervals
        for (i, u) in enumerate(U)
            u > 0 || continue
            (isempty(blocks) || blocks[end][2] != i - 1) ? push!(blocks, (i, i)) : (blocks[end] = (blocks[end][1], i))
        end
        dd(i) = @sprintf("%g", (i - 1) * Δ0)
        kind = isempty(blocks) ? "no treatment" :
               length(U) == blocks[1][2] - blocks[1][1] + 1 ? "treats throughout" :
               "treats on days " * join(["$(dd(a))–$(dd(b + 1))" for (a, b) in blocks], " and ") * " ($(@sprintf("%.1f", mp.o.cum_dose)) vs $(@sprintf("%.1f", c1.o.cum_dose)) µM·d under constant 1 µM)"
        extra = String[]
        # stated against the goal (fewer cells in total, fewer resistant cells); both directions are reported
        push!(extra, "$(mp.o.total_final <= c1.o.total_final ? "fewer" : "more") total cells ($(cells(mp.o.total_final)) vs $(cells(c1.o.total_final)))")
        push!(extra, "$(mp.o.C_final <= c1.o.C_final ? "fewer" : "more") resistant cells ($(cells(mp.o.C_final)) vs $(cells(c1.o.C_final)))")
        push!(parts, "$(env == :normal ? "normal medium" : "low resource") $scen: $kind" * (isempty(extra) ? "" : "; " * join(extra, ", ")))
    end
    "MPC schedule pattern (50:50, not hard-coded): " * join(parts, "; ") * "."
end
function write_control_js()
    js(x::Real) = isfinite(x) ? string(round(Float64(x); sigdigits=6)) : "null"
    jo(o) = "{\"total\":$(js(o.total_final)),\"rf\":$(js(o.resistant_fraction)),\"cum\":$(js(o.cum_dose)),\"maxTox\":$(js(o.max_tox))}"
    esc_js(t) = replace(t, "\\" => "\\\\", "\"" => "\\\"")
    rows = String[]
    for env in (:normal, :lowres), scen in (:A, :B)
        m = main_of(env, scen, "50:50"); P = m.policies
        run = "$(env == :normal ? "Normal medium" : "Low resource") · scenario $scen"
        push!(rows, "{\"run\":\"$(esc_js(run))\",\"schedule\":[$(join(js.(P["MPC, measured doses only"].U), ","))],\"mpc\":$(jo(P["MPC, measured doses only"].o)),\"const1\":$(jo(P["Constant 1 µM"].o)),\"mtd\":$(jo(P["Maximum measured dose (1.47 µM)"].o))}")
    end
    jsvar = LONG ? "CONTROL$(round(Int, DAYS))" : "CONTROL"
    out = "// written by control_report.jl — optimal-control summary for the Summary and Presentation\nwindow.$jsvar = {\"generated\":\"$(Dates.format(now(), "yyyy-mm-dd HH:MM"))\",\"days\":$(round(Int, DAYS)),\"dataDays\":$(round(Int, DATA_DAYS)),\"Btox\":$(CFG["toxicity"]["B_tox"]),\"finding\":\"$(esc_js(adaptive_finding()))\",\"rows\":[$(join(rows, ","))]};\n"
    write(joinpath(ROOT, "docs", "outputs", "results", LONG ? "control$(round(Int, DAYS)).js" : "control.js"), stage_names(out))
end
write_control_report()
write_control_js()
