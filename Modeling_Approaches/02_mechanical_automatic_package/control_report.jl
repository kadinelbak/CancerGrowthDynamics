# HTML report for the optimal-control experiments (run_control.jl checkpoint).
#   julia --project=. control_report.jl
# Writes outputs/control/report.html and docs/Modeling_Approaches/02_mechanical_automatic_package/outputs/control/report.html

include(joinpath(@__DIR__, "control_fasting.jl"))     # includes control_model.jl and control_mpc.jl; lr_plant lives there
const ZOO_DIR = joinpath(@__DIR__, "outputs", "zoo_fits")
const ZOO_DOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "zoo_fits")
const BOOT_B = ZOO_STATE.B
include(joinpath(@__DIR__, "zoo_report.jl"))           # hesc, tip(), INFO_TIPS_HEAD

# Run length: CONTROL_DAYS (default closed_loop_days = 14). Another value renders the dose-only report of
# outputs/control_<days>d/ (run_control.jl with the same setting) into docs/.../outputs/control_<days>d/.
# The long dose + fasting run (long_days) has its own report, control_fasting_report.jl.
# Headline drug-response scenario (control_config.toml [scenarios].headline; A since 2026-11-12, user): listed first
# and labelled "(main result)"; B is the labelled what-if (control_model.jl scenario_role). Both are computed and reported.
const HEADLINE = String(get(get(CFG, "scenarios", Dict{String,Any}()), "headline", "A"))
const SCENS = HEADLINE == "B" ? (:B, :A) : (:A, :B)
scen_label(scen) = "$scen ($(scenario_role(scen)))"
const DAYS = parse(Float64, get(ENV, "CONTROL_DAYS", string(CFG["schedule"]["closed_loop_days"])))
const LONG_D = round(Int, CFG["schedule"]["long_days"])        # the long dose + fasting run (control_<LONG_D>d)
const CL_NAME = "Closed-loop MPC (observes N, C with noise)"   # as in run_control.jl
const CV_OBS = Float64(CFG["observation"]["count_cv"])
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
dover(o) = hasproperty(o, :days_over) ? o.days_over : NaN     # days with N + C ≥ B_cancer
f3(x) = isfinite(x) ? @sprintf("%.3f", x) : "—"
d1(x) = isfinite(x) ? @sprintf("%.1f", x) : "—"
# time to progression (control_mpc.jl: time_to_progression); NaN = not progressed within the run
dayx(x) = isnan(x) ? "not progressed" : @sprintf("%.1f", x)
ttpkey(o) = isnan(o.ttp) ? Inf : o.ttp                        # for ranking: not progressed is best (latest)
const LOOK_TO_END = get(ENV, "CONTROL_LOOKAHEAD_TO_END", string(Int(get(CFG["schedule"], "lookahead_to_end", false)))) == "1"
const LOOKB = Int.(CFG["schedule"]["lookahead_blocks"]); const LOOK_DAYS = round(Int, sum(LOOKB) * Float64(CFG["schedule"]["decision_interval_days"]))
const T_TTP = "Time to progression = the start of the first spell in which N + C stays at or above B_prog = $(CFG["progression"]["threshold_fraction_of_K"]) K for at least $(CFG["progression"]["sustain_days"]) days, or until the end of the run. A brief excursion that treatment brings back under B_prog within that time does not count. 'not progressed' = no such spell within the run (the best outcome)."

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

const T_TOTAL = "Total cells N + C under each dosing policy. The dashed orange line is B_cancer ($(CFG["cancer"]["B_cancer_fraction_of_K"]) of the fitted capacity K), where the soft burden penalty rises steeply. Lower is less cancer burden."
# caveats specific to this page, shown behind ⚠ buttons (the full list is warnings_section)
const W_SIM = "Model simulation: no data were collected for these schedules."
const W_BANDS = "The bootstrap bands are approximate: each stage was bootstrapped with earlier-stage values held at their estimates, and scenario B's crowding (no bootstrap of its own) borrows the Stage 3a bootstrap ratios. The 10% count noise and the state reconstruction in the closed loop are assumptions."
const W_DEMO = "Above 1.47 µM the drug response is extrapolated: the demonstration mode is not supported by data."
const T_NC = "Sensitive (A2780, solid) and resistant (A2780cis, dashed) cells under the data-supported MPC schedule and under constant 1 µM. Whether MPC keeps sensitive cells to hold resistant cells back shows here: the resistant curve under MPC vs constant dosing."
const T_TOX = "Toxicity burden T = 1 − H (fraction of healthy host tissue lost; literature model of de Pillis & Radunskaya 2001). The dashed line is the bound B_tox (0.25: healthy tissue kept above 75%)."
const T_DOSE = "The cisplatin concentration (µM) each policy puts in the medium for every 2-day interval between media changes. Bar height is the dose (hover a bar for its value)."
const T_OUT = "Outcomes on $DL of the model simulation. Resistant fraction = C/(N+C). Max toxicity is the highest T reached; cumulative dose is ∫u dt in µM·days. 'Above B_cancer' counts the soft burden penalty level: whether, and for how many days, N + C was at or above B_cancer. 'Toxicity bound hit' means T touched B_tox at any time."

function outcome_table(m)
    rows = String[]
    for pn in pol_order(m)
        v = m.policies[pn]; o = v.o
        flag = o.tox_violation ? " <span class='flag'>breaks the toxicity limit</span>" : ""
        push!(rows, trow([td("<span style='color:$(polcol(pn))'>■</span> " * hesc(pn) * flag), tn(dayx(o.ttp)), tn(cells(o.total_final)), tn(d1(get(o, :days_half_K, NaN))),
                          tn(cells(o.C_final)), tn(cells(o.N_final)),
                          tn(@sprintf("%.2f", o.resistant_fraction)), tn(@sprintf("%.3f", o.max_tox)), tn(@sprintf("%.1f", o.cum_dose)), td(yn(o.cancer_violation)),
                          tn(d1(dover(o))), td(yn(o.tox_violation)), tn(dayx(o.ttp_first)), tn(d1(o.days_progressed))];
                     cls=startswith(pn, "MPC, measured") ? "win" : ""))
    end
    table([("Policy", "The dosing rule. MPC = model-predictive control; 'measured doses only' restricts it to 0, 0.62, 1 and 1.47 µM, the doses that were measured. A policy that breaks the toxicity limit is flagged in its row."),
           ("Time to progression (day)", T_TTP), ("N + C $DL", "Total cells at $DL.", :lower), ("Days above K/2", "Days with N + C at or above half the shared capacity K.", :lower),
           ("C $DL", "Resistant A2780cis cells at $DL."), ("N $DL", "Sensitive A2780 cells at $DL."),
           ("Resistant fraction", "C / (N + C) at $DL. Lower keeps the population drug-sensitive."), ("Max toxicity", T_TOX),
           ("Cumulative dose (µM·d)", "∫u dt over $(round(Int, DAYS)) days."), ("Above B_cancer", "Whether N + C reached B_cancer (where the soft burden penalty rises steeply) at any time."),
           ("Days above B_cancer", "How many days N + C spent at or above B_cancer."), ("Toxicity bound hit", "Whether T reached B_tox at any time."),
           ("First progression (day)", "The first day N + C reached B_prog (it may be brought back under it)."),
           ("Days progressed", "Days N + C spent at or above B_prog.")], rows) |> x -> withtip(x, T_OUT, "The goal is time to progression: later (or 'not progressed') is better; lower is better for cell counts, resistant fraction, toxicity, dose and days above the thresholds; No is better for the yes/no columns.")
end

# Pareto view: every policy as one point, final cells against cumulative dose (filled: N + C; open: resistant C).
# The dashed line joins the policies no other policy beats on both total cells and dose.
function pareto_svg(title, pts; W=760, H=340)
    xm = max(maximum(q.x for q in pts) * 1.1, 1.0); ym = ceil(maximum(max(q.tot, q.C) for q in pts) * 1.12; sigdigits=2)
    L, R, T, B = 64, 16, 28, 92
    sx(v) = L + v / xm * (W - L - R); sy(v) = H - B - clamp(v, 0, ym) / ym * (H - T - B)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(hesc(title))'><rect width='100%' height='100%' fill='white'/><text x='$L' y='17' font-size='13' font-weight='bold'>$(hesc(title))</text>",
         "<text x='$(W-R)' y='17' font-size='10' text-anchor='end' fill='#52636b'>● N + C   ○ resistant C   - - not beaten on both</text>"]
    for v in range(0, ym; length=5)
        push!(p, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#e5e9ec'/><text x='%d' y='%.1f' font-size='10' text-anchor='end'>%.0f</text>", L, W - R, sy(v), sy(v), L - 6, sy(v) + 3, v))
    end
    for v in range(0, xm; length=6)
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%.1f</text>", sx(v), H - B + 14, v))
    end
    push!(p, "<line x1='$L' y1='$(H-B)' x2='$(W-R)' y2='$(H-B)' stroke='#555'/><text x='$((L+W-R)÷2)' y='$(H-B+28)' font-size='11' text-anchor='middle'>Cumulative dose (µM·d)</text><text transform='translate(13 $((T+H-B)÷2)) rotate(-90)' font-size='11' text-anchor='middle'>Cells at $DL</text>")
    front = sort([q for q in pts if !any(r -> r.x <= q.x && r.tot <= q.tot && (r.x < q.x || r.tot < q.tot), pts)]; by=q -> q.x)
    length(front) > 1 && push!(p, "<polyline points='" * join([@sprintf("%.1f,%.1f", sx(q.x), sy(q.tot)) for q in front], " ") * "' fill='none' stroke='#94a3b8' stroke-dasharray='4,3'/>")
    for q in pts
        push!(p, @sprintf("<circle cx='%.1f' cy='%.1f' r='5.5' fill='%s'><title>%s: N + C = %.0f, dose %.1f µM·d</title></circle>", sx(q.x), sy(q.tot), q.color, hesc(q.label), q.tot, q.x))
        push!(p, @sprintf("<circle cx='%.1f' cy='%.1f' r='4.5' fill='white' stroke='%s' stroke-width='2'><title>%s: C = %.0f, dose %.1f µM·d</title></circle>", sx(q.x), sy(q.C), q.color, hesc(q.label), q.C, q.x))
    end
    for (i, q) in enumerate(pts)
        x = L + ((i - 1) % 3) * 230; y = H - B + 46 + ((i - 1) ÷ 3) * 13
        push!(p, "<circle cx='$(x+6)' cy='$(y-4)' r='4.5' fill='$(q.color)'/><text x='$(x+16)' y='$y' font-size='10'>$(hesc(first(q.label, 38)))</text>")
    end
    push!(p, "</svg>"); join(p)
end
const T_PARETO = "Each policy is one point at its cumulative dose (x) and its cells at $DL (y): filled = total N + C, open = resistant C (50:50 start). " *
    "A point lower and further left is better: fewer cells for less drug. The dashed line joins the policies that no other policy beats on both total cells and dose. " *
    "MPC maximises the time to progression and prices the resistant cells it leaves, traded against toxicity, dose and the soft penalties, so it need not have the fewest cells."

# the stage cost split into its terms over the run (control_mpc.jl: cost_terms), same weights for every policy
function cost_table(env, scen, m)
    P = plant_for_report(env, scen); W = default_weights(P)
    rows = String[]
    for pn in pol_order(m)
        ct = cost_terms(P, m.x0, m.policies[pn].U, Δ0, W)
        push!(rows, trow([td("<span style='color:$(polcol(pn))'>■</span> " * hesc(pn)), tn(f3(ct.progression)), tn(f3(ct.burden)), tn(f3(ct.resistance)), tn(f3(ct.toxicity)), tn(f3(ct.dose)),
                          tn(f3(ct.penalty_burden)), tn(f3(ct.barrier_tox)), tn(f3(ct.total))]; cls=startswith(pn, "MPC, measured") ? "win" : ""))
    end
    table([("Policy", "The dosing rule.", :none), ("Progression", "Σ q_P·σ((N + C − B_prog)/s): about q_P for every interval spent at or above B_prog = $(CFG["progression"]["threshold_fraction_of_K"]) K, the controller's main goal.", :lower), ("Burden", "Σ q_B·(N + C)/K over the $NRUN intervals.", :lower), ("Resistance premium", "Σ q_R·C/K.", :lower),
           ("Toxicity", "Σ q_T·(T/B_tox)².", :lower), ("Dose", "Σ r_u·(u/u_max)², with u_max = $(CFG["dose"]["u_max_data"]) µM for every policy (so a dose above it, in the demonstration mode, costs more than 0.1).", :lower),
           ("Burden penalty", "Σ of the soft burden penalty b_μ(1 − (N + C)/B_cancer): the log barrier below B_cancer and its steep quadratic continuation above it.", :lower),
           ("Toxicity barrier", "Σ b_μ(1 − T/B_tox).", :lower),
           ("Total", "Σ ℓ over the run, the stage cost MPC trades off (the terminal cost is a planning device and is left out).", :lower)], rows) |>
        x -> withtip(x, "Each policy's run scored with the controller's own stage cost, term by term (50:50 start). Large 'Burden penalty' values show where the soft penalty above B_cancer, not the cell count itself, drives the cost. MPC minimises the total over its $(LOOK_DAYS)-day look-ahead at every decision, not each term, and not this run total exactly.", "Lower is better (in every column, and above all in Total).")
end

# best fixed (non-MPC) policy on one outcome: (name, value)
function best_fixed(m, f; higher=false)
    fx = [pn for pn in pol_order(m) if !startswith(pn, "MPC")]
    vals = [f(m.policies[pn].o) for pn in fx]
    i = higher ? argmax(vals) : argmin(vals)
    (name=fx[i], v=vals[i])
end

function env_section(env)
    envname = env == :normal ? "normal medium (10% FBS), Stage 3 model" : "low resource (0.5% FBS), Stage 3 model + low-resource modules and the Stage 4 interaction"
    secs = String[]
    for scen in SCENS
        m = main_of(env, scen, "50:50")
        ts = m.policies["No treatment"].o.ts
        tot = [(label=pn, color=polcol(pn), dash=poldash(pn), t=ts, y=m.policies[pn].o.N .+ m.policies[pn].o.C) for pn in pol_order(m)]
        mp = m.policies["MPC, measured doses only"].o; cst = m.policies["Constant 1 µM"].o
        nc = [(label="MPC: sensitive N", color="#0b6e99", dash="", t=ts, y=mp.N), (label="MPC: resistant C", color="#0b6e99", dash="6,4", t=ts, y=mp.C),
              (label="Constant 1 µM: sensitive N", color="#d95f02", dash="", t=ts, y=cst.N), (label="Constant 1 µM: resistant C", color="#d95f02", dash="6,4", t=ts, y=cst.C)]
        tox = [(label=pn, color=polcol(pn), dash=poldash(pn), t=ts, y=m.policies[pn].o.T) for pn in pol_order(m)]
        W = default_weights(plant_for_report(env, scen))
        ratios = join(["<h4>Seeding ratio $(r)</h4>" * outcome_table(main_of(env, scen, r)) for r in ("25:75", "50:50", "75:25")])
        pts = [(label=pn, color=polcol(pn), x=m.policies[pn].o.cum_dose, tot=m.policies[pn].o.total_final, C=m.policies[pn].o.C_final) for pn in pol_order(m)]
        push!(secs, """<section class='panel' id='$(env)-$(scen)'><h3>Scenario $(hesc(m.name)) <small>· $(scenario_role(scen))</small></h3>$(string(scen) == "B" ? "<p class='muted'>" * hesc(B_NOTE) * "</p>" : "")
<div class='grid'><div>$(withtip(ts_svg("Total N + C, 50:50 start", tot; hline=(m.Bc, "B_cancer"), ymax=ceil(max(m.Bc * 1.15, maximum(maximum(s.y) for s in tot) * 1.05); sigdigits=2)), T_TOTAL, "Lower is better; time above the dashed line is penalised, not forbidden."; warn=[W_SIM, :soft_penalty]))</div>
<div>$(withtip(pareto_svg("Final cells vs cumulative dose, every policy (50:50)", pts), T_PARETO, "Lower and further left is better (fewer cells for less drug)."))</div>
<div>$(withtip(ts_svg("Sensitive and resistant cells", nc), T_NC, "Lower is better, above all for the resistant (dashed) curve."))</div>
<div>$(withtip(dose_svg("Dose schedules (µM per 2-day interval)", [(pn, m.policies[pn].U) for pn in pol_order(m)]), T_DOSE, "Neither on its own: less drug is better when the cell counts come out similar."))</div>
<div>$(withtip(ts_svg("Toxicity burden T", tox; hline=(W.Btox, "B_tox"), ymax=max(W.Btox * 1.3, maximum(maximum(s.y) for s in tox) * 1.1), ylabel="T (fraction)"), T_TOX, "Lower is better, and it must stay under the dashed line."; warn=:tox_proxy))</div></div>
<h4>Cost of each run, term by term (50:50)</h4>$(cost_table(env, scen, m))
$ratios</section>""")
    end
    notes = vcat(LONG ? [EXTRAP_NOTE] : String[], env == :lowres ? [LOWRES_NOTE] : String[], [KILL_LIMIT_NOTE])
    "<h2>Results: $(envname)</h2><div class='intro'><ul>" * join("<li>" .* hesc.(notes) .* "</li>") * "</ul></div>" * join(secs)
end
plant_for_report(env, scen) = env == :normal ? scenario_plant(scen) : lr_plant(scen)

mrange(v) = any(isnan, v) || isempty(v) ? "—" : "$(cells(median(v))) ($(cells(quantile(v, .025)))–$(cells(quantile(v, .975))))"
function robust_section()
    secs = String[]; rrows = String[]; prows = String[]
    for env in (:normal, :lowres), scen in SCENS
        R = C[:robust][(env, scen)]; ts = R.ts
        m = main_of(env, scen, "50:50")
        q(v, pr) = isempty(v) ? NaN : quantile(v, pr)
        bandof(M) = (t=ts, y=[median(M[i, :]) for i in axes(M, 1)], lo=[q(M[i, :], 0.025) for i in axes(M, 1)], hi=[q(M[i, :], 0.975) for i in axes(M, 1)])
        band(pn, f) = bandof(hcat([f(d[pn]) for d in R.draws if get(d, pn, nothing) !== nothing]...))
        ser = [merge(band(pn, o -> o.N .+ o.C), (label=pn * " (median, 95%)", color=polcol(pn), dash="")) for pn in ("Constant 1 µM", "MPC, measured doses only")]
        isempty(R.closed) || push!(ser, merge(bandof(hcat([c.N .+ c.C for c in R.closed]...)), (label="Closed-loop MPC (noisy counts)", color="#0f766e", dash="2,3")))
        serC = [merge(band(pn, o -> o.C), (label=pn * ": resistant C", color=polcol(pn), dash="6,4")) for pn in ("Constant 1 µM", "MPC, measured doses only")]
        push!(secs, "<div>" * withtip(ts_svg("$(env == :normal ? "Normal medium" : "Low resource") · scenario $(scen_label(scen)): N + C over bootstrap parameter sets", vcat(ser, serC); hline=(m.Bc, "B_cancer")),
            "Each schedule (fixed, planned with the best-fit model) re-simulated on $(length(R.draws)) bootstrap parameter sets. Line: median; band: 2.5–97.5 percentiles. A wide band means the fitted parameters leave the outcome uncertain; a band crossing the dashed line means the soft burden penalty level B_cancer is passed under plausible parameters.", "Lower and narrower is better; time above the dashed line is penalised.") * "</div>")
        rowfor(name, ds; cls="") = (tot = [d.tot for d in ds]; cf = [get(d, :Cf, NaN) for d in ds]; nf = [get(d, :Nf, NaN) for d in ds];
            trow([td("$(env == :normal ? "normal" : "low resource")"), td(string(scen)), td(name), tn(length(ds)),
                  tn(@sprintf("%.0f%%", 100 * mean(d.cv for d in ds))), tn(@sprintf("%.0f%%", 100 * mean(d.tv for d in ds))),
                  tn(mrange(nf)), tn(mrange(cf)),
                  tn("$(cells(median(tot))) ($(cells(q(tot, .025)))–$(cells(q(tot, .975))))"), tn(@sprintf("%.1f", median(d.cum for d in ds))),
                  tn(@sprintf("%.0f%%", 100 * mean(!isnan(get(d, :ttp, NaN)) for d in ds))),
                  tn((mt = median([isnan(get(d, :ttp, NaN)) ? Inf : d.ttp for d in ds]); isfinite(mt) ? @sprintf("%.1f", mt) : "not progressed"))]; cls=cls))
        for pn in ROB_POLS
            ds = [d[pn] for d in R.draws if get(d, pn, nothing) !== nothing]; isempty(ds) && continue
            push!(rrows, rowfor(hesc(pn), ds))
        end
        isempty(R.closed) || push!(rrows, rowfor("<b>Closed-loop MPC</b> (re-plans every 2 days from noisy counts; first $(length(R.closed)) draws)", R.closed; cls="win"))
        append!(prows, paired_rows(env, scen, R))
    end
    """<h2>Robustness over bootstrap parameter sets $(tip("The planned schedules replayed on bootstrap parameter sets, open loop and closed loop.", :none; warn=W_BANDS))</h2>
<div class='intro'><p>The schedules above are planned with the best-fit parameters. Here each one is replayed on bootstrap parameter sets: the growth, drug and coculture parameters of the same bootstrap replicate of their own stages. Scenario B's crowding was fitted to the untreated cocultures, which have no bootstrap of their own: each draw borrows the Stage 3a bootstrap ratio (replicate ÷ estimate) of the same parameter and applies it to the scenario-B value.</p>
<p><b>Closed loop.</b> MPC plans with the best-fit model and runs on the bootstrap plant. It sees only the cell counts N and C at each media change, with $(round(Int, 100CV_OBS))% multiplicative measurement noise (log-normal, mean 1; control_config.toml [observation]), and reconstructs every hidden state (drug-effect chains, systemic burden B, toxicity T) by simulating its own model under the doses it has applied. Threshold adaptive therapy re-decides on each plant from the same noisy counts. The comparison is <b>paired</b>: on the first $(get(first(values(C[:robust])), :ncl, 0)) draws every policy meets the same plant and the same noise, and the second table compares the closed loop with each policy draw by draw on exactly those draws.</p></div>
<div class='grid'>$(join(secs))</div>
$(table([("Environment", "Normal medium (Stage 3 model) or low resource (Stage 3 model with the fitted low-resource modules and the Stage 4 drug × low-resource interaction; it equals the Stage 4 fit at 1 µM)."), ("Scenario", "A: the coculture fit as selected (Stage 3); B: monoculture drug response for both lines with crowding from the untreated cocultures."),
         ("Policy", "The dosing rule; its schedule is planned with the best-fit model."), ("Draws", "Bootstrap parameter sets simulated."),
         ("Above B_cancer", "Share of parameter sets in which N + C reaches B_cancer (the soft burden penalty level) at any time."), ("Toxicity bound hit", "Share of parameter sets in which T reaches B_tox."),
         ("N $DL", "Median sensitive cells (95% range) across parameter sets."), ("C $DL", "Median resistant cells (95% range)."),
         ("N + C $DL", "Median total (95% range)."), ("Dose (µM·d)", "Median cumulative dose."),
         ("Progressed at the end", "Share of parameter sets in which the culture has lost control (N + C stays at or above B_prog to the end)."),
         ("Median time to progression (day)", "Median over parameter sets; a set that has not progressed counts as later than any that has. " * T_TTP)], rrows))
<h3>Paired comparison on the closed-loop draws $(tip("Closed-loop MPC minus each policy, computed draw by draw on the same bootstrap plants and the same measurement noise (paired), then summarised: the median difference and the share of draws in which the closed loop does better.", "Negative differences and shares above 50% favour the closed-loop MPC."))</h3>
$(table([("Environment", "Normal medium or low resource.", :none), ("Scen.", "Scenario.", :none), ("Compared with", "The policy the closed-loop MPC is compared with on the same draws.", :none),
         ("Pairs", "Draws where both runs finished.", :higher),
         ("Δ N + C", "Median of (closed-loop − policy) total cells at $DL.", "Negative favours the closed loop."), ("Closed loop fewer N + C", "Share of draws where the closed loop ends with fewer total cells.", :higher),
         ("Δ C", "Median of (closed-loop − policy) resistant cells at $DL.", "Negative favours the closed loop."), ("Closed loop fewer C", "Share of draws where the closed loop ends with fewer resistant cells.", :higher),
         ("Δ dose (µM·d)", "Median of (closed-loop − policy) cumulative dose.", "Negative means the closed loop uses less drug.")], prows))"""
end
const ROB_POLS = ["No treatment", "Constant 1 µM", "Maximum measured dose (1.47 µM)", "Periodic 1 µM, 2 days on / 2 off", "Threshold adaptive (Zhang-style)", "MPC, measured doses only"]
# closed loop minus each policy, draw by draw, on the draws that carry a closed-loop run
function paired_rows(env, scen, R)
    rows = String[]
    P = [d for d in R.draws if get(d, CL_NAME, nothing) !== nothing]
    isempty(P) && return rows
    for pn in ROB_POLS
        pr = [(d[CL_NAME], d[pn]) for d in P if get(d, pn, nothing) !== nothing]; isempty(pr) && continue
        dt = [a.tot - b.tot for (a, b) in pr]; dc = [a.Cf - b.Cf for (a, b) in pr]; dd = [a.cum - b.cum for (a, b) in pr]
        push!(rows, trow([td("$(env == :normal ? "normal" : "low resource")"), td(string(scen)), td(hesc(pn) * (pn == "MPC, measured doses only" ? " (open loop: the nominal schedule replayed)" : "")), tn(length(pr)),
                          tn(@sprintf("%+.0f", median(dt))), tn(@sprintf("%.0f%%", 100 * mean(dt .< 0))), tn(@sprintf("%+.0f", median(dc))), tn(@sprintf("%.0f%%", 100 * mean(dc .< 0))),
                          tn(@sprintf("%+.1f", median(dd)))]))
    end
    rows
end

function sweep_section()
    rows = String[]
    for s in C[:sweeps]
        s.kind == "burden penalty" && continue            # its own table (penalty_section)
        o = s.res.o
        push!(rows, trow([td(hesc(s.kind)), td(string(s.scen)), td(hesc(s.label)), td(join([@sprintf("%.3g", u) for u in s.res.U], " ")),
                          tn(cells(o.N_final)), tn(cells(o.C_final)), tn(cells(o.total_final)), tn(@sprintf("%.2f", o.resistant_fraction)), tn(@sprintf("%.3f", o.max_tox)), tn(@sprintf("%.1f", o.cum_dose)),
                          tn(dayx(o.ttp))]))
    end
    srows = [trow([td(string(s.scen)), tn(string(s.h)), tn(cells(s.N)), tn(cells(s.C)), tn(@sprintf("%.5f", s.T))]) for s in C[:steps]]
    tb = [s for s in C[:sweeps] if s.kind == "tie-break scale" && string(s.scen) == HEADLINE]
    tbnote = isempty(tb) ? "" : "<div class='intro'><p><b>Tie-break scale (scenario $HEADLINE, 50:50).</b> The second tier is scaled by tiebreak_scale = $(tiebreak_scale()) (default). " *
        (all(s -> isnan(s.res.o.ttp), tb) ? "The MPC does not progress within $(round(Int, DAYS)) days at any of these scales (only fixed comparators do), so the second tier ranks its schedules, and the schedule still moves with the scale: " :
                                             "The second tier ranks the MPC schedules that do not progress, and the schedule moves with the scale: ") *
        join(["$(hesc(s.label)): [$(join([@sprintf("%.3g", u) for u in s.res.U], ", "))] µM, N + C $(cells(s.res.o.total_final)), resistant $(cells(s.res.o.C_final))" for s in tb], "; ") *
        ". Between 0.1 and 0.3 the schedules differ by at most one dose level per interval and the totals by a few percent, but the resistant counts by up to about 10%: the 14-day schedule is a weighting choice, not a property of the data (user, 2026-11-12).</p></div>"
    """<h2>Sensitivity of the controller</h2>
$tbnote
<div class='intro'><p>Data-supported MPC in normal medium from the 50:50 start, varying one setting at a time: the cost weights, the toxicity model (host recovery rate r_H, dose mapping s_B, bound B_tox), the planning horizon, the decision interval and the simulation length (the soft burden penalty and the terminal weight are in their own table, under the objective). If the schedule barely changes across plausible weights, the recommendation does not hinge on an arbitrary weight choice.</p></div>
$(table([("Setting", "Which assumption is varied."), ("Scen.", "Drug-response scenario (A or B)."), ("Value", "The value used; everything else is at its default."),
         ("Schedule (µM)", "The dose applied in each 2-day interval (1-day intervals for the 'every 1 day' row" * (LONG ? ")." : "; 14 intervals for the 28-day row).")),
         ("N $DL", "Sensitive cells at the end."), ("C $DL", "Resistant cells at the end."), ("N + C", "Total at the end."),
         ("Resistant fraction", "C/(N + C) at the end."), ("Max toxicity", "Highest T."), ("Dose (µM·d)", "Cumulative dose."), ("Time to progression (day)", T_TTP)], rows))
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
    "The report varies s_B, the limit, and r_H down to 0.2/day (about a 5-day recovery, the order of the " *
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
             prow("\\alpha_{NC}", get(v, "α_NC", 1.0), "—", haskey(v, "α_NC") ? "DATA · Stage 3" : "fixed", "crowding of N by C"), prow("\\alpha_{CN}", get(v, "α_CN", 1.0), "—", haskey(v, "α_CN") ? "DATA · Stage 3" : "fixed", "crowding of C by N"),
             (haskey(v, "d_N") ? prow("d_N", v["d_N"], "/day", "DATA · Stage 3", "sensitive death in coculture") : ""), (haskey(v, "β_N") ? prow("\\beta_N", v["β_N"], "/day", "DATA · Stage 3", "max sensitive kill (scenario A; B uses Stage 2)") :
              prow("\\beta_N", c3.sys.lines[1].s.β.v, "/day", "DATA · Stage 2", "max sensitive kill, kept from the monoculture in coculture (user, 2026-11-03)")),
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
<div class='intro'><p>The cell equations are exactly the carried compartment models, with the same parameter values: the control plant simulates a flattened, type-stable copy of the fits' right-hand side (tests check it against the zoo's own right-hand side to 10⁻¹², and that the plant reproduces the fitted Stage 3 curves at 1 µM to 10⁻⁵). Every run starts from the nominal seeding (100 cells split by the ratio, e.g. 25:75 → N = 25, C = 75), not from a measured day-0 count. The controller only adds the dose input \\(u\\) and a toxicity proxy. Drug → effect → kill follows \\(u \\to A_{X,1} \\to A_{X,2} \\to A_{X,3}=E_X \\to\\) kill, and separately \\(u \\to B \\to T\\).</p>
$(join(["<div class='eq'>$(i == length(eqs) ? tip(TOX_TIP, "Lower T is better, and it must stay under the 0.25 bound (at least 75% of healthy tissue kept)."; warn=[:tox_proxy, "It is a generic host model, not calibrated to cisplatin or to these cells; s_B (1 µM in the medium = 1 unit of burden) is an assumption."]) : "")\\[$e\\]</div>" for (i, e) in enumerate(eqs)]))
<p><b>Low resource (D = u, u<sub>R</sub> = $(U_R)).</b> $(tip("The low-resource equations as the zoo writes them for the carried Stage 4 model: the cisplatin module and the low-resource module of each line. The drug-effect compartments and the toxicity states B, T are as above.", :none)) $(isempty(fast_plant(:A).iε) ? "Stage 4 carries a drug × low-resource interaction: each line's cisplatin kill is multiplied by \\((1+\\chi\\,u_R/0.955)\\), χ = $(fmt(fast_plant(:A).spec[1].chi)), so with no drug the model is the low-resource fit alone and at 1 µM it is exactly the Stage 4 fit." : "The low-resource strengths depend on the dose, \\(\\varepsilon_X(u)=1+(\\varepsilon_X^{S4}-1)\\min(u/1\\,\\mu\\mathrm{M},1)\\): with no drug the model is the low-resource fit alone (\\(\\varepsilon=1\\)), and at 1 µM it is exactly the Stage 4 fit.") This is the same plant as the $(LONG_D)-day dose + fasting controller held in low resource throughout, so the two reports agree.</p>
$(join(["<div class='eq'>\\[$e\\]</div>" for e in eqlr]))
<p>Discrete transition over one decision interval (Δ = $(Δ0) days, u constant): \\(x_{k+1}=F(x_k,u_k)\\), fixed-step RK4 with step 0.05 day. Local linearisation \\(A_k=\\partial F/\\partial x\\), \\(B_k=\\partial F/\\partial u\\) by forward-mode automatic differentiation through the RK4 steps (one Jacobian of F over [x; u]: nx + 1 dual directions, about the cost of reverse mode at this size, with no tape).</p></div>
<div class='grid'><div>$(table([("State", "Symbol."), ("Units", "Units."), ("Meaning", "What it represents.")], srows))</div>
<div>$(withtip(flow_svg(), "How parameters flow between stages: each stage fits only its new parameters and inherits the rest. The control plant inherits everything and adds the dose and the toxicity proxy.", :none))</div></div>
$(table([("Parameter", "Symbol in the equations."), ("Value", "Value used by the controller (scenario A, normal medium)."), ("Units", "Units."),
         ("Source", "DATA = fitted to these experiments (stage named); LITERATURE = taken from a cited model; ASSUMPTION = exploratory, swept below, to be replaced by measurements."),
         ("Meaning", "What the parameter does.")], prows))"""
end

function objective_section()
    """<h2>Objective, soft penalties and barriers</h2>
<div class='intro'>
<div class='eq'>\\[\\min_{u_1,\\dots,u_H}\\ \\sum_{k=1}^{H}\\ell(x_{k+1},u_k)+\\phi(x_{H+1}),\\qquad 0\\le u_k\\le u_{max}\\]</div>
<div class='eq'>\\[\\ell=q_P\\,\\sigma\\!\\Big(\\tfrac{N+C-B_{prog}}{s}\\Big)+q_B\\,\\tfrac{N+C}{K}+q_R\\,\\tfrac{C}{K}+q_T\\Big(\\tfrac{T}{B_{tox}}\\Big)^2+r_u\\Big(\\tfrac{u}{u_{max}}\\Big)^2+b_{\\mu_C}\\!\\Big(1-\\tfrac{N+C}{B_{cancer}}\\Big)+b_{\\mu_T}\\!\\Big(1-\\tfrac{T}{B_{tox}}\\Big),\\qquad \\sigma(z)=\\tfrac{1}{1+e^{-z}}\\]</div>
<div class='eq'>\\[b_\\mu(z)=\\begin{cases}-\\mu\\log z, & z\\ge\\delta\\\\ \\mu\\left[\\tfrac12\\left(\\tfrac{z-2\\delta}{\\delta}\\right)^2-\\tfrac12-\\log\\delta\\right], & z<\\delta\\end{cases}\\qquad \\phi=$(CFG["weights"]["terminal_factor"])\\left[q_P\\,\\sigma(\\cdot)+q_B\\,\\tfrac{N+C}{K}+q_R\\,\\tfrac{C}{K}+q_T\\Big(\\tfrac{T}{B_{tox}}\\Big)^2\\right]+q_{RT}\\,\\tfrac{C}{K}\\]</div>
<p><b>Goal: time to progression</b> (user, 2026-10-07). Progression is N + C at or above B_prog = $(CFG["progression"]["threshold_fraction_of_K"]) K. The step σ (width s = $(CFG["progression"]["step_width_fraction_of_K"]) K) is ≈ 0 below B_prog and ≈ 1 above it, so with q_P = $(CFG["progression"]["qP"]) the summed stage cost counts the intervals spent progressed: the controller holds the culture under B_prog as long as it can, treating hard when it nears the threshold and holding off when it does not, the quantity adaptive therapy maximises (Zhang et al. 2017; Gatenby et al. 2009). At the end of the look-ahead \\(q_{RT}\\,C/K\\) (q_RT = $(CFG["progression"]["qRT"])) prices the resistant cells left behind: they drive the next progression, which no dose can reverse. Progression is either line over its own threshold (rule "$(CFG["progression"]["rule"])", user 2026-11-03): resistant cells at $(CFG["progression"]["resistant_threshold_fraction_of_K"]) K or sensitive cells at $(CFG["progression"]["sensitive_threshold_fraction"]) × their untreated plateau. Everything else forms a second tier at $(CFG["progression"]["tiebreak_scale"]) × its weight (user, 2026-11-12): the total burden (q_B = $(CFG["weights"]["qB"])), the resistant cells (q_R = $(CFG["weights"]["qR"])), toxicity, dose and the inside of the host barriers; crossing a host limit keeps the full penalty. The objective without progression (q_P = q_RT = 0) and the 2026-11-01 weights are rows of the weight sweep. Time to progression is the start of the first spell above B_prog that lasts at least $(CFG["progression"]["sustain_days"]) days (or to the end of the run): a brief excursion brought back under the threshold does not count, and one dip after a long spell above does not restart the clock.</p>
$(stopping_note())
<p><b>Previous goal (2026-10-02 to 2026-10-07).</b> Keep the total tumour burden low <i>and</i> hold resistant cells back. The burden term \\(q_B(N+C)/K\\) is linear: every cell costs the same however large the tumour is, so the controller still treats a small culture when the drug removes cells. The resistance premium \\(q_R\\,C/K\\) makes a resistant cell cost more than a sensitive one (with \\(q_R=q_B\\) it counts as two cells), because further cisplatin cannot remove it. Sensitive cells are never rewarded: a schedule that leaves more cells in total is worse unless it buys a larger drop in resistant cells. On top of the linear burden there is a <b>soft burden penalty</b> (a relaxed barrier above B_cancer = $(CFG["cancer"]["B_cancer_fraction_of_K"]) K), not a hard cap. MPC therefore minimises this weighted cost (burden, resistance premium, toxicity, dose and the soft penalties), not the burden alone, so a fixed policy can end with fewer cells; each results section compares MPC with the best fixed policy on every metric and plots all policies on final cells against dose. Precedents: tumour burden as a linear cost (Schättler &amp; Ledzewicz 2015) and in-vitro chemotherapy optimised to avoid resistant growth (Carrère 2017).</p>
<p>States are normalised (by K and B_tox) so that numbers of cells do not swamp the toxicity terms. The burden penalty grows as N + C nears B_cancer and favours treating; the toxicity barrier grows as T nears B_tox and favours stopping; the schedule is the balance of the two, not a hard-coded rule. The relaxed barrier (Feller &amp; Ebenbauer 2017) is the exact log barrier inside the limit and a steep quadratic below z = δ, so a state beyond the limit still gives a finite, informative cost; exceedances are counted in the tables (whether and for how many days), never hidden. Defaults: q_P = $(CFG["progression"]["qP"]), q_RT = $(CFG["progression"]["qRT"]), q_B = $(CFG["weights"]["qB"]), q_R = $(CFG["weights"]["qR"]), q_T = $(CFG["weights"]["qT"]), r_u = $(CFG["weights"]["ru"]), μ = 0.05, δ = 0.02, terminal factor $(CFG["weights"]["terminal_factor"]); the weights (including the previous objective, no terminal resistant cost and a 3× one) are swept below.</p>
<p><b>Why a penalty and not a hard cap.</b> In normal medium the culture grows towards its capacity K under the measured doses: $(breach_summary()). As a hard state constraint, every schedule that cannot hold the culture below B_cancer would be infeasible and the optimiser would have no cost to rank them by; as a penalty, they are still ranked, and the exceedance is reported. A relaxed log barrier is, by construction, a soft-constraint (penalty) formulation (Feller &amp; Ebenbauer 2017), and in cancer optimal control a tumour-burden limit is usually a penalised state constraint or a progression threshold rather than a hard bound (e.g. adaptive-therapy thresholds set relative to the burden, Zhang et al. 2017). B_cancer is therefore reported as the level where a soft penalty starts to dominate. Because its quadratic branch above B_cancer can outweigh the burden term itself (see the cost tables in the results), the table below re-plans the 50:50 MPC with B_cancer at $(join(filter(>(0), CFG["cancer"]["B_cancer_sweep"]), ", ")) K and with no burden penalty at all, each with terminal factor 1, 5 and 10; the outcomes are always scored against the default B_cancer.</p>
$(penalty_table())
<p><b>Controller: sequential (SQP) MPC.</b> At each decision time: (1) simulate the nonlinear model over the horizon with the current plan; (2) linearise F along it (A_k, B_k by automatic differentiation); (3) gradient by the adjoint recursion \\(\\mu_k=\\nabla\\ell_k+A_k^\\top\\mu_{k+1}\\) (vector–Jacobian products); (4) Gauss–Newton Hessian from the forward sensitivities \\(S_{k+1}=A_kS_k+B_ke_k^\\top\\); (5) solve the box-constrained quadratic subproblem; (6) line search on the true nonlinear cost; repeat to convergence; apply the first dose; re-plan at the next media change. The local subproblem has the LQR structure \\(\\delta x_{k+1}=A_k\\delta x_k+B_k\\delta u_k\\), but it is re-linearised every iteration, bounded, and has non-quadratic barriers, so this is successive-linearisation MPC, not a fixed LQR.</p>
<p><b>Two modes.</b> <i>Data-supported</i>: doses restricted to 0, 0.62, 1.0 and 1.47 µM (the measured doses); the finite-horizon problem is then solved exactly by searching every dose sequence (prefix sharing, a warm-start bound and pruning) at each of the $NRUN media changes. The $(LOOK_DAYS)-day look-ahead is split into $(length(LOOKB)) blocks of $(join(LOOKB, ", ")) intervals with one dose each (<i>move blocking</i>, Cagienard et al. 2007), so $(4^length(LOOKB)) sequences are searched; the first block is a single interval, so the dose actually applied is decided at full resolution. <i>Demonstration</i>: continuous doses up to $(CFG["dose"]["u_max_demo"]) µM with SQP; above 1.47 µM $DRUG_CURVE is extrapolated $(warn(W_DEMO))</p>
$(LOOK_TO_END ? "<p><b>Look-ahead ends with the run</b> (user, 2026-10-27). Each decision plans up to $(LOOK_DAYS) days ahead but never past day $(round(Int, DAYS)), so the last decisions know the run is ending and can spend dose on a final kill; the resistant cells left at the end are still priced (q_RT). From 2026-10-07 to 2026-10-27 every decision planned the full $(LOOK_DAYS) days, past the end of the run; both are in the control sweep." : "<p><b>Fixed look-ahead.</b> Every decision, including the last ones, plans $(LOOK_DAYS) days ahead (14 before 2026-10-07; both are in the horizon sweep); near the end of the run the plan simply runs past day $(round(Int, DAYS)) (the model is simulated beyond it). The schedule shown is therefore the first $NRUN decisions of a continuing policy, read off at day $(round(Int, DAYS)), not a plan that stops caring about what happens after its last day (a shrinking horizon stops dosing near the end, because later regrowth no longer counts).")</p></div>"""
end

# how often the soft burden penalty level is passed in normal medium (every ratio, scenario and policy)
function breach_summary()
    ms = [m for m in C[:main] if m.env == :normal]
    n = sum(length(m.policies) for m in ms); k = sum(count(v -> v.o.cancer_violation, collect(values(m.policies))) for m in ms)
    mp = [m.policies["MPC, measured doses only"].o for m in ms]
    "$k of $n normal-medium runs (3 ratios × 2 scenarios × $(length(first(ms).policies)) policies) pass B_cancer = $(CFG["cancer"]["B_cancer_fraction_of_K"]) K at some time, " *
    "the data-supported MPC in $(count(o -> o.cancer_violation, mp)) of $(length(mp)) (for $(d1(median(dover.(mp)))) days, median)"
end

# the soft burden penalty × terminal weight sweep (run_control.jl sweeps, kind "burden penalty")
function penalty_table()
    S = [s for s in C[:sweeps] if s.kind == "burden penalty"]
    isempty(S) && return "<p class='muted'>Burden-penalty sweep not in this checkpoint (re-run run_control.jl).</p>"
    rows = [trow([td(string(s.scen)), td(hesc(s.label)), td(join([@sprintf("%.3g", u) for u in s.res.U], " ")), tn(cells(s.res.o.total_final)), tn(cells(s.res.o.C_final)),
                  tn(@sprintf("%.1f", s.res.o.cum_dose)), tn(d1(dover(s.res.o))), tn(dayx(s.res.o.ttp))];
                 cls=(occursin("B_cancer = $(CFG["cancer"]["B_cancer_fraction_of_K"]) K", s.label) && endswith(s.label, "×$(Float64(CFG["weights"]["terminal_factor"]))")) ? "win" : "") for s in S]
    table([("Scen.", "Drug-response scenario.", :none), ("Setting", "Where the soft burden penalty starts (fraction of K; 'no burden penalty' removes it) and the terminal factor. The highlighted row is the default.", :none),
           ("Schedule (µM)", "Dose in each 2-day interval.", :none), ("N + C $DL", "Total cells at the end.", :lower), ("C $DL", "Resistant cells at the end.", :lower),
           ("Dose (µM·d)", "Cumulative dose."), ("Time to progression (day)", T_TTP), ("Days above B_cancer", "Days with N + C at or above the DEFAULT B_cancer ($(CFG["cancer"]["B_cancer_fraction_of_K"]) K), whatever the setting.")], rows) |>
        x -> withtip(x, "Data-supported MPC from the 50:50 start in normal medium, re-planned with each setting. If the schedule and outcomes barely change, the results do not hinge on where the soft penalty sits or on the terminal weight.", "Lower cell counts are better; rows that agree mean the conclusion is robust to these choices.")
end

function tests_section()
    f = joinpath(CTRL14, "tests.txt")
    body = isfile(f) ? "<pre style='white-space:pre-wrap;background:#fbfdfd;border:1px solid #e4ecef;border-radius:8px;padding:10px'>$(hesc(read(f, String)))</pre>" : "<p class='muted'>Run test/test_control.jl to fill this in.</p>"
    """<h2>Tests $(tip("From test/test_control.jl. Each check guards a step the controller relies on: that it simulates the fitted models, that the integrator is accurate, that the automatic derivatives are right, that the linearisation is locally exact, that the adjoint gradient is right, and that the exact search really finds the optimum.", "Better: every check passes."))</h2>
<div class='intro'><ul><li>The control plant at 1 µM reproduces the fitted Stage 3 coculture model (scenario A) and the held-out construction of scenario B (Ce0 crowding with the Stage 2 kills) on the 50:50 treated coculture (relative difference &lt; 10⁻⁵); the flattened right-hand side equals the zoo's own to 10⁻¹².</li>
<li>Fixed-step RK4 F(x, u) matches an adaptive Tsit5 solve (&lt; 10⁻⁵ relative; the Norton–Simon factor has a kink at capacity, where RK4 loses its smooth-case order) at random biologically valid states.</li>
<li>Automatic-differentiation Jacobians A, B match central finite differences (&lt; 10⁻⁴ relative, as the differences straddle the same kink); finite differences are used only here, never in the optimiser.</li>
<li>Linearisation error falls 4× when the perturbation halves (second order).</li>
<li>The adjoint (vector–Jacobian) gradient of the cost matches finite differences (&lt; 10⁻⁴ relative).</li>
<li>Every stage and terminal cost is ≥ 0 at random states (what makes pruning in the exact search safe), and the cost split by term adds up to the cost.</li>
<li>The exact discrete search equals brute force; SQP never breaks the dose bounds and never increases the cost; at the SQP solution the projected gradient is near zero, and SQP over [0, 1.47] µM is not worse than the exact optimum over the four measured doses (a subset).</li>
<li>The nominal closed loop with the fixed look-ahead equals replaying its own schedule, and each of its decisions equals the first move of the exact plan from the replayed state.</li></ul>$body</div>"""
end

function assumptions_section()
    I = get(AN, :intermittent, nothing)
    wash = I === nothing ? "" : join(["$(l): ΔBIC $(@sprintf("%+.2f", I[l].bic_koff - I[l].bic_equal))" for l in ("A2780", "A2780cis")], "; ")
    rr = joinpath(CTRL14, "resource_reset_test.csv")
    rrows = isfile(rr) ? [trow([td(r.line), td(hesc(r.model)), tn(r.k), tn(@sprintf("%.2f", r.bic)), tn(@sprintf("%+.2f", r.delta_bic)), td(hesc(r.params))]) for r in eachrow(CSV.read(rr, DataFrame))] : String[]
    """<h2>Tests of two model assumptions</h2>
<div class='intro'><p><b>Drug washout.</b> The drug-effect compartments assume the effect decays after the drug is removed at the same rate it built up. The intermittent experiments (1 µM on then off, or off then on, at day 7) test this: a separate washout rate is <b>not</b> supported ($wash; a negative value beyond −2 would be needed). The equal-rate model is kept, and its no-refit predictions of those experiments are on the <a href='../zoo_analyses/validation/report.html'>Validation</a> page.</p>
<p><b>Low resource with media changes.</b> Media are changed every 2 days, so a resource state R that is depleted between changes and reset at each change was compared with the carried constant effect (both a simple exponential depletion and a cell-driven depletion \\(dR/dt=-\\lambda X R/1000\\)). Neither passes the 2-BIC rule, so the constant low-resource effect is kept.</p>
$(isempty(rrows) ? "" : table([("Line", "Cell line (low-resource monoculture)."), ("Model", "Resource-reset candidate."), ("k", "Free parameters."), ("BIC", "BIC of the fit."),
                                ("ΔBIC vs carried", "BIC minus the BIC of the carried Stage 1b model; the rule keeps the carried model unless this is below −2, or within 2 with fewer parameters."), ("Fitted values", "λ (and d where present).")], rrows))</div>"""
end

const BEYOND_LONG = "<li><b>Days 15–$(round(Int, DAYS)) are extrapolation.</b> Every stage was fitted on 14 days of data or fewer, so everything after day 14 follows the fitted equations " *
    "without any data to check them, like the continuous-dose demonstration mode. Most normal-medium runs are near the carrying capacity K by about day 10, so the later days are " *
    "mostly a plateau plus competitive release; resistant takeover under constant dosing is much larger by day $(round(Int, DAYS)) than by day 14. The model has no passaging or " *
    "confluence handling, the low-resource level is held constant, and the bootstrap bands widen after day 14.</li>"
function warnings_section()
    P = scenario_plant(:A); Ln, Lc = P.sys.lines
    ecN, ecC = ec50_of(Ln.s, P.p, Float64), ec50_of(Lc.s, P.p, Float64)
    """<h2 id='limits'>What the data do not identify $(tip("Every caveat on this page in full; the ⚠ buttons next to the figures and tables show the one that applies there.", :none))</h2>
<div class='intro'><ul>
$(LONG ? BEYOND_LONG : "")<li><b>The toxicity limit is inactive at the default settings.</b> Since \\(dT/dt=(1-T)\\big[a_H(1-e^{-B})-r_H\\,T\\big]\\), T can never exceed \\(a_H/r_H\\) = $(fmt(default_tox().a_H / default_tox().r_H)), below the bound B_tox = $(CFG["toxicity"]["B_tox"]), whatever the dose. The toxicity term in the cost still discourages dosing, but the toxicity barrier only comes close to binding in the sweep with slower host recovery (r_H = $(minimum(CFG["toxicity"]["r_H_sweep"])) /day, bound \\(a_H/r_H\\) = $(fmt(default_tox().a_H / minimum(CFG["toxicity"]["r_H_sweep"])))). The schedules are therefore not toxicity-constrained in the main runs.</li>
<li><b>How each scenario is validated.</b> Scenario A (the Stage 3a coculture fit, the main result) reproduces the treated cocultures, since it was fitted to them; leave-one-ratio-out validation tests it on held-out ratios, and the plant-vs-data table above checks the controller's plant against the measured day-14 counts. Scenario B was built without the treated cocultures (Ce0 crowding, monoculture kills) and predicts them poorly: $(hesc(B_NOTE)) (<a href='../zoo_analyses/validation/report.html'>Validation</a>).</li>
<li><b>Dose response outside 0.62–1.47 µM.</b> The monocultures were treated at three doses and the cocultures at 1 µM only. $(DRUG_HILL ? "EC50 is fixed at each line's measured IC50 ($(fmt(ecN)) µM for A2780, $(fmt(ecC)) µM for A2780cis), not estimated." : "The carried kill is linear in the drug effect, so it never saturates: it cannot show a maximum kill at high doses.") Continuous-dose control above 1.47 µM (demonstration mode) extrapolates $DRUG_CURVE.</li>
<li><b>Resistant kill in coculture.</b> $(hesc(KILL_LIMIT_NOTE)) The coculture fit finds no detectable resistant kill at 1 µM (φ_C ≈ 0), so in scenario A the controller cannot kill resistant cells at any dose; scenario B shows what the controller would do if the resistant cells kept their monoculture kill. Both are reported; A is the headline (user, 2026-11-12; B was the headline 2026-11-02 to 2026-11-12).</li>
<li><b>Toxicity.</b> Nothing in these in-vitro data measures host toxicity. Its structure and values come from de Pillis &amp; Radunskaya (2001); the mapping from medium concentration to systemic burden (s_B) and B_cancer are assumptions. Treat toxicity-limited results as illustrative.</li>
$(LONG ? "" : "<li><b>Beyond 14 days.</b> The 28-day runs extrapolate every stage beyond the data.</li>")
<li><b>Low resource + drug.</b> Stage 4 is fitted to a single 1:1 culture at one dose; its low-resource strengths ε_N, ε_C are identified only for that condition. Its first observation is one day after exposure started; the hidden drug and stress states at that point are simulated once from the candidates' starting values and held fixed during fitting (a small approximation).</li>
<li><b>Parameter uncertainty propagation.</b> The bootstrap combines replicates across stages, but each stage was bootstrapped with earlier-stage values held fixed, so the bands are approximate. Scenario B's crowding (fitted to the untreated cocultures) has no bootstrap of its own: its draws borrow the Stage 3a bootstrap ratios (replicate ÷ estimate) of the same parameters.</li>
<li><b>Measurement and state estimation.</b> The closed-loop runs assume each media change gives a count of N and C with $(round(Int, 100CV_OBS))% noise and that the hidden drug-effect and toxicity states follow the planner's own model; a real loop needs the two lines counted separately at each media change, as the coculture experiments did.</li></ul></div>
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
<li>Strobl MAR et al. (2021). Turnover modulates the need for a cost of resistance in adaptive therapy. Cancer Res 81:1135–1147.</li></ul></div>"""
end

# MPC against the BEST fixed policy on each metric (each metric's own best), and whether any fixed policy
# beats MPC on total cells, resistant cells and dose at once
function key_findings()
    out = String[]
    for env in (:normal, :lowres), scen in SCENS
        m = main_of(env, scen, "50:50")
        mp = m.policies["MPC, measured doses only"]; o = mp.o
        sched = join([@sprintf("%.3g", u) for u in mp.U], ", ")
        bt = best_fixed(m, x -> x.total_final); bc = best_fixed(m, x -> x.C_final); bd = best_fixed(m, x -> x.cum_dose); bx = best_fixed(m, x -> x.max_tox)
        dom = [pn for pn in pol_order(m) if !startswith(pn, "MPC") && (q = m.policies[pn].o; q.total_final <= o.total_final && q.C_final <= o.C_final && q.cum_dose <= o.cum_dose)]
        push!(out, "<li><b>$(env == :normal ? "Normal medium" : "Low resource") · scenario $(scen_label(scen)).</b> MPC (measured doses) chooses [$sched] µM for the $NWORD 2-day intervals. " *
                   "Against the best fixed policy on each metric at $DL: total cells $(cells(o.total_final)) vs $(cells(bt.v)) ($(hesc(bt.name))); resistant cells $(cells(o.C_final)) vs $(cells(bc.v)) ($(hesc(bc.name))); " *
                   "cumulative dose $(@sprintf("%.1f", o.cum_dose)) vs $(@sprintf("%.1f", bd.v)) µM·d ($(hesc(bd.name))); peak toxicity $(@sprintf("%.3f", o.max_tox)) vs $(@sprintf("%.3f", bx.v)) ($(hesc(bx.name))); " *
                   (bp = best_fixed(m, ttpkey; higher=true); "time to progression $(dayx(o.ttp)) vs $(isfinite(bp.v) ? @sprintf("%.1f", bp.v) : "not progressed") ($(hesc(bp.name))). ") *
                   (isempty(dom) ? "No fixed policy has at once no more total cells, no more resistant cells and no more drug than MPC." :
                                   "Fixed policies with no more total cells, resistant cells and drug than MPC: $(hesc(join(dom, ", "))).") * "</li>")
    end
    "<ul>" * join(out) * "</ul><p>MPC maximises the time to progression over its $(LOOK_DAYS)-day look-ahead (and prices the resistant cells it leaves), traded against toxicity and dose, so a fixed policy can beat it on a single end-of-run metric; the cost tables and the final-cells-vs-dose plots in each results section show the trade-off.</p>"
end

# Drug-response shape: the main runs (the carried Stage 2b drug model) next to a Hill alternative run with
# ZOO_STATE / CONTROL_OUT overrides into outputs/sensitivity/hill/control (when it exists)
function hill_section()
    f = joinpath(@__DIR__, "outputs", "sensitivity", "hill", "control", "control_outcomes.csv")
    fm = joinpath(CTRL14, "control_outcomes.csv")
    shape = DRUG_HILL ? "Hill" : "linear" * (DRUG_NS ? " (Norton–Simon)" : "")
    head = "<h2>Drug-response shape (Hill vs linear) $(tip("The carried drug model is Stage 2b's selection ($shape kill). Here the 14-day results of a run whose zoo state uses the Hill (saturating) kill instead are set next to the main ones, to show how much the control conclusions rest on the shape of the dose response, which the data (three monoculture doses, one coculture dose) barely identify.", :none))</h2>"
    (isfile(f) && isfile(fm)) || return head * "<div class='intro'><p class='muted'>Not run: no Hill alternative at outputs/sensitivity/hill/control/control_outcomes.csv (run run_control.jl with ZOO_STATE pointing to the Hill zoo state and CONTROL_OUT=outputs/sensitivity/hill).</p></div>"
    A = CSV.read(f, DataFrame); M = CSV.read(fm, DataFrame); rows = String[]
    pick(D, env, scen, pn) = (r = D[(string.(D.environment) .== env) .& (string.(D.scenario) .== scen) .& (D.ratio .== "50:50") .& (D.policy .== pn), :]; nrow(r) == 1 ? r[1, :] : nothing)
    for env in ("normal", "lowres"), scen in string.(SCENS), pn in ("MPC, measured doses only", "Constant 1 µM")
        a = pick(M, env, scen, pn); b = pick(A, env, scen, pn)
        (a === nothing || b === nothing) && continue
        push!(rows, trow([td(env == "normal" ? "normal" : "low resource"), td(scen), td(hesc(pn)), tn(cells(a.total_final)), tn(cells(b.total_final)), tn(cells(a.C_final)), tn(cells(b.C_final)),
                          tn(@sprintf("%.1f", a.cumulative_dose_uM_days)), tn(@sprintf("%.1f", b.cumulative_dose_uM_days)), td(hesc(a.schedule_uM)), td(hesc(b.schedule_uM))]; cls=startswith(pn, "MPC") ? "win" : ""))
    end
    head * "<div class='intro'><p>50:50 start, day 14. 'Main' is the carried $shape drug model; 'Hill' is the same pipeline with the Hill (saturating) kill.</p>" *
        table([("Environment", "Normal medium or low resource.", :none), ("Scen.", "Scenario.", :none), ("Policy", "Data-supported MPC or constant 1 µM.", :none),
               ("N + C, main", "Total cells at day 14, main drug model.", :lower), ("N + C, Hill", "Total cells at day 14, Hill alternative.", :lower),
               ("C, main", "Resistant cells, main.", :lower), ("C, Hill", "Resistant cells, Hill alternative.", :lower),
               ("Dose, main (µM·d)", "Cumulative dose, main.", "Lower is better when the cell counts are similar."), ("Dose, Hill (µM·d)", "Cumulative dose, Hill alternative.", "Lower is better when the cell counts are similar."),
               ("Schedule, main", "MPC doses per 2-day interval (µM), main.", :none), ("Schedule, Hill", "The same with the Hill alternative.", :none)], rows) * "</div>"
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
    "Above B_cancer" => "Better: No / 0% (the soft burden penalty level is never reached).",
    "Days above B_cancer" => :lower,
    "First progression (day)" => "Later is better; 'not progressed' (never reached B_prog) is best.",
    "Time to progression (day)" => "Later is better; 'not progressed' is best: the goal of the controller.",
    "Median time to progression (day)" => "Later is better; 'not progressed' is best.",
    "Days progressed" => :lower,
    "Progressed at the end" => :lower,
    "Toxicity bound hit" => "Better: No / 0% (the bound is never reached).",
    "Draws" => "Higher gives more reliable percentages.",
    "T $DL" => :lower,
    "k" => "Lower is better (fewer parameters), as long as the fit stays within 2 BIC.",
    "BIC" => :lower,
    "ΔBIC vs carried" => "Below −2 would favour the resource-reset model; otherwise the carried model stays.")
const AGREE = "Better: the rows agree; the smallest step is the reference."
hcell(h::AbstractString) = "<th>$h</th>"; hcell(h::Tuple) = "<th>$(h[1])$(tip(h[2], header_better(BETTER, h)))</th>"
table(head, rows) = "<div class='tscroll'><table><tr>" * join(hcell.(head)) * "</tr>" * join(rows) * "</table></div>"
trow(cells; cls="") = "<tr$(isempty(cls) ? "" : " class='$cls'")>" * join(cells) * "</tr>"
td(x) = "<td>$x</td>"; tn(x) = "<td class='num'>$x</td>"
# withtip(html, text, better; warn=...) lives in info_tips.jl

const TITLE = LONG ? "Optimal control of cisplatin dosing over $(round(Int, DAYS)) days" : "Optimal control of cisplatin dosing"
# the long report says it runs past the data; the 14-day report links to the long one when it exists
function intro_extra()
    if LONG
        return "</p><p class='note'><b>Runs to day $(round(Int, DAYS)); the data end at day $(round(Int, DATA_DAYS)).</b> The controller decides at each of the $NRUN media changes, " *
               "looking $(LOOK_DAYS) days ahead each time. Days $(round(Int, DATA_DAYS) + 1)–$(round(Int, DAYS)) are model extrapolation (marked on every time plot); see " *
               "'What the data do not identify' below. The <a href='../control/report.html'>14-day report</a> stays within the data."
    end
    long = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "control_$(LONG_D)d", "report.html")
    isfile(long) ? " A <a href='../control_$(LONG_D)d/report.html'>$(LONG_D)-day version</a> runs the same controller past the end of the data, with fasting (low resource) as a second control." : ""
end
const INTRO_EXTRA = intro_extra()
const PVD = plant_vs_data()
const CLOSER = "Closer to the measured value (the next column) is better; within about 2 × the fit RMSE is within the data's scatter."
# How the 14-day schedule is ranked (user, 2026-11-12), stated from the runs: when no policy progresses within the
# run, time to progression cannot separate them and the second tier decides
function tier_note()
    prog(scen) = [pn for r in ("25:75", "50:50", "75:25") for (pn, v) in main_of(:normal, Symbol(HEADLINE), r).policies if !isnan(v.o.ttp)]
    p = prog(HEADLINE)
    txt = isempty(p) ? "Scenario $HEADLINE (main result): nothing progresses within $(round(Int, DAYS)) days from these starts, under any policy, so the $(round(Int, DAYS))-day schedule is ranked by the second tier of the objective (total cells, resistant cells, dose and host costs); time to progression separates the policies only in the longer runs (30, 60 and $(LONG_D) days)." :
          "Scenario $HEADLINE (main result): $(length(unique(p))) policies progress within $(round(Int, DAYS)) days ($(join(unique(p), ", "))); time to progression ranks the others against them, the second tier (total cells, resistant cells, dose and host costs) breaks the ties."
    "<div class='intro'><p><b>How to read the $(round(Int, DAYS))-day results.</b> $(hesc(txt))</p></div>"
end
const B_NOTE = b_whatif_note(PVD)
LONG || CSV.write(joinpath(CTRL14, "plant_vs_data.csv"), DataFrame([merge(r, (scenario=string(r.scenario),)) for r in PVD]))
# The controller's plant against the coculture data (user, 2026-11-12): each condition simulated from its own
# fitted start, untreated against Ce0 and at a constant 1 µM against Ce1, day-14 N and C
function plant_data_section()
    rows = [trow([td("$(r.scenario) ($(hesc(scenario_role(r.scenario))))"), td(r.ratio), td(hesc(r.condition)), tn(@sprintf("%.0f / %.0f", r.N0, r.C0)),
                  tn(cells(r.N)), tn(cells(r.N_data)), tn(cells(r.C)), tn(cells(r.C_data)), tn(cells(r.N + r.C)), tn(cells(r.N_data + r.C_data)),
                  tn(@sprintf("%.0f / %.0f", r.rmseN, r.rmseC))]) for r in PVD]
    """<h2 id='plantdata'>The controller's plant against the data (day 14)</h2>
<div class='intro'><p>Every control run starts from the Stage 3a fitted starting state of its ratio's treated (Ce1) series (control_config [initial] source = "$(initial_source())"), so the controller simulates the same wells as the fit. Here each condition is simulated from its own fitted start: untreated against the untreated coculture (Ce0) and at a constant 1 µM against the treated coculture (Ce1). $(hesc(B_NOTE))</p>
$(table([("Scenario", "A: the Stage 3a coculture fit; B: the monoculture drug response in coculture (what-if).", :none), ("Ratio", "Seeded sensitive : resistant.", :none),
         ("Condition", "Ce0: no drug; Ce1: a constant 1 µM.", :none), ("Start N / C", "The fitted starting counts of that series.", :none),
         ("N plant", "Sensitive cells at day 14, control plant.", CLOSER), ("N data", "Measured sensitive cells at day 14.", :none),
         ("C plant", "Resistant cells at day 14, control plant.", CLOSER), ("C data", "Measured resistant cells at day 14.", :none),
         ("N + C plant", "Total at day 14, control plant.", CLOSER), ("N + C data", "Measured total at day 14.", :none),
         ("Fit RMSE N / C", "The Stage 3a fit's residual RMSE over that series (cells): the scale of the data's scatter around the fit.", :none)], rows))</div>"""
end
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
<div class='intro'><p>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")) by <code>run_control.jl</code> and <code>control_report.jl</code>. Sequential model-predictive control (MPC) chooses the cisplatin concentration for every 2-day interval between media changes, using the carried compartment models: normal medium (Stage 3) and low resource (Stage 4), starting from the nominal 25:75, 50:50 and 75:25 seedings (100 cells split by the ratio). Two drug-response scenarios are run throughout: <b>A</b> the coculture fit as selected (Stage 3, resistant kill β_C = 0) and <b>B</b> the monoculture drug response for both lines. <b>Scenario $HEADLINE is the headline (main result)</b> (user, 2026-11-02) because resistant cells being completely inactive to the drug is not biologically plausible; $(HEADLINE == "B" ? "scenario A is kept as a sensitivity case: it brackets the weaker resistant kill (if the resistant kill is as weak as the coculture fit suggests, MPC waits and doses late)." : "scenario B is kept as a sensitivity case.") The control goal is the longest time to progression: keep N + C under B_prog = $(CFG["progression"]["threshold_fraction_of_K"]) K for as long as possible, with the resistant cells left at the end of each look-ahead priced, a literature-based toxicity barrier and a soft burden penalty above B_cancer. $(LOOK_TO_END ? "Every decision looks up to $(LOOK_DAYS) days ahead, never past the end of the run." : "Every decision looks $(LOOK_DAYS) days ahead, also the last ones, so each schedule is the first $NRUN decisions of a continuing policy.") Everything the controller adds to the fitted biology is listed with its source below. $(warn([W_SIM, "None of these schedules has been tested experimentally yet."]))$(INTRO_EXTRA)</p>
<h3>Key findings (50:50 start) $(tip("Auto-generated from the runs below. The schedule lists the dose (µM) chosen for each of the $NWORD 2-day intervals; each metric is compared with the fixed policy that is best on that metric.", :none))</h3>$(key_findings())<p>$(hesc(adaptive_finding()))</p></div>
$(LONG ? "" : tier_note())
$(model_section())
$(LONG ? "" : plant_data_section())
$(objective_section())
$(env_section(:normal))
$(env_section(:lowres))
$(robust_section())
$(sweep_section())
$(hill_section())
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
    for env in (:normal, :lowres), scen in SCENS
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
        push!(extra, "time to progression $(dayx(mp.o.ttp)) vs $(dayx(c1.o.ttp)) days")
        push!(parts, "$(env == :normal ? "normal medium" : "low resource") $scen: $kind" * (isempty(extra) ? "" : "; " * join(extra, ", ")))
    end
    "MPC schedule pattern (50:50, not hard-coded): " * join(parts, "; ") * "."
end
function write_control_js()
    js(x::Real) = isfinite(x) ? string(round(Float64(x); sigdigits=6)) : "null"
    jo(o) = "{\"total\":$(js(o.total_final)),\"rf\":$(js(o.resistant_fraction)),\"cum\":$(js(o.cum_dose)),\"maxTox\":$(js(o.max_tox)),\"ttp\":$(js(get(o, :ttp, NaN))),\"ttpFirst\":$(js(get(o, :ttp_first, NaN)))}"
    esc_js(t) = replace(t, "\\" => "\\\\", "\"" => "\\\"")
    rows = String[]
    for env in (:normal, :lowres), scen in SCENS
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
