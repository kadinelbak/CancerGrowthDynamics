# HTML report for the long-run optimal control with two inputs, cisplatin and fasting (run_control_fasting.jl).
#   julia --project=. control_fasting_report.jl          (CONTROL_DAYS overrides long_days, as in the runner)
# Writes outputs/control_<days>d/report.html, docs/.../outputs/control_<days>d/report.html and
# docs/outputs/results/controlLong.js (window.CONTROL_LONG, for the Summary and Presentation).
# One continuing run is read off at the snapshot days (control_config.toml): the look-ahead is fixed (lookahead_blocks, 28 days)
# at every decision (it runs past the end of the run), so a shorter run would make exactly the same decisions
# up to its end; there is no end-of-run effect.

include(joinpath(@__DIR__, "control_fasting.jl"))
const ZOO_DIR = joinpath(@__DIR__, "outputs", "zoo_fits")
const ZOO_DOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "zoo_fits")
const BOOT_B = ZOO_STATE.B
include(joinpath(@__DIR__, "zoo_report.jl"))           # hesc, fmt, tip(text, better), INFO_TIPS_HEAD

const DAYS = parse(Float64, get(ENV, "CONTROL_DAYS", string(CFG["schedule"]["long_days"]))); const DATA_DAYS = Float64(CFG["schedule"]["data_days"])
const DSTR = string(round(Int, DAYS))
const SNAP = [d for d in Float64.(CFG["schedule"]["snapshot_days"]) if d <= DAYS]
const TICK = DAYS <= 30 ? 5 : DAYS <= 60 ? 10 : 20
const FOUT = joinpath(@__DIR__, "outputs", control_subdir(DAYS, NaN))
const FDOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "control_$(DSTR)d")
FC = deserialize(joinpath(FOUT, "fasting_state.jls"))
const Δ0 = Float64(CFG["schedule"]["decision_interval_days"]); const NRUN = round(Int, DAYS / Δ0)
const UMAX_DEMO = Float64(CFG["dose"]["u_max_demo"])
const MPC2 = "MPC, dose + fasting (measured levels)"; const MPC1 = "MPC, dose only (fed)"
const MPCD = "MPC demonstration (continuous u ≤ $(UMAX_DEMO) µM, v ∈ [0, 1])"
const B_NOTE = b_whatif_note()
const POL_ORDER = ["No treatment", "Constant 1 µM", "Maximum measured dose (1.47 µM)", "Periodic 1 µM, 2 days on / 2 off", "Threshold adaptive (Zhang-style)",
                   "Fasting only, 48 h every 6 days", "Periodic 1 µM with a 48-h fast on each dosing interval", "Constant 1 µM + 48-h fast every 6 days", "Constant 1 µM + 48-h fast every 8 days", MPC1, MPC2, MPCD]
const POL_COL = Dict("No treatment" => "#6b7280", "Constant 1 µM" => "#d95f02", "Maximum measured dose (1.47 µM)" => "#b91c1c",
                     "Periodic 1 µM, 2 days on / 2 off" => "#8b5cf6", "Threshold adaptive (Zhang-style)" => "#3b8b4f",
                     "Fasting only, 48 h every 6 days" => "#a16207", "Periodic 1 µM with a 48-h fast on each dosing interval" => "#db2777",
                     "Constant 1 µM + 48-h fast every 6 days" => "#9a3412", "Constant 1 µM + 48-h fast every 8 days" => "#c2410c", MPC1 => "#64748b", MPC2 => "#0b6e99", MPCD => "#0f766e")
polcol(n) = get(POL_COL, n, "#0f766e"); poldash(n) = n == MPCD ? "6,4" : n == MPC1 ? "2,3" : ""
cells(x) = isfinite(x) ? @sprintf("%.0f", x) : "—"
yn(b) = b ? "<span class='flag'>yes</span>" : "no"
pct(x) = @sprintf("%.0f%%", 100x)
dover(o) = hasproperty(o, :days_over) ? o.days_over : NaN     # days with N + C ≥ B_cancer
# time to progression (control_mpc.jl: time_to_progression); NaN = not progressed within the run
dayx(x) = isnan(x) ? "not progressed" : @sprintf("%.1f", x)
ttpkey(o) = isnan(o.ttp) ? Inf : o.ttp
const LOOKB = Int.(CFG["schedule"]["lookahead_blocks"]); const LOOK_DAYS = round(Int, sum(LOOKB) * Float64(CFG["schedule"]["decision_interval_days"]))
const T_TTP = "Time to progression = the start of the first spell in which N + C stays at or above B_prog = $(CFG["progression"]["threshold_fraction_of_K"]) K (80% confluence) for at least $(CFG["progression"]["sustain_days"]) days, or until the end of the run; a brief excursion that treatment reverses within that time does not count. 'not progressed' = no such spell within the run (the best outcome)."
f3(x) = isfinite(x) ? @sprintf("%.3f", x) : "—"
d0(x) = isfinite(x) ? @sprintf("%.0f", x) : "—"
const CL_NAME = "Closed-loop MPC, dose + fasting (observes N, C with noise)"   # as in run_control_fasting.jl
const CV_OBS = Float64(CFG["observation"]["count_cv"])
const BC_FRAC = CFG["cancer"]["B_cancer_fraction_of_K"]
# the carried drug model (Stage 2b), for the drug-response-shape section
const DRUG_HILL = zgroup("s2-a2780").cands[1].sys.lines[1].action == :hill
const DRUG_SHAPE = DRUG_HILL ? "Hill" : "linear" * (ns_on(zgroup("s2-a2780").cands[1].sys.lines[1]) ? " (Norton–Simon)" : "")
# headline drug-response scenario (control_config.toml [scenarios].headline, user 2026-11-02): listed first, marked "(main result)"
const HEADLINE = String(get(get(CFG, "scenarios", Dict{String,Any}()), "headline", "A"))
const SCENS = HEADLINE == "B" ? (:B, :A) : (:A, :B)
scen_tag(scen) = scenario_role(scen)       # control_model.jl (A main result since 2026-11-12, user; B the what-if)
main_of(scen, ratio) = only(filter(m -> m.scen == scen && m.ratio == ratio, FC[:main]))
pol_order(m) = [p for p in POL_ORDER if haskey(m.policies, p)]
# withtip(html, text, better; warn=...) lives in info_tips.jl
# caveats specific to this page, shown behind ⚠ buttons (the full list is limits_section at the end)
const W_UPPER = "Fasting's benefit is an upper bound: the tumour-side effect is the in-vitro effect of 0.5% FBS at full strength, with no in-vivo attenuation (a patient's fast changes glucose, insulin and IGF-1 far less than a 20-fold serum cut)."
const W_HOST = "The host states (toxicity, weight, protection) are literature-based proxies: nothing in these experiments measured them."
const W_UNTESTED = "These are model simulations: no schedule on this page has been tested in an experiment."
const W_WMIN = "The weight limit W_min is an assumption informed by the > 5% cachexia criterion (Fearon et al. 2011), which is defined over 6 months, not for an acute fast; acute fasting loss is largely water and glycogen and is regained. It is swept (0.90–0.98)."
const W_PH = "The protection strength p_H is an assumption: the randomized DIRECT trial (de Groot et al. 2020) found less DNA damage in T cells but no difference in toxicity grades, so p_H = 0 (no protection) is in the sweep."
const W_BORROW = "Scenario B's crowding has no bootstrap of its own: each draw borrows the Stage 3a bootstrap ratio of the same parameter, so its bands are approximate. The 10% count noise and the state reconstruction in the closed loop are assumptions."
hcell(h::AbstractString) = "<th>$h</th>"; hcell(h::Tuple) = "<th>$(h[1])$(tip(h[2], h[3]; warn=length(h) > 3 ? h[4] : nothing))</th>"   # (label, tip, better[, warning])
table(head, rows) = "<div class='tscroll'><table><tr>" * join(hcell.(head)) * "</tr>" * join(rows) * "</table></div>"
trow(cs; cls="") = "<tr$(isempty(cls) ? "" : " class='$cls'")>" * join(cs) * "</tr>"
td(x) = "<td>$x</td>"; tn(x) = "<td class='num'>$x</td>"
eq(t) = "<div class='eq'>\\[$t\\]</div>"
eqt(t, text; warn=nothing) = "<div class='eqrow'>$(eq(t))$(tip(text, :none; warn=warn))</div>"

# ── plots ──────────────────────────────────────────────────────────────────────
function ts_svg(title, series; ymax=nothing, hline=nothing, W=760, H=330, ylabel="Mean Cells", ymin=0.0)
    ym = ymax === nothing ? ceil(maximum(maximum(s.y) for s in series) * 1.12; sigdigits=2) : ymax
    ym <= ymin && (ym = ymin + 1.0)
    L, R, T, B = 64, 16, 28, 72
    sx(v) = L + v / DAYS * (W - L - R); sy(v) = H - B - (clamp(v, ymin, ym) - ymin) / (ym - ymin) * (H - T - B)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(hesc(title))'><rect width='100%' height='100%' fill='white'/><text x='$L' y='17' font-size='13' font-weight='bold'>$(hesc(title))</text>"]
    for v in range(ymin, ym; length=5)
        push!(p, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#e5e9ec'/><text x='%d' y='%.1f' font-size='10' text-anchor='end'>%s</text>", L, W - R, sy(v), sy(v), L - 6, sy(v) + 3, ym < 5 ? @sprintf("%.2f", v) : @sprintf("%.0f", v)))
    end
    for v in 0:TICK:floor(Int, DAYS)
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%d</text>", sx(v), H - B + 14, v))
    end
    for d in SNAP          # snapshot days (where the tables read the run off)
        d < DAYS && d != DATA_DAYS && push!(p, @sprintf("<line x1='%.1f' x2='%.1f' y1='%d' y2='%d' stroke='#94a3b8' stroke-dasharray='2,3'/>", sx(d), sx(d), T, H - B))
    end
    push!(p, @sprintf("<rect x='%.1f' y='%d' width='%.1f' height='%d' fill='#f5f3ff'/><line x1='%.1f' x2='%.1f' y1='%d' y2='%d' stroke='#7c3aed' stroke-dasharray='4,4'/><text x='%.1f' y='%d' font-size='10' fill='#7c3aed'>beyond the data (day > %d)</text>",
                      sx(DATA_DAYS), T, sx(DAYS) - sx(DATA_DAYS), H - B - T, sx(DATA_DAYS), sx(DATA_DAYS), T, H - B, sx(DATA_DAYS) + 4, T + 11, round(Int, DATA_DAYS)))
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
        x = L + ((i - 1) % 3) * 230; y = H - 30 + ((i - 1) ÷ 3) * 13
        push!(p, "<line x1='$x' x2='$(x+18)' y1='$(y-4)' y2='$(y-4)' stroke='$(s.color)' stroke-width='2'$(isempty(s.dash) ? "" : " stroke-dasharray='$(s.dash)'")/><text x='$(x+22)' y='$y' font-size='10'>$(hesc(first(s.label, 40)))</text>")
    end
    push!(p, "</svg>"); join(p)
end
# schedules: dose as bar height, fasting intervals shaded purple with an F
function sched_svg(title, pols; W=760, umax=3.0)
    rowH = 28; H = 34 + rowH * length(pols) + 22; L = 260; R = 16
    sx(t) = L + t / DAYS * (W - L - R)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(hesc(title))'><rect width='100%' height='100%' fill='white'/><text x='8' y='17' font-size='13' font-weight='bold'>$(hesc(title))</text>",
         @sprintf("<rect x='%.1f' y='24' width='%.1f' height='%d' fill='#f5f3ff'/>", sx(DATA_DAYS), sx(DAYS) - sx(DATA_DAYS), rowH * length(pols) + 6)]
    for (i, (name, U)) in enumerate(pols)
        y0 = 28 + rowH * (i - 1)
        push!(p, "<text x='8' y='$(y0+17)' font-size='10'>$(hesc(first(name, 46)))</text>")
        for k in axes(U, 1)
            x0, x1 = sx((k - 1) * Δ0), sx(k * Δ0) - 1
            U[k, 2] > 0 && push!(p, @sprintf("<rect x='%.1f' y='%d' width='%.1f' height='%d' fill='#c4b5fd' opacity='%.2f'><title>days %g–%g: fasting (v = %.2f)</title></rect><text x='%.1f' y='%d' font-size='9' fill='#5b21b6' text-anchor='middle'>F</text>",
                                          x0, y0 + 2, x1 - x0, rowH - 6, 0.35 + 0.65U[k, 2], (k - 1) * Δ0, k * Δ0, U[k, 2], (x0 + x1) / 2, y0 + 11))
            U[k, 1] > 0 || continue
            h = (rowH - 8) * min(U[k, 1] / umax, 1.0)
            push!(p, @sprintf("<rect x='%.1f' y='%.1f' width='%.1f' height='%.1f' fill='%s'><title>days %g–%g: %.2f µM</title></rect>", x0, y0 + rowH - 4 - h, x1 - x0, h, polcol(name), (k - 1) * Δ0, k * Δ0, U[k, 1]))
        end
        push!(p, "<line x1='$L' x2='$(W-R)' y1='$(y0+rowH-4)' y2='$(y0+rowH-4)' stroke='#d3dde2'/>")
    end
    for v in 0:TICK:round(Int, DAYS)
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%d</text>", sx(v), H - 6, v))
    end
    push!(p, "</svg>"); join(p)
end

# ── model section ─────────────────────────────────────────────────────────────
ctex(a, X, Y, nm) = (a.i == 0 && a.v == 0) ? X : (a.i == 0 && a.v == 1) ? "$X+$Y" : "$X+\\alpha_{$nm}$Y"
# the de Pillis host equations (B and the healthy-tissue loss T), shown when [toxicity] model = "depillis"
host_tox_eqs_depillis() = """$(eqt(raw"\frac{dB}{dt}=s_B\,u-d_B\,B", "Systemic drug burden (de Pillis & Radunskaya 2001). s_B maps 1 µM in the medium to 1 unit of burden."; warn=:tox_proxy))
$(eqt(raw"\frac{dT}{dt}=a_H\,\left(1-p_H\,P\right)\left(1-e^{-B}\right)(1-T)-r_H\,T\,(1-T)", "Toxicity burden T = fraction of healthy tissue lost (de Pillis & Radunskaya 2001), with the fasting protection factor (1 − p_H·P): while protected, the same drug burden does less damage to healthy tissue. This is differential stress resistance (Raffaghello et al. 2008; de Groot et al. 2020 saw less chemotherapy-induced DNA damage in T cells)."; warn=W_PH))"""
function model_section()
    FP = fast_plant(:A); S = FP.base.sys
    celltex = vcat(line_tex(S.lines[1], "N", ctex(S.alpha[1][2], "N", "C", "NC"), "K"), line_tex(S.lines[2], "C", ctex(S.alpha[2][1], "C", "N", "CN"), "K"))
    h = default_host(); t = default_tox(); W = weights2(FP)
    spec = FP.spec
    mN = spec[1].mods; mC = spec[2].mods
    # each line's cisplatin module (Hill or linear kill, whichever Stage 2b selected) and its low-resource modules
    drug(ms) = only(filter(m -> m.input_D, ms))
    hillN = drug(mN); hillC = drug(mC); hill = hillN.action == :hill
    PB = fast_plant(:B)
    # each line's low-resource module, whatever Stage 2a selected: one row per parameter it uses
    kerneldesc = Dict(:none => "immediate", :comp1 => "effect compartment", :chain2 => "2-compartment chain", :chain3 => "3-compartment chain",
                      :second => "second-order delay", :auc => "cumulative exposure (no fading)")
    function lrrows(ms, X, line)
        rows = String[]
        for m in filter(m -> !m.input_D, ms)
            what = "$(line) low-resource module: $(kerneldesc[m.kernel]) → $(m.action == :capacity ? "lowers K" : "extra death")"
            m.action == :capacity && push!(rows, prow("c_$X", m.c, "—", "DATA · Stage 1b", what * " (capacity reduction per unit effect)"))
            m.action == :loss && m.d0 != 0 && push!(rows, prow("d_$X^{LR}", m.d0, "/day", "DATA · Stage 1b", what * " (immediate death per unit u_R)"))
            m.action == :loss && m.β != 0 && push!(rows, prow("d_$X^{acc}", m.β, "/day", "DATA · Stage 1b", what * " (death per unit accumulated effect)"))
            m.kernel_chain > 0 && push!(rows, prow("k_$X^R", m.k, "/day", "DATA · Stage 1b", what * " (transit rate)"))
            m.kernel == :second && push!(rows, prow("\\omega_$X^R,\\ \\zeta_$X^R", "$(fmt(m.ω)), $(fmt(m.ζ))", "/day, —", "DATA · Stage 1b", what * " (natural frequency, damping)"))
        end
        rows
    end
    chiN, chiC = spec[1].chi, spec[2].chi
    kname = hill ? "maximum" : "per unit drug effect"
    nsnote = spec[1].ns ? "; scaled by the Norton–Simon growth factor" : ""
    prow(sym, val, unit, src, note) = trow([td("\\($sym\\)"), tn(val isa Number ? fmt(val) : val), td(unit), td(src), td(hesc(note))])
    prows = [prow("r_N", spec[1].r, "/day", "DATA · Stage 1", "A2780 growth rate"), prow("r_C", spec[2].r, "/day", "DATA · Stage 1", "A2780cis growth rate"),
             prow("\\theta_C", spec[2].θ, "—", "DATA · Stage 1", "A2780cis crowding shape"),
             prow("K", spec[1].K, "cells", "DATA · Stage 3 (A) / untreated cocultures (B: $(cells(PB.spec[1].K)))", "shared carrying capacity"),
             prow("\\alpha_{NC}", spec[1].α, "—", "DATA · Stage 3", "crowding of A2780 by A2780cis"), prow("\\alpha_{CN}", spec[2].α, "—", "DATA · Stage 3", "crowding of A2780cis by A2780"),
             (spec[1].d != 0 ? prow("d_N", spec[1].d, "/day", "DATA · Stage 3", "A2780 death in coculture") : ""),
             prow("\\beta_N", hillN.β, "/day", "DATA · Stage 3 (A) / Stage 2 (B: $(fmt(drug(PB.spec[1].mods).β)))", "$kname A2780 kill$nsnote"),
             prow("\\beta_C", hillC.β, "/day", "DATA · $(stage3_betaC().source) (A) / Stage 2 (B: $(fmt(drug(PB.spec[2].mods).β)))", "$kname A2780cis kill$nsnote"),
             (hill ? prow("h_N,\\ h_C", "$(fmt(hillN.h)), $(fmt(hillC.h))", "—", "DATA · Stage 2", "Hill slopes") : ""),
             prow("k_N,\\ k_C", "$(fmt(hillN.k)), $(fmt(hillC.k))", "/day", "DATA · Stage 2", "drug-effect transit rates"),
             (hill ? prow("EC_{50,N},\\ EC_{50,C}", "$(fmt(hillN.ec50)), $(fmt(hillC.ec50))", "µM", "DATA (fixed)", "measured IC50 of each line") : ""),
             lrrows(mN, "N", "A2780")..., lrrows(mC, "C", "A2780cis")...,
             (length(FP.ε) == 2 ? prow("\\varepsilon_N^{S4},\\ \\varepsilon_C^{S4}", "$(fmt(FP.ε[1])), $(fmt(FP.ε[2]))", "—", "DATA · Stage 4", "low-resource strength with 1 µM cisplatin") : ""),
             (chiN != 0 || chiC != 0 ? prow("\\chi", chiN == chiC ? fmt(chiN) : "$(fmt(chiN)), $(fmt(chiC))", "—", "DATA · Stage 4",
                                             "drug × low-resource interaction: the cisplatin kill is × (1 + χ) in 0.5% FBS") : ""),
             prow("u_R", U_R, "—", "DATA", "low-resource level at 0.5% FBS"),
             (TOX_MODEL == "xie2020" ?
              prow("k_N,\\ E_{max},\\ EC_{50}", "$(fmt(t.kN)), $(fmt(t.Emax)), $(fmt(t.EC50))", "/day, —, µM", "LITERATURE", "kidney injury response (Xie et al. 2020)") *
              prow("\\gamma", t.γ, "per µM", "ASSUMPTION (calibration)", "medium µM → surrogate plasma level: constant 1 µM gives the clinical mean excess Ē = $(fmt(XIE.Ebar))") :
              prow("s_B,\\ d_B", "$(fmt(t.s_B)), $(fmt(t.d_B))", "—, /day", "ASSUMPTION, LITERATURE", "drug burden (de Pillis & Radunskaya 2001)") *
              prow("a_H,\\ r_H", "$(fmt(t.a_H)), $(fmt(t.r_H))", "—, /day", "LITERATURE", "healthy-tissue kill and regrowth (de Pillis & Radunskaya 2001)")),
             prow("p_H", h.p_H, "—", "ASSUMPTION (swept 0 / 0.5 / 0.8)", "largest fraction of $(TOX_MODEL == \"xie2020\" ? \"kidney\" : \"healthy-tissue\") damage prevented by fasting (Raffaghello 2008; de Groot 2020)"),
             prow("k_P", h.k_P, "/day", "ASSUMPTION", "protection builds over ≈ 1 day of fasting"),
             prow("k_W", h.k_W, "/day", "LITERATURE-informed (Safdie 2009; swept)", "weight lost per fasting day (fraction of baseline)"),
             prow("r_W", h.r_W, "/day", "LITERATURE-informed (Safdie 2009; swept)", "weight regained after refeeding"),
             prow("W_{min}", W.Wmin, "—", "ASSUMPTION informed by Fearon 2011 (swept)", "weight-loss limit; Fearon's > 5% criterion defines cachexia over 6 months, not an acute fast (acute fasting loss is largely water and glycogen and is regained)"),
             (TOX_MODEL == "xie2020" ? prow("E_{lim}", W.Btox, "—", "ASSUMPTION (calibration)", "peak kidney excess of 75 mg/m² every 21 days (one dose: ΔE = $(fmt(XIE.ΔE)))") : prow("B_{tox}", W.Btox, "—", "LITERATURE", "healthy tissue ≥ 75% (de Pillis & Radunskaya 2001)")),
             prow("B_{cancer}", "$(BC_FRAC) K = $(cells(W.Bc))", "cells", "ASSUMPTION (swept)", "where the soft burden penalty (relaxed barrier) rises steeply; a penalty, not a hard cap")]
    """<h2>The model the controller uses $(tip("Every equation the controller simulates, in the order the states appear. The cell equations are generated from the plant itself, the same way the zoo pages write them, so what is shown is what is simulated.", :none))</h2>
<div class='intro'><p>Two controls, each held constant over a 2-day interval between media changes: the cisplatin concentration \\(u\\) (µM) and fasting \\(v\\) (1 = the interval is run in 0.5% FBS, 0 = fed, 10% FBS). Fasting enters the cells through the fitted low-resource level \\(u_R = $(fmt(U_R))\\,v\\), and the host through the weight and protection states. When \\(v = 0\\) the model is exactly the 14-day control model; with \\(u = 1\\) µM and \\(v = 1\\) it is exactly the fitted Stage 4 model (both checked below).</p>
<h3>Cells (scenario A; generated from the plant, with \\(D = u\\))</h3>
$(join([eqt(e, i == 1 ? "A2780 (sensitive). The cisplatin kill is the drug term (as selected in Stage 2b$(ns_on(S.lines[1]) ? ", scaled by the Norton–Simon growth factor" : "")), the low-resource module is the one Stage 2a selected for A2780 (its equations and values below): fasting acts through u_R." : "A2780cis (resistant). The cisplatin kill is the drug term (β_C from $(stage3_betaC().source) in scenario A). The low-resource module is the one Stage 2a selected for A2780cis (its equations and values below).") for (i, e) in enumerate(celltex[1:2])]))
$(length(celltex) > 2 ? join([eq(e) for e in celltex[3:end]]) : "")
$(eqt(raw"\frac{dA_{X,1}}{dt}=k_X\,(u-A_{X,1}),\quad \frac{dA_{X,j}}{dt}=k_X\,(A_{X,j-1}-A_{X,j}),\ j=2,3,\qquad E_X=A_{X,3}", "Drug-effect chain of each line (Stage 2): the kill follows the dose with a delay of about 3/k_X days."))
$(isempty(FP.iε) ? eqt("\\text{drug losses}\\times\\left(1+\\chi\\,u_R/0.955\\right),\\quad \\chi=$(fmt(FP.spec[1].chi))", "Drug × low-resource interaction (Stage 4): the carried model multiplies each line's cisplatin kill by (1 + χ) in 0.5% FBS and leaves it unchanged in 10% FBS, so it reproduces both measured conditions; no interpolation in the dose is needed.") : eqt(raw"\varepsilon_X(u)=1+\left(\varepsilon_X^{S4}-1\right)\min\!\left(\frac{u}{1\,\mu\mathrm{M}},\,1\right)", "Drug × low-resource interaction. Low resource alone was fitted with strength 1 (Stages 1b/3b); low resource with 1 µM cisplatin with ε^S4 (Stage 4). The plant reproduces both measured conditions and interpolates linearly in the dose between them (capped above 1 µM). The alternative, ε^S4 at every dose, predicts the low-resource cocultures much worse (checks below)."))
<h3>Host</h3>
$(TOX_MODEL == "xie2020" ? eqt(TX.eqf, "Kidney toxicity (Xie et al. 2020), with the fasting protection factor (1 − p_H·P). " * TX.what; warn=:tox_proxy) : host_tox_eqs_depillis())
$(eqt(raw"\frac{dP}{dt}=k_P\,(v-P)", "Protection state: builds towards 1 while fasting (about 1 day to establish, which is why the trials fast 24–48 h before the drug) and fades after refeeding. Its timing matters: fasting protects best when it starts before the dose."))
$(eqt(raw"\frac{dW}{dt}=-k_W\,v+r_W\,(1-W)(1-v)", "Body weight as a fraction of baseline: falls by k_W per fasting day and recovers towards 1 at rate r_W after refeeding (Safdie et al. 2009: a patient lost ≈ 3 kg over a multi-day fast, regained part in days and the rest within ≈ 3 weeks). This is the host cost of fasting."))
<h3>Objective (minimised over the look-ahead at every decision)</h3>
$(eqt(tox_tex(raw"\ell=q_P\,\sigma\!\Big(\tfrac{N+C-B_{prog}}{s}\Big)+q_B\,\tfrac{N+C}{K}+q_R\,\tfrac{C}{K}+q_T\Big(\tfrac{T}{B_{tox}}\Big)^2+q_W\Big(\tfrac{1-W}{1-W_{min}}\Big)^2+r_u\Big(\tfrac{u}{u_{max}}\Big)^2+r_v\,v^2"), "Stage cost after each interval. Goal since 2026-10-07: TIME TO PROGRESSION. σ is a smooth step (logistic, width s = $(CFG["progression"]["step_width_fraction_of_K"]) K) at B_prog = $(CFG["progression"]["threshold_fraction_of_K"]) K (80% confluence), so with q_P = $(CFG["progression"]["qP"]) the summed cost counts the intervals spent progressed: the controller keeps the culture under B_prog as long as it can, as adaptive therapy does (Zhang et al. 2017; Gatenby et al. 2009). The burden and resistance terms (q_B = $(CFG["weights"]["qB"]), q_R = $(CFG["weights"]["qR"])) are now tie-breakers. Earlier goal: a low total tumour burden with resistant cells held back. The burden (N + C)/K is linear, so every cell costs the same at any tumour size; the resistance premium q_R·C/K makes a resistant cell count as two (q_R = q_B), because further cisplatin cannot remove it. Then toxicity, weight loss, dose and fasting. Every term is scaled by its limit so they are comparable. Defaults: q_B = q_R = q_T = q_W = 1, r_u = 0.1, r_v = 0.05; the weights (including q_R = 0 and 3) are swept below."))
$(eqt(tox_tex(raw"+\,b_{\mu}\!\Big(1-\tfrac{N+C}{B_{cancer}}\Big)+b_{\mu}\!\Big(1-\tfrac{T}{B_{tox}}\Big)+b_{\mu}\!\Big(1-\tfrac{1-W}{1-W_{min}}\Big),\qquad b_\mu(z)=\begin{cases}-\mu\log z,&z\ge\delta\\ \mu\big[\tfrac12\big(\tfrac{z-2\delta}{\delta}\big)^2-\tfrac12-\log\delta\big],&z<\delta\end{cases}"), "Relaxed log barriers (Feller & Ebenbauer 2017) at three levels: the cancer burden B_cancer, healthy-tissue loss B_tox and weight loss W_min. A relaxed barrier is a soft constraint (a penalty): it grows steeply as the level is approached and continues as a steep quadratic beyond it, so passing a level has a finite cost. At B_cancer it is the soft burden penalty, which growing cultures pass; the tables count whether and for how many days. Its position and the terminal weight are swept below."))
$(eqt(tox_tex(raw"J=\sum_{k=1}^{H}\ell(x_{k+1},u_k,v_k)+5\Big[q_P\,\sigma(\cdot)+q_B\,\tfrac{N+C}{K}+q_R\,\tfrac{C}{K}+q_T\Big(\tfrac{T}{B_{tox}}\Big)^2+q_W\Big(\tfrac{1-W}{1-W_{min}}\Big)^2\Big]_{H+1}+q_{RT}\,\tfrac{C_{H+1}}{K}"), "Total cost over the look-ahead of $(sum(LOOKB)) intervals ($(LOOK_DAYS) days), planned as $(length(LOOKB)) blocks of $(join(LOOKB, ", ")) intervals with one dose each (move blocking, Cagienard et al. 2007; a fast applies to a block's first interval), plus q_RT·C/K (q_RT = $(CFG["progression"]["qRT"])) for the resistant cells left at its end, which drive the next progression; and with a terminal cost (×$(CFG["weights"]["terminal_factor"])) at the end. The first interval of the best plan is applied, then the plan is recomputed at the next media change (receding horizon). The look-ahead stays 14 days for every decision, also near day $DSTR (the plan runs past it), so the schedule is the first $NRUN decisions of a continuing policy with no end-of-run effect. MPC minimises this weighted cost, not the burden alone."))
$(eqt(raw"u_k\in\{0,\,0.62,\,1,\,1.47\}\ \mu\mathrm{M},\quad v_k\in\{0,1\},\quad \text{at most 1 fasting interval in a row (48 h)},\quad \ge 1\ \text{fed interval between fasts}", "Measured-levels mode: only the doses and media that were measured, and fasts no longer than 48 h (fasts of 24–72 h around platinum chemotherapy were feasible: Dorff et al. 2016) with at least 2 days of refeeding between them. The problem is then solved exactly by searching every admissible schedule. The demonstration mode relaxes this to u ∈ [0, 3] µM and v ∈ [0, 1] (partial serum reduction) and uses SQP with the same automatic-differentiation linearisation as the 14-day report."; warn="The 2 days of refeeding between fasts is an assumption, not a measured requirement."))
</div>
$(table([("Parameter", "Symbol in the equations.", :none), ("Value", "Value used by the controller (scenario A unless noted).", :none), ("Units", "Units.", :none),
         ("Source", "DATA = fitted to these experiments (stage named); LITERATURE = taken or derived from a cited study; ASSUMPTION = exploratory, swept below, to be replaced by measurements.", :none),
         ("Meaning", "What the parameter does.", :none)], prows))"""
end

# the controller's plant against the dedicated Stage 3b fit, on the same points, stated from the numbers
s3b_verdict(rp, r3) = !(isfinite(rp) && isfinite(r3)) ? "" : rp > 1.15r3 ? "The controller's model is therefore clearly less accurate than the dedicated low-resource fit for fasting without drug (it keeps the Stage 3a crowding coefficients)." :
    rp < r3 / 1.15 ? "The controller's model is therefore more accurate than the Stage 3b fit on these points." : "The two are about equally accurate (within 15%), so keeping the Stage 3a crowding coefficients costs little accuracy for fasting without drug."
function checks_section()
    ck = FC[:checks]
    r3(p) = hasproperty(p, :rmse_3b) ? cells(p.rmse_3b) : "—"
    rows = [trow([td(string(p.run)), td(p.ratio), td(p.line), tn(cells(p.rmse)), tn(r3(p)), tn(cells(p.last_obs)), tn(cells(p.last_pred))]) for p in ck[:interp].per]
    n3 = get(ck, :s3b_n, 0)
    """<h2>Checks against the measured low-resource conditions $(tip("The two conditions in which low resource was measured, and how well the controller's model reproduces them.", :none))</h2>
<div class='intro'><ul>
<li><b>Low resource + 1 µM cisplatin (Stage 4):</b> the plant with u = 1 µM, v = 1 reproduces the fitted Stage 4 model to $(@sprintf("%.1e", ck[:stage4_maxrel])) (relative).</li>
<li><b>Fed:</b> with v = 0 the plant is exactly the 14-day control model (tested to 10⁻¹²).</li>
<li><b>Low resource alone</b> (six low-resource cocultures, no drug, predicted with no refit from the nominal seeding): RMSE $(cells(ck[:interp].rmse)) cells with the dose-dependent strength ε(u) used here, versus $(cells(ck[:stage4].rmse)) cells if the Stage 4 strength were used at every dose, over all $(ck[:interp].n) observations after day 0 (QC-excluded points left out). $(n3 > 0 ? "On the $n3 of them that the Stage 3b fit also used, the controller's plant has an RMSE of $(cells(get(ck[:interp], :rmse_common, NaN))) cells and the Stage 3b fit $(cells(ck[:s3b_rmse])) cells. " * s3b_verdict(get(ck[:interp], :rmse_common, NaN), ck[:s3b_rmse]) : "The Stage 3b comparison is not in this checkpoint (re-run run_control_fasting.jl).") Every number is an unweighted root-mean-square error in cells; the Stage 3b fit itself minimised a weighted error (it down-weights the noisier A2780 counts), which is why its error in cells is not small. Using ε^S4 at every dose is far worse, so the drug × low-resource interaction matters.</li></ul></div>
$(table([("Run", "Low-resource coculture run (1 or 2).", :none), ("Ratio", "Seeding ratio sensitive:resistant.", :none), ("Line", "Cell line.", :none),
         ("RMSE (cells)", "Error of the controller's no-refit prediction of this trajectory under constant low resource, no drug.", :lower),
         ("Stage 3b RMSE (cells)", "Error of the Stage 3b fit on the same observations (in-sample, for reference).", :lower),
         ("Last observed", "Last measured count.", :none), ("Last predicted", "Predicted count on that day. Closer to the observed value is better.", "Closer to the observed value is better.")], rows))"""
end

const T_OUT = "Model simulation to day $DSTR. Resistant fraction = C/(N + C). Max toxicity is the highest $(TX.short); min weight is the lowest W. 'Above B_cancer' counts the soft burden penalty level: whether, and for how many days, N + C was at or above B_cancer. The toxicity and weight limits are 'hit' if the trajectory reaches them at any time."
function outcome_table(m)
    rows = String[]
    for pn in pol_order(m)
        v = m.policies[pn]; o = v.o
        flags = join(["<span class='flag'>breaks the $l limit</span>" for (l, b) in (("weight", o.weight_violation), ("toxicity", o.tox_violation)) if b], " ")
        push!(rows, trow([td("<span style='color:$(polcol(pn))'>■</span> " * hesc(pn) * (isempty(flags) ? "" : " " * flags)), tn(dayx(o.ttp)), tn(cells(o.total_final)),
                          tn(d0(get(o, :days_half_K, NaN))), tn(cells(o.C_final)), tn(cells(o.N_final)),
                          tn(@sprintf("%.2f", o.resistant_fraction)), tn(@sprintf("%.3f", o.max_tox)), tn(@sprintf("%.3f", o.min_weight)),
                          tn(@sprintf("%.1f", o.cum_dose)), tn(@sprintf("%.0f", o.fast_days)), td(yn(o.cancer_violation)), tn(d0(dover(o))), td(yn(o.tox_violation)), td(yn(o.weight_violation)),
                          tn(d0(o.days_progressed))];
                     cls=pn == MPC2 ? "win" : ""))
    end
    table([("Policy", "The dosing and fasting rule. MPC = model-predictive control. 'Dose only' may not fast; 'dose + fasting' chooses both from the measured levels. A policy that breaks the weight or toxicity limit is flagged in its row.", :none),
           ("Time to progression (day)", T_TTP, "Later is better; 'not progressed' is best: the goal of the controller."),
           ("N + C day $DSTR", "Total cells at day $DSTR.", :lower), ("Days above K/2", "Days with N + C at or above half the shared capacity K.", :lower),
           ("C day $DSTR", "Resistant A2780cis cells at day $DSTR.", :lower),
           ("N day $DSTR", "Sensitive A2780 cells at day $DSTR.", "Lower is better, but a surviving sensitive population also holds resistant cells back."),
           ("Resistant fraction", "C/(N + C) at day $DSTR.", :lower), ("Max toxicity", "Highest healthy-tissue loss T (limit 0.25).", :lower),
           ("Min weight", "Lowest body weight W (fraction of baseline; limit 0.95).", :higher),
           ("Dose (µM·d)", "Cumulative cisplatin ∫u dt.", :lower), ("Fasting days", "Days spent fasting (v = 1).", :lower),
           ("Above B_cancer", "N + C reached B_cancer (where the soft burden penalty rises steeply) at some time.", "'no' is better."),
           ("Days above B_cancer", "Days with N + C at or above B_cancer.", :lower), ("Toxicity limit hit", "$(TX.sym) reached $(TX.lim).", "'no' is better."),
           ("Weight limit hit", "W fell to W_min (5% loss).", "'no' is better.", W_WMIN),
           ("Days progressed", "Days with N + C at or above B_prog.", :lower)], rows) |> x -> withtip(x, T_OUT, :none; warn=[:extrapolation, :not_patient_days])
end

# Pareto view: every policy as one point, final cells against cumulative dose (filled: N + C; open: resistant C).
# The dashed line joins the policies no other policy beats on both total cells and dose.
function pareto_svg(title, pts; W=760, H=360)
    xm = max(maximum(q.x for q in pts) * 1.1, 1.0); ym = ceil(maximum(max(q.tot, q.C) for q in pts) * 1.12; sigdigits=2)
    L, R, T, B = 64, 16, 28, 112
    sx(v) = L + v / xm * (W - L - R); sy(v) = H - B - clamp(v, 0, ym) / ym * (H - T - B)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(hesc(title))'><rect width='100%' height='100%' fill='white'/><text x='$L' y='17' font-size='13' font-weight='bold'>$(hesc(title))</text>",
         "<text x='$(W-R)' y='17' font-size='10' text-anchor='end' fill='#52636b'>● N + C   ○ resistant C   - - not beaten on both</text>"]
    for v in range(0, ym; length=5)
        push!(p, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#e5e9ec'/><text x='%d' y='%.1f' font-size='10' text-anchor='end'>%.0f</text>", L, W - R, sy(v), sy(v), L - 6, sy(v) + 3, v))
    end
    for v in range(0, xm; length=6)
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%.0f</text>", sx(v), H - B + 14, v))
    end
    push!(p, "<line x1='$L' y1='$(H-B)' x2='$(W-R)' y2='$(H-B)' stroke='#555'/><text x='$((L+W-R)÷2)' y='$(H-B+28)' font-size='11' text-anchor='middle'>Cumulative dose (µM·d)</text><text transform='translate(13 $((T+H-B)÷2)) rotate(-90)' font-size='11' text-anchor='middle'>Cells at day $DSTR</text>")
    front = sort([q for q in pts if !any(r -> r.x <= q.x && r.tot <= q.tot && (r.x < q.x || r.tot < q.tot), pts)]; by=q -> q.x)
    length(front) > 1 && push!(p, "<polyline points='" * join([@sprintf("%.1f,%.1f", sx(q.x), sy(q.tot)) for q in front], " ") * "' fill='none' stroke='#94a3b8' stroke-dasharray='4,3'/>")
    for q in pts
        push!(p, @sprintf("<circle cx='%.1f' cy='%.1f' r='5.5' fill='%s'><title>%s: N + C = %.0f, dose %.1f µM·d</title></circle>", sx(q.x), sy(q.tot), q.color, hesc(q.label), q.tot, q.x))
        push!(p, @sprintf("<circle cx='%.1f' cy='%.1f' r='4.5' fill='white' stroke='%s' stroke-width='2'><title>%s: C = %.0f, dose %.1f µM·d</title></circle>", sx(q.x), sy(q.C), q.color, hesc(q.label), q.C, q.x))
    end
    for (i, q) in enumerate(pts)
        x = L + ((i - 1) % 3) * 230; y = H - B + 46 + ((i - 1) ÷ 3) * 13
        push!(p, "<circle cx='$(x+6)' cy='$(y-4)' r='4.5' fill='$(q.color)'/><text x='$(x+16)' y='$y' font-size='10'>$(hesc(first(q.label, 40)))</text>")
    end
    push!(p, "</svg>"); join(p)
end
const T_PARETO = "Each policy is one point at its cumulative dose (x) and its cells at day $DSTR (y): filled = total N + C, open = resistant C (50:50 start). " *
    "A point lower and further left is better: fewer cells for less drug. The dashed line joins the policies that no other policy beats on both total cells and dose. " *
    "Fasting days and weight are not on this plot (see the tables). MPC maximises the time to progression and prices the resistant cells it leaves, traded against toxicity, weight, dose, fasting and the soft penalties, so it need not have the fewest cells."

# the stage cost split into its terms over the run (control_fasting.jl: cost_terms2), same weights for every policy
function cost_table(scen, m)
    FP = fast_plant(scen); W = weights2(FP); rows = String[]
    for pn in pol_order(m)
        ct = cost_terms2(FP, m.x0, m.policies[pn].U, Δ0, W)
        push!(rows, trow([td("<span style='color:$(polcol(pn))'>■</span> " * hesc(pn)), tn(f3(ct.progression)), tn(f3(ct.burden)), tn(f3(ct.resistance)), tn(f3(ct.toxicity)), tn(f3(ct.dose)),
                          tn(f3(ct.penalty_burden)), tn(f3(ct.barrier_tox)), tn(f3(ct.weight + ct.fasting + ct.barrier_weight)), tn(f3(ct.total))]; cls=pn == MPC2 ? "win" : ""))
    end
    table([("Policy", "The dosing and fasting rule.", :none), ("Progression", "Σ q_P·σ((N + C − B_prog)/s): about q_P for every interval spent at or above B_prog = $(CFG["progression"]["threshold_fraction_of_K"]) K, the controller's main goal.", :lower), ("Burden", "Σ q_B·(N + C)/K over the $NRUN intervals.", :lower), ("Resistance premium", "Σ q_R·C/K.", :lower),
           ("Toxicity", "Σ q_T·($(TX.sym)/$(TX.lim))².", :lower), ("Dose", "Σ r_u·(u/u_max)², u_max = $(CFG["dose"]["u_max_data"]) µM for every policy.", :lower),
           ("Burden penalty", "Σ of the soft burden penalty b_μ(1 − (N + C)/B_cancer): log barrier below B_cancer, steep quadratic above it.", :lower),
           ("Toxicity barrier", "Σ b_μ(1 − $(TX.sym)/$(TX.lim)).", :lower), ("Fasting / weight", "Σ of q_W·((1 − W)/(1 − W_min))², r_v·v² and the weight barrier.", :lower),
           ("Total", "Σ ℓ over the run, the stage cost MPC trades off (its terminal cost is a planning device and is left out).", :lower)], rows) |>
        x -> withtip(x, "Each policy's $DSTR-day run scored with the controller's own stage cost, term by term (50:50 start). Large 'Burden penalty' values show where the soft penalty above B_cancer, not the cell count itself, drives the cost. MPC minimises the total over its $(LOOK_DAYS)-day look-ahead at every decision, not each term.", "Lower is better (in every column, and above all in Total).")
end

function scen_section(scen)
    m = main_of(scen, "50:50"); ts = m.policies["No treatment"].o.ts; W = m.W
    shown = [p for p in ("No treatment", "Constant 1 µM", "Periodic 1 µM with a 48-h fast on each dosing interval", MPC1, MPC2, MPCD) if haskey(m.policies, p)]
    tot = [(label=pn, color=polcol(pn), dash=poldash(pn), t=ts, y=m.policies[pn].o.N .+ m.policies[pn].o.C) for pn in shown]
    nc = vcat([[(label="$pn: N", color=polcol(pn), dash="", t=ts, y=m.policies[pn].o.N), (label="$pn: C", color=polcol(pn), dash="6,4", t=ts, y=m.policies[pn].o.C)] for pn in (MPC2, MPC1, "Constant 1 µM")]...)
    tox = [(label=pn, color=polcol(pn), dash=poldash(pn), t=ts, y=m.policies[pn].o.T) for pn in shown]
    wt = [(label=pn, color=polcol(pn), dash=poldash(pn), t=ts, y=m.policies[pn].o.W) for pn in shown]
    pr = [(label=pn, color=polcol(pn), dash=poldash(pn), t=ts, y=m.policies[pn].o.P) for pn in shown]
    ratios = join(["<h4>Seeding ratio $r</h4>" * outcome_table(main_of(scen, r)) for r in ("25:75", "50:50", "75:25")])
    pts = [(label=pn, color=polcol(pn), x=m.policies[pn].o.cum_dose, tot=m.policies[pn].o.total_final, C=m.policies[pn].o.C_final) for pn in pol_order(m)]
    """<section class='panel' id='scen-$scen'><h3>Scenario $(hesc(m.name)) <small>· $(scen_tag(scen))</small></h3>$(string(scen) == "B" ? "<p class='muted'>" * hesc(B_NOTE) * "</p>" : "")
<div class='grid'>
<div>$(withtip(ts_svg("Total N + C, 50:50 start", tot; hline=(W.Bc, "B_cancer"), ymax=ceil(max(W.Bc * 1.15, maximum(maximum(s.y) for s in tot) * 1.05); sigdigits=2)), "Total cells under each policy. The dashed orange line is B_cancer ($(BC_FRAC) of K), where the soft burden penalty rises steeply. The shaded region is after day 14.", :lower; warn=[:extrapolation, :soft_penalty]))</div>
<div>$(withtip(pareto_svg("Final cells vs cumulative dose, every policy (50:50)", pts), T_PARETO, "Lower and further left is better (fewer cells for less drug)."))</div>
<div>$(withtip(ts_svg("Sensitive (solid) and resistant (dashed) cells", nc), "Sensitive N and resistant C under MPC with fasting, MPC without fasting and constant 1 µM. Fasting kills the starvation-sensitive A2780cis cells in the fitted models, so compare the dashed curves.", :lower))</div>
<div>$(withtip(sched_svg("Schedules: dose (bars) and fasting (purple, F)", [(pn, m.policies[pn].U) for pn in pol_order(m)]), "Each row is a policy over $NRUN two-day intervals. Bar height is the cisplatin dose (hover for the value); purple shading marks fasting intervals (its strength for the continuous demonstration).", :none))</div>
<div>$(withtip(ts_svg(TX.name, tox; hline=(W.Btox, TX.lim), ymax=max(W.Btox * 1.3, maximum(maximum(s.y) for s in tox) * 1.1), ylabel=TX.ylabel), "Host toxicity with fasting protection. " * TX.limit_note, :lower; warn=:tox_proxy))</div>
<div>$(withtip(ts_svg("Body weight W", wt; hline=(W.Wmin, "W_min (assumed limit)"), ymin=min(0.9, minimum(minimum(s.y) for s in wt) - 0.005), ymax=1.005, ylabel="W (fraction of baseline)"), "Host cost of fasting: weight falls by k_W per fasting day and recovers after refeeding. The dashed line is the 5%-loss limit W_min.", :higher; warn=[W_WMIN, :not_patient_days]))</div>
<div>$(withtip(ts_svg("Fasting protection P", pr; ymax=1.0, ylabel="P (0 = none, 1 = full)"), "Protection of healthy tissue builds while fasting and fades after. A dose given while P is high does less damage to healthy tissue.", "Higher when the drug is given is better."))</div>
</div><h4>Cost of each run, term by term (50:50)</h4>$(cost_table(scen, m))
$ratios</section>"""
end

const ROB_POLS = ["No treatment", "Constant 1 µM", "Maximum measured dose (1.47 µM)", "Periodic 1 µM, 2 days on / 2 off", "Periodic 1 µM with a 48-h fast on each dosing interval",
                  "Constant 1 µM + 48-h fast every 6 days", "Constant 1 µM + 48-h fast every 8 days", MPC1, MPC2]
mrange(v) = any(isnan, v) || isempty(v) ? "—" : "$(cells(median(v))) ($(cells(quantile(v, .025)))–$(cells(quantile(v, .975))))"
function robust_section()
    secs = String[]; rows = String[]; prows = String[]
    q(v, p) = isempty(v) ? NaN : quantile(v, p)
    for scen in SCENS
        R = FC[:robust][scen]; ts = R.ts; m = main_of(scen, "50:50")
        bandof(M) = (t=ts, y=[median(M[i, :]) for i in axes(M, 1)], lo=[q(M[i, :], 0.025) for i in axes(M, 1)], hi=[q(M[i, :], 0.975) for i in axes(M, 1)])
        band(pn, f) = bandof(hcat([f(d[pn]) for d in R.draws if get(d, pn, nothing) !== nothing]...))
        ser = vcat([merge(band(pn, o -> o.N .+ o.C), (label="$pn: N + C", color=polcol(pn), dash="")) for pn in ("Constant 1 µM", MPC2)],
                   [merge(band(pn, o -> o.C), (label="$pn: C", color=polcol(pn), dash="6,4")) for pn in ("Constant 1 µM", MPC2)])
        isempty(R.closed) || push!(ser, merge(bandof(hcat([c.N .+ c.C for c in R.closed]...)), (label="Closed-loop MPC (noisy counts): N + C", color="#0f766e", dash="2,3")))
        push!(secs, "<div>" * withtip(ts_svg("Scenario $scen: N + C and C over bootstrap parameter sets", ser; hline=(m.W.Bc, "B_cancer")),
            "Each schedule (planned once with the best-fit model) replayed on $(length(R.draws)) bootstrap parameter sets: line = median, band = 2.5–97.5 percentiles. Wide bands = outcome uncertain under plausible parameters. The dotted band is the closed-loop MPC on its $(length(R.closed)) draws.", :lower) * "</div>")
        rowfor(name, ds; cls="") = (tot = [d.tot for d in ds]; cf = [get(d, :Cf, NaN) for d in ds];
            trow([td(string(scen)), td(name), tn(length(ds)), tn(pct(mean(d.cv for d in ds))), tn(pct(mean(d.tv for d in ds))), tn(pct(mean(d.wv for d in ds))),
                  tn(mrange(cf)), tn("$(cells(median(tot))) ($(cells(q(tot, .025)))–$(cells(q(tot, .975))))"),
                  tn(@sprintf("%.1f", median(d.cum for d in ds))), tn(@sprintf("%.0f", median(d.fast for d in ds)))]; cls=cls))
        for pn in ROB_POLS
            ds = [d[pn] for d in R.draws if get(d, pn, nothing) !== nothing]; isempty(ds) && continue
            push!(rows, rowfor(hesc(pn), ds; cls=pn == MPC2 ? "win" : ""))
        end
        isempty(R.closed) || push!(rows, rowfor("<b>Closed-loop MPC, dose + fasting</b> (re-plans every 2 days from noisy counts; first $(length(R.closed)) draws)", R.closed; cls="win"))
        # paired: closed loop minus each policy, draw by draw, on the draws that carry a closed-loop run
        P = [d for d in R.draws if get(d, CL_NAME, nothing) !== nothing]
        for pn in ROB_POLS
            pr = [(d[CL_NAME], d[pn]) for d in P if get(d, pn, nothing) !== nothing]; isempty(pr) && continue
            dt = [a.tot - b.tot for (a, b) in pr]; dc = [a.Cf - b.Cf for (a, b) in pr]; dd = [a.cum - b.cum for (a, b) in pr]
            push!(prows, trow([td(string(scen)), td(hesc(pn) * (pn == MPC2 ? " (open loop: the nominal schedule replayed)" : "")), tn(length(pr)),
                               tn(@sprintf("%+.0f", median(dt))), tn(pct(mean(dt .< 0))), tn(@sprintf("%+.0f", median(dc))), tn(pct(mean(dc .< 0))), tn(@sprintf("%+.1f", median(dd)))]))
        end
    end
    ncl = get(first(values(FC[:robust])), :ncl, 0)
    """<h2>Robustness over bootstrap parameter sets $(tip("How the 50:50 results change when the fitted parameters are replaced by bootstrap replicates (growth, drug, coculture and Stage 4 low-resource strengths). The host parameters are held at their defaults here; their effect is in the sweeps.", :none; warn=W_BORROW))</h2>
<div class='intro'><p>Each fixed schedule is replayed on the bootstrap parameter sets. Scenario B's crowding was fitted to the untreated cocultures, which have no bootstrap of their own: each draw borrows the Stage 3a bootstrap ratio (replicate ÷ estimate) of the same parameter. <b>Closed loop:</b> on the first $ncl draws MPC plans with the best-fit model, sees only the cell counts N and C at each media change with $(round(Int, 100CV_OBS))% multiplicative log-normal noise, and reconstructs the hidden states (drug-effect and stress chains, B, T, P, W) by simulating its own model under the doses and fasting it has applied. The comparison with each fixed policy is paired: the same plants (and noise) for every policy, on exactly those draws.</p></div>
<div class='grid'>$(join(secs))</div>
$(table([("Scen.", "Drug-response scenario.", :none), ("Policy", "Schedule planned with the best-fit model.", :none), ("Draws", "Parameter sets simulated.", :higher),
         ("Above B_cancer", "Share of parameter sets where N + C reaches B_cancer (the soft burden penalty level).", :lower), ("Toxicity limit hit", "Share where $(TX.sym) reaches $(TX.lim).", :lower),
         ("Weight limit hit", "Share where W reaches W_min.", :lower), ("C day $DSTR", "Median resistant cells (95% range).", :lower),
         ("N + C day $DSTR", "Median total (95% range).", :lower), ("Dose (µM·d)", "Median cumulative dose.", :lower), ("Fasting days", "Median days fasting.", :lower)], rows))
<h3>Paired comparison on the closed-loop draws $(tip("Closed-loop MPC minus each policy, draw by draw on the same bootstrap plants (paired), then summarised: the median difference and the share of draws in which the closed loop does better.", "Negative differences and shares above 50% favour the closed-loop MPC."))</h3>
$(table([("Scen.", "Scenario.", :none), ("Compared with", "The policy the closed-loop MPC is compared with on the same draws.", :none), ("Pairs", "Draws where both runs finished.", :higher),
         ("Δ N + C", "Median of (closed-loop − policy) total cells at day $DSTR.", "Negative favours the closed loop."), ("Closed loop fewer N + C", "Share of draws where the closed loop ends with fewer total cells.", :higher),
         ("Δ C", "Median of (closed-loop − policy) resistant cells at day $DSTR.", "Negative favours the closed loop."), ("Closed loop fewer C", "Share of draws where the closed loop ends with fewer resistant cells.", :higher),
         ("Δ dose (µM·d)", "Median of (closed-loop − policy) cumulative dose.", "Negative means the closed loop uses less drug.")], prows))"""
end

function washout_section()
    haskey(FC, :washout) || return ""
    # the fitted fade of A2780cis low-resource stress after refeeding: the chain's own rate, or none (cumulative exposure)
    mC0 = only(filter(M -> !M.input_D, fast_plant(:A).spec[2].mods))
    k0 = mC0.kernel_chain > 0 ? mC0.k : 0.0
    rows = String[]
    for scen in SCENS
        m = main_of(scen, "50:50")
        base = (koff="$(fmt(k0)) (fitted)", fo=m.policies["Fasting only, 48 h every 6 days"].o, mp=m.policies[MPC2].o, re=m.policies[MPC2].o, U=m.policies[MPC2].U)
        variants = vcat([base], [(koff=fmt(w.koff), fo=w.replay["Fasting only, 48 h every 6 days"], mp=w.replay[MPC2], re=w.o, U=w.U) for w in FC[:washout] if w.scen == scen])
        for v in variants
            push!(rows, trow([td(string(scen)), td(v.koff), tn(cells(v.fo.C_final)), tn(cells(v.mp.C_final)), tn(cells(v.mp.total_final)),
                              tn(cells(v.re.C_final)), tn(cells(v.re.total_final)), tn(@sprintf("%.0f", v.re.fast_days)), tn(@sprintf("%.1f", v.re.cum_dose))]))
        end
    end
    """<h2 id='washout'>How much rests on recovery after refeeding $(tip("How the results change if starvation stress fades after refeeding, which the fitted model does not allow for.", :none; warn=:refeeding))</h2>
<div class='intro'><p>$(k0 > 0 ? "In the fitted low-resource model (Stage 2a), starvation stress in A2780cis builds up through a chain of compartments at rate \\(k_C^R = $(fmt(k0))\\)/day, and after refeeding it fades at the same slow rate (about $(round(Int, 3 / k0)) days). One 48-h fast therefore keeps killing resistant cells for weeks." :
"In the fitted low-resource model (Stage 2a), starvation stress in A2780cis is a cumulative exposure: it builds up while the culture is in 0.5% FBS and never fades after refeeding, so every fast adds permanently to the resistant cells' death rate.") This is why the controller needs so little fasting. No experiment measured recovery after refeeding: every low-resource culture stayed in 0.5% FBS for all 14 days. Here the stress is allowed to fade once fed, at rate \\(k_{off}\\) ($(k0 > 0 ? "the chain rate becomes \\(k = k_C^R\\,v + k_{off}(1-v)\\)" : "\\(dS_C/dt = \\varepsilon_C u_R - k_{off}(1-v)\\,S_C\\)")), and the same schedules are replayed and re-optimised.</p>
<p><b>Experiment that would settle it:</b> A2780cis (and a 1:1 coculture) in 0.5% FBS for 48 h, then back to 10% FBS, counted daily for 10–14 days.</p></div>
$(table([("Scen.", "Scenario.", :none), ("k_off (/day)", "Fade rate of starvation stress after refeeding. The fitted model uses the build-up rate.", :none),
         ("C day $DSTR, fasting only", "Resistant cells after the fixed schedule 'Fasting only, 48 h every 6 days' (no drug).", :lower),
         ("C day $DSTR, MPC plan", "Resistant cells when the MPC schedule planned with the fitted model is run on this variant.", :lower),
         ("N + C, MPC plan", "Total cells for that schedule.", :lower),
         ("C day $DSTR, re-optimised", "Resistant cells when MPC plans with this variant.", :lower), ("N + C, re-optimised", "Total cells, re-optimised.", :lower),
         ("Fasting days", "Days fasting in the re-optimised schedule.", :lower), ("Dose (µM·d)", "Cumulative dose, re-optimised.", :lower)], rows))"""
end

function sweep_section()
    rows = String[]
    for s in FC[:sweeps]
        o = s.res.o; U = s.res.U
        sched = join([(U[k, 2] > 0 ? "<b>F</b>" : "") * @sprintf("%.3g", U[k, 1]) for k in axes(U, 1)], " ")
        push!(rows, trow([td(hesc(s.kind)), td(string(s.scen)), td(hesc(s.label)), td("<span style='font-size:11px'>$sched</span>"), tn(cells(o.C_final)), tn(cells(o.total_final)),
                          tn(@sprintf("%.2f", o.resistant_fraction)), tn(@sprintf("%.3f", o.max_tox)), tn(@sprintf("%.3f", o.min_weight)), tn(@sprintf("%.1f", o.cum_dose)), tn(@sprintf("%.0f", o.fast_days)),
                          tn(d0(dover(o))), tn(dayx(o.ttp))]))
    end
    """<h2>Sensitivity of the controller $(tip("MPC with dose + fasting from the 50:50 start, changing one assumption at a time. If the schedule and outcomes hold across the plausible range, the conclusion does not rest on that assumption.", :none))</h2>
<div class='intro'><p>The 'burden penalty' rows move the soft burden penalty (B_cancer at 0.9 or 1.0 K, or none) and the terminal factor (1 or 10) one at a time; they are scored against the default B_cancer = $(BC_FRAC) K, so 'Days above B_cancer' is comparable across rows. The full B_cancer × terminal-factor grid is in the 14-day report.</p></div>
$(table([("Setting", "Which assumption is varied.", :none), ("Scen.", "Scenario.", :none), ("Value", "The value used; everything else at its default.", :none),
         ("Schedule (µM; F = fasting)", "Dose in each 2-day interval; F marks a fasting interval.", :none), ("C day $DSTR", "Resistant cells at the end.", :lower),
         ("N + C day $DSTR", "Total at the end.", :lower), ("Resistant fraction", "C/(N + C) at the end.", :lower), ("Max toxicity", "Highest T.", :lower),
         ("Min weight", "Lowest W.", :higher), ("Dose (µM·d)", "Cumulative dose.", :lower), ("Fasting days", "Days fasting.", :lower),
         ("Days above B_cancer", "Days with N + C at or above the default B_cancer ($(BC_FRAC) K).", :lower),
         ("Time to progression (day)", T_TTP, "Later is better; 'not progressed' is best.")], rows))"""
end

function tests_section()
    f = joinpath(FOUT, "tests.txt")
    body = isfile(f) ? "<pre style='white-space:pre-wrap;background:#fbfdfd;border:1px solid #e4ecef;border-radius:8px;padding:10px'>$(hesc(read(f, String)))</pre>" : "<p class='muted'>Run test/test_control_fasting.jl.</p>"
    """<h2>Tests $(tip("From test/test_control_fasting.jl: the checks the two-input controller relies on.", "All tests passing is better."))</h2>
<div class='intro'><ul><li>Fed (v = 0) reproduces the 14-day control model to 10⁻¹² (relative); u = 1 µM with v = 1 reproduces the fitted Stage 4 model (&lt; 10⁻⁵ relative).</li>
<li>The fast, type-stable right-hand side used by the controller equals the zoo's generic right-hand side to 10⁻¹².</li>
<li>Fixed-step RK4 matches an adaptive Tsit5 solve (&lt; 10⁻⁵ relative); automatic-differentiation Jacobians A = ∂F/∂x and B = ∂F/∂(u, v) match central finite differences (&lt; 10⁻⁴ relative); the adjoint gradient matches finite differences of the cost (&lt; 10⁻⁵ of its largest component).</li>
<li>Every stage and terminal cost is ≥ 0 at random states (what makes pruning in the exact search safe), and the cost split by term adds up to the cost.</li>
<li>The exact search equals brute force with the fasting rules, and a warm start never changes its answer; SQP keeps both controls in bounds and never raises the cost.</li>
<li>The nominal closed loop with the fixed look-ahead equals replaying its own schedule.</li>
<li>Weight, protection and toxicity follow their equations; with p_H = 0 fasting gives no protection.</li></ul>$body</div>"""
end

function limits_section()
    """<h2 id='limits'>What the data do not support $(tip("Every caveat on this page in full; the ⚠ buttons next to the figures and tables show the one that applies there.", :none))</h2>
<div class='intro'><ul>
<li><b>No in-vivo attenuation yet.</b> The tumour-side effect of fasting is the in-vitro effect of 0.5% FBS at full strength. A patient's fast changes serum glucose, insulin and IGF-1 far less than a 20-fold serum cut, so the benefit shown here is an upper bound. Scaling it down is the planned next step.</li>
<li><b>Recovery after refeeding</b> of starvation damage was never measured, and the result depends on it strongly (section above).</li>
<li><b>Days 14–$DSTR</b> extrapolate every stage beyond the data. The model has no passaging: a culture near its capacity K stays there, which no real culture does for $DSTR days, so read the late days as the long-run tendency of the fitted equations, not as a prediction for a flask.</li>
<li><b>In-vitro days are not patient days.</b> Every rate here (growth, drug effect, starvation stress) is a culture rate in a dish; a 2-day fasting interval is 2 days in 0.5% FBS, and day $DSTR is a culture day. Mapping the schedule onto a patient's calendar would need pharmacokinetics and tumour growth rates in vivo, which none of these data provide.</li>
<li><b>The weight limit is an assumption.</b> W_min = $(HCFG["W_min"]) is informed by the > 5% weight-loss criterion for cancer cachexia (Fearon et al. 2011), but that criterion is defined over 6 months, not for an acute fast; weight lost in a short fast is largely water and glycogen and is regained after refeeding. It is swept (0.90–0.98).</li>
<li><b>The burden limit is a soft penalty.</b> B_cancer = $(BC_FRAC) K is where a relaxed barrier (a penalty) starts to dominate the cost, not a hard cap; growing cultures pass it, the tables count whether and for how long, and its position and the terminal weight are swept.</li>
<li><b>Fasting without drug</b> is predicted with an RMSE of $(cells(FC[:checks][:interp].rmse)) cells (the Stage 3b fit: $(cells(FC[:checks][:s3b_rmse])) cells on the points it used; checks above), and the drug × low-resource interaction was measured at one dose only (1 µM).</li>
<li><b>Host states are literature-informed proxies</b>: nothing in these experiments measured toxicity, weight or protection. The weight rates come from a case series, the weight limit is an assumption informed by a consensus cachexia definition, and the protection strength is an assumption (the randomized DIRECT trial found less DNA damage in T cells but no difference in toxicity grades, so p_H = 0 is in the sweep).</li>
<li><b>Resistant kill in coculture</b> (scenario A vs B) remains the largest structural uncertainty, as in the 14-day report. $(hesc(KILL_LIMIT_NOTE)) $(hesc(B_NOTE))</li><li><b>Beyond the data.</b> $(hesc(EXTRAP_NOTE)) $(hesc(LOWRES_NOTE)) Fasting is modelled as low resource, so every fasting result rests on that single treated low-resource series.</li><li><b>The toxicity limit is inactive at the default settings.</b> T can never exceed a_H/r_H = $(round(default_tox().a_H / default_tox().r_H; sigdigits=3)) (less while fasting protects healthy tissue), below the bound B_tox = $(CFG["toxicity"]["B_tox"]); the toxicity cost still discourages dosing, but the toxicity barrier never decides a schedule here.</li>
<li><b>Measurement and state estimation.</b> The closed-loop robustness runs assume each media change gives separate counts of N and C with $(round(Int, 100CV_OBS))% noise and that the hidden drug-effect, stress and host states follow the planner's own model. Scenario B's bootstrap draws borrow the Stage 3a bootstrap ratios for its crowding parameters.</li></ul></div>
<h2>References</h2>
<div class='intro'><ul class='muted'>
<li>Raffaghello L et al. (2008). Starvation-dependent differential stress resistance protects normal but not cancer cells against high-dose chemotherapy. PNAS 105:8215–8220.</li>
<li>Lee C et al. (2012). Fasting cycles retard growth of tumors and sensitize a range of cancer cell types to chemotherapy. Sci Transl Med 4:124ra27.</li>
<li>Safdie FM et al. (2009). Fasting and cancer treatment in humans: a case series report. Aging 1:988–1007.</li>
<li>Dorff TB et al. (2016). Safety and feasibility of fasting in combination with platinum-based chemotherapy. BMC Cancer 16:360.</li>
<li>de Groot S et al. (2020). Fasting mimicking diet as an adjunct to neoadjuvant chemotherapy for breast cancer in the multicentre randomized phase 2 DIRECT trial. Nat Commun 11:3083.</li>
<li>Obrist F et al. (2018). Metabolic vulnerability of cisplatin-resistant cancers. EMBO J 37:e98597.</li>
<li>Fearon K et al. (2011). Definition and classification of cancer cachexia: an international consensus. Lancet Oncol 12:489–495.</li>
$(join(["<li>$(hesc(r))</li>" for r in TX.refs]))
<li>de Pillis LG, Radunskaya AE (2001). A mathematical tumor model with immune resistance and drug therapy: an optimal control approach. J Theor Med 3:79–100 (values as tabulated in DTR-Bench, Luo et al. 2024, arXiv:2405.18610).</li>
<li>Feller C, Ebenbauer C (2017). Relaxed logarithmic barrier function based model predictive control of linear systems. IEEE Trans Autom Control 62:1223–1238.</li>
<li>Zhang J et al. (2017). Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. Nat Commun 8:1816.</li></ul></div>"""
end

# best fixed (non-MPC) policy on one outcome: (name, value)
function best_fixed(m, f; higher=false)
    fx = [pn for pn in pol_order(m) if !startswith(pn, "MPC")]
    vals = [f(m.policies[pn].o) for pn in fx]
    i = higher ? argmax(vals) : argmin(vals)
    (name=fx[i], v=vals[i])
end
function key_findings()
    out = String[]
    for scen in SCENS
        m = main_of(scen, "50:50"); P = m.policies
        mp = P[MPC2].o; m1 = P[MPC1].o; c1 = P["Constant 1 µM"].o
        U = P[MPC2].U; fdays = [@sprintf("%g–%g", (k - 1) * Δ0, k * Δ0) for k in axes(U, 1) if U[k, 2] > 0]
        fd = length(fdays) > 8 ? join(fdays[1:8], ", ") * ", … ($(length(fdays)) fasting intervals)" : (isempty(fdays) ? "none" : join(fdays, ", "))
        bt = best_fixed(m, x -> x.total_final); bc = best_fixed(m, x -> x.C_final); bd = best_fixed(m, x -> x.cum_dose)
        bx = best_fixed(m, x -> x.max_tox); bw = best_fixed(m, x -> x.min_weight; higher=true)
        dom = [pn for pn in pol_order(m) if !startswith(pn, "MPC") && (q = P[pn].o; q.total_final <= mp.total_final && q.C_final <= mp.C_final && q.cum_dose <= mp.cum_dose)]
        push!(out, "<li><b>Scenario $scen.</b> MPC with fasting ends at $(cells(mp.total_final)) cells ($(cells(mp.C_final)) resistant) vs $(cells(m1.total_final)) ($(cells(m1.C_final))) for MPC without fasting and $(cells(c1.total_final)) ($(cells(c1.C_final))) for constant 1 µM; " *
                   "it fasts on days $fd ($(@sprintf("%.0f", mp.fast_days)) days in total), gives $(@sprintf("%.1f", mp.cum_dose)) µM·d of cisplatin, " *
                   "and keeps weight ≥ $(@sprintf("%.3f", mp.min_weight)) and toxicity ≤ $(@sprintf("%.3f", mp.max_tox)). " *
                   "Against the best fixed policy on each metric at day $DSTR: total cells $(cells(mp.total_final)) vs $(cells(bt.v)) ($(hesc(bt.name))); resistant cells $(cells(mp.C_final)) vs $(cells(bc.v)) ($(hesc(bc.name))); " *
                   "cumulative dose $(@sprintf("%.1f", mp.cum_dose)) vs $(@sprintf("%.1f", bd.v)) µM·d ($(hesc(bd.name))); peak toxicity $(@sprintf("%.3f", mp.max_tox)) vs $(@sprintf("%.3f", bx.v)) ($(hesc(bx.name))); " *
                   "lowest weight $(@sprintf("%.3f", mp.min_weight)) vs $(@sprintf("%.3f", bw.v)) ($(hesc(bw.name))); " *
                   (bp = best_fixed(m, ttpkey; higher=true); "time to progression $(dayx(mp.ttp)) (without fasting $(dayx(m1.ttp))) vs $(isfinite(bp.v) ? @sprintf("%.1f", bp.v) : "not progressed") ($(hesc(bp.name))). ") *
                   (isempty(dom) ? "No fixed policy has at once no more total cells, no more resistant cells and no more drug than MPC." :
                                   "Fixed policies with no more total cells, resistant cells and drug than MPC: $(hesc(join(dom, ", "))).") * "</li>")
    end
    push!(out, "<li>MPC maximises the time to progression over its $(LOOK_DAYS)-day look-ahead and prices the resistant cells it leaves, traded against toxicity, weight, dose and fasting, so a fixed policy can beat it on a single end-of-run metric; the cost tables and the final-cells-vs-dose plots in each scenario show the trade-off.</li>")
    note = "<p class='flag'>$(warn(:refeeding)) These results rest on how slowly starvation damage in A2780cis fades after refeeding (see <a href='#washout'>How much rests on recovery after refeeding</a>)"
    if haskey(FC, :washout)
        w = only(filter(w -> w.scen == :A && w.koff == 2.0, FC[:washout]))
        note *= ": if it fades in about a day and a half (k_off = 2/day), the scenario-A MPC plan ends with $(cells(w.replay[MPC2].C_final)) resistant cells instead of $(cells(main_of(:A, "50:50").policies[MPC2].o.C_final)), and re-optimising gives $(cells(w.o.C_final)) with $(@sprintf("%.0f", w.o.fast_days)) fasting days"
    end
    "<ul>" * join(out) * "</ul>" * note * ".</p>"
end

const CSS = """body{font:14px/1.5 system-ui,sans-serif;margin:0;background:#f4f7f8;color:#17212b}main{max-width:1250px;margin:auto;padding:24px}h1{color:#124b6e;margin:6px 0}h2{color:#124b6e;margin-top:30px}h3{margin:14px 0 6px}h4{margin:14px 0 4px}
.topnav{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:0 0 16px}.nav{display:inline-flex;align-items:center;min-height:34px;padding:0 11px;border-radius:8px;color:#176b87;background:#eef5f6;text-decoration:none;font:650 14px system-ui,sans-serif;white-space:nowrap}.nav-home{padding:6px 9px}.nav-active{color:#fff;background:#176b87}
.intro,.panel{background:#fff;border:1px solid #d9e2e6;border-radius:12px;padding:16px 18px;margin:14px 0;min-width:0;overflow-wrap:anywhere}.intro{border-left:6px solid #176b87}.warnline{background:#fff6d9;border-radius:12px;padding:8px 14px;margin:14px 0}
.eq{overflow-x:auto;overflow-y:hidden;background:#fbfdfd;border:1px solid #e4ecef;border-radius:8px;padding:4px 12px;margin:8px 0;flex:1;min-width:0}.eqrow{display:flex;align-items:center;gap:4px}svg{display:block;width:100%;height:auto;margin:10px 0}.tscroll{overflow-x:auto}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:12px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}td,th{border:1px solid #d3dde2;padding:6px 8px;text-align:left;vertical-align:top}th{background:#e9f3f8}td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}tr.win td{background:#eef8f1;font-weight:600}.muted{color:#5b6b73}.flag{color:#b45309;font-weight:600}
code{overflow-wrap:anywhere}.tipwrap{position:relative}.tipwrap>.ii{position:absolute;top:8px;right:8px;z-index:2}mjx-container:not([display]){display:inline-block;max-width:100%;overflow-x:auto;vertical-align:middle}"""

# Drug-response shape: the main run (the carried Stage 2b drug model) next to a Hill alternative run with
# ZOO_STATE / CONTROL_OUT overrides into outputs/sensitivity/hill/control_<days>d (when it exists)
function hill_section()
    f = joinpath(@__DIR__, "outputs", "sensitivity", "hill", "control_$(DSTR)d", "fasting_outcomes.csv")
    fm = joinpath(FOUT, "fasting_outcomes.csv")
    head = "<h2>Drug-response shape (Hill vs linear) $(tip("The carried drug model is Stage 2b's selection ($DRUG_SHAPE kill). Here the day-$DSTR results of a run whose zoo state uses the Hill (saturating) kill instead are set next to the main ones, to show how much the conclusions rest on the shape of the dose response, which the data (three monoculture doses, one coculture dose) barely identify.", :none))</h2>"
    (isfile(f) && isfile(fm)) || return head * "<div class='intro'><p class='muted'>Not run: no Hill alternative at outputs/sensitivity/hill/control_$(DSTR)d/fasting_outcomes.csv (run run_control_fasting.jl with ZOO_STATE pointing to the Hill zoo state and CONTROL_OUT=outputs/sensitivity/hill).</p></div>"
    A = CSV.read(f, DataFrame); M = CSV.read(fm, DataFrame); rows = String[]
    pick(D, scen, pn) = (r = D[(string.(D.scenario) .== scen) .& (D.ratio .== "50:50") .& (D.policy .== pn), :]; nrow(r) == 1 ? r[1, :] : nothing)
    for scen in ("A", "B"), pn in (MPC2, MPC1, "Constant 1 µM")
        a = pick(M, scen, pn); b = pick(A, scen, pn)
        (a === nothing || b === nothing) && continue
        push!(rows, trow([td(scen), td(hesc(pn)), tn(cells(a.total_final)), tn(cells(b.total_final)), tn(cells(a.C_final)), tn(cells(b.C_final)),
                          tn(@sprintf("%.1f", a.cumulative_dose_uM_days)), tn(@sprintf("%.1f", b.cumulative_dose_uM_days)), tn(@sprintf("%.0f", a.fasting_days)), tn(@sprintf("%.0f", b.fasting_days))]; cls=pn == MPC2 ? "win" : ""))
    end
    head * "<div class='intro'><p>50:50 start, day $DSTR. 'Main' is the carried $DRUG_SHAPE drug model; 'Hill' is the same pipeline with the Hill (saturating) kill.</p>" *
        table([("Scen.", "Scenario.", :none), ("Policy", "MPC with or without fasting, or constant 1 µM.", :none),
               ("N + C, main", "Total cells at day $DSTR, main drug model.", :lower), ("N + C, Hill", "Total cells at day $DSTR, Hill alternative.", :lower),
               ("C, main", "Resistant cells, main.", :lower), ("C, Hill", "Resistant cells, Hill alternative.", :lower),
               ("Dose, main (µM·d)", "Cumulative dose, main.", "Lower is better when the cell counts are similar."), ("Dose, Hill (µM·d)", "Cumulative dose, Hill alternative.", "Lower is better when the cell counts are similar."),
               ("Fasting days, main", "Days fasting, main.", :lower), ("Fasting days, Hill", "Days fasting, Hill alternative.", :lower)], rows) * "</div>"
end

function write_report()
    home = "<svg width='22' height='22' viewBox='0 0 24 24' aria-hidden='true'><path fill='currentColor' d='M12 3.1 2.8 10.6v10.3h6.1v-6.2h6.2v6.2h6.1V10.6L12 3.1Zm7.1 15.8h-2v-6.2H6.9v6.2h-2V11.6L12 5.8l7.1 5.8v7.3Z'/></svg>"
    D = "../../../../"
    nl(l, h, a) = "<a href='$h' class='nav$(a ? " nav-active" : "")'>$l</a>"
    nav = "<nav class='topnav' aria-label='Report navigation'><a href='$(D)index.html' class='nav nav-home' aria-label='Reports home' title='Reports home'>$home</a>" *
          nl("Summary", "$(D)outputs/results/summary.html", false) * nl("Model zoo fits", "../zoo_fits/report.html", false) *
          nl("Optimal control (14 days)", "../control/report.html", false) * nl("Optimal control + fasting ($DSTR days)", "report.html", true) *
          nl("Final report", "$(D)outputs/final_report/report.html", false) * "</nav>"
    html = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Cisplatin and fasting control</title>
<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}};</script><script defer src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>
<style>$CSS</style>$INFO_TIPS_HEAD</head><body><main>$nav
<h1>Optimal control over $DSTR days: cisplatin and fasting</h1>
<div class='intro'><p>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")) by <code>run_control_fasting.jl</code> and <code>control_fasting_report.jl</code>. Model-predictive control now chooses two things at every media change (every 2 days) for $DSTR days: the cisplatin dose and whether the next 2 days are run in low resource ("fasting", 0.5% FBS). The goal is the same as in the 14-day report: the longest time to progression (keep N + C under B_prog = $(CFG["progression"]["threshold_fraction_of_K"]) K, 80% confluence), with the resistant cells left at the end of each look-ahead priced and a soft burden penalty above B_cancer; MPC minimises this weighted cost. Fasting has a host cost (body weight, with a 5% weight-loss limit) and a host benefit (protection of healthy tissue from the drug), both informed by the fasting-and-chemotherapy literature. Every run starts from the nominal seeding (100 cells split by the ratio). The look-ahead is $(LOOK_DAYS) days at every decision, also near day $DSTR (the plan runs past it), re-planned every 2 days, so the schedule is the first $NRUN decisions of a continuing policy. Scenarios A and B are as in the <a href='../control/report.html'>14-day report</a>.</p>
$(stopping_note())
<h3>Key findings (50:50 start) $(tip("Auto-generated from the runs below; each metric is compared with the fixed policy that is best on that metric.", :none))</h3>$(key_findings())</div>
<p class='warnline'><b>Read first</b> $(warn([:extrapolation, :not_patient_days, W_UPPER, W_HOST, W_UNTESTED])) <span class='muted'>click each ⚠; all caveats are listed in full under <a href='#limits'>What the data do not support</a>.</span></p>
$(model_section())
$(checks_section())
<div class='intro'><ul><li>$(hesc(EXTRAP_NOTE)) The tables give day $DSTR; days 15–$DSTR are model extrapolation.</li><li>$(hesc(LOWRES_NOTE))</li><li>$(hesc(KILL_LIMIT_NOTE))</li></ul></div>
<h2>Results $(tip("For each scenario: plots from the 50:50 start, then outcome tables for all three seeding ratios. Policies: five fed baselines (as in the 14-day report), four fixed fasting schedules in the style of the clinical studies (the 8-day one keeps the weight limit; the others are flagged where they break it), MPC without fasting, MPC with fasting (measured levels, solved exactly) and a continuous demonstration.", :none))</h2>
$(snapshot_section())
$(join([scen_section(sc) for sc in SCENS]))
$(washout_section())
$(haskey(FC, :robust) ? robust_section() : "<p class='muted'>Robustness runs not finished yet.</p>")
$(haskey(FC, :sweeps) ? sweep_section() : "<p class='muted'>Sweeps not finished yet.</p>")
$(hill_section())
$(tests_section())
$(limits_section())
<p class='muted'>Outputs: <code>fasting_outcomes.csv</code> in <code>Modeling_Approaches/02_mechanical_automatic_package/outputs/control_$(DSTR)d/</code>. Code: control_config.toml ([fasting], [host]), control_fasting.jl, run_control_fasting.jl, control_fasting_report.jl, test/test_control_fasting.jl.</p>
</main></body></html>"""
    write(joinpath(FOUT, "report.html"), stage_names(html)); if isempty(CONTROL_TAG)
        mkpath(FDOCS); cp(joinpath(FOUT, "report.html"), joinpath(FDOCS, "report.html"); force=true)
    end
    println("Wrote $(joinpath(isempty(CONTROL_TAG) ? FDOCS : FOUT, "report.html"))")
end

# One continuing run read off at the snapshot days: totals, resistant cells, dose and fasting up to each day
const SNAP_POLS = ["No treatment", "Constant 1 µM", MPC1, MPC2, MPCD]
function at_day(v, d)
    o = v.o; j = argmin(abs.(o.ts .- d)); k = round(Int, d / Δ0)
    (total=o.N[j] + o.C[j], C=o.C[j], cum=sum(v.U[1:k, 1]) * Δ0, fast=sum(v.U[1:k, 2]) * Δ0)
end
function snapshot_section()
    rows = String[]
    for scen in SCENS
        m = main_of(scen, "50:50")
        for pn in SNAP_POLS
            haskey(m.policies, pn) || continue
            v = m.policies[pn]
            cs = [td(string(scen)), td("<span style='color:$(polcol(pn))'>■</span> " * hesc(pn))]
            for d in SNAP
                a = at_day(v, d)
                push!(cs, tn("$(cells(a.total)) ($(cells(a.C)))"))
            end
            a = at_day(v, DAYS)
            push!(rows, trow(vcat(cs, [tn(@sprintf("%.1f", a.cum)), tn(@sprintf("%.0f", a.fast))]); cls=pn == MPC2 ? "win" : ""))
        end
    end
    head = vcat([("Scen.", "Drug-response scenario.", :none), ("Policy", "Dosing and fasting rule.", :none)],
                [("Day $(round(Int, d)): N + C (C)", "Total cells (resistant in brackets) at day $(round(Int, d)) of the same run.", :lower, d > DATA_DAYS ? :extrapolation : nothing) for d in SNAP],
                [("Dose to day $DSTR (µM·d)", "Cumulative cisplatin over the whole run.", :lower), ("Fasting days", "Days fasting over the whole run.", :lower)])
    """<h2>The run read off at days $(join(round.(Int, SNAP), ", ")) $(tip("One $DSTR-day closed-loop run, not separate runs. The controller always looks $(LOOK_DAYS) days ahead, also near the end (the plan runs past it), so a separate shorter run would make exactly the same decisions up to its end: reading the long run at several days is the same as running each length.", :none))</h2>
<div class='intro'><p>50:50 start. Each cell is the state of the same continuing run on that day, so the columns show how the outcome develops, not the best schedule for that horizon.</p>$(table(head, rows))</div>"""
end

function write_control_long_js()
    js(x::Real) = isfinite(x) ? string(round(Float64(x); sigdigits=6)) : "null"
    jo(o) = "{\"total\":$(js(o.total_final)),\"C\":$(js(o.C_final)),\"rf\":$(js(o.resistant_fraction)),\"cum\":$(js(o.cum_dose)),\"fast\":$(js(o.fast_days)),\"maxTox\":$(js(o.max_tox)),\"minW\":$(js(o.min_weight)),\"ttp\":$(js(get(o, :ttp, NaN)))}"
    rows = String[]
    for scen in (:A, :B)
        m = main_of(scen, "50:50"); P = m.policies; U = P[MPC2].U
        sj(v, d) = (a = at_day(v, d); "{\"total\":$(js(a.total)),\"C\":$(js(a.C))}")
        snap = join(["{\"day\":$(round(Int, d)),\"mpc2\":$(sj(P[MPC2], d)),\"mpc1\":$(sj(P[MPC1], d)),\"const1\":$(sj(P["Constant 1 µM"], d))}" for d in SNAP], ",")
        push!(rows, "{\"run\":\"Scenario $scen\",\"dose\":[$(join(js.(U[:, 1]), ","))],\"fast\":[$(join(js.(U[:, 2]), ","))],\"mpc2\":$(jo(P[MPC2].o)),\"mpc1\":$(jo(P[MPC1].o)),\"const1\":$(jo(P["Constant 1 µM"].o)),\"snap\":[$snap]}")
    end
    caveat = ""
    if haskey(FC, :washout)
        w = only(filter(w -> w.scen == :A && w.koff == 2.0, FC[:washout])); m1 = main_of(:A, "50:50").policies[MPC1].o
        caveat = "Caveat: the fitted model lets starvation damage in A2780cis fade slowly after refeeding (never measured). If it faded within about 1.5 days, the scenario-A plan would leave $(cells(w.replay[MPC2].C_final)) resistant cells instead of $(cells(main_of(:A, "50:50").policies[MPC2].o.C_final)); re-planning with that assumption still ends at $(cells(w.o.total_final)) total cells vs $(cells(m1.total_final)) without fasting."
    end
    out = "// written by control_fasting_report.jl — long-run dose + fasting control summary for the Summary and Presentation\nwindow.CONTROL_LONG = {\"generated\":\"$(Dates.format(now(), "yyyy-mm-dd HH:MM"))\",\"days\":$(round(Int, DAYS)),\"dataDays\":$(round(Int, DATA_DAYS)),\"Wmin\":$(HCFG["W_min"]),\"caveat\":\"$(replace(caveat, "\"" => "'"))\",\"rows\":[$(join(rows, ","))]};\n"
    write(joinpath(ROOT, "docs", "outputs", "results", "controlLong.js"), stage_names(out))
end
write_report()
isempty(CONTROL_TAG) && write_control_long_js()
