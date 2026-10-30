# HTML report of the control sweep (run_control_sweep.jl): one page with, per scenario and starting ratio,
# a schedule heat map (rows = variants, columns = 2-day intervals, colour = dose) and an outcome table of
# the data-supported MPC against the comparators under every variant.
#   julia --project=. control_sweep_report.jl
# Reads outputs/control_sweep/sweep_outcomes.csv; writes outputs/control_sweep/report.html and its docs copy.

using CSV, DataFrames, Dates, Printf, TOML

const PKG = @__DIR__
const ROOT = normpath(joinpath(PKG, "..", ".."))
const CFG = TOML.parsefile(joinpath(PKG, "control_config.toml"))
const SOUT = joinpath(PKG, "outputs", "control_sweep")
const SDOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "control_sweep")
const Δ0 = Float64(CFG["schedule"]["decision_interval_days"])
const DATA_DAYS = Float64(CFG["schedule"]["data_days"])
const MPC = "MPC, measured doses only"
const COMPARATORS = ["No treatment", "Constant 1 µM", "Maximum measured dose (1.47 µM)", "Periodic 1 µM, 2 days on / 2 off", "Threshold adaptive (Zhang-style)"]
const LEVELS = Float64.(CFG["dose"]["levels_data"])
const DOSE_COL = ["#f3f4f6", "#9ecae1", "#3182bd", "#08306b"]       # 0, 0.62, 1.0, 1.47 µM

hesc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "'" => "&#39;")
num(x; d=0) = x isa Number && isfinite(x) ? (d == 0 ? @sprintf("%.0f", x) : @sprintf("%.*f", d, x)) : "—"
ttp(x) = x isa Number && isfinite(x) ? @sprintf("day %.0f", x) : "not progressed"
ttpv(x) = x isa Number && isfinite(x) ? x : Inf                         # not progressed = best
doses(s) = ismissing(s) || isempty(strip(string(s))) ? Float64[] : parse.(Float64, split(string(s)))
dcol(u) = DOSE_COL[argmin(abs.(LEVELS .- u))]

function heatmap(sub)
    vs = unique(sub.variant); nmax = maximum(length(doses(r.schedule_uM)) for r in eachrow(sub))
    cw, rh, L = 14, 16, 470; W = L + nmax * cw + 10; H = 24 + rh * length(vs) + 26
    p = ["<svg viewBox='0 0 $W $H' width='$W' role='img' aria-label='MPC schedules by variant'><rect width='100%' height='100%' fill='white'/>"]
    for j in 0:5:nmax
        push!(p, "<text x='$(L + j * cw)' y='14' font-size='10' fill='#555'>d$(round(Int, j * Δ0))</text>")
    end
    dj = round(Int, DATA_DAYS / Δ0)
    push!(p, "<line x1='$(L + dj * cw)' x2='$(L + dj * cw)' y1='18' y2='$(H - 22)' stroke='#7c3aed' stroke-dasharray='3,3'/>")
    for (i, v) in enumerate(vs)
        r = sub[findfirst(==(v), sub.variant), :]; y = 20 + (i - 1) * rh
        r.beyond_data && push!(p, "<rect x='0' y='$y' width='$L' height='$rh' fill='#f5f3ff'/>")
        push!(p, "<text x='4' y='$(y + 12)' font-size='10.5'>$(hesc(r.axis)): $(hesc(v))</text>")
        for (j, u) in enumerate(doses(r.schedule_uM))
            push!(p, "<rect x='$(L + (j - 1) * cw)' y='$y' width='$(cw - 1)' height='$(rh - 2)' fill='$(dcol(u))'><title>$(hesc(v)), days $(round(Int, (j - 1) * Δ0))–$(round(Int, j * Δ0)): $(u) µM</title></rect>")
        end
    end
    ly = H - 10
    for (k, l) in enumerate(LEVELS)
        push!(p, "<rect x='$(L + (k - 1) * 70)' y='$(ly - 9)' width='12' height='10' fill='$(DOSE_COL[k])' stroke='#999' stroke-width='0.5'/><text x='$(L + (k - 1) * 70 + 16)' y='$ly' font-size='10'>$(l) µM</text>")
    end
    push!(p, "<text x='4' y='$ly' font-size='10' fill='#7c3aed'>shaded rows and right of the dashed line: beyond the data (day $(round(Int, DATA_DAYS)))</text></svg>")
    "<div class='tscroll'>" * join(p) * "</div>"
end

function outcome_table(sub)
    hdr = "<tr><th>Variant</th><th>MPC time to progression</th><th>MPC N + C end</th><th>MPC resistant end</th><th>MPC resistant fraction</th><th>MPC dose (µM·d)</th><th>Best comparator (time to progression)</th><th>No treatment N + C end</th></tr>"
    rows = String[]
    for v in unique(sub.variant)
        g = sub[sub.variant .== v, :]; m = g[g.policy .== MPC, :]
        isempty(m) && continue
        m = m[1, :]; c = g[in.(g.policy, Ref(COMPARATORS)), :]
        best = isempty(c) ? nothing : c[argmax([(ttpv(x.time_to_progression_day), -x.total_final) for x in eachrow(c)]), :]
        nt = g[g.policy .== "No treatment", :]
        beats = best === nothing ? "" : ttpv(m.time_to_progression_day) > ttpv(best.time_to_progression_day) ? " class='win'" : ""
        push!(rows, "<tr$(m.beyond_data ? " class='beyond'" : "")><td>$(hesc(m.axis)): $(hesc(v))</td><td$beats>$(ttp(m.time_to_progression_day))</td><td class='num'>$(num(m.total_final))</td><td class='num'>$(num(m.C_final))</td>" *
                     "<td class='num'>$(num(m.resistant_fraction; d=2))</td><td class='num'>$(num(m.cumulative_dose_uM_days; d=1))</td>" *
                     "<td>$(best === nothing ? "—" : hesc(best.policy) * " (" * ttp(best.time_to_progression_day) * ")")</td><td class='num'>$(isempty(nt) ? "—" : num(nt[1, :total_final]))</td></tr>")
    end
    "<div class='tscroll'><table>$hdr$(join(rows))</table></div>"
end

# per scenario: does the MPC schedule change across barriers (14 days) and across thresholds?
function verdict(df, data, scen)
    m = df[(df.data .== data) .& (df.scenario .== scen) .& (df.policy .== MPC), :]
    d14 = m[m.run_days .== Float64(CFG["schedule"]["closed_loop_days"]), :]
    def = m[m.axis .== "default", :]
    bar = d14[in.(d14.axis, Ref(("default", "barrier × run length"))), :]
    thr = m[in.(m.axis, Ref(("default", "threshold"))), :]
    ch(g) = isempty(g) ? "no runs" : (n = [length(unique(g[g.ratio .== r, :schedule_uM])) for r in unique(g.ratio)];
                                       all(==(1), n) ? "does not change" : "changes at $(count(>(1), n)) of $(length(n)) starting ratios")
    "<li>Scenario $(hesc(scen)): the 14-day MPC schedule <b>$(ch(bar))</b> across the $(length(unique(bar.barrier))) barrier shapes, and <b>$(ch(thr))</b> across thresholds $(join(sort(unique(thr.threshold_fraction_of_K)), ", ")) K.</li>"
end

f = joinpath(SOUT, "sweep_outcomes.csv")
df = isfile(f) ? CSV.read(f, DataFrame; stringtype=String, types=Dict(:schedule_uM => String, :ratio => String)) : DataFrame()
secs = String[]
if isempty(df)
    push!(secs, "<div class='intro'><p class='muted'>Not run yet: <code>julia -t 12 --project=. run_control_sweep.jl</code>.</p></div>")
else
    for data in unique(df.data)
        push!(secs, "<h2>Data: $(hesc(data))</h2><div class='intro'><ul>" * join([verdict(df, data, s) for s in unique(df.scenario)]) * "</ul></div>")
        for scen in unique(df.scenario), r in unique(df.ratio)
            sub = df[(df.data .== data) .& (df.scenario .== scen) .& (df.ratio .== r), :]
            isempty(sub) && continue
            push!(secs, "<details$(r == "50:50" ? " open" : "")><summary>Scenario $(hesc(scen)) · $(hesc(r)) start</summary><div class='panel'>" *
                        "<h3>MPC schedule per variant</h3>" * heatmap(sub[sub.policy .== MPC, :]) * "<h3>Outcomes</h3>" * outcome_table(sub) * "</div></details>")
        end
    end
end
nv = isempty(df) ? 0 : length(unique(df.variant))
html = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Control sweep</title>
<style>body{font-family:system-ui,-apple-system,Segoe UI,sans-serif;margin:0;background:#f4f7f8;color:#1f2a30;line-height:1.45}main{max-width:1200px;margin:0 auto;padding:16px}
.intro,.panel{background:#fff;border:1px solid #d9e2e6;border-radius:12px;padding:14px 18px;margin:12px 0;min-width:0}.intro{border-left:6px solid #176b87}
table{border-collapse:collapse;font-size:13px;width:100%}th,td{border-bottom:1px solid #e5eaec;padding:5px 7px;text-align:left;vertical-align:top}th{background:#eef3f5}
td.num{text-align:right;font-variant-numeric:tabular-nums}td.win{background:#e7f6ec;font-weight:600}tr.beyond td{background:#f5f3ff}.tscroll{overflow-x:auto}
details{margin:8px 0}summary{cursor:pointer;font-weight:600;padding:6px 0}.muted{color:#5b6b73}code{background:#eef3f5;padding:1px 4px;border-radius:4px}</style></head><body><main>
<h1>Control sweep: barriers, thresholds, run lengths</h1>
<div class='intro'><p>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")) by <code>run_control_sweep.jl</code> and <code>control_sweep_report.jl</code> (PRD, user 2026-10-27). Normal medium; $nv variants, each axis changed one at a time from the default ($(CFG["cancer"]["barrier"]) barrier, B_cancer = B_prog = $(CFG["cancer"]["B_cancer_fraction_of_K"]) K, $(round(Int, CFG["schedule"]["closed_loop_days"])) days, look-ahead $(CFG["schedule"]["lookahead_to_end"] ? "ending with the run" : "fixed"), time-to-progression objective), plus the full barrier × run-length cross. Every variant is scored with its own thresholds, and every comparator is rerun under it.</p>
<p><b>Barriers.</b> relaxed_log: the relaxed log barrier (Feller &amp; Ebenbauer 2017); quadratic: zero below B_cancer, a quadratic penalty above; hard: a near-wall that makes any interval at or above B_cancer cost about $(CFG["cancer"]["hard_weight"]), so a crossing plan is chosen only when every plan crosses; none: no burden penalty, the progression step only.</p>
<p><b>Read with care.</b> Runs longer than $(round(Int, DATA_DAYS)) days (shaded) follow the fitted equations beyond the data. A green cell marks an MPC that progresses later than every comparator. Scenario A (no resistant kill in coculture) cannot be rescued by any cost design: killing sensitive cells frees room for resistant ones.</p></div>
$(join(secs))
</main></body></html>"""
for d in (SOUT, SDOCS)
    mkpath(d); write(joinpath(d, "report.html"), html)
end
println("Wrote $(joinpath(SOUT, "report.html")) and its docs copy")
