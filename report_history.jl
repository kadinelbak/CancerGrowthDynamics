# Final Report history and modelling decisions (user, 2026-10-27: "I do not want to lose old reports anymore",
# "all decisions at the top of a final report", "compare reports and modelling decisions").
# Included by build_final_report.jl; standard library only (TOML, Dates, SHA), like that script.
#
#   docs/outputs/final_report/history/<id>/report.html   every Final Report version, never overwritten
#   docs/outputs/final_report/history/<id>/version.toml  its settings, key results, chapters and commit
#   docs/outputs/final_report/history/index.html         every version side by side, differences highlighted
#
# Settings are read from stages.toml and control_config.toml, key results from the committed output CSVs,
# and the dated decision log from decisions.toml. tools/seed_report_history.py made the versions that
# predate this file from git history, with the same keys and the same value formats as below.

using TOML, Dates, SHA

const HPKG = joinpath(@__DIR__, "Modeling_Approaches", "02_mechanical_automatic_package")
const HIST = joinpath(@__DIR__, "docs", "outputs", "final_report", "history")
include(joinpath(HPKG, "fit_set.jl"))      # which saved fit set the version uses (user, 2026-11-23)

hfmt(x::Bool) = x ? "yes" : "no"
hfmt(x::Integer) = string(x)
hfmt(x::AbstractFloat) = isinteger(x) ? string(Int(x)) : string(x)
hfmt(x::AbstractString) = String(x)
hfmt(x::AbstractVector) = join(hfmt.(x), ", ")
hfmt(::Nothing) = "not set"

# headline drug-response scenario (control_config.toml [scenarios].headline; user, 2026-11-02); "A" when the
# table is absent, so versions from before the decision still parse
headline_scenario(pkg=HPKG) = (f = joinpath(pkg, "control_config.toml"); isfile(f) ? String(get(get(TOML.parsefile(f), "scenarios", Dict()), "headline", "A")) : "A")

# the settings that define a version: (key, label, value), in display order
function report_settings(pkg=HPKG)
    sf, cf = joinpath(pkg, "stages.toml"), joinpath(pkg, "control_config.toml")
    S = isfile(sf) ? TOML.parsefile(sf) : Dict{String,Any}(); C = isfile(cf) ? TOML.parsefile(cf) : Dict{String,Any}()
    stg(id) = (i = findfirst(s -> get(s, "id", "") == id, get(S, "stage", Any[])); i === nothing ? Dict{String,Any}() : S["stage"][i])
    sel = get(S, "selection", Dict()); sch = get(C, "schedule", Dict()); can = get(C, "cancer", Dict())
    prog = get(C, "progression", Dict()); w = get(C, "weights", Dict())
    rule = get(sel, "rule_bic", nothing)
    hs = get(get(C, "scenarios", Dict()), "headline", nothing)
    wt = get(sel, "residual_weighting", "none")
    ex = get(stg("s2"), "exclude", Any[])
    srcf = joinpath(pkg, "data", "sources.toml")
    pts = isfile(srcf) ? [(get(s, "line", "?"), get(s, "exclude_days", Any[])) for s in get(TOML.parsefile(srcf), "source", Any[]) if !isempty(get(s, "exclude_days", Any[]))] : Any[]
    x0m = get(get(S, "starting_values", Dict()), "estimate", "auto")
    nst = get(sel, "extra_starts", 0)
    lines = get(S, "lines", Dict())
    Δ = Float64(get(sch, "decision_interval_days", 2.0))
    look = haskey(sch, "lookahead_blocks") ? "$(hfmt(sum(sch["lookahead_blocks"]) * Δ)) days (blocks $(hfmt(sch["lookahead_blocks"])))" :
           haskey(sch, "horizon_days") ? "$(hfmt(sch["horizon_days"])) days" : nothing
    days = sort(unique(Float64[x for x in vcat(get(sch, "closed_loop_days", Float64[]), get(sch, "run_days", Float64[]), get(sch, "long_days", Float64[]))]))
    bc = get(can, "B_cancer_fraction_of_K", nothing)
    obj = join(["$k = $(hfmt(v))" for (k, v) in (("qP", get(prog, "qP", 0.0)), ("qRT", get(prog, "qRT", 0.0)), ("qB", get(w, "qB", nothing)), ("qR", get(w, "qR", nothing)),
                                                  ("qT", get(w, "qT", nothing)), ("ru", get(w, "ru", nothing)), ("terminal factor", get(w, "terminal_factor", nothing))) if v !== nothing], ", ")
    [("rule_bic", "Model selection rule", rule === nothing ? "not recorded" : "fewest parameters within $(hfmt(rule)) BIC"),
     ("weighting", "Residual weighting", wt == "power" ? "power: max(y, 50)^$(hfmt(get(sel, "weight_power", 0.55)))" : hfmt(wt)),
     ("s2_exclude", "Stage 2b series left out", isempty(ex) ? "none" : join(["$(get(e, "line", "?")) $(hfmt(get(e, "dose", "?"))) µM" for e in ex], "; ")),
     ("points_excluded", "Single points left out (sources.toml exclude_days)", isempty(pts) ? "none" : join(["$(l) day $(hfmt(d))" for (l, d) in pts], "; ")),
     ("x0_estimate", "Starting values estimated", x0m == "mono" ? "every monoculture series" : x0m == "all" ? "every series" : "only when the first count is far from the seeding"),
     ("extra_starts", "Seeded extra starts per main fit", hfmt(nst)),
     ("kill_coupled", "Kill coupled to growth (Norton–Simon)", hfmt(get(stg("s2"), "kill_growth_coupled", false))),
     ("betaC_floor", "Resistant kill floor in coculture (× Stage 2b)", hfmt(get(stg("s3"), "betaC_floor", nothing))),
     ("ec50", "EC50, A2780 / A2780cis (µM)", isempty(lines) ? "not recorded" : "$(hfmt(get(get(lines, "A2780", Dict()), "ec50", nothing))) / $(hfmt(get(get(lines, "A2780cis", Dict()), "ec50", nothing)))"),
     ("headline_scenario", "Headline drug-response scenario", hs === nothing ? "not set" : hs == "B" ? "B (monoculture kills for both lines)" : hs == "A" ? "A (Stage 3a coculture fit as selected)" : hfmt(hs)),
     ("B_cancer", "B_cancer (burden penalty level)", bc === nothing ? "not set" : bc > 0 ? "$(hfmt(bc)) K" : "none"),
     ("B_prog", "B_prog (progression threshold)", haskey(prog, "threshold_fraction_of_K") ? "$(hfmt(prog["threshold_fraction_of_K"])) K" : "no progression goal"),
     ("barrier", "Burden barrier shape", hfmt(get(can, "barrier", "relaxed_log"))),
     ("lookahead", "Look-ahead", hfmt(look)),
     ("to_end", "Look-ahead ends with the run", hfmt(get(sch, "lookahead_to_end", false))),
     ("run_days", "Control run lengths (days)", isempty(days) ? "not set" : hfmt(days)),
     ("objective", "Objective weights", obj),
     ("B_tox", "Toxicity limit B_tox", hfmt(get(get(C, "toxicity", Dict()), "B_tox", nothing))),
     ("doses", "Measured doses (µM)", hfmt(get(get(C, "dose", Dict()), "levels_data", nothing)))]
end

# a minimal CSV reader (quoted fields with commas), so this file needs no packages
function hcsv(path)
    isfile(path) || return Dict{String,String}[]
    rows = Vector{String}[]
    for line in eachline(path)
        f = String[]; buf = IOBuffer(); q = false; i = 1; cs = collect(line)
        while i <= length(cs)
            c = cs[i]
            if q
                c == '"' ? (i < length(cs) && cs[i+1] == '"' ? (print(buf, '"'); i += 1) : (q = false)) : print(buf, c)
            else
                c == '"' ? (q = true) : c == ',' ? push!(f, String(take!(buf))) : print(buf, c)
            end
            i += 1
        end
        push!(f, String(take!(buf))); push!(rows, f)
    end
    isempty(rows) && return Dict{String,String}[]
    [Dict(zip(rows[1], r)) for r in rows[2:end] if length(r) == length(rows[1])]
end

# key results of a version: the carried model of every fitted group and the 50:50 MPC of each scenario
function report_results(pkg=HPKG)
    out = Tuple{String,String,String}[]
    for r in hcsv(joinpath(pkg, "outputs", "zoo_fits", "zoo_model_selection.csv"))
        get(r, "optimal", "") == "true" && push!(out, ("carried:" * r["group"], "Carried model · " * r["group"], r["candidate"]))
    end
    cc = hcsv(joinpath(pkg, "outputs", "control", "control_outcomes.csv"))
    hl = headline_scenario(pkg)
    for scen in (hl == "B" ? ("B", "A") : ("A", "B"))
        tag = scen == hl ? " (main)" : ""
        i = findfirst(r -> get(r, "environment", "") == "normal" && get(r, "scenario", "") == scen && get(r, "ratio", "") == "50:50" &&
                           get(r, "policy", "") == "MPC, measured doses only", cc)
        i === nothing && continue
        r = cc[i]
        push!(out, ("mpc:$scen", "14-day MPC schedule · scenario $scen$tag, 50:50 (µM)", r["schedule_uM"]))
        push!(out, ("mpc_total:$scen", "14-day MPC N + C at the end · scenario $scen$tag, 50:50", string(round(Int, parse(Float64, r["total_final"])))))
        haskey(r, "time_to_progression_day") && push!(out, ("mpc_ttp:$scen", "14-day MPC time to progression · scenario $scen$tag, 50:50",
                                                            (v = tryparse(Float64, r["time_to_progression_day"]); v === nothing || !isfinite(v) ? "not progressed" : "day $(hfmt(round(v; digits=1)))")))
    end
    append!(out, research_quality(pkg))
    out
end

# research-quality metrics (user, 2026-11-01: "what metrics across the project are the strongest to maximize ...
# can I measure them between reports"). Keys start with "q+:" (higher is better) or "q-:" (lower is better), and
# the value starts with the number that is compared, so the history page can mark a change better or worse.
# tools/seed_report_history.py --quality computes the same rows from git history for older versions.
hnum(s) = (v = tryparse(Float64, string(s)); v === nothing ? NaN : v)
hpct(a, b) = b == 0 ? "n/a" : "$(round(Int, 100a / b))% ($a of $b)"
function hmedian(x)
    isempty(x) && return NaN
    y = sort(x); n = length(y); isodd(n) ? y[(n+1)÷2] : (y[n÷2] + y[n÷2+1]) / 2
end
# MPC beats a baseline: it progresses later, or at the same time (or never) with fewer cells at the end
function mpc_beats(m, b, run_days)
    ttp(r) = (v = hnum(get(r, "time_to_progression_day", "")); isfinite(v) ? v : run_days + 1)
    tm, tb = ttp(m), ttp(b)
    tm > tb + 1e-6 || (abs(tm - tb) <= 1e-6 && hnum(m["total_final"]) < hnum(b["total_final"]))
end
function mpc_wins(rows, groupkeys, run_days_of)
    groups = Dict{Vector{String},Vector{Dict{String,String}}}()
    for r in rows
        push!(get!(groups, [get(r, k, "") for k in groupkeys], Dict{String,String}[]), r)
    end
    won = 0; n = 0
    for (_, g) in groups
        pol(name) = (i = findfirst(r -> get(r, "policy", "") == name, g); i === nothing ? nothing : g[i])
        m = pol("MPC, measured doses only"); a = pol("No treatment"); c = pol("Constant 1 µM")
        (m === nothing || a === nothing || c === nothing) && continue
        n += 1; rd = run_days_of(m)
        won += mpc_beats(m, a, rd) && mpc_beats(m, c, rd)
    end
    won, n
end
function research_quality(pkg=HPKG)
    out = Tuple{String,String,String}[]
    o = joinpath(pkg, "outputs")
    v = [hnum(r["nrmse_heldout"]) for r in hcsv(joinpath(o, "zoo_analyses", "validation_blocked.csv")) if get(r, "status", "") == "completed"]
    v = filter(isfinite, v)
    isempty(v) || push!(out, ("q-:heldout", "Quality · held-out error, median NRMSE over blocked checks (lower is better)", "$(round(hmedian(v); digits=3)) (n = $(length(v)), worst $(round(maximum(v); digits=3)))"))
    ce = hcsv(joinpath(o, "zoo_fits", "heldout_ce1.csv"))
    i = findfirst(r -> startswith(get(r, "row", ""), "prediction"), ce)
    i === nothing || push!(out, ("q-:heldout_ce1", "Quality · held-out coculture run (Ce1) NRMSE (lower is better)", string(round(hnum(ce[i]["nrmse"]); digits=3))))
    pr = hcsv(joinpath(o, "zoo_analyses", "sensitivity_profile.csv"))
    if !isempty(pr)
        ok = count(r -> get(r, "status", "") == "identifiable", pr)
        push!(out, ("q+:pinned", "Quality · parameters pinned by their profile (higher is better)", hpct(ok, length(pr)) *
              (ok < length(pr) ? "; not pinned: " * join(["$(r["param"]) ($(r["group"]))" for r in pr if get(r, "status", "") != "identifiable"], ", ") : "")))
    end
    bs = [r for r in hcsv(joinpath(o, "zoo_analyses", "bootstrap_selection.csv")) if get(r, "optimal", "") == "true"]
    if !isempty(bs)
        p = [hnum(r["pct_chosen"]) for r in bs]; j = argmin(p)
        push!(out, ("q+:stability", "Quality · selection stability, mean % of bootstrap refits that pick the carried model (higher is better)",
                    "$(round(Int, sum(p) / length(p)))% (lowest $(round(Int, p[j]))%, $(bs[j]["group"]))"))
    end
    cc = [r for r in hcsv(joinpath(o, "control", "control_outcomes.csv")) if get(r, "environment", "") == "normal"]
    if !isempty(cc)
        w, n = mpc_wins(cc, ["scenario", "ratio"], _ -> 14.0)
        n > 0 && push!(out, ("q+:robust", "Quality · 14-day MPC beats no treatment and constant 1 µM, scenarios × ratios (higher is better)", hpct(w, n)))
        hl = headline_scenario(pkg)
        wm, nm = mpc_wins([r for r in cc if get(r, "scenario", "") == hl], ["scenario", "ratio"], _ -> 14.0)
        nm > 0 && push!(out, ("q+:robust_main", "Quality · 14-day MPC beats no treatment and constant 1 µM, headline scenario $hl × ratios (higher is better)", hpct(wm, nm)))
    end
    sw = [r for r in hcsv(joinpath(o, "control_sweep", "sweep_outcomes.csv")) if get(r, "environment", "normal") == "normal"]
    if !isempty(sw)
        w, n = mpc_wins(sw, ["data", "variant", "scenario", "ratio"], r -> (d = hnum(get(r, "run_days", "14")); isfinite(d) ? d : 14.0))
        n > 0 && push!(out, ("q+:robust_sweep", "Quality · MPC beats both baselines across the control sweep (higher is better)", hpct(w, n)))
    end
    if !isempty(cc)
        same = 0; tot = 0
        for ratio in unique(get(r, "ratio", "") for r in cc)
            f(sc) = (i = findfirst(r -> r["scenario"] == sc && r["ratio"] == ratio && r["policy"] == "MPC, measured doses only", cc); i === nothing ? nothing : split(cc[i]["schedule_uM"]))
            a, b = f("A"), f("B")
            (a === nothing || b === nothing || length(a) != length(b)) && continue
            same += count(k -> hnum(a[k]) == hnum(b[k]), eachindex(a)); tot += length(a)
        end
        tot > 0 && push!(out, ("q+:agree", "Quality · scenario agreement, 14-day MPC dose blocks the same in A and B (higher is better)", hpct(same, tot)))
    end
    out
end

# data provenance (user, 2026-11-02, RP1: "If the data was reused from previous report please have that indicated.
# A new report assumes new settings or a new optimal control/qualification."). Read from outputs/run_manifest.toml
# (run_pipeline.jl) when a version is archived and stored in its version.toml as [provenance]:
#   recorded, started, finished, commit, continued_from, resume, summary, partly_reused, and one [[provenance.steps]]
#   per tracked step with status = "new" (run in this version) | "reused_unchanged" (control inputs unchanged; the step
#   carries reused_from, written by run_control*.jl) | "reused_earlier" (run_pipeline.jl --from/--only: the step finished
#   before this run started) | "resumed" (fits with --resume: unchanged fits/bootstraps reused) | "not_recorded".
const PROV_STEPS = [("fits", "Fits"), ("analyses", "Analyses"), ("control14", "Control"), ("control_long", "Long control"),
                    ("control30", "30-day control"), ("control60", "60-day control"), ("sweep", "Sweep"), ("sweep_all", "Sweep (all doses)")]
function hwhen(s)
    t = string(s)
    try Dates.format(DateTime(t[1:min(19, end)]), "d u HH:MM") catch; t end
end
function report_provenance(pkg=HPKG)
    none = Dict{String,Any}("recorded" => false, "partly_reused" => false, "summary" => "provenance not recorded", "steps" => Dict{String,Any}[])
    mf = joinpath(pkg, "outputs", "run_manifest.toml")
    isfile(mf) || return none
    M = try TOML.parsefile(mf) catch; nothing end
    M === nothing && return none
    steps = get(M, "steps", Dict()); started = string(get(M, "started", "")); cont = string(get(M, "continued_from", ""))
    resume = get(M, "resume", false) === true
    recs = Dict{String,Any}[]; parts = String[]; reused = false
    for (name, label) in PROV_STEPS
        st = get(steps, name, nothing)
        status = "not_recorded"; when = ""; text = "not recorded"
        if st !== nothing
            fin = string(get(st, "finished", ""))
            if haskey(st, "reused_from")
                status = "reused_unchanged"; when = string(st["reused_from"])
                text = "reused from the run of $(hwhen(when)) (inputs unchanged)"
            elseif !isempty(started) && !isempty(fin) && fin < started
                status = "reused_earlier"; when = fin
                text = "reused from an earlier run of $(hwhen(when)) (--from / --only)"
            elseif name == "fits" && resume
                status = "resumed"; when = fin
                text = "resumed (--resume: unchanged fits and bootstraps reused from an earlier run)"
            else
                status = "new"; when = fin; text = "new in this version"
                # saved fits and bootstraps reused unchanged inside the fitting step (since 2026-11-08 / 2026-11-09)
                if name == "fits" && haskey(st, "fits_total")
                    nb = length(get(st, "bootstraps_reused", []))
                    text *= " ($(get(st, "fits_reused", 0)) of $(st["fits_total"]) candidate fits and $nb of $(get(st, "bootstraps_total", 0)) " *
                            "bootstraps reused unchanged)"
                    (get(st, "fits_reused", 0) > 0 || nb > 0) && (reused = true)
                end
            end
        elseif name == "control60"
            status = "not_run"; text = "not run (only with --control60)"
        end
        status in ("reused_unchanged", "reused_earlier", "resumed") && (reused = true)
        push!(recs, Dict{String,Any}("step" => name, "label" => label, "status" => status, "when" => when, "detail" => text))
        push!(parts, "$label: $text")
    end
    # the saved fit set the version's fits are (fit_set.jl): its name heads the summary, so the version list shows it
    fs = fit_set_info(joinpath(pkg, "outputs", "zoo_fits", "zoo_state.jls"))
    fitset = Dict{String,Any}("name" => string(get(fs, "name", "")), "matched" => get(fs, "matched", false), "fit_run" => string(get(fs, "fit_run", "")),
                              "fit_commit" => string(get(fs, "fit_commit", "")), "text" => fits_used_text(fs))
    pushfirst!(parts, "Fits used: $(fitset["name"])")
    Dict{String,Any}("recorded" => true, "started" => started, "finished" => string(get(M, "finished", "")),
                     "commit" => string(get(get(M, "git", Dict()), "commit", "")), "continued_from" => cont, "resume" => resume,
                     "partly_reused" => reused, "summary" => join(parts, " · "), "steps" => recs, "fit_set" => fitset)
end
# the provenance line shown at the top of the Final Report
function provenance_html(prov)
    (prov === nothing || get(prov, "recorded", false) !== true) &&
        return "<p class='prov'><b>Data provenance:</b> <span class='mut'>provenance not recorded</span> (no run manifest when this version was built).</p>"
    tag = get(prov, "partly_reused", false) ? " <b class='tag'>partly reused</b>" : " <b class='tag'>all steps new</b>"
    "<p class='prov$(get(prov, "partly_reused", false) ? " warn" : "")'><b>Data provenance:</b>$tag $(hesc_(get(prov, "summary", ""))). " *
    "Run started $(hesc_(hwhen(get(prov, "started", ""))))$(isempty(get(prov, "finished", "")) ? "" : ", finished " * hesc_(hwhen(prov["finished"])))" *
    "$(isempty(get(prov, "commit", "")) ? "" : ", commit <code>" * hesc_(first(prov["commit"], 7)) * "</code>")." *
    "$(get(prov, "partly_reused", false) ? " A new report assumes new settings or a new optimal control; read the reused steps as results of the earlier run." : "")</p>" *
    (haskey(prov, "fit_set") ? "<p class='prov'><b>$(hesc_(prov["fit_set"]["text"]))</b></p>" : "")
end

decision_log(pkg=HPKG) = (f = joinpath(pkg, "decisions.toml"); isfile(f) ? get(TOML.parsefile(f), "decision", Any[]) : Any[])

# every archived version, oldest first
function load_versions()
    isdir(HIST) || return Dict{String,Any}[]
    vs = [TOML.parsefile(joinpath(HIST, d, "version.toml")) for d in readdir(HIST) if isfile(joinpath(HIST, d, "version.toml"))]
    sort(vs; by=v -> v["id"])
end
latest_version() = (v = load_versions(); isempty(v) ? nothing : v[end])
lookup(rows, key) = (i = findfirst(r -> r["key"] == key, rows); i === nothing ? nothing : rows[i]["value"])

# an archived copy lives two folders deeper than docs/outputs/final_report/report.html: every link from
# the page (chapter frames, their sources, the home link) gets two more "../"
function rebase_archived(page)
    head, sep, rest = let i = findfirst("const PARTS", page); i === nothing ? (page, "", "") : (page[1:prevind(page, first(i))], "const PARTS", page[first(i)+length("const PARTS"):end]) end
    head = replace(head, "href='../../index.html'" => "href='../../../../index.html'", "href='../archive" => "href='../../../archive",
                   "href='history/" => "href='../", "action='history/" => "action='../", "href='report.html'" => "href='../../report.html'")
    rest = replace(rest, "<base href=\\\"../../" => "<base href=\\\"../../../../", "src:\"../../" => "src:\"../../../../")
    head * sep * rest
end

function hjson(x)
    x isa AbstractDict && return "{" * join([hjson(string(k)) * ":" * hjson(v) for (k, v) in x], ",") * "}"
    x isa AbstractVector && return "[" * join(hjson.(x), ",") * "]"
    x isa Bool && return x ? "true" : "false"
    x isa Number && return isfinite(x) ? string(x) : "null"
    x === nothing && return "null"
    io = IOBuffer(); print(io, '"')
    for c in string(x)
        c == '"' ? print(io, "\\\"") : c == '\\' ? print(io, "\\\\") : c == '\n' ? print(io, "\\n") : c < ' ' ? print(io, "\\u", string(UInt16(c); base=16, pad=4)) : print(io, c)
    end
    print(io, '"'); replace(String(take!(io)), "</" => "<\\/")
end

# the compare page: every version as a column, every setting and key result as a row
function write_history_index()
    vs = load_versions()
    tpl = read(joinpath(@__DIR__, "report_history_template.html"), String)
    data = Dict("versions" => vs, "decisions" => decision_log(), "generated" => Dates.format(now(), "yyyy-mm-dd HH:MM"))
    mkpath(HIST); write(joinpath(HIST, "index.html"), replace(tpl, "/*HISTORY_JSON*/null" => hjson(data)))
end

hesc_(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "'" => "&#39;")

# the decisions panel at the top of the Final Report: this version's settings and results against the
# previous archived version, and the dated decision log
function decisions_panel(id, settings, results, prev; chapters=Tuple{String,String}[], provenance=nothing)
    previd = prev === nothing ? nothing : prev["id"]
    # chapters generated before the newest decision show results from older settings
    lastdec = maximum((get(d, "date", "") for d in decision_log()); init="")
    stale = [(t, st) for (t, st) in chapters if !isempty(lastdec) && st[1:min(10, end)] < lastdec]
    stalehtml = isempty(stale) ? "" : "<p class='warn'><b>$(length(stale)) of $(length(chapters)) chapters predate the newest decision ($(hesc_(lastdec)))</b>, so they still show results from earlier settings: " *
        join(["$(hesc_(t)) ($(hesc_(st)))" for (t, st) in stale], "; ") * ". Rerun the pipeline and rebuild this report to bring them up to date.</p>"
    function row(k, l, v, rows)
        pv = prev === nothing ? nothing : lookup(get(prev, rows, Any[]), k)
        ch = pv !== nothing && pv != v
        "<tr$(ch ? " class='chg'" : "")><td>$(hesc_(l))</td><td>$(hesc_(v))</td><td>$(pv === nothing ? "<span class='mut'>—</span>" : hesc_(pv))$(ch ? " <b class='tag'>changed</b>" : "")</td></tr>"
    end
    th = "<tr><th>Setting</th><th>This version ($(hesc_(id)))</th><th>Previous version$(previd === nothing ? "" : " ($(hesc_(previd)))")</th></tr>"
    log = reverse(decision_log())
    logrows = join(["<tr><td>$(hesc_(get(d, "date", "")))</td><td>$(hesc_(get(d, "area", "")))</td><td>$(hesc_(get(d, "decision", "")))</td><td>$(hesc_(get(d, "was", "")))</td><td>$(hesc_(get(d, "by", "")))</td><td><code>$(hesc_(get(d, "where", "")))</code></td></tr>" for d in log])
    nchg = prev === nothing ? 0 : count(r -> (pv = lookup(get(prev, "settings", Any[]), r[1]); pv !== nothing && pv != r[3]), settings)
    """<section class='decisions' id='decisions'><h2>Modelling decisions in this version</h2>
<p class='lead'>Version <b>$(hesc_(id))</b>. Every setting that decides what is fitted and controlled, read from <code>stages.toml</code> and <code>control_config.toml</code>, next to the previous Final Report version$(previd === nothing ? "" : " ($(hesc_(previd)); $nchg setting$(nchg == 1 ? "" : "s") changed)"). Every version is archived and never overwritten: <a href='history/index.html'>compare all versions</a>.</p>
$(provenance_html(provenance))
$stalehtml<div class='tscroll'><table>$th$(join([row(k, l, v, "settings") for (k, l, v) in settings]))</table></div>
<details><summary>Key results against the previous version ($(length(results)) rows)</summary><div class='tscroll'><table>$th$(join([row(k, l, v, "results") for (k, l, v) in results]))</table></div></details>
<details><summary>Decision log ($(length(log)) decisions, newest first; <code>decisions.toml</code>)</summary><div class='tscroll'><table><tr><th>Date</th><th>Area</th><th>Decision</th><th>Was</th><th>By</th><th>Where</th></tr>$logrows</table></div></details>
</section>"""
end

# the version history on the Final Report page itself (user, 2026-10-30): every archived version, newest
# first; a version opens to show its modelling decisions (changed ones highlighted against the version
# before), links to its archived report, and any two versions open side by side on history/index.html.
# Each version also says whether any of its steps were reused from an earlier run (provenance, user 2026-11-02).
function history_panel(id, settings; date=Dates.format(now(), "yyyy-mm-dd HH:MM"), commit=gitc("rev-parse", "--short", "HEAD"), provenance=nothing)
    vs = load_versions()
    cur = Dict{String,Any}("id" => id, "date" => date, "commit" => commit,
                           "settings" => [Dict("key" => k, "label" => l, "value" => v) for (k, l, v) in settings])
    provenance === nothing || (cur["provenance"] = provenance)
    any(v -> v["id"] == id, vs) || push!(vs, cur)
    sort!(vs; by=v -> v["id"])
    function changed(v, p)
        p === nothing && return "<span class='mut'>first recorded version</span>"
        a = get(v, "settings", Any[]); b = get(p, "settings", Any[])
        (isempty(a) || isempty(b)) && return "<span class='mut'>settings not recorded</span>"
        ch = [r["label"] for r in a if (pv = lookup(b, r["key"]); pv !== nothing && pv != r["value"])]
        new_ = [r["label"] for r in a if lookup(b, r["key"]) === nothing]
        parts = String[]
        isempty(ch) || push!(parts, "changed: " * join(hesc_.(ch), ", "))
        isempty(new_) || push!(parts, "new: " * join(hesc_.(new_), ", "))
        isempty(parts) ? "<span class='mut'>same settings (results or chapters differ)</span>" : join(parts, "; ")
    end
    function decisions(v, p; grp="settings")
        a = get(v, grp, Any[])
        isempty(a) && return "<p class='mut'>This version's $(grp == "settings" ? "settings were" : "key results were") not recorded.</p>"
        b = p === nothing ? Any[] : get(p, grp, Any[])
        trs = map(a) do r
            pv = lookup(b, r["key"]); ch = pv !== nothing && pv != r["value"]
            "<tr$(ch ? " class='chg'" : "")><td>$(hesc_(r["label"]))</td><td>$(hesc_(r["value"]))</td><td>$(pv === nothing ? "<span class='mut'>—</span>" : hesc_(pv))$(ch ? " <b class='tag'>changed</b>" : "")</td></tr>"
        end
        "<div class='tscroll'><table><tr><th>Setting</th><th>$(hesc_(v["id"]))</th><th>Version before$(p === nothing ? "" : " ($(hesc_(p["id"])))")</th></tr>$(join(trs))</table></div>"
    end
    rows = String[]
    for i in length(vs):-1:1
        v = vs[i]; p = i > 1 ? vs[i-1] : nothing; this = v["id"] == id
        open_ = this ? "<span class='tag'>this page</span>" : "<a href='history/$(hesc_(v["id"]))/report.html'>open this report</a>"
        pv_ = get(v, "provenance", nothing)
        ptag = pv_ === nothing || get(pv_, "recorded", false) !== true ? "<span class='mut'>provenance not recorded</span>" :
               get(pv_, "partly_reused", false) ? "<b class='tag'>partly reused</b> " * hesc_(get(pv_, "summary", "")) : "<span class='mut'>all steps new</span>"
        push!(rows, "<details class='ver$(this ? " this" : "")'><summary><b>$(hesc_(v["id"]))</b> <span class='mut'>$(hesc_(get(v, "date", ""))) · <code>$(hesc_(get(v, "commit", "")))</code></span> — $(changed(v, p)) · $ptag</summary>" *
                    "<div class='verbody'><p>$open_ · modelling decisions of this version:</p>$(decisions(v, p))" *
                    "<details><summary>Key results of this version</summary>$(decisions(v, p; grp="results"))</details></div></details>")
    end
    opts(sel) = join(["<option value='$(hesc_(v["id"]))'$(v["id"] == sel ? " selected" : "")>$(hesc_(v["id"]))</option>" for v in reverse(vs)])
    older = length(vs) > 1 ? vs[end-1]["id"] : id
    """<section class='decisions' id='history'><h2>Report history</h2>
<p class='lead'>Every Final Report version, newest first. Old versions are archived and never overwritten. Click a version to see its modelling decisions (highlighted where they changed from the version before) and open that report as it was. Each line ends with whether any step of that version was reused from an earlier run.</p>
<form class='cmp' action='history/index.html' method='get'><b>Compare side by side:</b> <select name='a'>$(opts(older))</select> with <select name='b'>$(opts(id))</select>
<input type='hidden' name='side' value='1'><button type='submit'>Open both</button> <a href='history/index.html'>all versions in one table</a></form>
$(join(rows, "\n"))
</section>"""
end

gitc(args...) = try readchomp(Cmd(`git -C $(@__DIR__) $args`)) catch; "" end

# archive a page under history/<id>/ with its version record
function archive_version(id, page; settings, results, chapters, source, commit=gitc("rev-parse", "--short", "HEAD"), note="", content_sha="", provenance=nothing)
    d = joinpath(HIST, id); mkpath(d)
    write(joinpath(d, "report.html"), rebase_archived(page))
    rec = Dict{String,Any}("id" => id, "date" => Dates.format(now(), "yyyy-mm-dd HH:MM"), "commit" => commit, "source" => source, "note" => note,
                           "content_sha" => content_sha, "report_sha" => bytes2hex(sha256(page)),
                           "settings" => [Dict("key" => k, "label" => l, "value" => v) for (k, l, v) in settings],
                           "results" => [Dict("key" => k, "label" => l, "value" => v) for (k, l, v) in results],
                           "chapters" => [Dict("title" => t, "stamp" => s) for (t, s) in chapters])
    provenance === nothing || (rec["provenance"] = provenance)
    open(io -> TOML.print(io, rec), joinpath(d, "version.toml"), "w")
    rec
end

# a file's modification time on the local clock (unix2datetime alone gives UTC, while version dates use now())
local_mtime(path) = unix2datetime(mtime(path)) + Minute(round(Int, (datetime2unix(now()) - time()) / 60))

# the site-nav script tag that add_site_nav puts into the live page after its hash was recorded
strip_site_nav(html) = replace(html, r"<script defer src=[\"'][^\"']*site-nav\.js[\"']></script>\n?" => "")

# a report.html that is in no version yet (built before the history existed, or by an older script) is
# archived as it is, before it is replaced
function archive_unrecorded(current)
    isfile(current) || return
    html = read(current, String)
    shas = (bytes2hex(sha256(html)), bytes2hex(sha256(strip_site_nav(html))))   # a page that only gained the nav bar is not new
    any(v -> any(sha -> get(v, "report_sha", "") == sha || sha in get(v, "report_shas", String[]), shas), load_versions()) && return
    sha = shas[1]
    id = Dates.format(local_mtime(current), "yyyy-mm-dd_HHMM") * "_unrecorded"
    page = read(current, String)
    archive_version(id, page; settings=Tuple{String,String,String}[], results=Tuple{String,String,String}[], chapters=Tuple{String,String}[],
                    source="archived before replacement; settings not recorded", commit="")
    # report_sha must be the live file's hash so it is recognised next time
    f = joinpath(HIST, id, "version.toml"); r = TOML.parsefile(f); r["report_sha"] = sha
    open(io -> TOML.print(io, r), f, "w")
end
