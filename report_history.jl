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

hfmt(x::Bool) = x ? "yes" : "no"
hfmt(x::Integer) = string(x)
hfmt(x::AbstractFloat) = isinteger(x) ? string(Int(x)) : string(x)
hfmt(x::AbstractString) = String(x)
hfmt(x::AbstractVector) = join(hfmt.(x), ", ")
hfmt(::Nothing) = "not set"

# the settings that define a version: (key, label, value), in display order
function report_settings(pkg=HPKG)
    sf, cf = joinpath(pkg, "stages.toml"), joinpath(pkg, "control_config.toml")
    S = isfile(sf) ? TOML.parsefile(sf) : Dict{String,Any}(); C = isfile(cf) ? TOML.parsefile(cf) : Dict{String,Any}()
    stg(id) = (i = findfirst(s -> get(s, "id", "") == id, get(S, "stage", Any[])); i === nothing ? Dict{String,Any}() : S["stage"][i])
    sel = get(S, "selection", Dict()); sch = get(C, "schedule", Dict()); can = get(C, "cancer", Dict())
    prog = get(C, "progression", Dict()); w = get(C, "weights", Dict())
    rule = get(sel, "rule_bic", nothing)
    wt = get(sel, "residual_weighting", "none")
    ex = get(stg("s2"), "exclude", Any[])
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
     ("kill_coupled", "Kill coupled to growth (Norton–Simon)", hfmt(get(stg("s2"), "kill_growth_coupled", false))),
     ("betaC_floor", "Resistant kill floor in coculture (× Stage 2b)", hfmt(get(stg("s3"), "betaC_floor", nothing))),
     ("ec50", "EC50, A2780 / A2780cis (µM)", isempty(lines) ? "not recorded" : "$(hfmt(get(get(lines, "A2780", Dict()), "ec50", nothing))) / $(hfmt(get(get(lines, "A2780cis", Dict()), "ec50", nothing)))"),
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
    for scen in ("A", "B")
        i = findfirst(r -> get(r, "environment", "") == "normal" && get(r, "scenario", "") == scen && get(r, "ratio", "") == "50:50" &&
                           get(r, "policy", "") == "MPC, measured doses only", cc)
        i === nothing && continue
        r = cc[i]
        push!(out, ("mpc:$scen", "14-day MPC schedule · scenario $scen, 50:50 (µM)", r["schedule_uM"]))
        push!(out, ("mpc_total:$scen", "14-day MPC N + C at the end · scenario $scen, 50:50", string(round(Int, parse(Float64, r["total_final"])))))
        haskey(r, "time_to_progression_day") && push!(out, ("mpc_ttp:$scen", "14-day MPC time to progression · scenario $scen, 50:50",
                                                            (v = tryparse(Float64, r["time_to_progression_day"]); v === nothing || !isfinite(v) ? "not progressed" : "day $(hfmt(round(v; digits=1)))")))
    end
    out
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
                   "href='history/" => "href='../", "href='report.html'" => "href='../../report.html'")
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
function decisions_panel(id, settings, results, prev; chapters=Tuple{String,String}[])
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
$stalehtml<div class='tscroll'><table>$th$(join([row(k, l, v, "settings") for (k, l, v) in settings]))</table></div>
<details><summary>Key results against the previous version ($(length(results)) rows)</summary><div class='tscroll'><table>$th$(join([row(k, l, v, "results") for (k, l, v) in results]))</table></div></details>
<details><summary>Decision log ($(length(log)) decisions, newest first; <code>decisions.toml</code>)</summary><div class='tscroll'><table><tr><th>Date</th><th>Area</th><th>Decision</th><th>Was</th><th>By</th><th>Where</th></tr>$logrows</table></div></details>
</section>"""
end

gitc(args...) = try readchomp(Cmd(`git -C $(@__DIR__) $args`)) catch; "" end

# archive a page under history/<id>/ with its version record
function archive_version(id, page; settings, results, chapters, source, commit=gitc("rev-parse", "--short", "HEAD"), note="", content_sha="")
    d = joinpath(HIST, id); mkpath(d)
    write(joinpath(d, "report.html"), rebase_archived(page))
    rec = Dict{String,Any}("id" => id, "date" => Dates.format(now(), "yyyy-mm-dd HH:MM"), "commit" => commit, "source" => source, "note" => note,
                           "content_sha" => content_sha, "report_sha" => bytes2hex(sha256(page)),
                           "settings" => [Dict("key" => k, "label" => l, "value" => v) for (k, l, v) in settings],
                           "results" => [Dict("key" => k, "label" => l, "value" => v) for (k, l, v) in results],
                           "chapters" => [Dict("title" => t, "stamp" => s) for (t, s) in chapters])
    open(io -> TOML.print(io, rec), joinpath(d, "version.toml"), "w")
    rec
end

# a report.html that is in no version yet (built before the history existed, or by an older script) is
# archived as it is, before it is replaced
function archive_unrecorded(current)
    isfile(current) || return
    sha = bytes2hex(open(sha256, current))
    any(v -> get(v, "report_sha", "") == sha, load_versions()) && return
    id = Dates.format(unix2datetime(mtime(current)), "yyyy-mm-dd_HHMM") * "_unrecorded"
    page = read(current, String)
    archive_version(id, page; settings=Tuple{String,String,String}[], results=Tuple{String,String,String}[], chapters=Tuple{String,String}[],
                    source="archived before replacement; settings not recorded", commit="")
    # report_sha must be the live file's hash so it is recognised next time
    f = joinpath(HIST, id, "version.toml"); r = TOML.parsefile(f); r["report_sha"] = sha
    open(io -> TOML.print(io, r), f, "w")
end
