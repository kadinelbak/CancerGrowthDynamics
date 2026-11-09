# Builds the tidy measurement table data/measurements.csv from the raw sources in data/sources.toml.
#   julia --project=. data/build_tidy.jl
# One row per (source, day): every tag of the source, the day as measured and the model time t
# (day − time_origin), the value and its spread, and whether the point is excluded and why.
# Excluded points (before the time origin, or flagged by a QC rule the source lists) are kept, marked.

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using CSV, DataFrames, XLSX, TOML, Statistics

const PKG_DIR  = normpath(joinpath(@__DIR__, ".."))
const ROOT     = normpath(joinpath(PKG_DIR, "..", ".."))
# SOURCES_TOML / QC_TOML / DATA_TABLE replace the sources, the QC rules and the output table
# (sensitivity runs, e.g. without the day-0 or outlier exclusions); defaults are the files here
const SOURCES  = TOML.parsefile(get(ENV, "SOURCES_TOML", joinpath(@__DIR__, "sources.toml")))["source"]
const QC       = TOML.parsefile(get(ENV, "QC_TOML", joinpath(@__DIR__, "qc_rules.toml")))
const SCHED    = TOML.parsefile(joinpath(@__DIR__, "schedules.toml"))
const OUT      = get(ENV, "DATA_TABLE", joinpath(@__DIR__, "measurements.csv"))

clean(s) = strip(replace(string(s), "﻿" => ""))

# raw columns of one source as a DataFrame with stripped header names
function read_source(src)
    path = joinpath(ROOT, src["file"])
    if endswith(lowercase(path), ".xlsx")
        raw = XLSX.readxlsx(path)[src["sheet"]][:]
        header = Symbol.(clean.(raw[1, :]))
        rows = [Any[raw[i, j] for j in eachindex(header)] for i in 2:size(raw, 1)]
        DataFrame([getindex.(rows, j) for j in eachindex(header)], header)
    else
        df = CSV.read(path, DataFrame; silencewarnings=true)
        rename!(df, Symbol.(clean.(names(df))))
    end
end

num(v) = v === missing || v === nothing ? missing : Float64(v)

# The QC rules (qc_rules.toml): duplicates first, then transient drops and spikes among the points that are
# not duplicates.
function qc_flags(y, rules)
    n = length(y); bad = falses(n); reason = fill("", n)
    if "exact_duplicate" in rules
        for i in 2:n
            if y[i] > 0 && any(y[j] == y[i] for j in 1:i-1)
                bad[i] = true; reason[i] = QC["exact_duplicate"]["description"]
            end
        end
    end
    if "transient_drop" in rules
        R = QC["transient_drop"]; frac, recover, floor = R["frac"], R["recover"], R["floor"]
        prev = y[1]
        for i in 2:n
            bad[i] && continue
            rest = [y[j] for j in i+1:n if !bad[j]]
            later = isempty(rest) ? -Inf : maximum(rest)
            if y[i] < frac * prev && prev - y[i] > floor && later >= recover * prev
                bad[i] = true; reason[i] = R["description"]
            else
                prev = y[i]
            end
        end
    end
    if "transient_spike" in rules
        R = QC["transient_spike"]; frac, recover, floor = R["frac"], R["recover"], R["floor"]
        prev = y[1]
        for i in 2:n
            bad[i] && continue
            nxt = something(findfirst(k -> k > i && !bad[k], 1:n), 0)
            if nxt > 0 && y[i] > prev / frac && y[i] - prev > floor && y[nxt] <= prev / recover
                bad[i] = true; reason[i] = R["description"]
            else
                prev = y[i]
            end
        end
    end
    bad, reason
end

# the rules of a source: those applied to every source plus any it lists itself
qc_rules_of(src) = unique(vcat(String.(get(QC, "apply_to_all", String[])), String.(get(src, "qc", String[]))))

function tidy_rows(src)
    for r in qc_rules_of(src)
        haskey(QC, r) || error("source $(src["id"]): unknown QC rule $r")
    end
    haskey(SCHED, src["schedule"]) || error("source $(src["id"]): unknown schedule $(src["schedule"])")
    df = read_source(src)
    vcol = Symbol(get(src, "value_col", "Mean Cells"))
    sdc, sec, nc = Symbol(get(src, "sd_col", "SD Cells")), Symbol(get(src, "sem_col", "SEM Cells")), Symbol(get(src, "n_col", "N Samples"))
    df = filter(r -> !ismissing(r.Day) && !ismissing(r[vcol]), df)
    day = Float64.(df.Day); y = Float64.(df[!, vcol])
    col(c) = hasproperty(df, c) ? num.(df[!, c]) : fill(missing, nrow(df))
    sd, sem, nw = col(sdc), col(sec), col(nc)
    o = sortperm(day); day, y, sd, sem, nw = day[o], y[o], sd[o], sem[o], nw[o]
    t0 = Float64(get(src, "time_origin", 0))
    keep = day .>= t0
    bad = falses(length(y)); why = fill("", length(y))
    for i in findall(.!keep)
        bad[i] = true; why[i] = "before the time origin (day $(Int(t0))): " * get(src, "time_origin_reason", "")
    end
    kb, kr = qc_flags(y[keep], qc_rules_of(src))
    bad[keep] .= kb; why[keep] .= kr
    # points left out by hand (sources.toml exclude_days + exclude_reason): kept in the table, marked
    # "left out (user ...)" so zoo_data.jl can put them back with ZOO_NO_EXCLUDE=1
    for d in Float64.(get(src, "exclude_days", Float64[])), i in findall(==(d), day)
        bad[i] && continue
        bad[i] = true; why[i] = "left out (user, $(get(src, "exclude_reason", "no reason given")))"
    end
    tags = (source_id=src["id"], experiment=src["experiment"], condition=src["condition"], line=src["line"], culture=src["culture"],
            ratio=string(src["ratio"]), run=src["run"], seeding=src["seeding"], fbs=Float64(src["fbs"]), schedule=src["schedule"],
            units=get(src, "units", "cells per image (mean of wells)"), file=src["file"], sheet=get(src, "sheet", ""))
    [merge(tags, (day=day[i], t=day[i] - t0, value=y[i], sd=sd[i], sem=sem[i], n_wells=nw[i], excluded=bad[i], exclude_reason=why[i]))
     for i in eachindex(y)]
end

# Terminal collapse (qc_rules.toml [terminal_collapse], fix review X4, user 2026-11-09): on one of the last
# `last_days` measured days of a coculture well, one line falls below `frac` × the median of its previous `window`
# kept points while the other line rises above that median / `frac` (each by more than `floor` cells). Opposite
# moves of the two lines of one well on the same late day are an image artefact (well detachment, focus), not
# biology: both points are excluded, kept in the table and marked. Monocultures are never flagged by this rule.
function terminal_collapse!(df, well)
    haskey(QC, "terminal_collapse") || return df
    R = QC["terminal_collapse"]; frac, floor, win, lastn = R["frac"], R["floor"], Int(R["window"]), Int(R["last_days"])
    for w in unique(well[i] for i in eachindex(well) if df.culture[i] != "mono")
        rows = [i for i in eachindex(well) if well[i] == w]
        lines = unique(df.line[rows]); length(lines) == 2 || continue
        days = sort(unique(df.day[rows])); late = days[max(end - lastn + 1, 1):end]
        ref(l, d) = (k = [i for i in rows if df.line[i] == l && df.day[i] < d && !df.excluded[i] && df.t[i] > 0];
                     sort!(k; by=i -> df.day[i]); length(k) < win ? NaN : median(df.value[k[end-win+1:end]]))
        for d in late
            at = Dict(l => findfirst(i -> df.line[i] == l && df.day[i] == d, rows) for l in lines)
            any(isnothing, values(at)) && continue
            i1, i2 = rows[at[lines[1]]], rows[at[lines[2]]]
            (df.excluded[i1] || df.excluded[i2]) && continue
            m1, m2 = ref(lines[1], d), ref(lines[2], d)
            (isnan(m1) || isnan(m2)) && continue
            drop(i, m) = df.value[i] < frac * m && m - df.value[i] > floor
            rise(i, m) = df.value[i] > m / frac && df.value[i] - m > floor
            if (drop(i1, m1) && rise(i2, m2)) || (rise(i1, m1) && drop(i2, m2))
                for i in (i1, i2)
                    df.excluded[i] = true; df.exclude_reason[i] = R["description"]
                end
            end
        end
    end
    df
end

function build_tidy()
    ids = [s["id"] for s in SOURCES]
    length(unique(ids)) == length(ids) || error("duplicate source ids in sources.toml")
    df = DataFrame(reduce(vcat, [tidy_rows(s) for s in SOURCES]))
    # a QC artefact in one line of a coculture well is an artefact of the image: flag that day in both lines
    qc_reason(r) = !r.excluded || startswith(r.exclude_reason, "before the time origin") ? "" : r.exclude_reason
    well = [(r.experiment, r.condition, r.culture, r.ratio, r.run, r.seeding, r.schedule, r.fbs, r.file) for r in eachrow(df)]
    for i in eachindex(well)
        df.culture[i] == "mono" && continue
        df.excluded[i] && continue
        why = [qc_reason(df[j, :]) for j in eachindex(well) if j != i && well[j] == well[i] && df.day[j] == df.day[i]]
        k = findfirst(!isempty, why)
        k === nothing && continue
        df.excluded[i] = true; df.exclude_reason[i] = why[k] * " (flagged in the other line of the same well)"
    end
    terminal_collapse!(df, well)
    CSV.write(OUT, df)
    println("Wrote $OUT: $(nrow(df)) rows from $(length(SOURCES)) sources, $(count(df.excluded)) excluded")
    df
end

abspath(PROGRAM_FILE) == (@__FILE__) && build_tidy()
