# Builds the tidy measurement table data/measurements.csv from the raw sources in data/sources.toml.
#   julia --project=. data/build_tidy.jl
# One row per (source, day): every tag of the source, the day as measured and the model time t
# (day − time_origin), the value and its spread, and whether the point is excluded and why.
# Excluded points (before the time origin, or flagged by a QC rule the source lists) are kept, marked.

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using CSV, DataFrames, XLSX, TOML

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

# The two QC rules (qc_rules.toml), exactly as the fitting scripts applied them: duplicates first,
# then transient drops among the points that are not duplicates.
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
    bad, reason
end

function tidy_rows(src)
    for r in get(src, "qc", String[])
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
    kb, kr = qc_flags(y[keep], get(src, "qc", String[]))
    bad[keep] .= kb; why[keep] .= kr
    tags = (source_id=src["id"], experiment=src["experiment"], condition=src["condition"], line=src["line"], culture=src["culture"],
            ratio=string(src["ratio"]), run=src["run"], seeding=src["seeding"], fbs=Float64(src["fbs"]), schedule=src["schedule"],
            units=get(src, "units", "cells per image (mean of wells)"), file=src["file"], sheet=get(src, "sheet", ""))
    [merge(tags, (day=day[i], t=day[i] - t0, value=y[i], sd=sd[i], sem=sem[i], n_wells=nw[i], excluded=bad[i], exclude_reason=why[i]))
     for i in eachindex(y)]
end

function build_tidy()
    ids = [s["id"] for s in SOURCES]
    length(unique(ids)) == length(ids) || error("duplicate source ids in sources.toml")
    df = DataFrame(reduce(vcat, [tidy_rows(s) for s in SOURCES]))
    CSV.write(OUT, df)
    println("Wrote $OUT: $(nrow(df)) rows from $(length(SOURCES)) sources, $(count(df.excluded)) excluded")
    df
end

abspath(PROGRAM_FILE) == (@__FILE__) && build_tidy()
