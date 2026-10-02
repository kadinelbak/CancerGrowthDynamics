# The hard-coded data loaders of zoo_common.jl as they were before the tidy table (commit 3d1bd72b).
# Kept only for test/test_data.jl, which checks that data/measurements.csv reproduces them exactly.
module LegacyData
using CSV, DataFrames, XLSX
const ROOT     = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DATA_DIR = joinpath(ROOT, "New Datasets")
const COCULTURE_DOSE = 1.0

# ── data ───────────────────────────────────────────────────────────────────────

function load_csv_xy(fname)
    df = CSV.read(joinpath(DATA_DIR, fname), DataFrame; silencewarnings=true)
    rename!(df, Symbol.(strip.(replace.(string.(names(df)), "﻿" => ""))))
    filter!(r -> !ismissing(r.Day) && !ismissing(r[Symbol("Mean Cells")]), df)
    x = Float64.(df.Day); y = Float64.(df[!, Symbol("Mean Cells")])
    o = sortperm(x); x[o], y[o]
end

function load_xlsx_xy(fname, sheet)
    raw = XLSX.readxlsx(joinpath(DATA_DIR, fname))[sheet][:]
    header = Symbol.(strip.(string.(raw[1, :])))
    rows   = [Any[raw[i, j] for j in eachindex(header)] for i in 2:size(raw, 1)]
    df     = DataFrame([getindex.(rows, j) for j in eachindex(header)], header)
    filter!(r -> !ismissing(r.Day) && !ismissing(r[Symbol("Mean Cells")]), df)
    x = Float64.(df.Day); y = Float64.(df[!, Symbol("Mean Cells")])
    o = sortperm(x); x[o], y[o]
end

# Transient drops (below half of the last kept value by > 50 cells, with a later recovery to
# >= 80%) and exact repeats of an earlier day's value are excluded, as in the fitting report.
function flag_drops(y; frac=0.5, recover=0.8, floor=50.0)
    n = length(y); bad = falses(n); reason = fill("", n)
    for i in 2:n
        if y[i] > 0 && any(y[j] == y[i] for j in 1:i-1)
            bad[i] = true; reason[i] = "exact duplicate of an earlier day"
        end
    end
    prev = y[1]
    for i in 2:n
        bad[i] && continue
        rest = [y[j] for j in i+1:n if !bad[j]]
        later = isempty(rest) ? -Inf : maximum(rest)
        if y[i] < frac*prev && prev - y[i] > floor && later >= recover*prev
            bad[i] = true; reason[i] = "transient drop"
        else
            prev = y[i]
        end
    end
    bad, reason
end

const LINES = [(line="A2780", sym="N", untreated="A2780_Untreated_30k 1.csv", lr_sheet="1 - mono A2780 30k",
                treated=[(1.47, "A2780_IC25_1.47um_30k.csv"), (1.00, "A2780_IC50_1um_30k.csv"), (0.62, "A2780_IC75_0.62um_30k.csv")]),
               (line="A2780cis", sym="C", untreated="A2780cis_Untreated_30k.csv", lr_sheet="1 - mono A2780cis 30k",
                treated=[(1.47, "A2780cis_IC25_1.47um_30k.csv"), (1.00, "A2780cis_IC50_1um_30k.csv"), (0.62, "A2780cis_IC75_0.62um_30k.csv")])]
const LR_FILE = "ALL A2780-cis LR data.xlsx"
# A2780cis low-resource day 0 is the nominal seeding count, 5x below day 1; its day 0 is dropped
# and day 1 becomes the start (all days shift by -1), as in run_new_data_staged_fits.jl.
const SHIFT_FROM_DAY1 = Set(["A2780cis"])

untreated = Dict(L.line => (xy = load_csv_xy(L.untreated); (x=xy[1], y=xy[2])) for L in LINES)

lrmono = Dict{String,Any}()
for L in LINES
    x, y = load_xlsx_xy(LR_FILE, L.lr_sheet)
    if L.line in SHIFT_FROM_DAY1
        keep = x .>= 1; x, y = x[keep] .- 1, y[keep]
    end
    bad, why = flag_drops(y)
    lrmono[L.line] = (x=x, y=y, bad=bad, why=why)
end

treated = Dict(L.line => (doses=[d for (d, _) in L.treated], data=[load_csv_xy(f) for (_, f) in L.treated],
                          files=[f for (_, f) in L.treated]) for L in LINES)

# Ce0 (untreated) and Ce1 (1 µM) cocultures; env e holds sensitive then resistant
const COCULTURE_LAYOUT = [
    (ratio="25S75R", cond="Ce0", file="ce0_25S75R.xlsx", ssheet="25Sensitive_Untreated", rsheet="75Resistant_Untreated", dose=0.0),
    (ratio="50S50R", cond="Ce0", file="ce0_50S50R.xlsx", ssheet="50Sensitive_Untreated", rsheet="50Resistant_Untreated", dose=0.0),
    (ratio="75S25R", cond="Ce0", file="ce0_75S25R.xlsx", ssheet="75Sensitive_Untreated", rsheet="25Resistant_Untreated", dose=0.0),
    (ratio="25S75R", cond="Ce1", file="ce1_25S75R.xlsx", ssheet="25Sensitive_Treated",   rsheet="75Resistant_Treated",   dose=COCULTURE_DOSE),
    (ratio="50S50R", cond="Ce1", file="ce1_50S50R.xlsx", ssheet="50Sensitive_Treated",   rsheet="50Resistant_Treated",   dose=COCULTURE_DOSE),
    (ratio="75S25R", cond="Ce1", file="ce1_75S25R.xlsx", ssheet="75Sensitive_Treated",   rsheet="25Resistant_Treated",   dose=COCULTURE_DOSE),
]
co_data = [(S=load_xlsx_xy(e.file, e.ssheet), R=load_xlsx_xy(e.file, e.rsheet)) for e in COCULTURE_LAYOUT]

# low-resource cocultures, runs 1 and 2; sheet ratio "a-b" = sensitive:resistant
const LR_CC_LAYOUT = [(run=run, ratio=ratio, label=Dict("1-3" => "25S75R", "1-1" => "50S50R", "3-1" => "75S25R")[ratio])
                      for run in (1, 2) for ratio in ("1-3", "1-1", "3-1")]
lr_cc = NamedTuple[]
for e in LR_CC_LAYOUT
    sN, sC = "$(e.run) - A2780 30k $(e.ratio)", "$(e.run) - A2780cis 30k $(e.ratio)"
    xN, yN = load_xlsx_xy(LR_FILE, sN); xC, yC = load_xlsx_xy(LR_FILE, sC)
    (bN, _), (bC, _) = flag_drops(yN), flag_drops(yC)
    push!(lr_cc, merge(e, (sN=sN, sC=sC, xN=xN, yN=yN, bN=bN, xC=xC, yC=yC, bC=bC)))
end

# Stage 4: the combination culture, low resource (0.5% FBS) + 1 µM cisplatin, A2780 : A2780cis 1:1 at
# 30k (New Datasets/Validation; fitted as a stage at the user's request, 2026-09-29). Day 0 is the
# nominal seeding count (SD 0), 5–20× below day 1, so, as for the A2780cis low-resource monoculture,
# day 0 is dropped and day 1 becomes the start (all days shift by −1). The drug and stress states are
# started as they stand after one day of exposure (zoo_stage4.jl).
const COMBO_FILES = (N=joinpath("Validation", "ALL A2780 LR + cis(A2780 30k 1-1).csv"),
                     C=joinpath("Validation", "ALL A2780 LR + cis(A2780cis 30k 1-1).csv"))
const COMBO_DOSE = 1.0
combo = let
    (xN, yN), (xC, yC) = load_csv_xy(COMBO_FILES.N), load_csv_xy(COMBO_FILES.C)
    kN, kC = xN .>= 1, xC .>= 1
    xN, yN, xC, yC = xN[kN] .- 1, yN[kN], xC[kC] .- 1, yC[kC]
    (bN, _), (bC, _) = flag_drops(yN), flag_drops(yC)
    (xN=xN, yN=yN, bN=bN, xC=xC, yC=yC, bC=bC)
end


end
