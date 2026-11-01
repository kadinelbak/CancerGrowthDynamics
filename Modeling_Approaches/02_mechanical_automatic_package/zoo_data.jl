# Data for the zoo pipeline, selected by tag from the tidy table data/measurements.csv
# (built from data/sources.toml by data/build_tidy.jl). Returns the structures the fitting code uses:
#   untreated[line] = (x, y)                       Stage 1
#   lrmono[line]    = (x, y, bad, why)             Stage 2a (low-resource monocultures)
#   treated[line]   = (doses, data=[(x, y)], files) Stage 2b, one series per dose
#   co_data[e]      = (S=(x, y), R=(x, y))         Stage 3a, environments in COCULTURE_LAYOUT order
#   lr_cc[e]        = (run, ratio, label, xN, yN, bN, xC, yC, bC, …)   Stage 3b
#   combo           = (xN, yN, bN, xC, yC, bC)     Stage 4
#   intermittent[line][schedule] = (x, y, se, schedule)   validation only
# x is the model time t (day − time origin); bad marks points a QC rule excluded.
# Include from zoo_common.jl (needs CSV, DataFrames, TOML).

# Configuration files; each can be replaced for a sensitivity run through an environment variable:
#   DATA_TABLE (tidy table), SOURCES_TOML (data/sources.toml), STAGES_TOML (stages.toml)
const DATA_TABLE   = get(ENV, "DATA_TABLE", joinpath(@__DIR__, "data", "measurements.csv"))
const SOURCES_FILE = get(ENV, "SOURCES_TOML", joinpath(@__DIR__, "data", "sources.toml"))
const STAGES_FILE  = get(ENV, "STAGES_TOML", joinpath(@__DIR__, "stages.toml"))
isfile(DATA_TABLE) || error("$DATA_TABLE is missing: run `julia --project=. data/build_tidy.jl` first")
const MEAS      = CSV.read(DATA_TABLE, DataFrame; ntasks=1, stringtype=String, types=Dict(:ratio => String, :exclude_reason => String))
const SOURCES   = TOML.parsefile(SOURCES_FILE)["source"]
const SCHEDULES = TOML.parsefile(joinpath(@__DIR__, "data", "schedules.toml"))

# source ids whose tags match every keyword (in sources.toml order)
select_sources(; tags...) = [s["id"] for s in SOURCES if all(string(get(s, string(k), "")) == string(v) for (k, v) in tags)]
source(id) = SOURCES[findfirst(s -> s["id"] == id, SOURCES)]
only_source(; tags...) = (ids = select_sources(; tags...); length(ids) == 1 ? ids[1] : error("expected one source for $tags, found $ids"))

# Nominal seeding of a source's line (user, 2026-10-07): 100 cells per well at 30k, 66.7 at 20k, times the
# line's share in a coculture (25S75R → 25 sensitive, 75 resistant). Day 0 in the files is this nominal value
# (99.8 with SD 0), not a measurement.
function nominal_seed(src)
    total = 100.0 * Float64(src["seeding"]) / 30000
    r = string(src["ratio"])
    r == "100" && return total
    m = match(r"^(\d+)S(\d+)R$", r); m === nothing && error("ratio $r of $(src["id"]) not understood")
    total * parse(Float64, src["line"] == "A2780" ? m[1] : m[2]) / 100
end

# one source's series from the time origin on, sorted by t; nom = its nominal seeding
function series(id)
    d = sort(MEAS[(MEAS.source_id .== id) .& (MEAS.t .>= 0), :], :t)
    nrow(d) > 0 || error("no rows for source $id in measurements.csv")
    (x=Vector{Float64}(d.t), y=Vector{Float64}(d.value), bad=Vector{Bool}(d.excluded),
     why=String[coalesce(w, "") for w in d.exclude_reason], se=Vector{Union{Missing,Float64}}(d.sem),
     sd=Vector{Union{Missing,Float64}}(d.sd), nom=nominal_seed(source(id)))
end

# Starting value of a series (user, 2026-10-07): the nominal seeding, unless the first measured count is
# more than X0_FOLD^t (t = its day, normally 1) times above or below it. No growth or kill moves a culture
# 3-fold in a day (A2780 doubles in about a day and the drug and low-resource effects build up over days), so
# such a series was seeded or imaged at a different effective density and its starting value is estimated
# (a free parameter in that stage's fits, counted in the BIC), between nom/10 and 3 × the larger of nom and
# the first count.
const X0_FOLD = 3.0
function x0_info(s)
    k = findfirst(i -> s.x[i] > 0 && !s.bad[i], eachindex(s.x))
    k === nothing && error("series without measurements after day 0")
    y1 = s.y[k]; nom = s.nom; f = X0_FOLD^s.x[k]
    est = y1 > f * nom || y1 < nom / f
    lo, hi = nom / 10, 3 * max(nom, y1)
    (nom=nom, est=est, y1=y1, day1=s.x[k], start=clamp(sqrt(nom * y1), 1.01lo, 0.99hi), lo=lo, hi=hi)
end

# a schedule as (breakpoint times, doses); a constant schedule is a single step at t = 0
schedule_steps(name) = (s = SCHEDULES[name]; (Float64[p[1] for p in s], Float64[p[2] for p in s]))
# the model input of a schedule: a number if constant, otherwise a Schedule (zoo_models.jl)
dose_input(name) = ((ts, ds) = schedule_steps(name); length(ds) == 1 ? ds[1] : Schedule(ts, ds))
function constant_dose(name)
    ts, ds = schedule_steps(name)
    length(ds) == 1 || error("schedule $name is not constant")
    ds[1]
end

const STAGE_CFG = TOML.parsefile(STAGES_FILE)
# the two lines, their symbols and measured IC50 (stages.toml [lines])
const LINE_CFG = STAGE_CFG["lines"]
const LINES = [(line=L, sym=LINE_CFG[L]["sym"], ec50=Float64(LINE_CFG[L]["ec50"])) for L in ("A2780", "A2780cis")]
const LINE_SYM = Dict(L.line => L.sym for L in LINES)
const EC50_OF_SYM = Dict(L.sym => L.ec50 for L in LINES)       # EC50 by line symbol (N, C)

untreated = Dict(L.line => (s = series(only_source(condition="untreated", line=L.line)); (x=s.x, y=s.y)) for L in LINES)

lrmono = Dict(L.line => (s = series(only_source(condition="low_resource", culture="mono", line=L.line)); (x=s.x, y=s.y, bad=s.bad, why=s.why))
              for L in LINES)

# Series a stage leaves out of its fit: the stage's `exclude` list in stages.toml, each entry a set of tags
# that must all match (`dose` = the constant dose in µM). ZOO_NO_EXCLUDE=1 ignores every list (the
# original fit, kept as a sensitivity variant: user, 2026-10-27, dropped the A2780cis 1 µM monoculture).
const USE_EXCLUDE = get(ENV, "ZOO_NO_EXCLUDE", "0") != "1"
function excluded(st, s)
    USE_EXCLUDE || return false
    any(get(st, "exclude", Any[])) do e
        all(collect(keys(e))) do k
            k == "dose" ? (string(get(s, "condition", "")) == "treated" && constant_dose(s["schedule"]) ≈ Float64(e[k])) : string(get(s, k, "")) == string(e[k])
        end
    end
end
# the excluded sources of every stage, for the reports: [(stage, id, line, file)]
excluded_sources() = [(stage=st["id"], id=s["id"], line=s["line"], file=basename(s["file"]))
                      for st in STAGE_CFG["stage"] for s in SOURCES if haskey(st, "data") && excluded(st, s)]

treated = Dict(L.line => begin
                   st2 = STAGE_CFG["stage"][findfirst(s -> s["id"] == "s2", STAGE_CFG["stage"])]
                   ids = [i for i in select_sources(condition="treated", line=L.line) if !excluded(st2, source(i))]
                   (doses=[constant_dose(source(i)["schedule"]) for i in ids], data=[(s = series(i); (s.x, s.y)) for i in ids],
                    files=[basename(source(i)["file"]) for i in ids])
               end for L in LINES)

# Ce0 (untreated) and Ce1 (1 µM) cocultures; environment e holds sensitive then resistant
const COCULTURE_LAYOUT = [begin
                              sS = only_source(condition=cond, culture="co", ratio=ratio, line="A2780")
                              sR = only_source(condition=cond, culture="co", ratio=ratio, line="A2780cis")
                              (ratio=ratio, cond=cond, file=basename(source(sS)["file"]), ssheet=source(sS)["sheet"], rsheet=source(sR)["sheet"],
                               dose=constant_dose(source(sS)["schedule"]), sources=(sS, sR))
                          end for cond in ("Ce0", "Ce1") for ratio in ("25S75R", "50S50R", "75S25R")]
co_data = [(S=(s = series(e.sources[1]); (s.x, s.y)), R=(s = series(e.sources[2]); (s.x, s.y))) for e in COCULTURE_LAYOUT]

# low-resource cocultures, runs 1 and 2; `ratio` is the sheet's sensitive:resistant notation
const RATIO_SHEET = Dict("25S75R" => "1-3", "50S50R" => "1-1", "75S25R" => "3-1")
lr_cc = NamedTuple[]
for run in (1, 2), label in ("25S75R", "50S50R", "75S25R")
    iN = only_source(condition="low_resource", culture="co", run=run, ratio=label, line="A2780")
    iC = only_source(condition="low_resource", culture="co", run=run, ratio=label, line="A2780cis")
    N, C = series(iN), series(iC)
    push!(lr_cc, (run=run, ratio=RATIO_SHEET[label], label=label, sN=source(iN)["sheet"], sC=source(iC)["sheet"],
                  xN=N.x, yN=N.y, bN=N.bad, xC=C.x, yC=C.y, bC=C.bad))
end

# Stage 4: low resource + 1 µM cisplatin, 1:1
const COMBO_SOURCES = (N=only_source(condition="combination", line="A2780"), C=only_source(condition="combination", line="A2780cis"))
const COMBO_DOSE = constant_dose(source(COMBO_SOURCES.N)["schedule"])
combo = let N = series(COMBO_SOURCES.N), C = series(COMBO_SOURCES.C)
    (xN=N.x, yN=N.y, bN=N.bad, xC=C.x, yC=C.y, bC=C.bad)
end

# intermittent dosing (validation only): per line and schedule
intermittent = Dict(L.line => Dict(source(i)["schedule"] => (s = series(i); (x=s.x, y=s.y, se=s.se, schedule=source(i)["schedule"]))
                                   for i in select_sources(condition="intermittent", line=L.line))
                    for L in LINES)

# ── data of a fitting stage (stages.toml) ──────────────────────────────────────

selection_cfg(k, default) = get(get(STAGE_CFG, "selection", Dict()), k, default)

# Per-point residual scale of a series under stages.toml [selection] residual_weighting:
#   "none"  (default) unweighted squared error; returns nothing
#   "sd"    the replicate SD of each point, floored at the median positive SD of the series
#           (day 0 has SD 0)
#   "power" max(y, 50)^weight_power: the spread grows like a power of the count
# Weighted BICs compare models within one scheme, not across schemes.
function point_scale(s)
    mode = selection_cfg("residual_weighting", "none")
    mode == "none" && return nothing
    if mode == "sd"
        sd = Float64[coalesce(v, NaN) for v in s.sd]
        pos = sort(filter(v -> isfinite(v) && v > 0, sd))
        fl = isempty(pos) ? 1.0 : pos[cld(length(pos), 2)]
        return [isfinite(v) ? max(v, fl) : fl for v in sd]
    elseif mode == "power"
        pw = Float64(selection_cfg("weight_power", 0.55))
        return [max(y, 50.0)^pw for y in s.y]
    end
    error("stages.toml [selection] residual_weighting must be none, sd or power, not $mode")
end
# the fitted points of a series: QC-excluded points and day 0 (the seeding, which is the initial condition,
# not a measurement) dropped; with their scales w when weighting is on
function keep(s)
    w = point_scale(s)
    k = .!s.bad .& (s.x .> 0)
    (x=s.x[k], y=s.y[k], w=w === nothing ? nothing : w[k])
end
# a dataset for the fits; carries w only when weighting is on, so unweighted runs are unchanged
dataset(k, idx) = k.w === nothing ? (x=k.x, y=k.y, state_index=idx) : (x=k.x, y=k.y, state_index=idx, w=k.w)
const STAGES = STAGE_CFG["stage"]
stage_cfg(id) = STAGES[findfirst(s -> s["id"] == id, STAGES)]
ur_of(fbs) = (10 - fbs) / (10 - 0.05)                      # low-resource level u_R from the serum %, as U_R

# sources matching a stage's tag selector (a list value matches any entry), in sources.toml order
function stage_sources(st; line=nothing)
    sel = st["data"]
    ok(s) = all(k -> (v = sel[k]; w = string(get(s, k, "")); v isa AbstractVector ? w in string.(v) : w == string(v)), collect(keys(sel))) &&
            (line === nothing || s["line"] == line)
    [s for s in SOURCES if ok(s) && !excluded(st, s)]
end
env_inputs(s) = (D=dose_input(s["schedule"]), R=ur_of(Float64(s["fbs"])))
# coculture environments: distinct (condition, ratio, run), each with its sensitive and resistant source
function cc_envs(st)
    srcs = stage_sources(st)
    ks = unique([(s["condition"], string(s["ratio"]), s["run"]) for s in srcs])
    [begin
         pick(L) = only(s for s in srcs if (s["condition"], string(s["ratio"]), s["run"]) == k && s["line"] == L)
         sN, sC = pick("A2780"), pick("A2780cis")
         env_inputs(sN) == env_inputs(sC) || error("environment $k: the two lines have different inputs")
         (key=k, inputs=env_inputs(sN), N=series(sN["id"]), C=series(sC["id"]))
     end for k in ks]
end

# the tags of each dataset of a fitted group, in the group's dataset order (one per line per environment)
function dataset_tags(gid)
    st = stage_cfg(split(gid, "-")[1])
    if st["kind"] in ("growth", "mechanism")
        line = only(L.line for L in LINES if endswith(gid, "-" * lowercase(L.line)))
        return stage_sources(st; line=line)
    end
    st["kind"] == "coculture" || error("dataset_tags: no per-dataset tags for $gid")
    srcs = stage_sources(st)
    reduce(vcat, [[only(s for s in srcs if (s["condition"], string(s["ratio"]), s["run"]) == e.key && s["line"] == L) for L in ("A2780", "A2780cis")]
                  for e in cc_envs(st)])
end

# SHA-256 of every raw file and of the tidy table, for the run manifest
function data_checksums()
    files = unique([s["file"] for s in SOURCES])
    Dict(vcat([f => bytes2hex(open(sha256, joinpath(ROOT, f))) for f in files],
              ["PKG/data/measurements.csv" => bytes2hex(open(sha256, DATA_TABLE))]))
end
