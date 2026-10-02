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

const DATA_TABLE = joinpath(@__DIR__, "data", "measurements.csv")
isfile(DATA_TABLE) || error("data/measurements.csv is missing: run `julia --project=. data/build_tidy.jl` first")
const MEAS      = CSV.read(DATA_TABLE, DataFrame; stringtype=String, types=Dict(:ratio => String, :exclude_reason => String))
const SOURCES   = TOML.parsefile(joinpath(@__DIR__, "data", "sources.toml"))["source"]
const SCHEDULES = TOML.parsefile(joinpath(@__DIR__, "data", "schedules.toml"))

# source ids whose tags match every keyword (in sources.toml order)
select_sources(; tags...) = [s["id"] for s in SOURCES if all(string(get(s, string(k), "")) == string(v) for (k, v) in tags)]
source(id) = SOURCES[findfirst(s -> s["id"] == id, SOURCES)]
only_source(; tags...) = (ids = select_sources(; tags...); length(ids) == 1 ? ids[1] : error("expected one source for $tags, found $ids"))

# one source's series from the time origin on, sorted by t
function series(id)
    d = sort(MEAS[(MEAS.source_id .== id) .& (MEAS.t .>= 0), :], :t)
    nrow(d) > 0 || error("no rows for source $id in measurements.csv")
    (x=Vector{Float64}(d.t), y=Vector{Float64}(d.value), bad=Vector{Bool}(d.excluded),
     why=String[coalesce(w, "") for w in d.exclude_reason], se=Vector{Union{Missing,Float64}}(d.sem))
end

# a schedule as (breakpoint times, doses); a constant schedule is a single step at t = 0
schedule_steps(name) = (s = SCHEDULES[name]; (Float64[p[1] for p in s], Float64[p[2] for p in s]))
function constant_dose(name)
    ts, ds = schedule_steps(name)
    length(ds) == 1 || error("schedule $name is not constant")
    ds[1]
end

const LINE_SYM = Dict("A2780" => "N", "A2780cis" => "C")
const LINES = [(line=L, sym=LINE_SYM[L]) for L in ("A2780", "A2780cis")]

untreated = Dict(L.line => (s = series(only_source(condition="untreated", line=L.line)); (x=s.x, y=s.y)) for L in LINES)

lrmono = Dict(L.line => (s = series(only_source(condition="low_resource", culture="mono", line=L.line)); (x=s.x, y=s.y, bad=s.bad, why=s.why))
              for L in LINES)

treated = Dict(L.line => begin
                   ids = select_sources(condition="treated", line=L.line)
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

# Stage 4: low resource + 1 µM cisplatin, 1:1 (day 1 is t = 0, see sources.toml)
const COMBO_SOURCES = (N=only_source(condition="combination", line="A2780"), C=only_source(condition="combination", line="A2780cis"))
const COMBO_DOSE = constant_dose(source(COMBO_SOURCES.N)["schedule"])
combo = let N = series(COMBO_SOURCES.N), C = series(COMBO_SOURCES.C)
    (xN=N.x, yN=N.y, bN=N.bad, xC=C.x, yC=C.y, bC=C.bad)
end

# intermittent dosing (validation only): per line and schedule
intermittent = Dict(L.line => Dict(source(i)["schedule"] => (s = series(i); (x=s.x, y=s.y, se=s.se, schedule=source(i)["schedule"]))
                                   for i in select_sources(condition="intermittent", line=L.line))
                    for L in LINES)

# SHA-256 of every raw file and of the tidy table, for the run manifest
function data_checksums()
    files = unique([s["file"] for s in SOURCES])
    Dict(vcat([f => bytes2hex(open(sha256, joinpath(ROOT, f))) for f in files],
              ["PKG/data/measurements.csv" => bytes2hex(open(sha256, DATA_TABLE))]))
end
