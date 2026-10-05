# The whole pipeline, data to Final Report, in one command. Each step runs in its own Julia process
# with its log in outputs/logs/<step>.log; a run manifest (outputs/run_manifest.toml) records the
# code version, Julia and package versions, threads, data checksums and every step's time and exit code.
#
#   julia run_pipeline.jl                       everything from the raw data (about 5–6 h on 20 threads)
#   julia run_pipeline.jl --threads 12          threads for the fitting steps (default: all cores − 2)
#   julia run_pipeline.jl --from analyses       start at a step (earlier outputs are reused)
#   julia run_pipeline.jl --only report14       one step (see --list)
#   julia run_pipeline.jl --resume              reuse unchanged fits and bootstraps (ZOO_RESUME, ZOO_REUSE_BOOT)
#                                               and the control checkpoints instead of recomputing
#   julia run_pipeline.jl --list                the steps
#
# Configuration lives in data/sources.toml (data and their tags), data/schedules.toml (dose schedules),
# data/qc_rules.toml (exclusions), stages.toml (stages, candidates, line properties), validation.toml
# (hold-outs) and control_config.toml (optimal control).

using Dates, TOML, SHA

const PKG = @__DIR__
const ROOT = normpath(joinpath(PKG, "..", ".."))
const LOGS = joinpath(PKG, "outputs", "logs")
const JL = Base.julia_cmd()

args = copy(ARGS)
flag(f) = (i = findfirst(==(f), args); i === nothing ? false : (deleteat!(args, i); true))
function opt(f, default)
    i = findfirst(==(f), args); i === nothing && return default
    v = args[i+1]; deleteat!(args, i:i+1); v
end
const RESUME = flag("--resume")
const LIST = flag("--list")
const THREADS = parse(Int, opt("--threads", string(max(Sys.CPU_THREADS - 2, 1))))
const FROM = opt("--from", "")
const ONLY = opt("--only", "")
isempty(args) || error("unknown arguments: $(join(args, " "))")

# (name, what it does, working directory, script and arguments, threads?, extra environment)
redo = RESUME ? Dict{String,String}() : Dict("REDO" => "1")
const STEPS = [
    ("data",       "tidy data table from the raw sources",                  PKG,  ["data/build_tidy.jl"],             false, Dict{String,String}()),
    ("checks",     "the tidy table reproduces the source files", PKG, ["test/test_data.jl"],             false, Dict{String,String}()),
    ("fits",       "staged model-zoo fits and bootstrap",                   PKG,  ["run_zoo_fits.jl"],                true,  RESUME ? Dict("ZOO_RESUME" => "1", "ZOO_REUSE_BOOT" => "1") : Dict{String,String}()),
    ("analyses",   "bootstrap comparison, sensitivity, validation",         PKG,  ["zoo_analyses.jl"],                true,  Dict{String,String}()),
    ("control14",  "14-day optimal control (dose)",                         PKG,  ["run_control.jl"],                 true,  redo),
    ("report14",   "14-day control report",                                 PKG,  ["control_report.jl"],              false, Dict{String,String}()),
    ("control_long", "long optimal control, dose + fasting (control_config long_days)", PKG, ["run_control_fasting.jl"], true, redo),
    ("report_long",  "long control report",                                 PKG,  ["control_fasting_report.jl"],      false, Dict{String,String}()),
    ("tests",      "control, stage and validation tests on the new results", PKG, ["test/run_tests.jl"],              true,  Dict{String,String}()),
    ("final",      "Final Report (all chapters)",                           ROOT, ["build_final_report.jl"],          false, Dict{String,String}()),
]
stepnames = [s[1] for s in STEPS]
if LIST
    foreach(s -> println(rpad(s[1], 12), s[2]), STEPS); exit(0)
end
for s in (FROM, ONLY)
    isempty(s) || s in stepnames || error("unknown step $s; steps: $(join(stepnames, ", "))")
end
todo = !isempty(ONLY) ? [ONLY] : isempty(FROM) ? stepnames : stepnames[findfirst(==(FROM), stepnames):end]

# ── run manifest ───────────────────────────────────────────────────────────────
git(cmd...) = try readchomp(Cmd(`git -C $ROOT $cmd`)) catch; "unknown" end
function gpe_pin()
    m = TOML.parsefile(joinpath(PKG, "Manifest.toml"))
    d = only(get(m, "deps", m)["GrowthParameterEstimation"])
    Dict(k => get(d, k, "") for k in ("version", "repo-url", "repo-rev", "git-tree-sha1"))
end
function data_files()
    srcs = TOML.parsefile(joinpath(PKG, "data", "sources.toml"))["source"]
    files = sort(unique(vcat([s["file"] for s in srcs])))
    Dict(vcat([f => (isfile(joinpath(ROOT, f)) ? bytes2hex(open(sha256, joinpath(ROOT, f))) : "missing") for f in files],
              ["PKG/data/measurements.csv" => bytes2hex(open(sha256, joinpath(PKG, "data", "measurements.csv")))]))
end
manifest = Dict{String,Any}(
    "started" => string(now()), "steps_requested" => todo, "threads" => THREADS, "resume" => RESUME,
    "git" => Dict("commit" => git("rev-parse", "HEAD"), "branch" => git("rev-parse", "--abbrev-ref", "HEAD"),
                  "uncommitted_changes" => !isempty(git("status", "--porcelain"))),
    "julia" => string(VERSION), "GrowthParameterEstimation" => gpe_pin(), "steps" => Dict{String,Any}())
const MANIFEST = joinpath(PKG, "outputs", "run_manifest.toml")
save_manifest() = open(io -> TOML.print(io, manifest), MANIFEST, "w")
# a run continued with --from / --only keeps the records of the earlier steps of the same run
if first(todo) != stepnames[1] && isfile(MANIFEST)
    old = TOML.parsefile(MANIFEST)
    merge!(manifest["steps"], get(old, "steps", Dict()))
    for k in ("data_sha256", "time_limit_hits")
        haskey(old, k) && (manifest[k] = old[k])
    end
    manifest["continued_from"] = get(old, "started", "")
end

# ── run ────────────────────────────────────────────────────────────────────────
mkpath(LOGS)
println("Pipeline: $(join(todo, " → ")) ($THREADS threads$(RESUME ? ", resume" : ""))")
for (name, what, dir, script, threaded, env) in STEPS
    name in todo || continue
    log = joinpath(LOGS, "$name.log")
    # hold point: the file outputs/hold_before_<step> pauses the run before that step until it is deleted
    hold = joinpath(PKG, "outputs", "hold_before_$name")
    if isfile(hold)
        println("[$(Dates.format(now(), "HH:MM"))] HOLD before $name: delete outputs/hold_before_$name to continue"); flush(stdout)
        while isfile(hold)
            sleep(20)
        end
        println("[$(Dates.format(now(), "HH:MM"))] hold released, continuing with $name"); flush(stdout)
    end
    cmd = `$JL --startup-file=no $(threaded ? ["-t", string(THREADS)] : String[]) --project=$PKG $script`
    cmd = setenv(Cmd(cmd; dir=dir), merge(Dict(ENV), env))
    t0 = time(); println("\n[$(Dates.format(now(), "HH:MM"))] $name: $what  (log: outputs/logs/$name.log)"); flush(stdout)
    ok = try
        run(pipeline(cmd; stdout=log, stderr=log)); true
    catch
        false
    end
    manifest["steps"][name] = Dict("ok" => ok, "seconds" => round(time() - t0; digits=1), "log" => "outputs/logs/$name.log",
                                   "env" => env, "finished" => string(now()))
    name == "data" && (manifest["data_sha256"] = data_files())
    if name == "fits" && ok      # fits that reached their safety time limit (0 = the run is reproducible)
        f = joinpath(PKG, "outputs", "zoo_fits", "time_limit_hits.csv")
        isfile(f) && (manifest["time_limit_hits"] = countlines(f) - 1)
    end
    save_manifest()
    @info "$name $(ok ? "done" : "FAILED") in $(round(Int, time() - t0)) s"
    ok || (println("Stopped: step $name failed; see $log"); exit(1))
end
haskey(manifest, "data_sha256") || (manifest["data_sha256"] = data_files())
manifest["finished"] = string(now()); save_manifest()
println("\nDone. Run manifest: $MANIFEST")
