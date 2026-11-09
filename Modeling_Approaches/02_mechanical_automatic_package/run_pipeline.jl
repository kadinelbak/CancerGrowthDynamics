# The whole pipeline, data to Final Report, in one command. Each step runs in its own Julia process
# with its log in outputs/logs/<step>.log; a run manifest (outputs/run_manifest.toml) records the
# code version, Julia and package versions, threads, data checksums and every step's time and exit code.
#
#   julia run_pipeline.jl                       everything from the raw data (about 3.5–4 h on 12 threads; less when
#                                               fits and bootstraps are unchanged and reused)
#   julia run_pipeline.jl --threads 12          threads for each fitting/control step (default: all cores − 2)
#   julia run_pipeline.jl --jobs 3              steps run at the same time when they do not depend on each other
#                                               (default 3; --jobs 1 = one step after another, the order before 2026-11-09)
#   julia run_pipeline.jl --control60           also run the 60-day control and its report (skipped by default since
#                                               2026-11-09, user; the Final Report marks that chapter as from an earlier run)
#   julia run_pipeline.jl --from analyses       start at a step (earlier outputs are reused)
#   julia run_pipeline.jl --only report14       one step (see --list)
#   julia run_pipeline.jl --resume              reuse unchanged fits and bootstraps (ZOO_RESUME, ZOO_REUSE_BOOT)
#                                               and the control checkpoints instead of recomputing. Since 2026-11-09 the
#                                               fits (fit_cache/) and bootstraps (boot_checkpoints/) are reused by
#                                               default whenever nothing they depend on changed; ZOO_REFIT=1 / ZOO_REBOOT=1
#                                               recompute them
#   julia run_pipeline.jl --recompute-control   rerun the control steps even when their inputs are unchanged
#                                               (by default a control step whose code, configs, data and
#                                               carried models + bootstrap draws match its last finished run
#                                               is skipped and its outputs reused; CONTROL_FORCE=1)
#   julia run_pipeline.jl --no-sysimage         start the steps without the precompiled image (tools/sysimage), which is
#                                               used when it exists and matches Manifest.toml
#   julia run_pipeline.jl --list                the steps and what each waits for
#
# Steps run as soon as the steps they need have finished, up to --jobs at a time (e.g. the all-doses variant runs
# next to the control steps). Each step's wall time, CPU time and average cores used go into the run manifest, so
# idle cores show up there. The results do not depend on the order: every step reads only the outputs of the steps
# it waits for.
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
const RECOMPUTE_CONTROL = flag("--recompute-control")
const CONTROL60 = flag("--control60")
const NO_SYSIMAGE = flag("--no-sysimage")
const LIST = flag("--list")
const THREADS = parse(Int, opt("--threads", string(max(Sys.CPU_THREADS - 2, 1))))
const FROM = opt("--from", "")
const ONLY = opt("--only", "")
const JOBS = max(parse(Int, opt("--jobs", "3")), 1)
isempty(args) || error("unknown arguments: $(join(args, " "))")

# (name, what it does, working directory, script and arguments, threads?, extra environment, steps it waits for)
redo = RESUME ? Dict{String,String}() : Dict("REDO" => "1")
RECOMPUTE_CONTROL && (redo["CONTROL_FORCE"] = "1")     # the control steps get `redo`
const REPORTS_DONE = ["report14", "report30", "report_long", "sweep_report"]
const STEPS = [
    ("data",       "tidy data table from the raw sources",                  PKG,  ["data/build_tidy.jl"],             false, Dict{String,String}(), String[]),
    ("checks",     "the tidy table reproduces the source files", PKG, ["test/test_data.jl"],             false, Dict{String,String}(), ["data"]),
    ("fits",       "staged model-zoo fits and bootstrap",                   PKG,  ["run_zoo_fits.jl"],                true,  RESUME ? Dict("ZOO_RESUME" => "1", "ZOO_REUSE_BOOT" => "1") : Dict{String,String}(), ["checks"]),
    ("analyses",   "bootstrap comparison, sensitivity, validation",         PKG,  ["zoo_analyses.jl"],                true,  redo, ["fits"]),     # a full run recomputes every analysis (its cache only tracks the carried models)
    ("control14",  "14-day optimal control (dose)",                         PKG,  ["run_control.jl"],                 true,  redo, ["analyses"]),
    ("report14",   "14-day control report",                                 PKG,  ["control_report.jl"],              false, Dict{String,String}(), ["control14"]),
    ("control_long", "long optimal control, dose + fasting (control_config long_days)", PKG, ["run_control_fasting.jl"], true, redo, ["analyses"]),
    ("report_long",  "long control report",                                 PKG,  ["control_fasting_report.jl"],      false, Dict{String,String}(), ["control_long"]),
    # PRD 2026-10-27: 30- and 60-day dose-only runs, the all-doses Stage 2b variant and the control sweep
    ("control30",  "30-day optimal control (dose; beyond the data)",        PKG,  ["run_control.jl"],                 true,  merge(redo, Dict("CONTROL_DAYS" => "30")), ["analyses"]),
    ("report30",   "30-day control report",                                 PKG,  ["control_report.jl"],              false, Dict("CONTROL_DAYS" => "30"), ["control30"]),
    # the 60-day run (27 min on 12 threads) only with --control60 (user, 2026-11-09)
    ("control60",  "60-day optimal control (dose; beyond the data; only with --control60)", PKG, ["run_control.jl"], true, merge(redo, Dict("CONTROL_DAYS" => "60")), ["analyses"]),
    ("report60",   "60-day control report (only with --control60)",         PKG,  ["control_report.jl"],              false, Dict("CONTROL_DAYS" => "60"), ["control60"]),
    # the variant reads the main run's saved fits (same data = same fits), so it waits for the main fits
    ("alldoses",   "Stage fits keeping the A2780cis 1 µM monoculture (variant, point fits)", PKG, ["run_zoo_fits.jl"], true, Dict("ZOO_NO_EXCLUDE" => "1", "ZOO_OUT" => joinpath(PKG, "outputs", "sensitivity", "all_doses"), "BOOT_B" => "0"), ["fits"]),
    ("sweep",      "control sweep: barriers, thresholds, run lengths",       PKG,  ["run_control_sweep.jl"],           true,  redo, ["analyses"]),
    # both sweeps write outputs/control_sweep/sweep_outcomes.csv, so they run one after the other
    ("sweep_all",  "control sweep on the all-doses variant",                PKG,  ["run_control_sweep.jl"],           true,  merge(redo, Dict("ZOO_STATE" => joinpath(PKG, "outputs", "sensitivity", "all_doses", "zoo_state.jls"), "SWEEP_DATA" => "all doses")), ["alldoses", "sweep"]),
    ("sweep_report", "control sweep report",                                PKG,  ["control_sweep_report.jl"],        false, Dict{String,String}(), ["sweep", "sweep_all"]),
    ("findings",   "key findings, assumptions register and limitations",   PKG,  ["build_findings_report.jl"],       false, Dict{String,String}(), vcat(REPORTS_DONE, ["report60"])),
    ("tests",      "control, stage and validation tests on the new results", PKG, ["test/run_tests.jl"],              true,  Dict{String,String}(), ["findings"]),
    ("final",      "Final Report (all chapters)",                           ROOT, ["build_final_report.jl"],          false, Dict{String,String}(), ["tests"]),
]
const OPTIONAL = Set(CONTROL60 ? String[] : ["control60", "report60"])     # left out unless asked for
stepnames = [s[1] for s in STEPS]
if LIST
    foreach(s -> println(rpad(s[1], 13), s[2], isempty(s[7]) ? "" : "  [after $(join(s[7], ", "))]"), STEPS); exit(0)
end
for s in (FROM, ONLY)
    isempty(s) || s in stepnames || error("unknown step $s; steps: $(join(stepnames, ", "))")
end
todo = !isempty(ONLY) ? [ONLY] : isempty(FROM) ? stepnames : stepnames[findfirst(==(FROM), stepnames):end]
isempty(ONLY) && (todo = [t for t in todo if !(t in OPTIONAL)])

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
# when several steps could start, the ones that took longest last time go first (their wall time from the last manifest)
const PREV_SECS = try
    Dict{String,Float64}(k => Float64(get(v, "seconds", 0.0)) for (k, v) in get(TOML.parsefile(MANIFEST), "steps", Dict()))
catch
    Dict{String,Float64}()
end
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
# CPU seconds of a running step (Linux: /proc/<pid>/stat utime + stime + the same for its finished child
# processes, e.g. the test files); nothing elsewhere. Sampled every 2 s while the step runs.
const CLK_TCK = try parse(Int, readchomp(`getconf CLK_TCK`)) catch; 100 end
function cpu_seconds(pid)
    f = "/proc/$pid/stat"; isfile(f) || return nothing
    s = try read(f, String) catch; return nothing end
    fields = split(s[findlast(')', s)+2:end])          # after "pid (comm) "
    sum(parse(Int, fields[i]) for i in 12:15) / CLK_TCK  # utime stime cutime cstime (fields 14–17)
end

# precompiled image (tools/sysimage/build_sysimage.jl): used while the package Manifest it was built from is unchanged
const SYSIMAGE = let d = joinpath(PKG, "tools", "sysimage"), img = joinpath(d, "a2780_sysimage." * (Sys.iswindows() ? "dll" : Sys.isapple() ? "dylib" : "so"))
    sh = joinpath(d, "manifest.sha256")
    if NO_SYSIMAGE || !isfile(img)
        ""
    elseif isfile(sh) && strip(read(sh, String)) == bytes2hex(open(sha256, joinpath(PKG, "Manifest.toml")))
        img
    else
        println("tools/sysimage: the image is older than Manifest.toml, so it is not used; rebuild it with " *
                "julia --project=tools/sysimage tools/sysimage/build_sysimage.jl"); ""
    end
end
manifest["sysimage"] = isempty(SYSIMAGE) ? "none" : relpath(SYSIMAGE, PKG)

mkpath(LOGS)
stamp() = Dates.format(now(), "HH:MM")
println("Pipeline: $(join(todo, ", ")) ($THREADS threads per step, up to $JOBS steps at a time$(RESUME ? ", resume" : ""))")
isempty(OPTIONAL) || isempty(ONLY) && println("Left out (add --control60 to run them): $(join(sort(collect(OPTIONAL)), ", "))")
const SPEC = Dict(s[1] => s for s in STEPS)
done = Set{String}(); failed = String[]
running = Dict{String,Any}()        # name => (proc, io, t0, cpu, started)
held = Set{String}()
ready(name) = all(d -> d in done || !(d in todo), SPEC[name][7])
function launch(name)
    _, what, dir, script, threaded, env, _ = SPEC[name]
    cmd = `$JL --startup-file=no $(isempty(SYSIMAGE) ? String[] : ["-J", SYSIMAGE]) $(threaded ? ["-t", string(THREADS)] : String[]) --project=$PKG $script`
    cmd = setenv(Cmd(cmd; dir=dir), merge(Dict(ENV), env, Dict("PIPELINE_STEP" => name)))
    rm(joinpath(PKG, "outputs", "reuse", "$name.toml"); force=true)      # written by run_control*.jl when the step skipped its work
    io = open(joinpath(LOGS, "$name.log"), "w")       # one handle for both streams
    println("\n[$(stamp())] $name: $what  (log: outputs/logs/$name.log)"); flush(stdout)
    p = run(pipeline(cmd; stdout=io, stderr=io); wait=false)
    running[name] = (proc=p, io=io, t0=time(), cpu=Ref(0.0), started=string(now()))
end
function finish(name)
    r = pop!(running, name); close(r.io)
    ok = success(r.proc); wall = time() - r.t0
    env = SPEC[name][6]
    # the code version of every step, so a run continued with --from / --only records what produced each output
    rec = Dict{String,Any}("ok" => ok, "seconds" => round(wall; digits=1), "log" => "outputs/logs/$name.log",
                           "env" => env, "started" => r.started, "finished" => string(now()), "commit" => git("rev-parse", "HEAD"),
                           "uncommitted_changes" => !isempty(git("status", "--porcelain")),
                           "threads" => SPEC[name][5] ? THREADS : 1)
    if r.cpu[] > 0           # CPU time and the average number of busy cores (Linux)
        rec["cpu_seconds"] = round(r.cpu[]; digits=1); rec["cores_used"] = round(r.cpu[] / max(wall, 1e-9); digits=2)
    end
    reuse_file = joinpath(PKG, "outputs", "reuse", "$name.toml")
    if isfile(reuse_file)      # the step reused its last finished outputs: say so, with the earlier run's time and fingerprint
        ru = try TOML.parsefile(reuse_file) catch; Dict{String,Any}() end
        for k in ("reused_from", "reused_fingerprint")
            haskey(ru, k) && (rec[k] = ru[k])
        end
    end
    # fits and bootstraps reused unchanged inside a fitting step (run_zoo_fits.jl reuse_summary.toml)
    if name in ("fits", "alldoses")
        rs = joinpath(get(env, "ZOO_OUT", joinpath(PKG, "outputs", "zoo_fits")), "reuse_summary.toml")
        if isfile(rs) && mtime(rs) >= r.t0
            for (k, v) in TOML.parsefile(rs)
                rec[k] = v
            end
        end
    end
    manifest["steps"][name] = rec
    name == "data" && (manifest["data_sha256"] = data_files())
    if name == "fits" && ok      # fits that reached their safety time limit (0 = the run is reproducible)
        f = joinpath(PKG, "outputs", "zoo_fits", "time_limit_hits.csv")
        isfile(f) && (manifest["time_limit_hits"] = countlines(f) - 1)
    end
    save_manifest()
    cores = haskey(rec, "cores_used") ? ", $(rec["cores_used"]) cores busy on average" : ""
    println("[$(stamp())] $name $(ok ? "done" : "FAILED") in $(round(Int, wall)) s$cores"); flush(stdout)
    ok ? push!(done, name) : push!(failed, name)
end
pending = copy(todo)
while !isempty(pending) || !isempty(running)
    # start every step whose inputs are finished, up to JOBS at a time (none after a failure)
    if isempty(failed)
        for name in sort(pending; by=n -> -get(PREV_SECS, n, 0.0))
            length(running) >= JOBS && break
            ready(name) || continue
            # hold point: the file outputs/hold_before_<step> pauses that step until it is deleted
            if isfile(joinpath(PKG, "outputs", "hold_before_$name"))
                name in held || (println("[$(stamp())] HOLD before $name: delete outputs/hold_before_$name to continue"); flush(stdout); push!(held, name))
                continue
            end
            name in held && (println("[$(stamp())] hold released, continuing with $name"); flush(stdout))
            filter!(!=(name), pending); launch(name)
        end
    elseif isempty(running)
        break
    end
    for (name, r) in collect(running)
        if process_exited(r.proc)
            finish(name)
        else
            c = cpu_seconds(getpid(r.proc)); c === nothing || (r.cpu[] = max(r.cpu[], c))
        end
    end
    isempty(running) && !isempty(pending) && !any(ready, pending) && isempty(failed) && isempty(held) &&
        error("steps cannot start (their inputs are not in this run): $(join(pending, ", "))")
    sleep(isempty(running) && isempty(held) ? 0.5 : 2)
end
if !isempty(failed)
    println("Stopped: step(s) $(join(failed, ", ")) failed; see outputs/logs/<step>.log")
    exit(1)
end
haskey(manifest, "data_sha256") || (manifest["data_sha256"] = data_files())
manifest["finished"] = string(now()); save_manifest()
println("\nDone. Run manifest: $MANIFEST")
