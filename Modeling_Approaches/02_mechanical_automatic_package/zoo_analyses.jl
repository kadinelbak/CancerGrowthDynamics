# Bootstrapping, sensitivity analysis and validation of the compartment-model (zoo) fits.
# Reads the checkpoint of run_zoo_fits.jl (outputs/zoo_fits/zoo_state.jls); the carried-forward
# (optimal) model of every fit group is the first candidate. Uses the GrowthParameterEstimation
# routines the earlier reports used:
#   bootstrapping   the zoo's 200-replicate moving-block residual bootstrap (every candidate refitted,
#                   selection rule replayed) + GrowthParameterEstimation.bootstrap_joint_fit (within-
#                   trajectory residual bootstrap of the optimal model) for comparison; residual
#                   autocorrelation shows why blocks are used
#   sensitivity     local sensitivity index (the formula of GrowthParameterEstimation.
#                   parameter_sensitivity_analysis: ±10% per parameter, max relative change of the
#                   prediction / 0.1, applied to every observed state), profile likelihood refits
#                   (run_joint_multistart) and the bound check of profile_joint_fit_bounds
#   validation      GrowthParameterEstimation.blocked_joint_validation (leave one dose / ratio / run
#                   out), a time-block forecast (fit the first 70% of days, predict the rest) and the
#                   treated (Ce1) cocultures held out: crowding fitted to Ce0 only, cisplatin from
#                   Stage 2, Ce1 predicted (the same held-out set as the neural-network track)
# Heavy results are checkpointed in outputs/zoo_analyses/analyses_state.jls (REDO=1 recomputes).
#   julia -t 12 --project=. zoo_analyses.jl

include(joinpath(@__DIR__, "zoo_common.jl"))
include(joinpath(@__DIR__, "zoo_models.jl"))
include(joinpath(@__DIR__, "bootstrap_tools.jl"))
using Serialization, Statistics, Random
import LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
const GPE = GrowthParameterEstimation

const ZOO_DIR  = joinpath(@__DIR__, "outputs", "zoo_fits")
const ZOO_DOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "zoo_fits")
const BLUE, ORANGE, GREEN, PURPLE = "#0b6e99", "#d95f02", "#3b8b4f", "#8b5cf6"
state = deserialize(joinpath(ZOO_DIR, "zoo_state.jls"))
const BOOT_B = state.B
include(joinpath(@__DIR__, "zoo_selection.jl"))
JOINT_INFO = hasproperty(state, :joint) ? state.joint : joint_info_from_groups(state.groups)
apply_joint_replay!(state.groups, state.results)
include(joinpath(@__DIR__, "export_zoo_bundle.jl"))   # zoo_series, zoo_panels, zoo_inherited
include(joinpath(@__DIR__, "zoo_report.jl"))          # zoo_band_svg, hesc
const groups = state.groups
const results = state.results
gid(id) = groups[findfirst(g -> g.id == id, groups)]

const AN_DIR  = joinpath(@__DIR__, "outputs", "zoo_analyses")
const AN_DOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "zoo_analyses")
mkpath(AN_DIR)
# SMOKE=1: tiny settings and a separate checkpoint, to test the script end to end in minutes
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const NBOOT = SMOKE ? 2 : BOOT_B
const CKPT = joinpath(AN_DIR, SMOKE ? "analyses_smoke.jls" : "analyses_state.jls")
ck = (isfile(CKPT) && get(ENV, "REDO", "0") != "1") ? deserialize(CKPT) : Dict{Symbol,Any}()
save_ck() = serialize(CKPT, ck)
# Groups whose carried model is new or has changed since the checkpoint are recomputed; the rest reuse it.
# A model is identified by its name and fitted values, and a group is also stale when a stage it inherits
# from is (e.g. Stage 4 keeps its name when Stage 3 is refitted underneath it).
model_key(c) = c.name * " | " * join(string.(c.params), ",")
current_models = Dict(g.id => model_key(g.cands[1]) for g in groups)
prev_models = get(ck, :models, haskey(ck, :pkgboot) ? Dict(id => current_models[id] for id in keys(ck[:pkgboot]) if haskey(current_models, id)) : Dict{String,String}())
# older checkpoints stored the name only: compare names for those
same_model(id, k) = (p = get(prev_models, id, ""); occursin(" | ", p) ? p == k : p == first(split(k, " | ")))
stale = Set(id for (id, k) in current_models if !same_model(id, k))
const INHERITS = Dict("s1b-a2780" => ["s1-a2780"], "s1b-a2780cis" => ["s1-a2780cis"], "s2-a2780" => ["s1-a2780"], "s2-a2780cis" => ["s1-a2780cis"],
                      "s3" => ["s2-a2780", "s2-a2780cis"], "s3b" => ["s1b-a2780", "s1b-a2780cis"], "s4" => ["s3", "s3b"])
for _ in 1:4, (id, srcs) in INHERITS
    haskey(current_models, id) && any(in(stale), srcs) && push!(stale, id)
end
if !isempty(stale) && !isempty(ck)
    haskey(ck, :pkgboot) && (ck[:pkgboot] = Dict(k => v for (k, v) in ck[:pkgboot] if !(k in stale)))
    haskey(ck, :profile) && (ck[:profile] = [q for q in ck[:profile] if !(q.gid in stale)])
    haskey(ck, :validation) && (ck[:validation] = merge(ck[:validation], (rows=[r for r in ck[:validation].rows if !(r.group in stale)],)))
    haskey(ck, :gsa) && (ck[:gsa] = merge(ck[:gsa], (stages=Dict(k => v for (k, v) in ck[:gsa].stages if !(k in stale)),)))
end
any(id -> startswith(id, "s1-") || startswith(id, "s2-"), stale) && delete!(ck, :intermittent)   # predicted from Stages 1 and 2b
ck[:models] = current_models
println("=== zoo analyses: $(Threads.nthreads()) threads; checkpoint has $(collect(keys(ck))); to (re)compute: $(sort(collect(stale))) ==="); flush(stdout)

# ── shared helpers ─────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "zoo_heldout.jl"))   # FIT_KW, evaluate, gfit, nobs, VAL and the held-out Ce1 test
sub_ds(d, keep) = (x=d.x[keep], y=d.y[keep], state_index=d.state_index)
rmse_of(ev) = ev.rmse
nrmse_of(ev, ds) = ev.rmse / mean(reduce(vcat, [d.y for d in ds]))
# value q with the i-th entry v inserted (keeps ForwardDiff / Dual element types)
function insert_at(q, i, v)
    T = promote_type(eltype(q), typeof(v)); r = Vector{T}(undef, length(q) + 1)
    for j in eachindex(r)
        r[j] = j < i ? q[j] : j == i ? T(v) : q[j-1]
    end
    r
end
drop(v, i) = [v[j] for j in eachindex(v) if j != i]
pl(n) = "\\($(ptex(n))\\)"
f3(x) = isfinite(x) ? fmt(x) : "—"

# ── 1. bootstrapping ───────────────────────────────────────────────────────────

todo_pk = [g for g in groups if !haskey(get(ck, :pkgboot, Dict()), g.id)]
if !isempty(todo_pk)
    println("\n--- package bootstrap_joint_fit (within-trajectory residuals): $(join([g.id for g in todo_pk], ", ")) ---"); flush(stdout)
    out = Vector{Any}(undef, length(todo_pk))
    Threads.@threads :dynamic for i in eachindex(todo_pk)
        g = todo_pk[i]; c = g.cands[1]; t0 = time()
        out[i] = try
            ev = evaluate(c.model, g.datasets, c.u0, c.params)
            fitted = (params=c.params, predictions=ev.predictions)
            b = GPE.bootstrap_joint_fit(c.model, g.datasets, c.u0, fitted; bounds=c.bounds, n_bootstrap=NBOOT, starts_per_bootstrap=1,
                                        seed=20260929 + i, optimizer=:nelder_mead, maxiters=600, maxtime=30.0, reltol=1e-7, abstol=1e-7, initial_time=0.0)
            (summary=b.parameter_summary, ok=b.successful_replicates, method=b.method)
        catch err
            @warn "package bootstrap failed" g.id err
            nothing
        end
        @printf("  %-12s package bootstrap done in %.0fs\n", g.id, time() - t0); flush(stdout)
    end
    ck[:pkgboot] = merge(get(ck, :pkgboot, Dict{String,Any}()), Dict(todo_pk[i].id => out[i] for i in eachindex(todo_pk))); save_ck()
end

# lag-1 autocorrelation of the optimal model's residuals, per trajectory
function residual_acf(g)
    c = g.cands[1]; ev = evaluate(c.model, g.datasets, c.u0, c.params)
    acs = Float64[]
    for (d, p) in zip(g.datasets, ev.predictions)
        e = d.y .- p; length(e) < 4 && continue
        e = e .- mean(e); den = sum(abs2, e); den > 0 || continue
        push!(acs, sum(e[1:end-1] .* e[2:end]) / den)
    end
    acs
end

# ── 2. sensitivity ─────────────────────────────────────────────────────────────

# Local sensitivity index as in GrowthParameterEstimation.parameter_sensitivity_analysis, over every
# observed state: SI = max relative change of the prediction under a ±10% change / 0.1. Relative
# changes are taken against max(prediction, 1 cell) so trajectories that die out do not divide by 0.
function local_si(g; pert=0.1)
    c = g.cands[1]; tmax = maximum(maximum(d.x) for d in g.datasets)
    idx = [d.state_index for d in g.datasets]
    pred(p) = (s = solve(ODEProblem(c.model, c.u0, (0.0, tmax), p), Tsit5(); saveat=0.1, reltol=1e-8, abstol=1e-8); [u[j] for u in s.u, j in idx])
    base = pred(c.params)
    rows = NamedTuple[]
    for i in eachindex(c.params)
        maxrel = 0.0; meanabs = 0.0; cnt = 0; ok = true
        for s in (1 + pert, 1 - pert)
            p = copy(c.params); p[i] *= s
            P = try pred(p) catch; ok = false; break end
            size(P) == size(base) || (ok = false; break)
            a = abs.(P .- base)
            maxrel = max(maxrel, maximum(a ./ max.(abs.(base), 1.0))); meanabs += sum(a); cnt += length(a)
        end
        push!(rows, (group=g.id, param=c.pnames[i], value=c.params[i], si=ok ? maxrel / pert : NaN,
                     max_rel_change=ok ? maxrel : NaN, mean_abs_change=ok ? meanabs / cnt : NaN))
    end
    rows
end

# Profile likelihood: fix one parameter on a grid (0.5×–2× the estimate), refit the others with
# run_joint_multistart, and record n·ln(SSE/SSE_min); a 95% interval is where this stays ≤ 3.84 (χ²₁).
const PROFILE_MULT = [0.5, 0.63, 0.79, 0.9, 1.0, 1.11, 1.26, 1.59, 2.0]
function profile_jobs(gs)
    jobs = NamedTuple[]
    for g in gs, (i, n) in enumerate(g.cands[1].pnames)
        c = g.cands[1]; lo, hi = c.bounds[i]; est = c.params[i]
        vals = abs(est) > 1e-6 ? est .* PROFILE_MULT : collect(range(lo, lo + 0.1 * (hi - lo); length=9))
        vals = unique(clamp.(vals, lo, hi))
        for v in vals
            push!(jobs, (gid=g.id, i=i, name=n, v=v))
        end
    end
    jobs
end
done_prof = Set(q.gid for q in get(ck, :profile, NamedTuple[]))
todo_prof = [g for g in groups if !(g.id in done_prof) && !isempty(g.cands[1].pnames)]
if !isempty(todo_prof)
    jobs = profile_jobs(todo_prof)
    SMOKE && (jobs = [j for (n, j) in enumerate(jobs) if n % 25 == 1])
    println("\n--- profile likelihood: $(length(jobs)) refits ---"); flush(stdout)
    out = Vector{Any}(undef, length(jobs)); done = Threads.Atomic{Int}(0)
    Threads.@threads :dynamic for j in eachindex(jobs)
        J = jobs[j]; g = gid(J.gid); c = g.cands[1]; k = length(c.params)
        sse = try
            if k == 1
                evaluate(c.model, g.datasets, c.u0, [J.v]).raw_sse
            else
                m = (du, u, q, t) -> c.model(du, u, insert_at(q, J.i, J.v), t)
                gfit(m, g.datasets, c.u0, drop(c.params, J.i), drop(c.bounds, J.i); starts=2, seed=j).sse
            end
        catch err
            NaN
        end
        out[j] = merge(J, (sse=sse,))
        n = Threads.atomic_add!(done, 1) + 1
        n % 20 == 0 && (println("    $n / $(length(jobs)) profile points"); flush(stdout))
    end
    ck[:profile] = vcat(get(ck, :profile, NamedTuple[]), out); save_ck()
end

function profile_summary(g)
    c = g.cands[1]; n = nobs(g.datasets)
    rows = NamedTuple[]
    for (i, name) in enumerate(c.pnames)
        pts = sort([p for p in ck[:profile] if p.gid == g.id && p.i == i && isfinite(p.sse)]; by=p -> p.v)
        isempty(pts) && continue
        smin = min(c.sse, minimum(p.sse for p in pts))
        prof = [(v=p.v, d=n * log(p.sse / smin)) for p in pts]
        est = c.params[i]
        lower = [p for p in prof if p.v < est]; upper = [p for p in prof if p.v > est]
        lo_ok = any(p -> p.d > 3.84, lower); hi_ok = any(p -> p.d > 3.84, upper)
        inside = [p.v for p in prof if p.d <= 3.84]
        status = lo_ok && hi_ok ? "identifiable" : (lo_ok || hi_ok) ? "one-sided" : "flat (not identifiable within 0.5–2×)"
        lo, hi = c.bounds[i]; span = hi - lo
        pos = span > 0 ? (est - lo) / span : NaN
        push!(rows, (group=g.id, param=name, estimate=est, ci_lo=isempty(inside) ? NaN : minimum(inside), ci_hi=isempty(inside) ? NaN : maximum(inside),
                     reached_lo=lo_ok, reached_hi=hi_ok, status=status, bound_position=pos,
                     near_bound=(pos >= 0.98 || est - lo <= 0.02 * abs(est) + 1e-8), profile=prof))
    end
    rows
end

# ── 2b. global sensitivity (Sobol indices) ────────────────────────────────────
# Variance-based global sensitivity: every parameter is drawn uniformly over its 95% block-
# bootstrap interval (its plausible range given the data), the model is solved for each draw, and
# the variance of an output is split between the parameters. First-order index S1 = share
# explained by that parameter alone; total index ST = share including its interactions.
# Saltelli sampling with the Saltelli (2010) S1 and Jansen ST estimators; 95% intervals from
# resampling the draws. Not in GrowthParameterEstimation, so implemented here.
using Random: Xoshiro
const GSA_N = SMOKE ? 64 : 16384

function sobol(f, lo, hi, m; N=GSA_N, seed=7, nboot=200)
    k = length(lo); rng = Xoshiro(seed)
    draw() = [lo[j] + rand(rng) * (hi[j] - lo[j]) for _ in 1:N, j in 1:k]
    A = draw(); B = draw()
    mats = vcat([A, B], [(M = copy(A); M[:, i] = B[:, i]; M) for i in 1:k])
    Y = [fill(NaN, N, m) for _ in mats]
    jobs = [(a, n) for a in eachindex(mats) for n in 1:N]
    Threads.@threads :dynamic for (a, n) in jobs
        y = try f(mats[a][n, :]) catch; fill(NaN, m) end
        Y[a][n, :] = y
    end
    YA, YB = Y[1], Y[2]
    ok = [all(isfinite, YA[n, :]) && all(isfinite, YB[n, :]) && all(i -> all(isfinite, Y[2+i][n, :]), 1:k) for n in 1:N]
    rows = NamedTuple[]
    brng = Xoshiro(seed + 1)
    for o in 1:m
        idx = findall(ok)
        est(ix) = begin
            ya, yb = YA[ix, o], YB[ix, o]; V = var(vcat(ya, yb)); yb = yb .- mean(yb)   # centring lowers the S1 estimator variance
            V > 0 || return (fill(NaN, k), fill(NaN, k))
            ([mean(yb .* (Y[2+i][ix, o] .- ya)) / V for i in 1:k], [0.5 * mean((ya .- Y[2+i][ix, o]) .^ 2) / V for i in 1:k])
        end
        S1, ST = est(idx)
        bs = [est(idx[rand(brng, 1:length(idx), length(idx))]) for _ in 1:nboot]
        for i in 1:k
            s1b = filter(isfinite, [b[1][i] for b in bs]); stb = filter(isfinite, [b[2][i] for b in bs])
            push!(rows, (output=o, param=i, S1=S1[i], S1_lo=isempty(s1b) ? NaN : quantile(s1b, 0.025), S1_hi=isempty(s1b) ? NaN : quantile(s1b, 0.975),
                         ST=ST[i], ST_lo=isempty(stb) ? NaN : quantile(stb, 0.025), ST_hi=isempty(stb) ? NaN : quantile(stb, 0.975),
                         n_ok=length(idx), mean_out=mean(vcat(YA[idx, o], YB[idx, o])), sd_out=std(vcat(YA[idx, o], YB[idx, o]))))
        end
    end
    rows
end

# 95% block-bootstrap interval of free parameter i of group g (widened slightly if degenerate)
function boot_range(g, i)
    l, _, h = qs(results[g.id].res.pars[1][:, i])
    (isfinite(l) && isfinite(h) && h > l) ? (l, h) : (0.99 * g.cands[1].params[i], 1.01 * g.cands[1].params[i])
end
last_day(g, k) = maximum(g.datasets[k].x)

# outputs per fit: (label, dataset index)
function gsa_outputs(g)
    g.kind in (:s1, :s1b) && return [("count on the last day", 1)]
    g.kind == :s2 && return [("count on the last day at 1.0 µM", findfirst(==(1.0), treated[g.line].doses))]
    g.kind == :s3 && return [("Ce1 50:50 sensitive, last day", 9), ("Ce1 50:50 resistant, last day", 10)]
    g.kind == :s4 && return [("sensitive, last day", 1), ("resistant, last day", 2)]
    [("run 1, 1:1 sensitive, last day", 3), ("run 1, 1:1 resistant, last day", 4)]
end
function gsa_stage(g)
    c = g.cands[1]; outs = gsa_outputs(g)
    T = maximum(last_day(g, k) for (_, k) in outs)
    lo = [boot_range(g, i)[1] for i in eachindex(c.params)]; hi = [boot_range(g, i)[2] for i in eachindex(c.params)]
    f(p) = (s = solve(ODEProblem(c.model, c.u0, (0.0, T), p), Tsit5(); saveat=[T], reltol=1e-6, abstol=1e-6);
            [s.u[end][g.datasets[k].state_index] for (_, k) in outs])
    (outputs=[o[1] for o in outs], names=c.pnames, lo=lo, hi=hi, rows=sobol(f, lo, hi, length(outs)))
end

# The whole treated-coculture model (Stage 3, Ce1 at 1 µM): every parameter free, including the
# growth and drug parameters carried over from Stages 1 and 2, each over its own bootstrap interval.
function gsa_full()
    LN, LC = c3.sys.lines
    LNf = merge(LN, (s=merge(LN.s, (r=Slot(6, 0.0), h=Slot(7, 0.0), k=Slot(8, 0.0))),))
    LCf = merge(LC, (s=merge(LC.s, (r=Slot(9, 0.0), θ=Slot(10, 0.0), h=Slot(11, 0.0), k=Slot(12, 0.0))),))
    sys = build_system([LNf, LCf], c3.sys.envs, c3.sys.alpha); m = make_system_model(sys)
    s1N, s1C, s2N, s2C = gid("s1-a2780"), gid("s1-a2780cis"), gid("s2-a2780"), gid("s2-a2780cis")
    src = vcat([(n, g3, i) for (i, n) in enumerate(c3.pnames)],
               [("r_N", s1N, 1), ("h_N", s2N, 2), ("k_N", s2N, 3), ("r_C", s1C, 1), ("θ_C", s1C, 3), ("h_C", s2C, 2), ("k_C", s2C, 3)])
    @assert [s[2].cands[1].pnames[s[3]] for s in src] == [s[1] for s in src]
    lo = [boot_range(g, i)[1] for (_, g, i) in src]; hi = [boot_range(g, i)[2] for (_, g, i) in src]
    iN, iC = g3.datasets[9].state_index, g3.datasets[10].state_index      # Ce1 50:50
    T = XMAX
    f(p) = (s = solve(ODEProblem(m, c3.u0, (0.0, T), p), Tsit5(); saveat=[T], reltol=1e-6, abstol=1e-6);
            (N = s.u[end][iN]; C = s.u[end][iC]; [N, C, N + C, C / (N + C)]))
    (outputs=["sensitive (A2780), day 14", "resistant (A2780cis), day 14", "total, day 14", "resistant fraction, day 14"],
     names=[s[1] for s in src], carried=[s[2].id for s in src], lo=lo, hi=hi, rows=sobol(f, lo, hi, 4))
end


# ── 3. validation ──────────────────────────────────────────────────────────────

# blocks: one label per trajectory; uses GrowthParameterEstimation.blocked_joint_validation
function blocked(g, labels, scheme)
    c = g.cands[1]
    df = GPE.blocked_joint_validation(c.model, g.datasets, c.u0, c.params, labels; bounds=c.bounds, starts_per_fold=3, seed=42,
                                      optimizer=:nelder_mead, maxiters=2000, maxtime=180.0, reltol=1e-8, abstol=1e-8, initial_time=0.0)
    rows = NamedTuple[]
    for r in eachrow(df)
        held = [k for k in eachindex(labels) if string(labels[k]) == r.block]
        insample = evaluate(c.model, g.datasets[held], c.u0, c.params)
        m = mean(reduce(vcat, [g.datasets[k].y for k in held]))
        push!(rows, (group=g.id, scheme=scheme, block=r.block, status=r.status, n_test=r.n_test_points,
                     rmse_heldout=r.rmse, rmse_insample=insample.rmse, nrmse_heldout=r.rmse / m, nrmse_insample=insample.rmse / m))
    end
    rows
end

# time-block forecast: fit days ≤ 70% of the last day, predict the later days

function forecast(g)
    c = g.cands[1]; tmax = maximum(maximum(d.x) for d in g.datasets); split = round(VAL["forecast"]["fit_fraction"] * tmax)
    train = [sub_ds(d, d.x .<= split) for d in g.datasets]; test = [sub_ds(d, d.x .> split) for d in g.datasets]
    test = [d for d in test if !isempty(d.x)]
    f = gfit(c.model, train, c.u0, c.params, c.bounds; starts=3)
    ho = evaluate(c.model, test, c.u0, f.params); ins = evaluate(c.model, test, c.u0, c.params)
    m = mean(reduce(vcat, [d.y for d in test]))
    (group=g.id, scheme="forecast", block="fit days ≤ $(Int(split)), predict days > $(Int(split))", status="completed", n_test=ho.n_points,
     rmse_heldout=ho.rmse, rmse_insample=ins.rmse, nrmse_heldout=ho.rmse / m, nrmse_insample=ins.rmse / m)
end

# blocked validations (validation.toml [[leave_out]]): (group, scheme, held-out label of each dataset)
tag_label(s, tag) = tag == "schedule" ? "$(constant_dose(s["schedule"])) µM" : tag == "run" ? "run $(s["run"])" : string(s[tag])
const BLOCKED_SPECS = [(lo["group"], lo["scheme"], () -> [tag_label(s, lo["tag"]) for s in dataset_tags(lo["group"])]) for lo in VAL["leave_out"]]
if !haskey(ck, :validation)
    println("\n--- validation ---"); flush(stdout)
    tasks = Any[]
    for g in groups
        push!(tasks, () -> [forecast(g)])
    end
    append!(tasks, [() -> blocked(gid(id), labs(), scheme) for (id, scheme, labs) in BLOCKED_SPECS])
    SMOKE && (tasks = tasks[[1, 9]])
    out = Vector{Any}(undef, length(tasks))
    Threads.@threads :dynamic for j in eachindex(tasks)
        out[j] = try tasks[j]() catch err; @warn "validation task failed" j err; NamedTuple[] end
        println("    validation task $j done"); flush(stdout)
    end
    rows = reduce(vcat, out)
    ck[:validation] = merge((rows=rows,), heldout_ce1()); save_ck()
elseif "s3" in stale
    println("\n--- held-out Ce1 for the new Stage 3 model ---"); flush(stdout)
    ck[:validation] = merge(ck[:validation], heldout_ce1()); save_ck()
end

# global sensitivity runs here because it needs the Stage 3 model (c3, g3) defined above
if haskey(ck, :gsa) && get(ck[:gsa], :N, 0) == GSA_N && any(g -> !haskey(ck[:gsa].stages, g.id), groups)
    miss = [g for g in groups if !haskey(ck[:gsa].stages, g.id)]
    println("\n--- global sensitivity for new groups: $(join([g.id for g in miss], ", ")) ---"); flush(stdout)
    st = copy(ck[:gsa].stages)
    for g in miss
        st[g.id] = gsa_stage(g)
    end
    ck[:gsa] = merge(ck[:gsa], (stages=st,)); save_ck()
end
if !haskey(ck, :gsa) || get(ck[:gsa], :N, 0) != GSA_N
    println("\n--- global sensitivity (Sobol, N = $GSA_N) ---"); flush(stdout)
    t0 = time()
    full = gsa_full()
    stages = Dict(g.id => gsa_stage(g) for g in groups)
    ck[:gsa] = (full=full, stages=stages, N=GSA_N); save_ck()
    @printf("  done in %.0fs\n", time() - t0); flush(stdout)
end


# forecasts for groups added since the checkpoint (e.g. Stage 4)
if haskey(ck, :validation)
    have = Set(r.group for r in ck[:validation].rows if r.scheme == "forecast")
    miss = [g for g in groups if !(g.id in have)]
    if !isempty(miss)
        println("\n--- forecasts for new groups: $(join([g.id for g in miss], ", ")) ---"); flush(stdout)
        ck[:validation] = merge(ck[:validation], (rows=vcat(ck[:validation].rows, [forecast(g) for g in miss]),)); save_ck()
    end
    # blocked validations dropped because their group was refitted
    haveb = Set((r.group, r.scheme) for r in ck[:validation].rows)
    missb = [s for s in BLOCKED_SPECS if !((s[1], s[2]) in haveb)]
    if !isempty(missb)
        println("\n--- blocked validation for refitted groups: $(join(["$(s[1]) $(s[2])" for s in missb], "; ")) ---"); flush(stdout)
        newrows = Vector{Any}(undef, length(missb))
        Threads.@threads :dynamic for j in eachindex(missb)
            id, scheme, labs = missb[j]
            newrows[j] = try blocked(gid(id), labs(), scheme) catch err; @warn "blocked validation failed" id scheme err; NamedTuple[] end
        end
        gorder = Dict(g.id => i for (i, g) in enumerate(groups))
        allrows = sort(vcat(ck[:validation].rows, reduce(vcat, newrows)); by=r -> (r.scheme == "forecast" ? 0 : 1, get(gorder, r.group, 99)))
        ck[:validation] = merge(ck[:validation], (rows=allrows,)); save_ck()
    end
end

# ── intermittent dosing: 1 µM on days 0–7 then off, or off then on ─────────────────────────────
# Predicted with no refit from each line's carried Stage 2b model (its growth is the carried
# Stage 1 law), driven by the schedule D(t) of the data (data/schedules.toml) through the model
# engine itself. Values are image areas normalised to cell size, so each trajectory has one scale,
# profiled out. A second question: does the drug effect decay after removal at the build-up rate k
# (the model's assumption)? A separate washout rate k_off of the first effect compartment while
# D = 0 is tested and kept only if the 3-BIC rule supports it (chain and compartment kernels only).
const SP = VAL["schedule_prediction"]
const INT_SCHED = Dict(Symbol(k) => v for (k, v) in SP["schedules"])
freeze_all(L, p) = merge(L, (s=NamedTuple{keys(L.s)}(Tuple((q = getfield(L.s, k); q.i > 0 ? Slot(0, p[q.i]) : q) for k in keys(L.s))),))
function int_system(line, sc)
    c = gid("$(SP["model_from"])-$(lowercase(line))").cands[1]
    L = freeze_all(c.sys.lines[1], c.params)
    sys = build_system([L], [(D=dose_input(INT_SCHED[sc]), R=0.0, lines=[1])], [[Slot(0, 1.0)]])
    sys, L
end
# the washout variant: the first effect compartment clears at koff instead of k while the dose is 0
function int_model(line, sc, koff)
    sys, L = int_system(line, sc); f! = make_system_model(sys)
    koff === nothing && return sys, f!
    (L.kernel in (:comp1, :chain2, :chain3) && ncell(L.action) == 0) || return sys, nothing
    i1 = sys.hb[1][1] + 1; k = L.s.k.v
    g!(du, u, p, t) = (f!(du, u, p, t); input_at(sys.envs[1].D, t) == 0 && (du[i1] += (k - koff) * u[i1]); nothing)
    sys, g!
end
# counts at days x from day 1 (= x0); days 0–1 run first from one cell so the effect states hold
# one day of the schedule, as in the experiment
function int_predict(line, sc, koff, x, x0)
    sys, m = int_model(line, sc, koff)
    m === nothing && return fill(NaN, length(x))
    u = sys_u0(sys, [[1.0]])
    pre = solve(ODEProblem(m, u, (0.0, 1.0), Float64[]), Tsit5(); reltol=1e-9, abstol=1e-9).u[end]
    pre[sys.xi[1][1]+1] = x0
    sol = solve(ODEProblem(m, pre, (1.0, maximum(x)), Float64[]), Tsit5(); saveat=x, tstops=sys_breakpoints(sys), reltol=1e-9, abstol=1e-9)
    [u[sys.xi[1][1]+1] for u in sol.u]
end
function int_fit_scale(line, sc, koff, x, y)
    f(ls) = (sc_ = exp(ls); sum(abs2, y .- sc_ .* int_predict(line, sc, koff, x, y[1] / sc_)))
    grid = range(log(0.05), log(20.0); length=40); v = [f(g) for g in grid]; i = argmin(v)
    a, b = grid[max(i - 1, 1)], grid[min(i + 1, end)]
    for _ in 1:40
        c = b - (b - a) / 1.618; d = a + (b - a) / 1.618
        f(c) < f(d) ? (b = d) : (a = c)
    end
    s = exp((a + b) / 2)
    (scale=s, sse=f((a + b) / 2), pred=s .* int_predict(line, sc, koff, x, y[1] / s))
end
if !haskey(ck, :intermittent)
    println("\n--- intermittent dosing ---"); flush(stdout)
    res = Dict{String,Any}()
    for line in ("A2780", "A2780cis")
        xy = Dict(sc => (d = intermittent[line][INT_SCHED[sc]]; (d.x, d.y, Float64.(coalesce.(d.se, NaN)))) for sc in (:onoff, :offon))
        n = sum(length(v[1]) for v in values(xy))
        tot(koff) = sum(int_fit_scale(line, sc, koff, xy[sc][1], xy[sc][2]).sse for sc in (:onoff, :offon))
        s_eq = tot(nothing)
        bic(sse, k) = n * log(sse / n) + k * log(n)
        grid = exp.(range(log(0.01), log(20.0); length=50)); vals = [tot(ko) for ko in grid]
        ok = SP["washout_test"] && all(isfinite, vals); ib = ok ? argmin(vals) : 1
        fits = Dict(sc => int_fit_scale(line, sc, nothing, xy[sc][1], xy[sc][2]) for sc in (:onoff, :offon))
        res[line] = (xy=xy, n=n, sse_equal=s_eq, bic_equal=bic(s_eq, 2), koff_best=ok ? grid[ib] : NaN, sse_koff=ok ? vals[ib] : NaN,
                     bic_koff=ok ? bic(vals[ib], 3) : NaN, profile=ok ? collect(zip(grid, vals)) : Tuple{Float64,Float64}[], fits=fits,
                     model=gid("$(SP["model_from"])-$(lowercase(line))").cands[1].name)
        @printf("  %-9s equal rates BIC %.2f; best k_off %.3g BIC %.2f (Δ %+.2f)\n", line, bic(s_eq, 2), res[line].koff_best, res[line].bic_koff,
                res[line].bic_koff - bic(s_eq, 2))
    end
    ck[:intermittent] = res; save_ck()
end

# ── HTML ───────────────────────────────────────────────────────────────────────
# Every table column and chart carries an "i" button (info_tips.jl, loaded by zoo_report.jl)
# that explains how to read it.

const CSS = """body{font:14px/1.5 system-ui,sans-serif;margin:0;background:#f4f7f8;color:#17212b}main{max-width:1250px;margin:auto;padding:24px}h1{color:#124b6e;margin:6px 0}h2{color:#124b6e;margin-top:30px}h3{margin:14px 0 6px}
.topnav{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:0 0 16px}.nav{display:inline-flex;align-items:center;min-height:34px;padding:0 11px;border-radius:8px;color:#176b87;background:#eef5f6;text-decoration:none;font:650 14px system-ui,sans-serif;white-space:nowrap}.nav-home{padding:6px 9px}.nav-active{color:#fff;background:#176b87}
.intro,.panel{background:#fff;border:1px solid #d9e2e6;border-radius:12px;padding:16px 18px;margin:14px 0;min-width:0;overflow-wrap:anywhere}.intro{border-left:6px solid #176b87}.note{background:#fff6d9;border-radius:8px;padding:10px 14px}
svg{display:block;width:100%;height:auto;margin:10px 0}.tscroll{overflow-x:auto}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:12px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}td,th{border:1px solid #d3dde2;padding:6px 8px;text-align:left;vertical-align:top}th{background:#e9f3f8}td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}tr.win td{background:#eef8f1;font-weight:600}.muted{color:#5b6b73}.flag{color:#b45309;font-weight:600}.ok{color:#2f7d4a;font-weight:600}
.tipwrap{position:relative}.tipwrap>.ii{position:absolute;top:8px;right:8px;z-index:2}
mjx-container:not([display]){display:inline-block;max-width:100%;overflow-x:auto;vertical-align:middle}"""

function page(file, title, active, body)
    home = "<svg width='22' height='22' viewBox='0 0 24 24' aria-hidden='true'><path fill='currentColor' d='M12 3.1 2.8 10.6v10.3h6.1v-6.2h6.2v6.2h6.1V10.6L12 3.1Zm7.1 15.8h-2v-6.2H6.9v6.2h-2V11.6L12 5.8l7.1 5.8v7.3Z'/></svg>"
    nl(label, href, on) = "<a href='$href' class='nav$(on ? " nav-active" : "")'>$label</a>"
    D = "../../../../../"
    nav = "<nav class='topnav' aria-label='Report navigation'><a href='$(D)index.html' class='nav nav-home' aria-label='Reports home' title='Reports home'>$home</a>" *
          nl("Summary", "$(D)outputs/results/summary.html", false) * nl("Model zoo fits", "../../zoo_fits/report.html", false) *
          nl("Bootstrapping", "../bootstrap/report.html", active == "bootstrap") * nl("Sensitivity", "../sensitivity/report.html", active == "sensitivity") *
          nl("Validation", "../validation/report.html", active == "validation") * nl("Final report", "$(D)outputs/final_report/report.html", false) * "</nav>"
    html = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>$title</title>
<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}};</script><script defer src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>
<style>$CSS</style>$INFO_TIPS_HEAD</head><body><main>$nav<h1>$title</h1>
<p class='muted'>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")) from the model zoo run of $(Dates.format(unix2datetime(mtime(joinpath(ZOO_DIR, "zoo_state.jls"))), "yyyy-mm-dd HH:MM")) by <code>zoo_analyses.jl</code>. Every model here is a compartment model with no time-based terms; the model analysed in each fit is the carried-forward optimal model (fewest free parameters within 3 BIC of the lowest). Buttons marked $(tip("Hover over, tap or focus an i button to see how to read the graph or column next to it. Click elsewhere or press Esc to close.", "Every tip ends with what is better: lower, higher, or closer to a target value.")) explain how to read the graph or column next to them.</p>
$body</main></body></html>"""
    dir = joinpath(AN_DIR, file); mkpath(dir); write(joinpath(dir, "report.html"), stage_names(html))
    ddir = joinpath(AN_DOCS, file); mkpath(ddir); cp(joinpath(dir, "report.html"), joinpath(ddir, "report.html"); force=true)
    println("Wrote $(joinpath(ddir, "report.html"))")
end
# head entries: "label" or ("label", "tip text")
# what is better for each table column (info_tips.jl: header_better); unknown columns stop the build
const BETTER = Dict{String,Any}(
    "Parameter" => :none,
    "Estimate" => :none,
    "Fit" => :none,
    "Optimal model" => :none,
    "Output" => :none,
    "Mean" => :none,
    "Range varied" => :none,
    "Fitted in" => :none,
    "Line" => :none,
    "Prediction" => :none,
    "Crowding fitted to Ce0" => :none,
    "Scheme" => :none,
    "Held-out block" => :none,
    "Best k_off (/day)" => :none,
    "Block median" => "Closer to the estimate is better; a large gap means a skewed or unstable fit.",
    "95% block bootstrap" => "Narrower is better: the parameter is more precisely determined.",
    "95% within-trajectory bootstrap (package)" => "Narrower means more precise, but this interval is over-confident; trust the block interval.",
    "Chosen by the rule" => :higher,
    "Lowest BIC" => :higher,
    "Sensitivity index" => "Neither: large means the prediction depends strongly on this parameter (so it must be known precisely); near 0 means it barely matters.",
    "Largest change of the prediction" => "Neither: larger means the prediction is more sensitive to this parameter.",
    "Profile 95% range (within 0.5–2×)" => "Narrower is better: the parameter is more precisely determined.",
    "Profile" => "Better: Identifiable (pinned on both sides); One-sided is weaker; Flat is worst.",
    "Bound check" => "Better: not at a bound.",
    "SD" => "Lower means the prediction is less uncertain across the plausible parameter ranges.",
    "First-order S1 (95%)" => "Neither: a large value marks a parameter that drives the prediction's uncertainty.",
    "Total ST (95%)" => "Neither: near 0 means the parameter could be fixed; large means it matters, alone or through interactions.",
    "RMSE (µm²)" => :lower,
    "RMSE (cells)" => :lower,
    "nRMSE" => :lower,
    "BIC, equal rates" => :lower,
    "BIC with k_off" => :lower,
    "ΔBIC" => "More negative favours adding the washout rate; above −3 the simpler equal-rate model is enough.",
    "Points" => "Higher gives a more reliable score.",
    "Held-out RMSE" => :lower,
    "In-sample RMSE" => :lower,
    "Held-out nRMSE" => :lower,
    "In-sample nRMSE" => :lower,
    "Ratio" => "Closer to 1 is better: unseen data are predicted about as well as the fitted data.")
hcell(h::AbstractString) = "<th>$h</th>"
hcell(h::Tuple) = "<th>$(h[1])$(tip(h[2], header_better(BETTER, h)))</th>"
table(head, rows) = "<div class='tscroll'><table><tr>" * join(hcell.(head)) * "</tr>" * join(rows) * "</table></div>"
trow(cells; cls="") = "<tr$(isempty(cls) ? "" : " class='$cls'")>" * join(cells) * "</tr>"
td(x) = "<td>$x</td>"; tn(x) = "<td class='num'>$x</td>"
pct(v::Vector{Int}, j) = isempty(v) ? NaN : 100 * count(==(j), v) / length(v)
pctfmt(x) = isfinite(x) ? @sprintf("%.0f%%", x) : "—"
withtip(html, text, better) = "<div class='tipwrap'>$(tip(text, better))$html</div>"

# shared tip texts
const T_RMSE = "Root-mean-square error between prediction and data, in cells. Lower is better."
const T_NRMSE = "RMSE divided by the mean of the data it is scored on, so fits of different sizes can be compared. 0.3 means the typical error is about 30% of the typical count."
const T_BAND = "Points: measured mean cells (grey × = excluded outlier). Line: the carried model fitted to the real data. Shaded band: where 95% of the curves refitted to the bootstrap data sets lie, i.e. the uncertainty of the fitted curve (not of single measurements)."

# horizontal bar chart of shares (0–100)
function bars_svg(title, items; W=760)
    H = 34 + 24 * length(items)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(hesc(title))'><rect width='100%' height='100%' fill='white'/><text x='8' y='18' font-size='13' font-weight='bold'>$(hesc(title))</text>"]
    for (i, (lab, v, hl)) in enumerate(items)
        y = 30 + 24 * (i - 1); w = 380 * (isfinite(v) ? v : 0) / 100
        push!(p, @sprintf("<text x='8' y='%d' font-size='11'>%s</text><rect x='330' y='%d' width='%.1f' height='16' rx='3' fill='%s'/><text x='%.1f' y='%d' font-size='11'>%s</text>",
                          y + 12, hesc(first(lab, 52)), y, max(w, 1), hl ? "#176b87" : "#9fb7c3", 336 + w, y + 12, pctfmt(v)))
    end
    push!(p, "</svg>"); join(p)
end

# ---- bootstrapping page
function write_bootstrap()
    secs = String[]; csvp = NamedTuple[]; csvs = NamedTuple[]
    for g in groups
        r = results[g.id].res; band = results[g.id].band; c = g.cands[1]; hasb = size(r.bic, 1) > 0
        order = sortperm([-pct(r.pick, j) for j in eachindex(g.cands)])
        items = [(g.cands[j].label * " (k=$(length(g.cands[j].params)))", pct(r.pick, j), j == 1) for j in order[1:min(6, end)]]
        pk = ck[:pkgboot][g.id]
        prow = String[]
        for (i, n) in enumerate(c.pnames)
            l, m, h = hasb ? qs(r.pars[1][:, i]) : (NaN, NaN, NaN)
            pl2, ph2 = NaN, NaN
            if pk !== nothing
                s = pk.summary[pk.summary.parameter_index .== i, :]
                nrow(s) == 1 && ((pl2, ph2) = (s.ci95_lower[1], s.ci95_upper[1]))
            end
            lo, hi = c.bounds[i]
            atb = (isfinite(h) && hi - h < 1e-3 * (hi - lo)) || (isfinite(l) && l - lo < 1e-3 * (hi - lo))
            push!(prow, trow([td(pl(n)), tn(fmt(c.params[i])), tn(f3(m)), tn(isfinite(l) ? "$(fmt(l)) – $(fmt(h))" : "—"),
                            tn(isfinite(pl2) ? "$(fmt(pl2)) – $(fmt(ph2))" : "—"), td(atb ? "<span class='flag'>block interval reaches a bound</span>" : "")]))
            push!(csvp, (group=g.id, model=c.name, param=n, estimate=c.params[i], block_median=m, block_lo=l, block_hi=h, iid_lo=pl2, iid_hi=ph2))
        end
        for j in eachindex(g.cands)
            push!(csvs, (group=g.id, candidate=g.cands[j].name, k=length(g.cands[j].params), bic=g.cands[j].bic,
                         pct_chosen=pct(r.pick, j), pct_lowest_bic=pct(r.winner, j), optimal=(j == 1)))
        end
        acf = residual_acf(g)
        tmax = g.kind == :s1b ? 20.0 : XMAX
        ts, sol = dense(c.model, c.u0, c.params, tmax; dt=0.1)
        series = zoo_series(g); panels = zoo_panels(g)
        ymax = g.kind == :s1b ? ceil(maximum(maximum(s.y) for s in series) * 1.2; sigdigits=2) : YMAX
        pn = zoo_main_panel(g)
        plot = zoo_band_svg("$(g.title): optimal model, 95% block-bootstrap band" * (isempty(panels[pn]) ? "" : " · " * panels[pn]),
                            [s for s in series if s.panel == pn], ts, sol.u, band, tmax, ymax)
        bars = withtip(bars_svg("How often each candidate is chosen (top 6)", items),
                       "Each bar is the share of the $(size(r.bic, 1)) bootstrap data sets in which the selection rule picks that candidate; the dark bar is the model carried forward. One dominant bar = a stable choice; several similar bars = the data cannot separate those models.", "Better: one dominant dark bar, i.e. a stable choice.")
        push!(secs, "<section class='panel' id='$(g.id)'><h3>$(hesc(g.title))</h3>" *
            "<p><b>Optimal model:</b> $(hesc(c.label)) (k = $(length(c.params))). Chosen by the selection rule in <b>$(pctfmt(pct(r.pick, 1)))</b> of $(size(r.bic, 1)) replicates; lowest BIC in $(pctfmt(pct(r.winner, 1))). Package bootstrap: $(pk === nothing ? "failed" : "$(pk.ok) of $NBOOT replicates completed").</p>" *
            "<div class='grid'><div>$bars</div><div>$(withtip(plot, T_BAND, "Better: the line runs through the points and the shaded band is narrow."))</div></div>" *
            table([("Parameter", "Symbol of the fitted parameter; values carried over from earlier stages are fixed and not bootstrapped here."),
                   ("Estimate", "Value fitted to the real data."),
                   ("Block median", "Median of the values refitted to the bootstrap data sets."),
                   ("95% block bootstrap", "95% of the bootstrap refits fall in this range. Residuals are resampled in 3-day blocks, which keeps day-to-day correlation: this is the interval to use."),
                   ("95% within-trajectory bootstrap (package)", "Same idea with GrowthParameterEstimation.bootstrap_joint_fit, which resamples each day's residual independently. It ignores correlation between neighbouring days, so it is usually narrower (over-confident)."),
                   ("", "Flags an interval that runs into the search limit: the data do not bound the parameter on that side.", "Better: no flag, i.e. the data bound the parameter.")], prow) *
            "<p class='muted'>Residual lag-1 autocorrelation of the optimal fit: $(isempty(acf) ? "—" : join([@sprintf("%.2f", a) for a in acf], ", ")) (mean $(isempty(acf) ? "—" : @sprintf("%.2f", mean(acf))))$(tip("Correlation between each day's fit error and the next day's, one value per trajectory. 0 = independent errors; positive = errors come in runs, so treating days as independent overstates certainty. That is why whole 3-day blocks are resampled.", "Closer to 0 is better (independent errors).")). Positive values mean neighbouring days are correlated, which is why whole 3-day blocks of residuals are resampled.</p></section>")
    end
    CSV.write(joinpath(AN_DIR, "bootstrap_parameters.csv"), DataFrame(csvp)); CSV.write(joinpath(AN_DIR, "bootstrap_selection.csv"), DataFrame(csvs))
    summ = [trow([td("<a href='#$(g.id)'>$(hesc(g.title))</a>"), td(hesc(g.cands[1].label)), tn(pctfmt(pct(results[g.id].res.pick, 1))), tn(pctfmt(pct(results[g.id].res.winner, 1)))]) for g in groups]
    page("bootstrap", "Bootstrapping — compartment models", "bootstrap", """
<div class='intro'><p><b>Method.</b> Each fit is resampled $BOOT_B times. The residuals of the most complex candidate are resampled in circular blocks of $BLOCK days within each trajectory (day 0 kept as observed) and added back to its fit; <b>every candidate</b> is refitted to every replicate and the selection rule (fewest free parameters within 3 BIC of the lowest; shared cisplatin structure across lines; Stage 1b shared-vs-separate test) is replayed on each. This gives how often the optimal model is chosen, 95% parameter intervals and the band around the fit. For comparison the optimal model is also passed through <code>GrowthParameterEstimation.bootstrap_joint_fit</code>, which resamples residuals one by one within each trajectory, as the earlier linked-treatment reports did; it ignores day-to-day correlation, so its intervals are usually narrower.</p></div>
<h2>Selection stability</h2>$(table([("Fit", "The data set and stage."), ("Optimal model", "The model carried forward: fewest free parameters within 3 BIC of the lowest."),
    ("Chosen by the rule", "Share of the bootstrap data sets in which the same selection rule picks this model. Above about 70% the choice is robust; low values mean another model is about as good."),
    ("Lowest BIC", "Share of the bootstrap data sets in which this model had the lowest BIC outright, before preferring fewer parameters. It can differ from 'chosen' when a model with more terms often fits slightly better.")], summ))
<h2>Fits</h2>$(join(secs, "\n"))
<p class='muted'>Outputs: <code>bootstrap_parameters.csv</code>, <code>bootstrap_selection.csv</code> in <code>outputs/zoo_analyses/</code>.</p>""")
end

# ---- sensitivity page
function profile_svg(title, prof, est; W=360, H=200)
    L, R, T, B = 44, 12, 24, 30
    vs = [p.v for p in prof]; ds = [min(p.d, 30.0) for p in prof]
    x0, x1 = minimum(vs), maximum(vs); x1 > x0 || (x1 = x0 + 1)
    sx(v) = L + (v - x0) / (x1 - x0) * (W - L - R); sy(d) = H - B - d / 30 * (H - T - B)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(hesc(title))'><rect width='100%' height='100%' fill='white'/><text x='$L' y='15' font-size='12' font-weight='bold'>$(hesc(title))</text>",
         @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#b45309' stroke-dasharray='4,3'/><text x='%d' y='%.1f' font-size='9' fill='#b45309'>3.84</text>", L, W - R, sy(3.84), sy(3.84), W - R - 22, sy(3.84) - 3),
         "<line x1='$L' y1='$(H-B)' x2='$(W-R)' y2='$(H-B)' stroke='#555'/><line x1='$L' y1='$T' x2='$L' y2='$(H-B)' stroke='#555'/>",
         @sprintf("<text x='%d' y='%.1f' font-size='9' text-anchor='end'>0</text><text x='%d' y='%.1f' font-size='9' text-anchor='end'>30</text>", L - 4, sy(0) + 3, L - 4, sy(30) + 3),
         @sprintf("<text x='%d' y='%d' font-size='9'>%s</text><text x='%d' y='%d' font-size='9' text-anchor='end'>%s</text>", L, H - B + 13, fmt(x0), W - R, H - B + 13, fmt(x1)),
         @sprintf("<line x1='%.1f' x2='%.1f' y1='%d' y2='%d' stroke='#9ca3af'/>", sx(est), sx(est), T, H - B),
         "<polyline points='" * join([@sprintf("%.1f,%.1f", sx(v), sy(d)) for (v, d) in zip(vs, ds)], " ") * "' fill='none' stroke='#176b87' stroke-width='2'/>",
         join([@sprintf("<circle cx='%.1f' cy='%.1f' r='2.5' fill='#176b87'/>", sx(v), sy(d)) for (v, d) in zip(vs, ds)]), "</svg>"]
    join(p)
end
const T_PROFILE = "Profile likelihood. x: the value this parameter is fixed at (0.5× to 2× the estimate). y: n·ln(SSE/SSE_min) after refitting the other parameters, 0 at the best fit and capped at 30. Values under the dashed 3.84 line fit the data equally well (95% level). A narrow valley = well determined; a flat curve = poorly determined. Grey line: the estimate."

# grouped bars of Sobol indices for one output: first-order (dark) and total (light), with 95% whiskers
function sobol_svg(title, names, rows; W=560)
    rs = sort(rows; by=r -> -(isfinite(r.ST) ? r.ST : 0.0))
    H = 46 + 26 * length(rs); L = 120; S = W - L - 60
    sx(v) = L + clamp(v, 0, 1) * S
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(hesc(title))'><rect width='100%' height='100%' fill='white'/><text x='8' y='17' font-size='13' font-weight='bold'>$(hesc(title))</text>",
         join([@sprintf("<line x1='%.1f' x2='%.1f' y1='26' y2='%d' stroke='#eef2f4'/><text x='%.1f' y='%d' font-size='9' text-anchor='middle'>%.1f</text>", sx(v), sx(v), H - 16, sx(v), H - 4, v) for v in 0:0.25:1])]
    for (i, r) in enumerate(rs)
        y = 30 + 26 * (i - 1)
        push!(p, @sprintf("<text x='%d' y='%d' font-size='11' text-anchor='end'>%s</text>", L - 6, y + 13, hesc(names[r.param])))
        for (v, lo, hi, dy, col) in ((r.S1, r.S1_lo, r.S1_hi, 0, "#176b87"), (r.ST, r.ST_lo, r.ST_hi, 10, "#9fc6d6"))
            isfinite(v) || continue
            push!(p, @sprintf("<rect x='%d' y='%d' width='%.1f' height='9' fill='%s'/>", L, y + dy, max(sx(v) - L, 0.5), col))
            isfinite(lo) && push!(p, @sprintf("<line x1='%.1f' x2='%.1f' y1='%.1f' y2='%.1f' stroke='#374151'/>", sx(lo), sx(hi), y + dy + 4.5, y + dy + 4.5))
        end
    end
    push!(p, "</svg>"); join(p)
end
const T_SOBOL = "Sobol indices for this output. Dark bar = first-order index S1: the share of the output's variance caused by this parameter alone. Light bar = total index ST: its share including interactions with other parameters. Whiskers: 95% intervals. ST near 0 = the parameter can be fixed without affecting this output; ST much larger than S1 = it matters through interactions."

function write_sensitivity()
    si_all = NamedTuple[]; pr_all = NamedTuple[]; secs = String[]
    for g in groups
        c = g.cands[1]; si = local_si(g); pr = profile_summary(g)
        append!(si_all, si); append!(pr_all, [Base.structdiff(p, NamedTuple{(:profile,)}) for p in pr])
        maxsi = maximum(filter(isfinite, [s.si for s in si]); init=1.0)
        rows = String[]
        for (s, p) in zip(si, pr)
            status = p.status == "identifiable" ? "<span class='ok'>identifiable</span>" : "<span class='flag'>$(p.status)</span>"
            push!(rows, trow([td(pl(s.param)), tn(fmt(s.value)), tn(@sprintf("%.2f", s.si)),
                            td(@sprintf("<svg viewBox='0 0 120 12' style='width:120px;margin:0'><rect width='%.1f' height='10' fill='#176b87'/></svg>", 120 * (isfinite(s.si) ? s.si : 0) / maxsi)),
                            tn(isfinite(p.ci_lo) ? "$(fmt(p.ci_lo)) – $(fmt(p.ci_hi))" : "—"), td(status),
                            td(p.near_bound ? "<span class='flag'>at the edge of its search range</span>" : "inside bounds")]))
        end
        plots = join([withtip(profile_svg("$(s.param)", p.profile, p.estimate), T_PROFILE, "Better: a deep, narrow valley, i.e. a well-determined parameter.") for (s, p) in zip(si, pr)])
        push!(secs, "<section class='panel' id='$(g.id)'><h3>$(hesc(g.title)) · $(hesc(c.label))</h3>" *
            table([("Parameter", "Fitted parameter of this stage's carried model."), ("Estimate", "Value fitted to the data."),
                   ("Sensitivity index", "Local sensitivity: how much the predicted curves change when this parameter is moved ±10% with the others held fixed, divided by 0.1. 1 = proportional response; above 1 = the prediction amplifies changes in it; near 0 = the prediction barely depends on it."),
                   ("", "Bar length = sensitivity index relative to the largest one in this fit.", "Neither: a long bar means the prediction depends strongly on this parameter."),
                   ("Profile 95% range (within 0.5–2×)", "Values of this parameter the data cannot distinguish from the best fit once the other parameters are refitted (profile likelihood within 3.84 of its minimum), searched from 0.5× to 2× the estimate. A range touching 0.5× or 2× may extend further."),
                   ("Profile", "Identifiable: the data pin the value on both sides. One-sided: only an upper or a lower limit is pinned. Flat: within 0.5–2× the data hardly constrain it, usually because it trades off against another parameter."),
                   ("Bound check", "Whether the estimate sits against its search limit (top 2% of the range, or pinned at the lower limit); if so the limit should be widened and the model refitted.")], rows) *
            "<div class='grid'>$plots</div></section>")
    end
    CSV.write(joinpath(AN_DIR, "sensitivity_local.csv"), DataFrame(si_all)); CSV.write(joinpath(AN_DIR, "sensitivity_profile.csv"), DataFrame(pr_all))
    top = sort(filter(s -> isfinite(s.si), si_all); by=s -> -s.si)
    toprows = [trow([td(hesc(gid(s.group).title)), td(pl(s.param)), tn(@sprintf("%.2f", s.si)), tn(@sprintf("%.1f%%", 100 * s.max_rel_change))]) for s in top[1:min(10, end)]]
    nflat = count(p -> p.status != "identifiable", pr_all)
    nb = count(p -> p.near_bound, pr_all)

    # global sensitivity
    G = ck[:gsa]; F = G.full
    gcsv = NamedTuple[]
    for r in F.rows
        push!(gcsv, (scope="treated coculture, all parameters", output=F.outputs[r.output], param=F.names[r.param], lo=F.lo[r.param], hi=F.hi[r.param],
                     S1=r.S1, S1_lo=r.S1_lo, S1_hi=r.S1_hi, ST=r.ST, ST_lo=r.ST_lo, ST_hi=r.ST_hi, n=r.n_ok))
    end
    fullcharts = join([withtip(sobol_svg(o, F.names, [r for r in F.rows if r.output == j]), T_SOBOL * " Output: " * o * ".", "Neither: long bars mark the parameters that drive the prediction's uncertainty, the ones worth measuring more precisely.") for (j, o) in enumerate(F.outputs)])
    src_title = Dict(g.id => g.title for g in groups)
    rangerows = [trow([td(pl(n)), tn("$(fmt(F.lo[i])) – $(fmt(F.hi[i]))"), td(hesc(src_title[F.carried[i]]))]) for (i, n) in enumerate(F.names)]
    outstat = [trow([td(hesc(o)), tn(fmt(first(r.mean_out for r in F.rows if r.output == j))), tn(fmt(first(r.sd_out for r in F.rows if r.output == j)))]) for (j, o) in enumerate(F.outputs)]
    stagerows = String[]
    for g in groups
        S = G.stages[g.id]
        for r in S.rows
            push!(gcsv, (scope=g.id, output=S.outputs[r.output], param=S.names[r.param], lo=S.lo[r.param], hi=S.hi[r.param],
                         S1=r.S1, S1_lo=r.S1_lo, S1_hi=r.S1_hi, ST=r.ST, ST_lo=r.ST_lo, ST_hi=r.ST_hi, n=r.n_ok))
            push!(stagerows, trow([td(hesc(g.title)), td(hesc(S.outputs[r.output])), td(pl(S.names[r.param])), tn("$(fmt(S.lo[r.param])) – $(fmt(S.hi[r.param]))"),
                                   tn(@sprintf("%.2f (%.2f – %.2f)", r.S1, r.S1_lo, r.S1_hi)), tn(@sprintf("%.2f (%.2f – %.2f)", r.ST, r.ST_lo, r.ST_hi))]))
        end
    end
    CSV.write(joinpath(AN_DIR, "sensitivity_global.csv"), DataFrame(gcsv))
    # the parameter with the largest total index for the resistant fraction
    rf = sort([r for r in F.rows if r.output == 4]; by=r -> -r.ST)
    lead = isempty(rf) ? "" : "For the resistant fraction at day 14 the largest total index is $(pl(F.names[rf[1].param])) (ST = $(@sprintf("%.2f", rf[1].ST)))" *
                              (length(rf) > 1 ? ", then $(pl(F.names[rf[2].param])) ($(@sprintf("%.2f", rf[2].ST)))" : "") * "."

    page("sensitivity", "Sensitivity analysis — compartment models", "sensitivity", """
<div class='intro'><p>Three views of how the fitted models depend on their parameters:</p><ul>
<li><b>Global sensitivity (Sobol indices)</b>: every parameter varied together across its whole plausible range (its 95% bootstrap interval), and the variance of an output split between parameters. This is the view that matters for optimal control: which uncertain parameters actually change the predicted outcome.</li>
<li><b>Local sensitivity</b>: the index of <code>GrowthParameterEstimation.parameter_sensitivity_analysis</code>, each fitted parameter moved ±10% with the others held: \\(SI = \\max|\\Delta \\hat y| / \\max(\\hat y, 1) / 0.1\\) over every observed trajectory. (The package function handles one-state models, so the same formula is applied here to the multi-state compartment models.)</li>
<li><b>Profile likelihood</b>: each parameter fixed at 0.5× to 2× its estimate and the others refitted with <code>run_joint_multistart</code>, to see whether the data pin it down; plus the bound check of <code>profile_joint_fit_bounds</code> ($(nb == 0 ? "no estimate sits at a search limit, so no bound needed widening" : "$(nb) estimate(s) sit at a search limit")).</li></ul>
<p class='note'>$(nflat) of $(length(pr_all)) fitted parameters are not fully identifiable from these data. $lead</p></div>

<h2>Global sensitivity: the treated coculture, all parameters $(tip("The Stage 3 model of the 50:50 coculture at 1 µM cisplatin, with every parameter free, including those carried over from Stages 1 and 2. Each is drawn uniformly across its 95% block-bootstrap interval; $(GSA_N) × (k + 2) model runs (Saltelli sampling). Outputs are the predicted counts on day 14.", :none))</h2>
<div class='grid'>$fullcharts</div>
<div class='grid'><div>$(table([("Parameter", "Every parameter of the treated-coculture model."), ("Range varied", "Uniform range used for the draws: the parameter's 95% block-bootstrap interval, i.e. its plausible range given the data."), ("Fitted in", "The stage whose data determine this parameter.")], rangerows))</div>
<div>$(table([("Output", "What is predicted: day-14 counts in the 50:50 treated coculture."), ("Mean", "Average of the output across all parameter draws."), ("SD", "Standard deviation across the draws: how much the output varies over the plausible parameter ranges. Sobol indices split this variance.")], outstat))</div></div>

<h2>Global sensitivity: each stage's own parameters</h2>
$(table([("Fit", "The stage and data set."), ("Output", "Predicted count on the last observed day of the named trajectory."), ("Parameter", "Free parameter of that stage's carried model; carried values are held fixed."),
         ("Range varied", "Uniform range for the draws: the 95% block-bootstrap interval."),
         ("First-order S1 (95%)", "Share of the output's variance caused by this parameter alone. Values near 0 can be slightly negative from sampling noise."),
         ("Total ST (95%)", "Share of the output's variance involving this parameter, including interactions. ST ≈ 0: the parameter could be fixed without changing this output.")], stagerows))

<h2>Local sensitivity: most influential parameters</h2>$(table([("Fit", "The stage and data set."), ("Parameter", "Fitted parameter."), ("Sensitivity index", "Relative change of the prediction per relative change of the parameter (±10%). 1 = proportional; above 1 = amplified."), ("Largest change of the prediction", "The biggest relative change of any predicted count on any day when this parameter is moved ±10%.")], toprows))
<h2>Local sensitivity and profile likelihood, every fit</h2>$(join(secs, "\n"))
<p class='muted'>Outputs: <code>sensitivity_global.csv</code>, <code>sensitivity_local.csv</code>, <code>sensitivity_profile.csv</code> in <code>outputs/zoo_analyses/</code>.</p>""")
end


function intermittent_html()
    haskey(ck, :intermittent) || return ""
    I = ck[:intermittent]
    rows = String[]; plots = String[]
    for line in ("A2780", "A2780cis")
        R = I[line]
        m = mean(vcat(R.xy[:onoff][2], R.xy[:offon][2]))
        rm = sqrt(R.sse_equal / R.n)
        d = R.bic_koff - R.bic_equal
        push!(rows, trow([td(line), tn(fmt(rm)), tn(@sprintf("%.2f", rm / m)), tn(@sprintf("%.2f", R.bic_equal)), tn(isfinite(R.koff_best) ? fmt(R.koff_best) : "—"), tn(isfinite(R.bic_koff) ? @sprintf("%.2f", R.bic_koff) : "—"),
                          td(!isfinite(d) ? "not tested (kernel has no effect compartment)" :
                             @sprintf("%+.2f", d) * (d > -3 ? " <span class='ok'>equal rates kept</span>" : " <span class='flag'>separate washout supported</span>"))]))
        for (sc, lab) in ((:onoff, "1 µM days 0–7, then off"), (:offon, "off, then 1 µM from day 7"))
            x, y, se = R.xy[sc]; f = R.fits[sc]
            ts = collect(1.0:0.1:maximum(x)); pred = f.scale .* int_predict(line, sc, nothing, ts, y[1] / f.scale)
            ymax = ceil(max(maximum(y .+ se), maximum(pred)) * 1.1; sigdigits=2)
            series = [(panel=1, name="observed (area)", color=line == "A2780" ? BLUE : ORANGE, filled=true, dashed=false, state=1, x=x, y=y, xbad=Float64[], ybad=Float64[])]
            svg = zoo_band_svg("$line · $lab: prediction (line) vs data (points)", series, ts, [[v] for v in pred], nothing, 14.0, ymax)
            push!(plots, withtip(svg, "Points: measured cell area (µm², mean of 3 wells) in the intermittent experiment. Line: the Stage 1 growth + Stage 2 drug model, not refitted to these data; only one area-per-cell scale is fitted per trajectory. The drug switch at day 7 tests whether the model's delay and recovery are right.", "Better: the line runs through the points and the shaded band is narrow."))
        end
    end
    """<h2>Intermittent dosing (drug switched on or off at day 7)</h2>
<div class='intro'><p>The intermittent monoculture experiments (30k; 1 µM cisplatin for days 0–7 then removed, or added from day 7) were never used for fitting. The Stage 1 growth and Stage 2 drug models predict them with no refit (values are image areas, so one area-per-cell scale is profiled per trajectory). They also test the one untested assumption that matters for scheduling: after the drug is removed, the drug effect decays at the same rate k as it built up. A separate washout rate \\(k_{off}\\) is fitted and kept only if it passes the 3-BIC rule.</p></div>
$(table([("Line", "Cell line; each row uses both intermittent trajectories of that line."), ("RMSE (µm²)", "Error of the no-refit prediction, in image-area units."),
         ("nRMSE", T_NRMSE), ("BIC, equal rates", "BIC of the prediction with decay rate = build-up rate (the model as fitted)."),
         ("Best k_off (/day)", "Separate decay rate after the drug is removed that fits best."), ("BIC with k_off", "BIC with that extra parameter."),
         ("ΔBIC", "BIC(with k_off) − BIC(equal rates). Negative by more than 3 would justify the extra washout parameter.")], rows))
<div class='grid'>$(join(plots))</div>"""
end

# ---- validation page
function write_validation()
    V = ck[:validation]
    CSV.write(joinpath(AN_DIR, "validation_blocked.csv"), DataFrame(V.rows))
    rows = [trow([td(hesc(gid(r.group).title)), td(hesc(r.scheme)), td(hesc(r.block)), tn(r.n_test),
                tn(isfinite(r.rmse_heldout) ? fmt(r.rmse_heldout) : "failed"), tn(fmt(r.rmse_insample)),
                tn(isfinite(r.nrmse_heldout) ? @sprintf("%.2f", r.nrmse_heldout) : "—"), tn(@sprintf("%.2f", r.nrmse_insample)),
                tn(isfinite(r.rmse_heldout) ? @sprintf("%.2f", r.rmse_heldout / r.rmse_insample) : "—")]) for r in V.rows]
    prim = V.ce1[1]; m = heldout_ce1_model(βN2, βC2)
    ts, sol = dense(m, c3.u0, prim.params, XMAX; dt=0.1)
    series = zoo_series(g3)
    T_CE1 = "Points: the measured treated (Ce1, 1 µM) coculture, never used to fit these lines. Lines: the prediction from crowding fitted to the untreated cocultures plus the monoculture drug fits. Blue = sensitive A2780, orange = resistant A2780cis. Close agreement = the model generalises to a condition it has not seen."
    plots = join([withtip(zoo_band_svg("Held-out Ce1 · $pt: prediction (lines) vs data (points)", [s for s in series if s.panel == pn && s.dashed], ts, sol.u, nothing, XMAX, YMAX), T_CE1, "Better: the lines run through the points.")
                  for (pn, pt) in enumerate(zoo_panels(g3))])
    fitplots = join([withtip(zoo_band_svg("Ce0 · $pt: crowding fitted here", [s for s in series if s.panel == pn && !s.dashed], ts, sol.u, nothing, XMAX, YMAX),
                             "The untreated (Ce0) cocultures that the crowding terms were fitted to for this validation. Points: data; lines: the fit.", "Better: the lines run through the points.")
                     for (pn, pt) in enumerate(zoo_panels(g3))])
    n1 = nobs(g3.datasets[CE1]); m1 = mean(reduce(vcat, [d.y for d in g3.datasets[CE1]]))
    ce1rows = vcat([trow([td(hesc(v.variant)), td(join(["$(pl(n)) = $(fmt(x))" for (n, x) in zip(c3.pnames[HO_RI], v.params)], ", ")),
                        tn(fmt(v.pred.rmse)), tn(@sprintf("%.2f", v.pred.rmse / m1))]; cls=i == 1 ? "win" : "") for (i, v) in enumerate(V.ce1)],
                   [trow([td("In-sample reference: the Stage 3 optimal model fitted with Ce1"), td("—"), tn(fmt(V.ce1_full.rmse)), tn(@sprintf("%.2f", V.ce1_full.rmse / m1))]),
                    trow([td(get(V, :ce1_baseline, "no_change") == "logistic" ? "Untreated-growth reference: the same logistic growth and crowding fitted to Ce0, no drug term" :
                              "No-change reference: each trajectory stays at its day-0 count"), td("—"), tn(fmt(V.ce1_baseline_rmse)), tn(@sprintf("%.2f", V.ce1_baseline_rmse / m1))])])
    page("validation", "Validation — compartment models", "validation", """
<div class='intro'><p>Three kinds of validation, each refitting the carried model without the held-out data and then predicting it:</p><ul>
<li><b>Held-out treated cocultures (Ce1).</b> The same held-out set as the neural-network track. The Stage 3 model's crowding terms ($(join(HO_NAMES, ", "))) are fitted to the untreated Ce0 cocultures only; cisplatin comes from the Stage 2 monoculture fits (\\(\\beta_N = $(fmt(βN2))\\), \\(\\beta_C = $(fmt(βC2))\\), with their delays$(any(L -> L.action == :hill, c3.sys.lines) ? " and Hill slopes" : "")$(any(ns_on, c3.sys.lines) ? "; each kill is scaled by its line's Norton–Simon growth factor" : "")). The Ce1 cocultures ($(n1) points) are then predicted without refitting.</li>
<li><b>Blocked validation</b> with <code>GrowthParameterEstimation.blocked_joint_validation</code>, as in the earlier linked-treatment reports: each whole block of trajectories (one dose, one seeding ratio, one experimental run) is held out in turn, the model is refitted to the rest from 3 starts, and the held-out block is predicted.</li>
<li><b>Forecast.</b> Each fit is refitted to the first 70% of days and predicts the remaining days.</li></ul>
<p>Scores: RMSE in cells; nRMSE = RMSE / mean of the held-out data. The ratio held-out / in-sample compares with the full fit on the same points (1 = no loss from holding the data out).</p></div>
<h2>Held-out treated cocultures (Ce1)</h2>
$(table([("Prediction", "First row: the real test (crowding fitted without Ce1, drug from the monocultures). In-sample reference: the Stage 3 fit that saw Ce1, the best case. No-change reference: counts that never change, a floor any useful model must beat."),
         ("Crowding fitted to Ce0", "Carrying capacity, crowding coefficients and sensitive death rate refitted to the untreated cocultures only."),
         ("RMSE (cells)", T_RMSE), ("nRMSE", T_NRMSE)], ce1rows))
<div class='grid'>$plots</div>
<details><summary>The Ce0 fits used for the prediction</summary><div class='grid'>$fitplots</div></details>
$(intermittent_html())
<h2>Blocked validation and forecasts</h2>
$(table([("Fit", "The stage and data set."), ("Scheme", "What is held out: one whole dose, seeding ratio or experimental run (blocked validation), or the last 30% of days (forecast)."),
         ("Held-out block", "The part predicted without being used for fitting."), ("Points", "Number of held-out measurements scored."),
         ("Held-out RMSE", "Error of the refitted model on the held-out data, in cells. " * T_RMSE), ("In-sample RMSE", "Error of the full fit (which saw these data) on the same points: the best case."),
         ("Held-out nRMSE", T_NRMSE), ("In-sample nRMSE", "In-sample RMSE divided by the mean of the same data."),
         ("Ratio", "Held-out ÷ in-sample RMSE. 1 = holding the data out costs nothing; 2 = errors double on unseen data.")], rows))
<p class='muted'>The low-resource + cisplatin coculture is now fitted as Stage 4 (at the user's request), so it is not used as held-out validation; its forecast (fit the first 70% of days, predict the rest) is in the table above. Outputs: <code>validation_blocked.csv</code> in <code>outputs/zoo_analyses/</code>.</p>""")
end

write_bootstrap()
write_sensitivity()
write_validation()
println("\nDone.")
