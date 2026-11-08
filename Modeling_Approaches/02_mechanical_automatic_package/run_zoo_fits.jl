# Staged fits with the autonomous model zoo (zoo_models.jl): no model uses time t directly,
# so every fitted model can be driven by a dose schedule D(t) for optimal control.
#
#   Stage 1   untreated monocultures          logistic / theta-logistic / Gompertz
#   Stage 1b  low-resource monocultures       Stage 1 growth + 16 low-resource mechanisms
#   Stage 2   treated monocultures (3 doses)  Stage 1 growth + 11 cisplatin mechanisms
#   Stage 3   Ce0 + Ce1 cocultures, shared K  Stage 1 growth + each line's Stage 2 mechanism
#             (the resistant line keeps a cisplatin kill: inherited from Stage 2, or refitted no lower
#             than BETA_C_FLOOR × its Stage 2 value; user decision 2026-10-01)
#   Stage 3b  low-resource cocultures         Stage 1 growth + each line's Stage 1b mechanism
#
# At every stage the optimal model is the candidate with the fewest free parameters within 2 BIC of the
# lowest BIC (lower BIC breaks ties); it is carried into later stages. At most 3 compartments per line.
# Every candidate is then refitted to B residual block-bootstrap replicates (bootstrap_tools.jl)
# for 95% parameter intervals, prediction bands and how often the optimal model stays lowest.
# Writes outputs/zoo_fits/ (CSVs + report, copied into docs/) and docs/outputs/results/bundle.js.
#
# The stages, their data (by tag, data/sources.toml), inheritance and candidates are configuration:
# stages.toml. This script runs them in order.
#
# Run:  julia -t auto run_zoo_fits.jl        (BOOT_B = replicates, default 200)
#   DRYRUN=1          build every stage's candidates without fitting (test/test_stages.jl)
#   ZOO_RESUME=1      reuse the Stage 1 / 1b / 2 fits from zoo_mono_stages.jls and finished group bootstraps
#                     from boot_checkpoints/ (when the group's candidates and fits are unchanged)
#   BOOT_WORK_FACTOR  bootstrap work budget, × the evaluations of one refit (default stages.toml, 50; 0 = none)
#   BOOT_SCREEN       two-pass bootstrap: every candidate gets the first BOOT_SCREEN replicates, the rest go only
#                     to candidates inside the selection band at least once (default stages.toml, 40; 0 = all get B)
#   ZOO_REFIT=1       refit every candidate instead of reusing its saved fit (fit_cache/, see cached_fit)
#   ZOO_REUSE_BOOT=1  reuse the bootstrap of every group whose candidates and fits are unchanged
#                     since the last zoo_state.jls (same names, BICs and parameters)

include(joinpath(@__DIR__, "zoo_common.jl"))
include(joinpath(@__DIR__, "zoo_models.jl"))
include(joinpath(@__DIR__, "bootstrap_tools.jl"))

const BOOT_B   = parse(Int, get(ENV, "BOOT_B", "200"))
# ZOO_OUT: a sensitivity run writes its fits, CSVs and report there and publishes nothing
# (no docs copy, site bundle or Final Report)
const ZOO_PUBLISH = !haskey(ENV, "ZOO_OUT")
const ZOO_DIR  = get(ENV, "ZOO_OUT", joinpath(@__DIR__, "outputs", "zoo_fits"))
const ZOO_DOCS = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "zoo_fits")
const BLUE, ORANGE, GREEN, PURPLE = "#0b6e99", "#d95f02", "#3b8b4f", "#8b5cf6"
mkpath(ZOO_DIR)
println("=== Autonomous model zoo: $(Threads.nthreads()) threads, B = $BOOT_B ==="); flush(stdout)

include(joinpath(@__DIR__, "zoo_selection.jl"))

# ── candidate construction ─────────────────────────────────────────────────────

mono_env(D, R) = (D=D, R=R, lines=[1])
const NOALPHA = [[Slot(0, 1.0)]]
inner(bounds, v) = [clamp(x, lo + 1e-3*(hi - lo), hi - 1e-3*(hi - lo)) for (x, (lo, hi)) in zip(v, bounds)]

# X0: nominal starting value per environment and line; x0: the x0_info record of each observed series
# (zoo_data.jl), which marks the series whose starting value is estimated
function make_cand(name, label, family, lines, envs, alpha, X0, pnames, bounds, p0s; x0=nothing, extra...)
    sys = build_system(lines, envs, alpha)
    xf = x0 === nothing ? (x0nom=Float64[], x0est=Int[], x0bounds=Tuple{Float64,Float64}[], x0start=Float64[]) : x0_fields(sys, x0)
    (name=name, label=label, family=family, sys=sys, model=make_system_model(sys), u0=sys_u0(sys, X0),
     pnames=pnames, bounds=bounds, p0s=[inner(bounds, p) for p in p0s], xf..., extra...)
end

function fit_zoo(cands, datasets)
    rows = Vector{Any}(undef, length(cands))
    Threads.@threads :dynamic for i in eachindex(cands)
        rows[i] = fit_cand(cands[i], datasets)          # bootstrap_tools.jl (estimated starting values included)
    end
    failed = [r.name for r in rows if !isfinite(r.bic)]
    isempty(failed) || @warn "no fit" failed
    sort!([r for r in rows if isfinite(r.bic)]; by=r -> r.bic)
end

function show_rows(tag, rows)
    for r in rows
        @printf("  %-9s %-20s k=%d BIC=%8.2f  %s\n", tag, r.name, length(r.params), r.bic,
                join(["$n = $(fmt(v))" for (n, v) in zip(r.pnames, r.params)], ", "))
    end
    flush(stdout)
end

# a single-line mechanism from a catalogue entry, on top of the carried-forward growth law
# ns: Norton–Simon coupling of the kill to growth (stages.toml kill_growth_coupled, drug stages only)
function zoo_line_cands(zoo, X, input, g, envs, X0; ns=false, x0=nothing)
    [begin
        free = Dict(k => i for (i, k) in enumerate(slots))
        L = line_spec(g.growth, kernel, action, input; fixed=Dict(:r => g.r, :K => g.K, :θ => g.θ, :EC50 => EC50_OF_SYM[X]), free=free)
        ns && (L = merge(L, (ns=true,)))
        bnds = [SLOT_BOUNDS[k] for k in slots]
        make_cand(name, label, fam, [L], envs, NOALPHA, X0, [slot_name(k, X, input; combo=(:d0 in slots)) for k in slots], bnds,
                  [[SLOT_STARTS[k][j] for k in slots] for j in 1:3]; slots=slots, x0=x0)
     end for (name, label, fam, kernel, action, slots) in zoo]
end

# ── stages from stages.toml, data from the tidy table by tag ──────────────────

using Serialization
STAGE_CFG["selection"]["rule_bic"] == RULE_BIC || error("stages.toml rule_bic differs from RULE_BIC in zoo_models.jl")
const CATALOGUES = Dict("DRUG_ZOO" => DRUG_ZOO, "LR_ZOO" => LR_ZOO)

# DRYRUN=1: build every stage's candidates without fitting; fitted rows are taken from the last
# zoo_state.jls so later stages inherit exactly what they did (used by test/test_stages.jl)
const DRYRUN = get(ENV, "DRYRUN", "0") == "1"
const DRY_PREV = DRYRUN ? deserialize(joinpath(ZOO_DIR, "zoo_state.jls")) : nothing
const DRY_BUILT = Dict{String,Any}()
function fit_stage(cands, ds, gid; fitter=fit_zoo)
    DRYRUN || return cached_fit(cands, ds, gid, fitter)
    DRY_BUILT[gid] = cands
    DRY_PREV.groups[findfirst(g -> g.id == gid, DRY_PREV.groups)].cands
end

# ── saved fits per candidate (user, 2026-11-08) ───────────────────────────────
# Adding a candidate must not refit the others: every fitted candidate is saved to fit_cache/<group>.jls,
# keyed by everything that decides its fit (name, label, free parameters, bounds, starts, fixed values and
# inherited mechanisms in its lines, crowding, inputs, starting states, the data, the multistart settings and
# the fitting code in bootstrap_tools.jl). A candidate whose key is saved is not refitted. The cache is seeded
# from the fits in the last zoo_state.jls, so the first run after this change already reuses them. A change
# upstream (e.g. a new Stage 2b kill) changes the inherited values in the key, so later stages refit as needed.
#   ZOO_REFIT=1   ignore the saved fits and refit every candidate (the saved fits are then replaced)
const FIT_CACHE_ON = get(ENV, "ZOO_REFIT", "0") != "1"
const FIT_CACHE_DIR = joinpath(ZOO_DIR, "fit_cache")
const FIT_FIELDS = (:params, :bic, :sse, :se, :u0, :x0fit, :nstarts, :nbest)
const FIT_CODE_HASH = bytes2hex(sha256(read(joinpath(@__DIR__, "bootstrap_tools.jl"))))
# the nominal starting states with the estimated ones masked (a fitted row holds its estimates there)
function masked_u0(c)
    u = Vector{Float64}(c.u0)
    hasproperty(c, :x0est) && (u[c.x0est] .= NaN)
    u
end
function cand_key(c, ds)
    parts = [c.name, c.label, c.family, repr(c.pnames), repr(c.bounds), repr(c.p0s), repr(c.sys.lines), repr(c.sys.alpha),
             repr(c.sys.envs), repr(masked_u0(c)), repr(get(c, :x0nom, nothing)), repr(get(c, :x0est, nothing)),
             repr(get(c, :x0bounds, nothing)), repr(get(c, :x0start, nothing)), repr(ds),
             string(EXTRA_STARTS), string(EXTRA_WORK), string(EXTRA_SPREAD), FIT_CODE_HASH]
    bytes2hex(sha256(join(parts, "\0")))
end
fit_record(r) = merge((name=r.name,), NamedTuple{FIT_FIELDS}(Tuple(getfield(r, k) for k in FIT_FIELDS)))
fit_values(rec) = NamedTuple{FIT_FIELDS}(Tuple(getfield(rec, k) for k in FIT_FIELDS))
const FIT_SEED = Ref{Any}(nothing)
# fits of the last finished run (zoo_state.jls), by group id
function seed_fits(gid)
    if FIT_SEED[] === nothing
        FIT_SEED[] = Dict{String,Any}()
        f = joinpath(ZOO_DIR, "zoo_state.jls")
        if isfile(f)
            st = try deserialize(f) catch err; @warn "saved fits: zoo_state.jls not readable" err; nothing end
            st === nothing || for g in st.groups
                FIT_SEED[][g.id] = (cands=g.cands, ds=g.datasets)
            end
        end
    end
    get(FIT_SEED[], gid, nothing)
end
# saved as plain TOML (outputs/zoo_fits/fit_cache/<group>.toml, committed with the run's outputs), so any
# machine with the repository reuses the fits
function read_fit_cache(f)
    isfile(f) || return Dict{String,Any}()
    d = try TOML.parsefile(f) catch err; @warn "saved fits not readable" f err; return Dict{String,Any}() end
    Dict{String,Any}(k => (name=string(get(v, "name", "")), params=Float64.(v["params"]), bic=Float64(v["bic"]), sse=Float64(v["sse"]), se=Float64.(v["se"]),
                           u0=Float64.(v["u0"]), x0fit=Float64.(v["x0fit"]), nstarts=Int(v["nstarts"]), nbest=Int(v["nbest"]))
                     for (k, v) in d if k != "about")
end
function load_fit_cache(gid)
    cache = read_fit_cache(joinpath(FIT_CACHE_DIR, "$gid.toml"))
    sd = seed_fits(gid)
    if sd !== nothing
        for r in sd.cands
            k = try cand_key(r, sd.ds) catch; continue end
            haskey(cache, k) || (cache[k] = fit_record(r))
        end
    end
    cache
end
function cached_fit(cands, ds, gid, fitter)
    FIT_CACHE_ON || isempty(cands) || return save_fits(gid, fitter(cands, ds), ds, Dict{String,Any}())
    cache = load_fit_cache(gid)
    ks = [cand_key(c, ds) for c in cands]
    hit = [haskey(cache, k) && isfinite(cache[k].bic) for k in ks]
    todo = cands[.!hit]
    @printf("  %s: %d of %d candidates reused from saved fits, %d fitted\n", gid, count(hit), length(cands), length(todo)); flush(stdout)
    fresh = isempty(todo) ? Any[] : fitter(todo, ds)
    save_fits(gid, fresh, ds, cache)
    rows = sort!(vcat(Any[merge(c, fit_values(cache[k])) for (c, k, h) in zip(cands, ks, hit) if h], collect(Any, fresh)); by=r -> r.bic)
    fitter === fit_zoo ? rows : order_rule(rows)
end
function save_fits(gid, rows, ds, cache)
    for r in rows
        cache[cand_key(r, ds)] = fit_record(r)
    end
    try
        mkpath(FIT_CACHE_DIR)
        d = Dict{String,Any}("about" => "Saved fits of group $gid (run_zoo_fits.jl cached_fit), keyed by a SHA-256 of each candidate's structure, data and fitting settings.")
        for (k, r) in cache
            d[k] = Dict{String,Any}("name" => r.name, "params" => collect(r.params), "bic" => r.bic, "sse" => r.sse, "se" => collect(r.se),
                                    "u0" => collect(r.u0), "x0fit" => collect(r.x0fit), "nstarts" => r.nstarts, "nbest" => r.nbest)
        end
        open(io -> TOML.print(io, d; sorted=true), joinpath(FIT_CACHE_DIR, "$gid.toml"), "w")
    catch err
        @warn "saved fits not written" gid err
    end
    rows
end

# ── per-line stages: growth laws, then mechanisms on top of the carried growth ─

function growth_rows(st, L)
    srcs = stage_sources(st; line=L.line); length(srcs) == 1 || error("$(st["id"]) $(L.line): expected one series")
    s = series(srcs[1]["id"]); d = keep(s); K0 = maximum(s.y)
    ds = [dataset(d, 1)]; env = [mono_env(env_inputs(srcs[1])...)]
    rel(k, v) = k == :K ? v * K0 : v
    cands = [begin
                 free = Symbol[Symbol(x) for x in c["free"]]
                 Lsp = line_spec(Symbol(c["growth"]), :none, :loss, :D; fixed=Dict(:EC50 => L.ec50), free=Dict(k => i for (i, k) in enumerate(free)))
                 make_cand(c["name"], c["label"], "Growth", [Lsp], env, NOALPHA, [[s.nom]], ["$(k)_$(L.sym)" for k in free],
                           [(rel(k, b[1]), rel(k, b[2])) for (k, b) in zip(free, c["bounds"])],
                           [[rel(k, v) for (k, v) in zip(free, p)] for p in c["starts"]]; slots=free, x0=[x0_info(s)])
             end for c in st["candidate"]]
    rows = order_rule(fit_stage(cands, ds, "$(st["id"])-$(lowercase(L.line))"))
    b = rows[1]; v = Dict(zip(b.slots, b.params))
    (rows=rows, ds=ds, growth=b.sys.lines[1].growth, r=v[:r], K=v[:K], θ=get(v, :θ, 1.0), X=L.sym)
end

function mechanism_rows(st, L, g)
    srcs = stage_sources(st; line=L.line); input = Symbol(st["input"])
    S = [series(s["id"]) for s in srcs]
    envs = [mono_env(env_inputs(s)...) for s in srcs]
    cands = zoo_line_cands(CATALOGUES[st["catalogue"]], L.sym, input, g, envs, [[s.nom] for s in S];
                           ns=get(st, "kill_growth_coupled", false) && input == :D, x0=[x0_info(s) for s in S])
    # optional restriction to named catalogue entries (stages.toml only = [...]; e.g. a Hill-only variant)
    haskey(st, "only") && (cands = [c for c in cands if c.name in st["only"]]; isempty(cands) && error("$(st["id"]): no candidate left by `only`"))
    obs = sys_obs(cands[1].sys)
    ds = [dataset(keep(s), obs[i]) for (i, s) in enumerate(S)]
    rows = order_rule(fit_stage(cands, ds, "$(st["id"])-$(lowercase(L.line))"))
    if input == :R          # long-run fate under low resource, shown on the pages
        rows = [begin
                    sol = solve(ODEProblem(r.model, r.u0, (0.0, 200.0), r.params), Tsit5(); reltol=1e-8, abstol=1e-8)
                    merge(r, (x30=sol(30.0)[1], x200=sol(200.0)[1], extinct=sol(200.0)[1] < 1.0))
                end for r in rows]
    end
    (rows=rows, ds=ds, best=rows[1])
end

# ── cocultures: two lines per environment, one shared K ───────────────────────

const CC_BOUNDS = Dict(Symbol(k) => (Float64(v[1]), Float64(v[2])) for (k, v) in STAGE_CFG["coculture"]["bounds"])
const CC_PNAME = Dict(:K => "K", :αNC => "α_NC", :αCN => "α_CN", :βN => "β_N", :βC => "β_C",
                      :dN => "d_N", :dC => "d_C", :εN => "ε_N", :εC => "ε_C",
                      # coculture modifiers (zoo_models.jl L.cc, 2026-11-08): kill factor φ (one shared or per line),
                      # protection ρ, cytostatic γ^cs, interference η, and a refitted delay rate k (cytostatic-only variants)
                      :φ => "φ", :φN => "φ_N", :φC => "φ_C", :ρN => "ρ_N", :ρC => "ρ_C", :γN => "γ_N^cs", :γC => "γ_C^cs",
                      :ηN => "η_N", :ηC => "η_C", :kN => "k_N", :kC => "k_C")
const CC_MOD_KEYS = (:φ, :φN, :φC, :ρN, :ρC, :γN, :γC, :ηN, :ηC)

# Rebuild a line from its fitted mechanism row: mechanism values fixed, growth from Stage 1,
# K replaced by the coculture K; `refit` names mechanism slots that are fitted again here.
function cc_line(row, g, Kslot, free::Dict{Symbol,Int}, fixed::Dict{Symbol,Float64}, ec50)
    L0 = row.sys.lines[1]
    fx = Dict{Symbol,Float64}(:r => g.r, :θ => g.θ, :K => 0.0, :EC50 => ec50)
    for (k, v) in zip(row.slots, row.params)
        fx[k] = v
    end
    merge!(fx, fixed)
    L = line_spec(L0.growth, L0.kernel, L0.action, L0.input; fixed=fx, free=free)
    merge(L, (s=merge(L.s, (K=Kslot,)),))
end

# fixed: values used for K / α when not fitted; noinput: switch the low-resource input off;
# betaC_floor: a refitted resistant kill stays ≥ this fraction of its inherited value
# cytostatic: "imm" (the cytostatic effect follows the dose) or "delay" (it follows each line's Stage 2b delay chain)
function make_cc(name, label, free, rows, growth, envs, X0, datasets_max; fixed=Dict{Symbol,Float64}(), noinput=false,
                 family="Coculture", betaC_floor=0.0, kill_growth_coupled=false, x0=nothing, cytostatic="imm")
    idx = Dict(k => i for (i, k) in enumerate(free))
    slot(k, v) = haskey(idx, k) ? Slot(idx[k], 0.0) : Slot(0, get(fixed, k, v))
    Kslot = slot(:K, NaN)
    lf(ks) = Dict(s => idx[k] for (k, s) in ks if haskey(idx, k))
    fN = lf([(:βN, :β), (:dN, :d), (:εN, :ε), (:kN, :k)]); fC = lf([(:βC, :β), (:dC, :d), (:εC, :ε), (:kC, :k)])
    # fixed line values from the candidate (e.g. βC = 0: no resistant cisplatin term) override the inherited ones
    fixline(ks) = Dict{Symbol,Float64}(s => Float64(fixed[k]) for (k, s) in ks if haskey(fixed, k) && !haskey(idx, k))
    xN = merge(noinput ? Dict(:ε => 0.0) : Dict{Symbol,Float64}(), fixline([(:βN, :β), (:dN, :d), (:εN, :ε)]))
    xC = merge(noinput ? Dict(:ε => 0.0) : Dict{Symbol,Float64}(), fixline([(:βC, :β), (:dC, :d), (:εC, :ε)]))
    LN = cc_line(rows[1], growth["A2780"], Kslot, fN, xN, LINES[1].ec50)
    LC = cc_line(rows[2], growth["A2780cis"], Kslot, fC, xC, LINES[2].ec50)
    # Norton–Simon: each line's drug kill scales with its growth factor (crowded cells die less)
    kill_growth_coupled && (LN = merge(LN, (ns=true,)); LC = merge(LC, (ns=true,)))
    # coculture modifiers (only on candidates that use them, so the earlier candidates are rebuilt exactly)
    if any(k -> k in free || haskey(fixed, k), CC_MOD_KEYS)
        cytostatic in ("imm", "delay") || error("$name: cytostatic must be \"imm\" or \"delay\"")
        kf(k) = haskey(idx, :φ) ? Slot(idx[:φ], 0.0) : slot(k, 1.0)
        mk(X) = (kf=kf(Symbol("φ$X")), pr=slot(Symbol("ρ$X"), 0.0), cs=slot(Symbol("γ$X"), 0.0), csimm=cytostatic == "imm",
                 di=slot(Symbol("η$X"), 0.0))
        LN = merge(LN, (cc=mk("N"),)); LC = merge(LC, (cc=mk("C"),))
    end
    alpha = [[Slot(0, 1.0), slot(:αNC, 1.0)], [slot(:αCN, 1.0), Slot(0, 1.0)]]
    inh(k, row) = (i = findfirst(==(k), row.slots); i === nothing ? (k == :k ? 1.0 : 0.5) : row.params[i])
    start = Dict(:K => 1.1datasets_max, :αNC => 1.0, :αCN => 1.0, :βN => inh(:β, rows[1]), :βC => inh(:β, rows[2]),
                 :dN => 0.05, :dC => 0.05, :εN => 1.0, :εC => 1.0,
                 :φ => 1.0, :φN => 1.0, :φC => 1.0, :ρN => 0.5, :ρC => 0.5, :γN => 0.1, :γC => 0.1, :ηN => 0.05, :ηC => 0.05,
                 :kN => inh(:k, rows[1]), :kC => inh(:k, rows[2]))
    alt = Dict(:K => 1.6datasets_max, :αNC => 0.5, :αCN => 1.5, :βN => 0.5start[:βN], :βC => 0.75start[:βC],
               :dN => 0.2, :dC => 0.01, :εN => 0.5, :εC => 2.0,
               :φ => 0.5, :φN => 0.5, :φC => 0.5, :ρN => 2.0, :ρC => 2.0, :γN => 0.4, :γC => 0.4, :ηN => 0.2, :ηC => 0.2,
               :kN => 0.5start[:kN], :kC => 0.5start[:kC])
    alt2 = Dict(:K => 0.8datasets_max, :αNC => 1.5, :αCN => 0.5, :βN => 2start[:βN] + 0.05, :βC => 2start[:βC] + 0.05,
                :dN => 0.01, :dC => 0.2, :εN => 2.0, :εC => 0.5,
                :φ => 0.2, :φN => 0.2, :φC => 0.2, :ρN => 8.0, :ρC => 8.0, :γN => 0.8, :γC => 0.8, :ηN => 0.01, :ηC => 0.01,
                :kN => 2start[:kN], :kC => 2start[:kC])
    bounds = [k == :K ? (CC_BOUNDS[:K][1] * datasets_max, CC_BOUNDS[:K][2] * datasets_max) :
              k == :βC ? (betaC_floor * start[:βC], CC_BOUNDS[:βC][2]) : CC_BOUNDS[k] for k in free]
    make_cand(name, label, family, [LN, LC], envs, alpha, X0, [CC_PNAME[k] for k in free], bounds,
              [[s[k] for k in free] for s in (start, alt, alt2)]; slots=free, fixedcc=fixed, noinput=noinput, mech=rows, x0=x0)
end
ccval(r, k, default) = (i = findfirst(==(k), r.slots); i === nothing ? get(r.fixedcc, k, default) : r.params[i])

function coculture_rows(st)
    E = cc_envs(st)
    envs = [(D=e.inputs.D, R=e.inputs.R, lines=[1, 2]) for e in E]
    X0 = [[e.N.nom, e.C.nom] for e in E]
    x0 = reduce(vcat, [[x0_info(e.N), x0_info(e.C)] for e in E])
    ymax = maximum(maximum(vcat(e.N.y, e.C.y)) for e in E)
    mech = [PL[st["inherits"]["mechanism"]]["A2780"].best, PL[st["inherits"]["mechanism"]]["A2780cis"].best]
    growth = PL[st["inherits"]["growth"]]
    cands = [begin
                 fx = Dict{Symbol,Float64}(Symbol(k) => Float64(v) for (k, v) in get(c, "fixed", Dict()))
                 if haskey(c, "fixed_from")
                     b = CC[c["fixed_from"]].rows[1]
                     merge!(fx, Dict(:K => ccval(b, :K, NaN), :αNC => ccval(b, :αNC, 1.0), :αCN => ccval(b, :αCN, 1.0)))
                 end
                 make_cc(c["name"], c["label"], Symbol[Symbol(x) for x in c["free"]], mech, growth, envs, X0, ymax; fixed=fx,
                         noinput=get(c, "noinput", false), betaC_floor=Float64(get(st, "betaC_floor", 0.0)),
                         kill_growth_coupled=get(c, "kill_growth_coupled", get(st, "kill_growth_coupled", false)), x0=x0,
                         cytostatic=get(c, "cytostatic", "imm"), family=get(c, "family", "Coculture"))
             end for c in st["candidate"]]
    obs = sys_obs(cands[1].sys)
    ds = reduce(vcat, [[dataset(keep(e.N), obs[2i-1]), dataset(keep(e.C), obs[2i])] for (i, e) in enumerate(E)])
    rows = order_rule(fit_stage(cands, ds, st["id"]))
    (rows=rows, ds=ds)
end

# ── run the stages in order ───────────────────────────────────────────────────

const MONO_CKPT = joinpath(ZOO_DIR, "zoo_mono_stages.jls")
const PERLINE = [st["id"] for st in STAGES if st["kind"] in ("growth", "mechanism")]
PL = Dict{String,Any}()          # per-line stage results: PL[stage id][line]
CC = Dict{String,Any}()          # coculture stage results
JOINT_INFO = Dict{String,Any}()
groups = NamedTuple[]
resumed = false
if !DRYRUN && get(ENV, "ZOO_RESUME", "0") == "1" && isfile(MONO_CKPT)
    # reuse the per-line fits of an interrupted run (same catalogue); the joint rules are applied again
    ck = deserialize(MONO_CKPT)
    ck isa Tuple && (ck = Dict("s1" => ck[1], "s1b" => ck[2], "s2" => ck[3]))       # older checkpoint format
    merge!(PL, Dict(k => copy(v) for (k, v) in ck)); resumed = true
    println("Resumed the per-line stages $(join(sort(collect(keys(ck))), ", ")) from ", MONO_CKPT)
end

function apply_joint_rule!(st)
    haskey(st, "joint_rule") || return
    id = st["id"]; R = PL[id]
    sel = st["joint_rule"] == "shared" ? select_shared(R["A2780"].rows, R["A2780cis"].rows) : select_ab(R["A2780"].rows, R["A2780cis"].rows)
    names = st["joint_rule"] == "shared" ? (sel.name, sel.name) : sel.names
    for (L, n) in zip(("A2780", "A2780cis"), names)
        r = order_by_name(R[L].rows, n); R[L] = merge(R[L], (rows=r, best=r[1]))
    end
    JOINT_INFO[id] = sel.info
    if st["joint_rule"] == "shared"
        @printf("Stage %s shared structure: %s (joint BIC %.2f, k=%d)\n", id, sel.name, sel.info.pick.bic, sel.info.pick.k)
    else
        i = sel.info
        @printf("Stage %s: %s (A per line %.2f, k=%d; B shared %s %.2f, k=%d; ΔBIC %+.2f)\n", id, i.useB ? "shared (B)" : "per line (A)",
                i.A.bic, i.A.k, i.B.name, i.B.bic, i.B.k, i.dbic)
    end
end

for st in STAGES
    id = st["id"]; kind = st["kind"]
    println("\n=== Stage $id: $(st["stage_label"]) ==="); flush(stdout)
    if kind in ("growth", "mechanism")
        if !(resumed && haskey(PL, id))
            PL[id] = Dict{String,Any}()
            for L in LINES
                PL[id][L.line] = kind == "growth" ? growth_rows(st, L) : mechanism_rows(st, L, PL[st["inherits"]["growth"]][L.line])
                show_rows(L.line, PL[id][L.line].rows)
            end
        end
        if !DRYRUN && !resumed && id == PERLINE[end]
            serialize(MONO_CKPT, Dict(k => PL[k] for k in PERLINE))     # before the joint rules, as on resume
        end
        apply_joint_rule!(st)
        for L in LINES
            R = PL[id][L.line]
            push!(groups, (id="$id-$(lowercase(L.line))", stage=st["stage_label"], title=replace(st["title"], "{line}" => L.line),
                           datasets=R.ds, u0=R.rows[1].u0, cands=R.rows, line=L.line, kind=Symbol(id)))
        end
    elseif kind == "coculture"
        CC[id] = coculture_rows(st)
        show_rows(st["line"], CC[id].rows)
        push!(groups, (id=id, stage=st["stage_label"], title=st["title"], datasets=CC[id].ds, u0=CC[id].rows[1].u0,
                       cands=CC[id].rows, line=st["line"], kind=Symbol(id)))
    elseif kind == "combination"
        include(joinpath(@__DIR__, "zoo_stage4.jl"))
        inh = st["inherits"]
        rk = PL[inh["resistant_kill"]]["A2780cis"].rows[1]
        cands = stage4_cands(CC[inh["drug"]].rows[1], CC[inh["lowres"]].rows[1], rk.params[findfirst(==(:β), rk.slots)])
        ds = combo_datasets(sys_obs(cands[1].sys))
        rows = fit_stage(cands, ds, id; fitter=fit_stage4)
        show_rows(st["line"], rows)
        push!(groups, (id=id, stage=st["stage_label"], title=st["title"], datasets=ds, u0=rows[1].u0, cands=rows,
                       line=st["line"], kind=Symbol(id)))
    else
        error("stage $id: unknown kind $kind")
    end
end
const JOINT_RULES = [(st["id"], st["joint_rule"]) for st in STAGES if haskey(st, "joint_rule")]

# ── bootstrap, CSVs, checkpoint and pages (skipped in DRYRUN) ─────────────────
if !DRYRUN

tmax_of(g) = Float64(get(stage_cfg(split(g.id, "-")[1]), "plot_tmax", XMAX))
using Serialization
const PREV = get(ENV, "ZOO_REUSE_BOOT", "0") == "1" && isfile(joinpath(ZOO_DIR, "zoo_state.jls")) ? deserialize(joinpath(ZOO_DIR, "zoo_state.jls")) : nothing
# the previous bootstrap of a group, if its candidates and fits are exactly the same as now
function previous_boot(g)
    (PREV === nothing || PREV.B != BOOT_B || get(PREV, :screen, 0) != BOOT_SCREEN) && return nothing
    i = findfirst(h -> h.id == g.id, PREV.groups); i === nothing && return nothing
    h = PREV.groups[i]
    same = length(h.cands) == length(g.cands) && all(a.name == b.name && a.bic == b.bic && a.params == b.params for (a, b) in zip(h.cands, g.cands))
    same ? PREV.results[g.id] : nothing
end
# Per-group bootstrap checkpoints (2026-10-28): each finished group is saved to boot_checkpoints/<id>.jls, so a
# run stopped part-way resumes (ZOO_RESUME=1) without redoing the groups whose candidates and fits are unchanged
const BOOT_CKPT_DIR = joinpath(ZOO_DIR, "boot_checkpoints")
const BOOT_SCREEN = parse(Int, get(ENV, "BOOT_SCREEN", string(selection_cfg("bootstrap_screen", 40))))
const BOOT_KEY = (B=BOOT_B, work_factor=get(ENV, "BOOT_WORK_FACTOR", string(selection_cfg("bootstrap_work_factor", 50))),
                  generator=selection_cfg("bootstrap_generator", "most_complex"), screen=BOOT_SCREEN)
fit_key(g) = [(c.name, c.bic, c.params) for c in g.cands]
function checkpoint_boot(g)
    get(ENV, "ZOO_RESUME", "0") == "1" || return nothing
    f = joinpath(BOOT_CKPT_DIR, "$(g.id).jls"); isfile(f) || return nothing
    ck = try deserialize(f) catch; return nothing end
    ck.key == BOOT_KEY && ck.fits == fit_key(g) ? ck.res : nothing
end
results = Dict{String,Any}()
for g in groups
    t0 = time()
    prev = previous_boot(g)
    ckres = prev === nothing ? checkpoint_boot(g) : nothing
    if prev !== nothing
        res = prev.res
        println("  $(g.id): candidates and fits unchanged, bootstrap reused from the previous run")
    elseif ckres !== nothing
        res = ckres
        println("  $(g.id): candidates and fits unchanged, bootstrap reused from its checkpoint")
    elseif BOOT_B > 0
        res = bootstrap(g; B=BOOT_B)
        try
            mkpath(BOOT_CKPT_DIR)
            serialize(joinpath(BOOT_CKPT_DIR, "$(g.id).jls"), (key=BOOT_KEY, fits=fit_key(g), res=res))
        catch err
            @warn "bootstrap checkpoint not written" g.id err
        end
    else
        nc = length(g.cands); ks = [length(c.params) for c in g.cands]
        res = (cands=g.cands, ks=ks, ref=g.cands[1], bic=fill(NaN, 0, nc), pars=[fill(NaN, 0, k) for k in ks],
               winner=Int[], pick=Int[], selected=1, orig_pick=1)
    end
    # candidate 1 is the one the rule picks on the data; replay the rule on every replicate
    ks = [length(c.params) for c in g.cands]
    picks = size(res.bic, 1) > 0 ? [rule_pick(res.bic[b, :], ks) for b in 1:size(res.bic, 1)] : Int[]
    res = merge(res, (selected=1, pick=picks))
    band = prev !== nothing ? prev.band : BOOT_B > 0 ? bands(g, res; tmax=tmax_of(g), j=1) : nothing
    results[g.id] = (res=res, band=band)
    @printf("  %-12s optimal %-24s k=%d  chosen in %5.1f%% of replicates (lowest BIC in %5.1f%%)  [%.0fs]\n", g.id, g.cands[1].name,
            length(g.cands[1].params), BOOT_B > 0 ? pct(res.pick, 1) : NaN, BOOT_B > 0 ? pct(res.winner, 1) : NaN, time() - t0)
    flush(stdout)
end

apply_joint_replay!(groups, results; rules=JOINT_RULES)

# ── CSVs ──────────────────────────────────────────────────────────────────────

sel_rows = NamedTuple[]; par_rows = NamedTuple[]
for g in groups
    r = results[g.id].res
    for (j, c) in enumerate(g.cands)
        push!(sel_rows, (group=g.id, stage=g.stage, candidate=c.name, label=c.label, family=c.family, k=length(c.params),
                         bic=c.bic, delta_bic=c.bic - g.cands[1].bic, pct_lowest_bic=BOOT_B > 0 ? pct(r.winner, j) : NaN, pct_chosen=BOOT_B > 0 ? pct(r.pick, j) : NaN,
                         n_boot=BOOT_B > 0 ? (hasproperty(r, :nboot) ? r.nboot[j] : BOOT_B) : 0, optimal=(j == 1), params=join(["$n=$(v)" for (n, v) in zip(c.pnames, c.params)], "; ")))
    end
    c = g.cands[1]
    for (i, n) in enumerate(c.pnames)
        l, m, h = BOOT_B > 0 ? qs(r.pars[1][:, i]) : (NaN, NaN, NaN)
        push!(par_rows, (group=g.id, candidate=c.name, param=n, estimate=c.params[i], se=c.se[i], boot_median=m, boot_lo=l, boot_hi=h))
    end
end
for attempt in 1:3
    try
        CSV.write(joinpath(ZOO_DIR, "zoo_model_selection.csv"), DataFrame(sel_rows))
        CSV.write(joinpath(ZOO_DIR, "zoo_parameters.csv"), DataFrame(par_rows))
        break
    catch err
        attempt == 3 ? rethrow() : sleep(2)
    end
end

# starting values of the carried model of every group: nominal seeding, or estimated (zoo_data.jl x0_info)
x0_rows = NamedTuple[]
for g in groups
    c = g.cands[1]; hasproperty(c, :x0nom) && !isempty(c.x0nom) || continue
    obs = sys_obs(c.sys)
    names = g.id == "s4" ? [COMBO_SOURCES.N, COMBO_SOURCES.C] : [s["id"] for s in dataset_tags(g.id)]
    for (i, o) in enumerate(obs)
        push!(x0_rows, (group=g.id, series=i <= length(names) ? names[i] : "state $o", nominal=c.x0nom[i],
                        estimated=o in c.x0est, start=c.u0[o], fold=c.u0[o] / c.x0nom[i]))
    end
end
CSV.write(joinpath(ZOO_DIR, "starting_values.csv"), isempty(x0_rows) ? DataFrame(group=String[]) : DataFrame(x0_rows))

# multistart record of the main fits (decision F5): starts tried per candidate (fixed + seeded Latin hypercube,
# bootstrap_tools.jl fit_cand) and how many of them reached the best BIC (within START_BIC_TOL)
start_rows = NamedTuple[]
for g in groups, c in g.cands
    hasproperty(c, :nstarts) || continue
    push!(start_rows, (group=g.id, candidate=c.name, k=length(c.params), bic=c.bic, starts=c.nstarts, starts_at_best=c.nbest,
                       share_at_best=c.nstarts == 0 ? NaN : c.nbest / c.nstarts, carried=c.name == g.cands[1].name))
end
CSV.write(joinpath(ZOO_DIR, "fit_starts.csv"), isempty(start_rows) ? DataFrame(group=String[]) : DataFrame(start_rows))

# bootstrap optimiser runs that went over their work budget (bootstrap_tools.jl); limit = evaluations, not seconds
CSV.write(joinpath(ZOO_DIR, "time_limit_hits.csv"),
          isempty(LIMIT_HITS) ? DataFrame(fit=String[], optimizer=String[], seconds=Float64[], limit=Float64[]) : DataFrame(LIMIT_HITS))

# checkpoint: zoo_rebuild_pages.jl redraws the report and site from this without refitting
try
    serialize(joinpath(ZOO_DIR, "zoo_state.jls"), (groups=groups, results=results, B=BOOT_B, screen=BOOT_SCREEN, joint=JOINT_INFO,
              config=(stages=STAGES_FILE, sources=SOURCES_FILE, data_table=DATA_TABLE, selection=get(STAGE_CFG, "selection", Dict()))))
catch err
    @warn "checkpoint not written" err
end

include(joinpath(@__DIR__, "export_zoo_bundle.jl"))
include(joinpath(@__DIR__, "zoo_report.jl"))
write_zoo_report(groups, results)
if ZOO_PUBLISH
    export_zoo_bundle(groups, results)
    try
        run(`$(Base.julia_cmd()) --startup-file=no $(joinpath(ROOT, "build_final_report.jl"))`)
    catch err
        @warn "Final report rebuild failed" err
    end
end
println("\nDone.")

end   # !DRYRUN
