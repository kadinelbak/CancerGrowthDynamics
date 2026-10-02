# Resource-depletion model search. Each culture carries a nutrient level R (fresh 10% FBS medium = 1)
# that the cells consume; the medium is replaced every other day (days 2, 4, ..., 12), resetting R.
# Every combination of mechanisms is fitted once, jointly, to all training cultures, then predicts
# the held-out low-resource + cisplatin coculture with no refit. Fitting uses BlackBoxOptim
# (bounded differential evolution) because the medium changes are discrete events.
#
#   julia -t 20 run_new_data_resource_search.jl        (MAX_STRUCTS=n, EVALS=n to shorten)

import Pkg
Pkg.activate(@__DIR__)

using CSV, DataFrames, XLSX, Dates, Printf
using OrdinaryDiffEq: ODEProblem, solve, Tsit5, DiscreteCallback
using BlackBoxOptim
using Base.Threads

const ROOT     = normpath(joinpath(@__DIR__, "..", ".."))
const DATA_DIR = joinpath(ROOT, "New Datasets")
const VAL_DIR  = joinpath(DATA_DIR, "Validation")
const OUT_DIR  = joinpath(@__DIR__, "outputs", "new_data_resource_search")
const LR_FILE  = "ALL A2780-cis LR data.xlsx"
const T0       = 1.0                               # day 0 is the nominal seeding count
const CHANGES  = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]  # medium replaced every other day
const R0_LR_FIXED = 0.5 / 10                       # 0.5% vs 10% FBS
const EVALS    = parse(Int, get(ENV, "EVALS", "25000"))
mkpath(OUT_DIR)

# ── data (same cultures and exclusion rule as run_new_data_model_search.jl) ──────
clean(df) = (rename!(df, Symbol.(strip.(replace.(string.(names(df)), "﻿" => "")))); df)
function xy(df)
    df = filter(r -> !ismissing(r.Day) && !ismissing(r[Symbol("Mean Cells")]), df)
    x = Float64.(df.Day); y = Float64.(df[!, Symbol("Mean Cells")]); o = sortperm(x); x[o], y[o]
end
load_csv(path) = xy(clean(CSV.read(path, DataFrame; silencewarnings=true)))
function load_xlsx(fname, sheet)
    raw = XLSX.readxlsx(joinpath(DATA_DIR, fname))[sheet][:]
    header = Symbol.(strip.(string.(raw[1, :])))
    rows = [Any[raw[i, j] for j in eachindex(header)] for i in 2:size(raw, 1)]
    xy(DataFrame([getindex.(rows, j) for j in eachindex(header)], header))
end
function flag_drops(y; frac=0.5, recover=0.8, floor=50.0)
    n = length(y); bad = falses(n)
    for i in 2:n
        y[i] > 0 && any(y[j] == y[i] for j in 1:i-1) && (bad[i] = true)
    end
    prev = y[1]
    for i in 2:n
        bad[i] && continue
        rest = [y[j] for j in i+1:n if !bad[j]]
        later = isempty(rest) ? -Inf : maximum(rest)
        if y[i] < frac*prev && prev - y[i] > floor && later >= recover*prev
            bad[i] = true
        else
            prev = y[i]
        end
    end
    bad
end

envs = NamedTuple[]
function add_env!(name, group, lowres, D, N, C; lr=false)
    prep(tr) = tr === nothing ? nothing : (keep = lr ? .!flag_drops(tr[2]) : trues(length(tr[2])); (x=tr[1][keep], y=tr[2][keep]))
    Np, Cp = prep(N), prep(C)
    start(tr) = tr === nothing ? 0.0 : tr.y[findfirst(==(T0), tr.x)]
    push!(envs, (name=name, group=group, lowres=lowres, D=D, N=Np, C=Cp, N1=start(Np), C1=start(Cp)))
end
add_env!("A2780 untreated", "1 untreated monoculture", false, 0.0, load_csv(joinpath(DATA_DIR, "A2780_Untreated_30k 1.csv")), nothing)
add_env!("A2780cis untreated", "1 untreated monoculture", false, 0.0, nothing, load_csv(joinpath(DATA_DIR, "A2780cis_Untreated_30k.csv")))
add_env!("A2780 low resource", "1b low-resource monoculture", true, 0.0, load_xlsx(LR_FILE, "1 - mono A2780 30k"), nothing; lr=true)
add_env!("A2780cis low resource", "1b low-resource monoculture", true, 0.0, nothing, load_xlsx(LR_FILE, "1 - mono A2780cis 30k"); lr=true)
for (line, files) in (("A2780", [(1.47, "A2780_IC25_1.47um_30k.csv"), (1.0, "A2780_IC50_1um_30k.csv"), (0.62, "A2780_IC75_0.62um_30k.csv")]),
                      ("A2780cis", [(1.47, "A2780cis_IC25_1.47um_30k.csv"), (1.0, "A2780cis_IC50_1um_30k.csv"), (0.62, "A2780cis_IC75_0.62um_30k.csv")]))
    for (D, f) in files
        tr = load_csv(joinpath(DATA_DIR, f))
        add_env!("$line $(D) µM", "2 treated monoculture", false, D, line == "A2780" ? tr : nothing, line == "A2780" ? nothing : tr)
    end
end
for (f, s, r, D, cond) in [("ce0_25S75R.xlsx", "25Sensitive_Untreated", "75Resistant_Untreated", 0.0, "Ce0 25:75"),
                           ("ce0_50S50R.xlsx", "50Sensitive_Untreated", "50Resistant_Untreated", 0.0, "Ce0 50:50"),
                           ("ce0_75S25R.xlsx", "75Sensitive_Untreated", "25Resistant_Untreated", 0.0, "Ce0 75:25"),
                           ("ce1_25S75R.xlsx", "25Sensitive_Treated", "75Resistant_Treated", 1.0, "Ce1 25:75"),
                           ("ce1_50S50R.xlsx", "50Sensitive_Treated", "50Resistant_Treated", 1.0, "Ce1 50:50"),
                           ("ce1_75S25R.xlsx", "75Sensitive_Treated", "25Resistant_Treated", 1.0, "Ce1 75:25")]
    add_env!(cond, "3 Ce0 + Ce1 coculture", false, D, load_xlsx(f, s), load_xlsx(f, r))
end
for run in (1, 2), ratio in ("1-3", "1-1", "3-1")
    add_env!("LR run $run $ratio", "3b low-resource coculture", true, 0.0,
             load_xlsx(LR_FILE, "$run - A2780 30k $ratio"), load_xlsx(LR_FILE, "$run - A2780cis 30k $ratio"); lr=true)
end
const NENV = length(envs)
const ENV_LR = Bool[e.lowres for e in envs]
const ENV_D = Float64[e.D for e in envs]
# observations: (env, population 1 = N / 2 = C, day, value)
obs = [(j, 1, x, y) for (j, e) in enumerate(envs) if e.N !== nothing for (x, y) in zip(e.N.x, e.N.y) if x > T0]
append!(obs, [(j, 2, x, y) for (j, e) in enumerate(envs) if e.C !== nothing for (x, y) in zip(e.C.x, e.C.y) if x > T0])
const N_TRAIN = length(obs)
const SAVE_DAYS = collect(2.0:1.0:14.0)

vN = load_csv(joinpath(VAL_DIR, "ALL A2780 LR + cis(A2780 30k 1-1).csv"))
vC = load_csv(joinpath(VAL_DIR, "ALL A2780 LR + cis(A2780cis 30k 1-1).csv"))
vbN, vbC = flag_drops(vN[2]), flag_drops(vC[2])
v0 = [vN[2][findfirst(==(T0), vN[1])], vC[2][findfirst(==(T0), vC[1])]]

# ── model structures ──────────────────────────────────────────────────────────
structs = [(crowd=c, consume=q, starve=s, drug=d, inter=i, r0=r)
           for c in (:shared, :asym) for q in (:same, :line) for s in (:death, :growth, :both)
           for d in (:delay, :now) for i in (:indep, :amp) for r in (:fixed, :fitted)]
haskey(ENV, "MAX_STRUCTS") && (structs = structs[1:parse(Int, ENV["MAX_STRUCTS"])])

const PSPEC = Dict(  # name => (lower, upper)
    "r_N" => (0.05, 2.5), "r_C" => (0.05, 2.5), "K" => (1000.0, 15000.0), "α_NC" => (0.0, 4.0), "α_CN" => (0.0, 4.0),
    "K_R" => (0.005, 3.0), "q" => (1e-4, 3.0), "q_N" => (1e-4, 3.0), "q_C" => (1e-4, 3.0),
    "δ_N" => (0.0, 2.0), "δ_C" => (0.0, 2.0), "β_N" => (0.0, 4.0), "β_C" => (0.0, 4.0),
    "λ_N" => (0.02, 5.0), "λ_C" => (0.02, 5.0), "h_N" => (0.5, 8.0), "h_C" => (0.5, 8.0), "γ" => (0.0, 10.0), "R0_LR" => (0.005, 1.0))

function pnames(s)
    n = ["r_N", "r_C", "K"]
    s.crowd == :asym && append!(n, ["α_NC", "α_CN"])
    push!(n, "K_R")
    s.consume == :same ? push!(n, "q") : append!(n, ["q_N", "q_C"])
    s.starve != :growth && append!(n, ["δ_N", "δ_C"])
    append!(n, ["β_N", "β_C"])
    s.drug == :delay && append!(n, ["λ_N", "λ_C"])
    append!(n, ["h_N", "h_C"])
    s.inter == :amp && push!(n, "γ")
    s.r0 == :fitted && push!(n, "R0_LR")
    n
end
label(s) = join([s.crowd == :shared ? "shared K, one-for-one" : "shared K + asymmetric α",
                 s.consume == :same ? "same uptake per cell" : "uptake per line",
                 Dict(:death => "starvation kills", :growth => "starvation slows growth", :both => "starvation slows growth and kills")[s.starve],
                 s.drug == :delay ? "delayed kill" : "immediate kill",
                 s.inter == :amp ? "starvation amplifies drug" : "drug independent of nutrients",
                 s.r0 == :fixed ? "LR medium = 5% of normal" : "LR medium level fitted"], " · ")

# nutrients: Monod saturation normalised so fresh normal medium (R = 1) gives f = 1
const StructT = NamedTuple{(:crowd, :consume, :starve, :drug, :inter, :r0), NTuple{6,Symbol}}
struct Model
    s::StructT; ix::Dict{String,Int}; lr::Vector{Bool}; D::Vector{Float64}
end
pget(m::Model, p, n, default) = (i = get(m.ix, n, 0); i > 0 ? p[i] : default)
r0_of(m::Model, p, lowres) = lowres ? (m.s.r0 == :fitted ? p[m.ix["R0_LR"]] : R0_LR_FIXED) : 1.0

function (m::Model)(du, u, p, t)
    s = m.s
    rN, rC, K = p[m.ix["r_N"]], p[m.ix["r_C"]], p[m.ix["K"]]
    aNC, aCN = pget(m, p, "α_NC", 1.0), pget(m, p, "α_CN", 1.0)
    KR = p[m.ix["K_R"]]
    qN = s.consume == :same ? p[m.ix["q"]] : p[m.ix["q_N"]]; qC = s.consume == :same ? qN : p[m.ix["q_C"]]
    δN, δC = pget(m, p, "δ_N", 0.0), pget(m, p, "δ_C", 0.0)
    βN, βC = p[m.ix["β_N"]], p[m.ix["β_C"]]
    AN = s.drug == :delay ? 1 - exp(-p[m.ix["λ_N"]]*t) : 1.0
    AC = s.drug == :delay ? 1 - exp(-p[m.ix["λ_C"]]*t) : 1.0
    hN, hC = p[m.ix["h_N"]], p[m.ix["h_C"]]
    γ = pget(m, p, "γ", 0.0)
    for j in eachindex(m.D)
        N = max(u[3j-2], 0.0); C = max(u[3j-1], 0.0); R = max(u[3j], 0.0)
        f = R*(1 + KR)/(R + KR)                  # nutrient availability, 1 in fresh normal medium
        hunger = max(1 - f, 0.0)
        gf = s.starve == :death ? 1.0 : min(f, 1.0)
        dieN = s.starve == :growth ? 0.0 : δN*hunger
        dieC = s.starve == :growth ? 0.0 : δC*hunger
        D = m.D[j]
        amp = 1 + γ*hunger
        kN = D > 0 ? βN*AN*D^hN/(1 + D^hN)*amp : 0.0
        kC = D > 0 ? βC*AC*D^hC/(1 + D^hC)*amp : 0.0
        crowdN = 1 - (N + aNC*C)/K; crowdC = 1 - (C + aCN*N)/K
        du[3j-2] = rN*gf*N*crowdN - (dieN + kN)*N
        du[3j-1] = rC*gf*C*crowdC - (dieC + kC)*C
        du[3j]   = -(qN*N + qC*C)/1000 * R/(R + KR)   # uptake per 1000 cells per day
    end
end

function simulate(m::Model, p, u0N, u0C; saveat=SAVE_DAYS)
    u0 = reduce(vcat, [[u0N[j], u0C[j], r0_of(m, p, m.lr[j])] for j in eachindex(m.D)])
    reset!(integ) = (for j in eachindex(m.D); integ.u[3j] = r0_of(m, integ.p, m.lr[j]); end)
    cb = DiscreteCallback((u, t, integ) -> t in CHANGES, reset!; save_positions=(false, false))
    solve(ODEProblem(m, u0, (T0, 14.0), p), Tsit5(); callback=cb, tstops=CHANGES, saveat=saveat,
          reltol=1e-6, abstol=1e-6, maxiters=100_000, verbose=false)
end

const U0N = Float64[e.N1 for e in envs]
const U0C = Float64[e.C1 for e in envs]
function train_sse(m::Model, p)
    sol = simulate(m, p, U0N, U0C)
    (sol.retcode == :Success || string(sol.retcode) == "Success") || return 1e14
    S = 0.0
    for (j, pop, x, y) in obs
        S += (sol.u[Int(x) - 1][3j - 3 + pop] - y)^2
    end
    isfinite(S) ? S : 1e14
end

function fit_struct(s)
    names = pnames(s)
    m = Model(s, Dict(n => i for (i, n) in enumerate(names)), ENV_LR, ENV_D)
    res = bboptimize(p -> train_sse(m, p); SearchRange=[PSPEC[n] for n in names], NumDimensions=length(names),
                     Method=:adaptive_de_rand_1_bin_radiuslimited, MaxFuncEvals=EVALS, TraceMode=:silent)
    p = best_candidate(res); S = best_fitness(res)
    (names=names, params=p, sse=S, bic=N_TRAIN*log(S/N_TRAIN) + length(p)*log(N_TRAIN), model=m)
end

sse(v) = sum(abs2, v)
function validate(fit)
    s = fit.model.s
    m = Model(s, fit.model.ix, [true], [1.0])
    sol = simulate(m, fit.params, [v0[1]], [v0[2]]; saveat=collect(T0:0.1:14.0))
    at(t, k) = sol(t)[k]
    mN = .!vbN .& (vN[1] .> T0); mC = .!vbC .& (vC[1] .> T0)
    rN_ = [at(x, 1) - y for (x, y) in zip(vN[1][mN], vN[2][mN])]
    rC_ = [at(x, 2) - y for (x, y) in zip(vC[1][mC], vC[2][mC])]
    n = length(rN_) + length(rC_); S = sse(rN_) + sse(rC_)
    r2(r, y) = 1 - sse(r)/sse(y .- sum(y)/length(y))
    (bic=n*log(S/n), sse=S, rmseN=sqrt(sse(rN_)/length(rN_)), rmseC=sqrt(sse(rC_)/length(rC_)),
     r2N=r2(rN_, vN[2][mN]), r2C=r2(rC_, vC[2][mC]), sol=sol)
end

println("Training: $(NENV) cultures, $(N_TRAIN) points. Structures: $(length(structs)). Evaluations per fit: $(EVALS). Threads: $(nthreads()).")
flush(stdout)
if haskey(ENV, "DEBUG")
    @time f = fit_struct(structs[1]); println(f.bic, " ", f.names, " ", round.(f.params; sigdigits=4)); exit()
end
t_start = time()
results = Vector{Any}(undef, length(structs))
done = Atomic{Int}(0)
@threads for k in eachindex(structs)
    results[k] = try fit_struct(structs[k]) catch err; @warn "fit failed" k err; nothing end
    c = atomic_add!(done, 1) + 1
    c % 8 == 0 && (@printf("  %d/%d fitted (%.0f s)\n", c, length(structs), time() - t_start); flush(stdout))
end
@printf("Fitting done in %.0f s\n", time() - t_start)

rows = NamedTuple[]; curves = NamedTuple[]
ok = [k for k in eachindex(structs) if results[k] !== nothing && isfinite(results[k].bic)]
for k in ok
    s, r = structs[k], results[k]
    v = validate(r)
    push!(rows, (id=k, crowd=string(s.crowd), consume=string(s.consume), starve=string(s.starve), drug=string(s.drug),
                 interaction=string(s.inter), r0=string(s.r0), label=label(s), k=length(r.params), train_bic=r.bic, train_sse=r.sse,
                 val_bic=v.bic, val_sse=v.sse, val_rmse_sensitive=v.rmseN, val_rmse_resistant=v.rmseC,
                 val_r2_sensitive=v.r2N, val_r2_resistant=v.r2C,
                 params=join(["$(n)=$(r.params[i])" for (i, n) in enumerate(r.names)], ";")))
    for (t, u) in zip(v.sol.t, v.sol.u)
        push!(curves, (id=k, t=t, N=u[1], C=u[2], R=u[3]))
    end
end
res = sort(DataFrame(rows), :train_bic)
CSV.write(joinpath(OUT_DIR, "model_search_results.csv"), res)
CSV.write(joinpath(OUT_DIR, "model_search_validation_curves.csv"), DataFrame(curves))

tcurves = NamedTuple[]
for k in unique([res.id[1], res.id[argmin(res.val_bic)]])
    r = results[k]
    sol = simulate(r.model, r.params, U0N, U0C; saveat=collect(T0:0.25:14.0))
    for (i, e) in enumerate(envs), (t, u) in zip(sol.t, sol.u)
        push!(tcurves, (id=k, env=e.name, group=e.group, t=t, N=u[3i-2], C=u[3i-1], R=u[3i]))
    end
end
CSV.write(joinpath(OUT_DIR, "model_search_training_curves.csv"), DataFrame(tcurves))
CSV.write(joinpath(OUT_DIR, "model_search_training_data.csv"), DataFrame(vcat(
    [[(env=e.name, group=e.group, line="N", day=x, value=y) for (x, y) in zip(e.N.x, e.N.y)] for e in envs if e.N !== nothing]...,
    [[(env=e.name, group=e.group, line="C", day=x, value=y) for (x, y) in zip(e.C.x, e.C.y)] for e in envs if e.C !== nothing]...)))

println("\nTop 8 by training BIC:")
for r in eachrow(first(res, 8))
    @printf("  train %8.1f  val %7.2f  R² N %5.2f C %5.2f  k=%2d  %s\n", r.train_bic, r.val_bic, r.val_r2_sensitive, r.val_r2_resistant, r.k, r.label)
end
println("Top 8 by validation BIC:")
for r in eachrow(first(sort(res, :val_bic), 8))
    @printf("  val %7.2f  train %8.1f  R² N %5.2f C %5.2f  k=%2d  %s\n", r.val_bic, r.train_bic, r.val_r2_sensitive, r.val_r2_resistant, r.k, r.label)
end
