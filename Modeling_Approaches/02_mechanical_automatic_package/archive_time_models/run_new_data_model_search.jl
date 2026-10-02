# Global model search. Every combination of mechanisms is fitted ONCE, jointly, to all training
# data (untreated, low-resource and treated monocultures; Ce0 + Ce1 cocultures; low-resource
# cocultures) with one shared parameter set in which resource level, dose and co-culture are
# conditions. Each fitted model then predicts the held-out low-resource + cisplatin coculture
# (New Datasets/Validation) with no refit.
#
# Run with threads:  julia -t 20 run_new_data_model_search.jl   (MAX_STRUCTS=n limits the grid)

import Pkg
Pkg.activate(@__DIR__)

using GrowthParameterEstimation
using CSV, DataFrames, XLSX, Dates, Printf
using OrdinaryDiffEq: ODEProblem, solve, Tsit5
using Base.Threads

const ROOT     = normpath(joinpath(@__DIR__, "..", ".."))
const DATA_DIR = joinpath(ROOT, "New Datasets")
const VAL_DIR  = joinpath(DATA_DIR, "Validation")
const OUT_DIR  = joinpath(@__DIR__, "outputs", "new_data_model_search")
const LR_FILE  = "ALL A2780-cis LR data.xlsx"
const U_R      = (10 - 0.5) / (10 - 0.05)   # 0.5% FBS; normal medium (10%) is u_R = 0
const T0       = 1.0                        # day 0 is the nominal seeding count, so start at day 1
mkpath(OUT_DIR)

# ── data ───────────────────────────────────────────────────────────────────────
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

# an environment = one culture: resource level u, dose D, and its N and C trajectories
envs = NamedTuple[]
function add_env!(name, group, u, D, N, C; lr=false)
    function prep(tr)
        tr === nothing && return nothing
        x, y = tr
        keep = lr ? .!flag_drops(y) : trues(length(y))
        (x=x[keep], y=y[keep])
    end
    Np, Cp = prep(N), prep(C)
    start(tr) = tr === nothing ? 0.0 : tr.y[findfirst(==(T0), tr.x)]
    push!(envs, (name=name, group=group, u=u, D=D, N=Np, C=Cp, N1=start(Np), C1=start(Cp)))
end

add_env!("A2780 untreated", "1 untreated monoculture", 0.0, 0.0, load_csv(joinpath(DATA_DIR, "A2780_Untreated_30k 1.csv")), nothing)
add_env!("A2780cis untreated", "1 untreated monoculture", 0.0, 0.0, nothing, load_csv(joinpath(DATA_DIR, "A2780cis_Untreated_30k.csv")))
add_env!("A2780 low resource", "1b low-resource monoculture", U_R, 0.0, load_xlsx(LR_FILE, "1 - mono A2780 30k"), nothing; lr=true)
add_env!("A2780cis low resource", "1b low-resource monoculture", U_R, 0.0, nothing, load_xlsx(LR_FILE, "1 - mono A2780cis 30k"); lr=true)
for (line, files) in (("A2780", [(1.47, "A2780_IC25_1.47um_30k.csv"), (1.0, "A2780_IC50_1um_30k.csv"), (0.62, "A2780_IC75_0.62um_30k.csv")]),
                      ("A2780cis", [(1.47, "A2780cis_IC25_1.47um_30k.csv"), (1.0, "A2780cis_IC50_1um_30k.csv"), (0.62, "A2780cis_IC75_0.62um_30k.csv")]))
    for (D, f) in files
        tr = load_csv(joinpath(DATA_DIR, f))
        add_env!("$line $(D) µM", "2 treated monoculture", 0.0, D, line == "A2780" ? tr : nothing, line == "A2780" ? nothing : tr)
    end
end
for (f, s, r, D, cond) in [("ce0_25S75R.xlsx", "25Sensitive_Untreated", "75Resistant_Untreated", 0.0, "Ce0 25:75"),
                           ("ce0_50S50R.xlsx", "50Sensitive_Untreated", "50Resistant_Untreated", 0.0, "Ce0 50:50"),
                           ("ce0_75S25R.xlsx", "75Sensitive_Untreated", "25Resistant_Untreated", 0.0, "Ce0 75:25"),
                           ("ce1_25S75R.xlsx", "25Sensitive_Treated", "75Resistant_Treated", 1.0, "Ce1 25:75"),
                           ("ce1_50S50R.xlsx", "50Sensitive_Treated", "50Resistant_Treated", 1.0, "Ce1 50:50"),
                           ("ce1_75S25R.xlsx", "75Sensitive_Treated", "25Resistant_Treated", 1.0, "Ce1 75:25")]
    add_env!(cond, "3 Ce0 + Ce1 coculture", 0.0, D, load_xlsx(f, s), load_xlsx(f, r))
end
for run in (1, 2), ratio in ("1-3", "1-1", "3-1")
    add_env!("LR run $run $ratio", "3b low-resource coculture", U_R, 0.0,
             load_xlsx(LR_FILE, "$run - A2780 30k $ratio"), load_xlsx(LR_FILE, "$run - A2780cis 30k $ratio"); lr=true)
end
const NENV = length(envs)
const ENV_U = Float64[e.u for e in envs]
const ENV_D = Float64[e.D for e in envs]

datasets = NamedTuple[]
for (i, e) in enumerate(envs)
    e.N === nothing || (m = e.N.x .> T0; push!(datasets, (x=e.N.x[m], y=e.N.y[m], state_index=2i-1)))
    e.C === nothing || (m = e.C.x .> T0; push!(datasets, (x=e.C.x[m], y=e.C.y[m], state_index=2i)))
end
u0_train = reduce(vcat, [[e.N1, e.C1] for e in envs])
n_train = sum(length(d.x) for d in datasets)

# validation culture
vN = load_csv(joinpath(VAL_DIR, "ALL A2780 LR + cis(A2780 30k 1-1).csv"))
vC = load_csv(joinpath(VAL_DIR, "ALL A2780 LR + cis(A2780cis 30k 1-1).csv"))
vbN, vbC = flag_drops(vN[2]), flag_drops(vC[2])
const VAL_D = 1.0
v0 = [vN[2][findfirst(==(T0), vN[1])], vC[2][findfirst(==(T0), vC[1])]]

# ── model structures ──────────────────────────────────────────────────────────
const CROWD = (:shared, :asym, :lv)
const LOWRES = (:death, :rate, :cap, :ratecap)
const DRUG = (:delay, :now, :arrest)
const INTER = (:indep, :sens, :prolif)
const DEATH = (:none, :base)
structs = [(crowd=c, lr=l, drug=d, inter=i, death=b) for c in CROWD for l in LOWRES for d in DRUG for i in INTER for b in DEATH
           if !(d == :arrest && i == :prolif)]   # growth arrest has no kill term to scale
haskey(ENV, "MAX_STRUCTS") && (structs = structs[1:parse(Int, ENV["MAX_STRUCTS"])])

const PSPEC = Dict(  # name => (p0, lower, upper)
    "r_N" => (0.8, 0.01, 3.0), "r_C" => (0.6, 0.01, 3.0),
    "K" => (4000.0, 500.0, 20000.0), "K_N" => (4000.0, 500.0, 20000.0), "K_C" => (4000.0, 500.0, 20000.0),
    "α_NC" => (1.0, 0.0, 5.0), "α_CN" => (1.0, 0.0, 5.0),
    "d^LR_N" => (0.15, 0.0, 3.0), "d^LR_C" => (0.5, 0.0, 3.0),
    "ρ_N" => (0.3, 0.0, 1.0), "ρ_C" => (0.8, 0.0, 1.0), "κ" => (0.3, 0.0, 0.95),
    "β_N" => (0.8, 0.0, 5.0), "β_C" => (0.1, 0.0, 5.0), "E_N" => (0.6, 0.0, 1.0), "E_C" => (0.2, 0.0, 1.0),
    "λ_N" => (0.2, 0.01, 10.0), "λ_C" => (0.2, 0.01, 10.0), "h_N" => (3.0, 0.2, 8.0), "h_C" => (1.0, 0.2, 8.0),
    "γ" => (0.5, 0.0, 10.0), "d_N" => (0.05, 0.0, 2.0), "d_C" => (0.05, 0.0, 2.0))

function pnames(s)
    n = ["r_N", "r_C"]
    s.crowd == :shared && push!(n, "K")
    s.crowd == :asym && append!(n, ["K", "α_NC", "α_CN"])
    s.crowd == :lv && append!(n, ["K_N", "K_C", "α_NC", "α_CN"])
    s.lr == :death && append!(n, ["d^LR_N", "d^LR_C"])
    s.lr in (:rate, :ratecap) && append!(n, ["ρ_N", "ρ_C"])
    s.lr in (:cap, :ratecap) && push!(n, "κ")
    s.drug == :arrest ? append!(n, ["E_N", "E_C"]) : append!(n, ["β_N", "β_C"])
    s.drug != :now && append!(n, ["λ_N", "λ_C"])
    append!(n, ["h_N", "h_C"])
    s.inter == :sens && push!(n, "γ")
    s.death == :base && append!(n, ["d_N", "d_C"])
    n
end

label(s) = join([Dict(:shared => "shared K, one-for-one", :asym => "shared K + asymmetric α", :lv => "separate K per line + α")[s.crowd],
                 Dict(:death => "LR: extra death", :rate => "LR: slower growth", :cap => "LR: lower K", :ratecap => "LR: slower growth + lower K")[s.lr],
                 Dict(:delay => "drug: delayed kill", :now => "drug: immediate kill", :arrest => "drug: delayed growth arrest")[s.drug],
                 Dict(:indep => "independent", :sens => "LR amplifies drug", :prolif => "kills dividing cells only")[s.inter],
                 Dict(:none => "no baseline death", :base => "baseline death")[s.death]], " · ")

# One right-hand side for every culture; U and Dv hold each culture's resource level and dose.
function make_rhs(s, names, U::Vector{Float64}, Dv::Vector{Float64})
    ix(n) = something(findfirst(==(n), names), 0)
    i_rN, i_rC = ix("r_N"), ix("r_C")
    i_K, i_KN, i_KC, i_aNC, i_aCN = ix("K"), ix("K_N"), ix("K_C"), ix("α_NC"), ix("α_CN")
    i_dLN, i_dLC, i_ρN, i_ρC, i_κ = ix("d^LR_N"), ix("d^LR_C"), ix("ρ_N"), ix("ρ_C"), ix("κ")
    i_bN, i_bC = s.drug == :arrest ? (ix("E_N"), ix("E_C")) : (ix("β_N"), ix("β_C"))
    i_λN, i_λC, i_hN, i_hC, i_γ, i_dN, i_dC = ix("λ_N"), ix("λ_C"), ix("h_N"), ix("h_C"), ix("γ"), ix("d_N"), ix("d_C")
    arrest, prolif, delayed = s.drug == :arrest, s.inter == :prolif, s.drug != :now
    nenv = length(U)
    function f!(du, u, p, t)
        z = zero(p[1])
        rN, rC = p[i_rN], p[i_rC]
        KN0 = i_KN > 0 ? p[i_KN] : p[i_K]; KC0 = i_KC > 0 ? p[i_KC] : p[i_K]
        aNC = i_aNC > 0 ? p[i_aNC] : one(z); aCN = i_aCN > 0 ? p[i_aCN] : one(z)
        AN = delayed ? 1 - exp(-p[i_λN]*t) : one(z); AC = delayed ? 1 - exp(-p[i_λC]*t) : one(z)
        for j in 1:nenv
            N = max(u[2j-1], z); C = max(u[2j], z)
            ue, D = U[j], Dv[j]
            rNe = i_ρN > 0 ? rN*(1 - p[i_ρN]*ue) : rN
            rCe = i_ρC > 0 ? rC*(1 - p[i_ρC]*ue) : rC
            κf = i_κ > 0 ? (1 - p[i_κ]*ue) : one(z)
            gN = 1 - (N + aNC*C)/(KN0*κf); gC = 1 - (C + aCN*N)/(KC0*κf)
            if D > 0
                HN = D^p[i_hN]/(1 + D^p[i_hN]); HC = D^p[i_hC]/(1 + D^p[i_hC])
                amp = i_γ > 0 ? 1 + p[i_γ]*ue : one(z)
                effN = p[i_bN]*AN*HN*amp; effC = p[i_bC]*AC*HC*amp
            else
                effN = effC = z
            end
            deathN = (i_dN > 0 ? p[i_dN] : z) + (i_dLN > 0 ? p[i_dLN]*ue : z)
            deathC = (i_dC > 0 ? p[i_dC] : z) + (i_dLC > 0 ? p[i_dLC]*ue : z)
            if arrest
                du[2j-1] = rNe*max(1 - effN, z)*N*gN - deathN*N
                du[2j]   = rCe*max(1 - effC, z)*C*gC - deathC*C
            else
                sN = prolif ? clamp(gN, z, one(z)) : one(z); sC = prolif ? clamp(gC, z, one(z)) : one(z)
                du[2j-1] = rNe*N*gN - deathN*N - effN*sN*N
                du[2j]   = rCe*C*gC - deathC*C - effC*sC*C
            end
        end
    end
end

sse(v) = sum(abs2, v)
function validate(s, names, p)
    f = make_rhs(s, names, [U_R], [VAL_D])
    sol = solve(ODEProblem(f, v0, (T0, 14.0), p), Tsit5(); reltol=1e-8, abstol=1e-8)
    mN = .!vbN .& (vN[1] .> T0); mC = .!vbC .& (vC[1] .> T0)
    rN_ = [sol(x)[1] - y for (x, y) in zip(vN[1][mN], vN[2][mN])]
    rC_ = [sol(x)[2] - y for (x, y) in zip(vC[1][mC], vC[2][mC])]
    n = length(rN_) + length(rC_); S = sse(rN_) + sse(rC_)
    r2(r, y) = 1 - sse(r)/sse(y .- sum(y)/length(y))
    (bic=n*log(S/n), sse=S, rmseN=sqrt(sse(rN_)/length(rN_)), rmseC=sqrt(sse(rC_)/length(rC_)),
     r2N=r2(rN_, vN[2][mN]), r2C=r2(rC_, vC[2][mC]), sol=sol)
end

function fit_struct(s; p0=nothing, optimizer=:bfgs, maxiters=2000)
    names = pnames(s)
    p0 = p0 === nothing ? [PSPEC[n][1] for n in names] : p0
    bounds = [(PSPEC[n][2], PSPEC[n][3]) for n in names]
    f = make_rhs(s, names, ENV_U, ENV_D)
    r = run_joint_fit(f, datasets, u0_train, p0; bounds=bounds, initial_time=T0, optimizer=optimizer,
                      maxiters=maxiters, reltol=1e-7, abstol=1e-7)
    (names=names, params=r.params, bic=r.bic, sse=r.sse)
end

println("Training: $(NENV) cultures, $(length(datasets)) trajectories, $(n_train) points. Structures: $(length(structs)). Threads: $(nthreads()).")
if haskey(ENV, "DEBUG")
    @time r = fit_struct(structs[1])
    println(r.bic, " ", r.names, " ", r.params)
    exit()
end
t_start = time()
results = Vector{Any}(undef, length(structs))
done = Atomic{Int}(0)
@threads for k in eachindex(structs)
    s = structs[k]
    fit = try
        fit_struct(s)
    catch err
        try fit_struct(s; optimizer=:nelder_mead, maxiters=6000) catch; nothing end
    end
    results[k] = fit
    c = atomic_add!(done, 1) + 1
    c % 10 == 0 && @printf("  %d/%d fitted (%.0f s)\n", c, length(structs), time() - t_start)
end

# polish the 30 best by training BIC with bounded Nelder-Mead from their BFGS optimum
ok = [k for k in eachindex(structs) if results[k] !== nothing && isfinite(results[k].bic)]
order = sort(ok; by=k -> results[k].bic)
polish = order[1:min(30, end)]
@threads for k in polish
    try
        r2 = fit_struct(structs[k]; p0=results[k].params, optimizer=:nelder_mead, maxiters=8000)
        r2.bic < results[k].bic && (results[k] = r2)
    catch
    end
end
@printf("Fitting done in %.0f s\n", time() - t_start)

rows = NamedTuple[]
curves = NamedTuple[]
for k in ok
    s, r = structs[k], results[k]
    v = validate(s, r.names, r.params)
    push!(rows, (id=k, crowd=string(s.crowd), low_resource=string(s.lr), drug=string(s.drug), interaction=string(s.inter), death=string(s.death),
                 label=label(s), k=length(r.params), train_bic=r.bic, train_sse=r.sse, val_bic=v.bic, val_sse=v.sse,
                 val_rmse_sensitive=v.rmseN, val_rmse_resistant=v.rmseC, val_r2_sensitive=v.r2N, val_r2_resistant=v.r2C,
                 params=join(["$(n)=$(r.params[i])" for (i, n) in enumerate(r.names)], ";")))
    for t in T0:0.1:14
        push!(curves, (id=k, t=t, N=v.sol(t)[1], C=v.sol(t)[2]))
    end
end
res = sort(DataFrame(rows), :train_bic)
CSV.write(joinpath(OUT_DIR, "model_search_results.csv"), res)
CSV.write(joinpath(OUT_DIR, "model_search_validation_curves.csv"), DataFrame(curves))

# training fit curves for the best two by training BIC and the best by validation, for the report
tcurves = NamedTuple[]
for k in unique([res.id[1], res.id[2], res.id[argmin(res.val_bic)]])
    s, r = structs[k], results[k]
    f = make_rhs(s, r.names, ENV_U, ENV_D)
    sol = solve(ODEProblem(f, u0_train, (T0, 14.0), r.params), Tsit5(); reltol=1e-8, abstol=1e-8)
    for (i, e) in enumerate(envs), t in T0:0.25:14
        push!(tcurves, (id=k, env=e.name, group=e.group, t=t, N=sol(t)[2i-1], C=sol(t)[2i]))
    end
end
CSV.write(joinpath(OUT_DIR, "model_search_training_curves.csv"), DataFrame(tcurves))
CSV.write(joinpath(OUT_DIR, "model_search_training_data.csv"), DataFrame(vcat(
    [[(env=e.name, group=e.group, line="N", day=x, value=y) for (x, y) in zip(e.N.x, e.N.y)] for e in envs if e.N !== nothing]...,
    [[(env=e.name, group=e.group, line="C", day=x, value=y) for (x, y) in zip(e.C.x, e.C.y)] for e in envs if e.C !== nothing]...)))

println("\nTop 10 by training BIC (with validation):")
for r in eachrow(first(res, 10))
    @printf("  train %8.1f  val %7.2f  R² N %5.2f C %5.2f  k=%2d  %s\n", r.train_bic, r.val_bic, r.val_r2_sensitive, r.val_r2_resistant, r.k, r.label)
end
vs = sort(res, :val_bic)
println("Top 10 by validation BIC:")
for r in eachrow(first(vs, 10))
    @printf("  val %7.2f  train %8.1f  R² N %5.2f C %5.2f  k=%2d  %s\n", r.val_bic, r.train_bic, r.val_r2_sensitive, r.val_r2_resistant, r.k, r.label)
end
