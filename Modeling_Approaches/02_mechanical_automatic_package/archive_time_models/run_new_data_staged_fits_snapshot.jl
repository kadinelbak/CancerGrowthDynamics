# Staged fits of the new A2780 / A2780cis data, using the candidate models from
# docs/.../reports/a2780_low_resource_modeling_plan.html:
#   Stage 1  untreated monocultures      logistic / theta-logistic / Gompertz
#   Stage 1b low-resource monocultures   inherited logistic + subtractive loss / rho modifiers / refitted r (30k sheets)
#   Stage 2  treated monocultures        inherited logistic + time-delayed kill (Hill or linear), EC50 = 1 uM
#   Stage 3  Ce0 + Ce1 cocultures jointly asymmetric crowding with one shared K (+/- death) + delayed kill on Ce1
#   Stage 3b low-resource cocultures     inherited crowding + subtractive loss (30k, runs 1 and 2)
# Stages 1b, 2, 3 and 3b also fit a healthy / damaged two-population family (healthy_damaged.jl);
# the report ranks both families together and reports the optimal (lowest-BIC) model.
# Writes CSVs and the "New data fitting" report (also copied into docs/).

import Pkg
Pkg.activate(@__DIR__)

using GrowthParameterEstimation
using CSV, DataFrames, XLSX, Dates, Printf
using OrdinaryDiffEq: ODEProblem, solve, Tsit5
using LinearAlgebra: pinv, diag

const ROOT     = normpath(joinpath(@__DIR__, "..", ".."))
const DATA_DIR = joinpath(ROOT, "New Datasets")
const OUT_DIR  = joinpath(@__DIR__, "outputs", "new_data_fits")
const DOCS_DIR = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "new_data_fits")
const EC50     = 1.0
const COCULTURE_DOSE = 1.0

mkpath(OUT_DIR)

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

# ── fitting helpers ────────────────────────────────────────────────────────────

# Local BFGS from p0 can stall, so each candidate is also fitted with bounded
# Nelder-Mead and the lower-BIC result is kept.
function fit_candidate(model, datasets, u0, p0, bounds)
    best = nothing
    for opt in (:bfgs, :nelder_mead)
        try
            r = run_joint_fit(model, datasets, u0, p0; bounds=bounds, optimizer=opt,
                maxiters = opt == :bfgs ? 3_000 : 8_000, reltol=1e-8, abstol=1e-8)
            if isfinite(r.bic) && (best === nothing || r.bic < best.bic)
                best = (params=r.params, bic=r.bic, sse=r.sse, optimizer=opt)
            end
        catch err
            @warn "fit failed" opt err
        end
    end
    best === nothing && error("all optimizers failed")
    best
end

# Asymptotic standard errors from the least-squares Jacobian: Cov = σ² (JᵀJ)⁻¹, σ² = SSE/(n − k).
# Residuals of mean trajectories are autocorrelated, so these are lower bounds on the real uncertainty.
function param_se(model, datasets, u0, p)
    k = length(p)
    k == 0 && return Float64[]
    tmax = maximum(maximum(d.x) for d in datasets)
    resid(q) = (sol = solve(ODEProblem(model, u0, (0.0, tmax), q), Tsit5(); reltol=1e-10, abstol=1e-10);
                reduce(vcat, [[sol(xi)[d.state_index] - yi for (xi, yi) in zip(d.x, d.y)] for d in datasets]))
    r0 = resid(p); n = length(r0)
    J = zeros(n, k)
    for j in 1:k
        h = 1e-5 * max(abs(p[j]), 1e-2)
        qp = copy(p); qp[j] += h
        qm = copy(p); qm[j] -= h
        J[:, j] = (resid(qp) - resid(qm)) / (2h)
    end
    σ2 = sum(abs2, r0) / max(n - k, 1)
    sqrt.(max.(diag(σ2 * pinv(J' * J)), 0.0))
end

# Optional on-disk cache of finished fits (run with environment variable FIT_CACHE=1): an
# interrupted run then resumes instead of refitting. Off by default. The key covers the
# candidate name, parameter names, start values, bounds and the data, so editing any of them
# refits that candidate; editing the equations under an unchanged name does NOT, so delete
# outputs/new_data_fits/fit_cache.tsv after changing a model.
const FIT_CACHE_ON = get(ENV, "FIT_CACHE", "0") == "1"
const FIT_CACHE_FILE = joinpath(OUT_DIR, "fit_cache.tsv")
const FIT_CACHE = Dict{String,String}()
if FIT_CACHE_ON && isfile(FIT_CACHE_FILE)
    for l in eachline(FIT_CACHE_FILE)
        k, v = split(l, '\t'); FIT_CACHE[k] = v
    end
end
cache_key(c, datasets, u0) = string(c.name, "|", join(c.pnames, ","), "|", c.p0, "|", c.bounds, "|",
                                    hash((u0, [(d.x, d.y, d.state_index) for d in datasets])))
floats(s) = isempty(s) ? Float64[] : parse.(Float64, split(s, ','))

function fit_candidates(cands, datasets, u0)
    rows = NamedTuple[]
    for c in cands
        key = FIT_CACHE_ON ? cache_key(c, datasets, u0) : ""
        if FIT_CACHE_ON && haskey(FIT_CACHE, key)
            b, s_, o, ps, ss = split(FIT_CACHE[key], ';')
            f = (params=floats(ps), bic=parse(Float64, b), sse=parse(Float64, s_), optimizer=Symbol(o))
            se = floats(ss)
        else
            f = fit_candidate(c.model, datasets, u0, c.p0, c.bounds)
            se = try param_se(c.model, datasets, u0, f.params) catch; fill(NaN, length(f.params)) end
            if FIT_CACHE_ON
                v = join([repr(f.bic), repr(f.sse), string(f.optimizer), join(repr.(f.params), ","), join(repr.(se), ",")], ';')
                FIT_CACHE[key] = v
                open(io -> println(io, key, '\t', v), FIT_CACHE_FILE, "a")
            end
        end
        push!(rows, (name=c.name, label=c.label, model=c.model, pnames=c.pnames, params=f.params,
                     bic=f.bic, sse=f.sse, optimizer=f.optimizer, bounds=c.bounds, se=se))
    end
    sort!(rows; by=r -> r.bic)
end

function dense(model, u0, p, tmax; dt=0.1)
    ts = collect(0.0:dt:tmax)
    sol = solve(ODEProblem(model, u0, (0.0, tmax), p), Tsit5(); saveat=ts, reltol=1e-8, abstol=1e-8)
    ts, sol
end

# Scores a fully inherited model (no free parameters) with the same BIC formula as run_joint_fit.
function predict_candidate(name, label, model, datasets, u0)
    tmax = maximum(maximum(d.x) for d in datasets)
    sol = solve(ODEProblem(model, u0, (0.0, tmax), Float64[]), Tsit5(); reltol=1e-8, abstol=1e-8)
    sse = sum(sum((sol(xi)[d.state_index] - yi)^2 for (xi, yi) in zip(d.x, d.y)) for d in datasets)
    n = sum(length(d.x) for d in datasets)
    (name=name, label=label, model=model, pnames=String[], params=Float64[],
     bic=n*log(max(sse, 1e-12)/n), sse=sse, optimizer=:none, bounds=Tuple{Float64,Float64}[], se=Float64[])
end

# A transient drop: below half of the last kept value, by more than 50 cells, with a
# later point recovering to >= 80% of that value. Genuine declines never recover, so
# they are kept; spikes are not touched.
# A value exactly repeating an earlier day's value is treated as a copy/paste error and
# excluded first, so it cannot count as a "recovery". Returns (bad, reason per point).
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

fmt(v) = abs(v) >= 100 ? @sprintf("%.1f", v) : @sprintf("%.4g", v)
param_str(r) = join(["$(n) = $(fmt(v))" for (n, v) in zip(r.pnames, r.params)], ", ")

ranking_rows = NamedTuple[]
all_param_rows = NamedTuple[]   # full-precision values, for downstream scripts
function record!(stage, condition, rows, inherited; line="coculture")
    for r in rows
        push!(ranking_rows, (stage=stage, condition=condition, model=r.name, bic=r.bic,
            delta_bic=r.bic - rows[1].bic, sse=r.sse, params=param_str(r),
            std_errors=join(["$(n) ± $(fmt(s))" for (n, s) in zip(r.pnames, r.se)], ", "), inherited=inherited,
            optimizer=string(r.optimizer)))
        for (i, (n, v)) in enumerate(zip(r.pnames, r.params))
            push!(all_param_rows, (stage=first(split(stage)), line=line, model=r.name, label=r.label, bic=r.bic,
                                   k=length(r.params), param=n, value=v, se=isempty(r.se) ? NaN : r.se[i]))
        end
    end
end

# ── Stage 1: untreated monocultures ────────────────────────────────────────────

println("\n=== Stage 1: untreated monocultures ===")

logistic!(du, u, p, t) = (X = max(u[1], 0.0*p[1]); du[1] = p[1]*X*(1 - X/p[2]))
theta_logistic!(du, u, p, t) = (X = max(u[1], 0.0*p[1]); du[1] = p[1]*X*(1 - (X/p[2])^p[3]))
gompertz!(du, u, p, t) = (X = max(u[1], 1e-8 + 0.0*p[1]); du[1] = p[1]*X*log(p[2]/X))

stage1_lines = [("A2780", "N", "A2780_Untreated_30k 1.csv"), ("A2780cis", "C", "A2780cis_Untreated_30k.csv")]
stage1 = Dict{String,Any}()

for (line, sym, fname) in stage1_lines
    x, y = load_csv_xy(fname)
    K0 = maximum(y)
    kb = (0.5*K0, 20*K0)
    cands = [
        (name="logistic", label="Logistic ★", model=logistic!, pnames=["r_$sym", "K_$sym"],
         p0=[0.6, 1.05K0], bounds=[(1e-4, 5.0), kb]),
        (name="theta_logistic", label="Theta logistic", model=theta_logistic!, pnames=["r_$sym", "K_$sym", "θ_$sym"],
         p0=[0.6, 1.05K0, 1.0], bounds=[(1e-4, 10.0), kb, (0.05, 10.0)]),
        (name="gompertz", label="Gompertz", model=gompertz!, pnames=["r_$sym", "K_$sym"],
         p0=[0.3, 1.05K0], bounds=[(1e-4, 5.0), kb]),
    ]
    rows = fit_candidates(cands, [(x=x, y=y, state_index=1)], [y[1]])
    log_row = only(filter(r -> r.name == "logistic", rows))
    stage1[line] = (x=x, y=y, rows=rows, r=log_row.params[1], K=log_row.params[2], sym=sym, file=fname)
    record!("1 untreated monoculture", "$line ($fname)", rows, "none"; line=line)
    for r in rows
        @printf("  %-9s %-15s BIC=%8.2f  %s\n", line, r.name, r.bic, param_str(r))
    end
end

# ── Stage 1b: low-resource monocultures (30k, run 1) ──────────────────────────

println("\n=== Stage 1b: low-resource monocultures (30k) ===")

const LR_FILE = "ALL A2780-cis LR data.xlsx"
# Normal medium is 10% FBS (u_R = 0, so Stages 1-3 have no low-resource term); low resource is 0.5% FBS.
const FBS_LR = 0.5
const U_R = (10 - FBS_LR) / (10 - 0.05)

outlier_rows = NamedTuple[]
function record_outliers!(stage, trajectory, x, y, bad, reason, sol, idx)
    for i in findall(bad)
        push!(outlier_rows, (stage=stage, trajectory=trajectory, day=x[i], observed=y[i], replacement=sol(x[i])[idx], reason=reason[i]))
    end
end

# fraction of full effect reached after passing through n damage steps at rate λ: 1 - Σ_{k<n} e^{-x} x^k/k!, x = λt
erlang_onset(x, n) = 1 - exp(-x) * sum(x^k / factorial(k) for k in 0:n-1)

function make_lr_mono(kind, rX, KX)
    kind == :loss  && return (du, u, p, t) -> (X = max(u[1], zero(p[1])); du[1] = rX*X*(1 - X/KX) - p[1]*U_R*X)
    kind == :rho   && return (du, u, p, t) -> (X = max(u[1], zero(p[1])); du[1] = p[1]*rX*X*(1 - X/(p[2]*KX)))
    # damage accumulates through n steps before it kills (Erlang onset, mean delay n/λ)
    kind == :accum4 && return (du, u, p, t) -> (X = max(u[1], zero(p[1])); du[1] = rX*X*(1 - X/KX) - p[1]*U_R*erlang_onset(p[2]*t, 4)*X)
    # an immediate loss plus extra loss that accumulates through n damage steps
    kind == :both4  && return (du, u, p, t) -> (X = max(u[1], zero(p[1])); du[1] = rX*X*(1 - X/KX) - (p[1] + p[2]*erlang_onset(p[3]*t, 4))*U_R*X)
    kind == :both8  && return (du, u, p, t) -> (X = max(u[1], zero(p[1])); du[1] = rX*X*(1 - X/KX) - (p[1] + p[2]*erlang_onset(p[3]*t, 8))*U_R*X)
    (du, u, p, t) -> (X = max(u[1], zero(p[1])); du[1] = p[1]*X*(1 - X/KX))
end

lr_mono_sheets = [("A2780", "1 - mono A2780 30k"), ("A2780cis", "1 - mono A2780cis 30k")]
# The day-0 value is the nominal seeding count (SD 0). For A2780cis it is 5x below day 1, which no
# growth rate can bridge, so its day 0 is dropped and day 1 becomes the start (all days shift by -1).
const SHIFT_FROM_DAY1 = Set(["A2780cis"])
stage1b = Dict{String,Any}()
for (line, sheet) in lr_mono_sheets
    s1 = stage1[line]; X = s1.sym
    x, y = load_xlsx_xy(LR_FILE, sheet)
    if line in SHIFT_FROM_DAY1
        keep = x .>= 1
        x, y = x[keep] .- 1, y[keep]
    end
    bad, why = flag_drops(y)
    xf, yf = x[.!bad], y[.!bad]
    cands = [
        (name="lr_subtractive_loss", label="Lineage-specific subtractive loss ★", model=make_lr_mono(:loss, s1.r, s1.K),
         pnames=["d_$(X)^LR"], p0=[0.3], bounds=[(0.0, 5.0)]),
        (name="lr_rho_modifiers", label="Growth-rate and carrying-capacity modifiers", model=make_lr_mono(:rho, s1.r, s1.K),
         pnames=["ρ_r", "ρ_K"], p0=[1.0, 0.8], bounds=[(0.0, 20.0), (1e-3, 3.0)]),
        (name="lr_separate_r", label="Separately fitted growth rate", model=make_lr_mono(:r, s1.r, s1.K),
         pnames=["r_$(X)^LR"], p0=[0.5], bounds=[(1e-4, 5.0)]),
        (name="lr_accumulated_loss4", label="Accumulated low-resource damage, 4 steps", model=make_lr_mono(:accum4, s1.r, s1.K),
         pnames=["d_$(X)^LR", "λ_$(X)^LR"], p0=[1.0, 1.0], bounds=[(0.0, 10.0), (0.01, 10.0)]),
        (name="lr_loss_plus_accum4", label="Immediate loss + accumulated damage, 4 steps", model=make_lr_mono(:both4, s1.r, s1.K),
         pnames=["d_$(X)^LR", "d_$(X)^acc", "λ_$(X)^LR"], p0=[0.5, 1.0, 0.8], bounds=[(0.0, 10.0), (0.0, 10.0), (0.01, 10.0)]),
        (name="lr_loss_plus_accum8", label="Immediate loss + accumulated damage, 8 steps", model=make_lr_mono(:both8, s1.r, s1.K),
         pnames=["d_$(X)^LR", "d_$(X)^acc", "λ_$(X)^LR"], p0=[0.5, 1.0, 1.5], bounds=[(0.0, 10.0), (0.0, 10.0), (0.01, 10.0)]),
    ]
    rows = fit_candidates(cands, [(x=xf, y=yf, state_index=1)], [y[1]])
    star = only(filter(r -> r.name == "lr_subtractive_loss", rows))
    _, solw = dense(rows[1].model, [y[1]], rows[1].params, maximum(x))
    record_outliers!("1b LR monoculture", "$line ($sheet)", x, y, bad, why, solw, 1)
    stage1b[line] = (x=x, y=y, bad=bad, rows=rows, dLR=star.params[1], sym=X, sheet=sheet)
    record!("1b low-resource monoculture", "$line ($sheet); $(count(bad)) point(s) excluded", rows,
            @sprintf("r_%s = %.4f, K_%s = %.1f from Stage 1 logistic; u_R = %.4f (0.5%% FBS)", X, s1.r, X, s1.K, U_R); line=line)
    @printf("  %-9s excluded days: %s\n", line, isempty(findall(bad)) ? "none" : join(Int.(x[bad]), ", "))
    for r in rows
        @printf("  %-9s %-22s BIC=%8.2f  %s\n", line, r.name, r.bic, param_str(r))
    end
end

# ── Stage 2: treated monocultures, three doses jointly ────────────────────────

println("\n=== Stage 2: treated monocultures (EC50 fixed at 1 uM) ===")

stage2_files = Dict(
    "A2780"    => [(1.47, "A2780_IC25_1.47um_30k.csv"), (1.00, "A2780_IC50_1um_30k.csv"), (0.62, "A2780_IC75_0.62um_30k.csv")],
    "A2780cis" => [(1.47, "A2780cis_IC25_1.47um_30k.csv"), (1.00, "A2780cis_IC50_1um_30k.csv"), (0.62, "A2780cis_IC75_0.62um_30k.csv")],
)

hill(D, h) = D^h / (EC50^h + D^h)

function make_treated_hill(r, K, doses)
    function f!(du, u, p, t)
        β, λ, h = p
        A = 1 - exp(-λ*t)
        for i in eachindex(doses)
            X = max(u[i], zero(β))
            du[i] = r*X*(1 - X/K) - β*A*hill(doses[i], h)*X
        end
    end
end

function make_treated_linear(r, K, doses)
    function f!(du, u, p, t)
        β, λ = p
        A = 1 - exp(-λ*t)
        for i in eachindex(doses)
            X = max(u[i], zero(β))
            du[i] = r*X*(1 - X/K) - β*A*doses[i]*X
        end
    end
end

stage2 = Dict{String,Any}()
for (line, _, _) in stage1_lines
    s1 = stage1[line]; sym = s1.sym
    doses = [d for (d, _) in stage2_files[line]]
    data  = [load_csv_xy(f) for (_, f) in stage2_files[line]]
    datasets = [(x=x, y=y, state_index=i) for (i, (x, y)) in enumerate(data)]
    u0 = [y[1] for (_, y) in data]
    cands = [
        (name="delayed_hill_kill", label="Time-delayed Hill kill ★", model=make_treated_hill(s1.r, s1.K, doses),
         pnames=["β_$sym", "λ_$sym", "h_$sym"], p0=[0.5, 0.5, 2.0], bounds=[(0.0, 10.0), (1e-3, 10.0), (0.1, 10.0)]),
        (name="delayed_linear_kill", label="Time-delayed kill, no Hill saturation", model=make_treated_linear(s1.r, s1.K, doses),
         pnames=["β_$sym", "λ_$sym"], p0=[0.3, 0.5], bounds=[(0.0, 10.0), (1e-3, 10.0)]),
    ]
    rows = fit_candidates(cands, datasets, u0)
    best = rows[1]
    λ = best.params[2]
    h = best.name == "delayed_hill_kill" ? best.params[3] : NaN
    stage2[line] = (doses=doses, data=data, u0=u0, rows=rows, best=best, λ=λ, h=h,
                    files=[f for (_, f) in stage2_files[line]], sym=sym)
    record!("2 treated monoculture", "$line, 0.62 / 1.0 / 1.47 uM jointly", rows,
            @sprintf("r_%s = %.4f, K_%s = %.1f from Stage 1 logistic; EC50 = 1 uM", sym, s1.r, sym, s1.K); line=line)
    for r in rows
        @printf("  %-9s %-20s BIC=%8.2f  %s\n", line, r.name, r.bic, param_str(r))
    end
end

# ── Stage 3: Ce0 + Ce1 cocultures jointly, shared K ──────────────────────────

println("\n=== Stage 3: Ce0 + Ce1 cocultures, shared K ===")

# env i -> states (2i-1 = sensitive N, 2i = resistant C)
coculture_layout = [
    (ratio="25S75R", cond="Ce0", file="ce0_25S75R.xlsx", ssheet="25Sensitive_Untreated", rsheet="75Resistant_Untreated", dose=0.0),
    (ratio="50S50R", cond="Ce0", file="ce0_50S50R.xlsx", ssheet="50Sensitive_Untreated", rsheet="50Resistant_Untreated", dose=0.0),
    (ratio="75S25R", cond="Ce0", file="ce0_75S25R.xlsx", ssheet="75Sensitive_Untreated", rsheet="25Resistant_Untreated", dose=0.0),
    (ratio="25S75R", cond="Ce1", file="ce1_25S75R.xlsx", ssheet="25Sensitive_Treated",   rsheet="75Resistant_Treated",   dose=COCULTURE_DOSE),
    (ratio="50S50R", cond="Ce1", file="ce1_50S50R.xlsx", ssheet="50Sensitive_Treated",   rsheet="50Resistant_Treated",   dose=COCULTURE_DOSE),
    (ratio="75S25R", cond="Ce1", file="ce1_75S25R.xlsx", ssheet="75Sensitive_Treated",   rsheet="25Resistant_Treated",   dose=COCULTURE_DOSE),
]

co_data = [(S=load_xlsx_xy(e.file, e.ssheet), R=load_xlsx_xy(e.file, e.rsheet)) for e in coculture_layout]
co_datasets = NamedTuple[]
for (i, d) in enumerate(co_data)
    push!(co_datasets, (x=d.S[1], y=d.S[2], state_index=2i-1))
    push!(co_datasets, (x=d.R[1], y=d.R[2], state_index=2i))
end
u0_s3 = reduce(vcat, [[d.S[2][1], d.R[2][1]] for d in co_data])

rN, rC = stage1["A2780"].r, stage1["A2780cis"].r
λN, λC = stage2["A2780"].λ, stage2["A2780cis"].λ
# Drug effect multiplier at the coculture dose, following each line's Stage 2 winner.
effect(line) = stage2[line].best.name == "delayed_hill_kill" ? hill(COCULTURE_DOSE, stage2[line].h) : COCULTURE_DOSE
HN, HC = effect("A2780"), effect("A2780cis")
const ENV_DOSE = [e.dose for e in coculture_layout]

function make_coculture(with_death::Bool, rN::Float64, rC::Float64, λN::Float64, λC::Float64, HN::Float64, HC::Float64)
    function f!(du, u, p, t)
        K, aNC, aCN = p[1], p[2], p[3]
        dN, dC = with_death ? (p[4], p[5]) : (zero(K), zero(K))
        βN, βC = with_death ? (p[6], p[7]) : (p[4], p[5])
        AN = 1 - exp(-λN*t); AC = 1 - exp(-λC*t)
        for i in 1:6
            N = max(u[2i-1], zero(K)); C = max(u[2i], zero(K))
            killN = ENV_DOSE[i] > 0 ? βN*AN*HN : zero(K)
            killC = ENV_DOSE[i] > 0 ? βC*AC*HC : zero(K)
            du[2i-1] = rN*N*(1 - (N + aNC*C)/K) - dN*N - killN*N
            du[2i]   = rC*C*(1 - (C + aCN*N)/K) - dC*C - killC*C
        end
    end
end

# Treatment fully inherited from each line's Stage 2 winner; only K (and optionally α) is fitted.
# Without α the populations crowd each other one-for-one (α_NC = α_CN = 1).
βN2, βC2 = stage2["A2780"].best.params[1], stage2["A2780cis"].best.params[1]
function make_coculture_inherited_kill(with_alpha::Bool, rN::Float64, rC::Float64, λN::Float64, λC::Float64,
                                       HN::Float64, HC::Float64, βN::Float64, βC::Float64)
    function f!(du, u, p, t)
        K = p[1]
        aNC, aCN = with_alpha ? (p[2], p[3]) : (one(K), one(K))
        AN = 1 - exp(-λN*t); AC = 1 - exp(-λC*t)
        for i in 1:6
            N = max(u[2i-1], zero(K)); C = max(u[2i], zero(K))
            killN = ENV_DOSE[i] > 0 ? βN*AN*HN : 0.0
            killC = ENV_DOSE[i] > 0 ? βC*AC*HC : 0.0
            du[2i-1] = rN*N*(1 - (N + aNC*C)/K) - killN*N
            du[2i]   = rC*C*(1 - (C + aCN*N)/K) - killC*C
        end
    end
end

# Terms named in `pn` are fitted; the rest are off: α = 1 (one-for-one crowding), d = 0, β = 0.
function make_coculture_terms(pn::Vector{String}, rN::Float64, rC::Float64, λN::Float64, λC::Float64, HN::Float64, HC::Float64)
    ix(name) = something(findfirst(==(name), pn), 0)
    iK, iaNC, iaCN, idN, idC, iβN, iβC = ix("K"), ix("α_NC"), ix("α_CN"), ix("d_N"), ix("d_C"), ix("β_N"), ix("β_C")
    function f!(du, u, p, t)
        K = p[iK]
        aNC = iaNC > 0 ? p[iaNC] : one(K); aCN = iaCN > 0 ? p[iaCN] : one(K)
        dN = idN > 0 ? p[idN] : zero(K);   dC = idC > 0 ? p[idC] : zero(K)
        βN = iβN > 0 ? p[iβN] : zero(K);   βC = iβC > 0 ? p[iβC] : zero(K)
        AN = 1 - exp(-λN*t); AC = 1 - exp(-λC*t)
        for i in 1:6
            N = max(u[2i-1], zero(K)); C = max(u[2i], zero(K))
            killN = ENV_DOSE[i] > 0 ? βN*AN*HN : zero(K)
            killC = ENV_DOSE[i] > 0 ? βC*AC*HC : zero(K)
            du[2i-1] = rN*N*(1 - (N + aNC*C)/K) - dN*N - killN*N
            du[2i]   = rC*C*(1 - (C + aCN*N)/K) - dC*C - killC*C
        end
    end
end

Kmax = maximum(maximum(d.y) for d in co_datasets)
Kb = (0.5Kmax, 20Kmax)
s3_term_cand(name, label, pn, p0, bounds) =
    (name=name, label=label, model=make_coculture_terms(pn, rN, rC, λN, λC, HN, HC), pnames=pn, p0=p0, bounds=bounds)
s3_cands = [
    s3_term_cand("crowding_deathN_killN", "Crowding + sensitive death + sensitive kill (no resistant terms)",
        ["K", "α_NC", "α_CN", "d_N", "β_N"], [1.1Kmax, 1.0, 1.0, 0.1, 0.5], [Kb, (0.0, 5.0), (0.0, 5.0), (0.0, 3.0), (0.0, 10.0)]),
    s3_term_cand("crowding_killN", "Crowding + sensitive kill (no death, no resistant kill)",
        ["K", "α_NC", "α_CN", "β_N"], [1.1Kmax, 1.0, 1.0, 0.5], [Kb, (0.0, 5.0), (0.0, 5.0), (0.0, 10.0)]),
    s3_term_cand("sharedK_deathN_killN", "One-for-one crowding + sensitive death + sensitive kill",
        ["K", "d_N", "β_N"], [1.1Kmax, 0.1, 0.5], [Kb, (0.0, 3.0), (0.0, 10.0)]),
    s3_term_cand("sharedK_killN", "One-for-one crowding + sensitive kill (no death)",
        ["K", "β_N"], [1.1Kmax, 0.5], [Kb, (0.0, 10.0)]),
    (name="shared_K_only", label="Shared K only, treatment inherited from Stage 2", model=make_coculture_inherited_kill(false, rN, rC, λN, λC, HN, HC, βN2, βC2),
     pnames=["K"], p0=[1.1Kmax], bounds=[Kb]),
    (name="crowding_inherited_kill", label="Shared K + asymmetric crowding, treatment inherited", model=make_coculture_inherited_kill(true, rN, rC, λN, λC, HN, HC, βN2, βC2),
     pnames=["K", "α_NC", "α_CN"], p0=[1.1Kmax, 1.0, 1.0], bounds=[Kb, (0.0, 5.0), (0.0, 5.0)]),
    (name="crowding_delayed_hill", label="Asymmetric crowding + time-delayed Hill kill on Ce1 ★", model=make_coculture(false, rN, rC, λN, λC, HN, HC),
     pnames=["K", "α_NC", "α_CN", "β_N", "β_C"], p0=[1.1Kmax, 1.0, 1.0, 0.5, 0.1],
     bounds=[Kb, (0.0, 5.0), (0.0, 5.0), (0.0, 10.0), (0.0, 10.0)]),
    (name="crowding_death_delayed_hill", label="Asymmetric crowding + death + time-delayed Hill kill on Ce1", model=make_coculture(true, rN, rC, λN, λC, HN, HC),
     pnames=["K", "α_NC", "α_CN", "d_N", "d_C", "β_N", "β_C"], p0=[1.1Kmax, 1.0, 1.0, 0.02, 0.02, 0.5, 0.1],
     bounds=[Kb, (0.0, 5.0), (0.0, 5.0), (0.0, 3.0), (0.0, 3.0), (0.0, 10.0), (0.0, 10.0)]),
]
s3_rows = fit_candidates(s3_cands, co_datasets, u0_s3)
s3_best = s3_rows[1]
s3_inherited = @sprintf("r_N = %.4f, r_C = %.4f (Stage 1); λ_N = %.4f, λ_C = %.4f, H_N(1 µM) = %.3f, H_C(1 µM) = %.3f (Stage 2)", rN, rC, λN, λC, HN, HC)
record!("3 Ce0 + Ce1 coculture (joint)", "3 ratios × untreated + treated, 12 trajectories", s3_rows, s3_inherited)
for r in s3_rows
    @printf("  %-28s k=%d BIC=%8.2f  %s\n", r.name, length(r.params), r.bic, param_str(r))
    println("      ", join(["$(n) = $(fmt(v)) ± $(fmt(1.96s))" for (n, v, s) in zip(r.pnames, r.params, r.se)], ", "))
end

# ── Stage 3b: low-resource cocultures (30k, runs 1 and 2) ─────────────────────

println("\n=== Stage 3b: low-resource cocultures (30k) ===")

# sheet ratio "a-b" = sensitive:resistant
lr_cc_layout = [(run=run, ratio=ratio, label=Dict("1-3" => "25S75R", "1-1" => "50S50R", "3-1" => "75S25R")[ratio])
                for run in (1, 2) for ratio in ("1-3", "1-1", "3-1")]
lr_cc = NamedTuple[]
lr_cc_datasets = NamedTuple[]
for (i, e) in enumerate(lr_cc_layout)
    sN, sC = "$(e.run) - A2780 30k $(e.ratio)", "$(e.run) - A2780cis 30k $(e.ratio)"
    xN, yN = load_xlsx_xy(LR_FILE, sN); xC, yC = load_xlsx_xy(LR_FILE, sC)
    (bN, wN), (bC, wC) = flag_drops(yN), flag_drops(yC)
    push!(lr_cc, merge(e, (sN=sN, sC=sC, xN=xN, yN=yN, bN=bN, wN=wN, xC=xC, yC=yC, bC=bC, wC=wC)))
    push!(lr_cc_datasets, (x=xN[.!bN], y=yN[.!bN], state_index=2i-1))
    push!(lr_cc_datasets, (x=xC[.!bC], y=yC[.!bC], state_index=2i))
end
u0_3b = reduce(vcat, [[c.yN[1], c.yC[1]] for c in lr_cc])

s3v = Dict(zip(s3_best.pnames, s3_best.params))
K3, aNC3, aCN3 = s3v["K"], get(s3v, "α_NC", 1.0), get(s3v, "α_CN", 1.0)
dN1b, dC1b = stage1b["A2780"].dLR, stage1b["A2780cis"].dLR

# getp maps the free-parameter vector to (K, α_NC, α_CN, d_N^LR, d_C^LR)
function make_lr_cc(getp, rN::Float64, rC::Float64)
    function f!(du, u, p, t)
        K, aNC, aCN, dN, dC = getp(p)
        for i in 1:6
            N = max(u[2i-1], zero(u[2i-1])); C = max(u[2i], zero(u[2i]))
            du[2i-1] = rN*N*(1 - (N + aNC*C)/K) - dN*U_R*N
            du[2i]   = rC*C*(1 - (C + aCN*N)/K) - dC*U_R*C
        end
    end
end

lrmax = maximum(maximum(d.y) for d in lr_cc_datasets)
s3b_rows = fit_candidates([
    (name="lr_cc_sharedK", label="LR coculture: shared K only (one-for-one crowding, d^LR from Stage 1b)", model=make_lr_cc(p -> (p[1], 1.0, 1.0, dN1b, dC1b), rN, rC),
     pnames=["K"], p0=[K3], bounds=[(0.5lrmax, 20lrmax)]),
    (name="lr_cc_K_only", label="LR coculture: K refitted, α from Stage 3, d^LR from Stage 1b", model=make_lr_cc(p -> (p[1], aNC3, aCN3, dN1b, dC1b), rN, rC),
     pnames=["K"], p0=[K3], bounds=[(0.5lrmax, 20lrmax)]),
    (name="lr_cc_refit_d", label="LR coculture ★, K and α inherited, d^LR refitted", model=make_lr_cc(p -> (K3, aNC3, aCN3, p[1], p[2]), rN, rC),
     pnames=["d_N^LR", "d_C^LR"], p0=[dN1b, dC1b], bounds=[(0.0, 5.0), (0.0, 5.0)]),
    (name="lr_cc_refit_all", label="LR coculture ★, K, α and d^LR all refitted", model=make_lr_cc(p -> (p[1], p[2], p[3], p[4], p[5]), rN, rC),
     pnames=["K", "α_NC", "α_CN", "d_N^LR", "d_C^LR"], p0=[K3, aNC3, aCN3, dN1b, dC1b],
     bounds=[(0.5lrmax, 20lrmax), (0.0, 5.0), (0.0, 5.0), (0.0, 5.0), (0.0, 5.0)]),
], lr_cc_datasets, u0_3b)
push!(s3b_rows, predict_candidate("lr_cc_inherited", "LR coculture ★, all parameters inherited (no refit)",
    make_lr_cc(p -> (K3, aNC3, aCN3, dN1b, dC1b), rN, rC), lr_cc_datasets, u0_3b))
sort!(s3b_rows; by=r -> r.bic)
s3b_best = s3b_rows[1]
_, sol3bw = dense(s3b_best.model, u0_3b, s3b_best.params, maximum(maximum(c.xN) for c in lr_cc))
for (i, c) in enumerate(lr_cc)
    record_outliers!("3b LR coculture", "run $(c.run) $(c.ratio) sensitive ($(c.sN))", c.xN, c.yN, c.bN, c.wN, sol3bw, 2i-1)
    record_outliers!("3b LR coculture", "run $(c.run) $(c.ratio) resistant ($(c.sC))", c.xC, c.yC, c.bC, c.wC, sol3bw, 2i)
end
s3b_inherited = @sprintf("r_N = %.4f, r_C = %.4f (Stage 1); K = %.1f, α_NC = %.4f, α_CN = %.4f (Stage 3 winner); d_N^LR = %.4f, d_C^LR = %.4f (Stage 1b ★); u_R = %.4f (0.5%% FBS)",
    rN, rC, K3, aNC3, aCN3, dN1b, dC1b, U_R)
record!("3b low-resource coculture", "3 ratios × runs 1 and 2 (30k), 12 trajectories; $(sum(count(c.bN) + count(c.bC) for c in lr_cc)) point(s) excluded",
        s3b_rows, s3b_inherited)
for r in s3b_rows
    @printf("  %-18s BIC=%8.2f  %s\n", r.name, r.bic, param_str(r))
end
println("  excluded points (all LR): ", length(outlier_rows))

# ── Healthy / damaged track (Stages 1b, 2, 3, 3b) ─────────────────────────────
# A second model family fitted in parallel: each line has healthy (dividing) and damaged
# (non-dividing, still counted) cells; see healthy_damaged.jl. It carries its own parameters
# from stage to stage. These rows are kept apart from stage1b/stage2/s3_rows/s3b_rows, which
# run_bootstrap_fits.jl and run_lr_mechanism_fits.jl read, and are merged only in the report.

include(joinpath(@__DIR__, "healthy_damaged.jl"))
println("\n=== Healthy / damaged models ===")

const HD_FAMILY = "Healthy / damaged"
# Each candidate carries its own layout (u0, observed states), so fit one at a time.
function fit_hd(c)
    datasets = [merge(d, (state_index=c.obs[k],)) for (k, d) in enumerate(c.datasets)]
    r = isempty(c.p0) ? predict_candidate(c.name, c.label, c.model, datasets, c.u0) :
                        only(fit_candidates([c], datasets, c.u0))
    merge(r, (u0=c.u0, obs=c.obs, family=HD_FAMILY, info=c.info(r.params), texf=c.texf, extra=c.extra))
end
hd_rank(rows) = sort(rows; by=r -> r.bic)
hdb(lo=0.0, hi=10.0) = (lo, hi)

# TeX for one line: X observed, damage expression kX (string), repair ρ, exit δ, n compartments.
function hd_tex_line(X, r_, crowd, K, chans)
    D(c, j, n) = n == 1 ? "$(X)_{D$(c.tag)}" : "$(X)_{D$(c.tag)$j}"
    Hs = "$(X)_H"
    death = join(["-$(tn(c.δ))\\,$(D(c, c.n, c.n))" for c in chans if c.δ != 0])
    lines = ["\\frac{d$X}{dt}=$(tn(r_))\\,$Hs\\left(1-\\frac{$crowd}{$(tn(K))}\\right)$death"]
    for c in chans
        out = c.ρ == 0 ? tn(c.δ) : "($(tn(c.ρ))+$(tn(c.δ)))"
        s = "\\frac{d$(D(c, 1, c.n))}{dt}=$(c.k)\\,$Hs-$out\\,$(D(c, 1, c.n))"
        c.n == 2 && (s *= ",\\qquad \\frac{d$(D(c, 2, 2))}{dt}=$(tn(c.δ))\\left($(D(c, 1, 2))-$(D(c, 2, 2))\\right)")
        push!(lines, s)
    end
    push!(lines, "$Hs=$X-" * join([c.n == 1 ? D(c, 1, 1) : "$(D(c, 1, 2))-$(D(c, 2, 2))" for c in chans], "-"))
    lines
end
hill_tex(κ, h, D) = "$(tn(κ))\\,\\frac{$(D)^{$(tn(h))}}{1+$(D)^{$(tn(h))}}"

# damage that builds up after the drug is added: A(t) = 1 − e^{−λt}; λ = 0 means no build-up
hd_onset(λ, t) = λ > 0 ? 1 - exp(-λ*t) : one(t)
onset_tex(λ) = λ > 0 ? "\\left(1-e^{-$(tn(λ))\\,t}\\right)" : ""

# ---- Stage 1b: low resource damages healthy cells at rate σ·u_R (or σ·u_R·crowding) ----
hd1b = Dict{String,Any}()
for (line, sheet) in lr_mono_sheets
    s1 = stage1[line]; s = stage1b[line]; X = s1.sym
    ds = [(x=s.x[.!s.bad], y=s.y[.!s.bad], state_index=1)]
    # get(p) → (σ, ρ, δ); dens = damage rate scales with crowding X/K (Gyllenberg & Webb)
    mk(name, label, n, pn, p0, get, dens) = begin
        ns = [[n]]
        act = dens ? ((e, l, c, t, g) -> U_R*g) : ((e, l, c, t, g) -> U_R)
        (name=name, label=label, pnames=pn, p0=p0, bounds=[hdb() for _ in pn], datasets=ds,
         u0=hd_u0(1, ns, [[s.y[1]]]), obs=hd_obs(1, ns),
         model=make_hd(1, ns, p -> (v = get(p); (r=(s1.r,), K=(s1.K,), α=((1.0,),), ch=((v,),), act=act))),
         info=p -> (v = get(p); (n=n, σ=v[1], ρ=v[2], δ=v[3], dens=dens)),
         texf=p -> (v = get(p); hd_tex_line(X, s1.r, X, s1.K, [(tag="", n=n, ρ=v[2], δ=v[3],
                    k="$(tn(v[1]))\\cdot $(tn(U_R))" * (dens ? "\\,\\frac{$X}{$(tn(s1.K))}" : ""))])),
         extra=[("r_$X", s1.r, "inherited from Stage 1 (logistic)"), ("K_$X", s1.K, "inherited from Stage 1 (logistic)"),
                ("u_R", U_R, "fixed: (10 − 0.5) / (10 − 0.05) for 0.5% FBS")])
    end
    z3(a, b, c) = (a, b, c)
    cands = [
        mk("hd_lr_irreversible", "Healthy → damaged → death (irreversible starvation damage)", 1,
           ["σ_$X", "δ_$(X)^LR"], [0.3, 0.5], p -> z3(p[1], zero(p[1]), p[2]), false),
        mk("hd_lr_repair", "Healthy ⇄ damaged → death (damage can be repaired)", 1,
           ["σ_$X", "δ_$(X)^LR", "ρ_$(X)^LR"], [0.3, 0.5, 0.3], p -> z3(p[1], p[3], p[2]), false),
        mk("hd_lr_transit", "Healthy → damaged, 2 steps → death (delayed death)", 2,
           ["σ_$X", "δ_$(X)^LR"], [0.3, 1.0], p -> z3(p[1], zero(p[1]), p[2]), false),
        mk("hd_lr_transit4", "Healthy → damaged, 4 steps → death (longer delayed death)", 4,
           ["σ_$X", "δ_$(X)^LR"], [0.3, 2.0], p -> z3(p[1], zero(p[1]), p[2]), false),
        mk("hd_lr_arrest", "Healthy ⇄ arrested, no death (quiescence)", 1,
           ["σ_$X", "ρ_$(X)^LR"], [0.3, 0.3], p -> z3(p[1], p[2], zero(p[1])), false),
        mk("hd_lr_density_arrest", "Healthy ⇄ arrested, entry rises with crowding, no death (Gyllenberg–Webb)", 1,
           ["σ_$X", "ρ_$(X)^LR"], [1.0, 0.3], p -> z3(p[1], p[2], zero(p[1])), true),
        mk("hd_lr_density_death", "Healthy → damaged → death, damage rises with crowding", 1,
           ["σ_$X", "δ_$(X)^LR"], [1.0, 0.5], p -> z3(p[1], zero(p[1]), p[2]), true),
    ]
    rows = hd_rank([fit_hd(c) for c in cands])
    hd1b[line] = (rows=rows, best=rows[1])
    record!("1b low-resource monoculture (healthy/damaged)", "$line ($sheet)", rows,
            @sprintf("r_%s = %.4f, K_%s = %.1f from Stage 1 logistic; u_R = %.4f", X, s1.r, X, s1.K, U_R); line=line)
    for r in rows
        @printf("  %-9s %-22s BIC=%8.2f  %s\n", line, r.name, r.bic, param_str(r))
    end
    flush(stdout)
end

# ---- Stage 2: cisplatin damages healthy cells at rate κ·H(D)·A(t) ----
hd2 = Dict{String,Any}()
for (line, _, _) in stage1_lines
    s1 = stage1[line]; s = stage2[line]; X = s1.sym; doses = s.doses
    ds = [(x=x, y=y, state_index=i) for (i, (x, y)) in enumerate(s.data)]
    # get(p) → (κ, ρ, δ, h, λ); λ = 0: damage at full rate from day 0
    mk(name, label, n, pn, p0, bounds, get, usehill) = begin
        ns = [[n]]
        (name=name, label=label, pnames=pn, p0=p0, bounds=bounds, datasets=ds,
         u0=hd_u0(3, ns, [[y[1]] for (_, y) in s.data]), obs=hd_obs(3, ns),
         model=make_hd(3, ns, p -> (v = get(p); (r=(s1.r,), K=(s1.K,), α=((1.0,),), ch=(((v[1], v[2], v[3]),),),
                                               act=(e, l, c, t, g) -> (usehill ? hill(doses[e], v[4]) : doses[e]) * hd_onset(v[5], t)))),
         info=p -> (v = get(p); (n=n, hill=usehill, κ=v[1], ρ=v[2], δ=v[3], h=usehill ? v[4] : NaN, λ=v[5])),
         texf=p -> (v = get(p); hd_tex_line(X, s1.r, X, s1.K, [(tag="", n=n, ρ=v[2], δ=v[3],
                    k=(usehill ? hill_tex(v[1], v[4], "D") : "$(tn(v[1]))\\,D") * onset_tex(v[5]))])),
         extra=vcat([("r_$X", s1.r, "inherited from Stage 1 (logistic)"), ("K_$X", s1.K, "inherited from Stage 1 (logistic)")],
                    usehill ? [("EC50", EC50, "fixed")] : Tuple{String,Float64,String}[]))
    end
    bκ, bδ, bh, bρ, bλ = (0.0, 20.0), (0.0, 10.0), (0.1, 10.0), (0.0, 10.0), (1e-3, 10.0)
    Z(p) = zero(p[1])
    cands = [
        mk("hd_drug_irreversible", "Healthy → damaged → death, Hill dose response", 1,
           ["κ_$X", "δ_$X", "h_$X"], [0.5, 0.5, 2.0], [bκ, bδ, bh], p -> (p[1], Z(p), p[2], p[3], Z(p)), true),
        mk("hd_drug_repair", "Healthy ⇄ damaged → death (repair), Hill dose response", 1,
           ["κ_$X", "δ_$X", "h_$X", "ρ_$X"], [0.5, 0.5, 2.0, 0.2], [bκ, bδ, bh, bρ], p -> (p[1], p[4], p[2], p[3], Z(p)), true),
        mk("hd_drug_transit", "Healthy → damaged, 2 steps → death, Hill dose response", 2,
           ["κ_$X", "δ_$X", "h_$X"], [0.5, 1.0, 2.0], [bκ, bδ, bh], p -> (p[1], Z(p), p[2], p[3], Z(p)), true),
        mk("hd_drug_repair_linear", "Healthy ⇄ damaged → death (repair), damage ∝ dose", 1,
           ["κ_$X", "δ_$X", "ρ_$X"], [0.3, 0.5, 0.2], [bκ, bδ, bρ], p -> (p[1], p[3], p[2], Z(p), Z(p)), false),
        mk("hd_drug_onset", "Healthy → damaged → death, damage builds up over time (Hill)", 1,
           ["κ_$X", "δ_$X", "h_$X", "λ_$X"], [1.0, 0.5, 2.0, 0.3], [bκ, bδ, bh, bλ], p -> (p[1], Z(p), p[2], p[3], p[4]), true),
        mk("hd_drug_transit_onset", "Healthy → damaged, 2 steps → death, damage builds up over time (Hill)", 2,
           ["κ_$X", "δ_$X", "h_$X", "λ_$X"], [1.0, 1.0, 2.0, 0.3], [bκ, bδ, bh, bλ], p -> (p[1], Z(p), p[2], p[3], p[4]), true),
    ]
    rows = hd_rank([fit_hd(c) for c in cands])
    hd2[line] = (rows=rows, best=rows[1])
    record!("2 treated monoculture (healthy/damaged)", "$line, 0.62 / 1.0 / 1.47 uM jointly", rows,
            @sprintf("r_%s = %.4f, K_%s = %.1f from Stage 1 logistic; EC50 = 1 uM", X, s1.r, X, s1.K); line=line)
    for r in rows
        @printf("  %-9s %-22s BIC=%8.2f  %s\n", line, r.name, r.bic, param_str(r))
    end
    flush(stdout)
end

# ---- Stage 3: Ce0 + Ce1, shared K, drug damage on Ce1 inherited from Stage 2 ----
iN2, iC2 = hd2["A2780"].best.info, hd2["A2780cis"].best.info
hd_Hfac(i) = i.hill ? hill(COCULTURE_DOSE, i.h) : COCULTURE_DOSE
HfN, HfC = hd_Hfac(iN2), hd_Hfac(iC2)
ns3 = [[iN2.n], [iC2.n]]
u0_hd3 = hd_u0(6, ns3, [[d.S[2][1], d.R[2][1]] for d in co_data])
obs_hd3 = hd_obs(6, ns3)
# get(p) → (K, α_NC, α_CN, (κ_N, ρ_N, δ_N, b_N), (κ_C, ρ_C, δ_C, b_C)); b = damage rate present in
# every coculture, treated or not (the healthy/damaged counterpart of a baseline death rate)
function hd3_k(v, l, e, t)
    w, Hf, λ = l == 1 ? (v[4], HfN, iN2.λ) : (v[5], HfC, iC2.λ)
    w[4] + (ENV_DOSE[e] > 0 ? w[1]*Hf*hd_onset(λ, t) : zero(w[1]))
end
function hd3_ktex(w, Hf, λ)
    drug = "\\mathbb{1}_{\\text{Ce1}}\\,$(tn(w[1]))\\cdot $(tn(Hf))" * onset_tex(λ)
    w[4] > 1e-8 ? "\\left($(tn(w[4]))+$drug\\right)" : drug
end
function hd3_cand(name, label, pn, p0, bounds, get)
    (name=name, label=label, pnames=pn, p0=p0, bounds=bounds, datasets=co_datasets, u0=u0_hd3, obs=obs_hd3,
     model=make_hd(6, ns3, p -> (v = get(p); (r=(rN, rC), K=(v[1], v[1]), α=((1.0, v[2]), (v[3], 1.0)),
                                              ch=(((one(v[1]), v[4][2], v[4][3]),), ((one(v[1]), v[5][2], v[5][3]),)),
                                              act=(e, l, c, t, g) -> hd3_k(v, l, e, t)))),
     info=p -> (v = get(p); (K=v[1], aNC=v[2], aCN=v[3],
                             N=(n=iN2.n, κ=v[4][1], ρ=v[4][2], δ=v[4][3], b=v[4][4], Hf=HfN, λ=iN2.λ),
                             C=(n=iC2.n, κ=v[5][1], ρ=v[5][2], δ=v[5][3], b=v[5][4], Hf=HfC, λ=iC2.λ))),
     texf=p -> (v = get(p); vcat(
        hd_tex_line("N", rN, "N+$(tn(v[2]))\\,C", v[1], [(tag="", n=iN2.n, k=hd3_ktex(v[4], HfN, iN2.λ), ρ=v[4][2], δ=v[4][3])]),
        hd_tex_line("C", rC, "C+$(tn(v[3]))\\,N", v[1], [(tag="", n=iC2.n, k=hd3_ktex(v[5], HfC, iC2.λ), ρ=v[5][2], δ=v[5][3])]))),
     extra=[("r_N", rN, "inherited from Stage 1 (A2780 logistic)"), ("r_C", rC, "inherited from Stage 1 (A2780cis logistic)"),
            ("H_N(1 µM)", HfN, "from the Stage 2 healthy/damaged winner ($(hd2["A2780"].best.name))"),
            ("H_C(1 µM)", HfC, "from the Stage 2 healthy/damaged winner ($(hd2["A2780cis"].best.name))"),
            ("λ_N", iN2.λ, "damage build-up rate from Stage 2 (0 = none)"), ("λ_C", iC2.λ, "damage build-up rate from Stage 2 (0 = none)")])
end
Kα = (["K", "α_NC", "α_CN"], [1.1Kmax, 1.0, 1.0], [Kb, (0.0, 5.0), (0.0, 5.0)])
# w(κ, ρ, δ, b) per line; the inherited Stage 2 values fill whatever is not fitted
wN(p; κ=iN2.κ, ρ=iN2.ρ, δ=iN2.δ, b=zero(p[1])) = (κ, ρ, δ, b)
wC(p; κ=iC2.κ, ρ=iC2.ρ, δ=iC2.δ, b=zero(p[1])) = (κ, ρ, δ, b)
cat3(a) = (vcat(Kα[1], a[1]), vcat(Kα[2], a[2]), vcat(Kα[3], a[3]))
hd3_cands = [
    hd3_cand("hd_cc_inherited_damage", "Healthy/damaged: shared K + asymmetric crowding, damage kinetics from Stage 2",
        Kα..., p -> (p[1], p[2], p[3], wN(p), wC(p))),
    hd3_cand("hd_cc_kappaN", "Healthy/damaged: + sensitive damage rate refitted",
        cat3((["κ_N"], [iN2.κ], [(0.0, 20.0)]))..., p -> (p[1], p[2], p[3], wN(p; κ=p[4]), wC(p))),
    hd3_cand("hd_cc_kappa", "Healthy/damaged: + both damage rates refitted",
        cat3((["κ_N", "κ_C"], [iN2.κ, iC2.κ], [(0.0, 20.0), (0.0, 20.0)]))..., p -> (p[1], p[2], p[3], wN(p; κ=p[4]), wC(p; κ=p[5]))),
    hd3_cand("hd_cc_kappa_baseN", "Healthy/damaged: + both damage rates, and coculture damage of sensitive cells",
        cat3((["κ_N", "κ_C", "b_N"], [iN2.κ, iC2.κ, 0.1], [(0.0, 20.0), (0.0, 20.0), (0.0, 3.0)]))...,
        p -> (p[1], p[2], p[3], wN(p; κ=p[4], b=p[6]), wC(p; κ=p[5]))),
    hd3_cand("hd_cc_kappa_baseN_delta", "Healthy/damaged: + damage, coculture damage of sensitive cells, and death rates",
        cat3((["κ_N", "κ_C", "b_N", "δ_N", "δ_C"], [iN2.κ, iC2.κ, 0.1, iN2.δ, iC2.δ], [(0.0, 20.0), (0.0, 20.0), (0.0, 3.0), (0.0, 10.0), (0.0, 10.0)]))...,
        p -> (p[1], p[2], p[3], wN(p; κ=p[4], b=p[6], δ=p[7]), wC(p; κ=p[5], δ=p[8]))),
    hd3_cand("hd_cc_kappa_delta", "Healthy/damaged: + damage and death rates refitted",
        cat3((["κ_N", "κ_C", "δ_N", "δ_C"], [iN2.κ, iC2.κ, iN2.δ, iC2.δ], [(0.0, 20.0), (0.0, 20.0), (0.0, 10.0), (0.0, 10.0)]))...,
        p -> (p[1], p[2], p[3], wN(p; κ=p[4], δ=p[6]), wC(p; κ=p[5], δ=p[7]))),
    hd3_cand("hd_cc_full", "Healthy/damaged: + damage, repair and death rates refitted",
        cat3((["κ_N", "κ_C", "δ_N", "δ_C", "ρ_N", "ρ_C"], [iN2.κ, iC2.κ, iN2.δ, iC2.δ, max(iN2.ρ, 0.05), max(iC2.ρ, 0.05)],
              [(0.0, 20.0), (0.0, 20.0), (0.0, 10.0), (0.0, 10.0), (0.0, 10.0), (0.0, 10.0)]))...,
        p -> (p[1], p[2], p[3], wN(p; κ=p[4], δ=p[6], ρ=p[8]), wC(p; κ=p[5], δ=p[7], ρ=p[9]))),
]
hd3_rows = hd_rank([fit_hd(c) for c in hd3_cands])
hd3_best = hd3_rows[1]
hd3_inherited = @sprintf("r_N = %.4f, r_C = %.4f (Stage 1); damage structure and non-refitted κ, ρ, δ, λ from the Stage 2 healthy/damaged winners (%s, %s); H_N(1 µM) = %.3f, H_C(1 µM) = %.3f",
    rN, rC, hd2["A2780"].best.name, hd2["A2780cis"].best.name, HfN, HfC)
record!("3 Ce0 + Ce1 coculture (joint, healthy/damaged)", "3 ratios × untreated + treated, 12 trajectories", hd3_rows, hd3_inherited)
for r in hd3_rows
    @printf("  %-24s k=%d BIC=%8.2f  %s\n", r.name, length(r.params), r.bic, param_str(r))
end
flush(stdout)

# ---- Stage 3b: low-resource cocultures, starvation damage inherited from Stage 1b ----
iN1b, iC1b = hd1b["A2780"].best.info, hd1b["A2780cis"].best.info
i3 = hd3_best.info
ns3b = [[iN1b.n], [iC1b.n]]
u0_hd3b = hd_u0(6, ns3b, [[c.yN[1], c.yC[1]] for c in lr_cc])
obs_hd3b = hd_obs(6, ns3b)
sN1b = (iN1b.σ, iN1b.ρ, iN1b.δ); sC1b = (iC1b.σ, iC1b.ρ, iC1b.δ)
lr_ktex(σ, dens, X) = "$(tn(σ))\\cdot $(tn(U_R))" * (dens ? "\\,g_$X" : "")
function hd3b_cand(name, label, pn, p0, bounds, get)
    (name=name, label=label, pnames=pn, p0=p0, bounds=bounds, datasets=lr_cc_datasets, u0=u0_hd3b, obs=obs_hd3b,
     model=make_hd(6, ns3b, p -> (v = get(p); (r=(rN, rC), K=(v[1], v[1]), α=((1.0, v[2]), (v[3], 1.0)), ch=((v[4],), (v[5],)),
                                                act=(e, l, c, t, g) -> U_R * ((l == 1 ? iN1b.dens : iC1b.dens) ? g : one(g))))),
     info=p -> (v = get(p); (K=v[1], aNC=v[2], aCN=v[3], N=(n=iN1b.n, σ=v[4][1], ρ=v[4][2], δ=v[4][3], dens=iN1b.dens),
                                                         C=(n=iC1b.n, σ=v[5][1], ρ=v[5][2], δ=v[5][3], dens=iC1b.dens))),
     texf=p -> (v = get(p); vcat(
        hd_tex_line("N", rN, "N+$(tn(v[2]))\\,C", v[1], [(tag="", n=iN1b.n, k=lr_ktex(v[4][1], iN1b.dens, "N"), ρ=v[4][2], δ=v[4][3])]),
        hd_tex_line("C", rC, "C+$(tn(v[3]))\\,N", v[1], [(tag="", n=iC1b.n, k=lr_ktex(v[5][1], iC1b.dens, "C"), ρ=v[5][2], δ=v[5][3])]),
        iN1b.dens || iC1b.dens ? ["g_N=\\frac{N+$(tn(v[2]))\\,C}{$(tn(v[1]))},\\quad g_C=\\frac{C+$(tn(v[3]))\\,N}{$(tn(v[1]))}"] : String[])),
     extra=[("r_N", rN, "inherited from Stage 1 (A2780 logistic)"), ("r_C", rC, "inherited from Stage 1 (A2780cis logistic)"),
            ("u_R", U_R, "fixed: (10 − 0.5) / (10 − 0.05) for 0.5% FBS")])
end
Kα3 = (["K", "α_NC", "α_CN"], [i3.K, max(i3.aNC, 0.05), max(i3.aCN, 0.05)], [(0.5lrmax, 20lrmax), (0.0, 5.0), (0.0, 5.0)])
hd3b_rows = hd_rank([fit_hd(c) for c in [
    hd3b_cand("hd_lrcc_inherited", "Healthy/damaged: everything inherited (K, α from Stage 3, damage from Stage 1b)",
        String[], Float64[], Tuple{Float64,Float64}[], p -> (i3.K, i3.aNC, i3.aCN, sN1b, sC1b)),
    hd3b_cand("hd_lrcc_K_alpha", "Healthy/damaged: K and α refitted, starvation damage from Stage 1b",
        Kα3..., p -> (p[1], p[2], p[3], sN1b, sC1b)),
    hd3b_cand("hd_lrcc_sigma", "Healthy/damaged: starvation damage rates refitted, K and α from Stage 3",
        ["σ_N", "σ_C"], [max(sN1b[1], 0.05), max(sC1b[1], 0.05)], [(0.0, 10.0), (0.0, 10.0)],
        p -> (i3.K, i3.aNC, i3.aCN, (p[1], sN1b[2], sN1b[3]), (p[2], sC1b[2], sC1b[3]))),
    hd3b_cand("hd_lrcc_K_alpha_sigma", "Healthy/damaged: K, α and starvation damage rates refitted",
        vcat(Kα3[1], "σ_N", "σ_C"), vcat(Kα3[2], max(sN1b[1], 0.05), max(sC1b[1], 0.05)), vcat(Kα3[3], [(0.0, 10.0), (0.0, 10.0)]),
        p -> (p[1], p[2], p[3], (p[4], sN1b[2], sN1b[3]), (p[5], sC1b[2], sC1b[3]))),
    hd3b_cand("hd_lrcc_full", "Healthy/damaged: K, α, starvation damage and death rates refitted",
        vcat(Kα3[1], "σ_N", "σ_C", "δ_N^LR", "δ_C^LR"), vcat(Kα3[2], max(sN1b[1], 0.05), max(sC1b[1], 0.05), max(sN1b[3], 0.05), max(sC1b[3], 0.05)),
        vcat(Kα3[3], [(0.0, 10.0), (0.0, 10.0), (0.0, 10.0), (0.0, 10.0)]),
        p -> (p[1], p[2], p[3], (p[4], sN1b[2], p[6]), (p[5], sC1b[2], p[7]))),
]])
hd3b_best = hd3b_rows[1]
hd3b_inherited = @sprintf("r_N = %.4f, r_C = %.4f (Stage 1); starvation-damage structure and non-refitted σ, ρ, δ from the Stage 1b healthy/damaged winners (%s, %s); K, α from %s (Stage 3); u_R = %.4f",
    rN, rC, hd1b["A2780"].best.name, hd1b["A2780cis"].best.name, hd3_best.name, U_R)
record!("3b low-resource coculture (healthy/damaged)", "3 ratios × runs 1 and 2 (30k), 12 trajectories", hd3b_rows, hd3b_inherited)
for r in hd3b_rows
    @printf("  %-22s k=%d BIC=%8.2f  %s\n", r.name, length(r.params), r.bic, param_str(r))
end
flush(stdout)

# Fully resolved healthy/damaged parameters of each stage winner (for the validation script)
hd_resolved = NamedTuple[]
hdput!(stage, line, model, kv) = for (k, v) in pairs(kv); push!(hd_resolved, (stage=stage, line=line, model=model, param=String(k), value=Float64(v))); end
for (line, _, _) in stage1_lines
    b = hd1b[line].best; hdput!("1b", line, b.name, b.info)      # dens stored as 1/0
    b = hd2[line].best;  hdput!("2", line, b.name, b.info)       # hill stored as 1/0
end
for (st, b) in (("3", hd3_best), ("3b", hd3b_best))
    hdput!(st, "coculture", b.name, (K=b.info.K, aNC=b.info.aNC, aCN=b.info.aCN))
    hdput!(st, "A2780", b.name, b.info.N); hdput!(st, "A2780cis", b.name, b.info.C)
end
CSV.write(joinpath(OUT_DIR, "hd_resolved_parameters.csv"), DataFrame(hd_resolved))

# ── CSV outputs ────────────────────────────────────────────────────────────────

for stale in ["stage1_baseline.csv", "stage1_predictions.csv", "stage2_results.csv",
              "stage2_predictions.csv", "stage3_params.csv", "stage3_predictions.csv"]
    rm(joinpath(OUT_DIR, stale); force=true)
end
CSV.write(joinpath(OUT_DIR, "bic_model_ranking.csv"), DataFrame(ranking_rows))

carried = NamedTuple[]
for (line, _, _) in stage1_lines
    s1 = stage1[line]; s2 = stage2[line]
    push!(carried, (stage="1", line=line, param="r_$(s1.sym)", value=s1.r, source="logistic (Stage 1)"))
    push!(carried, (stage="1", line=line, param="K_$(s1.sym)", value=s1.K, source="logistic (Stage 1)"))
    push!(carried, (stage="2", line=line, param="λ_$(s2.sym)", value=s2.λ, source="$(s2.best.name) (Stage 2)"))
    isnan(s2.h) || push!(carried, (stage="2", line=line, param="h_$(s2.sym)", value=s2.h, source="$(s2.best.name) (Stage 2)"))
end
for (line, _, _) in stage1_lines
    s = stage1b[line]
    push!(carried, (stage="1b", line=line, param="d_$(s.sym)^LR", value=s.dLR, source="lr_subtractive_loss ★ (Stage 1b)"))
end
for (n, v) in zip(s3_best.pnames, s3_best.params)
    push!(carried, (stage="3", line="coculture", param=n, value=v, source="$(s3_best.name) (Stage 3)"))
end
for (n, v) in zip(s3b_best.pnames, s3b_best.params)
    push!(carried, (stage="3b", line="LR coculture", param=n, value=v, source="$(s3b_best.name) (Stage 3b)"))
end
CSV.write(joinpath(OUT_DIR, "carried_parameters.csv"), DataFrame(carried))
CSV.write(joinpath(OUT_DIR, "lr_excluded_points.csv"), DataFrame(outlier_rows))
CSV.write(joinpath(OUT_DIR, "all_fitted_parameters.csv"), DataFrame(all_param_rows))

# ── report ─────────────────────────────────────────────────────────────────────

esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

all_y = vcat([s.y for s in values(stage1)]..., [y for s in values(stage2) for (_, y) in s.data]..., [d.y for d in co_datasets]...,
             [s.y for s in values(stage1b)]..., [c.yN for c in lr_cc]..., [c.yC for c in lr_cc]...)
all_x = vcat([s.x for s in values(stage1)]..., [x for s in values(stage2) for (x, _) in s.data]..., [d.x for d in co_datasets]...,
             [s.x for s in values(stage1b)]..., [c.xN for c in lr_cc]...)
const XMAX = maximum(all_x)
const YMAX = ceil(maximum(all_y) * 1.10; sigdigits=2)

const BLUE, ORANGE, GREEN, PURPLE, GREY = "#0b6e99", "#d95f02", "#3b8b4f", "#8b5cf6", "#374151"
const CAND_COLORS = [BLUE, GREEN, PURPLE, ORANGE, "#be185d", "#0f766e", "#b45309", "#4d7c0f", "#7c3aed", "#0369a1", "#9f1239", "#374151", "#a16207"]
const CAND_DASHES = ["", "7,4", "2,3", "10,3,2,3", "1,4", "12,4", "4,2", "8,2,2,2", "3,6", "14,3", "6,3,1,3", "2,2", "9,5"]

# healthy/damaged rows carry their own state layout, equations and inherited values
oix(r, k) = hasproperty(r, :obs) ? r.obs[k] : k
ru0(r, u0) = hasproperty(r, :u0) ? r.u0 : u0
family(r) = hasproperty(r, :family) ? r.family : "One population per line"
all_rows(a, b) = sort(vcat(a, b); by=r -> r.bic)

tn(v) = abs(v) < 1e-8 ? "0" : fmt(v)
shown(v) = abs(v) < 1e-8 ? "≈ 0 ($(@sprintf("%.1e", v)))" : fmt(v)
const UNITS = Dict("r" => "per day", "K" => "cells (Mean Cells units)", "θ" => "dimensionless", "β" => "per day",
    "λ" => "per day", "h" => "dimensionless", "α" => "dimensionless", "d" => "per day", "EC50" => "µM",
    "H" => "dimensionless (0–1)", "D" => "µM", "ρ" => "per day (repair)", "u" => "dimensionless (0–1)",
    "κ" => "per day (damage rate at full effect)", "b" => "per day (coculture damage rate)", "σ" => "per day (starvation damage rate)", "δ" => "per day (damaged-cell exit)")
unit_of(n) = startswith(n, "ρ_r") || startswith(n, "ρ_K") ? "dimensionless multiplier" : get(UNITS, String(strip(first(split(n, r"[_( ]")))), "")
evidence(d) = d == 0 ? "BIC winner" : d < 2 ? "weak (ΔBIC < 2)" : d < 6 ? "positive (2–6)" : d < 10 ? "strong (6–10)" : "very strong (> 10)"

# obs: (name, x, y, color, filled); fits: (t, y, color, dash); lines: legend entries (label, color, dash)
# marks: (x, y, kind = :excluded | :replacement, color)
function svg_panel(title, obs, fits; lines=NamedTuple[], marks=NamedTuple[], W=760)
    row2 = !isempty(lines) || !isempty(marks)
    H = row2 ? 352 : 330
    L, R, T, B = 64, 20, 30, row2 ? 84 : 62
    sx(v) = L + v/XMAX*(W-L-R)
    sy(v) = H - B - v/YMAX*(H-T-B)
    p = String["<svg viewBox='0 0 $W $H' role='img' aria-label='$(esc(title))'><rect width='100%' height='100%' fill='white'/>",
               "<text x='$L' y='18' font-size='14' font-weight='bold'>$(esc(title))</text>"]
    for v in range(0, YMAX; length=5)
        push!(p, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#e5e9ec'/><text x='%d' y='%.1f' font-size='10' text-anchor='end'>%.0f</text>", L, W-R, sy(v), sy(v), L-6, sy(v)+3, v))
    end
    for v in 0:2:floor(Int, XMAX)
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%d</text>", sx(v), H-B+14, v))
    end
    push!(p, "<line x1='$L' y1='$(H-B)' x2='$(W-R)' y2='$(H-B)' stroke='#555'/><line x1='$L' y1='$T' x2='$L' y2='$(H-B)' stroke='#555'/>")
    push!(p, "<text x='$((L+W-R)÷2)' y='$(H-B+30)' font-size='11' text-anchor='middle'>Day</text><text transform='translate(14 $((T+H-B)÷2)) rotate(-90)' font-size='11' text-anchor='middle'>Mean Cells</text>")
    for f in fits
        pts = join([@sprintf("%.1f,%.1f", sx(t), sy(clamp(y, 0, YMAX))) for (t, y) in zip(f.t, f.y)], " ")
        push!(p, "<polyline points='$pts' fill='none' stroke='$(f.color)' stroke-width='2'$(isempty(f.dash) ? "" : " stroke-dasharray='$(f.dash)'")/>")
    end
    oy = row2 ? H - 34 : H - 12
    for (i, o) in enumerate(obs)
        fill = o.filled ? o.color : "white"
        push!(p, join([@sprintf("<circle cx='%.1f' cy='%.1f' r='3.2' fill='%s' stroke='%s' stroke-width='1.4'/>", sx(x), sy(y), fill, o.color) for (x, y) in zip(o.x, o.y)]))
        lx = L + (i-1)*170
        push!(p, "<circle cx='$(lx+4)' cy='$(oy-4)' r='4' fill='$fill' stroke='$(o.color)' stroke-width='1.4'/><text x='$(lx+12)' y='$oy' font-size='11' fill='$(o.color)'>$(esc(o.name))</text>")
    end
    xmark(x, y) = @sprintf("<path d='M%.1f %.1fL%.1f %.1fM%.1f %.1fL%.1f %.1f' stroke='#9ca3af' stroke-width='1.8'/>", x-4, y-4, x+4, y+4, x-4, y+4, x+4, y-4)
    diamond(x, y, c) = @sprintf("<polygon points='%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f' fill='white' stroke='%s' stroke-width='1.6'/>", x, y-5, x+5, y, x, y+5, x-5, y, c)
    for m in marks
        push!(p, m.kind == :excluded ? xmark(sx(m.x), sy(clamp(m.y, 0, YMAX))) : diamond(sx(m.x), sy(clamp(m.y, 0, YMAX)), m.color))
    end
    lx = L
    for l in lines
        push!(p, "<line x1='$lx' x2='$(lx+26)' y1='$(H-16)' y2='$(H-16)' stroke='$(l.color)' stroke-width='2.2'$(isempty(l.dash) ? "" : " stroke-dasharray='$(l.dash)'")/><text x='$(lx+32)' y='$(H-12)' font-size='11'>$(esc(l.label))</text>")
        lx += 235
    end
    if any(m -> m.kind == :excluded, marks)
        push!(p, xmark(lx+6, H-16) * "<text x='$(lx+16)' y='$(H-12)' font-size='11'>excluded outlier</text>"); lx += 130
    end
    if any(m -> m.kind == :replacement, marks)
        push!(p, diamond(lx+6, H-16, GREY) * "<text x='$(lx+16)' y='$(H-12)' font-size='11'>replacement on fitted line</text>")
    end
    push!(p, "</svg>")
    join(p)
end

function bic_svg(rows; W=760)
    H = 40 + 50*length(rows)
    base = rows[1].bic; span = max(rows[end].bic - base, 1.0)
    p = ["<svg viewBox='0 0 $W $H' class='bic' role='img' aria-label='BIC comparison'><text x='8' y='18' font-size='14' font-weight='bold'>BIC by candidate (lower is better; bar length grows with ΔBIC)</text>"]
    for (i, r) in enumerate(rows)
        y = 32 + (i-1)*50; d = r.bic - base; w = 60 + 440*d/span
        push!(p, @sprintf("<text x='8' y='%d' font-size='13'%s>%s</text><rect x='8' y='%d' width='%.1f' height='22' rx='4' fill='%s'/><text x='%.1f' y='%d' font-size='13'>BIC %.2f · ΔBIC %.2f · %s</text>",
            y+12, i == 1 ? " font-weight='bold'" : "", esc(r.label), y+18, w, i == 1 ? "#176b87" : "#9fb7c3", 16+w, y+34, r.bic, d, evidence(d)))
    end
    push!(p, "</svg>")
    join(p)
end

function comparison_table(rows)
    body = join(["<tr$(i == 1 ? " class='win'" : "")><td>$(esc(r.label))</td><td>$(esc(family(r)))</td><td class='num'>$(length(r.params))</td><td class='num'>$(@sprintf("%.2f", r.bic))</td><td class='num'>$(@sprintf("%.2f", r.bic - rows[1].bic))</td><td>$(evidence(r.bic - rows[1].bic))</td><td class='num'>$(@sprintf("%.4g", r.sse))</td><td>$(esc(join(["$(n) = $(shown(v))" for (n, v) in zip(r.pnames, r.params)], ", ")))</td></tr>" for (i, r) in enumerate(rows)])
    "<div class='tscroll'><table class='cmp-table'><tr><th>Candidate</th><th>Model family</th><th>Free params</th><th>BIC</th><th>ΔBIC</th><th>Evidence against vs winner</th><th>SSE</th><th>Fitted parameters</th></tr>$body</table></div>"
end

function param_table(r, extra)
    at_bound(v, (lo, hi)) = isapprox(v, hi; rtol=1e-3) ? " — at upper bound $(fmt(hi))" :
                            isapprox(v, lo; atol=1e-6, rtol=1e-3) ? " — at lower bound $(fmt(lo))" : ""
    se = isempty(r.se) ? fill(NaN, length(r.params)) : r.se
    rows = vcat([(n, v, "fitted" * at_bound(v, b), s) for (n, v, b, s) in zip(r.pnames, r.params, r.bounds, se)],
                [(n, v, s, NaN) for (n, v, s) in extra])
    ci(v, s) = isnan(s) ? "—" : "$(fmt(v - 1.96s)) to $(fmt(v + 1.96s))"
    body = join(["<tr$(startswith(st, "fitted") ? " class='fitted'" : "")><td>$(esc(n))</td><td class='num'>$(esc(shown(v)))</td><td class='num'>$(isnan(s) ? "—" : fmt(s))</td><td class='num'>$(ci(v, s))</td><td>$(esc(unit_of(n)))</td><td>$(esc(st))</td></tr>" for (n, v, st, s) in rows])
    "<div class='tscroll'><table class='params'><tr><th>Parameter</th><th>Value</th><th>SE</th><th>95% interval</th><th>Units</th><th>Status</th></tr>$body</table></div>"
end

const HD_NOTE = "<p class='muted'>\\(X_H\\) = healthy (dividing) cells, \\(X_D\\) = damaged cells: they no longer divide but still take up space and are still counted, so the plotted count is \\(X = X_H + X_D\\). Damaged cells die at rate δ when they leave the last damaged compartment, or return to the healthy pool at repair rate ρ.</p>"
function fit_box(r, rows, tex_lines; note="", extra=Tuple{String,Float64,String}[])
    if hasproperty(r, :texf)
        tex_lines, extra, note = r.texf(r.params), r.extra, HD_NOTE
    end
    d = r.bic - rows[1].bic
    "<div class='fitbox'><h4>$(esc(r.label))$(d == 0 ? " <span class='badge'>BIC winner</span>" : "")</h4>" *
    "<p class='muted'>BIC $(@sprintf("%.2f", r.bic)) · ΔBIC $(@sprintf("%.2f", d)) · SSE $(@sprintf("%.4g", r.sse)) · $(length(r.params)) free parameters</p>" *
    "<h5>Fitted equation</h5>" * join(["<div class='eq'>\\[$t\\]</div>" for t in tex_lines]) * note *
    "<h5>Parameter values</h5>" * param_table(r, extra) * "</div>"
end

function compare_block(id, rows, views; overlay=nothing)
    btns = join(["<button type='button' class='cbtn$(i == 1 ? " active" : "")' data-view='$id-$i'>$(esc(r.label))<span class='d'>$(i == 1 ? "BIC winner" : @sprintf("ΔBIC %.1f", r.bic - rows[1].bic))</span></button>" for (i, r) in enumerate(rows)])
    overlay === nothing || (btns *= "<button type='button' class='cbtn' data-view='$id-all'>All candidates overlaid<span class='d'>visual comparison</span></button>")
    vs = join(["<div class='cview$(i == 1 ? " active" : "")' id='$id-$i'>$(views[i])</div>" for i in eachindex(rows)])
    overlay === nothing || (vs *= "<div class='cview' id='$id-all'>$overlay</div>")
    "<div class='cmp'><p class='switch-label'>Show fit for:</p><div class='cbtns'>$btns</div>$vs</div>"
end

vgrid(plots, box) = "<div class='vgrid'><div class='plots'>$plots</div>$box</div>"

# fitted-equation TeX, with values substituted
function tex_stage1(name, X, p)
    name == "logistic" && return "\\frac{d$X}{dt}=$(tn(p[1]))\\,$X\\left(1-\\frac{$X}{$(tn(p[2]))}\\right)"
    name == "theta_logistic" && return "\\frac{d$X}{dt}=$(tn(p[1]))\\,$X\\left[1-\\left(\\frac{$X}{$(tn(p[2]))}\\right)^{$(tn(p[3]))}\\right]"
    "\\frac{d$X}{dt}=$(tn(p[1]))\\,$X\\,\\log\\left(\\frac{$(tn(p[2]))}{$X}\\right)"
end

function tex_stage2(name, X, p, rX, KX)
    growth = "$(tn(rX))\\,$X\\left(1-\\frac{$X}{$(tn(KX))}\\right)"
    delay  = "\\left(1-e^{-$(tn(p[2]))\\,t}\\right)"
    kill   = name == "delayed_hill_kill" ? "\\frac{D^{$(tn(p[3]))}}{1^{$(tn(p[3]))}+D^{$(tn(p[3]))}}" : "D"
    "\\frac{d$X}{dt}=$growth-$(tn(p[1]))\\,$delay\\,$kill\\,$X"
end

const S3_INHERITED_KILL = ("shared_K_only", "crowding_inherited_kill")

function tex_stage3(r)
    v = Dict(zip(r.pnames, r.params))
    inh = r.name in S3_INHERITED_KILL
    line(X, Y, a, r_, λ, H, βinh) = begin
        crowd = haskey(v, "α_$a") ? "$X+$(tn(v["α_$a"]))\\,$Y" : "$X+$Y"
        s = "\\frac{d$X}{dt}=$(tn(r_))\\,$X\\left(1-\\frac{$crowd}{$(tn(v["K"]))}\\right)"
        haskey(v, "d_$X") && (s *= "-$(tn(v["d_$X"]))\\,$X")
        β = inh ? βinh : get(v, "β_$X", nothing)
        β === nothing || (s *= "-\\mathbb{1}_{\\text{Ce1}}\\,$(tn(β))\\left(1-e^{-$(tn(λ))\\,t}\\right)($(tn(H)))\\,$X")
        s
    end
    [line("N", "C", "NC", rN, λN, HN, βN2), line("C", "N", "CN", rC, λC, HC, βC2)]
end

function s3_extra_rows(r)
    rows = Tuple{String,Float64,String}[]
    if r.name in S3_INHERITED_KILL
        push!(rows, ("β_N", βN2, "inherited from Stage 2"), ("β_C", βC2, "inherited from Stage 2"))
    else
        for n in ("β_N", "β_C", "d_N", "d_C")
            n in r.pnames || push!(rows, (n, 0.0, "fixed at 0 (term not in model)"))
        end
    end
    "α_NC" in r.pnames || push!(rows, ("α_NC = α_CN", 1.0, "fixed (one-for-one crowding)"))
    rows
end

# Stage 1 panels
s1_html = String[]
for (line, _, _) in stage1_lines
    s = stage1[line]; X = s.sym
    obs = [(name="$line observed", x=s.x, y=s.y, color=GREY, filled=true)]
    curves = [begin ts, sol = dense(r.model, [s.y[1]], r.params, XMAX); (t=ts, y=[u[1] for u in sol.u]) end for r in s.rows]
    views = [vgrid(svg_panel("$line untreated — $(r.label)", obs, [(t=c.t, y=c.y, color=BLUE, dash="")]),
                   fit_box(r, s.rows, [tex_stage1(r.name, X, r.params)];
                           note = r.name == "logistic" ? "<p class='muted'>These logistic r and K are carried into Stage 2 (★ on the Candidate models page).</p>" : ""))
             for (r, c) in zip(s.rows, curves)]
    overlay = svg_panel("$line untreated — all candidates", obs,
        [(t=c.t, y=c.y, color=CAND_COLORS[i], dash=CAND_DASHES[i]) for (i, c) in enumerate(curves)];
        lines=[(label=r.label, color=CAND_COLORS[i], dash=CAND_DASHES[i]) for (i, r) in enumerate(s.rows)])
    push!(s1_html, "<div class='panel'><h3>$line <span class='muted'>($(esc(s.file)))</span></h3>" *
        bic_svg(s.rows) * comparison_table(s.rows) * compare_block("s1-$(lowercase(line))", s.rows, views; overlay=overlay) * "</div>")
end

# Stage 2 panels
s2_html = String[]
dose_colors = [PURPLE, GREEN, ORANGE]
for (line, _, _) in stage1_lines
    s = stage2[line]; s1 = stage1[line]; X = s.sym
    rows2 = all_rows(s.rows, hd2[line].rows)
    obs = [(name="$(d) µM observed", x=x, y=y, color=dose_colors[i], filled=true) for (i, (d, (x, y))) in enumerate(zip(s.doses, s.data))]
    sols = [dense(r.model, ru0(r, s.u0), r.params, XMAX) for r in rows2]
    views = String[]
    for (r, (ts, sol)) in zip(rows2, sols)
        fits = [(t=ts, y=[u[oix(r, i)] for u in sol.u], color=dose_colors[i], dash="") for i in eachindex(s.doses)]
        hill_note = r.name == "delayed_hill_kill" ?
            "<p class='muted'>Hill factor at each dose: " * join([@sprintf("H(%.2f µM) = %.3f", d, hill(d, r.params[3])) for d in s.doses], ", ") * ".</p>" : ""
        extra = [("r_$X", s1.r, "inherited from Stage 1 (logistic)"), ("K_$X", s1.K, "inherited from Stage 1 (logistic)")]
        r.name == "delayed_hill_kill" && push!(extra, ("EC50", EC50, "fixed"))
        push!(views, vgrid(svg_panel("$line treated — $(r.label)", obs, fits),
                           fit_box(r, rows2, [tex_stage2(r.name, X, r.params, s1.r, s1.K)]; note=hill_note, extra=extra)))
    end
    overlay = svg_panel("$line treated — all candidates (colour = dose, line style = model)", obs,
        [(t=ts, y=[u[oix(rows2[k], i)] for u in sol.u], color=dose_colors[i], dash=CAND_DASHES[k]) for (k, (ts, sol)) in enumerate(sols) for i in eachindex(s.doses)];
        lines=[(label=r.label, color=GREY, dash=CAND_DASHES[k]) for (k, r) in enumerate(rows2)])
    push!(s2_html, "<div class='panel'><h3>$line — three doses fitted jointly</h3><p class='muted'>Files: $(esc(join(s.files, ", ")))</p>" *
        bic_svg(rows2) * comparison_table(rows2) * compare_block("s2-$(lowercase(line))", rows2, views; overlay=overlay) * "</div>")
end

# Stage 3: one comparison block; each candidate view shows all three ratios
s3_extra = [("r_N", rN, "inherited from Stage 1 (A2780 logistic)"), ("r_C", rC, "inherited from Stage 1 (A2780cis logistic)"),
            ("λ_N", λN, "inherited from Stage 2 ($(stage2["A2780"].best.name))"), ("λ_C", λC, "inherited from Stage 2 ($(stage2["A2780cis"].best.name))"),
            ("H_N(1 µM)", HN, stage2["A2780"].best.name == "delayed_hill_kill" ? "derived: Hill factor at D = EC50" : "no Hill saturation: equals D"),
            ("H_C(1 µM)", HC, stage2["A2780cis"].best.name == "delayed_hill_kill" ? "derived: Hill factor at D = EC50" : "no Hill saturation: equals D"),
            ("D (Ce1)", COCULTURE_DOSE, "fixed"), ("EC50", EC50, "fixed")]
s3_all = all_rows(s3_rows, hd3_rows)
s3_views = String[]
for r in s3_all
    ts3, sol3 = dense(r.model, ru0(r, u0_s3), r.params, XMAX)
    plots = String[]
    for ratio in ["25S75R", "50S50R", "75S25R"]
        i0 = findfirst(e -> e.ratio == ratio && e.cond == "Ce0", coculture_layout)
        i1 = findfirst(e -> e.ratio == ratio && e.cond == "Ce1", coculture_layout)
        d0, d1 = co_data[i0], co_data[i1]
        obs = [(name="Ce0 sensitive", x=d0.S[1], y=d0.S[2], color=BLUE, filled=true),
               (name="Ce0 resistant", x=d0.R[1], y=d0.R[2], color=ORANGE, filled=true),
               (name="Ce1 sensitive", x=d1.S[1], y=d1.S[2], color=BLUE, filled=false),
               (name="Ce1 resistant", x=d1.R[1], y=d1.R[2], color=ORANGE, filled=false)]
        fits = [(t=ts3, y=[u[oix(r, 2i0-1)] for u in sol3.u], color=BLUE,   dash=""),
                (t=ts3, y=[u[oix(r, 2i0)]   for u in sol3.u], color=ORANGE, dash=""),
                (t=ts3, y=[u[oix(r, 2i1-1)] for u in sol3.u], color=BLUE,   dash="7,4"),
                (t=ts3, y=[u[oix(r, 2i1)]   for u in sol3.u], color=ORANGE, dash="7,4")]
        push!(plots, svg_panel("$ratio — Ce0 solid / filled, Ce1 (1 µM) dashed / open", obs, fits))
    end
    push!(s3_views, vgrid(join(plots), fit_box(r, s3_all, tex_stage3(r);
        note="<p class='muted'>\\(\\mathbb{1}_{\\text{Ce1}}\\) = 1 for treated (Ce1) cultures and 0 for untreated (Ce0). K is shared by both lineages, all ratios and both conditions.</p>",
        extra=vcat(s3_extra, s3_extra_rows(r)))))
end
s3_html = "<div class='panel'><h3>All 3 ratios × Ce0 + Ce1 (12 trajectories, one joint fit)</h3>" *
    bic_svg(s3_all) * comparison_table(s3_all) * compare_block("s3", s3_all, s3_views) * "</div>"

# excluded points and their on-the-line replacements for one trajectory
lr_marks(x, y, bad, sol, idx, color) = vcat(
    [(x=x[i], y=y[i], kind=:excluded, color=GREY) for i in findall(bad)],
    [(x=x[i], y=sol(x[i])[idx], kind=:replacement, color=color) for i in findall(bad)])

function tex_stage1b(name, X, p, rX, KX)
    name == "lr_subtractive_loss" && return "\\frac{d$X}{dt}=$(tn(rX))\\,$X\\left(1-\\frac{$X}{$(tn(KX))}\\right)-$(tn(p[1]))\\cdot $(tn(U_R))\\,$X"
    name == "lr_rho_modifiers" && return "\\frac{d$X}{dt}=$(tn(p[1]))\\cdot $(tn(rX))\\,$X\\left(1-\\frac{$X}{$(tn(p[2]))\\cdot $(tn(KX))}\\right)"
    for n in (4, 8)
        name == "lr_loss_plus_accum$n" && return "\\frac{d$X}{dt}=$(tn(rX))\\,$X\\left(1-\\frac{$X}{$(tn(KX))}\\right)-\\left($(tn(p[1]))+$(tn(p[2]))\\left[1-e^{-x}\\sum_{k=0}^{$(n-1)}\\frac{x^k}{k!}\\right]\\right)$(tn(U_R))\\,$X,\\quad x=$(tn(p[3]))\\,t"
    end
    name == "lr_accumulated_loss4" && return "\\frac{d$X}{dt}=$(tn(rX))\\,$X\\left(1-\\frac{$X}{$(tn(KX))}\\right)-$(tn(p[1]))\\cdot $(tn(U_R))\\left[1-e^{-x}\\sum_{k=0}^{3}\\frac{x^k}{k!}\\right]$X,\\quad x=$(tn(p[2]))\\,t"
    "\\frac{d$X}{dt}=$(tn(p[1]))\\,$X\\left(1-\\frac{$X}{$(tn(KX))}\\right)"
end

function tex_stage3b(r)
    K, aNC, aCN, dN, dC = r.name == "lr_cc_refit_all" ? Tuple(r.params) :
                          r.name == "lr_cc_refit_d" ? (K3, aNC3, aCN3, r.params[1], r.params[2]) :
                          r.name == "lr_cc_sharedK" ? (r.params[1], 1.0, 1.0, dN1b, dC1b) :
                          r.name == "lr_cc_K_only" ? (r.params[1], aNC3, aCN3, dN1b, dC1b) : (K3, aNC3, aCN3, dN1b, dC1b)
    ["\\frac{dN}{dt}=$(tn(rN))\\,N\\left(1-\\frac{N+$(tn(aNC))\\,C}{$(tn(K))}\\right)-$(tn(dN))\\cdot $(tn(U_R))\\,N",
     "\\frac{dC}{dt}=$(tn(rC))\\,C\\left(1-\\frac{C+$(tn(aCN))\\,N}{$(tn(K))}\\right)-$(tn(dC))\\cdot $(tn(U_R))\\,C"]
end

# Stage 1b panels
s1b_html = String[]
for (line, _, _) in stage1_lines
    s = stage1b[line]; s1 = stage1[line]; X = s.sym
    rows1b = all_rows(s.rows, hd1b[line].rows)
    obs = [(name="$line LR observed", x=s.x[.!s.bad], y=s.y[.!s.bad], color=GREY, filled=true)]
    sols = [dense(r.model, ru0(r, [s.y[1]]), r.params, XMAX) for r in rows1b]
    extra = [("r_$X", s1.r, "inherited from Stage 1 (logistic)"), ("K_$X", s1.K, "inherited from Stage 1 (logistic)"), ("u_R", U_R, "fixed: (10 − 0.5) / (10 − 0.05) for 0.5% FBS")]
    views = [vgrid(svg_panel("$line low resource — $(r.label)", obs, [(t=ts, y=[u[oix(r, 1)] for u in sol.u], color=BLUE, dash="")];
                             marks=lr_marks(s.x, s.y, s.bad, sol, oix(r, 1), BLUE)),
                   fit_box(r, rows1b, [tex_stage1b(r.name, X, r.params, s1.r, s1.K)]; extra=extra,
                           note = r.name == "lr_subtractive_loss" ? "<p class='muted'>This d^LR (★ on the Candidate models page) is carried into Stage 3b.</p>" : ""))
             for (r, (ts, sol)) in zip(rows1b, sols)]
    overlay = svg_panel("$line low resource — all candidates", obs,
        [(t=ts, y=[u[oix(rows1b[i], 1)] for u in sol.u], color=CAND_COLORS[i], dash=CAND_DASHES[i]) for (i, (ts, sol)) in enumerate(sols)];
        lines=[(label=r.label, color=CAND_COLORS[i], dash=CAND_DASHES[i]) for (i, r) in enumerate(rows1b)],
        marks=[(x=s.x[i], y=s.y[i], kind=:excluded, color=GREY) for i in findall(s.bad)])
    excl = isempty(findall(s.bad)) ? "none" : join(["day $(Int(s.x[i])) ($(fmt(s.y[i])))" for i in findall(s.bad)], ", ")
    push!(s1b_html, "<div class='panel'><h3>$line low resource <span class='muted'>($(esc(LR_FILE)) / $(esc(s.sheet)))</span></h3>" *
        "<p class='muted'>Excluded from the fit: $excl.$(line in SHIFT_FROM_DAY1 ? " <b>Time shifted:</b> the day-0 seeding value was dropped and day 1 is treated as day 0 (every later day moves back by one)." : "")</p>" *
        bic_svg(rows1b) * comparison_table(rows1b) * compare_block("s1b-$(lowercase(line))", rows1b, views; overlay=overlay) * "</div>")
end

# Stage 3b: one comparison block; each candidate view shows the 6 cocultures
s3b_extra_base = [("r_N", rN, "inherited from Stage 1 (A2780 logistic)"), ("r_C", rC, "inherited from Stage 1 (A2780cis logistic)"), ("u_R", U_R, "fixed: (10 − 0.5) / (10 − 0.05) for 0.5% FBS")]
s3b_inh = Dict(
    "lr_cc_inherited" => [("K", K3, "inherited from Stage 3"), ("α_NC", aNC3, "inherited from Stage 3"), ("α_CN", aCN3, "inherited from Stage 3"),
                          ("d_N^LR", dN1b, "inherited from Stage 1b ★"), ("d_C^LR", dC1b, "inherited from Stage 1b ★")],
    "lr_cc_refit_d"   => [("K", K3, "inherited from Stage 3"), ("α_NC", aNC3, "inherited from Stage 3"), ("α_CN", aCN3, "inherited from Stage 3")],
    "lr_cc_refit_all" => Tuple{String,Float64,String}[],
    "lr_cc_sharedK"   => [("α_NC = α_CN", 1.0, "fixed (one-for-one crowding)"),
                          ("d_N^LR", dN1b, "inherited from Stage 1b ★"), ("d_C^LR", dC1b, "inherited from Stage 1b ★")],
    "lr_cc_K_only"    => [("α_NC", aNC3, "inherited from Stage 3"), ("α_CN", aCN3, "inherited from Stage 3"),
                          ("d_N^LR", dN1b, "inherited from Stage 1b ★"), ("d_C^LR", dC1b, "inherited from Stage 1b ★")])
s3b_all = all_rows(s3b_rows, hd3b_rows)
s3b_views = String[]
for r in s3b_all
    ts, sol = dense(r.model, ru0(r, u0_3b), r.params, XMAX)
    plots = String[]
    for (i, c) in enumerate(lr_cc)
        obs = [(name="sensitive (A2780)", x=c.xN[.!c.bN], y=c.yN[.!c.bN], color=BLUE, filled=true),
               (name="resistant (A2780cis)", x=c.xC[.!c.bC], y=c.yC[.!c.bC], color=ORANGE, filled=true)]
        fits = [(t=ts, y=[u[oix(r, 2i-1)] for u in sol.u], color=BLUE, dash=""), (t=ts, y=[u[oix(r, 2i)] for u in sol.u], color=ORANGE, dash="")]
        marks = vcat(lr_marks(c.xN, c.yN, c.bN, sol, oix(r, 2i-1), BLUE), lr_marks(c.xC, c.yC, c.bC, sol, oix(r, 2i), ORANGE))
        push!(plots, svg_panel("Run $(c.run), ratio $(c.ratio) S:R ($(c.label))", obs, fits; marks=marks))
    end
    push!(s3b_views, vgrid(join(plots), fit_box(r, s3b_all, tex_stage3b(r);
        note="<p class='muted'>Low resource, untreated. No ordinary coculture death term is carried from Stage 3 (as on the Candidate models page).</p>",
        extra=vcat(s3b_extra_base, get(s3b_inh, r.name, Tuple{String,Float64,String}[])))))
end
s3b_html = "<div class='panel'><h3>3 ratios × runs 1 and 2 (12 trajectories, one joint fit)</h3>" *
    bic_svg(s3b_all) * comparison_table(s3b_all) * compare_block("s3b", s3b_all, s3b_views) * "</div>"

outlier_table = "<div class='tscroll'><table><tr><th>Stage</th><th>Trajectory</th><th>Day</th><th>Reason</th><th>Observed (excluded)</th><th>Replacement on BIC-winning fit</th></tr>" *
    join(["<tr><td>$(esc(o.stage))</td><td>$(esc(o.trajectory))</td><td class='num'>$(Int(o.day))</td><td>$(esc(o.reason))</td><td class='num'>$(fmt(o.observed))</td><td class='num'>$(fmt(o.replacement))</td></tr>" for o in outlier_rows]) *
    "</table></div>"

eq(s) = "<div class='eq'>$s</div>"
nav_link(label, href, active) = active ?
    "<a href='$href' aria-current='page' class='nav nav-active'>$label</a>" :
    "<a href='$href' class='nav'>$label</a>"
home_svg = "<svg width='22' height='22' viewBox='0 0 24 24' aria-hidden='true'><path fill='currentColor' d='M12 3.1 2.8 10.6v10.3h6.1v-6.2h6.2v6.2h6.1V10.6L12 3.1Zm7.1 15.8h-2v-6.2H6.9v6.2h-2V11.6L12 5.8l7.1 5.8v7.3Z'/></svg>"
nav = "<nav aria-label='New data report navigation' class='topnav'><a href='../../../../index.html' class='nav nav-home' aria-label='Reports home' title='Reports home'>$home_svg</a>" *
    nav_link("Summary", "../../../../outputs/results/summary.html", false) *
    nav_link("New data", "../../../../outputs/new_data_report/report.html", false) *
    nav_link("Candidate models", "../reports/a2780_low_resource_modeling_plan.html", false) *
    nav_link("New data fitting", "report.html", true) *
    nav_link("Bootstrapping", "../bootstrap_fits/report.html", false) *
    nav_link("LR mechanisms", "../lr_mechanisms/report.html", false) *
    nav_link("Low-resource fitting", "../low_resource_staged_fits/report.html", false) * "</nav>"

# Optimal model at each stage = lowest BIC over both model families, with the best model of
# the other family for comparison.
opt_rows = String[]
pvals(r) = isempty(r.params) ? "no free parameters (all inherited)" :
    esc(join(["$(n) = $(shown(v))" * (isempty(r.se) || isnan(r.se[i]) ? "" : " ± $(fmt(1.96r.se[i]))") for (i, (n, v)) in enumerate(zip(r.pnames, r.params))], ", "))
function opt_row!(stage, cond, rows)
    o = rows[1]
    j = findfirst(r -> family(r) != family(o), rows)
    other = j === nothing ? "—" : @sprintf("%s: %s (ΔBIC %.2f)", family(rows[j]), esc(rows[j].label), rows[j].bic - o.bic)
    push!(opt_rows, @sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td class='num'>%d</td><td class='num'>%.2f</td><td>%s</td><td>%s</td></tr>",
        stage, cond, esc(o.label), esc(family(o)), length(o.params), o.bic, pvals(o), other))
end
for (line, _, _) in stage1_lines; opt_row!("1", line, stage1[line].rows); end
for (line, _, _) in stage1_lines; opt_row!("1b", line, all_rows(stage1b[line].rows, hd1b[line].rows)); end
for (line, _, _) in stage1_lines; opt_row!("2", line, all_rows(stage2[line].rows, hd2[line].rows)); end
opt_row!("3", "Ce0 + Ce1 coculture", s3_all)
opt_row!("3b", "LR coculture", s3b_all)
opt_table = "<div class='tscroll'><table><tr><th>Stage</th><th>Condition</th><th>Optimal model (lowest BIC)</th><th>Family</th><th>Free params</th><th>BIC</th><th>Fitted values (± 95%)</th><th>Best model of the other family</th></tr>$(join(opt_rows, "\n"))</table></div>"

html = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>New data staged fitting</title>
<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}};</script><script defer src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>
<style>body{font:14px/1.5 system-ui,sans-serif;margin:0;background:#f4f7f8;color:#17212b}main{max-width:1400px;margin:auto;padding:24px}h1{color:#124b6e;margin:6px 0}h2{color:#124b6e;margin-top:34px;border-top:2px solid #cfdde3;padding-top:14px}
.topnav{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:0 0 16px}.nav{display:inline-flex;align-items:center;min-height:34px;padding:0 11px;border-radius:8px;color:#176b87;background:#eef5f6;text-decoration:none;font:650 14px system-ui,sans-serif;white-space:nowrap}.nav-home{padding:6px 9px}.nav-active{color:#fff;background:#176b87}
.intro,.panel{background:#fff;border:1px solid #d9e2e6;border-radius:12px;padding:16px 18px;margin:14px 0}.intro{border-left:6px solid #176b87}.note{background:#fff6d9;border-radius:8px;padding:10px 14px}
.vgrid>*{min-width:0}code{overflow-wrap:anywhere}.intro,.panel{min-width:0;overflow-wrap:anywhere}mjx-container:not([display]){display:inline-block;max-width:100%;overflow-x:auto;vertical-align:middle}.tscroll{overflow-x:auto}.eq{overflow-x:auto;overflow-y:hidden;background:#fbfdfd;border:1px solid #e4ecef;border-radius:8px;padding:4px 12px;margin:8px 0}svg{display:block;width:100%;height:auto;margin:10px 0}svg.bic{max-width:900px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}td,th{border:1px solid #d3dde2;padding:6px 8px;text-align:left;vertical-align:top}th{background:#e9f3f8}td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}tr.win td{background:#eef8f1;font-weight:600}table.params tr.fitted td{background:#f3faf5}.muted{color:#5b6b73;font-weight:400}h3{margin:0 0 6px}h4{margin:0 0 4px;font-size:16px}h5{margin:14px 0 4px;font-size:13px;text-transform:uppercase;letter-spacing:.04em;color:#52636b}
.cmp{margin-top:16px;border-top:1px dashed #cfdde3;padding-top:12px}.switch-label{margin:0 0 6px;font-weight:700;color:#124b6e}.cbtns{display:flex;gap:8px;flex-wrap:wrap;margin-bottom:10px}.cbtn{display:flex;flex-direction:column;align-items:flex-start;border:1px solid #176b87;background:#f7fbfd;color:#124b6e;border-radius:8px;padding:7px 12px;font:650 13px system-ui,sans-serif;cursor:pointer;text-align:left}.cbtn .d{font-weight:400;font-size:11.5px;opacity:.85}.cbtn.active,.cbtn:hover{background:#176b87;color:#fff}
.cview{display:none}.cview.active{display:block}.vgrid{display:grid;grid-template-columns:minmax(0,1.25fr) minmax(0,1fr);gap:20px;align-items:start}.fitbox{background:#fbfdfd;border:1px solid #e0e9ed;border-radius:10px;padding:12px 14px}.badge{display:inline-block;margin-left:6px;padding:1px 8px;border-radius:999px;background:#176b87;color:#fff;font-size:11px;vertical-align:middle}
@media(max-width:1050px){.vgrid{grid-template-columns:1fr}}@media(max-width:640px){main{padding:14px}}</style></head><body><main>
$nav
<h1>New data — staged parameter estimation</h1>
<div class='intro'><p>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")) with <code>GrowthParameterEstimation.run_joint_fit</code> on the Mean Cells trajectories in <code>New Datasets</code>. Two model families are compared: one population per line (candidates from the <a href='../reports/a2780_low_resource_modeling_plan.html'>Candidate models</a> page), and <b>healthy / damaged</b> two-population models (see below). At each stage every candidate is fitted, all are ranked together by BIC, and the <b>optimal</b> (lowest-BIC) model is reported.</p>
<ul><li><b>Stage 1 · Untreated monocultures</b> (A2780 = sensitive N, A2780cis = resistant C) → growth rate r and carrying capacity K per line.</li>
<li><b>Stage 1b · Low-resource monocultures</b> (30k sheets): inherit r and K; fit the low-resource terms.</li>
<li><b>Stage 2 · Treated monocultures</b>: inherit r and K; the three doses of each line are fitted <i>jointly</i> with shared drug parameters; <b>EC50 is fixed at 1 µM</b>.</li>
<li><b>Stage 3 · Untreated (Ce0) and treated (Ce1) cocultures together</b>: all three seeding ratios and both conditions (12 trajectories) are fitted at once with <b>one shared carrying capacity K</b>; Ce1 is 1 µM cisplatin.</li>
<li><b>Stage 3b · Low-resource cocultures</b> (30k, runs 1 and 2): inherit K and α from Stage 3 and d^LR from Stage 1b.</li></ul>
<p><b>How to compare models:</b> each fit shows a BIC chart and a comparison table for all candidates. Below them, the <i>Show fit for</i> buttons switch the plot, the fitted equation (with the estimated values substituted) and the parameter table to any candidate; <i>All candidates overlaid</i> draws every candidate on one plot. ΔBIC evidence bands: &lt; 2 weak, 2–6 positive, 6–10 strong, &gt; 10 very strong.</p>
<p class='note'>Low-resource data: only A2780 / A2780cis 30k sheets from <code>$(LR_FILE)</code> are used (Tyknu is not). FBS is 10% in normal medium (\\(u_R = 0\\), so Stages 1–3 have no low-resource term) and 0.5% in low resource, so \\(u_R = (10-0.5)/(10-0.05) = $(round(U_R; digits=4))\\) in Stages 1b and 3b. Sudden transient drops are excluded from the fit (a point below half of the previous kept value, by more than 50 cells, after which the trajectory recovers to at least 80% of that value; steady declines and spikes are kept), as are values that exactly repeat an earlier day's value (copy/paste errors); each excluded point is drawn as a grey × and replaced by a hollow diamond placed on the fitted line afterwards. All excluded points are listed at the end. Doses are read from the file names: IC25 file = 1.47 µM, IC50 file = 1 µM, IC75 file = 0.62 µM. Points are observations; lines are model fits. All plots share the same axes (0–$(round(Int, XMAX)) days, 0–$(round(Int, YMAX)) cells).</p></div>

<h2>Optimal model at each stage</h2>
<div class='intro'><p><b>Optimal</b> = lowest BIC among every candidate of both model families. Asymmetric models (different parameters for the sensitive and resistant line, \\(\\alpha_{NC} \\ne \\alpha_{CN}\\)) are allowed throughout. ± values are 95% intervals (1.96 × asymptotic SE; residuals of mean trajectories are correlated, so treat them as lower bounds). The last column gives the best model of the other family and how much worse its BIC is.</p>$opt_table</div>

<h2>Healthy / damaged two-population models</h2>
<div class='intro'>
<p>A second model family, fitted alongside the one-population-per-line models at Stages 1b, 2, 3 and 3b. Each line X (N = A2780 naive/sensitive, C = A2780cis resistant) is split into <b>healthy</b> cells \\(X_H\\), which divide, and <b>damaged</b> cells \\(X_D\\), which no longer divide but still take up space and are still counted. A stress (cisplatin at Stages 2 and 3, low resource at Stages 1b and 3b) damages healthy cells at rate \\(k_X\\). Damaged cells either die at rate \\(\\delta_X\\) or are repaired back to the healthy pool at rate \\(\\rho_X\\):</p>
$(eq(raw"\[\frac{dX_H}{dt}=r_X X_H\left(1-\frac{X+\alpha_{XY}Y}{K}\right)-k_X X_H+\rho_X X_D,\qquad \frac{dX_D}{dt}=k_X X_H-(\rho_X+\delta_X)X_D,\qquad X=X_H+X_D\]"))
$(eq(raw"\[k_X=\kappa_X H_X(D)\ \text{(cisplatin; } H_X=\tfrac{D^{h_X}}{1+D^{h_X}}\text{ or } D\text{)},\qquad k_X=\sigma_X u_R\ \text{(low resource)}\]"))
<p>Variants: <b>irreversible</b> (\\(\\rho = 0\\)); <b>repair</b> (\\(\\rho > 0\\)); <b>two-step transit</b> (damaged cells pass \\(X_{D1}\\to X_{D2}\\) at rate \\(\\delta\\) before dying, a longer and more peaked delay); and, for low resource, <b>arrest without death</b> (\\(\\delta = 0\\), damaged = quiescent). The damaged compartment replaces the \\(A(t)=1-e^{-\\lambda t}\\) delay of the one-population models: the delay comes from the time cells spend damaged before they die. Every rate is separate for N and C (asymmetric).</p>
<p><b>Carried forward within this family:</b> the Stage 2 winner's damage structure and \\(\\kappa, h, \\delta, \\rho\\) go to Stage 3; the Stage 1b winner's \\(\\sigma, \\delta, \\rho\\) and the Stage 3 winner's \\(K, \\alpha\\) go to Stage 3b. The one-population family keeps its own chain as before.</p>
<p class='muted'><b>Literature basis.</b> Proliferating and damaged (non-dividing, dying) compartments with transit through damage stages: Simeoni et al. 2004 (<i>Cancer Res</i> 64:1094) and the in-vitro transit-compartment model of Lobo &amp; Balthasar 2002 (<i>AAPS PharmSci</i> 4:E42). Repair of cisplatin–DNA lesions, faster in cisplatin-resistant A2780 derivatives, which motivates a separate \\(\\rho_C\\): Masuda et al. 1988 (<i>Cancer Res</i> 48:5713) and Parker et al. 1991 (<i>J Clin Invest</i> 87:772), which also reports lower platinum uptake (a lower \\(\\kappa_C\\)). Reversible switching between proliferating and quiescent cells under nutrient limitation: Gyllenberg &amp; Webb 1989 (<i>Growth Dev Aging</i> 53:25).</p></div>

<h2>Stage 1 — Untreated monoculture</h2>
<div class='intro'>
$(eq(raw"\[\text{Logistic ★:}\quad \frac{dX}{dt}=r_X X\left(1-\frac{X}{K_X}\right)\]"))
$(eq(raw"\[\text{Theta logistic:}\quad \frac{dX}{dt}=r_X X\left[1-\left(\frac{X}{K_X}\right)^{\theta_X}\right]\]"))
$(eq(raw"\[\text{Gompertz:}\quad \frac{dX}{dt}=r_X X\log\left(\frac{K_X}{X}\right)\]"))
<p>Each cell line is fitted separately. Carried forward: the logistic \\(r_X, K_X\\).</p></div>
$(join(s1_html, "\n"))

<h2>Stage 1b — Low-resource monoculture</h2>
<div class='intro'>
$(eq(raw"\[\text{Lineage-specific subtractive loss ★:}\quad \frac{dX}{dt}=r_X X\left(1-\frac{X}{K_X}\right)-d_X^{LR}\,u_R(t)\,X,\qquad u_R=\frac{10-F}{10-0.05}\]"))
$(eq(raw"\[\text{Growth-rate and carrying-capacity modifiers:}\quad \frac{dX}{dt}=\rho_r r_X X\left(1-\frac{X}{\rho_K K_X}\right)\]"))
$(eq(raw"\[\text{Separately fitted growth rate:}\quad \frac{dX}{dt}=r_X^{LR} X\left(1-\frac{X}{K_X}\right)\]"))
$(eq(raw"\[\text{Accumulated damage, 4 steps (delayed loss):}\quad \frac{dX}{dt}=r_X X\left(1-\frac{X}{K_X}\right)-d_X^{LR}\,u_R\,A_4(t)\,X,\qquad A_n(t)=1-e^{-\lambda_X^{LR} t}\sum_{k=0}^{n-1}\frac{(\lambda_X^{LR} t)^k}{k!}\]"))
$(eq(raw"\[\text{Immediate loss + accumulated damage, } n = 4 \text{ or } 8:\quad \frac{dX}{dt}=r_X X\left(1-\frac{X}{K_X}\right)-\left(d_X^{LR}+d_X^{acc}A_n(t)\right)u_R\,X\]"))
<p class='muted'>Accumulated damage: cells pass through \\(n\\) damage steps at rate \\(\\lambda\\) before the extra loss acts (mean delay \\(n/\\lambda\\)), the transit-compartment idea used for delayed drug effects (Lobo &amp; Balthasar 2002; Simeoni et al. 2004).</p>
<p>\\(r_X, K_X\\) are inherited from Stage 1. Sheets: <code>1 - mono A2780 30k</code> and <code>1 - mono A2780cis 30k</code>. Carried forward: \\(d_X^{LR}\\) from the ★ model.</p></div>
$(join(s1b_html, "\n"))

<h2>Stage 2 — Treated monoculture</h2>
<div class='intro'>
$(eq(raw"\[A_X(t)=1-e^{-\lambda_X t},\qquad H_X(D)=\frac{D^{h_X}}{EC_{50}^{h_X}+D^{h_X}},\qquad EC_{50}=1\ \mu\text{M (fixed)}\]"))
$(eq(raw"\[\text{Time-delayed Hill ★:}\quad \frac{dX}{dt}=r_X X\left(1-\frac{X}{K_X}\right)-\beta_X A_X(t) H_X(D) X\]"))
$(eq(raw"\[\text{No Hill saturation:}\quad \frac{dX}{dt}=r_X X\left(1-\frac{X}{K_X}\right)-\beta_X A_X(t)\, D\, X\]"))
<p>\\(r_X, K_X\\) are fixed from Stage 1. Free: \\(\\beta_X, \\lambda_X\\) (and \\(h_X\\) for the Hill model), shared across 0.62, 1.0 and 1.47 µM. The "distributed" form on the Candidate models page, \\(r_X X - \\tfrac{r_X}{K_X}X^2 - \\beta_X A_X H_X X\\), is the same equation expanded, so it is not fitted separately. Carried forward: \\(\\lambda_X\\) and \\(h_X\\).</p></div>
$(join(s2_html, "\n"))

<h2>Stage 3 — Untreated + treated coculture, shared carrying capacity</h2>
<div class='intro'>
$(eq(raw"\[\frac{dN}{dt}=r_N N\left(1-\frac{N+\alpha_{NC}C}{K}\right)-d_N N-\mathbb{1}_{\text{Ce1}}\,\beta_N A_N(t) H_N(D) N\]"))
$(eq(raw"\[\frac{dC}{dt}=r_C C\left(1-\frac{C+\alpha_{CN}N}{K}\right)-d_C C-\mathbb{1}_{\text{Ce1}}\,\beta_C A_C(t) H_C(D) C\]"))
<p>N = A2780 (sensitive), C = A2780cis (resistant). <b>K is one parameter shared by both lineages, all three ratios and both conditions.</b> Candidates switch terms on and off: one-for-one (\\(\\alpha_{NC}=\\alpha_{CN}=1\\)) or asymmetric crowding; treatment \\(\\beta, \\lambda, h\\) inherited from Stage 2 or \\(\\beta\\) refitted; death terms on neither, the sensitive line only, or both lines; a kill term on the resistant line or not. The kill term acts only on Ce1 (D = 1 µM). Parameter tables give asymptotic standard errors and 95% intervals; because consecutive points of a mean trajectory are correlated, treat them as lower bounds on the true uncertainty.</p>
<p class='muted'>Inherited: $(esc(s3_inherited)). With D = EC50 = 1 µM the Hill factor is exactly 0.5 whatever h is, so \\(h\\) cannot be re-estimated from the coculture data.</p></div>
$s3_html

<h2>Stage 3b — Low-resource coculture</h2>
<div class='intro'>
$(eq(raw"\[\frac{dN}{dt}=r_N N\left(1-\frac{N+\alpha_{NC}C}{K}\right)-d_N^{LR}\,u_R\,N,\qquad \frac{dC}{dt}=r_C C\left(1-\frac{C+\alpha_{CN}N}{K}\right)-d_C^{LR}\,u_R\,C\]"))
<p>The ★ low-resource coculture model from the Candidate models page, fitted jointly to runs 1 and 2 at 30k (3 ratios each; sheet ratios are sensitive:resistant). Five versions are compared: every parameter inherited (\\(K, \\alpha\\) from Stage 3, \\(d^{LR}\\) from Stage 1b; nothing refitted); only the shared \\(K\\) refitted, with one-for-one crowding (\\(\\alpha = 1\\)) or with \\(\\alpha\\) from Stage 3; \\(d^{LR}\\) refitted; and \\(K, \\alpha, d^{LR}\\) all refitted. \\(r_N, r_C\\) always come from Stage 1.</p>
<p class='muted'>Inherited: $(esc(s3b_inherited)).</p></div>
$s3b_html

<h2>Excluded low-resource points</h2>
<div class='intro'><p>$(length(outlier_rows)) point(s) were excluded from the low-resource fits. The replacement value is read off the BIC-winning fit for that stage after fitting, so it does not influence the fit.</p>$outlier_table</div>

<p class='muted'>Outputs: <code>bic_model_ranking.csv</code>, <code>carried_parameters.csv</code> in <code>Modeling_Approaches/02_mechanical_automatic_package/outputs/new_data_fits/</code>. Generated by <code>run_new_data_staged_fits.jl</code>.</p>
</main>
<script>
document.querySelectorAll('.cmp').forEach(function (c) {
  c.querySelectorAll('.cbtn').forEach(function (b) {
    b.addEventListener('click', function () {
      c.querySelectorAll('.cbtn').forEach(function (x) { x.classList.toggle('active', x === b); });
      c.querySelectorAll(':scope > .cview').forEach(function (v) { v.classList.toggle('active', v.id === b.dataset.view); });
      var view = document.getElementById(b.dataset.view);
      if (window.MathJax && MathJax.typesetPromise) MathJax.typesetPromise([view]);
    });
  });
});
</script></body></html>"""

write(joinpath(OUT_DIR, "report.html"), html)
mkpath(DOCS_DIR)
cp(joinpath(OUT_DIR, "report.html"), joinpath(DOCS_DIR, "report.html"); force=true)
println("\nWrote $(joinpath(OUT_DIR, "report.html"))\n  and $(joinpath(DOCS_DIR, "report.html"))")
