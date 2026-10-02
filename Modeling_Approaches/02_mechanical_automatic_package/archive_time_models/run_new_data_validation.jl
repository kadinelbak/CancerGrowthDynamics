# Validation: predict the held-out low-resource + cisplatin 1:1 coculture (30k) using ONLY
# parameters fitted in run_new_data_staged_fits.jl. Nothing is fitted here. Every combination of
# capacity/crowding, low-resource loss and kill term is simulated from the first measured day
# (day 1) and scored against the data. Writes outputs/new_data_validation/ and the docs copy.

import Pkg
Pkg.activate(@__DIR__)

using CSV, DataFrames, Dates, Printf
using OrdinaryDiffEq: ODEProblem, solve, Tsit5

const ROOT     = normpath(joinpath(@__DIR__, "..", ".."))
const VAL_DIR  = joinpath(ROOT, "New Datasets", "Validation")
const FITS     = joinpath(@__DIR__, "outputs", "new_data_fits", "all_fitted_parameters.csv")
const OUT_DIR  = joinpath(@__DIR__, "outputs", "new_data_validation")
const DOCS_DIR = joinpath(ROOT, "docs", "Modeling_Approaches", "02_mechanical_automatic_package", "outputs", "new_data_validation")
const FILE_N   = "ALL A2780 LR + cis(A2780 30k 1-1).csv"
const FILE_C   = "ALL A2780 LR + cis(A2780cis 30k 1-1).csv"
const DOSE     = 1.0                       # cisplatin, µM (IC50, as in Ce1)
const EC50     = 1.0
const U_R      = (10 - 0.5) / (10 - 0.05)  # 0.5% FBS
mkpath(OUT_DIR); mkpath(DOCS_DIR)

# ── data ───────────────────────────────────────────────────────────────────────
function load_mean(fname)
    df = CSV.read(joinpath(VAL_DIR, fname), DataFrame; silencewarnings=true)
    rename!(df, Symbol.(strip.(replace.(string.(names(df)), "﻿" => ""))))
    df = filter(r -> !ismissing(r.Day) && !ismissing(r[Symbol("Mean Cells")]), df)
    x = Float64.(df.Day); y = Float64.(df[!, Symbol("Mean Cells")])
    o = sortperm(x); x[o], y[o]
end

# same exclusion rule as the fitting report (exact duplicates, then transient drops)
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

xN, yN = load_mean(FILE_N); xC, yC = load_mean(FILE_C)
bN, wN = flag_drops(yN); bC, wC = flag_drops(yC)
# Day 0 is the nominal seeding count (SD 0 across wells), not a measurement, so the primary
# score starts each prediction from the first measured day and scores the days after it.
const T0 = 1.0
iT0 = findfirst(==(T0), xN)
u0 = [yN[iT0], yC[iT0]]
u0_seed = [yN[1], yC[1]]

# ── parameters from the staged fits ───────────────────────────────────────────
P = CSV.read(FITS, DataFrame)
function params_of(stage, line, model)
    sub = filter(r -> string(r.stage) == stage && r.line == line && r.model == model, P)
    isempty(sub) && error("no parameters for $stage / $line / $model")
    Dict(String(r.param) => Float64(r.value) for r in eachrow(sub))
end
# the combination grid uses the one-population family; healthy/damaged (hd_*) models are predicted separately below
best_model(stage, line) = (sub = filter(r -> string(r.stage) == stage && r.line == line && !startswith(r.model, "hd_"), P); String(sub[argmin(sub.bic), :model]))

s1N, s1C = params_of("1", "A2780", "logistic"), params_of("1", "A2780cis", "logistic")
rN, rC = s1N["r_N"], s1C["r_C"]
dNLR = params_of("1b", "A2780", "lr_subtractive_loss")["d_N^LR"]
dCLR = params_of("1b", "A2780cis", "lr_subtractive_loss")["d_C^LR"]
hill(D, h) = D^h / (EC50^h + D^h)
m2N, m2C = best_model("2", "A2780"), best_model("2", "A2780cis")
p2N, p2C = params_of("2", "A2780", m2N), params_of("2", "A2780cis", m2C)
kill2 = (βN=p2N["β_N"], λN=p2N["λ_N"], HN=m2N == "delayed_hill_kill" ? hill(DOSE, p2N["h_N"]) : DOSE,
         βC=p2C["β_C"], λC=p2C["λ_C"], HC=m2C == "delayed_hill_kill" ? hill(DOSE, p2C["h_C"]) : DOSE)
m3 = best_model("3", "coculture"); p3 = params_of("3", "coculture", m3)
p3s = params_of("3", "coculture", "shared_K_only")
m3b = best_model("3b", "coculture"); p3b = params_of("3b", "coculture", m3b)
p3bs = params_of("3b", "coculture", "lr_cc_sharedK")
kill3 = (βN=get(p3, "β_N", 0.0), λN=kill2.λN, HN=kill2.HN, βC=get(p3, "β_C", 0.0), λC=kill2.λC, HC=kill2.HC)

# ── the combination grid ──────────────────────────────────────────────────────
# Four ingredients, every one taken from an earlier fit. Only how they are combined changes.
p1bρN, p1bρC = params_of("1b", "A2780", "lr_rho_modifiers"), params_of("1b", "A2780cis", "lr_rho_modifiers")
rLRN = params_of("1b", "A2780", "lr_separate_r")["r_N^LR"]
rLRC = params_of("1b", "A2780cis", "lr_separate_r")["r_C^LR"]

capacities = [
    (key="S3opt+d", label="Stage 3 optimal, with sensitive death", K=p3["K"], aNC=get(p3, "α_NC", 1.0), aCN=get(p3, "α_CN", 1.0), dN=get(p3, "d_N", 0.0), dC=get(p3, "d_C", 0.0)),
    (key="S3opt",   label="Stage 3 optimal K and α (no death)",    K=p3["K"], aNC=get(p3, "α_NC", 1.0), aCN=get(p3, "α_CN", 1.0), dN=0.0, dC=0.0),
    (key="S3shK",   label="Stage 3 shared K only (α = 1)",           K=p3s["K"], aNC=1.0, aCN=1.0, dN=0.0, dC=0.0),
    (key="S3bopt",  label="Stage 3b optimal K and α",                K=p3b["K"], aNC=get(p3b, "α_NC", 1.0), aCN=get(p3b, "α_CN", 1.0), dN=0.0, dC=0.0),
    (key="S3bshK",  label="Stage 3b shared K only (α = 1)",          K=p3bs["K"], aNC=1.0, aCN=1.0, dN=0.0, dC=0.0),
]
# how low resource acts (all Stage 1b candidates, applied per line)
# immediate + accumulated (8-step) loss; A2780cis's Stage 1b fit counted time from day 1
p1bAN, p1bAC = params_of("1b", "A2780", "lr_loss_plus_accum8"), params_of("1b", "A2780cis", "lr_loss_plus_accum8")
noacc = (accN=0.0, laN=1.0, t0N=0.0, accC=0.0, laC=1.0, t0C=0.0)
lowres = [
    merge((key="none", label="no low-resource effect",                 short="none",      lossN=0.0, lossC=0.0, grN=1.0, grC=1.0, kN=1.0, kC=1.0), noacc),
    merge((key="loss", label="extra death (Stage 1b subtractive loss ★)", short="extra death", lossN=dNLR * U_R, lossC=dCLR * U_R, grN=1.0, grC=1.0, kN=1.0, kC=1.0), noacc),
    merge((key="rho",  label="lower growth rate and capacity (Stage 1b multipliers)", short="lower r and K", lossN=0.0, lossC=0.0,
                 grN=p1bρN["ρ_r"], grC=p1bρC["ρ_r"], kN=p1bρN["ρ_K"], kC=p1bρC["ρ_K"]), noacc),
    merge((key="rLR",  label="lower growth rate only (Stage 1b refitted r)", short="lower r", lossN=0.0, lossC=0.0, grN=rLRN / rN, grC=rLRC / rC, kN=1.0, kC=1.0), noacc),
    (key="accum", label="immediate + accumulated loss (Stage 1b, 8 steps)", short="delayed death", lossN=p1bAN["d_N^LR"] * U_R, lossC=p1bAC["d_C^LR"] * U_R,
     grN=1.0, grC=1.0, kN=1.0, kC=1.0, accN=p1bAN["d_N^acc"] * U_R, laN=p1bAN["λ_N^LR"], t0N=0.0, accC=p1bAC["d_C^acc"] * U_R, laC=p1bAC["λ_C^LR"], t0C=1.0),
]
erlang_onset(x, n) = x <= 0 ? 0.0 : 1 - exp(-x) * sum(x^k / factorial(k) for k in 0:n-1)
nokill = (βN=0.0, λN=kill2.λN, HN=kill2.HN, βC=0.0, λC=kill2.λC, HC=kill2.HC)
# how cisplatin acts
drugs = [
    (key="none",   label="no drug effect",                                  short="none",              k=nokill,  mode=:kill, delay=true),
    (key="S2kill", label="kills cells after a delay (Stage 2 fit ★)",      short="delayed kill (S2)", k=kill2,   mode=:kill, delay=true),
    (key="S2now",  label="kills cells from day 0 (Stage 2 β, no delay)",   short="immediate kill (S2)", k=kill2, mode=:kill, delay=false),
    (key="S2stat", label="stops growth instead of killing (Stage 2 β)",    short="growth arrest (S2)", k=kill2,  mode=:static, delay=true),
    (key="S3kill", label="kills cells after a delay (Stage 3 coculture fit)", short="delayed kill (S3)", k=kill3, mode=:kill, delay=true),
]
# how drug and growth conditions interact
interactions = [
    (key="indep",  label="independent: kill is the same whatever the growth conditions", short="independent"),
    (key="prolif", label="dividing cells only: kill scales with the room left to grow", short="dividing cells only"),
]

function make_model(cap, lr, drug, inter)
    k = drug.k
    rNe, rCe = rN * lr.grN, rC * lr.grC
    KN, KC = cap.K * lr.kN, cap.K * lr.kC
    function f!(du, u, p, t)
        N = max(u[1], 0.0); C = max(u[2], 0.0)
        gN = 1 - (N + cap.aNC*C)/KN
        gC = 1 - (C + cap.aCN*N)/KC
        AN = drug.delay ? 1 - exp(-k.λN*t) : 1.0
        AC = drug.delay ? 1 - exp(-k.λC*t) : 1.0
        effN, effC = k.βN*AN*k.HN, k.βC*AC*k.HC
        if drug.mode == :static
            # the drug removes up to its kill rate from the per-capita growth rate, never below zero
            growN = rNe*max(1 - effN/rNe, 0.0)*N*gN
            growC = rCe*max(1 - effC/rCe, 0.0)*C*gC
            killN = killC = 0.0
        else
            growN, growC = rNe*N*gN, rCe*C*gC
            sN = inter.key == "prolif" ? clamp(gN, 0.0, 1.0) : 1.0
            sC = inter.key == "prolif" ? clamp(gC, 0.0, 1.0) : 1.0
            killN, killC = effN*sN*N, effC*sC*C
        end
        accN = lr.accN > 0 ? lr.accN*erlang_onset(lr.laN*(t - lr.t0N), 8) : 0.0
        accC = lr.accC > 0 ? lr.accC*erlang_onset(lr.laC*(t - lr.t0C), 8) : 0.0
        du[1] = growN - (cap.dN + lr.lossN + accN)*N - killN
        du[2] = growC - (cap.dC + lr.lossC + accC)*C - killC
    end
end

simulate(model, t0, u) = solve(ODEProblem(model, u, (t0, 14.0)), Tsit5(); reltol=1e-9, abstol=1e-9)
sse(v) = sum(abs2, v)
function score(sol, keep_only, t0)
    mN = (keep_only ? .!bN : trues(length(xN))) .& (xN .> t0)
    mC = (keep_only ? .!bC : trues(length(xC))) .& (xC .> t0)
    rNv = [sol(x)[1] - y for (x, y) in zip(xN[mN], yN[mN])]
    rCv = [sol(x)[2] - y for (x, y) in zip(xC[mC], yC[mC])]
    n = length(rNv) + length(rCv); S = sse(rNv) + sse(rCv)
    r2(r, y) = 1 - sse(r) / sse(y .- sum(y)/length(y))
    (bic=n*log(S/n), sse=S, n=n, rmseN=sqrt(sse(rNv)/length(rNv)), rmseC=sqrt(sse(rCv)/length(rCv)),
     r2N=r2(rNv, yN[mN]), r2C=r2(rCv, yC[mC]))
end

combos = NamedTuple[]
for cap in capacities, lr in lowres, drug in drugs, inter in interactions
    # the interaction only changes kill terms; skip duplicates where there is no kill
    inter.key == "prolif" && (drug.key == "none" || drug.mode == :static) && continue
    m = make_model(cap, lr, drug, inter)
    sol, sol0 = simulate(m, T0, u0), simulate(m, 0.0, u0_seed)
    sc, sa, s0 = score(sol, true, T0), score(sol, false, T0), score(sol0, true, 0.0)
    star = cap.key == "S3opt" && lr.key == "loss" && drug.key == "S2kill" && inter.key == "indep"
    push!(combos, (cap=cap, lr=lr, drug=drug, inter=inter, sol=sol, sc=sc, sa=sa, s0=s0, star=star,
                   name="$(cap.label) · low resource: $(lr.short) · drug: $(drug.short) · $(inter.short)"))
end
sort!(combos; by=c -> c.sc.bic)
best = combos[1]
starc = combos[findfirst(c -> c.star, combos)]
star_rank = findfirst(c -> c.star, combos)

struct Flat end
(::Flat)(t) = u0
baseline = score(Flat(), true, T0)
n_beat = count(c -> c.sc.bic < baseline.bic, combos)

println("Validation: $(length(combos)) combinations, n = $(best.sc.n) points (excluded N days: $(Int.(xN[bN])), C days: $(Int.(xC[bC])))")
@printf("  baseline (no change after day 1): BIC=%.2f; %d combinations beat it\n", baseline.bic, n_beat)
for (i, c) in enumerate(combos[1:min(end, 12)])
    @printf("  %3d  BIC=%8.2f  RMSE N=%7.1f C=%7.1f  R² N=%6.2f C=%6.2f  %s%s\n", i, c.sc.bic, c.sc.rmseN, c.sc.rmseC, c.sc.r2N, c.sc.r2C, c.name, c.star ? "  ★" : "")
end
@printf("  ★ rank %d\n", star_rank)

CSV.write(joinpath(OUT_DIR, "validation_ranking.csv"), DataFrame([
    (rank=i, capacity=c.cap.label, low_resource=c.lr.label, drug=c.drug.label, interaction=c.inter.label, candidate_models_star=c.star,
     bic=c.sc.bic, delta_bic=c.sc.bic - best.sc.bic, beats_no_change_baseline=c.sc.bic < baseline.bic, sse=c.sc.sse, n=c.sc.n,
     rmse_sensitive=c.sc.rmseN, rmse_resistant=c.sc.rmseC, r2_sensitive=c.sc.r2N, r2_resistant=c.sc.r2C,
     bic_all_points=c.sa.bic, bic_from_day0_seeding=c.s0.bic) for (i, c) in enumerate(combos)]))

# ── healthy / damaged predictions ─────────────────────────────────────────────
# Both stresses act at once here, so each line gets a cisplatin damage channel and a
# low-resource damage channel, with parameters from the healthy/damaged stage winners.
include(joinpath(@__DIR__, "healthy_damaged.jl"))
const HD_FILE = joinpath(@__DIR__, "outputs", "new_data_fits", "hd_resolved_parameters.csv")
struct ObsSol{S}; sol::S; idx::Vector{Int}; end
(o::ObsSol)(t) = o.sol(t)[o.idx]

hd_combos = NamedTuple[]
if isfile(HD_FILE)
    HDP = CSV.read(HD_FILE, DataFrame)
    hdv(stage, line) = Dict(String(r.param) => Float64(r.value) for r in eachrow(filter(r -> string(r.stage) == stage && r.line == line, HDP)))
    hdm(stage, line) = String(first(filter(r -> string(r.stage) == stage && r.line == line, HDP)).model)
    hfac(d) = d["hill"] == 1.0 ? hill(DOSE, d["h"]) : DOSE
    lines_hd = ("A2780", "A2780cis")
    drug_src = [(key="S2", label="cisplatin damage from Stage 2 (monoculture)",
                 v=[(d = hdv("2", l); (n=Int(d["n"]), κ=d["κ"], ρ=d["ρ"], δ=d["δ"], Hf=hfac(d), λ=d["λ"], b=0.0)) for l in lines_hd], model=hdm("2", "A2780") * " / " * hdm("2", "A2780cis")),
                (key="S3", label="cisplatin damage from Stage 3 (coculture)",
                 v=[(d = hdv("3", l); (n=Int(d["n"]), κ=d["κ"], ρ=d["ρ"], δ=d["δ"], Hf=d["Hf"], λ=d["λ"], b=d["b"])) for l in lines_hd], model=hdm("3", "coculture"))]
    lr_src = [(key="S1b", label="starvation damage from Stage 1b (monoculture)",
               v=[(d = hdv("1b", l); (n=Int(d["n"]), σ=d["σ"], ρ=d["ρ"], δ=d["δ"], dens=d["dens"] == 1.0)) for l in lines_hd], model=hdm("1b", "A2780") * " / " * hdm("1b", "A2780cis")),
              (key="S3b", label="starvation damage from Stage 3b (coculture)",
               v=[(d = hdv("3b", l); (n=Int(d["n"]), σ=d["σ"], ρ=d["ρ"], δ=d["δ"], dens=d["dens"] == 1.0)) for l in lines_hd], model=hdm("3b", "coculture"))]
    cap_src = [(key="S3", label="K, α from Stage 3", v=hdv("3", "coculture")), (key="S3b", label="K, α from Stage 3b", v=hdv("3b", "coculture"))]
    pools = [(key="separate", label="separate damage pools (each stress keeps its own repair and death rates)"),
             (key="one", label="one damage pool (both stresses feed it; cisplatin repair and death rates)")]
    for cap in cap_src, dr in drug_src, lr in lr_src, pool in pools
        sep = pool.key == "separate"
        ns = sep ? [[dr.v[l].n, lr.v[l].n] for l in 1:2] : [[dr.v[l].n] for l in 1:2]
        ch = sep ? Tuple(((1.0, dr.v[l].ρ, dr.v[l].δ), (1.0, lr.v[l].ρ, lr.v[l].δ)) for l in 1:2) :
                   Tuple(((1.0, dr.v[l].ρ, dr.v[l].δ),) for l in 1:2)
        # cisplatin: coculture damage b + κ·H(1 µM)·A(t); low resource: σ·u_R (× crowding g if density-dependent)
        kdrug(l, t) = dr.v[l].b + dr.v[l].κ * dr.v[l].Hf * (dr.v[l].λ > 0 ? 1 - exp(-dr.v[l].λ * t) : 1.0)
        klr(l, g) = lr.v[l].σ * U_R * (lr.v[l].dens ? g : 1.0)
        act = sep ? ((e, l, c, t, g) -> c == 1 ? kdrug(l, t) : klr(l, g)) :
                    ((e, l, c, t, g) -> kdrug(l, t) + klr(l, g))
        q = (r=(rN, rC), K=(cap.v["K"], cap.v["K"]), α=((1.0, cap.v["aNC"]), (cap.v["aCN"], 1.0)), ch=ch, act=act)
        m = make_hd(1, ns, p -> q)
        idx = hd_obs(1, ns)
        sol = ObsSol(solve(ODEProblem(m, hd_u0(1, ns, [[u0[1], u0[2]]]), (T0, 14.0), Float64[]), Tsit5(); reltol=1e-9, abstol=1e-9), idx)
        push!(hd_combos, (cap=cap, dr=dr, lr=lr, pool=pool, sol=sol, sc=score(sol, true, T0), sa=score(sol, false, T0),
                          name="$(cap.label) · $(dr.label) · $(lr.label) · $(pool.key == "separate" ? "separate pools" : "one pool")"))
    end
    sort!(hd_combos; by=c -> c.sc.bic)
    CSV.write(joinpath(OUT_DIR, "hd_validation_ranking.csv"), DataFrame([
        (rank=i, capacity=c.cap.label, cisplatin=c.dr.label, low_resource=c.lr.label, pools=c.pool.label, bic=c.sc.bic, sse=c.sc.sse, n=c.sc.n,
         rmse_sensitive=c.sc.rmseN, rmse_resistant=c.sc.rmseC, r2_sensitive=c.sc.r2N, r2_resistant=c.sc.r2C, bic_all_points=c.sa.bic) for (i, c) in enumerate(hd_combos)]))
    println("Healthy/damaged predictions: $(length(hd_combos)); best BIC $(round(hd_combos[1].sc.bic; digits=2)): $(hd_combos[1].name)")
end

# ── summaries: each ingredient alone, and every pair together ─────────────────
opts_of(field) = field == :cap ? capacities : field == :lr ? lowres : field == :drug ? drugs : interactions
factor_title = Dict(:cap => "Capacity and crowding", :lr => "How low resource acts", :drug => "How cisplatin acts", :inter => "Drug × growth interaction")
function marginal(field)
    rows = NamedTuple[]
    for o in opts_of(field)
        cs = [c for c in combos if getfield(c, field).key == o.key]
        isempty(cs) && continue
        push!(rows, (label=o.label, best_d=cs[1].sc.bic - best.sc.bic, mean_d=sum(c.sc.bic - best.sc.bic for c in cs) / length(cs),
                     best_rank=findfirst(c -> getfield(c, field).key == o.key, combos), n=length(cs)))
    end
    sort(rows; by=r -> r.best_d)
end

esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
fmt(v) = abs(v) >= 100 ? @sprintf("%.1f", v) : abs(v) < 1e-8 ? "0" : @sprintf("%.4g", v)
cellcolor(d) = (x = clamp(d / 60, 0, 1); @sprintf("hsl(%d, 60%%, %d%%)", round(Int, 130 - 130x), round(Int, 88 - 6x)))

function pair_table(f1, f2)
    o1, o2 = opts_of(f1), opts_of(f2)
    head = "<tr><th>$(esc(factor_title[f1])) ↓ / $(esc(factor_title[f2])) →</th>" * join(["<th>$(esc(o.short))</th>" for o in o2]) * "</tr>"
    body = String[]
    for a in o1
        cells = String[]
        for b in o2
            cs = [c for c in combos if getfield(c, f1).key == a.key && getfield(c, f2).key == b.key]
            push!(cells, isempty(cs) ? "<td class='num'>—</td>" :
                  (d = cs[1].sc.bic - best.sc.bic; "<td class='num' style='background:$(cellcolor(d))'>$(@sprintf("%.1f", d))$(cs[1].sc.bic < baseline.bic ? " ✓" : "")</td>"))
        end
        push!(body, "<tr><th>$(esc(a.short))</th>" * join(cells) * "</tr>")
    end
    "<div class='tscroll'><table class='pairs'>$head$(join(body))</table></div>"
end
short_cap =Dict("S3opt+d" => "S3 optimal + death", "S3opt" => "S3 optimal K, α", "S3shK" => "S3 shared K", "S3bopt" => "S3b optimal K, α", "S3bshK" => "S3b shared K")
capacities = [merge(c, (short=short_cap[c.key],)) for c in capacities]
combos = [merge(c, (cap=merge(c.cap, (short=short_cap[c.cap.key],)),)) for c in combos]

# ── report ─────────────────────────────────────────────────────────────────────
const BLUE, ORANGE = "#0b6e99", "#d95f02"
function panel(title, sol; W=760, H=340)
    ts = T0:0.1:14
    pmax = maximum(maximum(sol(t)[i] for t in ts) for i in 1:2)
    ymax = ceil(max(maximum(yN), maximum(yC), min(pmax, 6000)) * 1.08; sigdigits=2)
    L, R, T, B = 64, 20, 30, 70
    sx(v) = L + v/14*(W-L-R); sy(v) = H - B - clamp(v, 0, ymax)/ymax*(H-T-B)
    p = ["<svg viewBox='0 0 $W $H' role='img' aria-label='$(esc(title))'><rect width='100%' height='100%' fill='white'/><text x='$L' y='18' font-size='13' font-weight='bold'>$(esc(title))</text>"]
    for v in range(0, ymax; length=5)
        push!(p, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#e5e9ec'/><text x='%d' y='%.1f' font-size='10' text-anchor='end'>%.0f</text>", L, W-R, sy(v), sy(v), L-6, sy(v)+3, v))
    end
    for d in 0:2:14
        push!(p, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%d</text>", sx(d), H-B+14, d))
    end
    push!(p, "<line x1='$L' y1='$(H-B)' x2='$(W-R)' y2='$(H-B)' stroke='#555'/><line x1='$L' y1='$T' x2='$L' y2='$(H-B)' stroke='#555'/><text x='$((L+W-R)÷2)' y='$(H-B+30)' font-size='11' text-anchor='middle'>Day</text>")
    for (i, c) in ((1, BLUE), (2, ORANGE))
        pts = join([@sprintf("%.1f,%.1f", sx(t), sy(sol(t)[i])) for t in ts], " ")
        push!(p, "<polyline points='$pts' fill='none' stroke='$c' stroke-width='2.2'/>")
    end
    pmax > ymax && push!(p, "<text x='$(W-R)' y='$(T+12)' font-size='11' text-anchor='end' fill='#b45309'>prediction exceeds the axis (peak $(round(Int, pmax)))</text>")
    for (x, y, b, c) in ((xN, yN, bN, BLUE), (xC, yC, bC, ORANGE))
        for j in eachindex(x)
            push!(p, b[j] ? @sprintf("<path d='M%.1f %.1fL%.1f %.1fM%.1f %.1fL%.1f %.1f' stroke='#9ca3af' stroke-width='1.8'/>", sx(x[j])-4, sy(y[j])-4, sx(x[j])+4, sy(y[j])+4, sx(x[j])-4, sy(y[j])+4, sx(x[j])+4, sy(y[j])-4) :
                            @sprintf("<circle cx='%.1f' cy='%.1f' r='3.4' fill='%s'/>", sx(x[j]), sy(y[j]), c))
        end
    end
    lx = L
    for (lab, c) in (("A2780 (sensitive) data + prediction", BLUE), ("A2780cis (resistant) data + prediction", ORANGE))
        push!(p, "<circle cx='$(lx+4)' cy='$(H-16)' r='4' fill='$c'/><text x='$(lx+12)' y='$(H-12)' font-size='11' fill='$c'>$lab</text>"); lx += 250
    end
    push!(p, "<path d='M$(lx) $(H-20)L$(lx+8) $(H-12)M$(lx) $(H-12)L$(lx+8) $(H-20)' stroke='#9ca3af' stroke-width='1.8'/><text x='$(lx+14)' y='$(H-12)' font-size='11'>excluded (transient drop)</text></svg>")
    join(p)
end

TOPN = 25
rank_rows = join([@sprintf("<tr%s><td class='num'>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s%s</td><td class='num'>%.2f</td><td class='num'>%.2f</td><td class='num'>%.1f</td><td class='num'>%.1f</td><td class='num'>%.2f</td><td class='num'>%.2f</td><td class='num'>%.2f</td><td class='num'>%.2f</td></tr>",
    i == 1 ? " class='win'" : c.star ? " class='star'" : "", i, esc(c.cap.short), esc(c.lr.short), esc(c.drug.short), esc(c.inter.short), c.star ? " ★" : "",
    c.sc.bic, c.sc.bic - best.sc.bic, c.sc.rmseN, c.sc.rmseC, c.sc.r2N, c.sc.r2C, c.sa.bic, c.s0.bic)
    for (i, c) in enumerate(combos) if i <= TOPN || c.star], "\n") *
    @sprintf("\n<tr class='ref'><td></td><td colspan='4'>Reference: no model — both populations stay at their day-1 counts</td><td class='num'>%.2f</td><td class='num'>%.2f</td><td class='num'>%.1f</td><td class='num'>%.1f</td><td class='num'>%.2f</td><td class='num'>%.2f</td><td></td><td></td></tr>",
        baseline.bic, baseline.bic - best.sc.bic, baseline.rmseN, baseline.rmseC, baseline.r2N, baseline.r2C)

marg_html = join(["<h3>$(esc(factor_title[f]))</h3><div class='tscroll'><table><tr><th>Option</th><th>Best ΔBIC reached</th><th>Best rank</th><th>Mean ΔBIC</th><th>Combinations</th></tr>" *
    join(["<tr><td>$(esc(m.label))</td><td class='num'>$(@sprintf("%.1f", m.best_d))</td><td class='num'>$(m.best_rank)</td><td class='num'>$(@sprintf("%.1f", m.mean_d))</td><td class='num'>$(m.n)</td></tr>" for m in marginal(f)]) *
    "</table></div>" for f in (:lr, :drug, :inter, :cap)])

src_rows = join(["<tr><td>$(esc(a))</td><td>$(esc(b))</td></tr>" for (a, b) in [
    ("Growth rates", @sprintf("r_N = %.4f, r_C = %.4f /day (Stage 1 logistic)", rN, rC)),
    ("Stage 3 optimal ($(m3))", join(["$(k) = $(fmt(v))" for (k, v) in sort(collect(p3))], ", ")),
    ("Stage 3 shared K only", @sprintf("K = %.1f, α_NC = α_CN = 1", p3s["K"])),
    ("Stage 3b optimal ($(m3b))", join(["$(k) = $(fmt(v))" for (k, v) in sort(collect(p3b))], ", ")),
    ("Stage 3b shared K only", @sprintf("K = %.1f, α_NC = α_CN = 1", p3bs["K"])),
    ("Low resource: extra death (Stage 1b ★)", @sprintf("d_N^LR = %.4f, d_C^LR = %.4f /day, × u_R = %.4f", dNLR, dCLR, U_R)),
    ("Low resource: lower r and K (Stage 1b multipliers)", @sprintf("A2780: ρ_r = %.3f, ρ_K = %.3f; A2780cis: ρ_r = %.3g (at its bound), ρ_K = %.4f", p1bρN["ρ_r"], p1bρN["ρ_K"], p1bρC["ρ_r"], p1bρC["ρ_K"])),
    ("Low resource: lower r only (Stage 1b)", @sprintf("r_N^LR = %.4f, r_C^LR = %.4f /day", rLRN, rLRC)),
    ("Low resource: delayed death (Stage 1b, 8 steps)", @sprintf("A2780: d = %.4f + %.4f after mean delay %.1f d; A2780cis: d = %.4f + %.4f after mean delay %.1f d from day 1 (all × u_R)",
        p1bAN["d_N^LR"], p1bAN["d_N^acc"], 8 / p1bAN["λ_N^LR"], p1bAC["d_C^LR"], p1bAC["d_C^acc"], 8 / p1bAC["λ_C^LR"])),
    ("Drug: Stage 2 monoculture fit", @sprintf("A2780 (%s): β_N = %.4f, λ_N = %.4f, H_N(1 µM) = %.3f; A2780cis (%s): β_C = %.4f, λ_C = %.4f, H_C(1 µM) = %.3f", m2N, kill2.βN, kill2.λN, kill2.HN, m2C, kill2.βC, kill2.λC, kill2.HC)),
    ("Drug: Stage 3 coculture fit", @sprintf("β_N = %.4f (λ_N, H_N from Stage 2), β_C = %.4f", kill3.βN, kill3.βC)),
]])
data_rows = join(["<tr><td class='num'>$(Int(xN[i]))</td><td class='num'>$(fmt(yN[i]))$(bN[i] ? " (excluded: $(wN[i]))" : "")</td><td class='num'>$(fmt(yC[i]))$(bC[i] ? " (excluded: $(wC[i]))" : "")</td></tr>" for i in eachindex(xN)])

home_svg = "<svg width='22' height='22' viewBox='0 0 24 24' aria-hidden='true'><path fill='currentColor' d='M12 3.1 2.8 10.6v10.3h6.1v-6.2h6.2v6.2h6.1V10.6L12 3.1Zm7.1 15.8h-2v-6.2H6.9v6.2h-2V11.6L12 5.8l7.1 5.8v7.3Z'/></svg>"
navl(l, h, a) = a ? "<a href='$h' aria-current='page' class='nav nav-active'>$l</a>" : "<a href='$h' class='nav'>$l</a>"
nav = "<nav aria-label='New data report navigation' class='topnav'><a href='../../../../index.html' class='nav nav-home' aria-label='Reports home' title='Reports home'>$home_svg</a>" *
      navl("Summary", "../../../../outputs/results/summary.html", false) *
      navl("New data", "../../../../outputs/new_data_report/report.html", false) * navl("Candidate models", "../reports/a2780_low_resource_modeling_plan.html", false) *
      navl("New data fitting", "../new_data_fits/report.html", false) * navl("Validation", "report.html", true) *
      navl("Bootstrapping", "../bootstrap_fits/report.html", false) * navl("LR mechanisms", "../lr_mechanisms/report.html", false) *
      navl("Low-resource fitting", "../low_resource_staged_fits/report.html", false) * "</nav>"

# ── global model search (run_new_data_model_search.jl), if it has been run ────
struct Curve; t::Vector{Float64}; N::Vector{Float64}; C::Vector{Float64}; end
function (c::Curve)(t)
    i = clamp(searchsortedlast(c.t, t), 1, length(c.t) - 1)
    w = (t - c.t[i]) / (c.t[i+1] - c.t[i])
    [c.N[i] + w*(c.N[i+1] - c.N[i]), c.C[i] + w*(c.C[i+1] - c.C[i])]
end

function search_section(dir, heading, intro, factors)
    isfile(joinpath(dir, "model_search_results.csv")) || return ""
    R = CSV.read(joinpath(dir, "model_search_results.csv"), DataFrame)
    VC = CSV.read(joinpath(dir, "model_search_validation_curves.csv"), DataFrame)
    TC = CSV.read(joinpath(dir, "model_search_training_curves.csv"), DataFrame)
    TD = CSV.read(joinpath(dir, "model_search_training_data.csv"), DataFrame)
    sort!(R, :train_bic)
    curve(id) = (s = sort(filter(r -> r.id == id, VC), :t); Curve(Float64.(s.t), Float64.(s.N), Float64.(s.C)))
    tw = R[1, :]; vw = R[argmin(R.val_bic), :]
    vrank(id) = findfirst(==(id), sort(R, :val_bic).id)
    nbeat = count(<(baseline.bic), R.val_bic)
    # Spearman rank correlation between training and validation ranking
    rk(v) = invperm(sortperm(v))
    a, b = rk(R.train_bic), rk(R.val_bic); n = length(a)
    ρs = 1 - 6*sum(abs2, a .- b) / (n*(n^2 - 1))

    # scatter: training ΔBIC vs validation BIC
    W, H, L, Rm, T, B = 760, 360, 70, 20, 20, 50
    dx = R.train_bic .- R.train_bic[1]; xmax = quantile_cap(dx)
    ylo, yhi = min(minimum(R.val_bic), baseline.bic) - 5, min(maximum(R.val_bic), baseline.bic + 120)
    sx(v) = L + clamp(v, 0, xmax)/xmax*(W-L-Rm); sy(v) = H - B - (clamp(v, ylo, yhi) - ylo)/(yhi - ylo)*(H-T-B)
    sc = ["<svg viewBox='0 0 $W $H' role='img' aria-label='training versus validation'><rect width='100%' height='100%' fill='white'/>"]
    for v in range(0, xmax; length=5)
        push!(sc, @sprintf("<line x1='%.1f' x2='%.1f' y1='%d' y2='%d' stroke='#eef1f3'/><text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%.0f</text>", sx(v), sx(v), T, H-B, sx(v), H-B+14, v))
    end
    for v in range(ylo, yhi; length=5)
        push!(sc, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#eef1f3'/><text x='%d' y='%.1f' font-size='10' text-anchor='end'>%.0f</text>", L, W-Rm, sy(v), sy(v), L-6, sy(v)+3, v))
    end
    push!(sc, @sprintf("<line x1='%d' x2='%d' y1='%.1f' y2='%.1f' stroke='#b45309' stroke-dasharray='6,4'/><text x='%d' y='%.1f' font-size='11' fill='#b45309' text-anchor='end'>no-change reference</text>", L, W-Rm, sy(baseline.bic), sy(baseline.bic), W-Rm, sy(baseline.bic)-4))
    for r in eachrow(R)
        c = r.id == tw.id ? "#16a34a" : r.id == vw.id ? "#7c3aed" : "#0b6e99"
        push!(sc, @sprintf("<circle cx='%.1f' cy='%.1f' r='%.1f' fill='%s' fill-opacity='%.2f'/>", sx(r.train_bic - R.train_bic[1]), sy(r.val_bic), r.id in (tw.id, vw.id) ? 6 : 3.5, c, r.id in (tw.id, vw.id) ? 1 : 0.45))
    end
    push!(sc, "<text x='$((L+W-Rm)÷2)' y='$(H-12)' font-size='11' text-anchor='middle'>Training ΔBIC (0 = best fit to all training data; capped at $(round(Int, xmax)))</text><text transform='translate(14 $((T+H-B)÷2)) rotate(-90)' font-size='11' text-anchor='middle'>Validation BIC (lower is better)</text></svg>")

    row(r, i) = @sprintf("<tr%s><td class='num'>%d</td><td>%s</td><td class='num'>%d</td><td class='num'>%.1f</td><td class='num'>%.2f</td><td class='num'>%.2f</td><td class='num'>%.2f</td><td class='num'>%d</td></tr>",
        r.val_bic < baseline.bic ? " class='win'" : "", i, esc(r.label), r.k, r.train_bic - R.train_bic[1], r.val_bic, r.val_r2_sensitive, r.val_r2_resistant, vrank(r.id))
    head = "<tr><th>Training rank</th><th>Mechanisms</th><th>Parameters</th><th>Training ΔBIC</th><th>Validation BIC</th><th>R² sensitive</th><th>R² resistant</th><th>Validation rank</th></tr>"
    top_train = join([row(r, i) for (i, r) in enumerate(eachrow(first(R, 15)))])
    VS = sort(R, :val_bic)
    top_val = join([row(r, findfirst(==(r.id), R.id)) for r in eachrow(first(VS, 15))])

    factor_rows = String[]
    for (col, title, names) in factors
        for (key, lab) in sort(collect(names); by=last)
            sub = filter(r -> string(r[col]) == key, R)
            isempty(sub) && continue
            push!(factor_rows, @sprintf("<tr><td>%s</td><td>%s</td><td class='num'>%.1f</td><td class='num'>%.2f</td><td class='num'>%d / %d</td></tr>",
                title, esc(lab), minimum(sub.train_bic) - R.train_bic[1], minimum(sub.val_bic), count(<(baseline.bic), sub.val_bic), nrow(sub)))
        end
    end

    function tpanel(id, env)
        d = filter(r -> r.env == env, TD); c = sort(filter(r -> r.id == id && r.env == env, TC), :t)
        ym = ceil(max(maximum(d.value), maximum(c.N), maximum(c.C)) * 1.1; sigdigits=2)
        w2, h2, l2, t2, b2 = 380, 220, 44, 22, 26
        px(v) = l2 + v/14*(w2-l2-8); py(v) = h2 - b2 - clamp(v, 0, ym)/ym*(h2-t2-b2)
        q = ["<svg viewBox='0 0 $w2 $h2' role='img' aria-label='$(esc(env))'><rect width='100%' height='100%' fill='white'/><text x='$l2' y='14' font-size='12' font-weight='bold'>$(esc(env))</text>",
             "<line x1='$l2' y1='$(h2-b2)' x2='$(w2-8)' y2='$(h2-b2)' stroke='#888'/><line x1='$l2' y1='$t2' x2='$l2' y2='$(h2-b2)' stroke='#888'/>",
             @sprintf("<text x='%d' y='%.1f' font-size='9' text-anchor='end'>%.0f</text><text x='%d' y='%.1f' font-size='9' text-anchor='end'>0</text>", l2-4, py(ym)+3, ym, l2-4, py(0)+3)]
        for (ln, col, ys) in (("N", BLUE, c.N), ("C", ORANGE, c.C))
            any(==(ln), d.line) || continue
            push!(q, "<polyline points='" * join([@sprintf("%.1f,%.1f", px(tt), py(v)) for (tt, v) in zip(c.t, ys)], " ") * "' fill='none' stroke='$col' stroke-width='1.8'/>")
            for rr in eachrow(filter(x -> x.line == ln, d))
                push!(q, @sprintf("<circle cx='%.1f' cy='%.1f' r='2.6' fill='%s'/>", px(rr.day), py(rr.value), col))
            end
        end
        push!(q, "</svg>"); join(q)
    end
    envs = unique(TD.env)
    train_grid = join(["<div class='panel small'>$(tpanel(tw.id, e))</div>" for e in envs])
    params_list(r) = esc(replace(r.params, ";" => " · ", r"=(-?\d+\.\d{0,4})\d*" => s"=\1"))

    # nutrient level in the validation prediction (resource models only)
    nutrient_html = ""
    if "R" in names(VC)
        function rpanel(id, title)
            c = sort(filter(r -> r.id == id, VC), :t)
            w2, h2, l2, t2, b2 = 760, 200, 64, 26, 40
            px(v) = l2 + v/14*(w2-l2-20); py(v) = h2 - b2 - clamp(v, 0, 1.05)/1.05*(h2-t2-b2)
            q = ["<svg viewBox='0 0 $w2 $h2' role='img' aria-label='nutrient level'><rect width='100%' height='100%' fill='white'/><text x='$l2' y='16' font-size='13' font-weight='bold'>$(esc(title))</text>",
                 "<line x1='$l2' y1='$(h2-b2)' x2='$(w2-20)' y2='$(h2-b2)' stroke='#888'/><line x1='$l2' y1='$t2' x2='$l2' y2='$(h2-b2)' stroke='#888'/>",
                 @sprintf("<text x='%d' y='%.1f' font-size='10' text-anchor='end'>1 (fresh 10%% FBS)</text><text x='%d' y='%.1f' font-size='10' text-anchor='end'>0</text>", l2-4, py(1)+3, l2-4, py(0)+3)]
            for d in 0:2:14
                push!(q, @sprintf("<text x='%.1f' y='%d' font-size='10' text-anchor='middle'>%d</text>", px(d), h2-b2+14, d))
            end
            push!(q, "<polyline points='" * join([@sprintf("%.1f,%.1f", px(tt), py(v)) for (tt, v) in zip(c.t, c.R)], " ") * "' fill='none' stroke='#15803d' stroke-width='2'/>")
            push!(q, "<text x='$((l2+w2-20)÷2)' y='$(h2-8)' font-size='11' text-anchor='middle'>Day (medium replaced on days 2, 4, 6, 8, 10, 12)</text></svg>")
            join(q)
        end
        nutrient_html = "<div class='grid'><div class='panel'>$(rpanel(tw.id, "Nutrient level R, best training fit"))</div><div class='panel'>$(rpanel(vw.id, "Nutrient level R, best validation prediction"))</div></div>"
    end

    """
$heading
<div class='intro'><p>Each of $(nrow(R)) combinations of mechanisms is <b>fitted once, jointly, to all 22 training cultures</b> (untreated, low-resource and treated monocultures, Ce0 + Ce1 cocultures and low-resource cocultures; $(sum(TD.day .> 1)) points after day 1) with one shared parameter set in which resource level, dose and co-culture are conditions. Each fitted model then predicts this held-out culture with no refit. Training data and prediction both start at day 1.</p>
$intro
<p class='note'><b>Result:</b> the best fit to the training data ($(esc(tw.label)); $(tw.k) parameters) scores validation BIC $(@sprintf("%.2f", tw.val_bic)) (R² $(@sprintf("%.2f", tw.val_r2_sensitive)) sensitive, $(@sprintf("%.2f", tw.val_r2_resistant)) resistant; validation rank $(vrank(tw.id)) of $(nrow(R))). The best validation score is $(@sprintf("%.2f", vw.val_bic)) ($(esc(vw.label)); training rank $(findfirst(==(vw.id), R.id))). $(nbeat) of $(nrow(R)) fitted models beat the no-change reference (BIC $(@sprintf("%.2f", baseline.bic))). Rank correlation between training and validation performance: ρ = $(@sprintf("%.2f", ρs)).</p></div>
<div class='grid'><div class='panel'>$(panel("Best training fit → validation prediction", curve(tw.id)))</div><div class='panel'>$(panel("Best validation prediction (training rank $(findfirst(==(vw.id), R.id)))", curve(vw.id)))</div></div>
$nutrient_html
<div class='panel'><h3>Training fit vs validation, all $(nrow(R)) models</h3>$(join(sc))<p class='muted'>Green = best training fit; purple = best validation score; dashed line = no-change reference. A model generalises well if it is low on both axes.</p></div>
<h3>Best 15 by training fit</h3><div class='intro tscroll'><table>$head$top_train</table></div>
<h3>Best 15 by validation</h3><div class='intro tscroll'><table>$head$top_val</table><p class='muted'>Green rows beat the no-change reference. All results, including every fitted parameter value, are in <code>model_search_results.csv</code>.</p></div>
<h3>Each mechanism option</h3><div class='intro tscroll'><table><tr><th>Factor</th><th>Option</th><th>Best training ΔBIC</th><th>Best validation BIC</th><th>Models beating reference</th></tr>$(join(factor_rows))</table></div>
<h3>Best training model on every training culture</h3><div class='intro'><p class='muted'>Parameters: $(params_list(tw)).</p><div class='grid4'>$train_grid</div></div>
"""
end
quantile_cap(v) = (s = sort(v); max(s[clamp(round(Int, 0.9length(s)), 1, length(s))], 1.0))
search_html =
    search_section(joinpath(@__DIR__, "outputs", "new_data_resource_search"),
        "<h2>Resource-depletion model search (medium changed every other day)</h2>",
        """<p><b>Why:</b> the constant-effect models below cannot produce a rise followed by a crash. Here each culture's medium carries a nutrient level \\(R\\) (fresh 10% FBS medium = 1) that both lines consume; the medium is replaced on days 2, 4, 6, 8, 10 and 12, resetting \\(R\\). Low resource simply starts each fresh medium lower. Cisplatin is assumed re-added with each change (constant dose). Supported by work showing batch-culture nutrients run out within days and that serum starvation sensitises cancer cells to cisplatin.</p>
<p class='note'><b>Why no model can predict this culture from the training data:</b> no training culture combines low resource <i>and</i> cisplatin. Any term for how the two interact, such as the starvation–drug amplification \\(\\gamma\\), leaves the training error unchanged, so the training data cannot pin it down: across the fitted models its value spans the whole allowed range (0 to 10). The validation culture is the first time the combination appears, and there A2780 is suppressed far more than by either condition alone. Predicting it needs at least one training culture with both, or a calibration of \\(\\gamma\\) on part of the validation data.</p>
<div class='eq'>\\[\\frac{dN}{dt}=r_N\\,g(f)\\,N\\left(1-\\frac{N+\\alpha_{NC}C}{K}\\right)-\\delta_N(1-f)\\,N-\\beta_N A_N(t)H_N(D)\\left(1+\\gamma(1-f)\\right)N,\\qquad \\frac{dR}{dt}=-\\frac{q_N N+q_C C}{1000}\\,\\frac{R}{R+K_R}\\]</div>
<p class='muted'>\\(f=R(1+K_R)/(R+K_R)\\) is nutrient availability (1 in fresh normal medium). <b>Options:</b> crowding one-for-one or asymmetric × uptake the same per cell or per line × starvation kills (\\(g=1\\)), slows growth (\\(g=f\\), no \\(\\delta\\)) or both × delayed or immediate kill × drug independent of nutrients or amplified by starvation (\\(\\gamma\\)) × low-resource medium fixed at 5% of normal (0.5/10% FBS) or fitted. Fitting: bounded differential evolution (BlackBoxOptim), $(get(ENV, "EVALS", "50000")) evaluations per model.</p>""",
        [(:crowd, "Crowding", Dict("shared" => "shared K, one-for-one", "asym" => "shared K + asymmetric α")),
         (:consume, "Nutrient uptake", Dict("same" => "same per cell", "line" => "per line")),
         (:starve, "What starvation does", Dict("death" => "kills", "growth" => "slows growth", "both" => "slows growth and kills")),
         (:drug, "Cisplatin timing", Dict("delay" => "delayed kill", "now" => "immediate kill")),
         (:interaction, "Drug × nutrients", Dict("indep" => "independent", "amp" => "starvation amplifies drug")),
         (:r0, "Low-resource medium", Dict("fixed" => "5% of normal (fixed)", "fitted" => "fitted"))]) *
    search_section(joinpath(@__DIR__, "outputs", "new_data_model_search"),
        "<h2>Constant-effect model search, fitted on all training data</h2>",
        """<p><b>Mechanism options:</b> crowding (shared K one-for-one, shared K with asymmetric α, or separate K per line with α) × low resource (extra death, slower growth \\(r(1-\\rho u_R)\\), lower capacity \\(K(1-\\kappa u_R)\\), or both) × cisplatin (delayed kill, immediate kill, or delayed growth arrest; Hill dose response with EC50 = 1 µM) × drug–resource interaction (independent, low resource amplifies the drug by \\(1+\\gamma u_R\\), or kills dividing cells only) × baseline death (none or \\(d_N, d_C\\)). Low resource acts as a constant condition. Fitting: bounded BFGS for all, then Nelder–Mead polishing of the 30 best by training BIC.</p>""",
        [(:crowd, "Crowding", Dict("shared" => "shared K, one-for-one", "asym" => "shared K + asymmetric α", "lv" => "separate K per line + α")),
         (:low_resource, "How low resource acts", Dict("death" => "extra death", "rate" => "slower growth", "cap" => "lower K", "ratecap" => "slower growth + lower K")),
         (:drug, "How cisplatin acts", Dict("delay" => "delayed kill", "now" => "immediate kill", "arrest" => "delayed growth arrest")),
         (:interaction, "Drug × resource", Dict("indep" => "independent", "sens" => "low resource amplifies drug", "prolif" => "kills dividing cells only")),
         (:death, "Baseline death", Dict("none" => "none", "base" => "d_N, d_C"))])

eq(s) = "<div class='eq'>$s</div>"

hd_html = isempty(hd_combos) ? "" : begin
    hb = hd_combos[1]
    rows = join([@sprintf("<tr%s><td class='num'>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td class='num'>%.2f%s</td><td class='num'>%.2f</td><td class='num'>%.1f</td><td class='num'>%.1f</td><td class='num'>%.2f</td><td class='num'>%.2f</td></tr>",
        i == 1 ? " class='win'" : "", i, esc(c.cap.label), esc(c.dr.label), esc(c.lr.label), esc(c.pool.label), c.sc.bic, c.sc.bic < baseline.bic ? " ✓" : "",
        c.sc.bic - best.sc.bic, c.sc.rmseN, c.sc.rmseC, c.sc.r2N, c.sc.r2C) for (i, c) in enumerate(hd_combos)], "\n")
    src = join(["<li><b>$(esc(s.label))</b>: $(esc(s.model))</li>" for s in vcat(drug_src, lr_src)])
    """<h2>Healthy / damaged prediction</h2>
<div class='intro'><p>The healthy / damaged family from the fitting report, again with <b>nothing refitted</b>. In this culture both stresses act at once, so each line has a cisplatin damage channel (rate \\(\\kappa_X H_X(1\\,\\mu\\text{M})\\)) and a low-resource damage channel (rate \\(\\sigma_X u_R\\)), each with the repair and death rates it was fitted with. \\(H_X\\) here also carries any damage build-up \\(A(t)\\) and coculture damage \\(b_X\\) of the source fit, and \\(\\sigma_X u_R\\) is multiplied by crowding when the source model is density-dependent:</p>
$(eq(raw"\[\frac{dX_H}{dt}=r_X X_H\left(1-\frac{X+\alpha_{XY}Y}{K}\right)-\left(\kappa_X H_X+\sigma_X u_R\right)X_H+\rho_X^{\text{cis}}X_D^{\text{cis}}+\rho_X^{LR}X_D^{LR}\]"))
$(eq(raw"\[\frac{dX_D^{\text{cis}}}{dt}=\kappa_X H_X X_H-(\rho_X^{\text{cis}}+\delta_X^{\text{cis}})X_D^{\text{cis}},\qquad \frac{dX_D^{LR}}{dt}=\sigma_X u_R X_H-(\rho_X^{LR}+\delta_X^{LR})X_D^{LR},\qquad X=X_H+X_D^{\text{cis}}+X_D^{LR}\]"))
<p>Options: capacity K, α from the Stage 3 or 3b winner; cisplatin damage from the Stage 2 or Stage 3 winner; starvation damage from the Stage 1b or 3b winner; and either <b>separate damage pools</b> (above) or <b>one pool</b>, where both stresses damage cells into a single compartment with the cisplatin repair and death rates. Source models:</p><ul>$src</ul>
<p><b>Best healthy/damaged prediction:</b> $(esc(hb.name)) — BIC $(@sprintf("%.2f", hb.sc.bic)) (ΔBIC $(@sprintf("%.2f", hb.sc.bic - best.sc.bic)) against the best one-population combination above; no-change reference $(@sprintf("%.2f", baseline.bic))).</p></div>
<div class='grid'><div class='panel'>$(panel("Best healthy/damaged: " * hb.name, hb.sol))</div><div class='panel'>$(panel("Best one-population combination: " * best.name, best.sol))</div></div>
<div class='intro tscroll'><table><tr><th>Rank</th><th>Capacity</th><th>Cisplatin damage</th><th>Starvation damage</th><th>Damage pools</th><th>BIC</th><th>ΔBIC vs best one-population</th><th>RMSE sensitive</th><th>RMSE resistant</th><th>R² sensitive</th><th>R² resistant</th></tr>$rows</table>
<p class='muted'>✓ = beats the no-change reference. ΔBIC is against the best of the $(length(combos)) one-population combinations (negative = better).</p></div>"""
end

best_line =@sprintf("%s — BIC %.2f, RMSE %.0f (sensitive) and %.0f (resistant) cells, R² %.2f and %.2f", best.name, best.sc.bic, best.sc.rmseN, best.sc.rmseC, best.sc.r2N, best.sc.r2C)
html = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Validation: low-resource + cisplatin coculture</title>
<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}};</script><script defer src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>
<style>body{font:14px/1.5 system-ui,sans-serif;margin:0;background:#f4f7f8;color:#17212b}main{max-width:1400px;margin:auto;padding:24px}h1{color:#124b6e;margin:6px 0}h2{color:#124b6e;margin-top:34px;border-top:2px solid #cfdde3;padding-top:14px}h3{margin:14px 0 6px}
.topnav{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:0 0 16px}.nav{display:inline-flex;align-items:center;min-height:34px;padding:0 11px;border-radius:8px;color:#176b87;background:#eef5f6;text-decoration:none;font:650 14px system-ui,sans-serif;white-space:nowrap}.nav-home{padding:6px 9px}.nav-active{color:#fff;background:#176b87}
.intro,.panel{background:#fff;border:1px solid #d9e2e6;border-radius:12px;padding:16px 18px;margin:14px 0;min-width:0;overflow-wrap:anywhere}.intro{border-left:6px solid #176b87}.note{background:#fff6d9;border-radius:8px;padding:10px 14px}.eq{overflow-x:auto;overflow-y:hidden;background:#fbfdfd;border:1px solid #e4ecef;border-radius:8px;padding:4px 12px;margin:8px 0}
.tscroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}td,th{border:1px solid #d3dde2;padding:6px 8px;text-align:left;vertical-align:top}th{background:#e9f3f8}td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}tr.win td{background:#eef8f1;font-weight:600}tr.star td{background:#fff8e6}tr.ref td{background:#f1f3f5;font-style:italic}
table.pairs th:first-child{white-space:nowrap}table.pairs td{text-align:center;font-weight:600}
.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}.grid4{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px}.panel.small{padding:6px;margin:0}@media(max-width:1100px){.grid4{grid-template-columns:repeat(2,minmax(0,1fr))}}svg{display:block;width:100%;height:auto}code{overflow-wrap:anywhere}.muted{color:#5b6b73}dl{display:grid;grid-template-columns:max-content 1fr;gap:6px 14px}dt{font-weight:700}
@media(max-width:900px){.grid{grid-template-columns:1fr}}@media(max-width:640px){main{padding:14px}}</style></head><body><main>
$nav
<h1>Validation — low-resource + cisplatin coculture</h1>
<div class='intro'><p>Generated $(Dates.format(now(), "yyyy-mm-dd HH:MM")). <b>No parameters are fitted here.</b> A2780 and A2780cis were co-cultured 1:1 at 30k in low resource (0.5% FBS, \\(u_R = $(round(U_R; digits=3))\\)) with cisplatin, assumed at the 1 µM IC50 used for Ce1. This condition was never used for fitting. Every prediction combines parameters already estimated in the <a href='../new_data_fits/report.html'>New data fitting</a> report; what varies is <i>how</i> the ingredients are assumed to act and interact. Each prediction starts from the <b>day-1</b> counts (N = $(fmt(u0[1])), C = $(fmt(u0[2]))) and is scored on days 2–14: the day-0 values ($(fmt(u0_seed[1])) each, SD 0) are the nominal seeding count, and the measured counts jump $(round(Int, u0[1] / u0_seed[1]))-fold (A2780) and $(round(Int, u0[2] / u0_seed[2]))-fold (A2780cis) by day 1, beyond any fitted growth rate.</p>
<p>Files: <code>New Datasets/Validation/$(esc(FILE_N))</code> and <code>$(esc(FILE_C))</code> (mean of 3 wells per day).</p>
<p class='note'>Scoring: with no free parameters, BIC = n ln(SSE/n), so ranking by BIC is ranking by squared error. R² compares with a flat line at each population's mean (below 0 = worse than that line). The naive reference "both populations stay at their day-1 counts" scores BIC $(@sprintf("%.2f", baseline.bic)); a ✓ marks combinations that beat it. The same exclusion rule as the fitting report removes $(count(bN) + count(bC)) transient drops (marked ×).</p></div>

<h2>How the terms can interact</h2>
<div class='intro'>
$(eq(raw"\[\frac{dN}{dt}=\underbrace{\rho_r\,r_N\,g_{\text{drug}}\,N\left(1-\frac{N+\alpha_{NC}C}{\rho_K K}\right)}_{\text{growth}}-\underbrace{\left(d_N+d_N^{LR}u_R\right)N}_{\text{death}}-\underbrace{\beta_N A_N(t)\,H_N(D)\,s_N\,N}_{\text{drug kill}}\]"))
<p>The same form holds for C with its own parameters. Four ingredients are switched between options, each option using values from an earlier stage:</p>
<dl>
<dt>Capacity and crowding</dt><dd>K and α from the Stage 3 or Stage 3b fits (optimal or shared-K-only); the Stage 3 optimal also brings its sensitive death rate d_N.</dd>
<dt>How low resource acts</dt><dd>not at all · as an extra death rate \\(d^{LR}u_R\\) · as a lower growth rate and capacity (\\(\\rho_r, \\rho_K\\)) · as a lower growth rate only · as an immediate loss plus extra loss that accumulates through 8 damage steps (<b>delayed death</b>; the best Stage 1b model for A2780cis). These are Stage 1b candidates.</dd>
<dt>How cisplatin acts</dt><dd>not at all · kills after a delay \\(A(t)=1-e^{-\\lambda t}\\) (Stage 2 or Stage 3 fit) · kills from day 0 (\\(A=1\\)) · <b>stops growth</b> instead of killing: \\(g_{\\text{drug}}=\\max(0,\\,1-\\beta A H/r)\\) and no kill term.</dd>
<dt>Drug × growth interaction</dt><dd><b>independent</b> (\\(s=1\\)): the kill rate is the same whatever the conditions · <b>dividing cells only</b> (\\(s = \\max(0, 1-\\text{crowding})\\)): the drug only kills cells with room to divide, so crowding and starvation protect cells. This is the usual picture for cisplatin.</dd>
</dl></div>

<h2>Result</h2>
<div class='intro'><p><b>Best of $(length(combos)) combinations:</b> $(esc(best_line)).</p>
<p>$(n_beat) of $(length(combos)) combinations beat the no-change reference (BIC $(@sprintf("%.2f", baseline.bic))). The Candidate models page's low-resource treated coculture model (★: Stage 3 K and α, extra low-resource death, delayed Stage 2 kill, independent) ranks $(star_rank) (ΔBIC $(@sprintf("%.1f", starc.sc.bic - best.sc.bic))).</p></div>
<div class='grid'><div class='panel'>$(panel("Best: " * best.name, best.sol))</div><div class='panel'>$(panel("★ Candidate-models combination (rank $(star_rank))", starc.sol))</div></div>

$hd_html

<h2>Which pairs of mechanisms work together?</h2>
<div class='intro'><p>Each cell is the best ΔBIC reached by any combination using that pair (green = close to the best, red = 60 or more worse; ✓ = beats the no-change reference).</p>
<h3>Low resource × cisplatin</h3>$(pair_table(:lr, :drug))
<h3>Cisplatin × interaction</h3>$(pair_table(:drug, :inter))
<h3>Capacity × low resource</h3>$(pair_table(:cap, :lr))</div>

<h2>Each ingredient on its own</h2>
<div class='intro'><p>Best and mean ΔBIC over every combination that uses each option.</p>$marg_html</div>

<h2>Top $(TOPN) combinations</h2>
<div class='intro tscroll'><table><tr><th>Rank</th><th>Capacity</th><th>Low resource</th><th>Cisplatin</th><th>Interaction</th><th>BIC</th><th>ΔBIC</th><th>RMSE sensitive</th><th>RMSE resistant</th><th>R² sensitive</th><th>R² resistant</th><th>BIC, all points</th><th>BIC, from day-0 seeding</th></tr>$rank_rows</table>
<p class='muted'>All $(length(combos)) combinations are in <code>validation_ranking.csv</code>.</p></div>

<h2>Next five combinations</h2>
<div class='grid'>$(join(["<div class='panel'>$(panel(@sprintf("#%d · ΔBIC %.1f · %s", i + 1, c.sc.bic - best.sc.bic, c.name), c.sol))</div>" for (i, c) in enumerate(combos[2:6])]))</div>

$search_html

<h2>Parameter sources</h2>
<div class='intro tscroll'><table><tr><th>Ingredient</th><th>Values (from the staged fits)</th></tr>$src_rows</table>
<p class='muted'>The Candidate models page also lists a resource-driven sensitive-to-resistant transfer rate γ. It was never fitted, so it cannot be tested without fitting and is left out.</p></div>

<h2>Validation data</h2>
<div class='intro tscroll'><table><tr><th>Day</th><th>A2780 mean cells</th><th>A2780cis mean cells</th></tr>$data_rows</table></div>
<p class='muted'>Outputs: <code>validation_ranking.csv</code> in <code>Modeling_Approaches/02_mechanical_automatic_package/outputs/new_data_validation/</code>. Generated by <code>run_new_data_validation.jl</code> from <code>all_fitted_parameters.csv</code>.</p>
</main></body></html>"""

write(joinpath(OUT_DIR, "report.html"), html)
cp(joinpath(OUT_DIR, "report.html"), joinpath(DOCS_DIR, "report.html"); force=true)
println("Wrote $(joinpath(DOCS_DIR, "report.html"))")
