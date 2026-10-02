# Optimal-control experiments on the selected compartment models (control_model.jl, control_mpc.jl).
# Writes outputs/control/control_state.jls (checkpoint; REDO=1 recomputes) and CSVs; control_report.jl
# renders the report.   julia -t 12 --project=. run_control.jl
#
#   environments  normal medium (Stage 3 model), low resource (Stage 4 model)
#   scenarios     A (selected coculture drug response), B (monoculture drug response for both lines)
#   ratios        25:75, 50:50, 75:25 sensitive:resistant, started from the day-0 counts
#   policies      no treatment, constant 1 µM, maximum measured dose, periodic, threshold adaptive,
#                 MPC data-supported (measured doses only), MPC demonstration (continuous, up to u_max_demo)
#   sweeps        weights, toxicity (r_H, s_B, B_tox), horizon, decision interval, integration step,
#                 28-day extrapolation
#   robustness    every schedule re-simulated on 200 bootstrap parameter sets; closed-loop MPC on 40

include(joinpath(@__DIR__, "control_model.jl"))
include(joinpath(@__DIR__, "control_mpc.jl"))

# Run length: CONTROL_DAYS (default closed_loop_days = 14). Longer runs (e.g. CONTROL_DAYS=30) keep the
# 14-day look-ahead, re-planned at every media change, and write to outputs/control_<days>d/.
const DAYS = parse(Float64, get(ENV, "CONTROL_DAYS", string(CFG["schedule"]["closed_loop_days"])))
const COUT = DAYS == Float64(CFG["schedule"]["closed_loop_days"]) ? joinpath(@__DIR__, "outputs", "control") :
             joinpath(@__DIR__, "outputs", "control_$(round(Int, DAYS))d"); mkpath(COUT)
const CKP = joinpath(COUT, "control_state.jls")
cck = (isfile(CKP) && get(ENV, "REDO", "0") != "1") ? deserialize(CKP) : Dict{Symbol,Any}()
savec() = serialize(CKP, cck)

const Δ0 = Float64(CFG["schedule"]["decision_interval_days"])
const NDEC = round(Int, CFG["schedule"]["horizon_days"] / Δ0)    # look-ahead (decisions)
const NRUN = round(Int, DAYS / Δ0)                                  # decisions over the whole run
const LEVELS = Float64.(CFG["dose"]["levels_data"])
const UMAX_DATA = Float64(CFG["dose"]["u_max_data"]); const UMAX_DEMO = Float64(CFG["dose"]["u_max_demo"])
const BOOT_N = min(200, size(ZOO_STATE.results["s3"].res.pars[1], 1))

ratios_normal = [(ratio=("25:75", "50:50", "75:25")[i], N0=co_data[i].S[2][1], C0=co_data[i].R[2][1]) for i in 1:3]      # Ce0 day-0 counts
ratios_lr = [(ratio=Dict("25S75R" => "25:75", "50S50R" => "50:50", "75S25R" => "75:25")[c.label], N0=c.yN[1], C0=c.yC[1]) for c in lr_cc if c.run == 1]
plant_for(env, scen; b=nothing, tox=default_tox()) = env == :normal ? scenario_plant(scen; b=b, tox=tox) : lr_plant(scen; b=b, tox=tox)
ratios_for(env) = env == :normal ? ratios_normal : ratios_lr
println("=== optimal control: $(Threads.nthreads()) threads; Δ = $Δ0 d, $NRUN decisions ($DAYS days, look-ahead $NDEC); checkpoint has $(collect(keys(cck))) ==="); flush(stdout)

function run_all_policies(P, x0; Δ=Δ0, n=NRUN, H=NDEC)
    W = default_weights(P); Wd = default_weights(P; umax=UMAX_DEMO)
    out = Dict{String,Any}()
    for (name, pol) in baseline_policies(P, W)
        U = run_policy(pol, P, x0, Δ, n); out[name] = (U=U, o=outcomes(P, x0, U, Δ, W))
    end
    t = time(); cl = closed_loop(P, P, x0, Δ, n, H, W; mode=:discrete, levels=LEVELS)
    out["MPC, measured doses only"] = (U=cl.U, o=outcomes(P, x0, cl.U, Δ, W), secs=time() - t)
    t = time(); cl = closed_loop(P, P, x0, Δ, n, H, Wd; mode=:continuous)
    out["MPC demonstration (continuous, ≤ $(UMAX_DEMO) µM)"] = (U=cl.U, o=outcomes(P, x0, cl.U, Δ, W), secs=time() - t, iters=cl.iters)
    out
end

# ── main runs: environment × scenario × ratio ────────────────────────────────
if !haskey(cck, :main)
    tasks = [(env, scen, r) for env in (:normal, :lowres) for scen in (:A, :B) for r in ratios_for(env)]
    out = Vector{Any}(undef, length(tasks))
    Threads.@threads :dynamic for i in eachindex(tasks)
        env, scen, r = tasks[i]; t = time()
        P = plant_for(env, scen); x0 = x0_for(P, r.N0, r.C0)
        out[i] = (env=env, scen=scen, ratio=r.ratio, x0=x0, name=P.name, K=P.K, Bc=default_weights(P).Bc, policies=run_all_policies(P, x0))
        @printf("  %-7s %s %s done [%.0fs]\n", env, scen, r.ratio, time() - t); flush(stdout)
    end
    cck[:main] = out; savec()
end

# ── sweeps (normal medium, 50:50, data-supported MPC) ─────────────────────────
function mpc_outcome(P, x0, W; Δ=Δ0, n=NRUN, H=NDEC, mode=:discrete)
    cl = closed_loop(P, P, x0, Δ, n, H, W; mode=mode, levels=LEVELS)
    (U=cl.U, o=outcomes(P, x0, cl.U, Δ, W))
end
if !haskey(cck, :sweeps)
    r = ratios_normal[2]
    tasks = Any[]
    for scen in (:A, :B)
        P = plant_for(:normal, scen); x0 = x0_for(P, r.N0, r.C0)
        for w in CFG["weights"]["sweep"]
            push!(tasks, (kind="weights", scen=scen, label="qN=$(w[1]) qC=$(w[2]) qT=$(w[3]) ru=$(w[4])", f=() -> mpc_outcome(P, x0, default_weights(P; w=Float64.(w)))))
        end
        base = default_tox()
        for v in CFG["toxicity"]["r_H_sweep"]
            Pt = plant_for(:normal, scen; tox=merge(base, (r_H=Float64(v),)))
            push!(tasks, (kind="toxicity", scen=scen, label="r_H = $v /day", f=() -> mpc_outcome(Pt, x0, default_weights(Pt))))
        end
        for v in CFG["toxicity"]["s_B_sweep"]
            Pt = plant_for(:normal, scen; tox=merge(base, (s_B=Float64(v),)))
            push!(tasks, (kind="toxicity", scen=scen, label="s_B = $v per µM", f=() -> mpc_outcome(Pt, x0, default_weights(Pt))))
        end
        for v in CFG["toxicity"]["B_tox_sweep"]
            push!(tasks, (kind="toxicity", scen=scen, label="B_tox = $v", f=() -> mpc_outcome(P, x0, default_weights(P; Btox=Float64(v)))))
        end
        for Hh in (1, 2, 4, NDEC)
            push!(tasks, (kind="horizon", scen=scen, label="$(Hh) decisions ($(Hh * Δ0) days)", f=() -> mpc_outcome(P, x0, default_weights(P); H=Hh)))
        end
        push!(tasks, (kind="interval", scen=scen, label="decide every 1 day (7-day horizon)", f=() -> mpc_outcome(P, x0, default_weights(P); Δ=1.0, n=2NRUN, H=7)))
        push!(tasks, (kind="interval", scen=scen, label="decide every 2 days (14-day horizon)", f=() -> mpc_outcome(P, x0, default_weights(P))))
        # the 28-day extrapolation row belongs to the 14-day report; longer reports are extrapolations throughout
        NRUN == NDEC && push!(tasks, (kind="extrapolation", scen=scen, label="28 days, receding 14-day horizon (beyond the data)", f=() -> mpc_outcome(P, x0, default_weights(P); n=14, H=NDEC)))
    end
    out = Vector{Any}(undef, length(tasks))
    Threads.@threads :dynamic for i in eachindex(tasks)
        T = tasks[i]; out[i] = (kind=T.kind, scen=T.scen, label=T.label, res=T.f())
    end
    # integration-step check: the nominal data-supported MPC schedule re-simulated with smaller RK4 steps
    steps = NamedTuple[]
    for scen in (:A, :B)
        m = only(filter(x -> x.env == :normal && x.scen == scen && x.ratio == "50:50", cck[:main]))
        P = plant_for(:normal, scen); U = m.policies["MPC, measured doses only"].U
        ref = nothing
        for h in (0.1, 0.05, 0.02, 0.01, 0.005)
            x = copy(m.x0)
            for u in U; x = F(x, u, P, Δ0; h=h); end
            ref === nothing && (ref = x)
            push!(steps, (scen=scen, h=h, N=x[P.iN], C=x[P.iC], T=x[P.iT]))
        end
    end
    cck[:sweeps] = out; cck[:steps] = steps; savec()
    println("  sweeps done"); flush(stdout)
end

# ── robustness over bootstrap parameter sets ──────────────────────────────────
if !haskey(cck, :robust)
    rob = Dict{Any,Any}()
    for env in (:normal, :lowres), scen in (:A, :B)
        m = only(filter(x -> x.env == env && x.scen == scen && x.ratio == "50:50", cck[:main]))
        P0 = plant_for(env, scen); W0 = default_weights(P0)
        pols = ["No treatment", "Constant 1 µM", "Periodic 1 µM, 2 days on / 2 off", "Threshold adaptive (Zhang-style)", "MPC, measured doses only"]
        nb = env == :normal ? BOOT_N : min(BOOT_N, size(ZOO_STATE.results["s4"].res.pars[1], 1))
        draws = Vector{Any}(undef, nb)
        Threads.@threads :dynamic for b in 1:nb
            Pb = try plant_for(env, scen; b=b) catch; nothing end
            if Pb === nothing
                draws[b] = nothing; continue
            end
            Wb = default_weights(Pb); Wb = Weights(Wb.qN, Wb.qC, Wb.qT, Wb.ru, Wb.tf, Wb.μC, Wb.μT, Wb.δ, W0.Bc, Wb.Btox, Wb.umax, Wb.K)   # same bound as the nominal plan
            d = Dict{String,Any}()
            for pn in pols
                U = pn == "Threshold adaptive (Zhang-style)" ? run_policy(threshold_policy(Pb, Wb), Pb, m.x0, Δ0, NRUN) : m.policies[pn].U
                o = try outcomes(Pb, m.x0, U, Δ0, Wb) catch; nothing end
                d[pn] = o === nothing ? nothing : (N=o.N, C=o.C, T=o.T, cum=o.cum_dose, cv=o.cancer_violation, tv=o.tox_violation, Nf=o.N_final, Cf=o.C_final, tot=o.total_final, rf=o.resistant_fraction)
            end
            draws[b] = d
        end
        # closed-loop MPC (plans with the nominal model, runs on a bootstrap plant) on a subset
        ncl = min(40, nb); cl = Vector{Any}(undef, ncl)
        Threads.@threads :dynamic for b in 1:ncl
            Pb = try plant_for(env, scen; b=b) catch; nothing end
            cl[b] = Pb === nothing ? nothing : (c = closed_loop(P0, Pb, m.x0, Δ0, NRUN, NDEC, W0; mode=:discrete, levels=LEVELS);
                                              o = outcomes(Pb, m.x0, c.U, Δ0, W0); (U=c.U, N=o.N, C=o.C, T=o.T, cum=o.cum_dose, cv=o.cancer_violation, tv=o.tox_violation, tot=o.total_final, rf=o.resistant_fraction))
        end
        rob[(env, scen)] = (ts=m.policies["No treatment"].o.ts, draws=[d for d in draws if d !== nothing], closed=[c for c in cl if c !== nothing])
        @printf("  robustness %-7s %s: %d draws, %d closed-loop\n", env, scen, length(rob[(env, scen)].draws), length(rob[(env, scen)].closed)); flush(stdout)
    end
    cck[:robust] = rob; savec()
end

# ── CSV of every main-run outcome ─────────────────────────────────────────────
rows = NamedTuple[]
for m in cck[:main], (pn, v) in m.policies
    o = v.o
    push!(rows, (environment=string(m.env), scenario=string(m.scen), ratio=m.ratio, policy=pn, schedule_uM=join(round.(v.U; digits=3), " "),
                 N_final=o.N_final, C_final=o.C_final, total_final=o.total_final, resistant_fraction=o.resistant_fraction,
                 max_total=o.max_total, max_toxicity=o.max_tox, cumulative_dose_uM_days=o.cum_dose,
                 cancer_bound_violated=o.cancer_violation, toxicity_bound_violated=o.tox_violation))
end
CSV.write(joinpath(COUT, "control_outcomes.csv"), DataFrame(rows))
println("Wrote $(joinpath(COUT, "control_outcomes.csv"))")
