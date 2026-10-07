# Optimal-control experiments on the selected compartment models (control_model.jl, control_mpc.jl).
# Writes outputs/control/control_state.jls (checkpoint; REDO=1 recomputes) and CSVs; control_report.jl
# renders the report.   julia -t 12 --project=. run_control.jl
#
#   environments  normal medium (Stage 3 model), low resource (Stage 4 model)
#   scenarios     A (selected coculture drug response), B (monoculture drug response for both lines)
#   ratios        25:75, 50:50, 75:25 sensitive:resistant, started from the NOMINAL seeding
#                 (control_config.toml [initial]: 100 cells split by the ratio), not a measured day-0 count
#   policies      no treatment, constant 1 µM, maximum measured dose, periodic, threshold adaptive,
#                 MPC data-supported (measured doses only), MPC demonstration (continuous, up to u_max_demo)
#                 MPC uses a FIXED 14-day look-ahead at every decision, also at the end of the run
#   sweeps        weights, soft burden penalty (B_cancer) × terminal factor, toxicity (r_H, s_B, B_tox),
#                 horizon, decision interval, integration step, 28-day extrapolation
#   robustness    every schedule re-simulated on the open-loop bootstrap draws; closed-loop MPC (observes
#                 N and C with measurement noise, reconstructs the hidden states with its own model) on the
#                 first closed-loop draws, where every policy is compared on the same draws and noise

include(joinpath(@__DIR__, "control_fasting.jl"))     # includes control_model.jl and control_mpc.jl; lr_plant lives there

# Run length: CONTROL_DAYS (default closed_loop_days = 14). Longer runs keep the 14-day look-ahead,
# re-planned at every media change, and write to outputs/control_<days>d/ (the long dose + fasting run is
# run_control_fasting.jl).
const DAYS = parse(Float64, get(ENV, "CONTROL_DAYS", string(CFG["schedule"]["closed_loop_days"])))
# CONTROL_OUT replaces the base folder "outputs" (a side run that must not overwrite the main results)
const COUT = joinpath(get(ENV, "CONTROL_OUT", joinpath(@__DIR__, "outputs")),
                      DAYS == Float64(CFG["schedule"]["closed_loop_days"]) ? "control" : "control_$(round(Int, DAYS))d"); mkpath(COUT)
const CKP = joinpath(COUT, "control_state.jls")
cck = (isfile(CKP) && get(ENV, "REDO", "0") != "1") ? deserialize(CKP) : Dict{Symbol,Any}()
savec() = serialize(CKP, cck)

const Δ0 = Float64(CFG["schedule"]["decision_interval_days"])
const NDEC = round(Int, CFG["schedule"]["horizon_days"] / Δ0)    # look-ahead (decisions), the same at every decision
const NRUN = round(Int, DAYS / Δ0)                                  # decisions over the whole run
const LEVELS = Float64.(CFG["dose"]["levels_data"])
const UMAX_DATA = Float64(CFG["dose"]["u_max_data"]); const UMAX_DEMO = Float64(CFG["dose"]["u_max_demo"])
const BOOT_N = min(Int(CFG["robustness"]["open_loop_draws"]), size(ZOO_STATE.results["s3"].res.pars[1], 1))
const NCL = min(Int(CFG["robustness"]["closed_loop_draws"]), BOOT_N)   # closed-loop draws (paired with every policy)
const CL_NAME = "Closed-loop MPC (observes N, C with noise)"

# starting counts: the nominal seeding split by the ratio (control_model.jl: nominal_seeding), the same in
# both environments; the data's day-0 values are not used
ratios_normal = [merge((ratio=r,), nominal_seeding(r)) for r in ("25:75", "50:50", "75:25")]
ratios_lr = ratios_normal
plant_for(env, scen; b=nothing, tox=default_tox()) = env == :normal ? scenario_plant(scen; b=b, tox=tox) : lr_plant(scen; b=b, tox=tox)
ratios_for(env) = env == :normal ? ratios_normal : ratios_lr
println("=== optimal control: $(Threads.nthreads()) threads; Δ = $Δ0 d, $NRUN decisions ($DAYS days, fixed look-ahead $NDEC); checkpoint has $(collect(keys(cck))) ==="); flush(stdout)

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
# Wo: the weights the outcomes are scored with (bound counts always at the default B_cancer)
function mpc_outcome(P, x0, W; Wo=W, Δ=Δ0, n=NRUN, H=NDEC, mode=:discrete)
    cl = closed_loop(P, P, x0, Δ, n, H, W; mode=mode, levels=LEVELS)
    (U=cl.U, o=outcomes(P, x0, cl.U, Δ, Wo))
end
if !haskey(cck, :sweeps)
    r = ratios_normal[2]
    tasks = Any[]
    for scen in (:A, :B)
        P = plant_for(:normal, scen); x0 = x0_for(P, r.N0, r.C0); W0 = default_weights(P)
        for w in CFG["weights"]["sweep"]
            push!(tasks, (kind="weights", scen=scen, label="qB=$(w[1]) qR=$(w[2]) qT=$(w[3]) ru=$(w[4])", f=() -> mpc_outcome(P, x0, default_weights(P; w=Float64.(w)))))
        end
        # soft burden penalty: where it starts (0 = none) × the terminal weight; scored against the default B_cancer
        for bc in CFG["cancer"]["B_cancer_sweep"], tf in CFG["cancer"]["terminal_factor_sweep"]
            lab = (bc > 0 ? "B_cancer = $(bc) K" : "no burden penalty") * ", terminal ×$(tf)"
            push!(tasks, (kind="burden penalty", scen=scen, label=lab, f=() -> mpc_outcome(P, x0, default_weights(P; bc_frac=Float64(bc), tf=Float64(tf)); Wo=W0)))
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
# Every draw b: the fixed schedules (planned once with the best-fit model) are replayed on plant b; the
# threshold rule re-decides on plant b from noisy counts; on the first NCL draws the closed-loop MPC plans
# with the best-fit model, observes plant b's N and C with the same noise, and reconstructs the hidden
# states with its own model. Paired: on those draws every policy meets the same plant and the same noise.
if !haskey(cck, :robust)
    rob = Dict{Any,Any}()
    for env in (:normal, :lowres), scen in (:A, :B)
        m = only(filter(x -> x.env == env && x.scen == scen && x.ratio == "50:50", cck[:main]))
        P0 = plant_for(env, scen); W0 = default_weights(P0)
        pols = ["No treatment", "Constant 1 µM", "Maximum measured dose (1.47 µM)", "Periodic 1 µM, 2 days on / 2 off", "Threshold adaptive (Zhang-style)", "MPC, measured doses only"]
        nb = env == :normal ? BOOT_N : min(BOOT_N, size(ZOO_STATE.results["s4"].res.pars[1], 1))
        draws = Vector{Any}(undef, nb)
        summ(o) = (N=o.N, C=o.C, T=o.T, cum=o.cum_dose, cv=o.cancer_violation, tv=o.tox_violation, dover=o.days_over, Nf=o.N_final, Cf=o.C_final, tot=o.total_final, rf=o.resistant_fraction, maxT=o.max_tox)
        Threads.@threads :dynamic for b in 1:nb
            Pb = try plant_for(env, scen; b=b) catch; nothing end
            if Pb === nothing
                draws[b] = nothing; continue
            end
            Wb = default_weights(Pb); Wb = Weights(Wb.qB, Wb.qR, Wb.qT, Wb.ru, Wb.tf, Wb.μC, Wb.μT, Wb.δ, W0.Bc, Wb.Btox, Wb.umax, Wb.K)   # same bound as the nominal plan
            noise = obs_noise(b, NRUN)
            d = Dict{String,Any}("b" => b)
            for pn in pols
                o = try
                    U = pn == "Threshold adaptive (Zhang-style)" ? run_policy(threshold_policy(Pb, Wb), Pb, m.x0, Δ0, NRUN; noise=noise) : m.policies[pn].U
                    outcomes(Pb, m.x0, U, Δ0, Wb)
                catch; nothing end
                d[pn] = o === nothing ? nothing : summ(o)
            end
            if b <= NCL
                c = try closed_loop(P0, Pb, m.x0, Δ0, NRUN, NDEC, W0; mode=:discrete, levels=LEVELS, noise=noise) catch; nothing end
                d[CL_NAME] = c === nothing ? nothing : merge(summ(outcomes(Pb, m.x0, c.U, Δ0, Wb)), (U=c.U,))
            end
            draws[b] = d
        end
        ds = [d for d in draws if d !== nothing]
        rob[(env, scen)] = (ts=m.policies["No treatment"].o.ts, draws=ds,
                            closed=[d[CL_NAME] for d in ds if get(d, CL_NAME, nothing) !== nothing], ncl=NCL, cv=Float64(CFG["observation"]["count_cv"]))
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
                 cancer_bound_violated=o.cancer_violation, days_above_B_cancer=o.days_over, toxicity_bound_violated=o.tox_violation))
end
CSV.write(joinpath(COUT, "control_outcomes.csv"), DataFrame(rows))
println("Wrote $(joinpath(COUT, "control_outcomes.csv"))")
