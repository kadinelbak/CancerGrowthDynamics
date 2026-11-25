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
#                 MPC looks ahead lookahead_blocks (move blocking; [1,1,2,2,4,4] = 28 days in 6 decisions);
#                 with control_config lookahead_to_end (default since 2026-10-27) the look-ahead stops at the
#                 end of the run, otherwise it is fixed at every decision (CONTROL_LOOKAHEAD_TO_END=0/1 overrides)
#   objective     time to progression (control_config [progression]): time spent at N + C ≥ B_prog (K/2) costs qP
#                 per interval (averaged over substeps); everything else (resistant cells at the end of the
#                 look-ahead, burden, dose, toxicity) is a tie-breaker scaled by tiebreak_scale (user, 2026-11-01)
#   stopping      (user, 2026-11-01) every comparator runs the whole run (cheap); each MPC run stops computing
#                 once progression is confirmed at confirm_decisions successive decisions, or once it has
#                 outlasted every comparator by settle_margin_days (control_mpc.jl: stop_now); its time to
#                 progression is read up to the stop, so the run length is only a cap
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
# skip when nothing this run reads has changed since its last finished run (control_model.jl: control_fingerprint)
const FPRINT = control_fingerprint(@__FILE__, "CONTROL_DAYS" => DAYS,
    "CONTROL_LOOKAHEAD_TO_END" => get(ENV, "CONTROL_LOOKAHEAD_TO_END", ""))
const FPFILE = joinpath(COUT, "inputs.sha256")
if control_unchanged(FPFILE, FPRINT, joinpath(COUT, "control_outcomes.csv"))
    println("Control inputs unchanged since the last finished run; reusing $(joinpath(COUT, "control_outcomes.csv")) (CONTROL_FORCE=1 recomputes)")
    # provenance (user, 2026-11-02): tell run_pipeline.jl this step reused its last finished outputs (it copies this into run_manifest.toml)
    let step = get(ENV, "PIPELINE_STEP", "")
        isempty(step) || (mkpath(joinpath(@__DIR__, "outputs", "reuse")); write(joinpath(@__DIR__, "outputs", "reuse", step * ".toml"),
            "reused_from = \"$(Libc.strftime("%Y-%m-%dT%H:%M:%S", mtime(FPFILE)))\"\nreused_fingerprint = \"$FPRINT\"\n"))
    end
    exit(0)
end
cck = (isfile(CKP) && get(ENV, "REDO", "0") != "1") ? deserialize(CKP) : Dict{Symbol,Any}()
savec() = serialize(CKP, cck)

const Δ0 = Float64(CFG["schedule"]["decision_interval_days"])
const LOOK = Int.(CFG["schedule"]["lookahead_blocks"])            # look-ahead blocks (intervals per decision), the same at every decision
const LOOK_DAYS = round(Int, sum(LOOK) * Δ0)
const NRUN = round(Int, DAYS / Δ0)                                  # decisions over the whole run
# look-ahead ends with the run (shrinking horizon) or is fixed at every decision (PRD 2026-10-27)
const SHRINK = get(ENV, "CONTROL_LOOKAHEAD_TO_END", string(Int(get(CFG["schedule"], "lookahead_to_end", false)))) == "1"
const LEVELS = Float64.(CFG["dose"]["levels_data"])
const UMAX_DATA = Float64(CFG["dose"]["u_max_data"]); const UMAX_DEMO = Float64(CFG["dose"]["u_max_demo"])
const BOOT_N = min(Int(CFG["robustness"]["open_loop_draws"]), size(ZOO_STATE.results["s3"].res.pars[1], 1))
const NCL = min(Int(CFG["robustness"]["closed_loop_draws"]), BOOT_N)   # closed-loop draws (paired with every policy)
const CL_NAME = "Closed-loop MPC (observes N, C with noise)"

# starting counts: control_config [initial] source (control_model.jl: start_counts; by default the Stage 3a fitted
# start of each ratio's treated Ce1 series), the same in both environments
ratios_normal = [merge((ratio=r,), start_counts(r)) for r in ("25:75", "50:50", "75:25")]
ratios_lr = ratios_normal
plant_for(env, scen; b=nothing, tox=default_tox()) = env == :normal ? scenario_plant(scen; b=b, tox=tox) : lr_plant(scen; b=b, tox=tox)
ratios_for(env) = env == :normal ? ratios_normal : ratios_lr
println("=== optimal control: $(Threads.nthreads()) threads; Δ = $Δ0 d, $NRUN decisions ($DAYS days, $(SHRINK ? "look-ahead $LOOK_DAYS days, ending with the run" : "fixed look-ahead $LOOK_DAYS days"), blocks $LOOK); checkpoint has $(collect(keys(cck))) ==="); flush(stdout)
# label of an objective-weight vector [qB, qR, qT, ru, qP, qRT] (control_config.toml [weights] sweep)
wlabel(w) = length(w) == 6 && w[5] == 0 && w[6] == 0 ? "previous objective (burden + resistance): qB=$(w[1]) qR=$(w[2])" :
    "qP=$(w[5]) qRT=$(w[6]) qB=$(w[1]) qR=$(w[2]) qT=$(w[3]) ru=$(w[4])"
blabel(b) = "$(round(Int, sum(b) * Δ0)) days ($(length(b)) decisions$(all(==(1), b) ? "" : ", blocks $(join(b, "/"))"))"


function run_all_policies(P, x0; Δ=Δ0, n=NRUN, H=LOOK)
    W = default_weights(P); Wd = default_weights(P; umax=UMAX_DEMO)
    out = Dict{String,Any}()
    for (name, pol) in baseline_policies(P, W)
        U = run_policy(pol, P, x0, Δ, n); out[name] = (U=U, o=outcomes(P, x0, U, Δ, W))
    end
    sk = stop_kw(W, Δ, [v.o.ttp for v in values(out)])          # stop once progressed, or once past every comparator
    t = time(); cl = closed_loop(P, P, x0, Δ, n, H, W; mode=:discrete, levels=LEVELS, shrink=SHRINK, sk...)
    out["MPC, measured doses only"] = (U=cl.U, o=outcomes(P, x0, cl.U, Δ, W; t_stop=tstop(cl)), secs=time() - t, stop=cl.stop_reason)
    t = time(); cl = closed_loop(P, P, x0, Δ, n, H, Wd; mode=:continuous, shrink=SHRINK, sk...)
    out["MPC demonstration (continuous, ≤ $(UMAX_DEMO) µM)"] = (U=cl.U, o=outcomes(P, x0, cl.U, Δ, W; t_stop=tstop(cl)), secs=time() - t, iters=cl.iters, stop=cl.stop_reason)
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
# (sweep runs stop once progression is confirmed; they have no comparators to outlast)
function mpc_outcome(P, x0, W; Wo=W, Δ=Δ0, n=NRUN, H=LOOK, mode=:discrete)
    cl = closed_loop(P, P, x0, Δ, n, H, W; mode=mode, levels=LEVELS, shrink=SHRINK, stop_kw(Wo, Δ; settle=false)...)
    (U=cl.U, o=outcomes(P, x0, cl.U, Δ, Wo; t_stop=tstop(cl)))
end
if !haskey(cck, :sweeps)
    r = ratios_normal[2]
    tasks = Any[]
    for scen in (:A, :B)
        P = plant_for(:normal, scen); x0 = x0_for(P, r.N0, r.C0); W0 = default_weights(P)
        for w in CFG["weights"]["sweep"]
            push!(tasks, (kind="weights", scen=scen, label=wlabel(w), f=() -> mpc_outcome(P, x0, default_weights(P; w=Float64.(w)); Wo=W0)))
        end
        # the tie-break scale of the second tier (user, 2026-11-12): the 14-day schedule is ranked by that tier
        for ts in CFG["progression"]["tiebreak_sweep"]
            push!(tasks, (kind="tie-break scale", scen=scen, label="tie-break scale = $ts", f=() -> mpc_outcome(P, x0, default_weights(P; ts=Float64(ts)); Wo=W0)))
        end
        # soft burden penalty: where it starts (0 = none) × the terminal weight; scored against the default B_cancer
        for bc in CFG["cancer"]["B_cancer_sweep"], tf in CFG["cancer"]["terminal_factor_sweep"]
            lab = (bc > 0 ? "B_cancer = $(bc) K" : "no burden penalty") * ", terminal ×$(tf)"
            # the relaxed log barrier explicitly: the default barrier is "none" since 2026-11-01
            push!(tasks, (kind="burden penalty", scen=scen, label=lab, f=() -> mpc_outcome(P, x0, default_weights(P; bc_frac=Float64(bc), tf=Float64(tf), barrier=:relaxed_log); Wo=W0)))
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
        for b in CFG["schedule"]["lookahead_sweep"]
            bb = Int.(b)
            push!(tasks, (kind="horizon", scen=scen, label=blabel(bb), f=() -> mpc_outcome(P, x0, default_weights(P); H=bb)))
        end
        # decide every day with the same 28-day look-ahead (blocks in 1-day intervals)
        push!(tasks, (kind="interval", scen=scen, label="decide every 1 day ($(LOOK_DAYS)-day look-ahead)", f=() -> mpc_outcome(P, x0, default_weights(P); Δ=1.0, n=2NRUN, H=[1, 1, 2, 4, 4, 8, 8])))
        push!(tasks, (kind="interval", scen=scen, label="decide every 2 days ($(LOOK_DAYS)-day look-ahead)", f=() -> mpc_outcome(P, x0, default_weights(P))))
        # the 28-day extrapolation row belongs to the 14-day report; longer reports are extrapolations throughout
        DAYS == Float64(CFG["schedule"]["closed_loop_days"]) && push!(tasks, (kind="extrapolation", scen=scen, label="28 days, receding $(LOOK_DAYS)-day look-ahead (beyond the data)", f=() -> mpc_outcome(P, x0, default_weights(P); n=14)))
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
        summ(o) = (N=o.N, C=o.C, T=o.T, cum=o.cum_dose, cv=o.cancer_violation, tv=o.tox_violation, dover=o.days_over, Nf=o.N_final, Cf=o.C_final, tot=o.total_final, rf=o.resistant_fraction, maxT=o.max_tox,
                   ttp=o.ttp, ttp_first=o.ttp_first, dprog=o.days_progressed)
        Threads.@threads :dynamic for b in 1:nb
            Pb = try plant_for(env, scen; b=b) catch; nothing end
            if Pb === nothing
                draws[b] = nothing; continue
            end
            Wb = with_fields(default_weights(Pb); Bc=W0.Bc, Bp=W0.Bp, BpN=W0.BpN, BpC=W0.BpC)   # the same thresholds as the nominal plan
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
                # stops once progressed on plant b, or once past every comparator on the same draw
                sk = stop_kw(Wb, Δ0, [d[pn] === nothing ? NaN : d[pn].ttp for pn in pols if pn != "MPC, measured doses only"])
                c = try closed_loop(P0, Pb, m.x0, Δ0, NRUN, LOOK, W0; mode=:discrete, levels=LEVELS, noise=noise, shrink=SHRINK, sk...) catch; nothing end
                d[CL_NAME] = c === nothing ? nothing : merge(summ(outcomes(Pb, m.x0, c.U, Δ0, Wb; t_stop=tstop(c))), (U=c.U, stop=c.stop_reason))
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
    o = v.o; sr = get(v, :stop, "")
    push!(rows, (environment=string(m.env), scenario=string(m.scen), ratio=m.ratio, policy=pn, schedule_uM=join(round.(v.U; digits=3), " "),
                 N_final=o.N_final, C_final=o.C_final, total_final=o.total_final, resistant_fraction=o.resistant_fraction,
                 max_total=o.max_total, max_toxicity=o.max_tox, cumulative_dose_uM_days=o.cum_dose,
                 cancer_bound_violated=o.cancer_violation, days_above_B_cancer=o.days_over, days_above_half_K=get(o, :days_half_K, NaN), toxicity_bound_violated=o.tox_violation,
                 first_progression_day=o.ttp_first, time_to_progression_day=o.ttp, days_progressed=o.days_progressed,
                 N_at_progression=o.N_prog, C_at_progression=o.C_prog, control_stopped_day=o.stop_day, stop_reason=sr))
end
CSV.write(joinpath(COUT, "control_outcomes.csv"), DataFrame(rows))
println("Wrote $(joinpath(COUT, "control_outcomes.csv"))")
write(FPFILE, FPRINT)
