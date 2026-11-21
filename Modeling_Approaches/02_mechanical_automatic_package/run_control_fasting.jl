# Long-run optimal control (long_days, default 120 days; CONTROL_DAYS overrides) with two inputs: cisplatin dose u and low resource / fasting v (control_fasting.jl).
#   julia -t 12 --project=. run_control_fasting.jl        (REDO=1 recomputes; the checkpoint is per section)
# Writes outputs/control_<days>d/fasting_state.jls and fasting_outcomes.csv; control_fasting_report.jl renders the page.
#
#   scenarios   A (selected coculture drug response), B (monoculture drug response for both lines)
#   ratios      25:75, 50:50, 75:25 from the NOMINAL seeding (control_config.toml [initial]: 100 cells split
#               by the ratio), not a measured day-0 count
#   policies    8 fixed schedules (5 fed, 3 with fasting), MPC dose only (fed), MPC dose + fasting (measured
#               levels, exact search), MPC demonstration (continuous u ≤ u_max_demo, v ∈ [0, 1], SQP);
#               every MPC looks a FIXED look-ahead ahead (control_config lookahead_blocks: 28 days in 6 blocked
#               decisions) at every decision, also near the end of the run; objective: time to progression
#   checks      the plant against both measured low-resource conditions (RMSE in cells)
#   sweeps      protection p_H, weight-loss rate k_W, weight limit W_min, fasting cost, interaction mode,
#               look-ahead, weights, soft burden penalty (B_cancer) and terminal factor
#   robustness  fixed schedules on the open-loop bootstrap draws; closed-loop MPC (observes N and C with
#               measurement noise, reconstructs the hidden states with its own model) on the first
#               closed-loop draws, where every policy is compared on the same draws and noise

include(joinpath(@__DIR__, "control_fasting.jl"))

const DAYS = parse(Float64, get(ENV, "CONTROL_DAYS", string(CFG["schedule"]["long_days"])))
const Δ0 = Float64(CFG["schedule"]["decision_interval_days"])
const NRUN = round(Int, DAYS / Δ0)                                   # one decision per 2-day interval
const LOOK = Int.(CFG["schedule"]["lookahead_blocks"])              # fixed look-ahead blocks (move blocking)
const LOOK_DAYS = round(Int, sum(LOOK) * Δ0)
const UMAX_DEMO = Float64(CFG["dose"]["u_max_demo"])
const FOUT = joinpath(get(ENV, "CONTROL_OUT", joinpath(@__DIR__, "outputs")), control_subdir(DAYS, NaN)); mkpath(FOUT)   # CONTROL_OUT: side runs, as in run_control.jl
const FCK = joinpath(FOUT, "fasting_state.jls")
# skip when nothing this run reads has changed since its last finished run (control_model.jl: control_fingerprint)
const FPRINT = control_fingerprint(@__FILE__, "CONTROL_DAYS" => DAYS)
const FPFILE = joinpath(FOUT, "inputs.sha256")
if control_unchanged(FPFILE, FPRINT, joinpath(FOUT, "fasting_outcomes.csv"))
    println("Control inputs unchanged since the last finished run; reusing $(joinpath(FOUT, "fasting_outcomes.csv")) (CONTROL_FORCE=1 recomputes)")
    # provenance (user, 2026-11-02): tell run_pipeline.jl this step reused its last finished outputs (it copies this into run_manifest.toml)
    let step = get(ENV, "PIPELINE_STEP", "")
        isempty(step) || (mkpath(joinpath(@__DIR__, "outputs", "reuse")); write(joinpath(@__DIR__, "outputs", "reuse", step * ".toml"),
            "reused_from = \"$(Libc.strftime("%Y-%m-%dT%H:%M:%S", mtime(FPFILE)))\"\nreused_fingerprint = \"$FPRINT\"\n"))
    end
    exit(0)
end
fck = (isfile(FCK) && get(ENV, "REDO", "0") != "1") ? deserialize(FCK) : Dict{Symbol,Any}()
savef() = serialize(FCK, fck)
const BOOT_N = min(Int(CFG["robustness"]["open_loop_draws"]), size(ZOO_STATE.results["s3"].res.pars[1], 1), size(ZOO_STATE.results["s4"].res.pars[1], 1))
const NCL = min(Int(CFG["robustness"]["closed_loop_draws"]), BOOT_N)   # closed-loop robustness runs per scenario (paired with every policy)
const CL_NAME = "Closed-loop MPC, dose + fasting (observes N, C with noise)"
# starting counts: control_config [initial] source (control_model.jl: start_counts; the fitted Ce1 starts by default)
ratios = [merge((ratio=r,), start_counts(r)) for r in ("25:75", "50:50", "75:25")]
const MPC2 = "MPC, dose + fasting (measured levels)"
const MPC1 = "MPC, dose only (fed)"
const MPCD = "MPC demonstration (continuous u ≤ $(UMAX_DEMO) µM, v ∈ [0, 1])"
println("=== $(round(Int, DAYS))-day control with fasting: $(Threads.nthreads()) threads; $NRUN decisions, fixed look-ahead $LOOK_DAYS days (blocks $LOOK); checkpoint has $(collect(keys(fck))) ==="); flush(stdout)
# label of an objective-weight vector [qB, qR, qT, ru, qP, qRT] (control_config.toml [weights] sweep)
wlabel(w) = length(w) == 6 && w[5] == 0 && w[6] == 0 ? "previous objective (burden + resistance): qB=$(w[1]) qR=$(w[2])" :
    "qP=$(w[5]) qRT=$(w[6]) qB=$(w[1]) qR=$(w[2]) qT=$(w[3]) ru=$(w[4])"
blabel(b) = "$(round(Int, sum(b) * Δ0)) days ($(length(b)) decisions$(all(==(1), b) ? "" : ", blocks $(join(b, "/"))"))"


# ── checks against the two measured low-resource conditions ──────────────────
# Every RMSE is in cells (unweighted). The controller's plant (fed the nominal seeding, constant low resource,
# no drug) is scored on every observation after day 0 (QC-excluded points left out) of the six low-resource
# cocultures; the plant and the Stage 3b fit (its own model and starting state) are also compared on exactly
# the points both use (the Stage 3b data can differ slightly in QC), so the two numbers are like for like.
function compute_checks()
    out = Dict{Symbol,Any}()
    g3b = zgroup("s3b"); c3b = g3b.cands[1]
    sol3b = solve(ODEProblem(c3b.model, c3b.u0, (0.0, maximum(maximum(d.x) for d in g3b.datasets)), c3b.params), Tsit5(); reltol=1e-10, abstol=1e-10)
    # the Stage 3b dataset holding most of these observations (matched on the (day, count) pairs) and the shared points
    function match3b(xk, yk)
        best = (j=0, idx=Int[])
        for (j, d) in enumerate(g3b.datasets)
            idx = [i for i in eachindex(xk) if any(m -> d.x[m] == xk[i] && isapprox(d.y[m], yk[i]; rtol=1e-9), eachindex(d.x))]
            length(idx) > length(best.idx) && (best = (j=j, idx=idx))
        end
        length(best.idx) >= 3 ? best : nothing
    end
    ratio_of = Dict("25S75R" => "25:75", "50S50R" => "50:50", "75S25R" => "75:25")
    for mode in (:interp, :stage4)
        FP = fast_plant(:A; eps_mode=mode); se = 0.0; n = 0; sec = 0.0; se3 = 0.0; n3 = 0; per = NamedTuple[]
        for c in lr_cc
            s = nominal_seeding(ratio_of[c.label])
            ts, xs = trajectory2(x0_fast(FP, s.N0, s.C0), repeat([0.0 1.0], 7), FP, 2.0; dt=0.5)
            for (lab, xd, yd, idx, bad) in (("A2780", c.xN, c.yN, FP.iN, c.bN), ("A2780cis", c.xC, c.yC, FP.iC, c.bC))
                keep = [i for i in eachindex(xd) if xd[i] > 0 && !bad[i]]
                xk, yk = xd[keep], yd[keep]
                pm = [xs[argmin(abs.(ts .- t))][idx] for t in xk]
                s2 = sum(abs2, pm .- yk); se += s2; n += length(keep)
                mt = match3b(xk, yk); r3 = NaN; rpc = NaN
                if mt !== nothing
                    I = mt.idx; p3 = [sol3b(xk[i])[g3b.datasets[mt.j].state_index] for i in I]
                    e3 = sum(abs2, p3 .- yk[I]); ep = sum(abs2, pm[I] .- yk[I])
                    se3 += e3; sec += ep; n3 += length(I); r3 = sqrt(e3 / length(I)); rpc = sqrt(ep / length(I))
                end
                push!(per, (run=c.run, ratio=c.label, line=lab, rmse=sqrt(s2 / length(keep)), rmse_3b=r3, rmse_common=rpc,
                            last_obs=yd[end], last_pred=xs[argmin(abs.(ts .- xd[end]))][idx]))
            end
        end
        out[mode] = (rmse=sqrt(se / n), n=n, rmse_common=n3 > 0 ? sqrt(sec / n3) : NaN, per=per)
        out[:s3b_rmse] = n3 > 0 ? sqrt(se3 / n3) : NaN; out[:s3b_n] = n3
    end
    FP = fast_plant(:A); c4 = zgroup("s4").cands[1]
    sol = solve(ODEProblem(c4.model, c4.u0, (0.0, 12.0), c4.params), Tsit5(); reltol=1e-10, abstol=1e-10)
    x = vcat(c4.u0, zeros(3), 1.0); mr = 0.0
    for k in 1:6
        x = F2(x, [1.0, 1.0], FP, 2.0); s = sol(2.0k)
        mr = max(mr, abs(x[FP.iN] - s[FP.iN]) / max(s[FP.iN], 50), abs(x[FP.iC] - s[FP.iC]) / max(s[FP.iC], 50))
    end
    out[:stage4_maxrel] = mr
    out
end
if !haskey(fck, :checks)
    out = compute_checks()
    fck[:checks] = out; savef()
    @printf("  checks: low resource alone RMSE (cells) interp %.0f, stage4 %.0f on %d points; on the %d points shared with Stage 3b: plant %.0f, Stage 3b fit %.0f; Stage 4 reproduced to %.1e\n",
            out[:interp].rmse, out[:stage4].rmse, out[:interp].n, out[:s3b_n], out[:interp].rmse_common, out[:s3b_rmse], out[:stage4_maxrel]); flush(stdout)
end

function run_all2(FP, x0; n=NRUN, H=LOOK)
    W = weights2(FP); Wd = weights2(FP; umax=UMAX_DEMO)
    out = Dict{String,Any}()
    for (name, pol) in baseline_policies2(FP, W)
        U = run_policy2(pol, FP, x0, Δ0, n); out[name] = (U=U, o=outcomes2(FP, x0, U, Δ0, W))
    end
    # each MPC stops once its progression is confirmed or once it has outlasted every comparator (control_mpc.jl)
    sk = stop_kw(W, Δ0, [v.o.ttp for v in values(out)])
    t = time(); cl = closed_loop2(FP, FP, x0, Δ0, n, H, W; vlevels=[0.0], sk...)
    out[MPC1] = (U=cl.U, o=outcomes2(FP, x0, cl.U, Δ0, W; t_stop=tstop(cl)), secs=time() - t, stop=cl.stop_reason)
    t = time(); cl = closed_loop2(FP, FP, x0, Δ0, n, H, W; sk...)
    out[MPC2] = (U=cl.U, o=outcomes2(FP, x0, cl.U, Δ0, W; t_stop=tstop(cl)), secs=time() - t, stop=cl.stop_reason)
    t = time(); cl = closed_loop2(FP, FP, x0, Δ0, n, H, Wd; mode=:continuous, sk...)
    out[MPCD] = (U=cl.U, o=outcomes2(FP, x0, cl.U, Δ0, W; t_stop=tstop(cl)), secs=time() - t, stop=cl.stop_reason)
    out
end

if !haskey(fck, :main)
    tasks = [(scen, r) for scen in (:A, :B) for r in ratios]
    out = Vector{Any}(undef, length(tasks))
    Threads.@threads :dynamic for i in eachindex(tasks)
        scen, r = tasks[i]; t = time()
        FP = fast_plant(scen); x0 = x0_fast(FP, r.N0, r.C0)
        out[i] = (scen=scen, ratio=r.ratio, x0=x0, name=FP.name, K=FP.K, W=weights2(FP), policies=run_all2(FP, x0))
        @printf("  main %s %s done [%.0fs]\n", scen, r.ratio, time() - t); flush(stdout)
    end
    fck[:main] = out; savef()
end
main_of(scen, ratio) = only(filter(m -> m.scen == scen && m.ratio == ratio, fck[:main]))

# ── sweeps (50:50, MPC dose + fasting) ─────────────────────────────────────────
if !haskey(fck, :sweeps)
    r = ratios[2]; tasks = Any[]
    # Wo: the weights the outcomes are scored with (bound counts always at the default B_cancer)
    # sweep runs stop once progression is confirmed (no comparators to outlast)
    mpc(FP, W; Wo=W, H=LOOK, kw...) = (x0 = x0_fast(FP, r.N0, r.C0); cl = closed_loop2(FP, FP, x0, Δ0, NRUN, H, W; stop_kw(Wo, Δ0; settle=false)..., kw...);
                                      (U=cl.U, o=outcomes2(FP, x0, cl.U, Δ0, Wo; t_stop=tstop(cl))))
    h0 = default_host()
    bc0 = Float64(CFG["cancer"]["B_cancer_fraction_of_K"]); tf0 = Float64(CFG["weights"]["terminal_factor"])
    # soft burden penalty, one setting at a time (the full B_cancer × terminal grid is in the 14-day report)
    pen = vcat([(Float64(b), tf0) for b in CFG["cancer"]["B_cancer_sweep"] if b != bc0], [(bc0, Float64(t)) for t in CFG["cancer"]["terminal_factor_sweep"] if t != tf0])
    for scen in (:A, :B)
        for v in HCFG["p_H_sweep"]
            push!(tasks, (kind="protection", scen=scen, label="p_H = $v", f=() -> (FP = fast_plant(scen; host=merge(h0, (p_H=Float64(v),))); mpc(FP, weights2(FP)))))
        end
        for v in HCFG["k_W_sweep"]
            push!(tasks, (kind="weight loss rate", scen=scen, label="k_W = $v /day", f=() -> (FP = fast_plant(scen; host=merge(h0, (k_W=Float64(v),))); mpc(FP, weights2(FP)))))
        end
        for v in HCFG["W_min_sweep"]
            push!(tasks, (kind="weight limit", scen=scen, label="W_min = $v", f=() -> (FP = fast_plant(scen); mpc(FP, weights2(FP; Wmin=Float64(v))))))
        end
        for v in (0.01, 0.05, 0.2)
            push!(tasks, (kind="fasting cost", scen=scen, label="r_v = $v", f=() -> (FP = fast_plant(scen); mpc(FP, weights2(FP; rv=v)))))
        end
        push!(tasks, (kind="drug × low resource", scen=scen, label="ε = Stage 4 value at every dose", f=() -> (FP = fast_plant(scen; eps_mode=:stage4); mpc(FP, weights2(FP)))))
        push!(tasks, (kind="fasting rules", scen=scen, label="fasts up to 4 days (2 intervals)", f=() -> (FP = fast_plant(scen); mpc(FP, weights2(FP); maxrun=2))))
        for b in CFG["schedule"]["lookahead_sweep"]
            bb = Int.(b); bb == LOOK && continue
            push!(tasks, (kind="look-ahead", scen=scen, label=blabel(bb), f=() -> (FP = fast_plant(scen); mpc(FP, weights2(FP); H=bb))))
        end
        for w in CFG["weights"]["sweep"][2:end]          # the first entry is the default
            push!(tasks, (kind="weights", scen=scen, label=wlabel(w), f=() -> (FP = fast_plant(scen); mpc(FP, weights2(FP; w=Float64.(w)); Wo=weights2(FP)))))
        end
        for (bc, tf) in pen
            lab = (bc > 0 ? "B_cancer = $(bc) K" : "no burden penalty") * ", terminal ×$(tf)"
            push!(tasks, (kind="burden penalty", scen=scen, label=lab, f=() -> (FP = fast_plant(scen); mpc(FP, weights2(FP; bc_frac=bc, tf=tf); Wo=weights2(FP)))))   # (the barrier follows [cancer] barrier)
        end
    end
    out = Vector{Any}(undef, length(tasks))
    Threads.@threads :dynamic for i in eachindex(tasks)
        T = tasks[i]; t = time(); out[i] = (kind=T.kind, scen=T.scen, label=T.label, res=T.f())
        @printf("  sweep %s %s %s [%.0fs]\n", T.scen, T.kind, T.label, time() - t); flush(stdout)
    end
    fck[:sweeps] = out; savef()
end

# ── robustness over bootstrap parameter sets ───────────────────────────────────
# Every draw b: the fixed schedules (planned once with the best-fit model) are replayed on plant b; on the
# first NCL draws the closed-loop MPC plans with the best-fit model, observes plant b's N and C with
# measurement noise and reconstructs the hidden states with its own model. Paired: every policy meets the
# same plants on those draws.
if !haskey(fck, :robust)
    rob = Dict{Any,Any}()
    pols = ["No treatment", "Constant 1 µM", "Maximum measured dose (1.47 µM)", "Periodic 1 µM, 2 days on / 2 off",
            "Periodic 1 µM with a 48-h fast on each dosing interval", "Constant 1 µM + 48-h fast every 6 days", "Constant 1 µM + 48-h fast every 8 days", MPC1, MPC2]
    summ(o) = (N=o.N, C=o.C, T=o.T, W=o.W, cum=o.cum_dose, fast=o.fast_days, cv=o.cancer_violation, tv=o.tox_violation, wv=o.weight_violation, dover=o.days_over,
               Nf=o.N_final, Cf=o.C_final, tot=o.total_final, rf=o.resistant_fraction, minW=o.min_weight, maxT=o.max_tox,
               ttp=o.ttp, ttp_first=o.ttp_first, dprog=o.days_progressed)
    for scen in (:A, :B)
        m = main_of(scen, "50:50"); FP0 = fast_plant(scen); W0 = weights2(FP0)
        draws = Vector{Any}(undef, BOOT_N)
        Threads.@threads :dynamic for b in 1:BOOT_N
            FPb = try fast_plant(scen; b=b) catch; nothing end
            if FPb === nothing
                draws[b] = nothing; continue
            end
            d = Dict{String,Any}("b" => b)
            for pn in pols
                o = try outcomes2(FPb, m.x0, m.policies[pn].U, Δ0, W0) catch; nothing end
                d[pn] = o === nothing ? nothing : summ(o)
            end
            if b <= NCL
                sk = stop_kw(W0, Δ0, [d[pn] === nothing ? NaN : d[pn].ttp for pn in pols if !(pn in (MPC1, MPC2))])
                c = try closed_loop2(FP0, FPb, m.x0, Δ0, NRUN, LOOK, W0; noise=obs_noise(b, NRUN), sk...) catch; nothing end
                d[CL_NAME] = c === nothing ? nothing : merge(summ(outcomes2(FPb, m.x0, c.U, Δ0, W0; t_stop=tstop(c))), (U=c.U, stop=c.stop_reason))
            end
            draws[b] = d
        end
        ds = [d for d in draws if d !== nothing]
        rob[scen] = (ts=m.policies["No treatment"].o.ts, draws=ds, closed=[d[CL_NAME] for d in ds if get(d, CL_NAME, nothing) !== nothing],
                     ncl=NCL, cv=Float64(CFG["observation"]["count_cv"]))
        @printf("  robustness %s: %d draws, %d closed-loop\n", scen, length(rob[scen].draws), length(rob[scen].closed)); flush(stdout)
    end
    fck[:robust] = rob; savef()
end

# ── recovery of starvation damage after refeeding (not measured; the fitted model lets it fade at k_C^R) ──
if !haskey(fck, :washout)
    tasks = Any[]
    for scen in (:A, :B), ko in (0.5, 2.0)
        push!(tasks, (scen=scen, koff=ko))
    end
    out = Vector{Any}(undef, length(tasks))
    Threads.@threads :dynamic for i in eachindex(tasks)
        T = tasks[i]; t = time(); m = main_of(T.scen, "50:50")
        FP = fast_plant(T.scen; koff=T.koff); W = weights2(FP)
        replay = Dict(pn => outcomes2(FP, m.x0, m.policies[pn].U, Δ0, W) for pn in ("Fasting only, 48 h every 6 days", "Constant 1 µM", MPC1, MPC2))
        cl = closed_loop2(FP, FP, m.x0, Δ0, NRUN, LOOK, W; stop_kw(W, Δ0; settle=false)...)
        out[i] = (scen=T.scen, koff=T.koff, replay=replay, U=cl.U, o=outcomes2(FP, m.x0, cl.U, Δ0, W; t_stop=tstop(cl)))
        @printf("  washout %s k_off=%.1f [%.0fs]\n", T.scen, T.koff, time() - t); flush(stdout)
    end
    fck[:washout] = out; savef()
end

# ── CSV ────────────────────────────────────────────────────────────────────────
rows = NamedTuple[]
for m in fck[:main], (pn, v) in m.policies
    o = v.o; sr = get(v, :stop, "")
    push!(rows, (scenario=string(m.scen), ratio=m.ratio, policy=pn, dose_uM=join(round.(v.U[:, 1]; digits=3), " "), fasting=join(round.(v.U[:, 2]; digits=3), " "),
                 N_final=o.N_final, C_final=o.C_final, total_final=o.total_final, resistant_fraction=o.resistant_fraction, max_total=o.max_total,
                 max_toxicity=o.max_tox, min_weight=o.min_weight, cumulative_dose_uM_days=o.cum_dose, fasting_days=o.fast_days,
                 cancer_bound_violated=o.cancer_violation, days_above_B_cancer=o.days_over, days_above_half_K=get(o, :days_half_K, NaN), toxicity_bound_violated=o.tox_violation, weight_bound_violated=o.weight_violation,
                 first_progression_day=o.ttp_first, time_to_progression_day=o.ttp, days_progressed=o.days_progressed,
                 N_at_progression=o.N_prog, C_at_progression=o.C_prog, control_stopped_day=o.stop_day, stop_reason=sr))
end
CSV.write(joinpath(FOUT, "fasting_outcomes.csv"), DataFrame(rows))
println("Wrote $(joinpath(FOUT, "fasting_outcomes.csv"))")
write(FPFILE, FPRINT)
