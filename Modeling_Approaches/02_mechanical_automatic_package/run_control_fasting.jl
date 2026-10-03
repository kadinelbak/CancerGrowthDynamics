# 30-day optimal control with two inputs: cisplatin dose u and low resource / fasting v (control_fasting.jl).
#   julia -t 12 --project=. run_control_fasting.jl        (REDO=1 recomputes; the checkpoint is per section)
# Writes outputs/control_30d/fasting_state.jls and fasting_outcomes.csv; control_fasting_report.jl renders the page.
#
#   scenarios   A (selected coculture drug response), B (monoculture drug response for both lines)
#   ratios      25:75, 50:50, 75:25 from the untreated-coculture day-0 counts
#   policies    8 fixed schedules (5 fed, 3 with fasting), MPC dose only (fed), MPC dose + fasting (measured
#               levels, exact search), MPC demonstration (continuous u ≤ u_max_demo, v ∈ [0, 1], SQP)
#   checks      the plant against both measured low-resource conditions
#   sweeps      protection p_H, weight-loss rate k_W, weight limit W_min, fasting cost, interaction mode,
#               look-ahead, weights
#   robustness  fixed schedules on 200 bootstrap parameter sets; closed-loop MPC on a subset

include(joinpath(@__DIR__, "control_fasting.jl"))

const DAYS = Float64(CFG["schedule"]["long_days"])
const Δ0 = Float64(CFG["schedule"]["decision_interval_days"])
const NRUN = round(Int, DAYS / Δ0)                                   # 15 decisions over 30 days
const NDEC = round(Int, CFG["schedule"]["horizon_days"] / Δ0)       # 7-interval (14-day) look-ahead
const UMAX_DEMO = Float64(CFG["dose"]["u_max_demo"])
const FOUT = joinpath(get(ENV, "CONTROL_OUT", joinpath(@__DIR__, "outputs")), "control_30d"); mkpath(FOUT)   # CONTROL_OUT: side runs, as in run_control.jl
const FCK = joinpath(FOUT, "fasting_state.jls")
fck = (isfile(FCK) && get(ENV, "REDO", "0") != "1") ? deserialize(FCK) : Dict{Symbol,Any}()
savef() = serialize(FCK, fck)
const BOOT_N = min(200, size(ZOO_STATE.results["s3"].res.pars[1], 1), size(ZOO_STATE.results["s4"].res.pars[1], 1))
const NCL = 24                                                       # closed-loop robustness runs per scenario
ratios = [(ratio=("25:75", "50:50", "75:25")[i], N0=co_data[i].S[2][1], C0=co_data[i].R[2][1]) for i in 1:3]
const MPC2 = "MPC, dose + fasting (measured levels)"
const MPC1 = "MPC, dose only (fed)"
const MPCD = "MPC demonstration (continuous u ≤ $(UMAX_DEMO) µM, v ∈ [0, 1])"
println("=== 30-day control with fasting: $(Threads.nthreads()) threads; $NRUN decisions, look-ahead $NDEC; checkpoint has $(collect(keys(fck))) ==="); flush(stdout)

# ── checks against the two measured low-resource conditions ──────────────────
function compute_checks()
    out = Dict{Symbol,Any}()
    for mode in (:interp, :stage4)
        FP = fast_plant(:A; eps_mode=mode); se = 0.0; n = 0; per = NamedTuple[]
        for c in lr_cc
            ts, xs = trajectory2(x0_fast(FP, c.yN[1], c.yC[1]), repeat([0.0 1.0], 7), FP, 2.0; dt=0.5)
            for (lab, xd, yd, idx, bad) in (("A2780", c.xN, c.yN, FP.iN, c.bN), ("A2780cis", c.xC, c.yC, FP.iC, c.bC))
                s2 = 0.0; m = 0
                for (i, t) in enumerate(xd)
                    (i == 1 || bad[i]) && continue
                    j = argmin(abs.(ts .- t)); s2 += (xs[j][idx] - yd[i])^2; m += 1
                end
                se += s2; n += m; push!(per, (run=c.run, ratio=c.label, line=lab, rmse=sqrt(s2 / m), last_obs=yd[end], last_pred=xs[argmin(abs.(ts .- xd[end]))][idx]))
            end
        end
        out[mode] = (rmse=sqrt(se / n), n=n, per=per)
    end
    g3b = zgroup("s3b"); c3b = g3b.cands[1]
    out[:s3b_rmse] = sqrt(c3b.sse / sum(length(d.x) for d in g3b.datasets))
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
    @printf("  checks: low resource alone RMSE interp %.0f, stage4 %.0f (Stage 3b fit %.0f); Stage 4 reproduced to %.1e\n", out[:interp].rmse, out[:stage4].rmse, out[:s3b_rmse], out[:stage4_maxrel]); flush(stdout)
end

function run_all2(FP, x0; n=NRUN, H=NDEC)
    W = weights2(FP); Wd = weights2(FP; umax=UMAX_DEMO)
    out = Dict{String,Any}()
    for (name, pol) in baseline_policies2(FP, W)
        U = run_policy2(pol, FP, x0, Δ0, n); out[name] = (U=U, o=outcomes2(FP, x0, U, Δ0, W))
    end
    t = time(); cl = closed_loop2(FP, FP, x0, Δ0, n, H, W; vlevels=[0.0])
    out[MPC1] = (U=cl.U, o=outcomes2(FP, x0, cl.U, Δ0, W), secs=time() - t)
    t = time(); cl = closed_loop2(FP, FP, x0, Δ0, n, H, W)
    out[MPC2] = (U=cl.U, o=outcomes2(FP, x0, cl.U, Δ0, W), secs=time() - t)
    t = time(); cl = closed_loop2(FP, FP, x0, Δ0, n, H, Wd; mode=:continuous)
    out[MPCD] = (U=cl.U, o=outcomes2(FP, x0, cl.U, Δ0, W), secs=time() - t)
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
    mpc(FP, W; H=NDEC, kw...) = (x0 = x0_fast(FP, r.N0, r.C0); cl = closed_loop2(FP, FP, x0, Δ0, NRUN, H, W; kw...); (U=cl.U, o=outcomes2(FP, x0, cl.U, Δ0, W)))
    h0 = default_host()
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
        for Hh in (3, 5)
            push!(tasks, (kind="look-ahead", scen=scen, label="$(Hh) decisions ($(round(Int, Hh * Δ0)) days)", f=() -> (FP = fast_plant(scen); mpc(FP, weights2(FP); H=Hh))))
        end
        for w in ([1.0, 3.0, 1.0, 0.1], [3.0, 1.0, 1.0, 0.1], [1.0, 1.0, 3.0, 0.1])
            push!(tasks, (kind="weights", scen=scen, label="qB=$(w[1]) qR=$(w[2]) qT=$(w[3]) ru=$(w[4])", f=() -> (FP = fast_plant(scen); mpc(FP, weights2(FP; w=w)))))
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
if !haskey(fck, :robust)
    rob = Dict{Any,Any}()
    pols = ["No treatment", "Constant 1 µM", "Periodic 1 µM, 2 days on / 2 off", "Periodic 1 µM with a 48-h fast on each dosing interval", MPC1, MPC2]
    for scen in (:A, :B)
        m = main_of(scen, "50:50"); FP0 = fast_plant(scen); W0 = weights2(FP0)
        draws = Vector{Any}(undef, BOOT_N)
        Threads.@threads :dynamic for b in 1:BOOT_N
            FPb = try fast_plant(scen; b=b) catch; nothing end
            if FPb === nothing
                draws[b] = nothing; continue
            end
            d = Dict{String,Any}()
            for pn in pols
                o = try outcomes2(FPb, m.x0, m.policies[pn].U, Δ0, W0) catch; nothing end
                d[pn] = o === nothing ? nothing : (N=o.N, C=o.C, T=o.T, W=o.W, cum=o.cum_dose, fast=o.fast_days, cv=o.cancer_violation, tv=o.tox_violation,
                                                   wv=o.weight_violation, Nf=o.N_final, Cf=o.C_final, tot=o.total_final, rf=o.resistant_fraction, minW=o.min_weight)
            end
            draws[b] = d
        end
        cl = Vector{Any}(undef, NCL)
        Threads.@threads :dynamic for b in 1:NCL
            FPb = try fast_plant(scen; b=b) catch; nothing end
            cl[b] = FPb === nothing ? nothing : (c = closed_loop2(FP0, FPb, m.x0, Δ0, NRUN, NDEC, W0);
                                               o = outcomes2(FPb, m.x0, c.U, Δ0, W0);
                                               (U=c.U, N=o.N, C=o.C, T=o.T, W=o.W, cum=o.cum_dose, fast=o.fast_days, cv=o.cancer_violation, tv=o.tox_violation,
                                                wv=o.weight_violation, tot=o.total_final, rf=o.resistant_fraction, minW=o.min_weight))
        end
        rob[scen] = (ts=m.policies["No treatment"].o.ts, draws=[d for d in draws if d !== nothing], closed=[c for c in cl if c !== nothing])
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
        cl = closed_loop2(FP, FP, m.x0, Δ0, NRUN, NDEC, W)
        out[i] = (scen=T.scen, koff=T.koff, replay=replay, U=cl.U, o=outcomes2(FP, m.x0, cl.U, Δ0, W))
        @printf("  washout %s k_off=%.1f [%.0fs]
", T.scen, T.koff, time() - t); flush(stdout)
    end
    fck[:washout] = out; savef()
end

# ── CSV ────────────────────────────────────────────────────────────────────────
rows = NamedTuple[]
for m in fck[:main], (pn, v) in m.policies
    o = v.o
    push!(rows, (scenario=string(m.scen), ratio=m.ratio, policy=pn, dose_uM=join(round.(v.U[:, 1]; digits=3), " "), fasting=join(round.(v.U[:, 2]; digits=3), " "),
                 N_final=o.N_final, C_final=o.C_final, total_final=o.total_final, resistant_fraction=o.resistant_fraction, max_total=o.max_total,
                 max_toxicity=o.max_tox, min_weight=o.min_weight, cumulative_dose_uM_days=o.cum_dose, fasting_days=o.fast_days,
                 cancer_bound_violated=o.cancer_violation, toxicity_bound_violated=o.tox_violation, weight_bound_violated=o.weight_violation))
end
CSV.write(joinpath(FOUT, "fasting_outcomes.csv"), DataFrame(rows))
println("Wrote $(joinpath(FOUT, "fasting_outcomes.csv"))")
