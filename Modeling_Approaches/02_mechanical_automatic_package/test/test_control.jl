# Tests for the optimal-control layer:  julia --project=. test/test_control.jl
include(joinpath(@__DIR__, "..", "control_fasting.jl"))     # includes control_model.jl and control_mpc.jl; lr_plant lives there
using Test
using Random: Xoshiro

const Δ = Float64(CFG["schedule"]["decision_interval_days"])
const LEVELS = Float64.(CFG["dose"]["levels_data"])
rng = Xoshiro(20260929)

# a random state the plant can actually reach: random counts, then 1–3 intervals of random doses
# (hidden drug and stress compartments then sit in their biological ranges)
function random_state(P, rng)
    x = zeros(P.nx)
    x[P.iN] = 50 + rand(rng) * 0.6P.K; x[P.iC] = 50 + rand(rng) * 0.6P.K
    P.fp === nothing || (x[P.fp.iW] = 1.0)
    for _ in 1:rand(rng, 1:3)
        x = F(x, rand(rng) * 1.47, P, Δ; h=0.005)
    end
    x
end

# the held-out construction scenario B is built from (zoo_heldout.jl: heldout_ce1_model): the Stage 3 lines
# with the Stage 2 kills fixed, parameters = the crowding fitted to the untreated cocultures
function heldout_model_B()
    c3 = zgroup("s3").cands[1]; LN, LC = c3.sys.lines
    βof(c) = c.params[findfirst(==(:β), c.slots)]
    kof(c) = (i = findfirst(==(:k), c.slots); i === nothing ? nothing : c.params[i])
    # as zoo_heldout.jl mono_drug_line: Stage 2b kill and delay rate, no coculture drug modifier (φ, ρ, γ)
    function mono(L, c)
        s = merge(L.s, (β=Slot(0, βof(c)),))
        k = kof(c); k !== nothing && L.s.k.i > 0 && (s = merge(s, (k=Slot(0, k),)))
        L2 = merge(L, (s=s,))
        cc_on(L2) ? merge(L2, (cc=merge(L2.cc, (kf=Slot(0, 1.0), pr=Slot(0, 0.0), cs=Slot(0, 0.0))),)) : L2
    end
    LN2 = mono(LN, zgroup("s2-a2780").cands[1]); LC2 = mono(LC, zgroup("s2-a2780cis").cands[1])
    make_system_model(build_system([LN2, LC2], c3.sys.envs, c3.sys.alpha)), AN_STATE[:validation].ce1[1].params
end

@testset verbose = true "control plant and MPC" begin
    for (scen, P) in ((:A, scenario_plant(:A)), (:B, scenario_plant(:B)), (Symbol("A, low resource"), lr_plant(:A)), (Symbol("B, low resource"), lr_plant(:B)))
        W = default_weights(P)
        normal = scen in (:A, :B)

        @testset "scenario $scen: plant reproduces the fitted zoo model at 1 µM" begin
            g3 = zgroup("s3"); c3 = g3.cands[1]
            if normal
                # 50:50 Ce1 coculture (environment 5): zoo model vs control plant with u = 1 µM.
                # A: the Stage 3 fit itself; B: its held-out construction (Ce0 crowding, Stage 2 kills)
                e = 5; obs = [d.state_index for d in g3.datasets]
                model, p = scen == :A ? (c3.model, c3.params) : heldout_model_B()
                u0z = c3.u0
                solz = solve(ODEProblem(model, u0z, (0.0, 14.0), p), Tsit5(); reltol=1e-10, abstol=1e-10)
                x0 = x0_for(P, u0z[obs[2*e-1]], u0z[obs[2*e]])
                ts, xs = trajectory(x0, fill(1.0, 7), P, 2.0; dt=0.5)
                maxrel = maximum(max(abs(xs[i][P.iN] - solz(t)[obs[2*e-1]]) / max(solz(t)[obs[2*e-1]], 50),
                                     abs(xs[i][P.iC] - solz(t)[obs[2*e]]) / max(solz(t)[obs[2*e]], 50)) for (i, t) in enumerate(ts))
                @test maxrel < 1e-5
            elseif scen == Symbol("A, low resource")
                # Stage 4 fit (low resource + 1 µM, from the pre-exposed state): zoo model vs control plant
                c4 = zgroup("s4").cands[1]
                solz = solve(ODEProblem(c4.model, c4.u0, (0.0, 12.0), c4.params), Tsit5(); reltol=1e-10, abstol=1e-10)
                x0 = vcat(c4.u0, zeros(3), 1.0)            # + B, T, protection P, weight W of the fasting plant held at v = 1
                ts, xs = trajectory(x0, fill(1.0, 6), P, 2.0; dt=0.5)
                maxrel = maximum(max(abs(xs[i][P.iN] - solz(t)[P.iN]) / max(solz(t)[P.iN], 50),
                                     abs(xs[i][P.iC] - solz(t)[P.iC]) / max(solz(t)[P.iC], 50)) for (i, t) in enumerate(ts))
                @test maxrel < 1e-5
            end
        end

        normal && @testset "scenario $scen: flattened right-hand side equals the zoo right-hand side" begin
            @test P.spec !== nothing                      # the normal-medium plant uses the flattened equations
            e = 0.0
            for _ in 1:100
                x = random_state(P, rng); u = rand(rng) * 1.6
                d1 = zeros(P.nx); d2 = zeros(P.nx); plant_rhs!(d1, x, u, P); plant_rhs!(d2, x, u, P; ref=true)
                e = max(e, maximum(abs.(d1 .- d2) ./ max.(abs.(d2), 1.0)))
            end
            @test e < 1e-12
        end

        @testset "scenario $scen: RK4 transition vs adaptive Tsit5" begin
            for _ in 1:5
                x = random_state(P, rng); u = rand(rng) * 1.47
                xr = F(x, u, P, Δ)
                sol = solve(ODEProblem((du, y, p, t) -> plant_rhs!(du, y, u, P), x, (0.0, Δ)), Tsit5(); reltol=1e-11, abstol=1e-11)
                @test maximum(abs.(xr .- sol.u[end]) ./ max.(abs.(sol.u[end]), 1.0)) < 1e-5   # 1e-5: the Norton–Simon factor max(1 − x, 0) has a kink at capacity, where RK4 loses its smooth-case order
            end
        end

        @testset "scenario $scen: AD Jacobians vs central finite differences" begin
            for _ in 1:5
                x = random_state(P, rng); u = rand(rng) * 1.47
                _, A, B = linearize(x, u, P, Δ)
                z = vcat(x, u); nx = length(x)
                Jfd = zeros(nx, nx + 1)
                for j in 1:nx+1
                    h = 1e-6 * max(abs(z[j]), 1.0); zp = copy(z); zm = copy(z); zp[j] += h; zm[j] -= h
                    Jfd[:, j] = (F(zp[1:nx], zp[end], P, Δ) .- F(zm[1:nx], zm[end], P, Δ)) ./ (2h)
                end
                scale = max.(abs.(hcat(A, B)), 1e-3 * maximum(abs.(hcat(A, B))))
                @test maximum(abs.(hcat(A, B) .- Jfd) ./ scale) < 1e-4   # finite differences straddle the Norton–Simon kink at capacity
            end
        end

        @testset "scenario $scen: linearisation error is second order" begin
            x = random_state(P, rng); u = 0.7
            Fx, A, B = linearize(x, u, P, Δ)
            dx = randn(rng, length(x)) .* abs.(x) .* 0.05; du = 0.1       # relative: states stay in the valid (non-negative) domain
            err(ε) = maximum(abs.(F(x .+ ε .* dx, u + ε * du, P, Δ) .- (Fx .+ A * (ε .* dx) .+ B .* (ε * du))))
            e1, e2 = err(0.2), err(0.1)
            @test 3.0 < e1 / e2 < 5.0              # halving the step quarters the error
        end

        @testset "scenario $scen: adjoint gradient vs finite differences" begin
            x0 = x0_for(P, 300.0, 300.0); U = [0.2, 0.9, 0.5, 1.2]      # off u = 1 µM, the kink of the low-resource ε(u)
            g = cost_gradient(x0, U, P, Δ, W)
            gfd = [(Up = copy(U); Um = copy(U); Up[k] += 1e-6; Um[k] -= 1e-6;
                    (total_cost(x0, Up, P, Δ, W) - total_cost(x0, Um, P, Δ, W)) / 2e-6) for k in eachindex(U)]
            @test maximum(abs.(g .- gfd) ./ max.(abs.(gfd), 1e-6)) < 1e-4
        end

        @testset "scenario $scen: cost is linear in the burden, with the resistance premium (no progression term)" begin
            x = x0_for(P, 1000.0, 500.0); x2 = x0_for(P, 2000.0, 1000.0)
            Wc = default_weights(P; w=[1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
            c(x) = stage_cost(x, 0.0, P, Wc) - burden_penalty(x[P.iN] + x[P.iC], Wc) - relaxed_barrier(1.0, W.μT, W.δ)
            @test isapprox(c(x2), 2c(x); rtol=1e-12)                                   # doubling every count doubles the burden cost
            @test isapprox(c(x), (1000 + 2 * 500) / P.K; rtol=1e-12)                    # qR = 1: a resistant cell counts twice
        end

        @testset "scenario $scen: every stage and terminal cost is ≥ 0 (exact pruning relies on it)" begin
            Ws = [W, default_weights(P; bc_frac=1.0, tf=10.0), default_weights(P; bc_frac=0.0, tf=1.0), default_weights(P; umax=3.0)]
            ok = true
            for _ in 1:200
                x = random_state(P, rng); x[P.iN] *= 2rand(rng); x[P.iC] *= 2rand(rng)     # also far above B_cancer (quadratic branch)
                u = rand(rng) * 3.0
                for Wi in Ws
                    ok &= stage_cost(x, u, P, Wi) >= 0 && terminal_cost(x, P, Wi) >= 0
                end
            end
            @test ok
            # the per-term split used in the reports adds up to the stage cost along a schedule
            x0 = x0_for(P, 50.0, 50.0); U = [0.0, 1.47, 0.62, 1.0, 0.0]
            xs = rollout(x0, U, P, Δ)
            @test isapprox(cost_terms(P, x0, U, Δ, W).total, sum(stage_cost(xs[k+1], U[k], P, W) for k in eachindex(U)); rtol=1e-12)
        end

        @testset "scenario $scen: time-to-progression cost: a smooth step at B_prog, resistant cells priced at the end" begin
            # placed in step widths from B_prog, so the check holds for any threshold (K/2 since 2026-10-27)
            above = x0_for(P, W.Bp, W.Bp)
            if uses_lines(W)
                # "lines" rule (user, 2026-11-03): resistant at its threshold OR sensitive at its threshold
                lN = W.BpN - 15W.sP; lC = W.BpC - 15W.sP; hN = W.BpN + 15W.sP; hC = W.BpC + 15W.sP
                @test progression_cost(lN, lC, W) < 1e-5 * W.qP            # both lines under: ≈ 0
                @test progression_cost(lN, hC, W) > (1 - 1e-5) * W.qP      # resistant over: progressed
                @test progression_cost(hN, lC, W) > (1 - 1e-5) * W.qP      # sensitive over: progressed
                @test isapprox(progression_cost(W.BpN - 40W.sP, W.BpC, W), W.qP / 2; rtol=1e-6)   # half way at the resistant threshold
                @test W.BpC == CFG["progression"]["resistant_threshold_fraction_of_K"] * P.K
                @test 0 < W.BpN < sensitive_plateau(P)                   # reachable without treatment
                @test prog_index(W.BpN, 0.0, W) == 1.0 && prog_index(0.0, W.BpC, W) == 1.0
            else
                lo = (W.Bp - 15W.sP) / 2; hi = (W.Bp + 15W.sP) / 2
                @test progression_cost(lo, lo, W) < 1e-6 * W.qP              # ≈ 0 well under the threshold
                @test progression_cost(hi, hi, W) > (1 - 1e-6) * W.qP        # ≈ qP well above it
                at = x0_for(P, 0.5W.Bp, 0.5W.Bp)
                @test isapprox(progression_cost(at[P.iN], at[P.iC], W), W.qP / 2; rtol=1e-12)   # half way at the threshold
            end
            W0 = with_fields(W; qRT=0.0)
            @test isapprox(terminal_cost(above, P, W) - terminal_cost(above, P, W0), W.qRT * above[P.iC] / W.K; rtol=1e-12)
            @test W.Bp == CFG["progression"]["threshold_fraction_of_K"] * P.K
        end

        @testset "scenario $scen: blocked (move-blocking) search equals brute force over blocked plans" begin
            x0 = x0_for(P, 300.0, 300.0); blk = [1, 2]
            sol = discrete_solve(x0, P, Δ, blk, LEVELS, W)
            brute = minimum(plan_cost(x0, expand_blocks([a, b], blk), P, Δ, W) for a in LEVELS, b in LEVELS)
            @test isapprox(sol.J, brute; rtol=1e-10)
            @test isapprox(plan_cost(x0, expand_blocks(sol.U, blk), P, Δ, W), sol.J; rtol=1e-10)   # the cost counts every interval
            @test discrete_solve(x0, P, Δ, [1, 1], LEVELS, W).J ≈ discrete_solve(x0, P, Δ, 2, LEVELS, W).J   # single blocks = no blocking
        end

        @testset "scenario $scen: discrete search equals brute force" begin
            x0 = x0_for(P, 300.0, 300.0); H = 3
            sol = discrete_solve(x0, P, Δ, H, LEVELS, W)
            brute = minimum(plan_cost(x0, [a, b, c], P, Δ, W) for a in LEVELS, b in LEVELS, c in LEVELS)
            @test isapprox(sol.J, brute; rtol=1e-10)
            @test isapprox(discrete_solve(x0, P, Δ, H, LEVELS, W; Uinit=[1.47, 0.0, 0.62]).J, brute; rtol=1e-10)   # a warm-start bound never changes the optimum
        end

        @testset "scenario $scen: SQP respects bounds, does not increase the cost, and is optimal" begin
            x0 = x0_for(P, 300.0, 300.0); U0 = fill(0.5, 5)
            sol = sqp_solve(x0, U0, P, Δ, W)
            @test all(0 .<= sol.U .<= W.umax + 1e-12)
            @test sol.J <= total_cost(x0, U0, P, Δ, W) + 1e-9
            cl = closed_loop(P, P, x0, Δ, 4, 4, W; mode=:continuous)
            @test all(0 .<= cl.U .<= W.umax + 1e-12)
            # over [0, 1.47] µM the measured doses are a subset; the time-to-progression step makes the cost non-convex,
            # so SQP (local) is started at the exact discrete optimum, as the continuous closed loop does, and must not end above it
            H = 4; d4 = discrete_solve(x0, P, Δ, H, LEVELS, W); s4 = sqp_solve(x0, d4.U, P, Δ, W)
            @test s4.J <= total_cost(x0, d4.U, P, Δ, W) * (1 + 1e-6) + 1e-9   # SQP's own (end-of-interval) cost
            # first-order optimality: the projected gradient vanishes at the SQP solution (normal medium only: the
            # low-resource plant has a kink in ε(u) at u = 1 µM, where SQP may stop at a non-smooth point)
            if normal
                g = cost_gradient(x0, s4.U, P, Δ, W)
                @test maximum(abs.(clamp.(s4.U .- g, 0.0, W.umax) .- s4.U)) < 2e-3
            end
        end
    end

    @testset "nominal closed loop (fixed look-ahead) equals replaying its own schedule" begin
        P = scenario_plant(:A); W = default_weights(P); s = nominal_seeding("50:50"); x0 = x0_for(P, s.N0, s.C0)
        n, H = 4, 3
        cl = closed_loop(P, P, x0, Δ, n, H, W; mode=:discrete, levels=LEVELS)
        xs = rollout(x0, cl.U, P, Δ)
        @test xs[end] == cl.x                                              # the planner's reconstructed state is the plant's
        @test all(discrete_solve(xs[k], P, Δ, H, LEVELS, W).U[1] == cl.U[k] for k in 1:n)   # each decision = first move of the exact plan
        @test closed_loop(P, P, x0, Δ, n, H, W; mode=:discrete, levels=LEVELS, noise=ones(n, 2)).U == cl.U   # noise factors of 1 = exact counts
        # the fixed look-ahead keeps planning past the end: the last decision sees H intervals, the old shrinking one only 1
        @test closed_loop(P, P, x0, Δ, n, H, W; mode=:discrete, levels=LEVELS, shrink=true).U[1:n-H+1] == cl.U[1:n-H+1]
    end

    @testset "blocked look-ahead: the closed loop equals replaying its own schedule" begin
        P = scenario_plant(:A); W = default_weights(P); s = nominal_seeding("50:50"); x0 = x0_for(P, s.N0, s.C0)
        n, blk = 3, [1, 2]
        cl = closed_loop(P, P, x0, Δ, n, blk, W; mode=:discrete, levels=LEVELS)
        xs = rollout(x0, cl.U, P, Δ)
        @test xs[end] == cl.x
        @test all(discrete_solve(xs[k], P, Δ, blk, LEVELS, W).U[1] == cl.U[k] for k in 1:n)   # each decision = first move of the blocked plan
        @test cut_blocks([1, 1, 2, 2, 4, 4], 5) == [1, 1, 2, 1] && cut_blocks([1, 2], 10) == [1, 2]
        @test Int.(CFG["schedule"]["lookahead_blocks"])[1] == 1                         # the applied dose is decided at full resolution
    end

    @testset "time to progression on a trajectory" begin
        ts = collect(0.0:1.0:10.0); hi, lo = 9.0, 1.0
        # above on days 1–2 only (reversed within 3 days): not progression; above from day 4 to the end: progression on day 4
        tp = time_to_progression(ts, [lo, hi, hi, lo, hi, hi, hi, hi, hi, hi, hi], 4.0; sustain=3.0)
        @test tp.first == 1.0 && tp.loss == 4.0 && tp.days == 9.0
        # above on days 1–6, then one dip: progression on day 1 (the spell lasted ≥ 3 days), the dip does not restart it
        @test time_to_progression(ts, [lo, hi, hi, hi, hi, hi, hi, lo, hi, hi, hi], 4.0; sustain=3.0).loss == 1.0
        @test isnan(time_to_progression(ts, [lo, hi, hi, lo, lo, lo, lo, lo, lo, lo, lo], 4.0; sustain=3.0).loss)   # reversed: not progressed
        @test time_to_progression(ts, [lo, lo, lo, lo, lo, lo, lo, lo, lo, hi, hi], 4.0; sustain=3.0).loss == 9.0  # runs to the end
        @test isnan(time_to_progression(ts, fill(lo, 11), 4.0).first)
        @test time_to_progression(ts, fill(hi, 11), 4.0).loss == 0.0
    end

    @testset "measurement noise: log-normal, mean 1, the configured CV, paired by draw" begin
        z = obs_noise(1, 20000)
        cv = Float64(CFG["observation"]["count_cv"])
        @test isapprox(mean(z), 1.0; atol=0.005) && isapprox(std(z) / mean(z), cv; rtol=0.05)
        @test obs_noise(7, 5) == obs_noise(7, 5) && obs_noise(7, 5) != obs_noise(8, 5)
    end

    @testset "nominal seeding: 100 cells split by the ratio" begin
        s = nominal_seeding("25:75")
        @test s.N0 ≈ 0.25 * CFG["initial"]["seeding_total"] && s.C0 ≈ 0.75 * CFG["initial"]["seeding_total"]
    end

    @testset "long runs: 30 days, 14-day look-ahead re-planned every 2 days" begin
        P = scenario_plant(:A); W = default_weights(P)
        days = 30.0; n = round(Int, days / Δ); H = round(Int, CFG["schedule"]["horizon_days"] / Δ)   # a run longer than the look-ahead (the 120-day run uses the same code, 4× longer)
        @test n == 15 && H == 7
        x0 = x0_for(P, 300.0, 300.0)
        cl = closed_loop(P, P, x0, Δ, n, H, W; mode=:discrete, levels=LEVELS)
        @test length(cl.U) == n                                # one dose per 2-day interval up to day 30
        @test all(u -> any(isapprox(u, l) for l in LEVELS), cl.U)   # only measured doses
        o = outcomes(P, x0, cl.U, Δ, W)
        @test isapprox(o.ts[end], days; atol=1e-9)              # the trajectory reaches day 30
        @test all(isfinite, o.N) && all(isfinite, o.C) && all(>=(0), o.N) && all(>=(0), o.C)
        @test isapprox(o.cum_dose, sum(cl.U) * Δ)
        @test 0 <= o.days_over <= days
    end

    @testset "relaxed barrier" begin
        μ, δ = 0.05, 0.02
        @test relaxed_barrier(0.5, μ, δ) ≈ -μ * log(0.5)
        @test isapprox(relaxed_barrier(δ + 1e-9, μ, δ), relaxed_barrier(δ - 1e-9, μ, δ); atol=1e-6)     # continuous
        d(z) = (relaxed_barrier(z + 1e-7, μ, δ) - relaxed_barrier(z - 1e-7, μ, δ)) / 2e-7
        @test isapprox(d(δ + 1e-5), d(δ - 1e-5); rtol=1e-2)                                             # C¹ at δ
        @test relaxed_barrier(-0.5, μ, δ) > relaxed_barrier(0.0, μ, δ) > relaxed_barrier(δ, μ, δ)         # steep, finite outside
        @test isfinite(relaxed_barrier(-1.0, μ, δ))
        @test relaxed_barrier(1.0, μ, δ) == 0                                                            # no burden penalty: B_cancer = Inf gives z = 1
    end

    # host limits (user, 2026-11-12): the inside of the barrier is a tie-breaker, a violation keeps the full penalty
    @testset "host barrier: tie-break inside, full penalty outside" begin
        μ, δ, s = 0.05, 0.02, 0.01
        for z in (0.9, 0.5, 0.03, 0.01, 0.0)
            @test host_barrier(z, μ, δ, s) ≈ s * relaxed_barrier(z, μ, δ)                     # inside the limit: scaled
        end
        @test host_barrier(-0.1, μ, δ, s) ≈ s * relaxed_barrier(-0.1, μ, δ) + μ / (2δ^2) * 0.01   # outside: + the full quadratic
        @test host_barrier(-0.1, μ, δ, s) > 0.5 * μ / (2δ^2) * 0.01                               # a violation is never cheap
        @test host_barrier(1.0, μ, δ, s) == 0 && all(host_barrier(z, μ, δ, s) >= 0 for z in -1:0.01:1)
        @test isapprox(host_barrier(1e-9, μ, δ, s), host_barrier(-1e-9, μ, δ, s); atol=1e-8)       # continuous at the limit
        P = scenario_plant(:A); W = default_weights(P)
        @test W.sB == CFG["progression"]["tiebreak_scale"] && default_weights(P; w=[1.0, 1.0, 1.0, 0.1, 0.0, 0.0]).sB == 1.0
    end

    # starting counts (user, 2026-11-12): the Stage 3a fitted start of each ratio's treated (Ce1) series by default
    @testset "starting counts: the fitted start of the treated coculture" begin
        c3 = zgroup("s3").cands[1]
        for ratio in ("25:75", "50:50", "75:25")
            s = fitted_start(ratio, true); e = coculture_env(ratio, true)
            @test s.N0 == c3.u0[e.iN] && s.C0 == c3.u0[e.iC] && s.N0 > 0 && s.C0 > 0
            @test start_counts(ratio) == (initial_source() == "fitted_treated" ? s : nominal_seeding(ratio))
            u = fitted_start(ratio, false); n = nominal_seeding(ratio)
            @test u.C0 ≈ n.C0                                    # the resistant start is the nominal seeding in both series
        end
    end

    # The 14-day plant against the treated cocultures (user, 2026-11-12). Scenario A is the Stage 3a fit itself, so
    # from the fitted Ce1 start at a constant 1 µM it must (1) equal the fitted model at day 14 (fixed-step RK4 vs
    # Tsit5: 0.5%) and (2) sit within 2 × the fit's residual RMSE of that series of the measured day-14 counts
    # (a ~95% band for a residual of that size; the residuals are the only measure of the data's scatter here).
    @testset "plant vs data: scenario A at 1 µM reproduces the treated cocultures at day 14" begin
        rows = [r for r in plant_vs_data() if r.scenario == :A && startswith(r.condition, "Ce1")]
        @test length(rows) == 3
        for r in rows
            d = coculture_day14(r.ratio, true)
            @test isapprox(r.N, d.fitN; rtol=5e-3) && isapprox(r.C, d.fitC; rtol=5e-3)
            @test abs(r.N - r.N_data) <= 2 * r.rmseN && abs(r.C - r.C_data) <= 2 * r.rmseC
        end
    end

    # PRD 2026-10-27: the four burden-penalty shapes at B_cancer ([cancer] barrier, swept by run_control_sweep.jl)
    @testset "burden barrier shapes" begin
        P = scenario_plant(:A)
        for b in BARRIERS
            W = default_weights(P; barrier=b)
            @test W.barrier == b && W.Bc ≈ CFG["cancer"]["B_cancer_fraction_of_K"] * P.K
            lo, at, hi = burden_penalty(0.1W.Bc, W), burden_penalty(W.Bc, W), burden_penalty(1.1W.Bc, W)
            @test 0 <= lo < 0.01                                          # ≈ 0 well below the threshold
            @test all(isfinite, (lo, at, hi)) && all(>=(0), (lo, at, hi))  # ≥ 0: the exact search can still prune
            b == :none ? (@test hi == 0) : (@test hi > lo && hi > at)      # rises above it
            g = ForwardDiff.derivative(t -> burden_penalty(t, W), 1.05W.Bc)
            @test isfinite(g) && (b == :none || g > 0)                    # differentiable for SQP
        end
        Wl = default_weights(P; barrier=:relaxed_log)
        @test burden_penalty(0.4P.K, Wl) == relaxed_barrier(1 - 0.4P.K / Wl.Bc, CFG["weights"]["mu_cancer"], CFG["weights"]["barrier_delta"])
        @test default_weights(P).barrier == Symbol(CFG["cancer"]["barrier"])   # "none" since 2026-11-01
        @test_throws ErrorException default_weights(P; barrier=:cubic)
        # hard: any interval at or above B_cancer costs more than a whole look-ahead of every other term
        Wh = default_weights(P; barrier=:hard)
        @test burden_penalty(Wh.Bc, Wh) > 20 * (Wh.qP + Wh.qRT + Wh.qB + Wh.ru) && burden_penalty(Wh.Bc - 0.02P.K, Wh) < 1
        # so the exact search avoids a crossing plan when one exists: from just under B_cancer, no drug crosses
        # and the plan that stays under (if any) is chosen
        x0 = x0_for(P, 0.25Wh.Bc, 0.2Wh.Bc)
        sol = discrete_solve(x0, P, Δ, 3, LEVELS, Wh)
        cross(U) = any(x -> x[P.iN] + x[P.iC] >= Wh.Bc, rollout(x0, U, P, Δ)[2:end])
        plans = vec([[a, b, c] for a in LEVELS, b in LEVELS, c in LEVELS])
        any(!cross(U) for U in plans) && @test !cross(sol.U)
    end

    # user, 2026-11-01: time to progression first, and no MPC decisions once the result is known
    @testset "time-to-progression objective and early stop" begin
        P = scenario_plant(:A); W = default_weights(P); x0 = x0_for(P, 300.0, 300.0); LOOK = Int.(CFG["schedule"]["lookahead_blocks"])
        # the progression cost of an interval is the time-averaged step: between 0 and qP, and with one
        # piece it is the end-of-interval form
        y, pc = step_progressed(x0, 0.0, P, Δ, W)
        @test isapprox(y, F(x0, 0.0, P, Δ); rtol=1e-10)
        @test 0 <= pc <= W.qP
        @test plan_cost(x0, [0.0, 1.0], P, Δ, W) >= 0
        @test CFG["progression"]["tiebreak_scale"] < 1 && W.qB ≈ CFG["weights"]["qB"] * CFG["progression"]["tiebreak_scale"]
        @test default_weights(P; w=[1.0, 1.0, 1.0, 0.1, 0.0, 0.0]).qB == 1.0      # no progression term: not scaled
        # progression confirmed: stops after `confirm` successive decisions above B_prog; arrays keep the run length
        n = round(Int, 30 / Δ); s = nominal_seeding("50:50"); x1 = x0_for(P, s.N0, s.C0)
        cl = closed_loop(P, P, x1, Δ, n, LOOK, W; mode=:discrete, levels=LEVELS, stop_kw(W, Δ)...)
        @test length(cl.U) == n
        if cl.stop_reason == "progression confirmed"
            o = outcomes(P, x1, cl.U, Δ, W; t_stop=tstop(cl))
            @test o.ttp <= cl.stop_day && o.stop_day == cl.stop_day && isfinite(o.N_prog)
            # the last stop_confirm decisions were all at or above B_prog
            xs = rollout(x1, cl.U, P, Δ); k = round(Int, cl.stop_day / Δ)
            @test all(prog_index(xs[j+1][P.iN], xs[j+1][P.iC], W) >= 1 for j in k-stop_confirm(Δ)+1:k)
        end
        # outlasted: with every comparator progressed by day 2, a run under B_prog on day 12 stops there
        cl2 = closed_loop(P, P, x1, Δ, n, LOOK, W; mode=:discrete, levels=LEVELS, Bp=Inf, confirm=0, settle_day=12.0)
        @test cl2.stop_reason == "outlasted every comparator" && cl2.stop_day == 12.0 && length(cl2.U) == n
        @test settle_after([5.0, NaN]) == Inf && settle_after([5.0, 7.0]) == 7.0 + CFG["progression"]["settle_margin_days"]
        # no stopping rule: the old behaviour
        @test isnan(closed_loop(P, P, x1, Δ, 3, 3, W; mode=:discrete, levels=LEVELS).stop_day)
    end

    @testset "look-ahead that ends with the run" begin
        LOOK = Int.(CFG["schedule"]["lookahead_blocks"])
        @test cut_blocks(LOOK, 1) == [1]                                   # the last decision plans one interval
        @test all(sum(cut_blocks(LOOK, m)) == min(m, sum(LOOK)) for m in 1:20)
        P = scenario_plant(:A); W = default_weights(P); n = round(Int, CFG["schedule"]["closed_loop_days"] / Δ)
        x0 = x0_for(P, 50.0, 50.0)
        cl = closed_loop(P, P, x0, Δ, n, LOOK, W; mode=:discrete, levels=LEVELS, shrink=true)
        @test length(cl.U) == n && all(u -> any(isapprox(u, l) for l in LEVELS), cl.U)
        # the last decision equals the exact one-interval problem from the state it starts in
        xl = foldl((x, u) -> F(x, u, P, Δ), cl.U[1:end-1]; init=copy(x0))
        @test cl.U[end] == discrete_solve(xl, P, Δ, [1], LEVELS, W).U[1]
    end
    @testset "kidney toxicity (xie2020): calibration and dynamics" begin
        if TOX_MODEL == "xie2020"
            t = default_tox(); P = scenario_plant(:A)
            # γ solves Emax·γ/(EC50 + γ) = Ē, and the periodic-steady-state peak sits above the mean
            @test t.Emax * t.γ / (t.EC50 + t.γ) ≈ XIE.Ebar rtol = 1e-10
            @test XIE.Elim > XIE.Ebar > 0 && XIE.ΔE > 0 && tox_limit() == XIE.Elim
            # a constant 1 µM drives E to Ē (k_N ≈ 1/6 per day: 60 days ≈ 10 time constants), with no dose E decays
            x = x0_for(P, 50.0, 50.0)
            for _ in 1:30; x = F(x, 1.0, P, 2.0); end
            @test x[P.iT] ≈ XIE.Ebar rtol = 1e-3
            for _ in 1:10; x = F(x, 0.0, P, 2.0); end
            @test x[P.iT] ≈ XIE.Ebar * exp(-t.kN * 20) rtol = 1e-3
            # the drug burden B no longer feeds the toxicity, and the limit is the weights' limit
            @test default_weights(P).Btox == XIE.Elim
            @test xie_drive(t, 0.0) == 0 && xie_drive(merge(t, (γ=3t.γ,)), 1.0) > xie_drive(t, 1.0)
        else
            @test default_tox().Elim == CFG["toxicity"]["B_tox"]
        end
    end
end
