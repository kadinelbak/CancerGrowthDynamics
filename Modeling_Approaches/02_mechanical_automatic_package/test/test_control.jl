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
    LN2 = merge(LN, (s=merge(LN.s, (β=Slot(0, βof(zgroup("s2-a2780").cands[1])),)),))
    LC2 = merge(LC, (s=merge(LC.s, (β=Slot(0, βof(zgroup("s2-a2780cis").cands[1])), λ=Slot(0, 1.0))),))   # λ_C = 1: the Stage 2b kill
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
            c(x) = stage_cost(x, 0.0, P, default_weights(P; w=[1.0, 1.0, 0.0, 0.0, 0.0, 0.0])) - relaxed_barrier(1 - (x[P.iN] + x[P.iC]) / W.Bc, W.μC, W.δ) - relaxed_barrier(1.0, W.μT, W.δ)
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
            below = x0_for(P, 0.25W.Bp, 0.25W.Bp); above = x0_for(P, W.Bp, W.Bp)   # 2 B_prog: 40 step widths above
            @test progression_cost(below[P.iN], below[P.iC], W) < 1e-6 * W.qP              # ≈ 0 well under the threshold
            @test progression_cost(above[P.iN], above[P.iC], W) > (1 - 1e-6) * W.qP        # ≈ qP well above it
            at = x0_for(P, 0.5W.Bp, 0.5W.Bp)
            @test isapprox(progression_cost(at[P.iN], at[P.iC], W), W.qP / 2; rtol=1e-12)   # half way at the threshold
            W0 = with_fields(W; qRT=0.0)
            @test isapprox(terminal_cost(above, P, W) - terminal_cost(above, P, W0), W.qRT * above[P.iC] / W.K; rtol=1e-12)
            @test W.Bp == CFG["progression"]["threshold_fraction_of_K"] * P.K
        end

        @testset "scenario $scen: blocked (move-blocking) search equals brute force over blocked plans" begin
            x0 = x0_for(P, 300.0, 300.0); blk = [1, 2]
            sol = discrete_solve(x0, P, Δ, blk, LEVELS, W)
            brute = minimum(total_cost(x0, expand_blocks([a, b], blk), P, Δ, W) for a in LEVELS, b in LEVELS)
            @test isapprox(sol.J, brute; rtol=1e-10)
            @test isapprox(total_cost(x0, expand_blocks(sol.U, blk), P, Δ, W), sol.J; rtol=1e-10)   # the cost counts every interval
            @test discrete_solve(x0, P, Δ, [1, 1], LEVELS, W).J ≈ discrete_solve(x0, P, Δ, 2, LEVELS, W).J   # single blocks = no blocking
        end

        @testset "scenario $scen: discrete search equals brute force" begin
            x0 = x0_for(P, 300.0, 300.0); H = 3
            sol = discrete_solve(x0, P, Δ, H, LEVELS, W)
            brute = minimum(total_cost(x0, [a, b, c], P, Δ, W) for a in LEVELS, b in LEVELS, c in LEVELS)
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
            @test s4.J <= d4.J * (1 + 1e-6) + 1e-9
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
end
