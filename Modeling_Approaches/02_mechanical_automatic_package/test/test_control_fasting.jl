# Tests for the two-input (dose + fasting) controller:  julia --project=. test/test_control_fasting.jl
include(joinpath(@__DIR__, "..", "control_fasting.jl"))
using Test
using Random: Xoshiro

const Δ = Float64(CFG["schedule"]["decision_interval_days"])
rng = Xoshiro(20260930)

# a reachable state: random counts, then 1–3 intervals of random doses and fasting
function random_state2(FP, rng)
    x = x0_fast(FP, 50 + rand(rng) * 0.6FP.K, 50 + rand(rng) * 0.6FP.K)
    x = F2(x, [0.1 + rand(rng) * 1.37, 1.0], FP, Δ; h=0.005)   # one dosed, fasting interval: every state > 0 (no clamp at 0)
    for _ in 1:rand(rng, 0:2)
        x = F2(x, [rand(rng) * 1.47, Float64(rand(rng, 0:1))], FP, Δ; h=0.005)
    end
    x
end

@testset verbose = true "dose + fasting control" begin
    for scen in (:A, :B)
        FP = fast_plant(scen); P = scenario_plant(scen); W = weights2(FP)

        @testset "scenario $scen: fed (v = 0) is exactly the normal-medium plant" begin
            x = x0_for(P, 300.0, 300.0); xf = x0_fast(FP, 300.0, 300.0); e = 0.0
            for u in (0.0, 0.62, 1.0, 1.47, 0.3)
                x = F(x, u, P, Δ); xf = F2(xf, [u, 0.0], FP, Δ)
                e = max(e, abs(x[P.iN] - xf[FP.iN]) / x[P.iN], abs(x[P.iC] - xf[FP.iC]) / x[P.iC])
            end
            @test e < 1e-12
        end

        @testset "scenario $scen: fast right-hand side equals the zoo right-hand side" begin
            e = 0.0
            for mode in (:interp, :stage4)
                FPm = fast_plant(scen; eps_mode=mode)
                for _ in 1:100
                    x = random_state2(FPm, rng); w = [rand(rng) * 1.6, rand(rng)]
                    d1 = zeros(FPm.nx); d2 = zeros(FPm.nx); fast_rhs!(d1, x, w, FPm); fast_rhs_ref!(d2, x, w, FPm)
                    e = max(e, maximum(abs.(d1 .- d2) ./ max.(abs.(d2), 1.0)))
                end
            end
            @test e < 1e-12
        end

        @testset "scenario $scen: RK4 transition vs adaptive Tsit5" begin
            for _ in 1:5
                x = random_state2(FP, rng); w = [rand(rng) * 1.47, Float64(rand(rng, 0:1))]
                sol = solve(ODEProblem((du, y, p, t) -> fast_rhs!(du, y, w, FP), x, (0.0, Δ)), Tsit5(); reltol=1e-11, abstol=1e-11)
                @test maximum(abs.(F2(x, w, FP, Δ) .- sol.u[end]) ./ max.(abs.(sol.u[end]), 1.0)) < 1e-6
            end
            x = random_state2(FP, rng); y = zeros(FP.nx)
            @test F2!(y, x, (1.0, 1.0), FP, Δ, RK4Buf(FP.nx)) ≈ F2(x, [1.0, 1.0], FP, Δ)
        end

        @testset "scenario $scen: AD Jacobians vs central finite differences" begin
            for _ in 1:5
                x = random_state2(FP, rng); w = [0.1 + 0.8rand(rng), 0.2 + 0.6rand(rng)]   # u < u_ref: away from the kink of min(u/u_ref, 1)
                _, A, B = linearize2(x, w, FP, Δ)
                z = vcat(x, w); nx = length(x); Jfd = zeros(nx, nx + 2)
                for j in 1:nx+2
                    h = 1e-6 * max(abs(z[j]), 1.0); zp = copy(z); zm = copy(z); zp[j] += h; zm[j] -= h
                    Jfd[:, j] = (F2(zp[1:nx], zp[nx+1:end], FP, Δ) .- F2(zm[1:nx], zm[nx+1:end], FP, Δ)) ./ (2h)
                end
                scale = max.(abs.(hcat(A, B)), 1e-3 * maximum(abs.(hcat(A, B))))
                @test maximum(abs.(hcat(A, B) .- Jfd) ./ scale) < 1e-5
            end
        end

        @testset "scenario $scen: adjoint gradient vs finite differences" begin
            x0 = x0_fast(FP, 300.0, 300.0); U = [0.2 0.3; 1.0 0.0; 0.5 0.7; 1.2 0.1]
            g, _ = grad_hess2(x0, U, FP, Δ, W); z = vec(permutedims(U))
            Um(z) = permutedims(reshape(z, 2, size(U, 1)))
            fd = [(total_cost2(x0, Um(z .+ 1e-6 .* (1:length(z) .== i)), FP, Δ, W) - total_cost2(x0, Um(z .- 1e-6 .* (1:length(z) .== i)), FP, Δ, W)) / 2e-6 for i in eachindex(z)]
            @test maximum(abs.(fd .- g)) < 1e-5 * max(1.0, maximum(abs.(g)))
        end

        @testset "scenario $scen: exact search equals brute force and obeys the fasting rules" begin
            x0 = x0_fast(FP, 300.0, 300.0); H = 3; best = Inf; bestU = nothing
            levels = [(u, v) for v in (0.0, 1.0) for u in Float64.(CFG["dose"]["levels_data"])]
            for a in levels, b in levels, c in levels
                vs = (a[2], b[2], c[2])
                any(vs[i] > 0 && vs[i+1] > 0 for i in 1:2) && continue              # max 1 consecutive, min gap 1
                U = [a[1] a[2]; b[1] b[2]; c[1] c[2]]; J = total_cost2(x0, U, FP, Δ, W)
                J < best && (best = J; bestU = U)
            end
            s = discrete_solve2(x0, FP, Δ, H, W)
            @test s.J ≈ best rtol = 1e-10
            @test s.U == bestU
            s2 = discrete_solve2(x0, FP, Δ, H, W; Uinit=[1.0 0.0; 0.0 0.0; 1.0 0.0])      # a warm start never changes the optimum
            @test s2.J ≈ best rtol = 1e-10
        end

        @testset "scenario $scen: SQP respects bounds and does not increase the cost" begin
            x0 = x0_fast(FP, 300.0, 300.0); U0 = [0.5 0.2; 0.5 0.2; 0.5 0.2; 0.5 0.2]
            s = sqp_solve2(x0, U0, FP, Δ, W)
            @test all(0 .<= s.U[:, 1] .<= W.umax + 1e-12) && all(0 .<= s.U[:, 2] .<= 1 + 1e-12)
            @test s.J <= total_cost2(x0, U0, FP, Δ, W) + 1e-12
        end
    end

    @testset "host states behave as specified" begin
        FP = fast_plant(:A); h = FP.host
        x = x0_fast(FP, 300.0, 300.0); x1 = F2(x, [0.0, 1.0], FP, Δ)            # one 2-day fast
        @test x1[FP.iW] ≈ 1 - h.k_W * Δ rtol = 1e-9                              # dW/dt = −k_W while fasting
        @test x1[FP.iP] ≈ 1 - exp(-h.k_P * Δ) rtol = 1e-6                        # dP/dt = k_P(1 − P) (RK4 accuracy)
        x2 = F2(x1, [0.0, 0.0], FP, Δ)                                            # refeeding
        @test 1 - x2[FP.iW] ≈ (1 - x1[FP.iW]) * exp(-h.r_W * Δ) rtol = 1e-6
        # protection lowers the toxicity gained from the same dose
        xa = F2(F2(x, [0.0, 1.0], FP, Δ), [1.47, 1.0], FP, Δ); xb = F2(F2(x, [0.0, 0.0], FP, Δ), [1.47, 0.0], FP, Δ)
        @test xa[FP.iT] < xb[FP.iT]
        FP0 = fast_plant(:A; host=merge(h, (p_H=0.0,)))
        xc = F2(F2(x, [0.0, 1.0], FP0, Δ), [1.47, 1.0], FP0, Δ)
        @test xc[FP0.iT] ≈ xb[FP.iT] rtol = 1e-12                                 # p_H = 0: no protection
    end

    @testset "stress-recovery variant only changes the fed (refeeding) phase" begin
        FP = fast_plant(:A); FPk = fast_plant(:A; koff=2.0)
        x = random_state2(FP, rng); d1 = zeros(FP.nx); d2 = zeros(FP.nx)
        fast_rhs!(d1, x, [0.7, 1.0], FP); fast_rhs!(d2, x, [0.7, 1.0], FPk)
        @test d1 == d2                                                            # fasting: identical
        x1 = F2(x, [0.0, 0.0], FP, Δ); x2 = F2(x, [0.0, 0.0], FPk, Δ)
        S = [M for M in FP.spec[2].mods if !M.input_D && M.kernel_chain == 3]
        o = only(S).o
        @test sum(x2[o+1:o+3]) < sum(x1[o+1:o+3])                                 # refeeding: total stress (d/dt = −k·S_3) fades faster
    end

    @testset "Stage 4 is reproduced at 1 µM + low resource" begin
        FP = fast_plant(:A); c4 = zgroup("s4").cands[1]
        sol = solve(ODEProblem(c4.model, c4.u0, (0.0, 12.0), c4.params), Tsit5(); reltol=1e-10, abstol=1e-10)
        x = vcat(c4.u0, zeros(3), 1.0); mr = 0.0
        for k in 1:6
            x = F2(x, [1.0, 1.0], FP, Δ); s = sol(Δ * k)
            mr = max(mr, abs(x[FP.iN] - s[FP.iN]) / max(s[FP.iN], 50), abs(x[FP.iC] - s[FP.iC]) / max(s[FP.iC], 50))
        end
        @test mr < 1e-5
    end
end
