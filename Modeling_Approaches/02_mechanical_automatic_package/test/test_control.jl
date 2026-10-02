# Tests for the optimal-control layer:  julia --project=. test/test_control.jl
include(joinpath(@__DIR__, "..", "control_model.jl"))
include(joinpath(@__DIR__, "..", "control_mpc.jl"))
using Test
using Random: Xoshiro

const Δ = Float64(CFG["schedule"]["decision_interval_days"])
rng = Xoshiro(20260929)

# a random state the plant can actually reach: random counts, then 1–3 intervals of random doses
# (hidden drug and stress compartments then sit in their biological ranges)
function random_state(P, rng)
    x = zeros(P.nx)
    x[P.iN] = 50 + rand(rng) * 0.6P.K; x[P.iC] = 50 + rand(rng) * 0.6P.K
    for _ in 1:rand(rng, 1:3)
        x = F(x, rand(rng) * 1.47, P, Δ; h=0.005)
    end
    x
end

@testset verbose = true "control plant and MPC" begin
    for (scen, P) in ((:A, scenario_plant(:A)), (:B, scenario_plant(:B)), (Symbol("A, low resource"), lr_plant(:A)), (Symbol("B, low resource"), lr_plant(:B)))
        W = default_weights(P)

        @testset "scenario $scen: plant reproduces the fitted zoo model at 1 µM" begin
            g3 = zgroup("s3"); c3 = g3.cands[1]
            if scen == :A
                # 50:50 Ce1 coculture (environment 5): zoo model vs control plant with u = 1 µM
                e = 5; obs = [d.state_index for d in g3.datasets]
                u0z = c3.u0
                solz = solve(ODEProblem(c3.model, u0z, (0.0, 14.0), c3.params), Tsit5(); reltol=1e-10, abstol=1e-10)
                x0 = x0_for(P, u0z[obs[2*e-1]], u0z[obs[2*e]])
                ts, xs = trajectory(x0, fill(1.0, 7), P, 2.0; dt=0.5)
                maxrel = maximum(max(abs(xs[i][P.iN] - solz(t)[obs[2*e-1]]) / max(solz(t)[obs[2*e-1]], 50),
                                     abs(xs[i][P.iC] - solz(t)[obs[2*e]]) / max(solz(t)[obs[2*e]], 50)) for (i, t) in enumerate(ts))
                @test maxrel < 1e-5
            elseif scen == Symbol("A, low resource")
                # Stage 4 fit (low resource + 1 µM, from the pre-exposed state): zoo model vs control plant
                c4 = zgroup("s4").cands[1]
                solz = solve(ODEProblem(c4.model, c4.u0, (0.0, 12.0), c4.params), Tsit5(); reltol=1e-10, abstol=1e-10)
                x0 = vcat(c4.u0, 0.0, 0.0)
                ts, xs = trajectory(x0, fill(1.0, 6), P, 2.0; dt=0.5)
                maxrel = maximum(max(abs(xs[i][P.iN] - solz(t)[P.iN]) / max(solz(t)[P.iN], 50),
                                     abs(xs[i][P.iC] - solz(t)[P.iC]) / max(solz(t)[P.iC], 50)) for (i, t) in enumerate(ts))
                @test maxrel < 1e-5
            end
        end

        @testset "scenario $scen: RK4 transition vs adaptive Tsit5" begin
            for _ in 1:5
                x = random_state(P, rng); u = rand(rng) * 1.47
                xr = F(x, u, P, Δ)
                sol = solve(ODEProblem((du, y, p, t) -> plant_rhs!(du, y, u, P), x, (0.0, Δ)), Tsit5(); reltol=1e-11, abstol=1e-11)
                @test maximum(abs.(xr .- sol.u[end]) ./ max.(abs.(sol.u[end]), 1.0)) < 1e-6
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
                @test maximum(abs.(hcat(A, B) .- Jfd) ./ scale) < 1e-5
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
            x0 = x0_for(P, 300.0, 300.0); U = [0.2, 1.0, 0.5, 1.2]
            xs = [copy(x0)]; As = Matrix{Float64}[]; Bs = Vector{Float64}[]
            for k in eachindex(U)
                xn, A, B = linearize(xs[k], U[k], P, Δ); push!(xs, xn); push!(As, A); push!(Bs, B)
            end
            g = adjoint_gradient(xs, As, Bs, U, P, W)
            gfd = [(Up = copy(U); Um = copy(U); Up[k] += 1e-6; Um[k] -= 1e-6;
                    (total_cost(x0, Up, P, Δ, W) - total_cost(x0, Um, P, Δ, W)) / 2e-6) for k in eachindex(U)]
            @test maximum(abs.(g .- gfd) ./ max.(abs.(gfd), 1e-6)) < 1e-4
        end

        @testset "scenario $scen: discrete search equals brute force" begin
            x0 = x0_for(P, 300.0, 300.0); levels = [0.0, 0.62, 1.0, 1.47]; H = 3
            sol = discrete_solve(x0, P, Δ, H, levels, W)
            brute = minimum(total_cost(x0, [a, b, c], P, Δ, W) for a in levels, b in levels, c in levels)
            @test isapprox(sol.J, brute; rtol=1e-10)
        end

        @testset "scenario $scen: SQP respects bounds and does not increase the cost" begin
            x0 = x0_for(P, 300.0, 300.0); U0 = fill(0.5, 5)
            sol = sqp_solve(x0, U0, P, Δ, W)
            @test all(0 .<= sol.U .<= W.umax + 1e-12)
            @test sol.J <= total_cost(x0, U0, P, Δ, W) + 1e-9
            cl = closed_loop(P, P, x0, Δ, 4, 4, W; mode=:continuous)
            @test all(0 .<= cl.U .<= W.umax + 1e-12)
        end
    end

    @testset "long runs: 30 days, 14-day look-ahead re-planned every 2 days" begin
        P = scenario_plant(:A); W = default_weights(P); levels = Float64.(CFG["dose"]["levels_data"])
        days = Float64(CFG["schedule"]["long_days"]); n = round(Int, days / Δ); H = round(Int, CFG["schedule"]["horizon_days"] / Δ)
        @test n == 15 && H == 7
        x0 = x0_for(P, 300.0, 300.0)
        cl = closed_loop(P, P, x0, Δ, n, H, W; mode=:discrete, levels=levels)
        @test length(cl.U) == n                                # one dose per 2-day interval up to day 30
        @test all(u -> any(isapprox(u, l) for l in levels), cl.U)   # only measured doses
        o = outcomes(P, x0, cl.U, Δ, W)
        @test isapprox(o.ts[end], days; atol=1e-9)              # the trajectory reaches day 30
        @test all(isfinite, o.N) && all(isfinite, o.C) && all(>=(0), o.N) && all(>=(0), o.C)
        @test isapprox(o.cum_dose, sum(cl.U) * Δ)
    end

    @testset "relaxed barrier" begin
        μ, δ = 0.05, 0.02
        @test relaxed_barrier(0.5, μ, δ) ≈ -μ * log(0.5)
        @test isapprox(relaxed_barrier(δ + 1e-9, μ, δ), relaxed_barrier(δ - 1e-9, μ, δ); atol=1e-6)     # continuous
        d(z) = (relaxed_barrier(z + 1e-7, μ, δ) - relaxed_barrier(z - 1e-7, μ, δ)) / 2e-7
        @test isapprox(d(δ + 1e-5), d(δ - 1e-5); rtol=1e-2)                                             # C¹ at δ
        @test relaxed_barrier(-0.5, μ, δ) > relaxed_barrier(0.0, μ, δ) > relaxed_barrier(δ, μ, δ)         # steep, finite outside
        @test isfinite(relaxed_barrier(-1.0, μ, δ))
    end
end
