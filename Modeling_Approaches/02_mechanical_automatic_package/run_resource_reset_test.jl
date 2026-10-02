# Low-resource test (user spec §3): do resource-depletion models with a reset at every media change
# (every 2 days) explain the 0.5% FBS monocultures better than the carried constant-effect models?
#   B1  exogenous: R = e^(−λτ), τ = time since the last media change (dR/dt = −λR, reset to 1)
#   B2  mechanistic: dR/dt = −λ·X·R/1000, reset to 1 at each media change (solver events)
# Each is fitted with GrowthParameterEstimation.run_joint_fit (media-change events are passed through
# a solver wrapper) and compared with the Stage 1b carried model by BIC. Writes
# outputs/control/resource_reset_test.csv (read by the optimal-control report).
#   julia -t 4 --project=. run_resource_reset_test.jl
include(joinpath(@__DIR__, "zoo_common.jl"))
using OrdinaryDiffEq: SciMLBase
using OrdinaryDiffEq: DiscreteCallback
const GPE = GrowthParameterEstimation

# solver wrapper that adds media-change resets (callback + tstops) inside the package fitter
struct WithEvents{A,C} <: SciMLBase.AbstractODEAlgorithm
    alg::A; cb::C; tstops::Vector{Float64}
end
SciMLBase.__solve(prob::SciMLBase.AbstractODEProblem, a::WithEvents, args...; kwargs...) =
    solve(prob, a.alg, args...; callback=a.cb, tstops=a.tstops, kwargs...)

P = CSV.read(joinpath(@__DIR__, "outputs", "zoo_fits", "zoo_parameters.csv"), DataFrame)
sel = CSV.read(joinpath(@__DIR__, "outputs", "zoo_fits", "zoo_model_selection.csv"), DataFrame)
out_rows = NamedTuple[]
pv(g, n) = P[(P.group .== g) .& (P.param .== n), :estimate][1]
growth = Dict("A2780" => (r=pv("s1-a2780", "r_N"), K=pv("s1-a2780", "K_N"), θ=1.0),
              "A2780cis" => (r=pv("s1-a2780cis", "r_C"), K=pv("s1-a2780cis", "K_C"), θ=pv("s1-a2780cis", "θ_C")))
cur(g) = (r = sel[(sel.group .== g) .& (sel.optimal .== true), :][1, :]; (Float64(r.bic), Int(r.k), String(r.candidate)))
current = Dict("A2780" => cur("s1b-a2780"), "A2780cis" => cur("s1b-a2780cis"))

for L in LINES
    m = lrmono[L.line]; g = growth[L.line]
    x, y = m.x[.!m.bad], m.y[.!m.bad]; n = length(x)
    shift = L.line in SHIFT_FROM_DAY1 ? 1.0 : 0.0
    changes = [c - shift for c in 2.0:2.0:30.0 if c - shift > 0]          # media changes every 2 days (original frame)
    last_change(t) = (i = searchsortedlast(changes, t); i == 0 ? -shift : changes[i])
    G(X, reff) = reff * X * (1 - (max(X, 0) / g.K)^g.θ)
    ds = [(x=x, y=y, state_index=1)]
    cb = DiscreteCallback((u, t, integ) -> t in changes, integ -> (integ.u[2] = one(eltype(integ.u))); save_positions=(false, false))
    ev = WithEvents(Tsit5(), cb, changes)
    cands = [
        ("B1: growth × R, R = e^{-λτ} (exogenous)", 1, [0.3], [(1e-4, 20.0)], [x[1]], nothing,
         (du, u, p, t) -> (R = exp(-p[1] * (t - last_change(t))); du[1] = G(u[1], g.r * R))),
        ("B1: starvation death d(1-R), R = e^{-λτ}", 2, [0.3, 0.5], [(1e-4, 20.0), (0.0, 20.0)], [x[1]], nothing,
         (du, u, p, t) -> (R = exp(-p[1] * (t - last_change(t))); du[1] = G(u[1], g.r) - p[2] * (1 - R) * max(u[1], 0))),
        ("B1: growth × R and death d(1-R)", 2, [0.3, 0.5], [(1e-4, 20.0), (0.0, 20.0)], [x[1]], nothing,
         (du, u, p, t) -> (R = exp(-p[1] * (t - last_change(t))); du[1] = G(u[1], g.r * R) - p[2] * (1 - R) * max(u[1], 0))),
        ("B2: growth × R, dR/dt = -λ X R/1000, reset", 1, [0.3], [(1e-5, 50.0)], [y[1], 1.0], ev,
         (du, u, p, t) -> (X = max(u[1], 0); R = max(u[2], 0); du[1] = G(X, g.r * R); du[2] = -p[1] * X * R / 1000)),
        ("B2: death d(1-R), dR/dt = -λ X R/1000, reset", 2, [0.3, 0.5], [(1e-5, 50.0), (0.0, 20.0)], [y[1], 1.0], ev,
         (du, u, p, t) -> (X = max(u[1], 0); R = max(u[2], 0); du[1] = G(X, g.r) - p[2] * (1 - R) * X; du[2] = -p[1] * X * R / 1000)),
        ("B2: growth × R and death d(1-R), reset", 2, [0.3, 0.5], [(1e-5, 50.0), (0.0, 20.0)], [y[1], 1.0], ev,
         (du, u, p, t) -> (X = max(u[1], 0); R = max(u[2], 0); du[1] = G(X, g.r * R) - p[2] * (1 - R) * X; du[2] = -p[1] * X * R / 1000)),
    ]
    println("\n== $(L.line) low resource: n = $n, current pick BIC $(current[L.line][1]) (k=$(current[L.line][2])); media changes at $(changes[1:4])…"); flush(stdout)
    for (name, k, p0, bnds, u0, solver, model) in cands
        u0v = length(u0) == 1 ? [y[1]] : u0
        best = (bic=Inf, params=Float64[])
        for opt in (:nelder_mead,)
            st = GPE.generate_multistarts(p0, bnds; n_starts=4, seed=3)
            for s in st
                f = try
                    GPE.run_joint_fit(model, ds, u0v, s; bounds=bnds, optimizer=opt, maxiters=1500, reltol=1e-8, abstol=1e-8,
                                      initial_time=0.0, (solver === nothing ? (;) : (; solver=solver))...)
                catch err
                    nothing
                end
                f !== nothing && isfinite(f.bic) && f.bic < best.bic && (best = (bic=f.bic, params=f.params))
            end
        end
        @printf("  %-48s k=%d BIC=%8.2f  ΔBIC vs current %+7.2f  params %s\n", name, k, best.bic, best.bic - current[L.line][1], round.(best.params; sigdigits=4)); flush(stdout)
        push!(out_rows, (line=L.line, model=name, k=k, bic=best.bic, delta_bic=best.bic - current[L.line][1], params=join(round.(best.params; sigdigits=4), "; "),
                         carried=current[L.line][3], carried_bic=current[L.line][1], carried_k=current[L.line][2]))
    end
end

mkpath(joinpath(@__DIR__, "outputs", "control"))
CSV.write(joinpath(@__DIR__, "outputs", "control", "resource_reset_test.csv"), DataFrame(out_rows))
println("wrote outputs/control/resource_reset_test.csv")
