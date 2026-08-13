module Simulation

using DifferentialEquations
using DataFrames

using ..Exposure
using ..Registry

export SimulationResult, SweepGrid, SweepResult, simulate, classify_failure, run_sweep

struct SimulationResult{T}
    success::Bool
    reason::String
    times::Vector{Float64}
    states::Matrix{T}
    observed::Vector{T}
end

struct SweepGrid
    seed_totals::Vector{Float64}
    resistant_fractions::Vector{Float64}
    doses::Vector{Float64}
    times::Vector{Float64}
end

struct SweepResult
    summary::DataFrame
    simulations::Dict{NamedTuple,SimulationResult}
end

function classify_failure(error_obj)
    msg = lowercase(string(error_obj))
    if occursin("dt", msg) || occursin("unstable", msg) || occursin("instability", msg)
        return "solver_failure"
    elseif occursin("domain", msg) || occursin("nan", msg)
        return "biological_domain_failure"
    else
        return "unknown_failure"
    end
end

function _solver_from_symbol(kind::Symbol)
    if kind == :ode
        return Tsit5()
    elseif kind == :stiff_ode
        return Rodas5()
    else
        return Tsit5()
    end
end

function _to_state_matrix(sol)
    n_t = length(sol.u)
    n_s = length(sol.u[1])
    arr = Matrix{eltype(sol.u[1])}(undef, n_s, n_t)
    for i in 1:n_t
        for j in 1:n_s
            arr[j, i] = sol.u[i][j]
        end
    end
    return arr
end

function simulate(
    spec::Registry.ModelSpec,
    times::AbstractVector{<:Real},
    params::AbstractVector{<:Real};
    u0::AbstractVector{<:Real},
    exposure::Exposure.AbstractExposure = Exposure.ConstantExposure(0.0),
    reltol::Float64 = 1e-8,
    abstol::Float64 = 1e-8,
    enforce_nonnegative::Bool = true,
)
    t_obs = Float64.(collect(times))
    tspan = (minimum(t_obs), maximum(t_obs))
    param_vec = collect(params)
    value_type = promote_type(Float64, eltype(param_vec), eltype(u0))

    wrapped! = function (du, u, p, t)
        spec.dynamics!(du, u, p, t, exposure)
        if enforce_nonnegative
            for i in eachindex(du)
                if !isfinite(du[i])
                    du[i] = 0.0
                end
            end
        end
        return nothing
    end

    prob = ODEProblem(wrapped!, value_type.(u0), tspan, value_type.(param_vec))

    try
        solver = _solver_from_symbol(spec.solver_type)
        sol = solve(prob, solver; saveat=t_obs, reltol=reltol, abstol=abstol, isoutofdomain=(u, p, t) -> any(x -> x < -1e-10 || !isfinite(x), u))

        if sol.retcode != ReturnCode.Success
            return SimulationResult(false, "solver_failure: $(sol.retcode)", t_obs, zeros(Float64, length(u0), length(t_obs)), fill(NaN, length(t_obs)))
        end

        state_matrix = _to_state_matrix(sol)
        if enforce_nonnegative
            state_matrix .= max.(state_matrix, 0.0)
        end

        observed = [spec.observation(state_matrix[:, i], params, t_obs[i]) for i in eachindex(t_obs)]
        return SimulationResult(true, "ok", t_obs, state_matrix, observed)
    catch err
        return SimulationResult(false, classify_failure(err), t_obs, zeros(Float64, length(u0), length(t_obs)), fill(NaN, length(t_obs)))
    end
end

function _make_sweep_u0(spec::Registry.ModelSpec, seed_total::Real, resistant_fraction::Real)
    n_states = length(spec.state_names)
    u0 = zeros(Float64, n_states)
    total = max(Float64(seed_total), 0.0)
    frac_r = clamp(Float64(resistant_fraction), 0.0, 1.0)

    if n_states == 1
        u0[1] = total
    elseif n_states >= 2
        u0[1] = total * (1 - frac_r)
        u0[2] = total * frac_r
    end

    return u0
end

function run_sweep(
    spec::Registry.ModelSpec,
    params::AbstractVector{<:Real},
    grid::SweepGrid;
    exposure_builder::Function = dose -> Exposure.ConstantExposure(Float64(dose)),
    reltol::Float64 = 1e-8,
    abstol::Float64 = 1e-8,
    enforce_nonnegative::Bool = true,
)
    rows = NamedTuple[]
    sims = Dict{NamedTuple,SimulationResult}()

    for seed_total in grid.seed_totals, resistant_fraction in grid.resistant_fractions, dose in grid.doses
        u0 = _make_sweep_u0(spec, seed_total, resistant_fraction)
        exposure = exposure_builder(dose)
        sim = simulate(
            spec,
            grid.times,
            params;
            u0=u0,
            exposure=exposure,
            reltol=reltol,
            abstol=abstol,
            enforce_nonnegative=enforce_nonnegative,
        )

        key = (seed_total=Float64(seed_total), resistant_fraction=Float64(resistant_fraction), dose=Float64(dose))
        sims[key] = sim

        final_total = sim.success ? sum(sim.states[:, end]) : NaN
        final_resistant_fraction = sim.success && size(sim.states, 1) >= 2 ? sim.states[2, end] / max(sum(sim.states[:, end]), 1e-12) : NaN
        auc_total = sim.success ? sum(sim.observed) : NaN

        push!(rows, (
            seed_total=Float64(seed_total),
            resistant_fraction=Float64(resistant_fraction),
            dose=Float64(dose),
            success=sim.success,
            reason=sim.reason,
            final_total=final_total,
            final_observed=sim.success ? sim.observed[end] : NaN,
            final_resistant_fraction=final_resistant_fraction,
            auc_total=auc_total,
        ))
    end

    return SweepResult(DataFrame(rows), sims)
end

end # module Simulation
