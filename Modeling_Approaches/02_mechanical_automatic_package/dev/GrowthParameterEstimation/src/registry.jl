module Registry

using DifferentialEquations

using ..Models

export ModelSpec, register_model, get_model, list_models, clear_registry!, register_builtin_models!

struct ModelSpec
    name::String
    dynamics!::Function
    param_names::Vector{Symbol}
    bounds::Vector{Tuple{Float64,Float64}}
    state_names::Vector{Symbol}
    observation::Function
    solver_type::Symbol
    metadata::Dict{Symbol,Any}
end

const MODEL_REGISTRY = Dict{String,ModelSpec}()

function register_model(spec::ModelSpec; overwrite::Bool=false)
    if haskey(MODEL_REGISTRY, spec.name) && !overwrite
        error("Model $(spec.name) already registered. Use overwrite=true to replace it.")
    end
    MODEL_REGISTRY[spec.name] = spec
    return spec
end

get_model(name::String) = MODEL_REGISTRY[name]
list_models() = sort(collect(keys(MODEL_REGISTRY)))
clear_registry!() = empty!(MODEL_REGISTRY)

function _ode_adapter(model!::Function)
    return function (du, u, p, t, exposure)
        model!(du, u, p, t)
        return nothing
    end
end

function _theta_hill_inhibition!(du, u, p, t, exposure)
    r, K, theta, ic50, hill = p
    N = max(u[1], 0.0)
    drug = max(exposure(t), 0.0)
    inhibition = drug^hill / (ic50^hill + drug^hill + 1e-12)
    growth = r * N * (1 - (N / max(K, 1e-8))^theta)
    du[1] = max(1e-12, 1 - inhibition) * growth
end

function _theta_hill_kill!(du, u, p, t, exposure)
    r, K, theta, emax_kill, ic50, hill = p
    N = max(u[1], 0.0)
    drug = max(exposure(t), 0.0)
    kill = emax_kill * (drug^hill / (ic50^hill + drug^hill + 1e-12))
    growth = r * N * (1 - (N / max(K, 1e-8))^theta)
    du[1] = growth - kill * N
end

function _sensitive_resistant!(du, u, p, t, exposure)
    rS, rR, K, kSR, emax, ic50, hill = p
    S = max(u[1], 0.0)
    R = max(u[2], 0.0)
    N = max(S + R, 1e-12)
    drug = max(exposure(t), 0.0)
    killS = emax * (drug^hill / (ic50^hill + drug^hill + 1e-12))
    common = max(0.0, 1 - N / max(K, 1e-8))
    du[1] = rS * S * common - killS * S - kSR * S
    du[2] = rR * R * common + kSR * S
end

function _null_coculture!(du, u, p, t, exposure)
    rS, KS, rR, KR = p
    S = max(u[1], 0.0)
    R = max(u[2], 0.0)
    du[1] = rS * S * max(0.0, 1 - S / max(KS, 1e-8))
    du[2] = rR * R * max(0.0, 1 - R / max(KR, 1e-8))
end

function _lotka_volterra_competition!(du, u, p, t, exposure)
    rS, KS, alpha_SR, rR, KR, alpha_RS = p
    S = max(u[1], 0.0)
    R = max(u[2], 0.0)
    du[1] = rS * S * max(0.0, 1 - (S + alpha_SR * R) / max(KS, 1e-8))
    du[2] = rR * R * max(0.0, 1 - (R + alpha_RS * S) / max(KR, 1e-8))
end

function _lotka_volterra_hill_competition!(du, u, p, t, exposure)
    rS, KS, alpha_SR, rR, KR, alpha_RS, emaxS, ic50S, emaxR, ic50R, hill = p
    S = max(u[1], 0.0)
    R = max(u[2], 0.0)
    drug = max(exposure(t), 0.0)
    effectS = emaxS * (drug^hill / (ic50S^hill + drug^hill + 1e-12))
    effectR = emaxR * (drug^hill / (ic50R^hill + drug^hill + 1e-12))
    growthS = rS * S * max(0.0, 1 - (S + alpha_SR * R) / max(KS, 1e-8))
    growthR = rR * R * max(0.0, 1 - (R + alpha_RS * S) / max(KR, 1e-8))
    du[1] = growthS - effectS * S
    du[2] = growthR - effectR * R
end

function _damage_repair_arrest!(du, u, p, t, exposure)
    r, K, k_damage, k_repair, k_death, ic50, hill = p
    S = max(u[1], 0.0)
    D = max(u[2], 0.0)
    N = max(S + D, 1e-12)
    drug = max(exposure(t), 0.0)
    effect = drug^hill / (ic50^hill + drug^hill + 1e-12)
    growth = r * S * max(0.0, 1 - N / max(K, 1e-8))
    damage_flux = k_damage * effect * S
    repair_flux = k_repair * D
    du[1] = growth - damage_flux + repair_flux
    du[2] = damage_flux - repair_flux - k_death * D
end

function _adaptive_ic50!(du, u, p, t, exposure)
    r, K, emax, ic50_0, hill, k_adapt = p
    N = max(u[1], 0.0)
    A = max(u[2], 0.0)
    ic50_t = ic50_0 * (1 + A)
    drug = max(exposure(t), 0.0)
    inhibition = emax * (drug^hill / (ic50_t^hill + drug^hill + 1e-12))
    growth = r * N * max(0.0, 1 - N / max(K, 1e-8))
    du[1] = growth * max(1e-12, 1 - inhibition)
    du[2] = k_adapt * inhibition - 0.05 * A
end

function _pkpd_inhibition!(du, u, p, t, exposure)
    r, K, emax, ic50, hill, k_elim, k_in = p
    N = max(u[1], 0.0)
    C = max(u[2], 0.0)
    inhibition = emax * (C^hill / (ic50^hill + C^hill + 1e-12))
    du[1] = r * N * max(0.0, 1 - N / max(K, 1e-8)) * max(1e-12, 1 - inhibition)
    du[2] = k_in * max(exposure(t), 0.0) - k_elim * C
end

function _bi_exponential!(du, u, p, t, exposure)
    a, b = p
    du[1] = -a * u[1]
    du[2] = -b * u[2]
end

function _transit_chain_erlang!(du, u, p, t, exposure)
    r, K, ktr, emax, ic50, hill = p
    N = max(u[1], 0.0)
    X1 = max(u[2], 0.0)
    X2 = max(u[3], 0.0)
    X3 = max(u[4], 0.0)
    drug = max(exposure(t), 0.0)
    effect = emax * (drug^hill / (ic50^hill + drug^hill + 1e-12))
    growth = r * N * max(0.0, 1 - N / max(K, 1e-8))
    du[1] = growth - ktr * X3
    du[2] = effect * N - ktr * X1
    du[3] = ktr * (X1 - X2)
    du[4] = ktr * (X2 - X3)
end

function register_builtin_models!()
    clear_registry!()

    register_model(ModelSpec(
        "logistic_growth",
        _ode_adapter(Models.logistic_growth!),
        [:r, :K],
        [(1e-6, 5.0), (1e-3, 1e7)],
        [:N],
        (u, p, t) -> u[1],
        :ode,
        Dict(:family => :baseline),
    ))

    register_model(ModelSpec(
        "gompertz_growth",
        _ode_adapter(Models.gompertz_growth!),
        [:a, :b, :K],
        [(1e-6, 5.0), (1e-6, 10.0), (1e-3, 1e7)],
        [:N],
        (u, p, t) -> u[1],
        :ode,
        Dict(:family => :baseline),
    ))

    register_model(ModelSpec(
        "theta_logistic_hill_inhibition",
        _theta_hill_inhibition!,
        [:r, :K, :theta, :ic50, :hill],
        [(1e-6, 5.0), (1e-3, 1e7), (0.1, 5.0), (1e-8, 1e4), (0.1, 8.0)],
        [:N],
        (u, p, t) -> u[1],
        :ode,
        Dict(:family => :baseline_hill),
    ))

    register_model(ModelSpec(
        "theta_logistic_hill_kill",
        _theta_hill_kill!,
        [:r, :K, :theta, :emax_kill, :ic50, :hill],
        [(1e-6, 5.0), (1e-3, 1e7), (0.1, 5.0), (0.0, 20.0), (1e-8, 1e4), (0.1, 8.0)],
        [:N],
        (u, p, t) -> u[1],
        :ode,
        Dict(:family => :baseline_hill),
    ))

    register_model(ModelSpec(
        "null_coculture",
        _null_coculture!,
        [:rS, :KS, :rR, :KR],
        [(1e-6, 5.0), (1e-3, 1e7), (1e-6, 5.0), (1e-3, 1e7)],
        [:S, :R],
        (u, p, t) -> u[1] + u[2],
        :ode,
        Dict(:family => :coculture_null),
    ))

    register_model(ModelSpec(
        "lotka_volterra_competition",
        _lotka_volterra_competition!,
        [:rS, :KS, :alpha_SR, :rR, :KR, :alpha_RS],
        [(1e-6, 5.0), (1e-3, 1e7), (0.0, 5.0), (1e-6, 5.0), (1e-3, 1e7), (0.0, 5.0)],
        [:S, :R],
        (u, p, t) -> u[1] + u[2],
        :ode,
        Dict(:family => :coculture_competition),
    ))

    register_model(ModelSpec(
        "lotka_volterra_hill_competition",
        _lotka_volterra_hill_competition!,
        [:rS, :KS, :alpha_SR, :rR, :KR, :alpha_RS, :emaxS, :ic50S, :emaxR, :ic50R, :hill],
        [(1e-6, 5.0), (1e-3, 1e7), (0.0, 5.0), (1e-6, 5.0), (1e-3, 1e7), (0.0, 5.0), (0.0, 20.0), (1e-8, 1e4), (0.0, 20.0), (1e-8, 1e4), (0.1, 8.0)],
        [:S, :R],
        (u, p, t) -> u[1] + u[2],
        :ode,
        Dict(:family => :coculture_competition_hill),
    ))

    register_model(ModelSpec(
        "sensitive_resistant",
        _sensitive_resistant!,
        [:rS, :rR, :K, :kSR, :emax, :ic50, :hill],
        [(1e-6, 5.0), (1e-6, 5.0), (1e-3, 1e7), (0.0, 2.0), (0.0, 20.0), (1e-8, 1e4), (0.1, 8.0)],
        [:S, :R],
        (u, p, t) -> u[1] + u[2],
        :ode,
        Dict(:family => :mechanistic),
    ))

    register_model(ModelSpec(
        "damage_repair_arrest",
        _damage_repair_arrest!,
        [:r, :K, :k_damage, :k_repair, :k_death, :ic50, :hill],
        [(1e-6, 5.0), (1e-3, 1e7), (0.0, 10.0), (0.0, 10.0), (0.0, 10.0), (1e-8, 1e4), (0.1, 8.0)],
        [:S, :D],
        (u, p, t) -> u[1] + u[2],
        :ode,
        Dict(:family => :mechanistic),
    ))

    register_model(ModelSpec(
        "adaptive_ic50",
        _adaptive_ic50!,
        [:r, :K, :emax, :ic50_0, :hill, :k_adapt],
        [(1e-6, 5.0), (1e-3, 1e7), (0.0, 10.0), (1e-8, 1e4), (0.1, 8.0), (0.0, 5.0)],
        [:N, :A],
        (u, p, t) -> u[1],
        :ode,
        Dict(:family => :adaptive),
    ))

    register_model(ModelSpec(
        "pkpd_inhibition",
        _pkpd_inhibition!,
        [:r, :K, :emax, :ic50, :hill, :k_elim, :k_in],
        [(1e-6, 5.0), (1e-3, 1e7), (0.0, 10.0), (1e-8, 1e4), (0.1, 8.0), (1e-6, 10.0), (0.0, 20.0)],
        [:N, :C],
        (u, p, t) -> u[1],
        :ode,
        Dict(:family => :pkpd),
    ))

    register_model(ModelSpec(
        "transit_chain_erlang",
        _transit_chain_erlang!,
        [:r, :K, :ktr, :emax, :ic50, :hill],
        [(1e-6, 5.0), (1e-3, 1e7), (1e-6, 10.0), (0.0, 10.0), (1e-8, 1e4), (0.1, 8.0)],
        [:N, :X1, :X2, :X3],
        (u, p, t) -> u[1],
        :ode,
        Dict(:family => :delay_surrogate),
    ))

    register_model(ModelSpec(
        "bi_exponential_response",
        _bi_exponential!,
        [:a, :b],
        [(1e-6, 10.0), (1e-6, 10.0)],
        [:N1, :N2],
        (u, p, t) -> u[1] + u[2],
        :ode,
        Dict(:family => :phenomenological),
    ))

    return nothing
end

end # module Registry
