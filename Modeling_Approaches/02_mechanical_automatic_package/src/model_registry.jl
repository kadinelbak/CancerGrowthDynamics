module ModelRegistry

using GrowthParameterEstimation

export ensure_model_registry!, local_model_specs

const _registry_loaded = Ref(false)

# ---------------------------------------------------------------------------
# ODE definitions: 5-arg (du, u, p, t, exposure) where exposure(t) = dose
# ---------------------------------------------------------------------------

function _one_stage_transit_linear_kill!(du, u, p, t, exposure)
    r, K, k_kill, tau = p
    N1 = max(u[1], 0.0)
    N2 = max(u[2], 0.0)
    C = max(exposure(t), 0.0)
    growth = r * N1 * max(0.0, 1 - N1 / max(K, 1e-8))
    damage_flux = k_kill * C * N1
    clearance_flux = (1 / max(tau, 1e-8)) * N2
    du[1] = growth - damage_flux
    du[2] = damage_flux - clearance_flux
end

function _hill_death_instantaneous!(du, u, p, t, exposure)
    r, K, emax, ec50, hill_n = p
    N = max(u[1], 0.0)
    C = max(exposure(t), 0.0)
    kill = emax * (C^hill_n / (ec50^hill_n + C^hill_n + 1e-12))
    growth = r * N * max(0.0, 1 - N / max(K, 1e-8))
    du[1] = growth - kill * N
end

function _transit_hill_combined!(du, u, p, t, exposure)
    r, K, emax, ec50, hill_n, tau = p
    N1 = max(u[1], 0.0)
    N2 = max(u[2], 0.0)
    C = max(exposure(t), 0.0)
    k_hill = emax * (C^hill_n / (ec50^hill_n + C^hill_n + 1e-12))
    growth = r * N1 * max(0.0, 1 - N1 / max(K, 1e-8))
    damage_flux = k_hill * N1
    clearance_flux = (1 / max(tau, 1e-8)) * N2
    du[1] = growth - damage_flux
    du[2] = damage_flux - clearance_flux
end

# ---------------------------------------------------------------------------
# Local model spec: stores everything needed for fitting without using the
# package's (non-existent) register_model / ModelSpec API.
# ---------------------------------------------------------------------------

struct LocalModelSpec
    name::String
    ode!::Function                          # 5-arg: (du, u, p, t, exposure)
    param_names::Vector{Symbol}
    bounds::Vector{Tuple{Float64,Float64}}
    n_states::Int
    observable::Function                    # (u::Vector) -> scalar
    p0_factory::Function                    # (r0, K0, dose) -> Vector{Float64}
end

function local_model_specs()::Vector{LocalModelSpec}
    return [
        LocalModelSpec(
            "hill_death_instantaneous",
            _hill_death_instantaneous!,
            [:r, :K, :emax, :ec50, :hill_n],
            [(1e-6, 5.0), (1e-3, 1e7), (0.0, 1.0), (1e-3, 10.0), (0.1, 5.0)],
            1,
            u -> u[1],
            (r0, K0, dose) -> [r0, K0, 0.3, max(dose, 0.1), 1.0],
        ),
        LocalModelSpec(
            "one_stage_transit_linear_kill",
            _one_stage_transit_linear_kill!,
            [:r, :K, :k_kill, :tau],
            [(1e-6, 5.0), (1e-3, 1e7), (0.0, 10.0), (1e-3, 20.0)],
            2,
            u -> u[1] + u[2],
            (r0, K0, dose) -> [r0, K0, 0.5, 5.0],
        ),
        LocalModelSpec(
            "transit_hill_combined",
            _transit_hill_combined!,
            [:r, :K, :emax, :ec50, :hill_n, :tau],
            [(1e-6, 5.0), (1e-3, 1e7), (0.0, 1.0), (1e-3, 10.0), (0.1, 5.0), (1e-3, 20.0)],
            2,
            u -> u[1] + u[2],
            (r0, K0, dose) -> [r0, K0, 0.3, max(dose, 0.1), 1.0, 5.0],
        ),
    ]
end

function ensure_model_registry!()
    # No-op: specs are now stored locally in local_model_specs().
    # The GrowthParameterEstimation built-ins (logistic_growth, gompertz_growth)
    # are used directly for untreated conditions via run_single_fit.
    _registry_loaded[] = true
    return nothing
end

end
