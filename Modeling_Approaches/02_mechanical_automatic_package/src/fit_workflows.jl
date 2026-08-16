module FitWorkflows

using CSV
using DataFrames
using Statistics
using OrdinaryDiffEq
using BlackBoxOptim
using GrowthParameterEstimation
using Plots

using ..IOUtils
using ..ModelRegistry

export load_untreated_baseline, run_condition_fit!

const DENSITY_LOG_CONTRAST_BOUND = log(1.05)
const PRIMARY_UNTREATED_MODELS = Set(["logistic_growth", "gompertz_growth", "theta_logistic_growth"])

# ---------------------------------------------------------------------------
# Untreated baseline loading
# ---------------------------------------------------------------------------

function load_untreated_baseline(; start::AbstractString = pwd())
    root = IOUtils.package_root(start)
    base_dir = joinpath(root, "outputs", "csv", "monoculture_untreated")
    summary_path = joinpath(base_dir, "best_model_summary.csv")
    params_path = joinpath(base_dir, "params", "best_params.csv")

    if !isfile(summary_path) || !isfile(params_path)
        return nothing
    end

    summary_df = CSV.read(summary_path, DataFrame)
    params_df  = CSV.read(params_path,  DataFrame)
    summary_map = Dict(string(row.key) => row.value for row in eachrow(summary_df))

    model_name = get(summary_map, "best_model", missing)
    params     = Float64.(params_df.value)

    (ismissing(model_name) || isempty(params)) && return nothing

    return (
        model       = String(model_name),
        params      = params,
        summary_path = summary_path,
        params_path  = params_path,
    )
end

# ---------------------------------------------------------------------------
# Write staged baseline CSV (for traceability)
# ---------------------------------------------------------------------------

function _write_staged_baseline(out, untreated_baseline)
    untreated_baseline === nothing && return nothing

    rows = Tuple{String,String}[("baseline_model", untreated_baseline.model)]
    append!(rows, [("param_$(i)", string(untreated_baseline.params[i])) for i in eachindex(untreated_baseline.params)])

    staged_path = joinpath(out.csv, "staged_untreated_baseline.csv")
    CSV.write(staged_path, DataFrame(key = first.(rows), value = last.(rows)))
    return staged_path
end

_nonmissing_strings(vals) = String[v for v in vals if v !== nothing]
_safe_string(v; default::AbstractString = "") = (v === missing || v === nothing) ? String(default) : string(v)

function _load_untreated_monoculture_cellline_baselines(; start::AbstractString = pwd())
    root = IOUtils.find_repo_root(start)
    auto_path = joinpath(IOUtils.package_root(start), "outputs", "csv", "monoculture_untreated", "untreated_group_baselines.csv")
    if isfile(auto_path)
        df_auto = CSV.read(auto_path, DataFrame)
        required_auto = Set([:cell_line, :density, :r, :K])
        if all(c -> c in propertynames(df_auto), required_auto)
            by_cell_density = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
            spec_by_cell_density = Dict{Tuple{String,String},NamedTuple}()
            by_cell = Dict{String,Tuple{Float64,Float64}}()

            for row in eachrow(df_auto)
                cell = String(row.cell_line)
                density = String(row.density)
                r = Float64(row.r)
                K = Float64(row.K)
                by_cell_density[(cell, density)] = (r, K)
                spec_by_cell_density[(cell, density)] = (
                    model = :best_model in propertynames(df_auto) ? String(row.best_model) : "logistic_growth",
                    pooling_mode = :pooling_mode in propertynames(df_auto) ? String(row.pooling_mode) : "legacy",
                    r = r,
                    K = K,
                    shape_parameter = :shape_parameter in propertynames(df_auto) ? _safe_string(row.shape_parameter) : "",
                    shape_value = :shape_value in propertynames(df_auto) && !ismissing(row.shape_value) ? Float64(row.shape_value) : NaN,
                    inheritance_allowed = :inheritance_allowed in propertynames(df_auto) ? Bool(row.inheritance_allowed) : true,
                )
            end

            cell_groups = groupby(df_auto, :cell_line)
            for g in cell_groups
                cell = String(first(g.cell_line))
                by_cell[cell] = (mean(Float64.(g.r)), mean(Float64.(g.K)))
            end

            return (
                path = auto_path,
                by_cell_density = by_cell_density,
                spec_by_cell_density = spec_by_cell_density,
                by_cell = by_cell,
                density_aware = true,
            )
        end
    end

    candidates = [
        joinpath(root, "Modeling_Approaches", "01_mechanical_manual", "Untreated MonoCulture", "outputs", "csv", "untreated_monoculture_logistic_parameters.csv"),
        joinpath(root, "Modeling_Approaches", "01_mechanical_manual", "Untreated MonoCulture", "untreated_monoculture_logistic_parameters.csv"),
    ]
    existing = filter(isfile, candidates)
    isempty(existing) && return nothing
    path = first(existing)

    df = CSV.read(path, DataFrame)
    required = Set([:CellLine, :Density, :r, :K])
    all(c -> c in propertynames(df), required) || return nothing

    by_cell_density = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
    by_cell = Dict{String,Tuple{Float64,Float64}}()

    for row in eachrow(df)
        cell = String(row.CellLine)
        density = String(row.Density)
        r = Float64(row.r)
        K = Float64(row.K)
        by_cell_density[(cell, density)] = (r, K)
    end

    cell_groups = groupby(df, :CellLine)
    for g in cell_groups
        cell = String(first(g.CellLine))
        by_cell[cell] = (mean(Float64.(g.r)), mean(Float64.(g.K)))
    end

    spec_by_cell_density = Dict(
        key => (
            model = "logistic_growth",
            pooling_mode = "legacy",
            r = values[1],
            K = values[2],
            shape_parameter = "",
            shape_value = NaN,
            inheritance_allowed = true,
        ) for (key, values) in by_cell_density
    )
    return (
        path = path,
        by_cell_density = by_cell_density,
        spec_by_cell_density = spec_by_cell_density,
        by_cell = by_cell,
        density_aware = false,
    )
end

function _baseline_rk_for_group(
    baseline_map,
    cell_line::Union{Missing,AbstractString},
    density::Union{Missing,AbstractString},
    fallback_r::Float64,
    fallback_K::Float64,
)
    baseline_map === nothing && return (fallback_r, fallback_K)
    ismissing(cell_line) && return (fallback_r, fallback_K)

    cell = String(cell_line)
    den = ismissing(density) ? "" : String(density)

    if haskey(baseline_map.by_cell_density, (cell, den))
        return baseline_map.by_cell_density[(cell, den)]
    end
    if haskey(baseline_map.by_cell, cell)
        return baseline_map.by_cell[cell]
    end
    return (fallback_r, fallback_K)
end

function _fit_anchored_linear_kill(
    x::Vector{Float64},
    y::Vector{Float64},
    dose::Float64,
    r0::Float64,
    K0::Float64;
    solver = Rodas5(),
    max_time::Float64 = 12.0,
)
    function ode!(du, u, p, t)
        k_kill = p[1]
        N = max(u[1], 0.0)
        growth = r0 * N * max(0.0, 1 - N / max(K0, 1e-8))
        kill = k_kill * max(dose, 0.0) * N
        du[1] = growth - kill
    end

    u0 = [max(y[1], 1.0)]
    prob = ODEProblem(ode!, u0, (x[1], x[end]), [0.05])

    function loss(p_vec)
        pv = Vector{Float64}(p_vec)
        try
            sol = solve(remake(prob; p = pv), solver;
                        reltol = 1e-6, abstol = 1e-6,
                        saveat = x, maxiters = 50_000)
            sol.retcode == ReturnCode.Success || return 1e12
            yhat = first.(sol.u)
            return sum((y .- yhat) .^ 2)
        catch
            return 1e12
        end
    end

    result = bboptimize(
        loss;
        SearchRange = [(0.0, 8.0)],
        NumDimensions = 1,
        Method = :de_rand_1_bin,
        MaxTime = max_time,
        TraceMode = :silent,
    )
    p_opt = Vector{Float64}(result.archive_output.best_candidate)

    bic, ssr = try
        sol2 = solve(remake(prob; p = p_opt), solver;
                     reltol = 1e-10, abstol = 1e-10, saveat = x)
        yhat2 = first.(sol2.u)
        ssr2 = sum((y .- yhat2) .^ 2)
        n, k = length(x), length(p_opt)
        n * log(max(ssr2, 1e-20) / n) + k * log(n), ssr2
    catch
        1e12, 1e12
    end

    return (params = p_opt, bic = bic, ssr = ssr)
end

function _fit_anchored_hill_kill(
    x::Vector{Float64},
    y::Vector{Float64},
    dose::Float64,
    r0::Float64,
    K0::Float64;
    solver = Rodas5(),
    max_time::Float64 = 12.0,
)
    function ode!(du, u, p, t)
        emax, ec50, hill_n = p
        N = max(u[1], 0.0)
        C = max(dose, 0.0)
        kill = emax * (C^hill_n / (ec50^hill_n + C^hill_n + 1e-12))
        growth = r0 * N * max(0.0, 1 - N / max(K0, 1e-8))
        du[1] = growth - kill * N
    end

    u0 = [max(y[1], 1.0)]
    prob = ODEProblem(ode!, u0, (x[1], x[end]), [0.3, max(dose, 0.1), 1.0])

    function loss(p_vec)
        pv = Vector{Float64}(p_vec)
        try
            sol = solve(remake(prob; p = pv), solver;
                        reltol = 1e-6, abstol = 1e-6,
                        saveat = x, maxiters = 50_000)
            sol.retcode == ReturnCode.Success || return 1e12
            yhat = first.(sol.u)
            return sum((y .- yhat) .^ 2)
        catch
            return 1e12
        end
    end

    result = bboptimize(
        loss;
        SearchRange = [(0.0, 1.5), (1e-3, 10.0), (0.1, 6.0)],
        NumDimensions = 3,
        Method = :de_rand_1_bin,
        MaxTime = max_time,
        TraceMode = :silent,
    )
    p_opt = Vector{Float64}(result.archive_output.best_candidate)

    bic, ssr = try
        sol2 = solve(remake(prob; p = p_opt), solver;
                     reltol = 1e-10, abstol = 1e-10, saveat = x)
        yhat2 = first.(sol2.u)
        ssr2 = sum((y .- yhat2) .^ 2)
        n, k = length(x), length(p_opt)
        n * log(max(ssr2, 1e-20) / n) + k * log(n), ssr2
    catch
        1e12, 1e12
    end

    return (params = p_opt, bic = bic, ssr = ssr)
end

function _treated_dose_metadata(dose::Float64)
    if isapprox(dose, 0.67; atol = 0.08)
        return (effect_level = 0.25, ic_label = "IC25")
    elseif isapprox(dose, 1.0; atol = 0.08)
        return (effect_level = 0.50, ic_label = "IC50")
    elseif isapprox(dose, 1.47; atol = 0.08)
        return (effect_level = 0.75, ic_label = "IC75")
    end
    return (effect_level = dose, ic_label = "concentration")
end

function _hill_effect(signal, emax, ec50, hill_n)
    s = max(signal, zero(signal))
    h = max(hill_n, oftype(hill_n, 1e-8))
    midpoint = max(ec50, oftype(ec50, 1e-8))
    return emax * s^h / (midpoint^h + s^h + oftype(s, 1e-12))
end

function _ramp_activation(t, t0, lambda, onset = 0.0)
    elapsed = max(t - t0 - onset, zero(t))
    return 1 - exp(-max(lambda, oftype(lambda, 1e-8)) * elapsed)
end

function _joint_treated_model_specs(
    doses::Vector{Float64},
    effect_levels::Vector{Float64},
    density_indices::Vector{Int},
    baseline_specs::Vector{NamedTuple},
    t0::Float64,
    pooling_mode::String,
)
    length(doses) == length(baseline_specs) == length(density_indices) || error("Treated trajectory metadata is inconsistent")

    function untreated_growth(N, trajectory_index)
        baseline = baseline_specs[trajectory_index]
        r = baseline.r
        K = max(baseline.K, 1e-8)
        if baseline.model == "gompertz_growth"
            positive_N = max(N, oftype(N, 1e-8))
            return r * positive_N * log(K / positive_N)
        elseif baseline.model == "theta_logistic_growth"
            theta = max(baseline.shape_value, 1e-8)
            return r * N * max(0.0, 1 - (N / K)^theta)
        elseif baseline.model == "logistic_simple_death"
            return r * N * max(0.0, 1 - N / K) - max(baseline.shape_value, 0.0) * N
        elseif baseline.model == "allee_growth"
            threshold = max(baseline.shape_value, 1e-8)
            return r * N * max(0.0, 1 - N / K) * ((N / threshold) - 1)
        end
        return r * N * max(0.0, 1 - N / K)
    end

    function pooled_spec(base_p0, base_bounds, base_names, amplitude_indices)
        nbase = length(base_p0)
        if pooling_mode == "shared"
            p0 = copy(base_p0)
            bounds = copy(base_bounds)
            names = copy(base_names)
        elseif pooling_mode == "partial_5pct"
            p0 = vcat(base_p0, [0.0])
            bounds = vcat(base_bounds, [(-DENSITY_LOG_CONTRAST_BOUND, DENSITY_LOG_CONTRAST_BOUND)])
            names = vcat(base_names, [:log_contrast_effect])
        elseif pooling_mode == "independent_diagnostic"
            p0 = vcat(base_p0, base_p0)
            bounds = vcat(base_bounds, base_bounds)
            names = vcat(Symbol.(string.(base_names) .* "_20k"), Symbol.(string.(base_names) .* "_30k"))
        else
            error("Unsupported treated pooling mode: $(pooling_mode)")
        end

        local_parameters = function (p, trajectory_index)
            density_index = density_indices[trajectory_index]
            values = if pooling_mode == "independent_diagnostic"
                first_index = (density_index - 1) * nbase + 1
                collect(p[first_index:(first_index + nbase - 1)])
            else
                collect(p[1:nbase])
            end
            if pooling_mode == "partial_5pct"
                sign = density_index == 1 ? -1.0 : 1.0
                multiplier = exp(sign * p[nbase + 1])
                for parameter_index in amplitude_indices
                    values[parameter_index] *= multiplier
                end
            end
            return values
        end
        return (
            p0 = p0,
            bounds = bounds,
            param_names = names,
            base_param_names = copy(base_names),
            base_bounds = copy(base_bounds),
            amplitude_indices = copy(amplitude_indices),
            local_parameters = local_parameters,
            nbase = nbase,
        )
    end

    linear = pooled_spec([0.1], [(0.0, 8.0)], [:k_kill], [1])
    hill = pooled_spec([0.3, max(median(doses), 0.1), 1.0], [(0.0, 3.0), (1e-3, 10.0), (0.1, 6.0)], [:emax, :ec50_uM, :hill_n], [1])
    decay = pooled_spec([0.2, 0.1], [(0.0, 4.0), (0.0, 5.0)], [:k_kill, :lambda], [1])
    logistic_hill = pooled_spec([1.5, 0.5, 3.0], [(0.0, 4.0), (0.05, 1.0), (0.2, 10.0)], [:emax, :ec50_effect, :hill_n], [1])
    ramp = pooled_spec([1.8, 0.5, 4.0, 0.6], [(0.0, 5.0), (0.05, 1.0), (0.2, 12.0), (0.01, 5.0)], [:emax, :ec50_effect, :hill_n, :lambda], [1])
    ramp_onset = pooled_spec([1.8, 0.5, 4.0, 0.8, 2.0], [(0.0, 5.0), (0.05, 1.0), (0.2, 12.0), (0.01, 5.0), (0.0, 7.0)], [:emax, :ec50_effect, :hill_n, :lambda, :t_onset], [1])
    transit = pooled_spec([1.8, 0.8, 1.5, 0.5], [(0.0, 5.0), (0.01, 5.0), (0.0, 7.0), (0.01, 5.0)], [:emax, :lambda, :t_onset, :k_clear], [1])
    populations = pooled_spec([2.0, 0.5, 0.8, 1.5, 0.1], [(0.0, 5.0), (0.0, 3.0), (0.01, 5.0), (0.0, 7.0), (0.0, 0.95)], [:emax_sensitive, :emax_tolerant, :lambda, :t_onset, :f_tolerant0], [1, 2])
    platinum_pkpd = pooled_spec(
        [1.5, 0.5, 3.0, 1.0, 0.5],
        [(0.0, 5.0), (0.01, 5.0), (0.2, 10.0), (0.02, 5.0), (0.02, 5.0)],
        [:emax, :ec50_bound, :hill_n, :k_efflux, :k_repair],
        [1],
    )

    function linear_kill!(du, u, p, t)
        for i in eachindex(doses)
            k_kill = linear.local_parameters(p, i)[1]
            N = max(u[i], 0.0)
            growth = untreated_growth(N, i)
            du[i] = growth - k_kill * max(doses[i], 0.0) * N
        end
    end

    function hill_kill!(du, u, p, t)
        for i in eachindex(doses)
            emax, ec50, hill_n = hill.local_parameters(p, i)
            h = max(hill_n, 1e-8)
            N = max(u[i], 0.0)
            C = max(doses[i], 0.0)
            effect = emax * (C^h / (max(ec50, 1e-8)^h + C^h + 1e-12))
            growth = untreated_growth(N, i)
            du[i] = growth - effect * N
        end
    end

    function time_decay_dose_scaled!(du, u, p, t)
        elapsed = max(t - t0, 0.0)
        for i in eachindex(doses)
            k_kill, lambda = decay.local_parameters(p, i)
            N = max(u[i], 0.0)
            growth = untreated_growth(N, i)
            kill = k_kill * max(doses[i], 0.0) * exp(-max(lambda, 0.0) * elapsed)
            du[i] = growth - kill * N
        end
    end

    function ic_effect_logistic_hill!(du, u, p, t)
        for i in eachindex(effect_levels)
            emax, ec50, hill_n = logistic_hill.local_parameters(p, i)
            N = max(u[i], zero(u[i]))
            kill = _hill_effect(effect_levels[i], emax, ec50, hill_n)
            growth = untreated_growth(N, i)
            du[i] = growth - kill * N
        end
    end

    function ic_effect_hill_ramp!(du, u, p, t)
        for i in eachindex(effect_levels)
            emax, ec50, hill_n, lambda = ramp.local_parameters(p, i)
            activation = _ramp_activation(t, t0, lambda)
            N = max(u[i], zero(u[i]))
            kill = activation * _hill_effect(effect_levels[i], emax, ec50, hill_n)
            growth = untreated_growth(N, i)
            du[i] = growth - kill * N
        end
    end

    function ic_effect_hill_ramp_onset!(du, u, p, t)
        for i in eachindex(effect_levels)
            emax, ec50, hill_n, lambda, onset = ramp_onset.local_parameters(p, i)
            activation = _ramp_activation(t, t0, lambda, onset)
            N = max(u[i], zero(u[i]))
            kill = activation * _hill_effect(effect_levels[i], emax, ec50, hill_n)
            growth = untreated_growth(N, i)
            du[i] = growth - kill * N
        end
    end

    function ic_effect_transit_death!(du, u, p, t)
        ec50 = 0.5
        hill_n = 4.0
        for i in eachindex(effect_levels)
            emax, lambda, onset, k_clear = transit.local_parameters(p, i)
            activation = _ramp_activation(t, t0, lambda, onset)
            live_idx = 2i - 1
            damaged_idx = 2i
            live = max(u[live_idx], zero(u[live_idx]))
            damaged = max(u[damaged_idx], zero(u[damaged_idx]))
            damage_rate = activation * _hill_effect(effect_levels[i], emax, ec50, hill_n)
            growth = untreated_growth(live, i)
            damage_flux = damage_rate * live
            du[live_idx] = growth - damage_flux
            du[damaged_idx] = damage_flux - k_clear * damaged
        end
    end

    function ic_effect_two_population!(du, u, p, t)
        ec50 = 0.5
        hill_n = 4.0
        for i in eachindex(effect_levels)
            emax_sensitive, emax_tolerant, lambda, onset, _ = populations.local_parameters(p, i)
            activation = _ramp_activation(t, t0, lambda, onset)
            sensitive_idx = 2i - 1
            tolerant_idx = 2i
            sensitive = max(u[sensitive_idx], zero(u[sensitive_idx]))
            tolerant = max(u[tolerant_idx], zero(u[tolerant_idx]))
            total = sensitive + tolerant
            signal = effect_levels[i]
            kill_sensitive = activation * _hill_effect(signal, emax_sensitive, ec50, hill_n)
            kill_tolerant = activation * _hill_effect(signal, emax_tolerant, ec50, hill_n)
            baseline = baseline_specs[i]
            total_growth = untreated_growth(total, i)
            sensitive_share = total > 0 ? sensitive / total : zero(total)
            tolerant_share = total > 0 ? tolerant / total : zero(total)
            du[sensitive_idx] = total_growth * sensitive_share - kill_sensitive * sensitive
            du[tolerant_idx] = total_growth * tolerant_share - kill_tolerant * tolerant
        end
    end

    # Lumped adaptation of El-Kareh and Secomb (2003). Uptake and DNA binding
    # scales are fixed to one because cell-count trajectories alone cannot
    # separately identify both scales and the downstream kill amplitude.
    function intracellular_platinum_pkpd!(du, u, p, t)
        for i in eachindex(doses)
            emax, ec50_bound, hill_n, k_efflux, k_repair = platinum_pkpd.local_parameters(p, i)
            live_idx = 3i - 2
            intracellular_idx = live_idx + 1
            bound_idx = live_idx + 2
            live = max(u[live_idx], zero(u[live_idx]))
            intracellular = max(u[intracellular_idx], zero(u[intracellular_idx]))
            dna_bound = max(u[bound_idx], zero(u[bound_idx]))
            extracellular = max(doses[i], 0.0)
            kill = _hill_effect(dna_bound, emax, ec50_bound, hill_n)

            du[live_idx] = untreated_growth(live, i) - kill * live
            du[intracellular_idx] = extracellular - max(k_efflux, 1e-8) * intracellular
            du[bound_idx] = intracellular - max(k_repair, 1e-8) * dna_bound
        end
    end

    return Dict(
        "joint_anchored_linear_kill" => (
            model = linear_kill!,
            p0 = linear.p0,
            bounds = linear.bounds,
            layout = :single,
            param_names = linear.param_names,
            base_param_names = linear.base_param_names,
            base_bounds = linear.base_bounds,
            amplitude_indices = linear.amplitude_indices,
            local_parameters = linear.local_parameters,
            dose_basis = "concentration_uM",
            behavior = "immediate net growth or extinction",
        ),
        "joint_anchored_hill_kill" => (
            model = hill_kill!,
            p0 = hill.p0,
            bounds = hill.bounds,
            layout = :single,
            param_names = hill.param_names,
            base_param_names = hill.base_param_names,
            base_bounds = hill.base_bounds,
            amplitude_indices = hill.amplitude_indices,
            local_parameters = hill.local_parameters,
            dose_basis = "concentration_uM",
            behavior = "immediate logistic growth or extinction",
        ),
        "joint_time_decay_dose_scaled" => (
            model = time_decay_dose_scaled!,
            p0 = decay.p0,
            bounds = decay.bounds,
            layout = :single,
            param_names = decay.param_names,
            base_param_names = decay.base_param_names,
            base_bounds = decay.base_bounds,
            amplitude_indices = decay.amplitude_indices,
            local_parameters = decay.local_parameters,
            dose_basis = "concentration_uM",
            behavior = "waning drug effect and regrowth",
        ),
        "joint_ic_effect_logistic_hill" => (
            model = ic_effect_logistic_hill!,
            p0 = logistic_hill.p0,
            bounds = logistic_hill.bounds,
            layout = :single,
            param_names = logistic_hill.param_names,
            base_param_names = logistic_hill.base_param_names,
            base_bounds = logistic_hill.base_bounds,
            amplitude_indices = logistic_hill.amplitude_indices,
            local_parameters = logistic_hill.local_parameters,
            dose_basis = "IC_effect_level",
            behavior = "logistic survival below threshold and extinction above it",
        ),
        "joint_ic_effect_hill_ramp" => (
            model = ic_effect_hill_ramp!,
            p0 = ramp.p0,
            bounds = ramp.bounds,
            layout = :single,
            param_names = ramp.param_names,
            base_param_names = ramp.base_param_names,
            base_bounds = ramp.base_bounds,
            amplitude_indices = ramp.amplitude_indices,
            local_parameters = ramp.local_parameters,
            dose_basis = "IC_effect_level",
            behavior = "early logistic growth followed by accumulating kill",
        ),
        "joint_ic_effect_hill_ramp_onset" => (
            model = ic_effect_hill_ramp_onset!,
            p0 = ramp_onset.p0,
            bounds = ramp_onset.bounds,
            layout = :single,
            param_names = ramp_onset.param_names,
            base_param_names = ramp_onset.base_param_names,
            base_bounds = ramp_onset.base_bounds,
            amplitude_indices = ramp_onset.amplitude_indices,
            local_parameters = ramp_onset.local_parameters,
            dose_basis = "IC_effect_level",
            behavior = "delayed switch from logistic growth to extinction",
        ),
        "joint_ic_effect_transit_death" => (
            model = ic_effect_transit_death!,
            p0 = transit.p0,
            bounds = transit.bounds,
            layout = :transit,
            param_names = transit.param_names,
            base_param_names = transit.base_param_names,
            base_bounds = transit.base_bounds,
            amplitude_indices = transit.amplitude_indices,
            local_parameters = transit.local_parameters,
            dose_basis = "IC_effect_level",
            behavior = "delayed damage transit; EC50=0.5, Hill n=4, visible damaged fraction=0.5",
        ),
        "joint_ic_effect_two_population" => (
            model = ic_effect_two_population!,
            p0 = populations.p0,
            bounds = populations.bounds,
            layout = :two_population,
            param_names = populations.param_names,
            base_param_names = populations.base_param_names,
            base_bounds = populations.base_bounds,
            amplitude_indices = populations.amplitude_indices,
            local_parameters = populations.local_parameters,
            dose_basis = "IC_effect_level",
            behavior = "sensitive/tolerant populations; EC50=0.5 and Hill n=4 fixed",
        ),
        "joint_intracellular_platinum_pkpd" => (
            model = intracellular_platinum_pkpd!,
            p0 = platinum_pkpd.p0,
            bounds = platinum_pkpd.bounds,
            layout = :platinum_pkpd,
            param_names = platinum_pkpd.param_names,
            base_param_names = platinum_pkpd.base_param_names,
            base_bounds = platinum_pkpd.base_bounds,
            amplitude_indices = platinum_pkpd.amplitude_indices,
            local_parameters = platinum_pkpd.local_parameters,
            dose_basis = "extracellular_cisplatin_uM",
            behavior = "lumped cisplatin uptake, DNA binding/repair, and DNA-bound Hill kill",
        ),
    )
end

function _joint_model_inputs(spec, datasets, initial_counts, density_indices)
    if spec.layout == :single
        mapped = [(
            x = datasets[i].x,
            y = datasets[i].y,
            state_index = i,
            residual_scale = datasets[i].residual_scale,
        ) for i in eachindex(datasets)]
        return mapped, copy(initial_counts), nothing
    elseif spec.layout == :transit
        mapped = [(
            x = datasets[i].x,
            y = datasets[i].y,
            residual_scale = datasets[i].residual_scale,
            observable = let live_idx = 2i - 1, damaged_idx = 2i
                (u, p, t) -> u[live_idx] + 0.5 * u[damaged_idx]
            end,
        ) for i in eachindex(datasets)]
        u0 = reduce(vcat, ([n0, 0.0] for n0 in initial_counts))
        return mapped, u0, nothing
    elseif spec.layout == :two_population
        mapped = [(
            x = datasets[i].x,
            y = datasets[i].y,
            residual_scale = datasets[i].residual_scale,
            observable = let sensitive_idx = 2i - 1, tolerant_idx = 2i
                (u, p, t) -> u[sensitive_idx] + u[tolerant_idx]
            end,
        ) for i in eachindex(datasets)]
        u0 = reduce(vcat, ([0.9 * n0, 0.1 * n0] for n0 in initial_counts))
        u0_builder = p -> begin
            reduce(vcat, ([begin
                local_p = spec.local_parameters(p, i)
                f_tolerant = clamp(local_p[5], zero(local_p[5]), one(local_p[5]))
                [n0 * (1 - f_tolerant), n0 * f_tolerant]
            end for (i, n0) in enumerate(initial_counts)]))
        end
        return mapped, u0, u0_builder
    elseif spec.layout == :platinum_pkpd
        mapped = [(
            x = datasets[i].x,
            y = datasets[i].y,
            residual_scale = datasets[i].residual_scale,
            observable = let live_idx = 3i - 2
                (u, p, t) -> u[live_idx]
            end,
        ) for i in eachindex(datasets)]
        u0 = reduce(vcat, ([n0, 0.0, 0.0] for n0 in initial_counts))
        return mapped, u0, nothing
    end
    error("Unknown joint model layout: $(spec.layout)")
end

function _joint_ic75_long_term_fate(model_name::String, p, r_anchor::Float64)
    signal = 0.75
    net_rates = if model_name in (
        "joint_ic_effect_logistic_hill",
        "joint_ic_effect_hill_ramp",
        "joint_ic_effect_hill_ramp_onset",
    )
        [r_anchor - _hill_effect(signal, p[1], p[2], p[3])]
    elseif model_name == "joint_ic_effect_transit_death"
        [r_anchor - _hill_effect(signal, p[1], 0.5, 4.0)]
    elseif model_name == "joint_ic_effect_two_population"
        [
            r_anchor - _hill_effect(signal, p[1], 0.5, 4.0),
            r_anchor - _hill_effect(signal, p[2], 0.5, 4.0),
        ]
    elseif model_name == "joint_intracellular_platinum_pkpd"
        emax, ec50_bound, hill_n, k_efflux, k_repair = p[1:5]
        steady_bound = 1.47 / max(k_efflux * k_repair, 1e-12)
        [r_anchor - _hill_effect(steady_bound, emax, ec50_bound, hill_n)]
    elseif model_name == "joint_time_decay_dose_scaled"
        [r_anchor]
    else
        Float64[]
    end

    isempty(net_rates) && return (net_rate = NaN, fate = "not_assessed")
    maximum(net_rates) < 0 && return (net_rate = maximum(net_rates), fate = "eventual_extinction")
    return (net_rate = maximum(net_rates), fate = "survival_or_regrowth")
end

function _write_treated_trajectory_summary(decoded::DataFrame, out, time_col::Symbol, area_col::Symbol)
    rows = NamedTuple[]
    for grp in groupby(decoded, [:cell_line, :density, :dose])
        means = combine(groupby(grp, time_col), area_col => mean => :observed_mean)
        sort!(means, time_col)
        nrow(means) >= 2 || continue
        peak_idx = argmax(means.observed_mean)
        first_count = Float64(first(means.observed_mean))
        peak_count = Float64(means.observed_mean[peak_idx])
        last_count = Float64(last(means.observed_mean))
        last_peak_ratio = last_count / max(peak_count, 1e-12)
        behavior = if last_peak_ratio <= 0.15 && peak_idx < nrow(means)
            "extinction_like"
        elseif last_peak_ratio <= 0.70 && peak_idx < nrow(means)
            "post_peak_decline"
        else
            "logistic_or_plateau"
        end
        dose = Float64(first(grp.dose))
        metadata = _treated_dose_metadata(dose)
        push!(rows, (
            cell_line = _safe_string(first(grp.cell_line); default = "pooled"),
            density = _safe_string(first(grp.density); default = "pooled"),
            concentration_uM = dose,
            ic_label = metadata.ic_label,
            effect_level = metadata.effect_level,
            first_count = first_count,
            peak_time = Float64(means[peak_idx, time_col]),
            peak_count = peak_count,
            last_time = Float64(last(means[!, time_col])),
            last_count = last_count,
            last_to_peak = last_peak_ratio,
            observed_behavior = behavior,
        ))
    end
    summary = isempty(rows) ? DataFrame() : DataFrame(rows)
    CSV.write(joinpath(out.csv, "monoculture_treated_trajectory_behavior.csv"), summary)
    return summary
end

function _run_joint_treated_monoculture_fitting(
    decoded::DataFrame,
    out;
    start::AbstractString = pwd(),
    untreated_baseline = nothing,
    max_time_per_fit::Float64 = 45.0,
    render_plots::Bool = true,
)
    time_col = :time
    :day in propertynames(decoded) && !(:time in propertynames(decoded)) && (time_col = :day)
    area_col = :count in propertynames(decoded) ? :count : error("Joint treated fitting requires a count column")

    colset = Set(Symbol.(names(decoded)))
    (:cell_line in colset && :density in colset && :dose in colset) || return DataFrame()
    _write_treated_trajectory_summary(decoded, out, time_col, area_col)

    r0 = isnothing(untreated_baseline) ? 0.3 : untreated_baseline.params[1]
    K0 = isnothing(untreated_baseline) ? 3000.0 : untreated_baseline.params[2]
    baseline_map = _load_untreated_monoculture_cellline_baselines(; start = start)

    rows = NamedTuple[]
    parameter_rows = NamedTuple[]
    initial_condition_rows = NamedTuple[]
    overlay_parts = DataFrame[]
    gcols = [:cell_line]
    for grp in groupby(decoded, gcols)
        cell_label = _safe_string(first(grp.cell_line); default = "pooled")
        r_anchor, K_anchor = _baseline_rk_for_group(baseline_map, cell_label, missing, r0, K0)

        trajectory_groups = collect(groupby(grp, [:density, :dose]))
        length(trajectory_groups) < 4 && continue
        sort!(trajectory_groups, by = g -> (
            _safe_string(first(g.density); default = "pooled"),
            Float64(first(g.dose)),
        ))

        datasets = NamedTuple[]
        doses = Float64[]
        effect_levels = Float64[]
        ic_labels = String[]
        u0 = Float64[]
        density_labels = String[]
        for dgrp in trajectory_groups
            density_label = _safe_string(first(dgrp.density); default = "pooled")
            dose_val = Float64(first(dgrp.dose))
            means = combine(groupby(dgrp, time_col), area_col => mean => :y_mean)
            sort!(means, time_col)
            nrow(means) >= 3 || continue
            push!(doses, dose_val)
            push!(density_labels, density_label)
            metadata = _treated_dose_metadata(dose_val)
            push!(effect_levels, metadata.effect_level)
            push!(ic_labels, metadata.ic_label)
            first_observed_value = max(Float64(first(means.y_mean)), 1.0)
            fixed_u0 = _fixed_day0_total(density_label)
            residual_scale = max(maximum(Float64.(means.y_mean)), 1.0)
            push!(u0, fixed_u0)
            push!(datasets, (
                x = Float64.(means[!, time_col]),
                y = Float64.(means.y_mean),
                state_index = length(datasets) + 1,
                residual_scale = residual_scale,
            ))
            push!(initial_condition_rows, (
                cell_line = cell_label,
                density = density_label,
                dose = dose_val,
                ic_label = metadata.ic_label,
                first_observed_value = first_observed_value,
                fixed_u0 = fixed_u0,
                u0_time_day = 0.0,
                u0_over_K = fixed_u0 / max(K_anchor, 1e-12),
                residual_scale = residual_scale,
                u0_strategy = "fixed_nominal_day0_67_for_20k_100_for_30k",
            ))
        end

        length(datasets) >= 4 || continue
        all_times = sort(unique(vcat([ds.x for ds in datasets]...)))
        specs = _joint_treated_model_specs(doses, effect_levels, r_anchor, K_anchor, 0.0)
        fitted_densities = sort(unique(density_labels))
        fitted_doses = sort(unique(doses))
        den_label = join(fitted_densities, "/")
        joint_label = "densities=$(join(fitted_densities, "/"));doses=$(join(fitted_doses, "/"))"

        for model_name in sort(collect(keys(specs)))
            spec = specs[model_name]
            println("  Joint fitting $(model_name), cell=$(cell_label), density=$(den_label), doses=$(joint_label)...")
            model_datasets, model_u0, u0_builder = _joint_model_inputs(spec, datasets, u0)
            model_maxiters = spec.layout == :single ?
                max(40, Int(round(max_time_per_fit * 50))) :
                max(30, Int(round(max_time_per_fit * 35)))
            fit = try
                GrowthParameterEstimation.run_joint_fit(
                    spec.model,
                    model_datasets,
                    model_u0,
                    spec.p0;
                    solver = Tsit5(),
                    bounds = spec.bounds,
                    u0_builder = u0_builder,
                    maxiters = model_maxiters,
                    reltol = 1e-7,
                    abstol = 1e-7,
                    optimizer = :nelder_mead,
                )
            catch e
                @warn "Joint treated monoculture fit failed" model_name cell_label den_label exception=e
                nothing
            end
            fit === nothing && continue
            isfinite(Float64(fit.bic)) && abs(Float64(fit.sse)) < 9.99e11 || continue
            fate = _joint_ic75_long_term_fate(model_name, fit.params, r_anchor)

            for (parameter_index, parameter_name) in enumerate(spec.param_names)
                parameter_index <= length(fit.params) || continue
                lower_bound, upper_bound = spec.bounds[parameter_index]
                value = Float64(fit.params[parameter_index])
                bound_span = upper_bound - lower_bound
                push!(parameter_rows, (
                    model = model_name,
                    cell_line = cell_label,
                    density = den_label,
                    environment = "$(cell_label) | $(den_label)",
                    parameter = String(parameter_name),
                    value = value,
                    lower_bound = Float64(lower_bound),
                    upper_bound = Float64(upper_bound),
                    bound_position = bound_span > 0 ? (value - lower_bound) / bound_span : NaN,
                    r_anchor = r_anchor,
                    K_anchor = K_anchor,
                    parameter_scope = "shared_across_densities_and_doses",
                ))
            end


            push!(rows, (
                model = model_name,
                dose = NaN,
                bic = Float64(fit.bic),
                ssr = Float64(fit.raw_sse),
                scaled_ssr = Float64(fit.scaled_sse),
                params = string((names = spec.param_names, fit = fit.params, r_anchor = r_anchor, K_anchor = K_anchor, joint_doses_uM = fitted_doses, ic_effect_levels = sort(unique(effect_levels)), density_coupling = "shared_kinetics_fixed_nominal_day0_u0", residual_scaling = "trajectory_peak", package_api = "GrowthParameterEstimation.run_joint_fit")),
                cell_line = cell_label,
                density = den_label,
                mix = "",
                joint_group = joint_label,
                n_doses = length(fitted_doses),
                n_densities = length(fitted_densities),
                n_points = sum(length(ds.x) for ds in datasets),
                dose_basis = spec.dose_basis,
                model_behavior = spec.behavior,
                ic_mapping = join(unique(["$(doses[i])=$(ic_labels[i])" for i in eachindex(doses)]), ";"),
                initial_condition_strategy = "fixed nominal day-zero total: 67 for 20k, 100 for 30k",
                parameter_scope = "shared across 20k and 30k within cell line",
                residual_scaling = "divide each trajectory residual by its observed peak",
                ic75_low_density_net_rate = fate.net_rate,
                ic75_long_term_fate = fate.fate,
            ))

            for (ds_idx, ds) in enumerate(datasets)
                push!(overlay_parts, DataFrame(
                    time = ds.x,
                    observed = ds.y,
                    predicted = fit.predictions[ds_idx],
                    model = fill(model_name, length(ds.x)),
                    cell_line = fill(cell_label, length(ds.x)),
                    density = fill(density_labels[ds_idx], length(ds.x)),
                    dose = fill(doses[ds_idx], length(ds.x)),
                    ic_label = fill(ic_labels[ds_idx], length(ds.x)),
                    effect_level = fill(effect_levels[ds_idx], length(ds.x)),
                    joint_group = fill(joint_label, length(ds.x)),
                ))
            end
        end
    end

    joint_df = isempty(rows) ? DataFrame() : sort!(DataFrame(rows), :bic)
    joint_path = joinpath(out.csv, "monoculture_treated_joint_dose_model_ranking.csv")
    CSV.write(joint_path, joint_df)
    CSV.write(
        joinpath(out.csv, "monoculture_treated_joint_initial_condition_diagnostics.csv"),
        isempty(initial_condition_rows) ? DataFrame() : DataFrame(initial_condition_rows),
    )

    if !isempty(joint_df)
        bic_summary = GrowthParameterEstimation.summarize_joint_bic(joint_df)
        CSV.write(joinpath(out.csv, "monoculture_treated_joint_bic_long.csv"), bic_summary.long)
        CSV.write(joinpath(out.csv, "monoculture_treated_joint_bic_matrix.csv"), bic_summary.matrix)
        CSV.write(joinpath(out.csv, "monoculture_treated_joint_bic_summary.csv"), bic_summary.aggregate)
        CSV.write(joinpath(out.csv, "monoculture_treated_joint_environment_winners.csv"), bic_summary.winners)
        cell_line_top5 = GrowthParameterEstimation.summarize_joint_bic_by_group(
            joint_df;
            group_col = :cell_line,
            top_n = 5,
            environment_cols = [:density],
        )
        if :environments_fit in propertynames(cell_line_top5)
            density_counts = Dict(
                String(cell) => length(unique(String.(decoded[String.(decoded.cell_line) .== String(cell), :density])))
                for cell in unique(String.(decoded.cell_line))
            )
            cell_line_top5.environments_fit = [get(density_counts, String(cell), 0) for cell in cell_line_top5.cell_line]
        end
        CSV.write(joinpath(out.csv, "monoculture_treated_joint_cell_line_top5.csv"), cell_line_top5)

        parameter_df = DataFrame(parameter_rows)
        CSV.write(joinpath(out.csv, "monoculture_treated_joint_parameter_estimates.csv"), parameter_df)
        parameter_stability = GrowthParameterEstimation.summarize_joint_parameter_stability(parameter_df)
        CSV.write(joinpath(out.csv, "monoculture_treated_joint_parameter_stability.csv"), parameter_stability)

        model_stability = combine(
            groupby(parameter_stability, :model),
            nrow => :parameters_assessed,
            :stability_class => (values -> count(==("tight"), values)) => :tight_parameters,
            :stability_class => (values -> count(==("moderate"), values)) => :moderate_parameters,
            :stability_class => (values -> count(==("broad"), values)) => :broad_parameters,
            :near_optimization_bound => count => :parameters_near_bound,
            :bound_range_fraction => maximum => :maximum_bound_range_fraction,
            :coefficient_of_variation => (values -> begin
                finite_values = filter(isfinite, values)
                isempty(finite_values) ? NaN : median(finite_values)
            end) => :median_parameter_cv,
        )
        model_selection = leftjoin(bic_summary.aggregate, model_stability; on = :model)
        sort!(model_selection, [:delta_total_bic, :mean_environment_rank])
        CSV.write(joinpath(out.csv, "monoculture_treated_joint_model_selection_summary.csv"), model_selection)
    end
    if !isempty(overlay_parts)
        fig_dir = joinpath(out.csv, "figures")
        mkpath(fig_dir)
        overlay_df = vcat(overlay_parts...)
        CSV.write(joinpath(fig_dir, "monoculture_treated_joint_dose_overlays.csv"), overlay_df)
        render_plots || return joint_df

        @eval using Plots

        image_fig_dir = joinpath(out.images, "figures")
        mkpath(image_fig_dir)
        model_labels = Dict(
            "joint_anchored_linear_kill" => "Anchored linear kill",
            "joint_anchored_hill_kill" => "Anchored Hill kill",
            "joint_time_decay_dose_scaled" => "Time-decay kill",
            "joint_ic_effect_logistic_hill" => "IC-effect logistic Hill",
            "joint_ic_effect_hill_ramp" => "IC-effect Hill ramp",
            "joint_ic_effect_hill_ramp_onset" => "Delayed Hill ramp",
            "joint_ic_effect_transit_death" => "Transit damage/death",
            "joint_ic_effect_two_population" => "Sensitive/tolerant",
            "joint_intracellular_platinum_pkpd" => "Intracellular platinum PK/PD",
        )
        for grp in groupby(overlay_df, [:cell_line, :density])
            cell_label = _safe_string(first(grp.cell_line); default = "pooled")
            den_label = _safe_string(first(grp.density); default = "pooled")
            model_order = sort(unique(String.(grp.model)))
            color_values = palette(:tab10)
            model_colors = Dict(model => color_values[mod1(i, length(color_values))] for (i, model) in enumerate(model_order))
            dose_groups = sort(collect(groupby(grp, :dose)); by = dgrp -> first(dgrp.effect_level), rev = true)
            panels = Any[]
            for (dose_index, dgrp) in enumerate(dose_groups)
                dose_label = "$(first(dgrp.ic_label)) ($(first(dgrp.dose)) uM)"
                show_legend = dose_index == length(dose_groups)
                panel = plot(
                    xlabel = "Time (day)",
                    ylabel = dose_index == 1 ? "Cell count" : "",
                    title = dose_label,
                    legend = false,
                    legendfontsize = 7,
                    titlefontsize = 10,
                )
                obs = unique(dgrp[:, [:time, :observed]])
                sort!(obs, :time)
                scatter!(panel, obs.time, obs.observed;
                    label = "",
                    color = :black,
                    ms = 3.5,
                    markerstrokewidth = 0,
                )
                for model in model_order
                    mgrp = dgrp[String.(dgrp.model) .== model, :]
                    isempty(mgrp) && continue
                    pred = unique(mgrp[:, [:time, :predicted]])
                    sort!(pred, :time)
                    plot!(panel, pred.time, pred.predicted;
                        color = model_colors[model],
                        lw = 1.8,
                        alpha = 0.85,
                        label = "",
                    )
                end
                push!(panels, panel)
            end
            legend_panel = plot(; axis = false, grid = false, framestyle = :none, legend = :left, legendfontsize = 9)
            scatter!(legend_panel, [NaN], [NaN]; label = "Observed mean", color = :black, ms = 4, markerstrokewidth = 0)
            for model in model_order
                plot!(legend_panel, [NaN], [NaN];
                    color = model_colors[model],
                    lw = 1.8,
                    label = get(model_labels, model, model),
                )
            end
            push!(panels, legend_panel)
            p = plot(
                panels...;
                layout = (2, 2),
                size = (1400, 920),
                plot_title = "Joint dose fit: $(cell_label) $(den_label)",
                plot_titlefontsize = 13,
                margin = 4 * Plots.mm,
            )
            safe_name = IOUtils.sanitize_name("monoculture_treated_joint_dose_$(cell_label)_$(den_label)")
            savefig(p, joinpath(image_fig_dir, "$(safe_name).png"))
        end

        environment_groups = sort(
            collect(groupby(overlay_df, [:cell_line, :density]));
            by = grp -> begin
                cell = lowercase(_safe_string(first(grp.cell_line)))
                cell_order = occursin("naive", cell) ? 1 : occursin("cis", cell) ? 2 : 3
                density = _safe_string(first(grp.density))
                density_value = tryparse(Int, replace(lowercase(density), "k" => "000"))
                (cell_order, something(density_value, typemax(Int)), density)
            end,
        )
        best_panels = Any[]
        best_overlay_parts = DataFrame[]
        for environment_group in environment_groups
            cell_label = _safe_string(first(environment_group.cell_line); default = "pooled")
            den_label = _safe_string(first(environment_group.density); default = "pooled")
            ranking_group = joint_df[
                String.(joint_df.cell_line) .== cell_label,
                :,
            ]
            isempty(ranking_group) && continue
            winner = ranking_group[argmin(ranking_group.bic), :]
            winner_model = String(winner.model)
            winner_bic = Float64(winner.bic)
            best_overlay = environment_group[String.(environment_group.model) .== winner_model, :]
            isempty(best_overlay) && continue
            push!(best_overlay_parts, DataFrame(best_overlay))

            dose_groups = sort(
                collect(groupby(best_overlay, :dose));
                by = dose_group -> first(dose_group.effect_level),
            )
            for (dose_index, dose_group) in enumerate(dose_groups)
                dose_label = "$(first(dose_group.ic_label)) ($(first(dose_group.dose)) uM)"
                winner_label = get(model_labels, winner_model, winner_model)
                panel_title = dose_index == 1 ?
                    "$(cell_label) $(den_label): $(dose_label)\n$(winner_label), BIC=$(round(winner_bic; digits = 1))" :
                    dose_label
                panel = plot(
                    xlabel = "Time (day)",
                    ylabel = dose_index == 1 ? "Cell count" : "",
                    title = panel_title,
                    legend = false,
                    titlefontsize = 9,
                )
                observed = unique(dose_group[:, [:time, :observed]])
                sort!(observed, :time)
                scatter!(panel, observed.time, observed.observed;
                    label = "Observed",
                    color = :black,
                    ms = 3.8,
                    markerstrokewidth = 0,
                )
                predicted = unique(dose_group[:, [:time, :predicted]])
                sort!(predicted, :time)
                plot!(panel, predicted.time, predicted.predicted;
                    label = winner_label,
                    color = :steelblue,
                    lw = 2.8,
                )
                push!(best_panels, panel)
            end
        end
        if !isempty(best_panels)
            environment_count = div(length(best_panels), 3)
            best_figure = plot(
                best_panels...;
                layout = (environment_count, 3),
                size = (1500, 390 * environment_count),
                plot_title = "Best density-coupled joint-dose model in each treated environment",
                plot_titlefontsize = 14,
                margin = 4 * Plots.mm,
            )
            savefig(best_figure, joinpath(image_fig_dir, "monoculture_treated_best_joint_model_by_environment.png"))
            CSV.write(
                joinpath(fig_dir, "monoculture_treated_best_joint_model_by_environment_overlays.csv"),
                vcat(best_overlay_parts...),
            )
        end
    end
    return joint_df
end

function _density_nominal_value(label::AbstractString)
    parsed = tryparse(Float64, replace(lowercase(String(label)), "k" => "000"))
    return something(parsed, NaN)
end

function _fixed_day0_total(label::AbstractString)
    normalized = lowercase(strip(String(label)))
    occursin("20k", normalized) && return 67.0
    occursin("30k", normalized) && return 100.0
    error("No fixed day-zero seeding anchor is defined for density $(label)")
end

function _render_treated_pooling_graph_grid(overlay::DataFrame, status::DataFrame, out)
    isempty(overlay) && return nothing
    panels = Any[]
    best_parts = DataFrame[]
    cell_lines = sort(unique(String.(overlay.cell_line)); by = cell -> occursin("naive", lowercase(cell)) ? 1 : 2)
    for cell_line in cell_lines
        status_row = status[String.(status.cell_line) .== cell_line, :]
        isempty(status_row) && continue
        winner_model = String(first(status_row.winning_model))
        winner_pooling = String(first(status_row.winning_pooling_mode))
        cell_overlay = overlay[
            (String.(overlay.cell_line) .== cell_line) .&
            (String.(overlay.model) .== winner_model) .&
            (String.(overlay.pooling_mode) .== winner_pooling),
            :,
        ]
        push!(best_parts, DataFrame(cell_overlay))
        densities = sort(unique(String.(cell_overlay.density)); by = _density_nominal_value)
        for density in densities
            density_overlay = cell_overlay[String.(cell_overlay.density) .== density, :]
            dose_groups = sort(collect(groupby(density_overlay, :dose)); by = group -> first(group.effect_level))
            for (dose_index, dose_group) in enumerate(dose_groups)
                observed = unique(dose_group[:, [:time, :observed]])
                predicted = unique(dose_group[:, [:time, :predicted]])
                sort!(observed, :time)
                sort!(predicted, :time)
                title = "$(first(dose_group.ic_label)) ($(first(dose_group.dose)) uM)"
                dose_index == 1 && (title = "$(cell_line) $(density)\n$(winner_model) | $(winner_pooling)\n$(title)")
                panel = plot(
                    predicted.time,
                    predicted.predicted;
                    color = :steelblue,
                    lw = 2.6,
                    label = "",
                    title = title,
                    titlefontsize = 8,
                    xlabel = "Time (day)",
                    ylabel = dose_index == 1 ? "Cell count" : "",
                    legend = false,
                )
                scatter!(panel, observed.time, observed.observed; color = :black, ms = 3.5, markerstrokewidth = 0, label = "")
                push!(panels, panel)
            end
        end
    end
    isempty(panels) && return nothing
    image_dir = joinpath(out.images, "figures")
    csv_dir = joinpath(out.csv, "figures")
    mkpath(image_dir)
    mkpath(csv_dir)
    figure = plot(
        panels...;
        layout = (div(length(panels), 3), 3),
        size = (1500, 390 * div(length(panels), 3)),
        plot_title = "Best density-aware joint-dose model in each treated environment",
        plot_titlefontsize = 14,
        margin = 4 * Plots.mm,
    )
    image_path = joinpath(image_dir, "monoculture_treated_best_joint_model_by_environment.png")
    savefig(figure, image_path)
    CSV.write(joinpath(csv_dir, "monoculture_treated_best_joint_model_by_environment_overlays.csv"), vcat(best_parts...; cols = :union))
    return image_path
end

function _run_density_aware_treated_monoculture_fitting(
    decoded::DataFrame,
    out;
    start::AbstractString = pwd(),
    max_time_per_fit::Float64 = 45.0,
)
    time_col = :time
    :day in propertynames(decoded) && !(:time in propertynames(decoded)) && (time_col = :day)
    area_col = :count in propertynames(decoded) ? :count : error("Joint treated fitting requires a count column")
    baseline_map = _load_untreated_monoculture_cellline_baselines(; start = start)
    baseline_map === nothing && error("Density-aware untreated baseline artifact is required before treated fitting")

    _write_treated_trajectory_summary(decoded, out, time_col, area_col)
    rows = NamedTuple[]
    parameter_rows = NamedTuple[]
    initial_rows = NamedTuple[]
    overlay_parts = DataFrame[]
    profile_parts = DataFrame[]
    identifiability_parts = DataFrame[]
    inheritance_rows = NamedTuple[]

    for cell_group in groupby(decoded, :cell_line)
        cell_line = _safe_string(first(cell_group.cell_line); default = "pooled")
        density_labels = sort(unique(String.(cell_group.density)); by = _density_nominal_value)
        length(density_labels) == 2 || error("$(cell_line) requires exactly two starting densities; found $(density_labels)")
        density_index = Dict(label => index for (index, label) in enumerate(density_labels))

        baseline_by_density = Dict{String,NamedTuple}()
        for density in density_labels
            key = (cell_line, density)
            haskey(baseline_map.spec_by_cell_density, key) || error("Missing untreated baseline for $(cell_line) $(density)")
            baseline = baseline_map.spec_by_cell_density[key]
            baseline.inheritance_allowed || error("Untreated density pooling was inadequate for $(cell_line); treated inheritance stopped")
            baseline_by_density[density] = baseline
            push!(inheritance_rows, (
                treated_cell_line = cell_line,
                density = density,
                untreated_growth_family = String(baseline.model),
                untreated_r = Float64(baseline.r),
                untreated_K = Float64(baseline.K),
                untreated_shape_parameter = String(baseline.shape_parameter),
                untreated_shape_value = Float64(baseline.shape_value),
                baseline_source = String(baseline_map.path),
                inheritance_mode = "exact_fixed_family_and_density_specific_parameters",
            ))
        end
        growth_families = unique([baseline_by_density[density].model for density in density_labels])
        length(growth_families) == 1 || error("Untreated winner must use one growth family across densities for $(cell_line)")

        trajectory_groups = sort(
            collect(groupby(cell_group, [:density, :dose]));
            by = group -> (
                density_index[String(first(group.density))],
                Float64(first(group.dose)),
            ),
        )
        datasets = NamedTuple[]
        doses = Float64[]
        effect_levels = Float64[]
        ic_labels = String[]
        initial_counts = Float64[]
        trajectory_densities = String[]
        trajectory_density_indices = Int[]
        trajectory_baselines = NamedTuple[]

        for trajectory in trajectory_groups
            density = String(first(trajectory.density))
            dose = Float64(first(trajectory.dose))
            observed = combine(groupby(trajectory, time_col), area_col => mean => :y_mean)
            sort!(observed, time_col)
            nrow(observed) >= 3 || continue
            y = Float64.(observed.y_mean)
            first_observed_value = max(first(y), 1.0)
            fixed_u0 = _fixed_day0_total(density)
            residual_scale = max(maximum(y), 1.0)
            metadata = _treated_dose_metadata(dose)
            baseline = baseline_by_density[density]
            push!(datasets, (x = Float64.(observed[!, time_col]), y = y, state_index = length(datasets) + 1, residual_scale = residual_scale))
            push!(doses, dose)
            push!(effect_levels, metadata.effect_level)
            push!(ic_labels, metadata.ic_label)
            push!(initial_counts, fixed_u0)
            push!(trajectory_densities, density)
            push!(trajectory_density_indices, density_index[density])
            push!(trajectory_baselines, baseline)
            push!(initial_rows, (
                cell_line = cell_line,
                density = density,
                dose = dose,
                ic_label = metadata.ic_label,
                nominal_density = _density_nominal_value(density),
                first_observed_value = first_observed_value,
                fixed_u0 = fixed_u0,
                u0_time_day = 0.0,
                u0_over_K = fixed_u0 / baseline.K,
                residual_scale = residual_scale,
                untreated_growth_family = baseline.model,
                untreated_r = baseline.r,
                untreated_K = baseline.K,
                initial_condition_strategy = "fixed_nominal_day0_67_for_20k_100_for_30k",
            ))
        end
        length(datasets) == 6 || error("Expected six treated trajectories for $(cell_line), found $(length(datasets))")
        joint_group = "densities=$(join(density_labels, "/"));doses=$(join(sort(unique(doses)), "/"))"

        for pooling_mode in ("shared", "partial_5pct", "independent_diagnostic")
            specs = _joint_treated_model_specs(
                doses,
                effect_levels,
                trajectory_density_indices,
                trajectory_baselines,
                0.0,
                pooling_mode,
            )
            for model_name in sort(collect(keys(specs)))
                spec = specs[model_name]
                println("  Density-aware treated fit $(cell_line), $(model_name), $(pooling_mode)...")
                model_datasets, model_u0, u0_builder = _joint_model_inputs(spec, datasets, initial_counts, trajectory_density_indices)
                explicit_profiles = Dict{Symbol,Vector{Float64}}()
                for parameter_name in spec.param_names
                    occursin("ec50_effect", String(parameter_name)) && (explicit_profiles[parameter_name] = [1.0, 1.5, 2.0])
                    occursin("lambda", String(parameter_name)) && (explicit_profiles[parameter_name] = [5.0, 7.5, 10.0])
                end
                profile_parameters = filter(name -> !occursin("log_contrast", String(name)), spec.param_names)
                model_maxiters = pooling_mode == "independent_diagnostic" ?
                    max(250, Int(round(max_time_per_fit * 12))) :
                    max(350, Int(round(max_time_per_fit * 18)))
                profiled = try
                    GrowthParameterEstimation.profile_joint_fit_bounds(
                        spec.model,
                        model_datasets,
                        model_u0,
                        spec.p0;
                        bounds = spec.bounds,
                        parameter_names = spec.param_names,
                        explicit_upper_profiles = explicit_profiles,
                        profile_parameters = profile_parameters,
                        solver = Tsit5(),
                        u0_builder = u0_builder,
                        maxiters = model_maxiters,
                        reltol = 1e-7,
                        abstol = 1e-7,
                        optimizer = :nelder_mead,
                        initial_time = 0.0,
                    )
                catch error_value
                    @warn "Density-aware treated fit failed" cell_line model_name pooling_mode exception = error_value
                    nothing
                end
                profiled === nothing && continue
                fit = profiled.fit
                isfinite(Float64(fit.bic)) && isfinite(Float64(fit.raw_sse)) && Float64(fit.raw_sse) < 9.99e11 || continue

                profile_df = DataFrame(profiled.profile)
                if !isempty(profile_df)
                    profile_df.cell_line = fill(cell_line, nrow(profile_df))
                    profile_df.model = fill(model_name, nrow(profile_df))
                    profile_df.pooling_mode = fill(pooling_mode, nrow(profile_df))
                    push!(profile_parts, profile_df)
                end
                identifiability_df = DataFrame(profiled.identifiability)
                identifiability_df.cell_line = fill(cell_line, nrow(identifiability_df))
                identifiability_df.model = fill(model_name, nrow(identifiability_df))
                identifiability_df.pooling_mode = fill(pooling_mode, nrow(identifiability_df))
                push!(identifiability_parts, identifiability_df)
                boundary_issue = any(identifiability_df.identifiability .!= "interior")

                for density in density_labels
                    trajectory_index = findfirst(==(density), trajectory_densities)
                    effective = spec.local_parameters(fit.params, trajectory_index)
                    other_index = findfirst(!=(density), density_labels)
                    other_trajectory = findfirst(==(density_labels[other_index]), trajectory_densities)
                    other_effective = spec.local_parameters(fit.params, other_trajectory)
                    for parameter_index in eachindex(spec.base_param_names)
                        center = pooling_mode == "independent_diagnostic" ?
                            sqrt(max(effective[parameter_index], 0.0) * max(other_effective[parameter_index], 0.0)) :
                            Float64(fit.params[parameter_index])
                        deviation = abs(center) > 1e-12 ? 100 * (Float64(effective[parameter_index]) / center - 1) : 0.0
                        lower, upper = spec.base_bounds[parameter_index]
                        push!(parameter_rows, (
                            cell_line = cell_line,
                            model = model_name,
                            pooling_mode = pooling_mode,
                            density = density,
                            parameter = String(spec.base_param_names[parameter_index]),
                            center_value = center,
                            effective_value = Float64(effective[parameter_index]),
                            deviation_pct = deviation,
                            lower_bound = Float64(lower),
                            upper_bound = Float64(upper),
                            bound_position = (Float64(effective[parameter_index]) - lower) / max(upper - lower, eps()),
                            parameter_scope = pooling_mode == "partial_5pct" && parameter_index in spec.amplitude_indices ? "density_partial_5pct" : pooling_mode,
                        ))
                    end
                end

                representative = spec.local_parameters(fit.params, 1)
                fate = _joint_ic75_long_term_fate(model_name, representative, trajectory_baselines[1].r)
                push!(rows, (
                    model = model_name,
                    pooling_mode = pooling_mode,
                    eligible_for_inheritance = pooling_mode != "independent_diagnostic",
                    cell_line = cell_line,
                    density = join(density_labels, "/"),
                    dose = NaN,
                    bic = Float64(fit.bic),
                    ssr = Float64(fit.raw_sse),
                    scaled_ssr = Float64(fit.scaled_sse),
                    n_parameters = length(fit.params),
                    n_points = sum(length(dataset.x) for dataset in datasets),
                    n_densities = 2,
                    n_doses = length(unique(doses)),
                    ic_mapping = join(unique(["$(doses[i])=$(ic_labels[i])" for i in eachindex(doses)]), ";"),
                    joint_group = joint_group,
                    growth_family = first(growth_families),
                    growth_inheritance = "exact_density_specific_untreated_family_and_parameters",
                    residual_scaling = "trajectory_peak",
                    initial_condition_strategy = "fixed_nominal_day0_67_for_20k_100_for_30k",
                    effect_parameter_scope = pooling_mode == "partial_5pct" ? "amplitude_only_plus_or_minus_5pct" : pooling_mode,
                    boundary_issue = boundary_issue,
                    model_behavior = spec.behavior,
                    anchored_hill_limitation = model_name == "joint_anchored_hill_kill" ? "monotonic benchmark; cannot represent rise_then_decline" : "",
                    ic75_low_density_net_rate = fate.net_rate,
                    ic75_long_term_fate = fate.fate,
                    params = string((names = spec.param_names, values = fit.params, package_api = "GrowthParameterEstimation.profile_joint_fit_bounds/run_joint_fit")),
                ))

                for (trajectory_index, dataset) in enumerate(datasets)
                    fixed_u0 = initial_counts[trajectory_index]
                    overlay_length = length(dataset.x) + 1
                    push!(overlay_parts, DataFrame(
                        time = vcat(0.0, dataset.x),
                        observed = vcat(fixed_u0, dataset.y),
                        predicted = vcat(fixed_u0, fit.predictions[trajectory_index]),
                        fixed_day0_anchor = vcat(true, fill(false, length(dataset.x))),
                        model = fill(model_name, overlay_length),
                        pooling_mode = fill(pooling_mode, overlay_length),
                        cell_line = fill(cell_line, overlay_length),
                        density = fill(trajectory_densities[trajectory_index], overlay_length),
                        dose = fill(doses[trajectory_index], overlay_length),
                        ic_label = fill(ic_labels[trajectory_index], overlay_length),
                        effect_level = fill(effect_levels[trajectory_index], overlay_length),
                        growth_family = fill(first(growth_families), overlay_length),
                        joint_group = fill(joint_group, overlay_length),
                    ))
                end
            end
        end
    end

    ranking = isempty(rows) ? DataFrame() : sort!(DataFrame(rows), [:cell_line, :bic])
    isempty(ranking) && error("No finite density-aware treated fits were produced")
    pooling_summary = GrowthParameterEstimation.summarize_pooling_bic(ranking; top_n = 5)
    status = pooling_summary.status
    any(.!Bool.(status.inheritance_allowed)) && error("Independent treated fits beat eligible pooling by at least 10 BIC points")

    initial_df = DataFrame(initial_rows)
    ratio_rows = NamedTuple[]
    for group in groupby(initial_df, [:cell_line, :dose])
        sort!(group, :nominal_density)
        nrow(group) == 2 || continue
        push!(ratio_rows, (
            cell_line = String(first(group.cell_line)),
            dose = Float64(first(group.dose)),
            nominal_density_ratio = Float64(group.nominal_density[2] / group.nominal_density[1]),
            fixed_u0_ratio = Float64(group.fixed_u0[2] / group.fixed_u0[1]),
            first_observed_ratio = Float64(group.first_observed_value[2] / group.first_observed_value[1]),
        ))
    end
    initial_df = leftjoin(initial_df, DataFrame(ratio_rows); on = [:cell_line, :dose])

    figures_csv = joinpath(out.csv, "figures")
    mkpath(figures_csv)
    CSV.write(joinpath(out.csv, "monoculture_treated_joint_dose_model_ranking.csv"), ranking)
    CSV.write(joinpath(out.csv, "monoculture_treated_pooling_model_ranking.csv"), ranking)
    CSV.write(joinpath(out.csv, "monoculture_treated_pooling_top5.csv"), pooling_summary.ranking)
    CSV.write(joinpath(out.csv, "monoculture_treated_pooling_status.csv"), status)
    CSV.write(joinpath(out.csv, "monoculture_treated_joint_cell_line_top5.csv"), pooling_summary.ranking)
    CSV.write(joinpath(out.csv, "monoculture_treated_joint_parameter_estimates.csv"), DataFrame(parameter_rows))
    CSV.write(joinpath(out.csv, "monoculture_treated_joint_initial_condition_diagnostics.csv"), initial_df)
    CSV.write(joinpath(out.csv, "monoculture_treated_boundary_profiles.csv"), isempty(profile_parts) ? DataFrame() : vcat(profile_parts...; cols = :union))
    CSV.write(joinpath(out.csv, "monoculture_treated_identifiability.csv"), isempty(identifiability_parts) ? DataFrame() : vcat(identifiability_parts...; cols = :union))
    CSV.write(joinpath(out.csv, "monoculture_treated_inheritance_audit.csv"), unique(DataFrame(inheritance_rows)))
    overlay_df = vcat(overlay_parts...; cols = :union)
    CSV.write(joinpath(figures_csv, "monoculture_treated_joint_dose_overlays.csv"), overlay_df)
    _render_treated_pooling_graph_grid(overlay_df, status, out)
    return ranking
end

# ---------------------------------------------------------------------------
# Per-dose model fitting using local model specs
# ---------------------------------------------------------------------------

"""
Fit a single LocalModelSpec to one dose group (x, y time-series).
Returns (params, bic, ssr).
"""
function _fit_spec_to_dose(
    spec::ModelRegistry.LocalModelSpec,
    x::Vector{Float64},
    y::Vector{Float64},
    dose::Float64,
    r0::Float64,
    K0::Float64;
    solver   = Rodas5(),
    max_time::Float64 = 45.0,
)
    p0  = spec.p0_factory(r0, K0, dose)
    obs = spec.observable          # (u::Vector) -> scalar
    nst = spec.n_states

    # Constant-exposure closure, baked-in for this dose group
    exposure_fn = Returns(dose)
    ode4! = (du, u, p, t) -> spec.ode!(du, u, p, t, exposure_fn)

    u0    = zeros(nst)
    u0[1] = max(y[1], 1.0)
    tspan = (x[1], x[end])
    prob  = ODEProblem(ode4!, u0, tspan, p0)

    n_dims = length(p0)

    function loss(p_vec)
        pv = Vector{Float64}(p_vec)
        try
            sol = solve(remake(prob; p = pv), solver;
                        reltol = 1e-6, abstol = 1e-6,
                        saveat = x, maxiters = 50_000)
            sol.retcode == ReturnCode.Success || return 1e12
            yhat = [obs(u) for u in sol.u]
            length(yhat) == length(y) || return 1e12
            return sum((y .- yhat) .^ 2)
        catch
            return 1e12
        end
    end

    result = bboptimize(
        loss;
        SearchRange  = collect(zip(first.(spec.bounds), last.(spec.bounds))),
        NumDimensions = n_dims,
        Method       = :de_rand_1_bin,
        MaxTime      = max_time,
        TraceMode    = :silent,
    )
    p_opt = Vector{Float64}(result.archive_output.best_candidate)

    # Compute final BIC/SSR with tighter tolerances
    bic, ssr = try
        sol2 = solve(remake(prob; p = p_opt), solver;
                     reltol = 1e-10, abstol = 1e-10, saveat = x)
        yhat2 = [obs(u) for u in sol2.u]
        ssr2 = sum((y .- yhat2) .^ 2)
        n, k = length(x), length(p_opt)
        bic2 = n * log(max(ssr2, 1e-20) / n) + k * log(n)
        bic2, ssr2
    catch
        1e12, 1e12
    end

    return (params = p_opt, bic = bic, ssr = ssr)
end

"""
Fit all local model specs to each dose group in `decoded` and return a
ranked DataFrame with columns: model, dose, bic, ssr, params.
"""
function _run_treated_fitting(
    decoded::DataFrame,
    out;
    start::AbstractString = pwd(),
    untreated_baseline = nothing,
    max_time_per_fit::Float64 = 45.0,
)
    r0 = isnothing(untreated_baseline) ? 0.3 : untreated_baseline.params[1]
    K0 = isnothing(untreated_baseline) ? 3000.0 : untreated_baseline.params[2]
    baseline_map = _load_untreated_monoculture_cellline_baselines(; start = start)

    specs = ModelRegistry.local_model_specs()

    # Detect the measurement column (area or first numeric non-metadata col)
    meta_cols = Set([:dose, :time, :day, :replicate, :cell_line, :density, :sample_id, :file, :source_file, :condition, :mix])
    numeric_cols = [c for c in propertynames(decoded)
                    if !(c in meta_cols) && eltype(decoded[!, c]) <: Real]
    time_col = :time
    :day in propertynames(decoded) && !(:time in propertynames(decoded)) && (time_col = :day)
    area_col = isempty(numeric_cols) ? error("No numeric measurement column found in decoded data") :
                                       numeric_cols[1]

    rows = NamedTuple{(:model, :dose, :bic, :ssr, :params, :cell_line, :density, :mix),
                      Tuple{String, Float64, Float64, Float64, String, String, String, String}}[]

    colset = Set(Symbol.(names(decoded)))
    cell_col = :cell_line in colset ? :cell_line : (:cellline in colset ? :cellline : nothing)
    density_col = :density in colset ? :density : nothing
    mix_col = :mix in colset ? :mix : nothing
    has_cell = cell_col !== nothing
    has_density = density_col !== nothing
    has_mix = mix_col !== nothing

    gcols = Any[]
    has_cell && push!(gcols, cell_col)
    has_density && push!(gcols, density_col)
    has_mix && push!(gcols, mix_col)
    push!(gcols, :dose)

    for grp in groupby(decoded, gcols)
        dose_val = Float64(first(grp.dose))
        grouped = combine(groupby(grp, time_col), area_col => mean => :y_mean)
        sort!(grouped, time_col)
        x = Float64.(grouped[!, time_col])
        y = Float64.(grouped[!, :y_mean])

        cell_label = has_cell ? _safe_string(first(grp[!, cell_col]); default = "pooled") : "pooled"
        den_label = has_density ? _safe_string(first(grp[!, density_col]); default = "pooled") : "pooled"
        mix_label = has_mix ? _safe_string(first(grp[!, mix_col]); default = "pooled") : "pooled"
        r_group, K_group = _baseline_rk_for_group(baseline_map, cell_label, den_label, r0, K0)

        for spec in specs
            println("  Fitting $(spec.name) at dose=$(dose_val), cell=$(cell_label), density=$(den_label)...")
            fit = _fit_spec_to_dose(spec, x, y, dose_val, r_group, K_group; max_time = max_time_per_fit)
            push!(rows, (
                model  = spec.name,
                dose   = dose_val,
                bic    = fit.bic,
                ssr    = fit.ssr,
                params = string((fit = fit.params, r_anchor = r_group, K_anchor = K_group)),
                cell_line = cell_label,
                density = den_label,
                mix = mix_label,
            ))
        end
    end

    # Additional fast anchored models using untreated monoculture r/K per cell line and density.
    # These are fitted by (cell_line, density, dose) groups and only estimate drug-effect terms.
    if has_cell
        gcols = has_density ? [cell_col, density_col, :dose] : [cell_col, :dose]
        has_mix && push!(gcols, mix_col)
        for grp in groupby(decoded, gcols)
            dose_val = Float64(first(grp.dose))
            grouped = combine(groupby(grp, time_col), area_col => mean => :y_mean)
            sort!(grouped, time_col)
            x = Float64.(grouped[!, time_col])
            y = Float64.(grouped[!, :y_mean])

            cell_label = _safe_string(first(grp[!, cell_col]); default = "pooled")
            den_label = has_density ? _safe_string(first(grp[!, density_col])) : ""
            r_anchor, K_anchor = _baseline_rk_for_group(baseline_map, cell_label, den_label, r0, K0)
            anchored_budget = max(1.0, min(max_time_per_fit / 3, 4.0))

            println("  Fitting anchored_linear_kill at dose=$(dose_val), cell=$(cell_label), density=$(den_label)...")
            lin_fit = _fit_anchored_linear_kill(x, y, dose_val, r_anchor, K_anchor; max_time = anchored_budget)
            mix_label = :mix in propertynames(grp) ? _safe_string(first(grp.mix); default = "pooled") : "pooled"
            push!(rows, (
                model = "anchored_linear_kill",
                dose = dose_val,
                bic = lin_fit.bic,
                ssr = lin_fit.ssr,
                params = string((k_kill = lin_fit.params[1], r_anchor = r_anchor, K_anchor = K_anchor)),
                cell_line = cell_label,
                density = den_label,
                mix = mix_label,
            ))

            println("  Fitting anchored_hill_kill at dose=$(dose_val), cell=$(cell_label), density=$(den_label)...")
            hill_fit = _fit_anchored_hill_kill(x, y, dose_val, r_anchor, K_anchor; max_time = anchored_budget)
            push!(rows, (
                model = "anchored_hill_kill",
                dose = dose_val,
                bic = hill_fit.bic,
                ssr = hill_fit.ssr,
                params = string((emax = hill_fit.params[1], ec50 = hill_fit.params[2], hill_n = hill_fit.params[3], r_anchor = r_anchor, K_anchor = K_anchor)),
                cell_line = cell_label,
                density = den_label,
                mix = mix_label,
            ))
        end
    end

    ranking_df = sort!(DataFrame(rows), :bic)
    condition_name = :condition in propertynames(decoded) ? _safe_string(first(decoded.condition); default = "") : ""
    if condition_name == "monoculture_treated"
        joint_df = _run_joint_treated_monoculture_fitting(
            decoded,
            out;
            start = start,
            untreated_baseline = untreated_baseline,
            max_time_per_fit = max_time_per_fit,
        )
    end

    return ranking_df
end

# ---------------------------------------------------------------------------
# Untreated custom model fitting
# ---------------------------------------------------------------------------

function _fit_untreated_model(
    model_name::String,
    x::Vector{Float64},
    y::Vector{Float64};
    solver = Rodas5(),
    max_time::Float64 = 8.0,
)
    function model_spec(name::String)
        if name == "logistic_growth"
            ode! = function (du, u, p, t)
                r, K = p
                N = max(u[1], 0.0)
                du[1] = r * N * max(0.0, 1 - N / max(K, 1e-8))
            end
            return (ode! = ode!, p0 = [0.3, max(y[end], 100.0)], bounds = [(1e-6, 5.0), (1.0, 1e7)])
        elseif name == "gompertz_growth"
            ode! = function (du, u, p, t)
                r, K = p
                N = max(u[1], 1e-8)
                du[1] = r * N * log(max(K, 1e-8) / N)
            end
            return (ode! = ode!, p0 = [0.3, max(y[end], 100.0)], bounds = [(1e-6, 5.0), (1.0, 1e7)])
        elseif name == "logistic_simple_death"
            ode! = function (du, u, p, t)
                r, K, d = p
                N = max(u[1], 0.0)
                growth = r * N * max(0.0, 1 - N / max(K, 1e-8))
                du[1] = growth - d * N
            end
            return (ode! = ode!, p0 = [0.4, max(y[end], 100.0), 0.05], bounds = [(1e-6, 5.0), (1.0, 1e7), (0.0, 2.0)])
        elseif name == "allee_growth"
            ode! = function (du, u, p, t)
                r, K, A = p
                N = max(u[1], 0.0)
                allee_term = (N / max(A, 1e-8)) - 1
                du[1] = r * N * max(0.0, 1 - N / max(K, 1e-8)) * allee_term
            end
            return (ode! = ode!, p0 = [0.4, max(y[end], 100.0), max(0.2 * y[1], 1.0)], bounds = [(1e-6, 5.0), (1.0, 1e7), (1e-6, 1e6)])
        elseif name == "theta_logistic_growth"
            ode! = function (du, u, p, t)
                r, K, theta = p
                N = max(u[1], 0.0)
                du[1] = r * N * max(0.0, 1 - (N / max(K, 1e-8))^max(theta, 1e-8))
            end
            return (ode! = ode!, p0 = [0.3, max(y[end], 100.0), 1.0], bounds = [(1e-6, 5.0), (1.0, 1e7), (0.1, 4.0)])
        elseif name == "generalized_logistic_growth"
            ode! = function (du, u, p, t)
                r, K, nu = p
                N = max(u[1], 0.0)
                du[1] = r * N * max(0.0, 1 - (N / max(K, 1e-8))^max(nu, 1e-8))
            end
            return (ode! = ode!, p0 = [0.3, max(y[end], 100.0), 1.0], bounds = [(1e-6, 5.0), (1.0, 1e7), (0.05, 8.0)])
        else
            error("Unsupported untreated model: $(name)")
        end
    end

    spec = model_spec(model_name)
    u0 = [max(y[1], 1.0)]
    prob = ODEProblem(spec.ode!, u0, (x[1], x[end]), spec.p0)

    function loss(p_vec)
        pv = Vector{Float64}(p_vec)
        try
            sol = solve(remake(prob; p = pv), solver;
                        reltol = 1e-6, abstol = 1e-6,
                        saveat = x, maxiters = 50_000)
            sol.retcode == ReturnCode.Success || return 1e12
            yhat = first.(sol.u)
            return sum((y .- yhat) .^ 2)
        catch
            return 1e12
        end
    end

    result = bboptimize(
        loss;
        SearchRange = spec.bounds,
        NumDimensions = length(spec.p0),
        Method = :de_rand_1_bin,
        MaxTime = max_time,
        TraceMode = :silent,
    )
    p_opt = Vector{Float64}(result.archive_output.best_candidate)

    bic, ssr = try
        sol2 = solve(remake(prob; p = p_opt), solver;
                     reltol = 1e-10, abstol = 1e-10, saveat = x)
        yhat2 = first.(sol2.u)
        ssr2 = sum((y .- yhat2) .^ 2)
        n, k = length(x), length(p_opt)
        n * log(max(ssr2, 1e-20) / n) + k * log(n), ssr2
    catch
        1e12, 1e12
    end

    return (params = p_opt, bic = bic, ssr = ssr)
end

# ---------------------------------------------------------------------------
# Density-aware untreated monoculture fitting
# ---------------------------------------------------------------------------

function _untreated_base_spec(model_name::String, max_y::Float64, initial_counts::Vector{Float64})
    k0 = max(1.05 * max_y, 1.0)
    k_bounds = (max(0.75 * max_y, 1.0), max(3.0 * max_y, k0 * 1.1))
    if model_name == "logistic_growth" || model_name == "gompertz_growth"
        return (param_names = [:r, :K], p0 = [0.5, k0], bounds = [(1e-4, 2.0), k_bounds])
    elseif model_name == "theta_logistic_growth"
        return (param_names = [:r, :K, :theta], p0 = [0.5, k0, 1.0], bounds = [(1e-4, 2.0), k_bounds, (0.1, 8.0)])
    elseif model_name == "logistic_simple_death"
        return (param_names = [:r, :K, :death_rate], p0 = [0.5, k0, 0.05], bounds = [(1e-4, 2.0), k_bounds, (0.0, 2.0)])
    elseif model_name == "allee_growth"
        return (param_names = [:r, :K, :allee_threshold], p0 = [0.5, k0, max(0.2 * minimum(initial_counts), 1.0)], bounds = [(1e-4, 2.0), k_bounds, (1e-6, max_y)])
    end
    error("Unsupported density-aware untreated model: $(model_name)")
end

function _pooling_spec(base, pooling_mode::String, density_labels::AbstractVector{<:AbstractString})
    length(density_labels) == 2 || error("Density-aware monoculture fitting requires exactly two densities")
    if pooling_mode == "shared"
        return (p0 = copy(base.p0), bounds = copy(base.bounds), param_names = copy(base.param_names))
    elseif pooling_mode == "partial_5pct"
        return (
            p0 = vcat(base.p0, [0.0, 0.0]),
            bounds = vcat(base.bounds, [(-DENSITY_LOG_CONTRAST_BOUND, DENSITY_LOG_CONTRAST_BOUND), (-DENSITY_LOG_CONTRAST_BOUND, DENSITY_LOG_CONTRAST_BOUND)]),
            param_names = vcat(base.param_names, [:log_contrast_r, :log_contrast_K]),
        )
    elseif pooling_mode == "independent_diagnostic"
        names = Symbol[]
        for density in density_labels, name in base.param_names
            push!(names, Symbol("$(name)_$(density)"))
        end
        return (
            p0 = vcat(base.p0, base.p0),
            bounds = vcat(base.bounds, base.bounds),
            param_names = names,
        )
    end
    error("Unsupported pooling mode: $(pooling_mode)")
end

function _effective_untreated_params(p, base, pooling_mode::String, density_index::Int)
    nbase = length(base.p0)
    if pooling_mode == "shared"
        return collect(p[1:nbase])
    elseif pooling_mode == "partial_5pct"
        effective = collect(p[1:nbase])
        r_pair = GrowthParameterEstimation.symmetric_relative_pair(effective[1], p[nbase + 1])
        k_pair = GrowthParameterEstimation.symmetric_relative_pair(effective[2], p[nbase + 2])
        effective[1] = density_index == 1 ? r_pair.low : r_pair.high
        effective[2] = density_index == 1 ? k_pair.low : k_pair.high
        return effective
    elseif pooling_mode == "independent_diagnostic"
        first_index = (density_index - 1) * nbase + 1
        return collect(p[first_index:(first_index + nbase - 1)])
    end
    error("Unsupported pooling mode: $(pooling_mode)")
end

function _untreated_joint_ode(model_name::String, base, pooling_mode::String)
    return function (du, u, p, t)
        for i in eachindex(u)
            effective = _effective_untreated_params(p, base, pooling_mode, i)
            N = max(u[i], zero(u[i]))
            r, K = effective[1], max(effective[2], oftype(effective[2], 1e-8))
            if model_name == "logistic_growth"
                du[i] = r * N * max(0.0, 1 - N / K)
            elseif model_name == "gompertz_growth"
                du[i] = r * max(N, 1e-8) * log(K / max(N, 1e-8))
            elseif model_name == "theta_logistic_growth"
                theta = max(effective[3], oftype(effective[3], 1e-8))
                du[i] = r * N * max(0.0, 1 - (N / K)^theta)
            elseif model_name == "logistic_simple_death"
                du[i] = r * N * max(0.0, 1 - N / K) - max(effective[3], 0.0) * N
            elseif model_name == "allee_growth"
                threshold = max(effective[3], oftype(effective[3], 1e-8))
                du[i] = r * N * max(0.0, 1 - N / K) * ((N / threshold) - 1)
            end
        end
    end
end

function _effective_parameter_rows(cell_line, model_name, pooling_mode, density_labels, base, profiled)
    rows = NamedTuple[]
    effective = [_effective_untreated_params(profiled.fit.params, base, pooling_mode, i) for i in eachindex(density_labels)]
    for (parameter_index, parameter_name) in enumerate(base.param_names)
        values = [Float64(params[parameter_index]) for params in effective]
        center = pooling_mode == "independent_diagnostic" ? sqrt(max(values[1], 1e-12) * max(values[2], 1e-12)) : Float64(profiled.fit.params[parameter_index])
        for (density_index, density) in enumerate(density_labels)
            push!(rows, (
                cell_line = cell_line,
                model = model_name,
                pooling_mode = pooling_mode,
                density = density,
                parameter = String(parameter_name),
                center_value = center,
                effective_value = values[density_index],
                deviation_percent = 100 * (values[density_index] / max(center, 1e-12) - 1),
                parameter_scope = pooling_mode == "shared" ? "shared_across_densities" : pooling_mode,
            ))
        end
    end
    return rows
end

function _render_untreated_pooling_graph_grid(overlay::DataFrame, top5::DataFrame, out)
    isempty(overlay) && return nothing
    palette_values = Plots.palette(:tab10)
    panels = Any[]
    for cell_line in sort(unique(String.(overlay.cell_line)); by = cell -> occursin("naive", lowercase(cell)) ? 1 : 2)
        selected = top5[String.(top5.cell_line) .== cell_line, :]
        combinations = [(String(row.model), String(row.pooling_mode)) for row in eachrow(selected)]
        densities = sort(unique(String.(overlay[String.(overlay.cell_line) .== cell_line, :density])); by = _density_nominal_value)
        for density in densities
            environment = overlay[
                (String.(overlay.cell_line) .== cell_line) .&
                (String.(overlay.density) .== density),
                :,
            ]
            observed = unique(environment[:, [:time, :observed]])
            sort!(observed, :time)
            panel = plot(
                title = "$(cell_line) $(density)",
                titlefontsize = 10,
                xlabel = "Time (day)",
                ylabel = "Cell count",
                legend = :bottomright,
                legendfontsize = 7,
            )
            scatter!(panel, observed.time, observed.observed; color = :black, ms = 3.5, markerstrokewidth = 0, label = "Observed mean")
            for (index, (model, pooling_mode)) in enumerate(combinations)
                curve = environment[
                    (String.(environment.model) .== model) .&
                    (String.(environment.pooling_mode) .== pooling_mode),
                    :,
                ]
                isempty(curve) && continue
                predicted = unique(curve[:, [:time, :predicted]])
                sort!(predicted, :time)
                plot!(panel, predicted.time, predicted.predicted;
                    color = palette_values[mod1(index, length(palette_values))],
                    lw = index == 1 ? 2.8 : 1.7,
                    alpha = index == 1 ? 1.0 : 0.8,
                    label = "$(model) | $(pooling_mode)",
                )
            end
            push!(panels, panel)
        end
    end
    isempty(panels) && return nothing
    image_dir = joinpath(out.images, "figures")
    mkpath(image_dir)
    figure = plot(
        panels...;
        layout = (div(length(panels) + 1, 2), 2),
        size = (1450, 520 * div(length(panels) + 1, 2)),
        plot_title = "Untreated monoculture: top model and density-pooling comparisons",
        plot_titlefontsize = 14,
        margin = 5 * Plots.mm,
    )
    image_path = joinpath(image_dir, "monoculture_untreated_pooling_model_grid.png")
    savefig(figure, image_path)
    return image_path
end

function _run_density_aware_untreated_fitting(
    decoded::DataFrame,
    out;
    start::AbstractString = pwd(),
    max_time_per_fit::Float64 = 45.0,
)
    time_col = :time in propertynames(decoded) ? :time : :day
    area_col = :count in propertynames(decoded) ? :count : error("Untreated monoculture fitting requires count")
    model_names = ["logistic_growth", "gompertz_growth", "theta_logistic_growth", "logistic_simple_death", "allee_growth"]
    pooling_modes = ["shared", "partial_5pct", "independent_diagnostic"]
    ranking_rows = NamedTuple[]
    parameter_rows = NamedTuple[]
    profile_parts = DataFrame[]
    identifiability_parts = DataFrame[]
    initial_rows = NamedTuple[]
    overlay_parts = DataFrame[]

    for cell_group in groupby(decoded, :cell_line)
        cell_line = _safe_string(first(cell_group.cell_line); default = "pooled")
        density_groups = collect(groupby(cell_group, :density))
        sort!(density_groups, by = group -> _safe_string(first(group.density)))
        length(density_groups) == 2 || continue
        density_labels = [_safe_string(first(group.density)) for group in density_groups]
        datasets = NamedTuple[]
        initial_counts = Float64[]
        max_y = 0.0
        for (density_index, density_group) in enumerate(density_groups)
            means = combine(groupby(density_group, time_col), area_col => mean => :y_mean)
            sort!(means, time_col)
            x = Float64.(means[!, time_col])
            y = Float64.(means.y_mean)
            first_observed_value = max(first(y), 1.0)
            fixed_u0 = _fixed_day0_total(density_labels[density_index])
            residual_scale = max(maximum(y), 1.0)
            max_y = max(max_y, maximum(y))
            push!(initial_counts, fixed_u0)
            push!(datasets, (x = x, y = y, state_index = density_index, residual_scale = residual_scale))
            nominal = tryparse(Float64, replace(lowercase(density_labels[density_index]), "k" => "000"))
            push!(initial_rows, (
                cell_line = cell_line,
                density = density_labels[density_index],
                first_observed_value = first_observed_value,
                fixed_u0 = fixed_u0,
                u0_time_day = 0.0,
                residual_scale = residual_scale,
                nominal_density = something(nominal, NaN),
                first_observed_to_fixed_u0 = first_observed_value / fixed_u0,
                u0_strategy = "fixed_nominal_day0_67_for_20k_100_for_30k",
            ))
        end

        for model_name in model_names
            base = _untreated_base_spec(model_name, max_y, initial_counts)
            for pooling_mode in pooling_modes
                pooled = _pooling_spec(base, pooling_mode, density_labels)
                ode! = _untreated_joint_ode(model_name, base, pooling_mode)
                biological_parameters = [name for name in pooled.param_names if !startswith(String(name), "log_contrast")]
                println("  Density-aware untreated fit $(model_name), cell=$(cell_line), pooling=$(pooling_mode)...")
                profiled = GrowthParameterEstimation.profile_joint_fit_bounds(
                    ode!, datasets, initial_counts, pooled.p0;
                    bounds = pooled.bounds,
                    parameter_names = pooled.param_names,
                    profile_parameters = biological_parameters,
                    solver = Tsit5(),
                    optimizer = :nelder_mead,
                    maxiters = max(120, Int(round(max_time_per_fit * 45))),
                    reltol = 1e-7,
                    abstol = 1e-7,
                    initial_time = 0.0,
                )
                fit = profiled.fit
                valid = isfinite(Float64(fit.bic)) && abs(Float64(fit.sse)) < 9.99e11
                valid || continue
                boundary_issue = any(profiled.identifiability.identifiability .!= "interior")
                is_primary = model_name in PRIMARY_UNTREATED_MODELS
                push!(ranking_rows, (
                    model = model_name,
                    pooling_mode = pooling_mode,
                    eligible_for_inheritance = is_primary && pooling_mode != "independent_diagnostic",
                    diagnostic_model = !is_primary,
                    bic = Float64(fit.bic),
                    ssr = Float64(fit.raw_sse),
                    scaled_ssr = Float64(fit.scaled_sse),
                    params = string(fit.params),
                    cell_line = cell_line,
                    density = join(density_labels, "/"),
                    mix = "",
                    n_densities = length(density_labels),
                    n_points = sum(length(dataset.x) for dataset in datasets),
                    parameter_count = length(fit.params),
                    boundary_issue = boundary_issue,
                    accepted_boundary_expansions = nrow(profiled.profile) == 0 ? 0 : count(profiled.profile.accepted),
                    package_api = "GrowthParameterEstimation.profile_joint_fit_bounds",
                ))
                append!(parameter_rows, _effective_parameter_rows(cell_line, model_name, pooling_mode, density_labels, base, profiled))
                if nrow(profiled.profile) > 0
                    profile = copy(profiled.profile)
                    insertcols!(profile, 1, :cell_line => fill(cell_line, nrow(profile)), :model => fill(model_name, nrow(profile)), :pooling_mode => fill(pooling_mode, nrow(profile)))
                    push!(profile_parts, profile)
                end
                identifiability = copy(profiled.identifiability)
                insertcols!(identifiability, 1, :cell_line => fill(cell_line, nrow(identifiability)), :model => fill(model_name, nrow(identifiability)), :pooling_mode => fill(pooling_mode, nrow(identifiability)))
                push!(identifiability_parts, identifiability)
                for (density_index, dataset) in enumerate(datasets)
                    fixed_u0 = initial_counts[density_index]
                    overlay_length = length(dataset.x) + 1
                    push!(overlay_parts, DataFrame(
                        time = vcat(0.0, dataset.x),
                        observed = vcat(fixed_u0, dataset.y),
                        predicted = vcat(fixed_u0, fit.predictions[density_index]),
                        fixed_day0_anchor = vcat(true, fill(false, length(dataset.x))),
                        model = fill(model_name, overlay_length),
                        pooling_mode = fill(pooling_mode, overlay_length),
                        cell_line = fill(cell_line, overlay_length),
                        density = fill(density_labels[density_index], overlay_length),
                    ))
                end
            end
        end
    end

    ranking = isempty(ranking_rows) ? DataFrame() : DataFrame(ranking_rows)
    isempty(ranking) && error("No finite density-aware untreated fits were produced")
    for cell_line in unique(String.(ranking.cell_line))
        cell_mask = String.(ranking.cell_line) .== cell_line
        primary_mask = cell_mask .& Bool.(ranking.eligible_for_inheritance)
        primary_best = minimum(Float64.(ranking.bic[primary_mask]))
        for row_index in findall(cell_mask .& Bool.(ranking.diagnostic_model) .& (String.(ranking.pooling_mode) .!= "independent_diagnostic"))
            ranking.eligible_for_inheritance[row_index] = !ranking.boundary_issue[row_index] && Float64(ranking.bic[row_index]) <= primary_best - 10.0
        end
    end
    sort!(ranking, [:cell_line, :bic])
    pooling_summary = GrowthParameterEstimation.summarize_pooling_bic(ranking; top_n = 5)

    parameter_df = DataFrame(parameter_rows)
    status = pooling_summary.status
    baseline_rows = NamedTuple[]
    for status_row in eachrow(status)
        winner = ranking[
            (String.(ranking.cell_line) .== String(status_row.cell_line)) .&
            (String.(ranking.model) .== String(status_row.winning_model)) .&
            (String.(ranking.pooling_mode) .== String(status_row.winning_pooling_mode)),
            :,
        ]
        isempty(winner) && continue
        winner_parameters = parameter_df[
            (String.(parameter_df.cell_line) .== String(status_row.cell_line)) .&
            (String.(parameter_df.model) .== String(status_row.winning_model)) .&
            (String.(parameter_df.pooling_mode) .== String(status_row.winning_pooling_mode)),
            :,
        ]
        for density in unique(String.(winner_parameters.density))
            density_parameters = winner_parameters[String.(winner_parameters.density) .== density, :]
            r_row = density_parameters[String.(density_parameters.parameter) .== "r", :]
            k_row = density_parameters[String.(density_parameters.parameter) .== "K", :]
            shape_rows = density_parameters[.!in.(String.(density_parameters.parameter), Ref(["r", "K"])), :]
            push!(baseline_rows, (
                cell_line = String(status_row.cell_line),
                density = density,
                best_model = String(status_row.winning_model),
                pooling_mode = String(status_row.winning_pooling_mode),
                r = Float64(first(r_row.effective_value)),
                K = Float64(first(k_row.effective_value)),
                shape_parameter = isempty(shape_rows) ? "" : String(first(shape_rows.parameter)),
                shape_value = isempty(shape_rows) ? NaN : Float64(first(shape_rows.effective_value)),
                bic = Float64(first(winner.bic)),
                ssr = Float64(first(winner.ssr)),
                scaled_ssr = Float64(first(winner.scaled_ssr)),
                inheritance_allowed = Bool(status_row.inheritance_allowed),
                independent_bic_improvement = Float64(status_row.independent_bic_improvement),
            ))
        end
    end
    baselines = DataFrame(baseline_rows)
    if !isempty(initial_rows) && !isempty(baselines)
        initial_df = DataFrame(initial_rows)
        initial_df = leftjoin(initial_df, select(baselines, :cell_line, :density, :K); on = [:cell_line, :density])
        initial_df.u0_over_K = initial_df.fixed_u0 ./ initial_df.K
    else
        initial_df = DataFrame(initial_rows)
    end

    figures_csv = joinpath(out.csv, "figures")
    mkpath(figures_csv)
    CSV.write(joinpath(out.csv, "monoculture_untreated_pooling_model_ranking.csv"), ranking)
    CSV.write(joinpath(out.csv, "monoculture_untreated_pooling_top5.csv"), pooling_summary.ranking)
    CSV.write(joinpath(out.csv, "monoculture_untreated_pooling_status.csv"), status)
    CSV.write(joinpath(out.csv, "monoculture_untreated_pooling_parameter_estimates.csv"), parameter_df)
    CSV.write(joinpath(out.csv, "monoculture_untreated_initial_condition_diagnostics.csv"), initial_df)
    CSV.write(joinpath(out.csv, "untreated_group_baselines.csv"), baselines)
    CSV.write(joinpath(out.csv, "monoculture_untreated_boundary_profiles.csv"), isempty(profile_parts) ? DataFrame() : vcat(profile_parts...; cols = :union))
    CSV.write(joinpath(out.csv, "monoculture_untreated_identifiability.csv"), isempty(identifiability_parts) ? DataFrame() : vcat(identifiability_parts...; cols = :union))
    overlay_df = vcat(overlay_parts...; cols = :union)
    CSV.write(joinpath(figures_csv, "monoculture_untreated_pooling_overlays.csv"), overlay_df)
    _render_untreated_pooling_graph_grid(overlay_df, pooling_summary.ranking, out)

    ranking_path = joinpath(out.csv, "monoculture_untreated_automatic_model_ranking.csv")
    best_path = joinpath(out.csv, "monoculture_untreated_automatic_best_models_top10.csv")
    CSV.write(ranking_path, ranking)
    CSV.write(best_path, pooling_summary.ranking)
    summary_path = joinpath(out.csv, "best_model_summary.csv")
    CSV.write(summary_path, status)
    params_dir = joinpath(out.csv, "params")
    mkpath(params_dir)
    params_path = joinpath(params_dir, "best_params.csv")
    CSV.write(params_path, baselines)
    IOUtils.write_manifest_row(
        condition = "monoculture_untreated",
        step = "density_aware_fit",
        outputs = _nonmissing_strings([ranking_path, best_path, summary_path, params_path]),
        start = start,
    )
    inheritance_allowed = nrow(status) > 0 && all(Bool.(status.inheritance_allowed))
    return (
        result = nothing,
        ranking = ranking,
        ranking_path = ranking_path,
        best_path = best_path,
        untreated_baseline = inheritance_allowed ? (model = "density_aware", params = [mean(baselines.r), mean(baselines.K)], bic = minimum(ranking.bic), ssr = minimum(ranking.ssr), summary_path = summary_path, params_path = params_path) : nothing,
        staged_baseline_path = nothing,
        used_cached_results = false,
        fit_api_available = true,
        inheritance_allowed = inheritance_allowed,
    )
end

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

function run_condition_fit!(
    decoded::DataFrame,
    condition::AbstractString;
    start::AbstractString = pwd(),
    untreated_baseline = nothing,
    max_time_per_fit::Float64 = 45.0,
)
    ModelRegistry.ensure_model_registry!()
    out = IOUtils.condition_output_dirs(condition; start)

    staged_baseline_path = _write_staged_baseline(out, untreated_baseline)

    if condition == "monoculture_untreated"
        return _run_density_aware_untreated_fitting(
            decoded,
            out;
            start = start,
            max_time_per_fit = max_time_per_fit,
        )
    end

    # Monoculture treatment is fitted jointly across both densities and all doses.
    if condition == "monoculture_treated"
        ranking_df = _run_density_aware_treated_monoculture_fitting(
            decoded,
            out;
            start = start,
            max_time_per_fit = max_time_per_fit,
        )
        ranking_path = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
        best_path = joinpath(out.csv, "$(condition)_automatic_best_models_top10.csv")
        CSV.write(ranking_path, ranking_df)
        CSV.write(best_path, first(ranking_df, min(10, nrow(ranking_df))))
        IOUtils.write_manifest_row(
            condition = condition,
            step = "density_aware_joint_fit",
            outputs = _nonmissing_strings([ranking_path, best_path, staged_baseline_path]),
            start = start,
        )
        return (
            result = nothing,
            ranking = ranking_df,
            ranking_path = ranking_path,
            best_path = best_path,
            untreated_baseline = untreated_baseline,
            staged_baseline_path = staged_baseline_path,
            used_cached_results = false,
            fit_api_available = true,
        )
    end

    if condition in ("coculture_untreated", "coculture_treated")
        ranking_df = _run_density_aware_coculture_fitting(
            decoded,
            condition,
            out;
            start = start,
            max_time_per_fit = max_time_per_fit,
        )
        ranking_path = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
        best_path = joinpath(out.csv, "$(condition)_automatic_best_models_top10.csv")
        CSV.write(ranking_path, ranking_df)
        CSV.write(best_path, first(ranking_df, min(10, nrow(ranking_df))))
        IOUtils.write_manifest_row(
            condition = condition,
            step = "density_mix_aware_joint_fit",
            outputs = _nonmissing_strings([ranking_path, best_path, staged_baseline_path]),
            start = start,
        )
        return (
            result = nothing,
            ranking = ranking_df,
            ranking_path = ranking_path,
            best_path = best_path,
            untreated_baseline = untreated_baseline,
            staged_baseline_path = staged_baseline_path,
            used_cached_results = false,
            fit_api_available = true,
        )
    end

    # Legacy fallback retained for non-staged callers.
    if condition == "coculture_treated"
        ranking_df = _run_treated_fitting(decoded, out;
                                          start = start,
                                          untreated_baseline = untreated_baseline,
                                          max_time_per_fit = max_time_per_fit)
        ranking_path = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
        best_path    = joinpath(out.csv, "$(condition)_automatic_best_models_top10.csv")
        CSV.write(ranking_path, ranking_df)
        CSV.write(best_path, first(ranking_df, min(10, nrow(ranking_df))))

        IOUtils.write_manifest_row(
            condition = condition,
            step      = "fit",
            outputs   = _nonmissing_strings([ranking_path, best_path, staged_baseline_path]),
            start     = start,
        )
        return (
            result              = nothing,
            ranking             = ranking_df,
            ranking_path        = ranking_path,
            best_path           = best_path,
            untreated_baseline  = untreated_baseline,
            staged_baseline_path = staged_baseline_path,
            used_cached_results = false,
            fit_api_available   = true,
        )
    end

    # Untreated conditions: fit per experimental subgroup so r/K can carry over
    # to treated fits by cell line and density.
    model_map = Dict(
        "monoculture_untreated" => [
            "logistic_growth",
            "gompertz_growth",
            "logistic_simple_death",
            "allee_growth",
            "theta_logistic_growth",
            "generalized_logistic_growth",
        ],
        "coculture_untreated"   => ["logistic_growth"],
    )
    include_models = get(model_map, condition, ["logistic_growth"])

    # Detect columns
    time_col = :time
    :day in propertynames(decoded) && !(:time in propertynames(decoded)) && (time_col = :day)
    meta_cols = Set([:dose, :time, :day, :replicate, :cell_line, :density, :sample_id, :file, :source_file, :condition, :mix])
    numeric_cols = [c for c in propertynames(decoded) if !(c in meta_cols) && eltype(decoded[!, c]) <: Real]
    area_col = isempty(numeric_cols) ? error("No numeric measurement column") : numeric_cols[1]

    colset = Set(Symbol.(names(decoded)))
    cell_col = :cell_line in colset ? :cell_line : (:cellline in colset ? :cellline : nothing)
    density_col = :density in colset ? :density : nothing
    mix_col = :mix in colset ? :mix : nothing

    subgroup_cols = Symbol[]
    cell_col !== nothing && push!(subgroup_cols, cell_col)
    density_col !== nothing && push!(subgroup_cols, density_col)
    mix_col !== nothing && push!(subgroup_cols, mix_col)

    subgroups = isempty(subgroup_cols) ? [decoded] : collect(groupby(decoded, subgroup_cols))

    all_candidates = NamedTuple[]
    for (idx, subdf) in enumerate(subgroups)
        grouped = combine(groupby(subdf, time_col), area_col => mean => :y_mean)
        sort!(grouped, time_col)
        x = Float64.(grouped[!, time_col])
        y = Float64.(grouped[!, :y_mean])

        fits = Dict{String,NamedTuple}()
        for mname in include_models
            println("  Fitting $(mname) for untreated subgroup $(idx)...")
            fit = _fit_untreated_model(mname, x, y; max_time = min(max_time_per_fit, 8.0))
            fits[mname] = (bic = fit.bic, ssr = fit.ssr, params = fit.params)
        end

        subgroup_tmp = joinpath(out.csv, "tables", "subgroup_$(condition)_$(idx)_ranking.csv")
        mkpath(dirname(subgroup_tmp))
        subgroup_rank = DataFrame(
            model = String[k for k in keys(fits)],
            bic = Float64[v.bic for v in values(fits)],
            ssr = Float64[v.ssr for v in values(fits)],
            params = String[string(v.params) for v in values(fits)],
        )
        sort!(subgroup_rank, :bic)
        CSV.write(subgroup_tmp, subgroup_rank)

        cell_label = cell_col === nothing ? "pooled" : _safe_string(first(subdf[!, cell_col]); default = "pooled")
        den_label = density_col === nothing ? "pooled" : _safe_string(first(subdf[!, density_col]); default = "pooled")
        mix_label = mix_col === nothing ? "pooled" : _safe_string(first(subdf[!, mix_col]); default = "pooled")

        for (mname, v) in pairs(fits)
            push!(all_candidates, (
                model = String(mname),
                bic = Float64(v.bic),
                ssr = Float64(v.ssr),
                params_vec = Vector{Float64}(v.params),
                params = string(v.params),
                cell_line = cell_label,
                density = den_label,
                mix = mix_label,
            ))
        end
    end

    ranking_df = sort!(DataFrame(
        model = [c.model for c in all_candidates],
        bic = [c.bic for c in all_candidates],
        ssr = [c.ssr for c in all_candidates],
        params = [c.params for c in all_candidates],
        cell_line = [c.cell_line for c in all_candidates],
        density = [c.density for c in all_candidates],
        mix = [c.mix for c in all_candidates],
    ), :bic)

    # Build per-subgroup untreated baseline table used by treated models.
    subgroup_baselines = NamedTuple[]
    grouped_candidates = groupby(ranking_df, [:cell_line, :density])
    for g in grouped_candidates
        g_sorted = sort(g, :bic)
        best = first(g_sorted, 1)
        params_txt = String(best.params[1])
        # Parse out first two numeric entries from the parameter vector string.
        nums = [parse(Float64, m.match) for m in eachmatch(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", params_txt)]
        if length(nums) >= 2
            push!(subgroup_baselines, (
                cell_line = String(best.cell_line[1]),
                density = String(best.density[1]),
                best_model = String(best.model[1]),
                r = nums[1],
                K = nums[2],
                bic = Float64(best.bic[1]),
                ssr = Float64(best.ssr[1]),
            ))
        end
    end

    subgroup_baseline_df = DataFrame(subgroup_baselines)
    subgroup_baseline_path = joinpath(out.csv, "untreated_group_baselines.csv")
    CSV.write(subgroup_baseline_path, subgroup_baseline_df)

    if !isempty(all_candidates)
        best_idx = argmin([c.bic for c in all_candidates])
        best = all_candidates[best_idx]
        summary_path = joinpath(out.csv, "best_model_summary.csv")
        params_dir = joinpath(out.csv, "params")
        mkpath(params_dir)
        params_path = joinpath(params_dir, "best_params.csv")
        CSV.write(summary_path, DataFrame(key = ["best_model", "best_cell_line", "best_density"], value = [best.model, best.cell_line, best.density]))
        CSV.write(params_path, DataFrame(param = ["p$(i)" for i in eachindex(best.params_vec)], value = best.params_vec))
    end

    ranking_path = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
    best_path    = joinpath(out.csv, "$(condition)_automatic_best_models_top10.csv")
    CSV.write(ranking_path, ranking_df)
    CSV.write(best_path, first(ranking_df, min(10, nrow(ranking_df))))

    IOUtils.write_manifest_row(
        condition = condition,
        step      = "fit",
        outputs   = _nonmissing_strings([ranking_path, best_path, subgroup_baseline_path, staged_baseline_path]),
        start     = start,
    )
    return (
        result              = nothing,
        ranking             = ranking_df,
        ranking_path        = ranking_path,
        best_path           = best_path,
        subgroup_baseline_path = subgroup_baseline_path,
        untreated_baseline  = untreated_baseline,
        staged_baseline_path = staged_baseline_path,
        used_cached_results = false,
        fit_api_available   = true,
    )
end

include("coculture_joint.jl")
include("linked_treatment_joint.jl")
include("timing_hypotheses.jl")

end
