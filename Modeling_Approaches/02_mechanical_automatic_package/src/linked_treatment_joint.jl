# Linked treated-monoculture and treated-coculture fitting. Included in FitWorkflows.
#
# Scientific policy:
# - fit every linked-treatment hypothesis below;
# - keep the full ranking for audit and downstream reuse;
# - show only the top five hypotheses in report-facing summary artifacts.

using SHA

const LINKED_TREATMENT_VISIBLE_LIMIT = 5

const LINKED_TREATMENT_HYPOTHESES = [
    "strict_inheritance",
    "competitor_scaled",
    "load_scaled",
    "load_plus_context_amplitude",
    "tolerant_context_shift",
    "subpopulation_load_scaled",
    "load_plus_tolerant_context",
    "load_plus_tolerant_growth_context",
    "fully_free_context_diagnostic",
]

const LINKED_TREATMENT_CANONICAL_DATA_COMMIT = "41ee3c9d281b79634c441d481db4f7a8092d0b1d"

function _linked_treatment_seed(start)
    root = IOUtils.package_root(start)
    base = joinpath(root, "outputs", "csv", "monoculture_treated")
    status_path = joinpath(base, "monoculture_treated_pooling_status.csv")
    parameter_path = joinpath(base, "monoculture_treated_joint_parameter_estimates.csv")
    ranking_path = joinpath(base, "monoculture_treated_pooling_model_ranking.csv")
    isfile(status_path) || error("Treated-monoculture pooling status is required for linked treatment fitting")
    isfile(parameter_path) || error("Treated-monoculture parameter estimates are required for linked treatment fitting")
    isfile(ranking_path) || error("Treated-monoculture model ranking is required for linked treatment fitting")
    status = CSV.read(status_path, DataFrame)
    parameters = CSV.read(parameter_path, DataFrame)
    ranking = CSV.read(ranking_path, DataFrame)

    expected = Dict(
        "A2780Naive" => "joint_ic_effect_hill_ramp_onset",
        "A2780cis" => "joint_ic_effect_two_population",
    )
    values = Dict{Tuple{String,Symbol},Float64}()
    rows = NamedTuple[]
    for (cell_line, expected_model) in expected
        selected = status[String.(status.cell_line) .== cell_line, :]
        nrow(selected) == 1 || error("Expected one treated-monoculture winner for $(cell_line)")
        winner = first(selected)
        nominal_model = String(winner.winning_model)
        compatible = ranking[
            (String.(ranking.cell_line) .== cell_line) .&
            (String.(ranking.model) .== expected_model) .&
            Bool.(ranking.eligible_for_inheritance),
            :,
        ]
        isempty(compatible) && error("No eligible $(expected_model) fit is available for linked $(cell_line) treatment")
        sort!(compatible, :bic)
        compatible_winner = first(compatible)
        selected_pooling = String(compatible_winner.pooling_mode)
        nominal_bic = Float64(winner.winning_bic)
        compatible_bic = Float64(compatible_winner.bic)
        selection_reason = nominal_model == expected_model ?
            "nominal_stage2_winner" :
            "best_stage2_candidate_compatible_with_population_balance_stage4"
        winner_parameters = parameters[
            (String.(parameters.cell_line) .== cell_line) .&
            (String.(parameters.model) .== expected_model) .&
            (String.(parameters.pooling_mode) .== selected_pooling),
            :,
        ]
        for parameter in unique(String.(winner_parameters.parameter))
            parameter_rows = winner_parameters[String.(winner_parameters.parameter) .== parameter, :]
            center = Float64(first(parameter_rows.center_value))
            values[(cell_line, Symbol(parameter))] = center
            push!(rows, (
                cell_line = cell_line,
                treatment_family = expected_model,
                pooling_mode = selected_pooling,
                parameter = parameter,
                sequential_stage2_center = center,
                nominal_stage2_winner = nominal_model,
                nominal_stage2_bic = nominal_bic,
                compatible_stage2_bic = compatible_bic,
                compatible_delta_bic = compatible_bic - nominal_bic,
                selection_reason = selection_reason,
                source_path = parameter_path,
            ))
        end
        contrast_rows = parameters[
            (String.(parameters.cell_line) .== cell_line) .&
            (String.(parameters.model) .== expected_model) .&
            (String.(parameters.pooling_mode) .== selected_pooling) .&
            (String.(parameters.parameter) .== (cell_line == "A2780Naive" ? "emax" : "emax_sensitive")),
            :,
        ]
        low = Float64(first(contrast_rows.effective_value))
        high = Float64(last(contrast_rows.effective_value))
        values[(cell_line, :log_contrast_effect)] = 0.5 * log(max(high, 1e-12) / max(low, 1e-12))
    end
    timing_ranking_path = joinpath(base, "monoculture_treated_timing_hypothesis_ranking.csv")
    timing_parameter_path = joinpath(base, "monoculture_treated_timing_hypothesis_parameters.csv")
    timing_spec = nothing
    if isfile(timing_ranking_path) && isfile(timing_parameter_path)
        timing_ranking = sort!(CSV.read(timing_ranking_path, DataFrame), :bic)
        nrow(timing_ranking) > 0 || error("Treated timing ranking is empty")
        timing_winner = String(first(timing_ranking.model))
        timing_parameters = CSV.read(timing_parameter_path, DataFrame)
        selected = timing_parameters[String.(timing_parameters.timing_hypothesis) .== timing_winner, :]
        nrow(selected) > 0 || error("No parameters found for winning timing hypothesis $(timing_winner)")
        names = Symbol.(String.(selected.parameter))
        p0 = Float64.(selected.estimate)
        bounds = collect(zip(Float64.(selected.lower_bound), Float64.(selected.upper_bound)))
        all(isfinite, p0) || error("Winning timing hypothesis contains non-finite parameters")
        timing_spec = (
            hypothesis = timing_winner,
            names = names,
            p0 = p0,
            bounds = bounds,
            index = Dict(name => index for (index, name) in enumerate(names)),
        )
        for (name, estimate) in zip(names, p0)
            push!(rows, (
                cell_line = "joint treated timing",
                treatment_family = timing_winner,
                pooling_mode = "joint across both cell lines, densities, and doses",
                parameter = String(name),
                sequential_stage2_center = estimate,
                nominal_stage2_winner = timing_winner,
                nominal_stage2_bic = Float64(first(timing_ranking.bic)),
                compatible_stage2_bic = Float64(first(timing_ranking.bic)),
                compatible_delta_bic = 0.0,
                selection_reason = "timing_audit_winner",
                source_path = timing_parameter_path,
            ))
        end
    end
    return (
        values = values,
        audit = DataFrame(rows),
        status_path = status_path,
        parameter_path = parameter_path,
        timing_spec = timing_spec,
        timing_ranking_path = timing_ranking_path,
        timing_parameter_path = timing_parameter_path,
    )
end

function _linked_monoculture_environments(start)
    root = IOUtils.package_root(start)
    path = joinpath(root, "outputs", "csv", "monoculture_treated", "monoculture_treated_a2780_decoded.csv")
    isfile(path) || error("Decoded treated-monoculture data are required for linked fitting: $(path)")
    decoded = CSV.read(path, DataFrame)
    environments = NamedTuple[]
    for cell_line in ("A2780Naive", "A2780cis")
        cell_rows = decoded[String.(decoded.cell_line) .== cell_line, :]
        for density in ("20k", "30k")
            density_rows = cell_rows[String.(cell_rows.density) .== density, :]
            for dose in sort(unique(Float64.(density_rows.dose)))
                trajectory = density_rows[isapprox.(Float64.(density_rows.dose), dose; atol = 1e-10), :]
                observed = combine(groupby(trajectory, :time), :count => mean => :observed)
                sort!(observed, :time)
                metadata = _treated_dose_metadata(dose)
                sensitive_baseline, resistant_baseline = _coculture_monoculture_baselines(start, density)
                push!(environments, (
                    context = "monoculture",
                    cell_line = cell_line,
                    density = density,
                    dose = dose,
                    effect_level = Float64(metadata.effect_level),
                    ic_label = String(metadata.ic_label),
                    times = Float64.(observed.time),
                    observed = Float64.(observed.observed),
                    residual_scale = max(maximum(Float64.(observed.observed)), 1.0),
                    fixed_u0 = _fixed_day0_total(density),
                    baseline = cell_line == "A2780Naive" ? sensitive_baseline : resistant_baseline,
                ))
            end
        end
    end
    length(environments) == 12 || error("Expected 12 treated-monoculture trajectories, found $(length(environments))")
    return environments, path
end

function _linked_base_spec(seed)
    if seed.timing_spec !== nothing
        timing = seed.timing_spec
        inheritance_bounds = [
            (
                max(first(bound), value - max(0.05 * abs(value), 1e-4)),
                min(last(bound), value + max(0.05 * abs(value), 1e-4)),
            ) for (value, bound) in zip(timing.p0, timing.bounds)
        ]
        return (
            p0 = copy(timing.p0),
            bounds = inheritance_bounds,
            original_bounds = copy(timing.bounds),
            names = copy(timing.names),
            timing_spec = timing,
            nbase = length(timing.p0),
        )
    end
    values = seed.values
    p0 = [
        values[("A2780Naive", :emax)],
        values[("A2780Naive", :ec50_effect)],
        values[("A2780Naive", :hill_n)],
        values[("A2780Naive", :lambda)],
        values[("A2780Naive", :t_onset)],
        values[("A2780Naive", :log_contrast_effect)],
        values[("A2780cis", :emax_sensitive)],
        values[("A2780cis", :emax_tolerant)],
        values[("A2780cis", :lambda)],
        values[("A2780cis", :t_onset)],
        values[("A2780cis", :f_tolerant0)],
        values[("A2780cis", :log_contrast_effect)],
    ]
    bounds = [
        (0.0, 5.0), (0.05, 2.0), (0.2, 12.0), (0.01, 5.0), (0.0, 7.0), (-DENSITY_LOG_CONTRAST_BOUND, DENSITY_LOG_CONTRAST_BOUND),
        (0.0, 5.0), (0.0, 3.0), (0.01, 5.0), (0.0, 7.0), (0.0, 0.95), (-DENSITY_LOG_CONTRAST_BOUND, DENSITY_LOG_CONTRAST_BOUND),
    ]
    names = [
        :naive_emax, :naive_ec50_effect, :naive_hill_n, :naive_lambda, :naive_t_onset, :naive_log_contrast_density,
        :cis_emax_sensitive, :cis_emax_tolerant, :cis_lambda, :cis_t_onset, :cis_f_tolerant0, :cis_log_contrast_density,
    ]
    timing_spec = (
        hypothesis = "independent_onset_gradual",
        names = names,
        p0 = p0,
        bounds = bounds,
        index = Dict(name => index for (index, name) in enumerate(names)),
    )
    inheritance_bounds = [
        (
            max(first(bound), value - max(0.05 * abs(value), 1e-4)),
            min(last(bound), value + max(0.05 * abs(value), 1e-4)),
        ) for (value, bound) in zip(p0, bounds)
    ]
    return (
        p0 = p0,
        bounds = inheritance_bounds,
        original_bounds = bounds,
        names = names,
        timing_spec = timing_spec,
        nbase = length(p0),
    )
end

function _linked_hypothesis_spec(base, hypothesis)
    if hypothesis == "strict_inheritance"
        return (p0 = copy(base.p0), bounds = copy(base.bounds), names = copy(base.names))
    elseif hypothesis in ("competitor_scaled", "load_scaled")
        return (
            p0 = vcat(base.p0, [0.0, 0.0]),
            bounds = vcat(base.bounds, [(-3.0, 3.0), (-3.0, 3.0)]),
            names = vcat(base.names, [:beta_naive, :beta_cis]),
        )
    elseif hypothesis == "load_plus_context_amplitude"
        context_bound = log(1.20)
        return (
            p0 = vcat(base.p0, [0.0, 0.0, 0.0, 0.0]),
            bounds = vcat(base.bounds, [(-3.0, 3.0), (-3.0, 3.0), (-context_bound, context_bound), (-context_bound, context_bound)]),
            names = vcat(base.names, [:beta_naive, :beta_cis, :log_context_amplitude_naive, :log_context_amplitude_cis]),
        )
    elseif hypothesis == "tolerant_context_shift"
        return (
            p0 = vcat(base.p0, [-1.0, 0.5]),
            bounds = vcat(base.bounds, [(-4.0, 4.0), (-3.0, 3.0)]),
            names = vcat(base.names, [:log_cis_tolerant_context, :cis_tolerant_logit_shift]),
        )
    elseif hypothesis == "subpopulation_load_scaled"
        return (
            p0 = vcat(base.p0, [0.0, 0.0, 0.0]),
            bounds = vcat(base.bounds, [(-3.0, 3.0), (-3.0, 3.0), (-3.0, 3.0)]),
            names = vcat(base.names, [:beta_naive, :beta_cis_sensitive, :beta_cis_tolerant]),
        )
    elseif hypothesis == "load_plus_tolerant_context"
        return (
            p0 = vcat(base.p0, [0.0, 0.0, -1.0, 0.5]),
            bounds = vcat(base.bounds, [(-3.0, 3.0), (-3.0, 3.0), (-4.0, 4.0), (-3.0, 3.0)]),
            names = vcat(base.names, [:beta_naive, :beta_cis, :log_cis_tolerant_context, :cis_tolerant_logit_shift]),
        )
    elseif hypothesis == "load_plus_tolerant_growth_context"
        return (
            p0 = vcat(base.p0, [0.0, 0.0, -1.0, 0.5, 0.0]),
            bounds = vcat(base.bounds, [(-3.0, 3.0), (-3.0, 3.0), (-4.0, 4.0), (-3.0, 3.0), (-log(2.0), log(2.0))]),
            names = vcat(base.names, [:beta_naive, :beta_cis, :log_cis_tolerant_context, :cis_tolerant_logit_shift, :log_cis_tolerant_growth_context]),
        )
    elseif hypothesis == "fully_free_context_diagnostic"
        return (
            p0 = vcat(base.p0, base.p0),
            bounds = vcat(base.bounds, base.original_bounds),
            names = vcat(base.names, Symbol.("coculture_" .* string.(base.names))),
        )
    end
    error("Unknown linked treatment hypothesis $(hypothesis)")
end

function _linked_parameter_block(p, hypothesis, context)
    error("_linked_parameter_block requires the linked base specification")
end

function _linked_parameter_block(p, hypothesis, context, base)
    if context == :coculture && hypothesis == "fully_free_context_diagnostic"
        return @view p[(base.nbase + 1):(2base.nbase)]
    end
    return @view p[1:base.nbase]
end

function _linked_treatment_values(p, hypothesis, context, density, base)
    q = _linked_parameter_block(p, hypothesis, context, base)
    treatment = _timing_values(q, base.timing_spec, density)
    context == :coculture || return treatment
    if hypothesis == "tolerant_context_shift"
        gamma, shift = p[base.nbase + 1], p[base.nbase + 2]
    elseif hypothesis in ("load_plus_tolerant_context", "load_plus_tolerant_growth_context")
        gamma, shift = p[base.nbase + 3], p[base.nbase + 4]
    else
        return treatment
    end
    f = clamp(treatment.cis_f_tolerant0, 1e-8, 1 - 1e-8)
    shifted_fraction = inv(one(f) + exp(-(log(f / (1 - f)) + shift)))
    return merge(treatment, (
        cis_emax_tolerant = treatment.cis_emax_tolerant * exp(clamp(gamma, -4.0, 4.0)),
        cis_f_tolerant0 = shifted_fraction,
    ))
end

function _linked_tolerant_growth_scale(p, hypothesis, context, base)
    context == :coculture && hypothesis == "load_plus_tolerant_growth_context" || return one(eltype(p))
    return exp(clamp(p[base.nbase + 5], -log(4.0), log(4.0)))
end

function _linked_context_scales(p, hypothesis, sensitive_burden, resistant_burden, base)
    hypothesis in ("strict_inheritance", "tolerant_context_shift", "fully_free_context_diagnostic") &&
        return (one(sensitive_burden), one(resistant_burden), one(resistant_burden))
    if hypothesis == "subpopulation_load_scaled"
        beta_naive = p[base.nbase + 1]
        beta_sensitive = p[base.nbase + 2]
        beta_tolerant = p[base.nbase + 3]
        return (
            exp(clamp(beta_naive * sensitive_burden, -4.0, 4.0)),
            exp(clamp(beta_sensitive * resistant_burden, -4.0, 4.0)),
            exp(clamp(beta_tolerant * resistant_burden, -4.0, 4.0)),
        )
    end
    beta_naive, beta_cis = p[base.nbase + 1], p[base.nbase + 2]
    naive_log_scale = beta_naive * sensitive_burden
    cis_log_scale = beta_cis * resistant_burden
    if hypothesis == "load_plus_context_amplitude"
        naive_log_scale += p[base.nbase + 3]
        cis_log_scale += p[base.nbase + 4]
    end
    return (
        exp(clamp(naive_log_scale, -4.0, 4.0)),
        exp(clamp(cis_log_scale, -4.0, 4.0)),
        exp(clamp(cis_log_scale, -4.0, 4.0)),
    )
end

function _linked_intrinsic_growth(population, baseline, t)
    return _baseline_growth(population, population, baseline, t)
end

function _linked_problem(monoculture_environments, coculture_environments, start, seed, hypothesis)
    base = _linked_base_spec(seed)
    spec = _linked_hypothesis_spec(base, hypothesis)
    untreated_baseline = _load_coculture_untreated_baseline(; start = start)
    coculture_baselines = [_coculture_monoculture_baselines(start, environment.density) for environment in coculture_environments]
    untreated_specs = _untreated_coculture_specs(coculture_environments, coculture_baselines, untreated_baseline.pooling_mode)
    untreated_spec = untreated_specs[untreated_baseline.model]
    interaction_values = [untreated_spec.local_parameters(untreated_baseline.params, index) for index in eachindex(coculture_environments)]

    naive_mono = [environment for environment in monoculture_environments if environment.cell_line == "A2780Naive"]
    cis_mono = [environment for environment in monoculture_environments if environment.cell_line == "A2780cis"]
    n_naive = length(naive_mono)
    n_cis = length(cis_mono)
    mono_cis_offset = n_naive
    coculture_offset = n_naive + 2n_cis

    model! = function (du, u, p, t)
        for (index, environment) in enumerate(naive_mono)
            N = max(u[index], zero(u[index]))
            treatment = _linked_treatment_values(p, hypothesis, :monoculture, environment.density, base)
            activation = _ramp_activation(t, 0.0, treatment.naive_lambda, treatment.naive_onset)
            kill = activation * _hill_effect(environment.effect_level, treatment.naive_emax, treatment.naive_ec50, treatment.naive_hill)
            du[index] = _linked_intrinsic_growth(N, environment.baseline, t) - kill * N
        end
        for (index, environment) in enumerate(cis_mono)
            sensitive_index = mono_cis_offset + 2index - 1
            tolerant_index = sensitive_index + 1
            sensitive = max(u[sensitive_index], zero(u[sensitive_index]))
            tolerant = max(u[tolerant_index], zero(u[tolerant_index]))
            total = sensitive + tolerant
            treatment = _linked_treatment_values(p, hypothesis, :monoculture, environment.density, base)
            activation = _timing_activation(t, treatment.cis_lambda, treatment.cis_onset, treatment.cis_activation_mode)
            kill_sensitive = activation * _hill_effect(environment.effect_level, treatment.cis_emax_sensitive, 0.5, 4.0)
            kill_tolerant = activation * _hill_effect(environment.effect_level, treatment.cis_emax_tolerant, 0.5, 4.0)
            growth = _linked_intrinsic_growth(total, environment.baseline, t)
            sensitive_share = total > 0 ? sensitive / total : one(total)
            tolerant_share = total > 0 ? tolerant / total : zero(total)
            du[sensitive_index] = growth * sensitive_share - kill_sensitive * sensitive
            du[tolerant_index] = growth * tolerant_share - kill_tolerant * tolerant
        end
        for (index, environment) in enumerate(coculture_environments)
            naive_index = coculture_offset + 3index - 2
            cis_sensitive_index = naive_index + 1
            cis_tolerant_index = naive_index + 2
            naive = max(u[naive_index], zero(u[naive_index]))
            cis_sensitive = max(u[cis_sensitive_index], zero(u[cis_sensitive_index]))
            cis_tolerant = max(u[cis_tolerant_index], zero(u[cis_tolerant_index]))
            cis_total = cis_sensitive + cis_tolerant
            local_interactions = interaction_values[index]
            if untreated_baseline.model == "lv_symmetric_competition"
                alpha_sr, alpha_rs, death_naive, death_cis = local_interactions[1], local_interactions[1], 0.0, 0.0
            elseif untreated_baseline.model == "lv_asymmetric_competition"
                alpha_sr, alpha_rs, death_naive, death_cis = local_interactions[1], local_interactions[2], 0.0, 0.0
            else
                alpha_sr, alpha_rs, death_naive, death_cis = local_interactions
            end
            naive_baseline, cis_baseline = coculture_baselines[index]
            naive_load = naive + alpha_sr * cis_total
            cis_load = cis_total + alpha_rs * naive
            naive_growth = _anchored_component_growth(naive, naive_load, naive_baseline, t) - death_naive * naive
            cis_growth = _anchored_component_growth(cis_total, cis_load, cis_baseline, t) - death_cis * cis_total
            treatment = _linked_treatment_values(p, hypothesis, :coculture, environment.density, base)
            if hypothesis == "competitor_scaled"
                naive_burden = alpha_sr * cis_total / max(naive_baseline.K, 1e-8)
                cis_burden = alpha_rs * naive / max(cis_baseline.K, 1e-8)
            else
                naive_burden = naive_load / max(naive_baseline.K, 1e-8)
                cis_burden = cis_load / max(cis_baseline.K, 1e-8)
            end
            naive_scale, cis_sensitive_scale, cis_tolerant_scale =
                _linked_context_scales(p, hypothesis, naive_burden, cis_burden, base)
            naive_activation = _ramp_activation(t, 0.0, treatment.naive_lambda, treatment.naive_onset)
            cis_activation = _timing_activation(t, treatment.cis_lambda, treatment.cis_onset, treatment.cis_activation_mode)
            naive_kill = naive_scale * naive_activation * _hill_effect(0.5, treatment.naive_emax, treatment.naive_ec50, treatment.naive_hill)
            cis_kill_sensitive = cis_sensitive_scale * cis_activation * _hill_effect(0.5, treatment.cis_emax_sensitive, 0.5, 4.0)
            cis_kill_tolerant = cis_tolerant_scale * cis_activation * _hill_effect(0.5, treatment.cis_emax_tolerant, 0.5, 4.0)
            cis_sensitive_share = cis_total > 0 ? cis_sensitive / cis_total : one(cis_total)
            cis_tolerant_share = cis_total > 0 ? cis_tolerant / cis_total : zero(cis_total)
            tolerant_growth_scale = _linked_tolerant_growth_scale(p, hypothesis, :coculture, base)
            du[naive_index] = naive_growth - naive_kill * naive
            du[cis_sensitive_index] = cis_growth * cis_sensitive_share - cis_kill_sensitive * cis_sensitive
            du[cis_tolerant_index] = tolerant_growth_scale * cis_growth * cis_tolerant_share - cis_kill_tolerant * cis_tolerant
        end
    end

    u0_builder = function (p)
        values = eltype(p)[]
        for environment in naive_mono
            push!(values, environment.fixed_u0)
        end
        for environment in cis_mono
            treatment = _linked_treatment_values(p, hypothesis, :monoculture, environment.density, base)
            push!(values, environment.fixed_u0 * (1 - treatment.cis_f_tolerant0), environment.fixed_u0 * treatment.cis_f_tolerant0)
        end
        for environment in coculture_environments
            treatment = _linked_treatment_values(p, hypothesis, :coculture, environment.density, base)
            push!(values,
                environment.u0_sensitive,
                environment.u0_resistant * (1 - treatment.cis_f_tolerant0),
                environment.u0_resistant * treatment.cis_f_tolerant0,
            )
        end
        return values
    end

    datasets = NamedTuple[]
    metadata = NamedTuple[]
    for (index, environment) in enumerate(naive_mono)
        push!(datasets, (x = environment.times, y = environment.observed, state_index = index, residual_scale = environment.residual_scale))
        push!(metadata, (context = "monoculture", cell_line = environment.cell_line, density = environment.density, dose = environment.dose, mix = "", component = "total", fixed_u0 = environment.fixed_u0))
    end
    for (index, environment) in enumerate(cis_mono)
        sensitive_index = mono_cis_offset + 2index - 1
        tolerant_index = sensitive_index + 1
        observable = let sensitive_index = sensitive_index, tolerant_index = tolerant_index
            (u, p, t) -> u[sensitive_index] + u[tolerant_index]
        end
        push!(datasets, (x = environment.times, y = environment.observed, observable = observable, residual_scale = environment.residual_scale))
        push!(metadata, (context = "monoculture", cell_line = environment.cell_line, density = environment.density, dose = environment.dose, mix = "", component = "total", fixed_u0 = environment.fixed_u0))
    end
    for (index, environment) in enumerate(coculture_environments)
        naive_index = coculture_offset + 3index - 2
        cis_sensitive_index = naive_index + 1
        cis_tolerant_index = naive_index + 2
        resistant_observable = let cis_sensitive_index = cis_sensitive_index, cis_tolerant_index = cis_tolerant_index
            (u, p, t) -> u[cis_sensitive_index] + u[cis_tolerant_index]
        end
        push!(datasets, (x = environment.times, y = environment.sensitive, state_index = naive_index, residual_scale = environment.sensitive_scale))
        push!(metadata, (context = "coculture", cell_line = COCULTURE_SYSTEM_LABEL, density = environment.density, dose = 1.0, mix = environment.mix, component = "sensitive", fixed_u0 = environment.u0_sensitive))
        push!(datasets, (x = environment.times, y = environment.resistant, observable = resistant_observable, residual_scale = environment.resistant_scale))
        push!(metadata, (context = "coculture", cell_line = COCULTURE_SYSTEM_LABEL, density = environment.density, dose = 1.0, mix = environment.mix, component = "resistant", fixed_u0 = environment.u0_resistant))
    end
    return merge(spec, (
        model = model!,
        datasets = datasets,
        metadata = metadata,
        u0_builder = u0_builder,
        u0 = Float64.(u0_builder(spec.p0)),
        untreated_baseline = untreated_baseline,
        base = base,
        timing_hypothesis = base.timing_spec.hypothesis,
    ))
end

function _linked_overlay(fit, problem, hypothesis)
    parts = DataFrame[]
    for (index, dataset) in enumerate(problem.datasets)
        meta = problem.metadata[index]
        n = length(dataset.x) + 1
        push!(parts, DataFrame(
            time = vcat(0.0, dataset.x),
            observed = vcat(meta.fixed_u0, dataset.y),
            predicted = vcat(meta.fixed_u0, fit.predictions[index]),
            fixed_day0_anchor = vcat(true, fill(false, length(dataset.x))),
            model = fill(hypothesis, n),
            pooling_mode = fill("linked_global", n),
            context = fill(meta.context, n),
            cell_line = fill(meta.cell_line, n),
            density = fill(meta.density, n),
            dose = fill(meta.dose, n),
            mix = fill(meta.mix, n),
            component = fill(meta.component, n),
        ))
    end
    return vcat(parts...; cols = :union)
end

function _render_linked_monoculture_grid(overlay, winner, out)
    selected = overlay[(String.(overlay.model) .== String(winner.model)) .& (String.(overlay.context) .== "monoculture"), :]
    panels = Any[]
    for cell_line in ("A2780Naive", "A2780cis"), density in ("20k", "30k")
        rows = selected[(String.(selected.cell_line) .== cell_line) .& (String.(selected.density) .== density), :]
        for dose_group in sort(collect(groupby(rows, :dose)); by = group -> Float64(first(group.dose)))
            sort!(dose_group, :time)
            panel = plot(title = "$(cell_line) $(density)\n$(_treated_dose_metadata(Float64(first(dose_group.dose))).ic_label) ($(first(dose_group.dose)) uM)", titlefontsize = 8, xlabel = "Time (day)", ylabel = "Cell count", legend = false)
            scatter!(panel, dose_group.time, dose_group.observed; color = :black, ms = 3.0, markerstrokewidth = 0, label = "")
            plot!(panel, dose_group.time, dose_group.predicted; color = :steelblue, lw = 2.4, label = "")
            push!(panels, panel)
        end
    end
    figure = plot(panels...; layout = (4, 3), size = (1500, 1500), plot_title = "Linked global fit: treated monoculture", margin = 4 * Plots.mm)
    path = joinpath(out.images, "figures", "linked_treatment_monoculture_grid.png")
    mkpath(dirname(path)); savefig(figure, path); return path
end

function _render_linked_coculture_grid(overlay, winner, out)
    selected = overlay[(String.(overlay.model) .== String(winner.model)) .& (String.(overlay.context) .== "coculture"), :]
    panels = Any[]
    y_limits = _shared_y_limits(selected)
    for density in ("20k", "30k"), mix in ("25-75", "50-50", "75-25")
        environment = selected[(String.(selected.density) .== density) .& (String.(selected.mix) .== mix), :]
        panel = plot(title = "$(density), mix $(mix)", titlefontsize = 9, xlabel = "Time (day)", ylabel = mix == "25-75" ? "Measured population" : "", ylims = y_limits, legend = density == "20k" && mix == "75-25" ? :outertopright : false, legendfontsize = 8)
        for (component, color) in (("sensitive", :crimson), ("resistant", :steelblue))
            rows = environment[String.(environment.component) .== component, :]
            sort!(rows, :time)
            scatter!(panel, rows.time, rows.observed; color = color, ms = 3.2, markerstrokewidth = 0, label = "$(component) observed")
            plot!(panel, rows.time, rows.predicted; color = color, lw = 2.4, label = "$(component) linked fit")
        end
        push!(panels, panel)
    end
    figure = plot(panels...; layout = (2, 3), size = (1500, 820), plot_title = "Linked treatment inheritance: $(winner.model), BIC=$(round(Float64(winner.bic); digits = 1))", margin = 4 * Plots.mm)
    path = joinpath(out.images, "figures", "linked_treatment_coculture_grid.png")
    mkpath(dirname(path)); savefig(figure, path); return path
end

function _render_linked_hypothesis_grid(overlay, winner, out)
    hypotheses = unique(["strict_inheritance", String(winner.model)])
    styles = Dict("strict_inheritance" => :dash, String(winner.model) => :solid)
    panels = Any[]
    for density in ("20k", "30k"), mix in ("25-75", "50-50", "75-25")
        panel = plot(title = "$(density), mix $(mix)", titlefontsize = 9, xlabel = "Time (day)", ylabel = mix == "25-75" ? "Measured population" : "", legend = density == "20k" && mix == "75-25" ? :outertopright : false, legendfontsize = 7)
        reference = overlay[(String.(overlay.model) .== hypotheses[1]) .& (String.(overlay.context) .== "coculture") .& (String.(overlay.density) .== density) .& (String.(overlay.mix) .== mix), :]
        for (component, color, short) in (("sensitive", :crimson, "S"), ("resistant", :steelblue, "R"))
            observed = reference[String.(reference.component) .== component, :]
            sort!(observed, :time)
            scatter!(panel, observed.time, observed.observed; color = color, ms = 3.0, markerstrokewidth = 0, alpha = 0.7, label = "$(short) observed")
            for hypothesis in hypotheses
                rows = overlay[(String.(overlay.model) .== hypothesis) .& (String.(overlay.context) .== "coculture") .& (String.(overlay.density) .== density) .& (String.(overlay.mix) .== mix) .& (String.(overlay.component) .== component), :]
                sort!(rows, :time)
                plot!(panel, rows.time, rows.predicted; color = color, linestyle = styles[hypothesis], lw = 2.2, label = "$(short) $(replace(hypothesis, "_" => " "))")
            end
        end
        push!(panels, panel)
    end
    figure = plot(panels...; layout = (2, 3), size = (1600, 860), plot_title = "Strict drug inheritance versus selected coculture modifier", margin = 4 * Plots.mm)
    path = joinpath(out.images, "figures", "linked_treatment_hypothesis_comparison_grid.png")
    mkpath(dirname(path)); savefig(figure, path); return path
end

function _write_linked_provenance(paths, out; start = pwd())
    canonical_root = joinpath(IOUtils.find_repo_root(start), "Processed_Datasets")
    canonical_present = isdir(canonical_root)
    source_mode = canonical_present ?
        "canonical Processed_Datasets restored from Git history and decoded for this run" :
        "validated cached decoded artifact; canonical Processed_Datasets directory absent"
    rows = [(
        artifact = basename(path),
        path = path,
        exists = isfile(path),
        bytes = isfile(path) ? filesize(path) : 0,
        sha256 = isfile(path) ? bytes2hex(sha256(read(path))) : "",
        source_mode = source_mode,
        canonical_processed_data_present = canonical_present,
        canonical_source_commit = canonical_present ? LINKED_TREATMENT_CANONICAL_DATA_COMMIT : "",
        canonical_root = canonical_root,
        treated_monoculture_label_correction = "legacy IC25 folder -> 1.47 uM/IC75; legacy IC75 folder -> 0.67 uM/IC25",
    ) for path in paths]
    result = DataFrame(rows)
    CSV.write(joinpath(out.csv, "linked_treatment_data_provenance.csv"), result)
    return result
end

function _linked_resume_parameters(path, hypothesis, expected_length)
    isfile(path) || return nothing
    ranking = CSV.read(path, DataFrame)
    rows = ranking[String.(ranking.model) .== hypothesis, :]
    nrow(rows) == 1 || return nothing
    matched = match(r"values = \[([^\]]+)\]", String(first(rows.params)))
    matched === nothing && return nothing
    values = try
        parse.(Float64, strip.(split(matched.captures[1], ",")))
    catch
        return nothing
    end
    return length(values) == expected_length ? values : nothing
end

function _fit_linked_treatment_joint(
    coculture_environments,
    out;
    start,
    max_time_per_fit,
    resume_from_existing = true,
)
    seed = _linked_treatment_seed(start)
    monoculture_environments, monoculture_path = _linked_monoculture_environments(start)
    root = IOUtils.package_root(start)
    coculture_path = joinpath(root, "outputs", "csv", "coculture_treated", "coculture_treated_a2780_decoded.csv")
    _write_linked_provenance(
        [
            monoculture_path,
            coculture_path,
            seed.parameter_path,
            seed.status_path,
            seed.timing_ranking_path,
            seed.timing_parameter_path,
        ],
        out;
        start = start,
    )
    CSV.write(joinpath(out.csv, "linked_treatment_stage2_seed_audit.csv"), seed.audit)

    rows = NamedTuple[]
    fit_cache = Dict{String,Any}()
    overlay_parts = DataFrame[]
    resume_path = joinpath(out.csv, "linked_treatment_model_ranking.csv")
    for hypothesis in LINKED_TREATMENT_HYPOTHESES
        problem = _linked_problem(monoculture_environments, coculture_environments, start, seed, hypothesis)
        resumed = resume_from_existing ?
            _linked_resume_parameters(resume_path, hypothesis, length(problem.p0)) :
            nothing
        if resumed !== nothing && all(first(problem.bounds[index]) <= resumed[index] <= last(problem.bounds[index]) for index in eachindex(resumed))
            problem = merge(problem, (p0 = resumed, u0 = Float64.(problem.u0_builder(resumed))))
        end
        maxiters = hypothesis == "fully_free_context_diagnostic" ? max(400, Int(round(max_time_per_fit * 80))) : max(300, Int(round(max_time_per_fit * 60)))
        println("  Linked treatment fit $(hypothesis)...")
        fit = GrowthParameterEstimation.run_joint_fit(
            problem.model, problem.datasets, problem.u0, problem.p0;
            bounds = problem.bounds,
            u0_builder = problem.u0_builder,
            solver = Rodas5(),
            optimizer = :nelder_mead,
            maxiters = maxiters,
            reltol = 1e-7,
            abstol = 1e-7,
            initial_time = 0.0,
        )
        isfinite(fit.bic) && isfinite(fit.raw_sse) && fit.raw_sse < 9.99e11 || continue
        fit_cache[hypothesis] = (fit = fit, problem = problem)
        push!(rows, (
            cell_line = COCULTURE_SYSTEM_LABEL,
            model = hypothesis,
            pooling_mode = "linked_global",
            eligible_for_inheritance = hypothesis != "fully_free_context_diagnostic",
            diagnostic_model = hypothesis == "fully_free_context_diagnostic",
            bic = Float64(fit.bic),
            ssr = Float64(fit.raw_sse),
            scaled_ssr = Float64(fit.scaled_sse),
            n_parameters = length(fit.params),
            n_points = sum(length(dataset.x) for dataset in problem.datasets),
            n_densities = 2,
            n_mixes = 3,
            n_monoculture_trajectories = 12,
            n_coculture_component_trajectories = 12,
            inherited_timing_hypothesis = problem.timing_hypothesis,
            parameter_scope = hypothesis == "fully_free_context_diagnostic" ? "separate monoculture and coculture intrinsic drug vectors; diagnostic only" : "one intrinsic lineage-specific drug vector shared across monoculture and coculture",
            initial_condition_strategy = "fixed day-zero density totals; coculture split by nominal mix and inherited cis tolerant fraction",
            residual_scaling = "trajectory peak across combined 24-trajectory objective",
            boundary_issue = false,
            params = string((names = problem.names, values = fit.params, package_api = "GrowthParameterEstimation.run_joint_fit")),
        ))
    end
    strobl = _fit_strobl_linked_benchmarks(
        monoculture_environments, coculture_environments, start;
        max_time_per_fit = max_time_per_fit,
    )
    append!(rows, strobl.rows)
    ranking = sort!(DataFrame(rows), :bic)
    expected_model_count = length(LINKED_TREATMENT_HYPOTHESES) + length(STROBL_MODEL_VARIANTS)
    nrow(ranking) == expected_model_count || error("Not all linked treatment and Strobl benchmark hypotheses produced finite fits")
    eligible = ranking[Bool.(ranking.eligible_for_inheritance), :]
    winner = eligible[argmin(eligible.bic), :]
    cached = fit_cache[String(winner.model)]
    explicit_lower = Dict(name => [-4.5, -6.0] for name in cached.problem.names if startswith(String(name), "beta_"))
    for name in cached.problem.names
        name == :log_cis_tolerant_context && (explicit_lower[name] = [-6.0, -8.0])
        name == :cis_tolerant_logit_shift && (explicit_lower[name] = [-4.5, -6.0])
        name == :log_cis_tolerant_growth_context && (explicit_lower[name] = [-log(3.0), -log(4.0)])
    end
    physical_lower = Dict(name => -6.0 for name in keys(explicit_lower))
    haskey(physical_lower, :log_cis_tolerant_context) && (physical_lower[:log_cis_tolerant_context] = -10.0)
    explicit_upper = Dict{Symbol,Vector{Float64}}()
    for name in cached.problem.names
        occursin("lambda", String(name)) && (explicit_upper[name] = [7.5, 10.0])
        startswith(String(name), "beta_") && (explicit_upper[name] = [4.5, 6.0])
        occursin("emax", String(name)) && (explicit_upper[name] = [7.5, 10.0])
        name == :log_cis_tolerant_context && (explicit_upper[name] = [6.0, 8.0])
        name == :cis_tolerant_logit_shift && (explicit_upper[name] = [4.5, 6.0])
        name == :log_cis_tolerant_growth_context && (explicit_upper[name] = [log(3.0), log(4.0)])
    end
    modifier_indices = if cached.problem.base.nbase < length(cached.problem.names)
        collect((cached.problem.base.nbase + 1):length(cached.problem.names))
    else
        Int[]
    end
    expandable_parameters = [
        cached.problem.names[index]
        for index in modifier_indices
        if !occursin("f_tolerant0", String(cached.problem.names[index])) &&
            !occursin("log_contrast_density", String(cached.problem.names[index]))
    ]
    profiled = GrowthParameterEstimation.profile_joint_fit_bounds_two_sided(
        cached.problem.model, cached.problem.datasets, cached.problem.u0, cached.fit.params;
        bounds = cached.problem.bounds,
        parameter_names = cached.problem.names,
        explicit_lower_profiles = explicit_lower,
        physical_lower_limits = physical_lower,
        explicit_upper_profiles = explicit_upper,
        profile_parameters = expandable_parameters,
        u0_builder = cached.problem.u0_builder,
        solver = Rodas5(),
        optimizer = :nelder_mead,
        maxiters = max(400, Int(round(max_time_per_fit * 70))),
        reltol = 1e-7,
        abstol = 1e-7,
        initial_time = 0.0,
    )
    fit_cache[String(winner.model)] = (fit = profiled.fit, problem = merge(cached.problem, (bounds = profiled.bounds,)))
    winner_index = findfirst(String.(ranking.model) .== String(winner.model))
    ranking.bic[winner_index] = Float64(profiled.fit.bic)
    ranking.ssr[winner_index] = Float64(profiled.fit.raw_sse)
    ranking.scaled_ssr[winner_index] = Float64(profiled.fit.scaled_sse)
    ranking.boundary_issue[winner_index] = any(profiled.identifiability.identifiability .!= "interior")
    ranking.params[winner_index] = string((
        names = cached.problem.names,
        values = profiled.fit.params,
        package_api = "GrowthParameterEstimation.run_joint_fit + profile_joint_fit_bounds_two_sided",
    ))
    sort!(ranking, :bic)
    eligible = ranking[Bool.(ranking.eligible_for_inheritance), :]
    winner = eligible[argmin(eligible.bic), :]

    for hypothesis in LINKED_TREATMENT_HYPOTHESES
        cached_hypothesis = fit_cache[hypothesis]
        push!(overlay_parts, _linked_overlay(cached_hypothesis.fit, cached_hypothesis.problem, hypothesis))
    end
    append!(overlay_parts, strobl.overlays)
    overlay = vcat(overlay_parts...; cols = :union)
    free_row = first(ranking[String.(ranking.model) .== "fully_free_context_diagnostic", :])
    free_improvement = Float64(winner.bic) - Float64(free_row.bic)
    status = DataFrame(
        cell_line = [COCULTURE_SYSTEM_LABEL],
        winning_model = [String(winner.model)],
        winning_pooling_mode = ["linked_global"],
        winning_bic = [Float64(winner.bic)],
        fully_free_diagnostic_bic = [Float64(free_row.bic)],
        fully_free_bic_improvement = [free_improvement],
        inadequacy_delta = [10.0],
        inadequate_inheritance = [free_improvement >= 10.0],
        inheritance_allowed = [free_improvement < 10.0],
    )
    ranking.delta_bic = ranking.bic .- minimum(ranking.bic)
    ranking.rank = collect(1:nrow(ranking))

    winner_cached = fit_cache[String(winner.model)]
    parameter_rows = NamedTuple[]
    base = _linked_base_spec(seed)
    for (index, name) in enumerate(winner_cached.problem.names)
        seed_value = index <= base.nbase ? base.p0[index] : NaN
        push!(parameter_rows, (
            parameter = String(name),
            estimate = Float64(winner_cached.fit.params[index]),
            sequential_stage2_seed = seed_value,
            deviation_percent = isfinite(seed_value) && abs(seed_value) > 1e-12 ? 100 * (Float64(winner_cached.fit.params[index]) / seed_value - 1) : NaN,
            lower_bound = Float64(winner_cached.problem.bounds[index][1]),
            upper_bound = Float64(winner_cached.problem.bounds[index][2]),
            source = index <= base.nbase ? "shared intrinsic parameter jointly estimated from treated monoculture and coculture" : "coculture-specific modifier",
        ))
    end
    shared_audit = DataFrame(parameter_rows)

    effective_rows = NamedTuple[]
    for density in ("20k", "30k")
        mono = _linked_treatment_values(winner_cached.fit.params, String(winner.model), :monoculture, density, winner_cached.problem.base)
        co = _linked_treatment_values(winner_cached.fit.params, String(winner.model), :coculture, density, winner_cached.problem.base)
        for (lineage, names) in (("A2780Naive", (:naive_emax, :naive_ec50, :naive_hill, :naive_lambda, :naive_onset)), ("A2780cis", (:cis_emax_sensitive, :cis_emax_tolerant, :cis_lambda, :cis_onset, :cis_f_tolerant0)))
            for name in names
                mono_value = Float64(getproperty(mono, name))
                co_value = Float64(getproperty(co, name))
                push!(effective_rows, (
                    lineage = lineage,
                    density = density,
                    parameter = String(name),
                    treated_monoculture_value = mono_value,
                    treated_coculture_intrinsic_value = co_value,
                    ratio = abs(mono_value) > 1e-12 ? co_value / mono_value : NaN,
                    inheritance_mode = "same parameter index in one combined global objective",
                ))
            end
        end
    end

    CSV.write(joinpath(out.csv, "linked_treatment_model_ranking.csv"), ranking)
    CSV.write(joinpath(out.csv, "linked_treatment_top5.csv"), first(ranking, min(LINKED_TREATMENT_VISIBLE_LIMIT, nrow(ranking))))
    CSV.write(joinpath(out.csv, "linked_treatment_status.csv"), status)
    CSV.write(joinpath(out.csv, "linked_treatment_shared_parameter_audit.csv"), shared_audit)
    CSV.write(joinpath(out.csv, "linked_treatment_effective_parameter_inheritance.csv"), DataFrame(effective_rows))
    CSV.write(joinpath(out.csv, "linked_treatment_boundary_profiles.csv"), profiled.profile)
    CSV.write(joinpath(out.csv, "linked_treatment_identifiability.csv"), profiled.identifiability)
    figure_csv = joinpath(out.csv, "figures"); mkpath(figure_csv)
    CSV.write(joinpath(figure_csv, "linked_treatment_combined_overlays.csv"), overlay)
    _render_linked_monoculture_grid(overlay, winner, out)
    _render_linked_coculture_grid(overlay, winner, out)
    _render_linked_hypothesis_grid(overlay, winner, out)
    return ranking
end
