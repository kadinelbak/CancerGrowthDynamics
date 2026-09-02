# Joint A2780 coculture fitting. This file is included inside FitWorkflows.

const COCULTURE_SYSTEM_LABEL = "A2780Naive+A2780cis"

function _component_label(value, source_file)
    text = lowercase("$(_safe_string(value; default = "")) $(source_file)")
    return occursin("cis", text) ? "resistant" : "sensitive"
end

function _mix_fraction(label::AbstractString)
    matched = match(r"^(\d+)-(\d+)$", strip(label))
    matched === nothing && return NaN
    first_part = parse(Float64, matched.captures[1])
    second_part = parse(Float64, matched.captures[2])
    total = first_part + second_part
    return total > 0 ? first_part / total : NaN
end

"""
Recover density provenance that was lost when cached decoded files retained only
the basename of identically named 20k and 30k source files. Sorted source paths
placed the 20k replicate block first. Treated day-level files provide an
independent check of the recovered block means.
"""
function _recover_coculture_design(decoded::DataFrame, condition::AbstractString)
    required = Set([:time, :count, :source_file, :cell_line, :mix])
    all(column -> column in propertynames(decoded), required) ||
        error("Coculture decoding requires time, count, source_file, cell_line, and mix")

    treated = condition == "coculture_treated"
    block_size = treated ? 6 : 3
    source_names = lowercase.(String.(decoded.source_file))
    specific = decoded[
        startswith.(source_names, "measure_") .&
        endswith.(source_names, "_well_day_averages.csv"),
        :,
    ]
    isempty(specific) && error("No condition-specific well/day coculture rows were available")

    rows = NamedTuple[]
    provenance = NamedTuple[]
    for group in groupby(specific, [:source_file, :time])
        source_file = String(first(group.source_file))
        n = nrow(group)
        n > block_size || error("$(source_file) time=$(first(group.time)) has only $(n) rows; both density blocks are required")
        n <= 2block_size || error("$(source_file) time=$(first(group.time)) has $(n) rows; expected at most $(2block_size)")
        for (index, row) in enumerate(eachrow(group))
            density = index <= block_size ? "20k" : "30k"
            replicate = index <= block_size ? index : index - block_size
            push!(rows, (
                time = Float64(row.time),
                count = Float64(row.count),
                source_file = source_file,
                component = _component_label(row.cell_line, source_file),
                density = density,
                mix = _safe_string(row.mix; default = "unknown_mix"),
                replicate = replicate,
            ))
        end
    end
    recovered = sort!(DataFrame(rows), [:density, :mix, :component, :time, :replicate])

    for group in groupby(recovered, [:source_file, :time, :density])
        source_file = String(first(group.source_file))
        density = String(first(group.density))
        block_mean = mean(Float64.(group.count))
        day_source = replace(source_file, "_well_day_averages.csv" => "_day_averages.csv")
        day_rows = decoded[
            (String.(decoded.source_file) .== day_source) .&
            (Float64.(decoded.time) .== Float64(first(group.time))),
            :,
        ]
        day_value = NaN
        if treated && nrow(day_rows) >= 2
            day_value = Float64(day_rows[density == "20k" ? 1 : 2, :count])
        end
        push!(provenance, (
            source_file = source_file,
            time = Float64(first(group.time)),
            density = density,
            rows_in_block = nrow(group),
            recovered_mean = block_mean,
            day_aggregate = day_value,
            absolute_check_error = isfinite(day_value) ? abs(block_mean - day_value) : NaN,
            recovery_rule = "sorted source path block: 20k first, 30k second",
        ))
    end
    return recovered, DataFrame(provenance)
end

function _coculture_environments(recovered::DataFrame)
    environments = NamedTuple[]
    diagnostics = NamedTuple[]
    for environment in groupby(recovered, [:density, :mix])
        density = String(first(environment.density))
        mix = String(first(environment.mix))
        curves = Dict{String,DataFrame}()
        for component in ("sensitive", "resistant")
            component_rows = environment[String.(environment.component) .== component, :]
            isempty(component_rows) && error("Missing $(component) observations for $(density) $(mix)")
            curve = combine(groupby(component_rows, :time), :count => mean => :observed)
            sort!(curve, :time)
            nrow(curve) >= 3 || error("Insufficient observations for $(density) $(mix) $(component)")
            curves[component] = curve
        end
        sensitive = curves["sensitive"]
        resistant = curves["resistant"]
        sensitive.time == resistant.time || error("Component time grids differ for $(density) $(mix)")
        first_sensitive_value = max(Float64(first(sensitive.observed)), 1e-6)
        first_resistant_value = max(Float64(first(resistant.observed)), 1e-6)
        nominal_sensitive_fraction = _mix_fraction(mix)
        isfinite(nominal_sensitive_fraction) || error("Cannot initialize unknown coculture mix $(mix)")
        fixed_total_u0 = _fixed_day0_total(density)
        u0_sensitive = fixed_total_u0 * nominal_sensitive_fraction
        u0_resistant = fixed_total_u0 * (1 - nominal_sensitive_fraction)
        observed_sensitive_fraction = first_sensitive_value / (first_sensitive_value + first_resistant_value)
        push!(environments, (
            density = density,
            mix = mix,
            times = Float64.(sensitive.time),
            sensitive = Float64.(sensitive.observed),
            resistant = Float64.(resistant.observed),
            u0_sensitive = u0_sensitive,
            u0_resistant = u0_resistant,
            fixed_total_u0 = fixed_total_u0,
            first_sensitive_value = first_sensitive_value,
            first_resistant_value = first_resistant_value,
            sensitive_scale = max(maximum(Float64.(sensitive.observed)), 1.0),
            resistant_scale = max(maximum(Float64.(resistant.observed)), 1.0),
        ))
        push!(diagnostics, (
            density = density,
            mix = mix,
            nominal_sensitive_fraction = nominal_sensitive_fraction,
            fixed_total_u0 = fixed_total_u0,
            fixed_sensitive_u0 = u0_sensitive,
            fixed_resistant_u0 = u0_resistant,
            first_observed_sensitive = first_sensitive_value,
            first_observed_resistant = first_resistant_value,
            first_observed_total = first_sensitive_value + first_resistant_value,
            observed_sensitive_fraction = observed_sensitive_fraction,
            fraction_difference = observed_sensitive_fraction - nominal_sensitive_fraction,
            initial_condition_strategy = "fixed nominal day-zero total (67 for 20k, 100 for 30k) split by mix",
            mix_strategy = "first mix fraction initializes sensitive; second initializes cis-resistant",
        ))
    end
    sort!(environments, by = environment -> (_density_nominal_value(environment.density), environment.mix))
    return environments, DataFrame(diagnostics)
end

function _coculture_monoculture_baselines(start, density)
    baseline_map = _load_untreated_monoculture_cellline_baselines(; start = start)
    baseline_map === nothing && error("Monoculture baseline artifact is required for coculture fitting")
    sensitive_key = ("A2780Naive", density)
    resistant_key = ("A2780cis", density)
    haskey(baseline_map.spec_by_cell_density, sensitive_key) || error("Missing $(sensitive_key) baseline")
    haskey(baseline_map.spec_by_cell_density, resistant_key) || error("Missing $(resistant_key) baseline")
    sensitive = baseline_map.spec_by_cell_density[sensitive_key]
    resistant = baseline_map.spec_by_cell_density[resistant_key]
    sensitive.inheritance_allowed || error("Sensitive monoculture baseline is not eligible for inheritance")
    resistant.inheritance_allowed || error("Resistant monoculture baseline is not eligible for inheritance")
    return sensitive, resistant
end

function _anchored_component_growth(population, competitive_load, baseline, t)
    return _baseline_growth(population, competitive_load, baseline, t)
end

function _pooled_local_parameters(base_p0, base_bounds, base_names, pooling_mode, density_indices; contrast_indices = collect(eachindex(base_p0)))
    nbase = length(base_p0)
    if pooling_mode == "shared"
        p0 = copy(base_p0)
        bounds = copy(base_bounds)
        names = copy(base_names)
    elseif pooling_mode == "partial_5pct"
        p0 = vcat(base_p0, 0.0)
        bounds = vcat(base_bounds, [(-DENSITY_LOG_CONTRAST_BOUND, DENSITY_LOG_CONTRAST_BOUND)])
        names = vcat(base_names, [:log_contrast_density])
    elseif pooling_mode == "independent_diagnostic"
        p0 = vcat(base_p0, base_p0)
        bounds = vcat(base_bounds, base_bounds)
        names = vcat(Symbol.(string.(base_names) .* "_20k"), Symbol.(string.(base_names) .* "_30k"))
    else
        error("Unsupported coculture pooling mode: $(pooling_mode)")
    end
    local_parameters = function (p, environment_index)
        density_index = density_indices[environment_index]
        values = if pooling_mode == "independent_diagnostic"
            offset = (density_index - 1) * nbase
            collect(p[(offset + 1):(offset + nbase)])
        else
            collect(p[1:nbase])
        end
        if pooling_mode == "partial_5pct"
            multiplier = exp((density_index == 1 ? -1.0 : 1.0) * p[nbase + 1])
            for parameter_index in contrast_indices
                values[parameter_index] *= multiplier
            end
        end
        return values
    end
    return (p0 = p0, bounds = bounds, names = names, base_names = base_names,
            base_bounds = base_bounds, nbase = nbase, local_parameters = local_parameters)
end

function _untreated_coculture_specs(environments, baselines, pooling_mode)
    density_indices = [environment.density == "20k" ? 1 : 2 for environment in environments]
    definitions = Dict(
        "lv_symmetric_competition" => ([1.0], [(0.02, 4.0)], [:alpha], [1]),
        "lv_asymmetric_competition" => ([1.0, 1.0], [(0.02, 4.0), (0.02, 4.0)], [:alpha_sr, :alpha_rs], [1, 2]),
        "lv_asymmetric_competition_death" => ([1.0, 1.0, 0.02, 0.02], [(0.02, 4.0), (0.02, 4.0), (0.0, 1.0), (0.0, 1.0)], [:alpha_sr, :alpha_rs, :death_sensitive, :death_resistant], [1, 2]),
    )
    specs = Dict{String,NamedTuple}()
    for (model_name, definition) in definitions
        base_p0, base_bounds, base_names, contrast_indices = definition
        pooled = _pooled_local_parameters(base_p0, base_bounds, base_names, pooling_mode, density_indices; contrast_indices = contrast_indices)
        model! = function (du, u, p, t)
            for environment_index in eachindex(environments)
                sensitive_index = 2environment_index - 1
                resistant_index = 2environment_index
                sensitive = max(u[sensitive_index], zero(u[sensitive_index]))
                resistant = max(u[resistant_index], zero(u[resistant_index]))
                local_values = pooled.local_parameters(p, environment_index)
                if model_name == "lv_symmetric_competition"
                    alpha_sr, alpha_rs, death_sensitive, death_resistant = local_values[1], local_values[1], 0.0, 0.0
                elseif model_name == "lv_asymmetric_competition"
                    alpha_sr, alpha_rs, death_sensitive, death_resistant = local_values[1], local_values[2], 0.0, 0.0
                else
                    alpha_sr, alpha_rs, death_sensitive, death_resistant = local_values
                end
                sensitive_baseline, resistant_baseline = baselines[environment_index]
                sensitive_load = sensitive + alpha_sr * resistant
                resistant_load = resistant + alpha_rs * sensitive
                du[sensitive_index] = _anchored_component_growth(sensitive, sensitive_load, sensitive_baseline, t) - death_sensitive * sensitive
                du[resistant_index] = _anchored_component_growth(resistant, resistant_load, resistant_baseline, t) - death_resistant * resistant
            end
        end
        specs[model_name] = merge(pooled, (model = model!,))
    end
    return specs
end

function _untreated_coculture_inputs(environments)
    datasets = NamedTuple[]
    u0 = Float64[]
    metadata = NamedTuple[]
    for (environment_index, environment) in enumerate(environments)
        append!(u0, [environment.u0_sensitive, environment.u0_resistant])
        for (component, values, scale, state_index) in (
            ("sensitive", environment.sensitive, environment.sensitive_scale, 2environment_index - 1),
            ("resistant", environment.resistant, environment.resistant_scale, 2environment_index),
        )
            push!(datasets, (x = environment.times, y = values, state_index = state_index, residual_scale = scale))
            push!(metadata, (environment_index = environment_index, density = environment.density, mix = environment.mix, component = component))
        end
    end
    return datasets, u0, metadata
end

function _fit_untreated_coculture_joint(environments, out; start, max_time_per_fit)
    baselines = [_coculture_monoculture_baselines(start, environment.density) for environment in environments]
    inheritance_rows = NamedTuple[]
    for density in ("20k", "30k")
        sensitive, resistant = _coculture_monoculture_baselines(start, density)
        for (component, baseline) in (("sensitive", sensitive), ("resistant", resistant))
            push!(inheritance_rows, (
                density = density,
                component = component,
                monoculture_cell_line = component == "sensitive" ? "A2780Naive" : "A2780cis",
                growth_family = String(baseline.model),
                r = Float64(baseline.r),
                K = Float64(baseline.K),
                shape_parameter = String(baseline.shape_parameter),
                shape_value = Float64(baseline.shape_value),
                carrying_capacity_mode = "exact_fixed_untreated_K_no_rescaling",
                inheritance_mode = "exact_fixed_monoculture_family_and_density_specific_parameters",
            ))
        end
    end
    datasets, u0, metadata = _untreated_coculture_inputs(environments)
    rows = NamedTuple[]
    overlays = DataFrame[]
    parameter_rows = NamedTuple[]
    profile_parts = DataFrame[]
    identifiability_parts = DataFrame[]
    fit_cache = Dict{Tuple{String,String},Any}()
    pooling_modes = _tournament_smoke_test() ? ("shared", "partial_5pct") : ("shared", "partial_5pct", "independent_diagnostic")
    for pooling_mode in pooling_modes
        specs = _untreated_coculture_specs(environments, baselines, pooling_mode)
        model_names = sort(collect(keys(specs)))
        _tournament_smoke_test() && (model_names = filter(name -> name in ("lv_symmetric_competition", "lv_asymmetric_competition", "lv_asymmetric_competition_death"), model_names))
        for model_name in model_names
            spec = specs[model_name]
            println("  Joint untreated coculture fit $(model_name), $(pooling_mode)...")
            model_maxiters = max(500, Int(round(max_time_per_fit * 18)))
            _tournament_smoke_test() && (model_maxiters = 1)
            multistart = try
                starts = GrowthParameterEstimation.generate_multistarts(
                    spec.p0, spec.bounds;
                    n_starts = _configured_multistarts(1),
                    seed = _stable_multistart_seed("untreated_coculture", model_name, pooling_mode),
                )
                GrowthParameterEstimation.run_joint_multistart(
                    spec.model, datasets, u0, starts;
                    refine_optimizer = _configured_refine_optimizer(),
                    bounds = spec.bounds,
                    solver = Tsit5(),
                    maxiters = model_maxiters,
                    maxtime = _configured_fit_maxtime(),
                    reltol = 1e-7,
                    abstol = 1e-7,
                    optimizer = :nelder_mead,
                    initial_time = 0.0,
                )
            catch error_value
                @warn "Untreated coculture joint multistart failed" model_name pooling_mode exception = error_value
                nothing
            end
            multistart === nothing && continue
            profiled = try
                GrowthParameterEstimation.profile_joint_fit_bounds(
                    spec.model, datasets, u0, multistart.fit.params;
                    bounds = spec.bounds,
                    parameter_names = spec.names,
                    solver = Tsit5(),
                    maxiters = model_maxiters,
                    maxtime = _configured_fit_maxtime(),
                    reltol = 1e-7,
                    abstol = 1e-7,
                    optimizer = :nelder_mead,
                    initial_time = 0.0,
                )
            catch error_value
                @warn "Untreated coculture joint fit failed" model_name pooling_mode exception = error_value
                nothing
            end
            profiled === nothing && continue
            fit = profiled.fit
            isfinite(fit.bic) && isfinite(fit.raw_sse) && fit.raw_sse < 9.99e11 || continue
            profile = DataFrame(profiled.profile)
            if !isempty(profile)
                profile.model = fill(model_name, nrow(profile))
                profile.pooling_mode = fill(pooling_mode, nrow(profile))
                push!(profile_parts, profile)
            end
            identifiability = DataFrame(profiled.identifiability)
            identifiability.model = fill(model_name, nrow(identifiability))
            identifiability.pooling_mode = fill(pooling_mode, nrow(identifiability))
            push!(identifiability_parts, identifiability)
            boundary_issue = any(identifiability.identifiability .!= "interior")
            fit_cache[(model_name, pooling_mode)] = (fit = fit, spec = spec)
            push!(rows, (
                cell_line = COCULTURE_SYSTEM_LABEL,
                model = model_name,
                pooling_mode = pooling_mode,
                eligible_for_inheritance = pooling_mode != "independent_diagnostic",
                bic = Float64(fit.bic),
                ssr = Float64(fit.raw_sse),
                scaled_ssr = Float64(fit.scaled_sse),
                n_parameters = length(fit.params),
                n_points = sum(length(dataset.x) for dataset in datasets),
                n_densities = 2,
                n_mixes = length(unique(environment.mix for environment in environments)),
                parameter_scope = pooling_mode == "partial_5pct" ? "shared center with symmetric density effects bounded at +/-5%" : pooling_mode,
                growth_inheritance = "exact density-specific sensitive/resistant monoculture growth family/r/K; no K rescaling",
                initial_condition_strategy = "fixed nominal day-zero total split by mix",
                residual_scaling = "component-trajectory peak",
                boundary_issue = boundary_issue,
                params = string((names = spec.names, values = fit.params, package_api = "GrowthParameterEstimation.profile_joint_fit_bounds/run_joint_fit")),
            ))
            for (dataset_index, dataset) in enumerate(datasets)
                meta = metadata[dataset_index]
                component_u0 = metadata[dataset_index].component == "sensitive" ? environments[metadata[dataset_index].environment_index].u0_sensitive : environments[metadata[dataset_index].environment_index].u0_resistant
                overlay_length = length(dataset.x) + 1
                push!(overlays, DataFrame(
                    time = vcat(0.0, dataset.x),
                    observed = vcat(component_u0, dataset.y),
                    predicted = vcat(component_u0, fit.predictions[dataset_index]),
                    fixed_day0_anchor = vcat(true, fill(false, length(dataset.x))),
                    model = fill(model_name, overlay_length),
                    pooling_mode = fill(pooling_mode, overlay_length),
                    density = fill(meta.density, overlay_length),
                    mix = fill(meta.mix, overlay_length),
                    component = fill(meta.component, overlay_length),
                ))
            end
            for environment_index in eachindex(environments)
                local_values = spec.local_parameters(fit.params, environment_index)
                for parameter_index in eachindex(spec.base_names)
                    push!(parameter_rows, (
                        stage = "coculture_untreated",
                        model = model_name,
                        pooling_mode = pooling_mode,
                        density = environments[environment_index].density,
                        parameter = String(spec.base_names[parameter_index]),
                        effective_value = Float64(local_values[parameter_index]),
                        lower_bound = Float64(spec.base_bounds[parameter_index][1]),
                        upper_bound = Float64(spec.base_bounds[parameter_index][2]),
                        bound_position = (Float64(local_values[parameter_index]) - spec.base_bounds[parameter_index][1]) / max(spec.base_bounds[parameter_index][2] - spec.base_bounds[parameter_index][1], eps()),
                    ))
                end
            end
        end
    end
    strobl = _fit_strobl_untreated_benchmarks(
        environments, baselines, datasets, u0, metadata;
        max_time_per_fit = max_time_per_fit,
    )
    append!(rows, strobl.rows)
    append!(overlays, strobl.overlays)
    append!(parameter_rows, strobl.parameter_rows)
    ranking = sort!(DataFrame(rows), :bic)
    isempty(ranking) && error("No finite untreated coculture joint fits were produced")
    eligible = ranking[Bool.(ranking.eligible_for_inheritance), :]
    isempty(eligible) && error("No pooled untreated coculture fit is eligible for inheritance")
    winner = eligible[argmin(eligible.bic), :]
    cached = fit_cache[(String(winner.model), String(winner.pooling_mode))]
    pooling_summary = GrowthParameterEstimation.summarize_pooling_bic(ranking; top_n = 5)
    baseline = DataFrame(
        model = [String(winner.model)],
        pooling_mode = [String(winner.pooling_mode)],
        parameter_names = [join(string.(cached.spec.names), ";")],
        parameter_values = [join(string.(cached.fit.params), ";")],
        bic = [Float64(winner.bic)],
        ssr = [Float64(winner.ssr)],
        scaled_ssr = [Float64(winner.scaled_ssr)],
        monoculture_growth_inheritance = ["exact density-specific family/r/K; no carrying-capacity rescaling"],
    )
    overlay_df = vcat(overlays...; cols = :union)
    CSV.write(joinpath(out.csv, "coculture_untreated_pooling_model_ranking.csv"), ranking)
    CSV.write(joinpath(out.csv, "coculture_untreated_pooling_top5.csv"), pooling_summary.ranking)
    CSV.write(joinpath(out.csv, "coculture_untreated_pooling_status.csv"), pooling_summary.status)
    CSV.write(joinpath(out.csv, "coculture_untreated_joint_parameter_estimates.csv"), unique(DataFrame(parameter_rows)))
    CSV.write(joinpath(out.csv, "coculture_untreated_boundary_profiles.csv"), isempty(profile_parts) ? DataFrame() : vcat(profile_parts...; cols = :union))
    CSV.write(joinpath(out.csv, "coculture_untreated_identifiability.csv"), isempty(identifiability_parts) ? DataFrame() : vcat(identifiability_parts...; cols = :union))
    CSV.write(joinpath(out.csv, "coculture_untreated_joint_baseline.csv"), baseline)
    CSV.write(joinpath(out.csv, "coculture_untreated_monoculture_inheritance_audit.csv"), DataFrame(inheritance_rows))
    figure_csv = joinpath(out.csv, "figures")
    mkpath(figure_csv)
    CSV.write(joinpath(figure_csv, "coculture_untreated_joint_overlays.csv"), overlay_df)
    _render_coculture_graph_grid(overlay_df, winner, out, "coculture_untreated")
    return ranking
end

function _load_coculture_untreated_baseline(; start)
    root = IOUtils.package_root(start)
    path = joinpath(root, "outputs", "csv", "coculture_untreated", "coculture_untreated_joint_baseline.csv")
    isfile(path) || error("Untreated coculture baseline artifact is required before treated coculture fitting")
    row = first(CSV.read(path, DataFrame))
    values = parse.(Float64, split(String(row.parameter_values), ";"))
    names = Symbol.(split(String(row.parameter_names), ";"))
    return (model = String(row.model), pooling_mode = String(row.pooling_mode), params = values, names = names, path = path)
end

function _treated_coculture_specs(environments, baselines, untreated_baseline, untreated_specs, pooling_mode, t0)
    density_indices = [environment.density == "20k" ? 1 : 2 for environment in environments]
    definitions = Dict(
        "dual_constant_kill" => ([0.08, 0.03], [(0.0, 3.0), (0.0, 3.0)], [:kill_sensitive, :kill_resistant], [1, 2], :live),
        "dual_time_decay_kill" => ([0.4, 0.15, 0.2], [(0.0, 4.0), (0.0, 4.0), (0.0, 5.0)], [:kill_sensitive, :kill_resistant, :lambda], [1, 2], :live),
        "dual_delayed_ramp_kill" => ([0.4, 0.15, 0.8, 1.0], [(0.0, 4.0), (0.0, 4.0), (0.01, 5.0), (0.0, 7.0)], [:kill_sensitive, :kill_resistant, :lambda, :onset], [1, 2], :live),
        "dual_transit_damage" => ([0.5, 0.2, 0.8, 1.0, 0.4], [(0.0, 4.0), (0.0, 4.0), (0.01, 5.0), (0.0, 7.0), (0.01, 4.0)], [:kill_sensitive, :kill_resistant, :lambda, :onset, :clearance], [1, 2], :transit),
        "dual_transit_competitor_scaled" => ([0.5, 0.2, 0.8, 1.0, 0.4, 0.0, 0.0], [(0.0, 4.0), (0.0, 4.0), (0.01, 5.0), (0.0, 7.0), (0.01, 4.0), (-3.0, 3.0), (-3.0, 3.0)], [:kill_sensitive, :kill_resistant, :lambda, :onset, :clearance, :beta_sensitive, :beta_resistant], [1, 2], :transit),
        "dual_transit_load_scaled" => ([0.5, 0.2, 0.8, 1.0, 0.4, 0.0, 0.0], [(0.0, 4.0), (0.0, 4.0), (0.01, 5.0), (0.0, 7.0), (0.01, 4.0), (-3.0, 3.0), (-3.0, 3.0)], [:kill_sensitive, :kill_resistant, :lambda, :onset, :clearance, :beta_sensitive, :beta_resistant], [1, 2], :transit),
        "sensitive_tolerant_transition" => ([0.6, 0.08, 0.15, 0.8, 1.0, 0.08], [(0.0, 4.0), (0.0, 3.0), (0.0, 4.0), (0.01, 5.0), (0.0, 7.0), (0.0, 1.5)], [:kill_sensitive, :kill_tolerant, :kill_resistant, :lambda, :onset, :transition], [1, 2, 3], :tolerant),
    )
    specs = Dict{String,NamedTuple}()
    untreated_spec = untreated_specs[untreated_baseline.model]
    for (model_name, definition) in definitions
        base_p0, base_bounds, base_names, amplitude_indices, layout = definition
        pooled = _pooled_local_parameters(base_p0, base_bounds, base_names, pooling_mode, density_indices; contrast_indices = amplitude_indices)
        states_per_environment = layout == :live ? 2 : layout == :transit ? 4 : 3
        model! = function (du, u, p, t)
            elapsed = max(t - t0, 0.0)
            for environment_index in eachindex(environments)
                offset = (environment_index - 1) * states_per_environment
                values = pooled.local_parameters(p, environment_index)
                untreated_values = untreated_spec.local_parameters(untreated_baseline.params, environment_index)
                if untreated_baseline.model == "lv_symmetric_competition"
                    alpha_sr, alpha_rs, death_sensitive, death_resistant = untreated_values[1], untreated_values[1], 0.0, 0.0
                elseif untreated_baseline.model == "lv_asymmetric_competition"
                    alpha_sr, alpha_rs, death_sensitive, death_resistant = untreated_values[1], untreated_values[2], 0.0, 0.0
                else
                    alpha_sr, alpha_rs, death_sensitive, death_resistant = untreated_values
                end
                sensitive_baseline, resistant_baseline = baselines[environment_index]
                if layout == :live
                    sensitive = max(u[offset + 1], zero(u[offset + 1]))
                    resistant = max(u[offset + 2], zero(u[offset + 2]))
                    growth_sensitive = _anchored_component_growth(sensitive, sensitive + alpha_sr * resistant, sensitive_baseline, t) - death_sensitive * sensitive
                    growth_resistant = _anchored_component_growth(resistant, resistant + alpha_rs * sensitive, resistant_baseline, t) - death_resistant * resistant
                    if model_name == "dual_constant_kill"
                        kill_sensitive, kill_resistant = values
                        exposure = 1.0
                    elseif model_name == "dual_time_decay_kill"
                        kill_sensitive, kill_resistant, lambda = values
                        exposure = exp(-lambda * elapsed)
                    else
                        kill_sensitive, kill_resistant, lambda, onset = values
                        exposure = elapsed <= onset ? 0.0 : 1 - exp(-lambda * (elapsed - onset))
                    end
                    du[offset + 1] = growth_sensitive - exposure * kill_sensitive * sensitive
                    du[offset + 2] = growth_resistant - exposure * kill_resistant * resistant
                elseif layout == :transit
                    sensitive = max(u[offset + 1], zero(u[offset + 1]))
                    resistant = max(u[offset + 2], zero(u[offset + 2]))
                    damaged_sensitive = max(u[offset + 3], zero(u[offset + 3]))
                    damaged_resistant = max(u[offset + 4], zero(u[offset + 4]))
                    kill_sensitive, kill_resistant, lambda, onset, clearance = values[1:5]
                    exposure = elapsed <= onset ? 0.0 : 1 - exp(-lambda * (elapsed - onset))
                    growth_sensitive = _anchored_component_growth(sensitive, sensitive + alpha_sr * resistant, sensitive_baseline, t) - death_sensitive * sensitive
                    growth_resistant = _anchored_component_growth(resistant, resistant + alpha_rs * sensitive, resistant_baseline, t) - death_resistant * resistant
                    effect_scale_sensitive = 1.0
                    effect_scale_resistant = 1.0
                    if model_name != "dual_transit_damage"
                        beta_sensitive, beta_resistant = values[6], values[7]
                        if model_name == "dual_transit_competitor_scaled"
                            burden_sensitive = alpha_sr * resistant / sensitive_baseline.K
                            burden_resistant = alpha_rs * sensitive / resistant_baseline.K
                        else
                            burden_sensitive = (sensitive + alpha_sr * resistant) / sensitive_baseline.K
                            burden_resistant = (resistant + alpha_rs * sensitive) / resistant_baseline.K
                        end
                        effect_scale_sensitive = exp(clamp(beta_sensitive * burden_sensitive, -4.0, 4.0))
                        effect_scale_resistant = exp(clamp(beta_resistant * burden_resistant, -4.0, 4.0))
                    end
                    flux_sensitive = exposure * kill_sensitive * effect_scale_sensitive * sensitive
                    flux_resistant = exposure * kill_resistant * effect_scale_resistant * resistant
                    du[offset + 1] = growth_sensitive - flux_sensitive
                    du[offset + 2] = growth_resistant - flux_resistant
                    du[offset + 3] = flux_sensitive - clearance * damaged_sensitive
                    du[offset + 4] = flux_resistant - clearance * damaged_resistant
                else
                    sensitive = max(u[offset + 1], zero(u[offset + 1]))
                    tolerant = max(u[offset + 2], zero(u[offset + 2]))
                    resistant = max(u[offset + 3], zero(u[offset + 3]))
                    naive_total = sensitive + tolerant
                    kill_sensitive, kill_tolerant, kill_resistant, lambda, onset, transition = values
                    exposure = elapsed <= onset ? 0.0 : 1 - exp(-lambda * (elapsed - onset))
                    naive_growth = _anchored_component_growth(naive_total, naive_total + alpha_sr * resistant, sensitive_baseline, t) - death_sensitive * naive_total
                    resistant_growth = _anchored_component_growth(resistant, resistant + alpha_rs * naive_total, resistant_baseline, t) - death_resistant * resistant
                    sensitive_share = naive_total > 0 ? sensitive / naive_total : 1.0
                    tolerant_share = naive_total > 0 ? tolerant / naive_total : 0.0
                    du[offset + 1] = naive_growth * sensitive_share - exposure * kill_sensitive * sensitive - transition * sensitive
                    du[offset + 2] = naive_growth * tolerant_share - exposure * kill_tolerant * tolerant + transition * sensitive
                    du[offset + 3] = resistant_growth - exposure * kill_resistant * resistant
                end
            end
        end
        specs[model_name] = merge(pooled, (model = model!, layout = layout, states_per_environment = states_per_environment))
    end
    return specs
end

function _treated_coculture_inputs(environments, layout)
    states_per_environment = layout == :live ? 2 : layout == :transit ? 4 : 3
    u0 = Float64[]
    datasets = NamedTuple[]
    metadata = NamedTuple[]
    for (environment_index, environment) in enumerate(environments)
        offset = (environment_index - 1) * states_per_environment
        if layout == :live
            append!(u0, [environment.u0_sensitive, environment.u0_resistant])
            observables = (("sensitive", environment.sensitive, environment.sensitive_scale, offset + 1), ("resistant", environment.resistant, environment.resistant_scale, offset + 2))
            for (component, values, scale, state_index) in observables
                push!(datasets, (x = environment.times, y = values, state_index = state_index, residual_scale = scale))
                push!(metadata, (density = environment.density, mix = environment.mix, component = component))
            end
        elseif layout == :transit
            append!(u0, [environment.u0_sensitive, environment.u0_resistant, 0.0, 0.0])
            sensitive_observable = let live = offset + 1, damaged = offset + 3
                (u, p, t) -> u[live] + u[damaged]
            end
            resistant_observable = let live = offset + 2, damaged = offset + 4
                (u, p, t) -> u[live] + u[damaged]
            end
            for (component, values, scale, observable) in (("sensitive", environment.sensitive, environment.sensitive_scale, sensitive_observable), ("resistant", environment.resistant, environment.resistant_scale, resistant_observable))
                push!(datasets, (x = environment.times, y = values, observable = observable, residual_scale = scale))
                push!(metadata, (density = environment.density, mix = environment.mix, component = component))
            end
        else
            append!(u0, [environment.u0_sensitive, 0.0, environment.u0_resistant])
            sensitive_observable = let sensitive = offset + 1, tolerant = offset + 2
                (u, p, t) -> u[sensitive] + u[tolerant]
            end
            resistant_observable = let resistant = offset + 3
                (u, p, t) -> u[resistant]
            end
            for (component, values, scale, observable) in (("sensitive", environment.sensitive, environment.sensitive_scale, sensitive_observable), ("resistant", environment.resistant, environment.resistant_scale, resistant_observable))
                push!(datasets, (x = environment.times, y = values, observable = observable, residual_scale = scale))
                push!(metadata, (density = environment.density, mix = environment.mix, component = component))
            end
        end
    end
    return datasets, u0, metadata
end

function _fit_treated_coculture_joint(environments, out; start, max_time_per_fit)
    untreated_baseline = _load_coculture_untreated_baseline(; start = start)
    baselines = [_coculture_monoculture_baselines(start, environment.density) for environment in environments]
    untreated_specs = _untreated_coculture_specs(environments, baselines, untreated_baseline.pooling_mode)
    haskey(untreated_specs, untreated_baseline.model) || error("Unknown inherited coculture model $(untreated_baseline.model)")
    t0 = 0.0
    rows = NamedTuple[]
    overlays = DataFrame[]
    parameter_rows = NamedTuple[]
    profile_parts = DataFrame[]
    identifiability_parts = DataFrame[]
    pooling_modes = _tournament_smoke_test() ? ("shared",) : ("shared", "partial_5pct", "independent_diagnostic")
    for pooling_mode in pooling_modes
        specs = _treated_coculture_specs(environments, baselines, untreated_baseline, untreated_specs, pooling_mode, t0)
        model_names = sort(collect(keys(specs)))
        _tournament_smoke_test() && (model_names = filter(==("dual_constant_kill"), model_names))
        if lowercase(get(ENV, "A2780_CONDITIONAL_TOURNAMENT", "false")) == "true"
            model_names = filter(
                name -> name ∉ ("dual_transit_competitor_scaled", "dual_transit_load_scaled"),
                model_names,
            )
        end
        for model_name in model_names
            spec = specs[model_name]
            datasets, u0, metadata = _treated_coculture_inputs(environments, spec.layout)
            println("  Joint treated coculture fit $(model_name), $(pooling_mode)...")
            explicit_profiles = Dict{Symbol,Vector{Float64}}(name => [5.0, 7.5, 10.0] for name in spec.names if occursin("lambda", String(name)))
            model_maxiters = max(650, Int(round(max_time_per_fit * 22)))
            multistart = try
                requested_starts = _configured_multistarts(1)
                screening_starts = lowercase(get(ENV, "A2780_CONDITIONAL_TOURNAMENT", "false")) == "true" ?
                    min(requested_starts, 3) : requested_starts
                starts = GrowthParameterEstimation.generate_multistarts(
                    spec.p0, spec.bounds;
                    n_starts = screening_starts,
                    seed = _stable_multistart_seed("treated_coculture", model_name, pooling_mode),
                )
                GrowthParameterEstimation.run_joint_multistart(
                    spec.model, datasets, u0, starts;
                    refine_optimizer = _configured_refine_optimizer(),
                    bounds = spec.bounds,
                    solver = Rodas5(),
                    maxiters = model_maxiters,
                    maxtime = _configured_fit_maxtime(),
                    reltol = 1e-7,
                    abstol = 1e-7,
                    optimizer = :nelder_mead,
                    initial_time = 0.0,
                )
            catch error_value
                @warn "Treated coculture joint multistart failed" model_name pooling_mode exception = error_value
                nothing
            end
            multistart === nothing && continue
            profiled = try
                GrowthParameterEstimation.profile_joint_fit_bounds(
                    spec.model, datasets, u0, multistart.fit.params;
                    bounds = spec.bounds,
                    parameter_names = spec.names,
                    explicit_upper_profiles = explicit_profiles,
                    solver = Rodas5(),
                    maxiters = model_maxiters,
                    reltol = 1e-7,
                    abstol = 1e-7,
                    optimizer = :nelder_mead,
                    initial_time = 0.0,
                )
            catch error_value
                @warn "Treated coculture joint fit failed" model_name pooling_mode exception = error_value
                nothing
            end
            profiled === nothing && continue
            fit = profiled.fit
            isfinite(fit.bic) && isfinite(fit.raw_sse) && fit.raw_sse < 9.99e11 || continue
            profile = DataFrame(profiled.profile)
            if !isempty(profile)
                profile.model = fill(model_name, nrow(profile))
                profile.pooling_mode = fill(pooling_mode, nrow(profile))
                push!(profile_parts, profile)
            end
            identifiability = DataFrame(profiled.identifiability)
            identifiability.model = fill(model_name, nrow(identifiability))
            identifiability.pooling_mode = fill(pooling_mode, nrow(identifiability))
            push!(identifiability_parts, identifiability)
            boundary_issue = any(identifiability.identifiability .!= "interior")
            push!(rows, (
                cell_line = COCULTURE_SYSTEM_LABEL,
                model = model_name,
                pooling_mode = pooling_mode,
                eligible_for_inheritance = pooling_mode != "independent_diagnostic",
                bic = Float64(fit.bic),
                ssr = Float64(fit.raw_sse),
                scaled_ssr = Float64(fit.scaled_sse),
                n_parameters = length(fit.params),
                n_points = sum(length(dataset.x) for dataset in datasets),
                n_densities = 2,
                n_mixes = length(unique(environment.mix for environment in environments)),
                untreated_interaction_model = untreated_baseline.model,
                untreated_interaction_pooling = untreated_baseline.pooling_mode,
                treatment_interaction_mode = model_name == "dual_transit_competitor_scaled" ? "multiplicative drug scaling by inherited opposing-lineage competition burden" : model_name == "dual_transit_load_scaled" ? "multiplicative drug scaling by inherited total competitive load" : "additive treatment effect on inherited growth and interaction dynamics",
                parameter_scope = pooling_mode == "partial_5pct" ? "shared treatment center with all kill amplitudes scaled together within +/-5% by density" : pooling_mode,
                initial_condition_strategy = "fixed nominal day-zero total split by mix",
                residual_scaling = "component-trajectory peak",
                boundary_issue = boundary_issue,
                params = string((names = spec.names, values = fit.params, package_api = "GrowthParameterEstimation.profile_joint_fit_bounds/run_joint_fit")),
            ))
            for (dataset_index, dataset) in enumerate(datasets)
                meta = metadata[dataset_index]
                environment_index = findfirst(environment -> environment.density == meta.density && environment.mix == meta.mix, environments)
                component_u0 = meta.component == "sensitive" ? environments[environment_index].u0_sensitive : environments[environment_index].u0_resistant
                overlay_length = length(dataset.x) + 1
                push!(overlays, DataFrame(
                    time = vcat(0.0, dataset.x),
                    observed = vcat(component_u0, dataset.y),
                    predicted = vcat(component_u0, fit.predictions[dataset_index]),
                    fixed_day0_anchor = vcat(true, fill(false, length(dataset.x))),
                    model = fill(model_name, overlay_length),
                    pooling_mode = fill(pooling_mode, overlay_length),
                    density = fill(meta.density, overlay_length),
                    mix = fill(meta.mix, overlay_length),
                    component = fill(meta.component, overlay_length),
                ))
            end
            for environment_index in eachindex(environments)
                local_values = spec.local_parameters(fit.params, environment_index)
                for parameter_index in eachindex(spec.base_names)
                    push!(parameter_rows, (
                        stage = "coculture_treated",
                        model = model_name,
                        pooling_mode = pooling_mode,
                        density = environments[environment_index].density,
                        parameter = String(spec.base_names[parameter_index]),
                        effective_value = Float64(local_values[parameter_index]),
                        lower_bound = Float64(spec.base_bounds[parameter_index][1]),
                        upper_bound = Float64(spec.base_bounds[parameter_index][2]),
                        bound_position = (Float64(local_values[parameter_index]) - spec.base_bounds[parameter_index][1]) / max(spec.base_bounds[parameter_index][2] - spec.base_bounds[parameter_index][1], eps()),
                    ))
                end
            end
        end
    end
    ranking = sort!(DataFrame(rows), :bic)
    isempty(ranking) && error("No finite treated coculture joint fits were produced")
    pooling_summary = GrowthParameterEstimation.summarize_pooling_bic(ranking; top_n = 5)
    eligible = ranking[Bool.(ranking.eligible_for_inheritance), :]
    winner = eligible[argmin(eligible.bic), :]
    overlay_df = vcat(overlays...; cols = :union)
    CSV.write(joinpath(out.csv, "coculture_treated_pooling_model_ranking.csv"), ranking)
    CSV.write(joinpath(out.csv, "coculture_treated_pooling_top5.csv"), pooling_summary.ranking)
    CSV.write(joinpath(out.csv, "coculture_treated_pooling_status.csv"), pooling_summary.status)
    CSV.write(joinpath(out.csv, "coculture_treated_joint_parameter_estimates.csv"), unique(DataFrame(parameter_rows)))
    nonadditive_comparison = _coculture_treated_nonadditive_comparison(ranking)
    CSV.write(joinpath(out.csv, "coculture_treated_nonadditive_model_comparison.csv"), nonadditive_comparison)
    interaction_audit = DataFrame(
        inherited_model = fill(untreated_baseline.model, length(untreated_baseline.params)),
        inherited_pooling_mode = fill(untreated_baseline.pooling_mode, length(untreated_baseline.params)),
        parameter = String.(untreated_baseline.names),
        inherited_value = Float64.(untreated_baseline.params),
        baseline_source = fill(untreated_baseline.path, length(untreated_baseline.params)),
        inheritance_mode = fill("exact_fixed_untreated_coculture_interaction_model_and_parameters", length(untreated_baseline.params)),
    )
    CSV.write(joinpath(out.csv, "coculture_treated_inheritance_audit.csv"), interaction_audit)
    CSV.write(joinpath(out.csv, "coculture_treated_boundary_profiles.csv"), isempty(profile_parts) ? DataFrame() : vcat(profile_parts...; cols = :union))
    CSV.write(joinpath(out.csv, "coculture_treated_identifiability.csv"), isempty(identifiability_parts) ? DataFrame() : vcat(identifiability_parts...; cols = :union))
    figure_csv = joinpath(out.csv, "figures")
    mkpath(figure_csv)
    CSV.write(joinpath(figure_csv, "coculture_treated_joint_overlays.csv"), overlay_df)
    _render_coculture_graph_grid(overlay_df, winner, out, "coculture_treated")
    _render_coculture_nonadditive_grid(overlay_df, nonadditive_comparison, out)
    CSV.write(joinpath(out.csv, "coculture_treated_unlinked_model_ranking.csv"), ranking)
    return _fit_linked_treatment_joint(environments, out; start = start, max_time_per_fit = max_time_per_fit)
end

function _coculture_treated_nonadditive_comparison(ranking::DataFrame)
    hypotheses = Dict(
        "dual_transit_damage" => "additive transit-damage reference; equivalent to beta=0",
        "dual_transit_competitor_scaled" => "drug damage multiplied by exp(beta * opposing-lineage competition burden)",
        "dual_transit_load_scaled" => "drug damage multiplied by exp(beta * total competitive load)",
    )
    rows = NamedTuple[]
    for model_name in keys(hypotheses)
        candidates = ranking[
            (String.(ranking.model) .== model_name) .& Bool.(ranking.eligible_for_inheritance),
            :,
        ]
        isempty(candidates) && continue
        winner = candidates[argmin(candidates.bic), :]
        push!(rows, (
            model = model_name,
            pooling_mode = String(winner.pooling_mode),
            bic = Float64(winner.bic),
            scaled_ssr = Float64(winner.scaled_ssr),
            ssr = Float64(winner.ssr),
            n_parameters = Int(winner.n_parameters),
            n_points = Int(winner.n_points),
            scaling_hypothesis = hypotheses[model_name],
        ))
    end
    isempty(rows) && return DataFrame(
        model = String[], pooling_mode = String[], bic = Float64[],
        scaled_ssr = Float64[], ssr = Float64[], n_parameters = Int[],
        n_points = Int[], scaling_hypothesis = String[], delta_bic = Float64[],
        bic_support = String[],
    )
    comparison = sort!(DataFrame(rows), :bic)
    isempty(comparison) && return comparison
    comparison.delta_bic = comparison.bic .- minimum(comparison.bic)
    comparison.bic_support = [delta < 2 ? "substantial" : delta < 6 ? "moderate" : delta < 10 ? "weak" : "not supported" for delta in comparison.delta_bic]
    return comparison
end

function _render_coculture_nonadditive_grid(overlay::DataFrame, comparison::DataFrame, out)
    isempty(comparison) && return nothing
    model_styles = Dict(
        "dual_transit_damage" => (:solid, "additive"),
        "dual_transit_competitor_scaled" => (:dash, "competitor-scaled"),
        "dual_transit_load_scaled" => (:dot, "load-scaled"),
    )
    panels = Any[]
    for density in ("20k", "30k"), mix in ("25-75", "50-50", "75-25")
        panel = plot(
            title = "$(density), mix $(mix)",
            titlefontsize = 10,
            xlabel = "Time (day)",
            ylabel = mix == "25-75" ? "Measured population" : "",
            legend = density == "20k" && mix == "75-25" ? :outertopright : false,
            legendfontsize = 7,
        )
        first_choice = first(comparison)
        observed_environment = overlay[
            (String.(overlay.model) .== String(first_choice.model)) .&
            (String.(overlay.pooling_mode) .== String(first_choice.pooling_mode)) .&
            (String.(overlay.density) .== density) .&
            (String.(overlay.mix) .== mix),
            :,
        ]
        for (component, color, short_component) in (("sensitive", :crimson, "S"), ("resistant", :steelblue, "R"))
            observed_rows = observed_environment[String.(observed_environment.component) .== component, :]
            sort!(observed_rows, :time)
            scatter!(panel, observed_rows.time, observed_rows.observed; color = color, ms = 3.2, markerstrokewidth = 0, alpha = 0.75, label = "$(short_component) observed")
            for choice in eachrow(comparison)
                model_name = String(choice.model)
                style, short_model = model_styles[model_name]
                predicted_rows = overlay[
                    (String.(overlay.model) .== model_name) .&
                    (String.(overlay.pooling_mode) .== String(choice.pooling_mode)) .&
                    (String.(overlay.density) .== density) .&
                    (String.(overlay.mix) .== mix) .&
                    (String.(overlay.component) .== component),
                    :,
                ]
                sort!(predicted_rows, :time)
                plot!(panel, predicted_rows.time, predicted_rows.predicted; color = color, linestyle = style, lw = 2.1, alpha = 0.9, label = "$(short_component) $(short_model)")
            end
        end
        push!(panels, panel)
    end
    figure = plot(
        panels...;
        layout = (2, 3),
        size = (1600, 860),
        plot_title = "Treated coculture: additive versus interaction-scaled transit simulations",
        plot_titlefontsize = 14,
        margin = 4 * Plots.mm,
    )
    image_dir = joinpath(out.images, "figures")
    mkpath(image_dir)
    path = joinpath(image_dir, "coculture_treated_nonadditive_simulation_grid.png")
    savefig(figure, path)
    return path
end

function _render_coculture_graph_grid(overlay::DataFrame, winner, out, condition)
    selected = overlay[
        (String.(overlay.model) .== String(winner.model)) .&
        (String.(overlay.pooling_mode) .== String(winner.pooling_mode)),
        :,
    ]
    panels = Any[]
    y_limits = _shared_y_limits(selected)
    for density in ("20k", "30k"), mix in ("25-75", "50-50", "75-25")
        environment = selected[(String.(selected.density) .== density) .& (String.(selected.mix) .== mix), :]
        isempty(environment) && continue
        panel = plot(
            title = "$(density), mix $(mix)",
            titlefontsize = 10,
            xlabel = "Time (day)",
            ylabel = mix == "25-75" ? "Measured population" : "",
            ylims = y_limits,
            legend = density == "20k" && mix == "75-25" ? :outertopright : false,
            legendfontsize = 8,
        )
        for (component, color) in (("sensitive", :crimson), ("resistant", :steelblue))
            component_rows = environment[String.(environment.component) .== component, :]
            sort!(component_rows, :time)
            scatter!(panel, component_rows.time, component_rows.observed; color = color, ms = 3.5, markerstrokewidth = 0, label = "$(component) observed")
            plot!(panel, component_rows.time, component_rows.predicted; color = color, lw = 2.5, label = "$(component) fit")
        end
        push!(panels, panel)
    end
    isempty(panels) && return nothing
    figure = plot(
        panels...;
        layout = (2, 3),
        size = (1500, 820),
        plot_title = "$(replace(condition, "_" => " ")): $(winner.model) | $(winner.pooling_mode), BIC=$(round(Float64(winner.bic); digits = 1))",
        plot_titlefontsize = 14,
        margin = 4 * Plots.mm,
    )
    image_dir = joinpath(out.images, "figures")
    mkpath(image_dir)
    path = joinpath(image_dir, "$(condition)_best_mechanistic_fit_grid.png")
    savefig(figure, path)
    return path
end

function _run_density_aware_coculture_fitting(decoded::DataFrame, condition::AbstractString, out; start, max_time_per_fit)
    recovered, provenance = _recover_coculture_design(decoded, condition)
    environments, initial_diagnostics = _coculture_environments(recovered)
    length(environments) == 6 || error("Expected six density-by-mix coculture environments; found $(length(environments))")
    CSV.write(joinpath(out.csv, "$(condition)_recovered_design.csv"), recovered)
    CSV.write(joinpath(out.csv, "$(condition)_density_provenance_validation.csv"), provenance)
    CSV.write(joinpath(out.csv, "$(condition)_initial_mix_diagnostics.csv"), initial_diagnostics)
    if condition == "coculture_untreated"
        return _fit_untreated_coculture_joint(environments, out; start = start, max_time_per_fit = max_time_per_fit)
    end
    return _fit_treated_coculture_joint(environments, out; start = start, max_time_per_fit = max_time_per_fit)
end
