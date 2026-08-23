# Joint timing-hypothesis comparison for treated monoculture.

const TREATED_TIMING_HYPOTHESES = [
    "independent_onset_gradual",
    "shared_onset_gradual",
    "partial_onset_0_5day",
    "cis_gradual_only",
    "cis_onset_only",
]

function _timing_parameter_seed(seed, name)
    values = seed.values
    mapping = Dict(
        :naive_emax => ("A2780Naive", :emax),
        :naive_ec50_effect => ("A2780Naive", :ec50_effect),
        :naive_hill_n => ("A2780Naive", :hill_n),
        :naive_lambda => ("A2780Naive", :lambda),
        :naive_log_contrast_density => ("A2780Naive", :log_contrast_effect),
        :cis_emax_sensitive => ("A2780cis", :emax_sensitive),
        :cis_emax_tolerant => ("A2780cis", :emax_tolerant),
        :cis_f_tolerant0 => ("A2780cis", :f_tolerant0),
        :cis_log_contrast_density => ("A2780cis", :log_contrast_effect),
        :naive_t_onset => ("A2780Naive", :t_onset),
        :cis_lambda => ("A2780cis", :lambda),
        :cis_t_onset => ("A2780cis", :t_onset),
    )
    return Float64(values[mapping[name]])
end

function _timing_hypothesis_spec(seed, hypothesis::String)
    hypothesis in TREATED_TIMING_HYPOTHESES || error("Unknown timing hypothesis $(hypothesis)")
    names = Symbol[
        :naive_emax, :naive_ec50_effect, :naive_hill_n, :naive_lambda,
        :naive_log_contrast_density, :cis_emax_sensitive, :cis_emax_tolerant,
        :cis_f_tolerant0, :cis_log_contrast_density,
    ]
    bounds = Tuple{Float64,Float64}[
        (0.0, 5.0), (0.05, 2.0), (0.2, 12.0), (0.01, 5.0),
        (-DENSITY_LOG_CONTRAST_BOUND, DENSITY_LOG_CONTRAST_BOUND),
        (0.0, 5.0), (0.0, 3.0), (0.0, 0.95),
        (-DENSITY_LOG_CONTRAST_BOUND, DENSITY_LOG_CONTRAST_BOUND),
    ]
    p0 = [_timing_parameter_seed(seed, name) for name in names]
    naive_onset = _timing_parameter_seed(seed, :naive_t_onset)
    cis_onset = _timing_parameter_seed(seed, :cis_t_onset)
    cis_lambda = _timing_parameter_seed(seed, :cis_lambda)

    if hypothesis == "independent_onset_gradual"
        append!(names, [:naive_t_onset, :cis_lambda, :cis_t_onset])
        append!(bounds, [(0.0, 7.0), (0.01, 5.0), (0.0, 7.0)])
        append!(p0, [naive_onset, cis_lambda, cis_onset])
    elseif hypothesis == "shared_onset_gradual"
        append!(names, [:shared_t_onset, :cis_lambda])
        append!(bounds, [(0.0, 7.0), (0.01, 5.0)])
        append!(p0, [0.5 * (naive_onset + cis_onset), cis_lambda])
    elseif hypothesis == "partial_onset_0_5day"
        append!(names, [:onset_center, :onset_delta, :cis_lambda])
        append!(bounds, [(0.5, 6.5), (-0.5, 0.5), (0.01, 5.0)])
        append!(p0, [
            clamp(0.5 * (naive_onset + cis_onset), 0.5, 6.5),
            clamp(0.5 * (cis_onset - naive_onset), -0.5, 0.5),
            cis_lambda,
        ])
    elseif hypothesis == "cis_gradual_only"
        append!(names, [:naive_t_onset, :cis_lambda])
        append!(bounds, [(0.0, 7.0), (0.01, 5.0)])
        append!(p0, [naive_onset, cis_lambda])
    elseif hypothesis == "cis_onset_only"
        append!(names, [:naive_t_onset, :cis_t_onset])
        append!(bounds, [(0.0, 7.0), (0.0, 7.0)])
        append!(p0, [naive_onset, cis_onset])
    end
    index = Dict(name => position for (position, name) in enumerate(names))
    return (hypothesis = hypothesis, names = names, bounds = bounds, p0 = p0, index = index)
end

_timing_value(p, spec, name, default = NaN) =
    haskey(spec.index, name) ? p[spec.index[name]] : default

function _timing_values(p, spec, density::String)
    sign = density == "20k" ? -1.0 : 1.0
    naive_scale = exp(sign * _timing_value(p, spec, :naive_log_contrast_density))
    cis_scale = exp(sign * _timing_value(p, spec, :cis_log_contrast_density))
    hypothesis = spec.hypothesis
    naive_onset, cis_onset = if hypothesis == "shared_onset_gradual"
        onset = _timing_value(p, spec, :shared_t_onset)
        (onset, onset)
    elseif hypothesis == "partial_onset_0_5day"
        center = _timing_value(p, spec, :onset_center)
        delta = _timing_value(p, spec, :onset_delta)
        (center - delta, center + delta)
    elseif hypothesis == "cis_gradual_only"
        (_timing_value(p, spec, :naive_t_onset), 0.0)
    else
        (_timing_value(p, spec, :naive_t_onset), _timing_value(p, spec, :cis_t_onset))
    end
    return (
        naive_emax = _timing_value(p, spec, :naive_emax) * naive_scale,
        naive_ec50 = _timing_value(p, spec, :naive_ec50_effect),
        naive_hill = _timing_value(p, spec, :naive_hill_n),
        naive_lambda = _timing_value(p, spec, :naive_lambda),
        naive_onset = naive_onset,
        cis_emax_sensitive = _timing_value(p, spec, :cis_emax_sensitive) * cis_scale,
        cis_emax_tolerant = _timing_value(p, spec, :cis_emax_tolerant) * cis_scale,
        cis_lambda = _timing_value(p, spec, :cis_lambda, NaN),
        cis_onset = cis_onset,
        cis_f_tolerant0 = clamp(_timing_value(p, spec, :cis_f_tolerant0), 0.0, 1.0),
        cis_activation_mode = hypothesis == "cis_onset_only" ? :step : :ramp,
    )
end

function _timing_activation(t, lambda, onset, mode)
    mode == :step && return t > onset ? one(t) : zero(t)
    return _ramp_activation(t, 0.0, lambda, onset)
end

function _timing_starts(spec)
    starts = [copy(spec.p0)]
    for (onset_guess, lambda_guess) in ((3.0, 0.7), (4.5, 2.0))
        candidate = copy(spec.p0)
        for name in (:naive_t_onset, :cis_t_onset, :shared_t_onset, :onset_center)
            haskey(spec.index, name) && (candidate[spec.index[name]] = onset_guess)
        end
        haskey(spec.index, :onset_delta) && (candidate[spec.index[:onset_delta]] = 0.0)
        haskey(spec.index, :naive_lambda) && (candidate[spec.index[:naive_lambda]] = lambda_guess)
        haskey(spec.index, :cis_lambda) && (candidate[spec.index[:cis_lambda]] = lambda_guess)
        push!(starts, clamp.(candidate, first.(spec.bounds), last.(spec.bounds)))
    end
    return starts
end

function _treated_timing_problem(start, seed, hypothesis)
    environments, source_path = _linked_monoculture_environments(start)
    spec = _timing_hypothesis_spec(seed, hypothesis)
    naive_environments = [environment for environment in environments if environment.cell_line == "A2780Naive"]
    cis_environments = [environment for environment in environments if environment.cell_line == "A2780cis"]
    n_naive = length(naive_environments)

    model! = function (du, u, p, t)
        for (index, environment) in enumerate(naive_environments)
            N = max(u[index], zero(u[index]))
            treatment = _timing_values(p, spec, environment.density)
            activation = _timing_activation(t, treatment.naive_lambda, treatment.naive_onset, :ramp)
            kill = activation * _hill_effect(
                environment.effect_level,
                treatment.naive_emax,
                treatment.naive_ec50,
                treatment.naive_hill,
            )
            du[index] = _linked_intrinsic_growth(N, environment.baseline, t) - kill * N
        end
        for (index, environment) in enumerate(cis_environments)
            sensitive_index = n_naive + 2index - 1
            tolerant_index = sensitive_index + 1
            sensitive = max(u[sensitive_index], zero(u[sensitive_index]))
            tolerant = max(u[tolerant_index], zero(u[tolerant_index]))
            total = sensitive + tolerant
            treatment = _timing_values(p, spec, environment.density)
            activation = _timing_activation(
                t,
                treatment.cis_lambda,
                treatment.cis_onset,
                treatment.cis_activation_mode,
            )
            kill_sensitive = activation * _hill_effect(
                environment.effect_level,
                treatment.cis_emax_sensitive,
                0.5,
                4.0,
            )
            kill_tolerant = activation * _hill_effect(
                environment.effect_level,
                treatment.cis_emax_tolerant,
                0.5,
                4.0,
            )
            growth = _linked_intrinsic_growth(total, environment.baseline, t)
            sensitive_share = total > 0 ? sensitive / total : one(total)
            tolerant_share = total > 0 ? tolerant / total : zero(total)
            du[sensitive_index] = growth * sensitive_share - kill_sensitive * sensitive
            du[tolerant_index] = growth * tolerant_share - kill_tolerant * tolerant
        end
    end

    u0_builder = function (p)
        values = eltype(p)[]
        for environment in naive_environments
            push!(values, environment.fixed_u0)
        end
        for environment in cis_environments
            treatment = _timing_values(p, spec, environment.density)
            push!(values, environment.fixed_u0 * (1 - treatment.cis_f_tolerant0))
            push!(values, environment.fixed_u0 * treatment.cis_f_tolerant0)
        end
        return values
    end

    datasets = NamedTuple[]
    metadata = NamedTuple[]
    for (index, environment) in enumerate(naive_environments)
        push!(datasets, (
            x = environment.times,
            y = environment.observed,
            state_index = index,
            residual_scale = environment.residual_scale,
        ))
        push!(metadata, (
            cell_line = environment.cell_line,
            density = environment.density,
            dose = environment.dose,
            ic_label = environment.ic_label,
            fixed_u0 = environment.fixed_u0,
        ))
    end
    for (index, environment) in enumerate(cis_environments)
        sensitive_index = n_naive + 2index - 1
        tolerant_index = sensitive_index + 1
        observable = let sensitive_index = sensitive_index, tolerant_index = tolerant_index
            (u, p, t) -> u[sensitive_index] + u[tolerant_index]
        end
        push!(datasets, (
            x = environment.times,
            y = environment.observed,
            observable = observable,
            residual_scale = environment.residual_scale,
        ))
        push!(metadata, (
            cell_line = environment.cell_line,
            density = environment.density,
            dose = environment.dose,
            ic_label = environment.ic_label,
            fixed_u0 = environment.fixed_u0,
        ))
    end
    return (
        model = model!,
        datasets = datasets,
        metadata = metadata,
        u0_builder = u0_builder,
        u0 = Float64.(u0_builder(spec.p0)),
        spec = spec,
        p0 = spec.p0,
        bounds = spec.bounds,
        names = spec.names,
        source_path = source_path,
    )
end

function _timing_overlay(fit, problem, hypothesis)
    parts = DataFrame[]
    for (index, dataset) in enumerate(problem.datasets)
        meta = problem.metadata[index]
        push!(parts, DataFrame(
            time = vcat(0.0, Float64.(dataset.x)),
            observed = vcat(Float64(meta.fixed_u0), Float64.(dataset.y)),
            predicted = vcat(Float64(meta.fixed_u0), Float64.(fit.predictions[index])),
            cell_line = fill(String(meta.cell_line), length(dataset.x) + 1),
            density = fill(String(meta.density), length(dataset.x) + 1),
            dose = fill(Float64(meta.dose), length(dataset.x) + 1),
            ic_label = fill(String(meta.ic_label), length(dataset.x) + 1),
            timing_hypothesis = fill(hypothesis, length(dataset.x) + 1),
        ))
    end
    return vcat(parts...; cols = :union)
end

function _render_timing_hypothesis_grid(overlay, ranking, out)
    top_models = String.(first(ranking.model, min(3, nrow(ranking))))
    colors = [:royalblue3, :darkorange2, :seagreen4]
    panels = Any[]
    selected_overlay = overlay[in.(String.(overlay.timing_hypothesis), Ref(top_models)), :]
    y_limits = _shared_y_limits(selected_overlay)
    for cell_line in ("A2780Naive", "A2780cis"), density in ("20k", "30k"), dose in (0.67, 1.0, 1.47)
        panel_rows = overlay[
            (String.(overlay.cell_line) .== cell_line) .&
            (String.(overlay.density) .== density) .&
            isapprox.(Float64.(overlay.dose), dose; atol = 0.02),
            :,
        ]
        label = _treated_dose_metadata(dose).ic_label
        panel = plot(
            title = "$(cell_line) $(density), $(label)",
            xlabel = "Time (day)",
            ylabel = "Cell count",
            ylims = y_limits,
            legend = false,
        )
        observed = panel_rows[String.(panel_rows.timing_hypothesis) .== first(top_models), :]
        scatter!(panel, observed.time, observed.observed; color = :black, markersize = 2.5, label = "observed")
        for (model_index, model_name) in enumerate(top_models)
            rows = panel_rows[String.(panel_rows.timing_hypothesis) .== model_name, :]
            sort!(rows, :time)
            plot!(panel, rows.time, rows.predicted; color = colors[model_index], linewidth = 1.8, label = model_name)
        end
        push!(panels, panel)
    end
    path = joinpath(out.images, "figures", "monoculture_treated_timing_hypothesis_grid.png")
    mkpath(dirname(path))
    savefig(plot(
        panels...;
        layout = (4, 3),
        size = (1500, 1750),
        plot_title = "Treated monoculture timing hypotheses",
    ), path)
    return path
end

function _fit_treated_timing_hypotheses(; start, max_time_per_fit = 1.5)
    out = IOUtils.condition_output_dirs("monoculture_treated"; start = start)
    seed = _linked_treatment_seed(start)
    ranking_rows = NamedTuple[]
    start_parts = DataFrame[]
    overlay_parts = DataFrame[]
    cache = Dict{String,Any}()

    for hypothesis in TREATED_TIMING_HYPOTHESES
        println("  Treated timing fit $(hypothesis)...")
        problem = _treated_timing_problem(start, seed, hypothesis)
        multistart = GrowthParameterEstimation.run_joint_multistart(
            problem.model,
            problem.datasets,
            problem.u0,
            _timing_starts(problem.spec);
            bounds = problem.bounds,
            u0_builder = problem.u0_builder,
            solver = Rodas5(),
            optimizer = :nelder_mead,
            maxiters = max(180, Int(round(max_time_per_fit * 70))),
            reltol = 1e-7,
            abstol = 1e-7,
            initial_time = 0.0,
        )
        fit = multistart.fit
        cache[hypothesis] = (fit = fit, problem = problem)
        starts = copy(multistart.summary)
        starts.timing_hypothesis = fill(hypothesis, nrow(starts))
        push!(start_parts, starts)
        timing = _timing_values(fit.params, problem.spec, "20k")
        push!(ranking_rows, (
            model = hypothesis,
            bic = Float64(fit.bic),
            scaled_ssr = Float64(fit.scaled_sse),
            ssr = Float64(fit.raw_sse),
            n_parameters = length(fit.params),
            n_points = sum(length(dataset.x) for dataset in problem.datasets),
            naive_onset = Float64(timing.naive_onset),
            cis_onset = Float64(timing.cis_onset),
            naive_lambda = Float64(timing.naive_lambda),
            cis_lambda = isfinite(timing.cis_lambda) ? Float64(timing.cis_lambda) : NaN,
            cis_activation_mode = String(timing.cis_activation_mode),
            best_start_index = multistart.best_start_index,
            boundary_issue = false,
            params = string((
                names = problem.names,
                values = fit.params,
                package_api = "GrowthParameterEstimation.run_joint_multistart/run_joint_fit",
            )),
        ))
    end

    ranking = sort!(DataFrame(ranking_rows), :bic)
    winner_name = String(first(ranking.model))
    winner = cache[winner_name]
    explicit_lower = Dict(name => [0.0] for name in winner.problem.names if occursin("lambda", String(name)))
    physical_lower = Dict(name => 0.0 for name in keys(explicit_lower))
    explicit_upper = Dict{Symbol,Vector{Float64}}()
    for name in winner.problem.names
        occursin("lambda", String(name)) && (explicit_upper[name] = [7.5, 10.0])
        occursin("emax", String(name)) && (explicit_upper[name] = [7.5, 10.0])
        occursin("onset", String(name)) && !occursin("delta", String(name)) &&
            (explicit_upper[name] = [10.0, 14.0])
    end
    profile_parameters = [
        name for name in winner.problem.names
        if !occursin("contrast", String(name)) && name != :onset_delta && name != :cis_f_tolerant0
    ]
    profiled = GrowthParameterEstimation.profile_joint_fit_bounds_two_sided(
        winner.problem.model,
        winner.problem.datasets,
        winner.problem.u0,
        winner.fit.params;
        bounds = winner.problem.bounds,
        parameter_names = winner.problem.names,
        explicit_lower_profiles = explicit_lower,
        physical_lower_limits = physical_lower,
        explicit_upper_profiles = explicit_upper,
        profile_parameters = profile_parameters,
        u0_builder = winner.problem.u0_builder,
        solver = Rodas5(),
        optimizer = :nelder_mead,
        maxiters = max(240, Int(round(max_time_per_fit * 80))),
        reltol = 1e-7,
        abstol = 1e-7,
        initial_time = 0.0,
    )
    cache[winner_name] = (
        fit = profiled.fit,
        problem = merge(winner.problem, (bounds = profiled.bounds,)),
    )
    winner_index = findfirst(String.(ranking.model) .== winner_name)
    winner_timing = _timing_values(profiled.fit.params, winner.problem.spec, "20k")
    ranking.bic[winner_index] = Float64(profiled.fit.bic)
    ranking.scaled_ssr[winner_index] = Float64(profiled.fit.scaled_sse)
    ranking.ssr[winner_index] = Float64(profiled.fit.raw_sse)
    ranking.naive_onset[winner_index] = Float64(winner_timing.naive_onset)
    ranking.cis_onset[winner_index] = Float64(winner_timing.cis_onset)
    ranking.naive_lambda[winner_index] = Float64(winner_timing.naive_lambda)
    ranking.cis_lambda[winner_index] =
        isfinite(winner_timing.cis_lambda) ? Float64(winner_timing.cis_lambda) : NaN
    ranking.boundary_issue[winner_index] =
        any(profiled.identifiability.identifiability .!= "interior")
    ranking.params[winner_index] = string((
        names = winner.problem.names,
        values = profiled.fit.params,
        package_api = "GrowthParameterEstimation.run_joint_multistart/profile_joint_fit_bounds_two_sided",
    ))
    sort!(ranking, :bic)
    ranking.delta_bic = ranking.bic .- minimum(ranking.bic)
    ranking.rank = collect(1:nrow(ranking))

    parameter_parts = DataFrame[]
    for hypothesis in TREATED_TIMING_HYPOTHESES
        cached = cache[hypothesis]
        push!(overlay_parts, _timing_overlay(cached.fit, cached.problem, hypothesis))
        push!(parameter_parts, DataFrame(
            timing_hypothesis = fill(hypothesis, length(cached.problem.names)),
            parameter = String.(cached.problem.names),
            estimate = Float64.(cached.fit.params),
            lower_bound = first.(cached.problem.bounds),
            upper_bound = last.(cached.problem.bounds),
        ))
    end
    overlay = vcat(overlay_parts...; cols = :union)
    parameters = vcat(parameter_parts...; cols = :union)
    starts = vcat(start_parts...; cols = :union)
    CSV.write(joinpath(out.csv, "monoculture_treated_timing_hypothesis_ranking.csv"), ranking)
    CSV.write(joinpath(out.csv, "monoculture_treated_timing_hypothesis_parameters.csv"), parameters)
    CSV.write(joinpath(out.csv, "monoculture_treated_timing_hypothesis_multistart.csv"), starts)
    CSV.write(joinpath(out.csv, "monoculture_treated_timing_hypothesis_boundary_profiles.csv"), profiled.profile)
    CSV.write(joinpath(out.csv, "monoculture_treated_timing_hypothesis_identifiability.csv"), profiled.identifiability)
    figure_csv = joinpath(out.csv, "figures")
    mkpath(figure_csv)
    CSV.write(joinpath(figure_csv, "monoculture_treated_timing_hypothesis_overlays.csv"), overlay)
    figure_path = _render_timing_hypothesis_grid(overlay, ranking, out)
    return (ranking = ranking, cache = cache, figure_path = figure_path)
end
