# Reduced inheritance-route comparison for treated A2780 coculture.
# Included inside FitWorkflows so it can reuse the staged data and model loaders.

const REDUCED_STAGE_MODELS = (
    "stage2_direct_context_scaling",
    "stage2_plus_asymmetric_competition",
    "stage2_plus_asymmetric_competition_death",
)

function _reduced_treatment_first_problem(environments, start, model_name)
    seed = _linked_treatment_seed(start)
    base = _linked_base_spec(seed)
    baselines = [_coculture_monoculture_baselines(start, environment.density) for environment in environments]

    if model_name == "stage2_direct_context_scaling"
        p0 = [0.0, 0.0]
        bounds = [(-3.0, 3.0), (-3.0, 3.0)]
        names = [:log_naive_context_scale, :log_cis_context_scale]
    elseif model_name == "stage2_plus_asymmetric_competition"
        p0 = [0.5, 1.0]
        bounds = [(0.0, 4.0), (0.0, 4.0)]
        names = [:alpha_NC, :alpha_CN]
    elseif model_name == "stage2_plus_asymmetric_competition_death"
        p0 = [0.5, 1.0, 0.05, 0.01]
        bounds = [(0.0, 4.0), (0.0, 4.0), (0.0, 1.0), (0.0, 1.0)]
        names = [:alpha_NC, :alpha_CN, :death_naive, :death_cis]
    else
        error("Unknown reduced-stage model $(model_name)")
    end

    states_per_environment = 3
    model! = function (du, u, p, t)
        for (environment_index, environment) in enumerate(environments)
            offset = (environment_index - 1) * states_per_environment
            naive = max(u[offset + 1], zero(u[offset + 1]))
            cis_sensitive = max(u[offset + 2], zero(u[offset + 2]))
            cis_tolerant = max(u[offset + 3], zero(u[offset + 3]))
            cis_total = cis_sensitive + cis_tolerant
            naive_baseline, cis_baseline = baselines[environment_index]

            if model_name == "stage2_direct_context_scaling"
                alpha_NC, alpha_CN, death_naive, death_cis = 0.0, 0.0, 0.0, 0.0
                naive_context_scale = exp(clamp(p[1], -3.0, 3.0))
                cis_context_scale = exp(clamp(p[2], -3.0, 3.0))
            elseif model_name == "stage2_plus_asymmetric_competition"
                alpha_NC, alpha_CN = p
                death_naive, death_cis = 0.0, 0.0
                naive_context_scale, cis_context_scale = 1.0, 1.0
            else
                alpha_NC, alpha_CN, death_naive, death_cis = p
                naive_context_scale, cis_context_scale = 1.0, 1.0
            end

            naive_growth = _anchored_component_growth(
                naive, naive + alpha_NC * cis_total, naive_baseline, t,
            ) - death_naive * naive
            cis_growth = _anchored_component_growth(
                cis_total, cis_total + alpha_CN * naive, cis_baseline, t,
            ) - death_cis * cis_total

            treatment = _linked_treatment_values(
                base.p0, "strict_inheritance", :coculture, environment.density, base,
            )
            naive_activation = _ramp_activation(t, 0.0, treatment.naive_lambda, treatment.naive_onset)
            cis_activation = _timing_activation(
                t, treatment.cis_lambda, treatment.cis_onset, treatment.cis_activation_mode,
            )
            naive_kill = naive_context_scale * naive_activation *
                _hill_effect(0.5, treatment.naive_emax, treatment.naive_ec50, treatment.naive_hill)
            cis_sensitive_kill = cis_context_scale * cis_activation *
                _hill_effect(0.5, treatment.cis_emax_sensitive, 0.5, 4.0)
            cis_tolerant_kill = cis_context_scale * cis_activation *
                _hill_effect(0.5, treatment.cis_emax_tolerant, 0.5, 4.0)
            sensitive_share = cis_total > 0 ? cis_sensitive / cis_total : one(cis_total)
            tolerant_share = cis_total > 0 ? cis_tolerant / cis_total : zero(cis_total)

            du[offset + 1] = naive_growth - naive_kill * naive
            du[offset + 2] = cis_growth * sensitive_share - cis_sensitive_kill * cis_sensitive
            du[offset + 3] = cis_growth * tolerant_share - cis_tolerant_kill * cis_tolerant
        end
    end

    u0 = Float64[]
    datasets = NamedTuple[]
    metadata = NamedTuple[]
    for (environment_index, environment) in enumerate(environments)
        treatment = _linked_treatment_values(
            base.p0, "strict_inheritance", :coculture, environment.density, base,
        )
        append!(u0, [
            environment.u0_sensitive,
            environment.u0_resistant * (1 - treatment.cis_f_tolerant0),
            environment.u0_resistant * treatment.cis_f_tolerant0,
        ])
        offset = (environment_index - 1) * states_per_environment
        cis_observable = let sensitive_index = offset + 2, tolerant_index = offset + 3
            (state, parameters, time) -> state[sensitive_index] + state[tolerant_index]
        end
        push!(datasets, (
            x = environment.times,
            y = environment.sensitive,
            state_index = offset + 1,
            residual_scale = environment.sensitive_scale,
        ))
        push!(metadata, (density = environment.density, mix = environment.mix, component = "sensitive"))
        push!(datasets, (
            x = environment.times,
            y = environment.resistant,
            observable = cis_observable,
            residual_scale = environment.resistant_scale,
        ))
        push!(metadata, (density = environment.density, mix = environment.mix, component = "resistant"))
    end

    starts = Vector{Vector{Float64}}([copy(p0)])
    midpoint = [(lower + upper) / 2 for (lower, upper) in bounds]
    push!(starts, midpoint)
    if length(p0) == 2
        push!(starts, [0.15 * first(bound) + 0.85 * last(bound) for bound in bounds])
    else
        push!(starts, [0.75 * first(bound) + 0.25 * last(bound) for bound in bounds])
    end
    return (
        model = model!, datasets = datasets, metadata = metadata, u0 = u0,
        p0 = p0, starts = starts, bounds = bounds, names = names,
        seed = seed, base = base,
    )
end

function _reduced_overlay(fit, problem, environments, model_name)
    parts = DataFrame[]
    for (dataset_index, dataset) in enumerate(problem.datasets)
        metadata = problem.metadata[dataset_index]
        environment = only(filter(
            item -> item.density == metadata.density && item.mix == metadata.mix,
            environments,
        ))
        fixed_u0 = metadata.component == "sensitive" ? environment.u0_sensitive : environment.u0_resistant
        count = length(dataset.x) + 1
        push!(parts, DataFrame(
            time = vcat(0.0, dataset.x),
            observed = vcat(fixed_u0, dataset.y),
            predicted = vcat(fixed_u0, fit.predictions[dataset_index]),
            fixed_day0_anchor = vcat(true, fill(false, length(dataset.x))),
            model = fill(model_name, count),
            route = fill("treatment_first", count),
            density = fill(metadata.density, count),
            mix = fill(metadata.mix, count),
            component = fill(metadata.component, count),
        ))
    end
    return vcat(parts...)
end

"""
    run_reduced_stage_comparison!(; start=pwd(), maxiters=1200)

Fit the reciprocal Stage 1 -> Stage 2 -> Stage 4 conditional models to treated
coculture only, and combine their ranking with the existing Stage 1 -> Stage 3
-> Stage 4 fits. All BIC values in the returned comparison use the same 168
treated-coculture observations and trajectory-specific residual scaling.
"""
function run_reduced_stage_comparison!(; start::AbstractString = pwd(), maxiters::Int = 1200)
    root = IOUtils.package_root(start)
    treated_dir = joinpath(root, "outputs", "csv", "coculture_treated")
    recovered_path = joinpath(treated_dir, "coculture_treated_recovered_design.csv")
    isfile(recovered_path) || error("Recovered treated-coculture design is required: $(recovered_path)")
    environments, _ = _coculture_environments(CSV.read(recovered_path, DataFrame))

    comparison_dir = joinpath(root, "outputs", "csv", "reduced_stage_comparison")
    mkpath(comparison_dir)
    ranking_rows = NamedTuple[]
    overlays = DataFrame[]
    multistart_parts = DataFrame[]
    identifiability_parts = DataFrame[]

    for model_name in REDUCED_STAGE_MODELS
        problem = _reduced_treatment_first_problem(environments, start, model_name)
        result = GrowthParameterEstimation.run_joint_multistart(
            problem.model, problem.datasets, problem.u0, problem.starts;
            bounds = problem.bounds,
            solver = Rodas5(),
            optimizer = :nelder_mead,
            maxiters = maxiters,
            reltol = 1e-7,
            abstol = 1e-7,
            initial_time = 0.0,
            show_stats = false,
        )
        fit = result.fit
        summary = DataFrame(result.summary)
        summary.model = fill(model_name, nrow(summary))
        push!(multistart_parts, summary)
        stability_rows = NamedTuple[]
        for candidate in result.fits
            candidate === nothing && continue
            isfinite(candidate.bic) || continue
            for parameter_index in eachindex(problem.names)
                push!(stability_rows, (
                    model = model_name,
                    parameter = String(problem.names[parameter_index]),
                    value = Float64(candidate.params[parameter_index]),
                    lower_bound = Float64(first(problem.bounds[parameter_index])),
                    upper_bound = Float64(last(problem.bounds[parameter_index])),
                ))
            end
        end
        identifiability = GrowthParameterEstimation.summarize_joint_parameter_stability(
            DataFrame(stability_rows),
        )
        identifiability.model = fill(model_name, nrow(identifiability))
        push!(identifiability_parts, identifiability)
        push!(ranking_rows, (
            route = "Stage 1 -> Stage 2 -> Stage 4",
            model = model_name,
            pooling_mode = "shared",
            bic = Float64(fit.bic),
            raw_sse = Float64(fit.raw_sse),
            scaled_sse = Float64(fit.scaled_sse),
            n_parameters = length(fit.params),
            n_points = sum(length(dataset.x) for dataset in problem.datasets),
            inherited_stage = "Stage 2 treatment response fixed from treated monoculture",
            fitted_stage = model_name == "stage2_direct_context_scaling" ?
                "Stage 4 lineage-specific coculture treatment scales" :
                "Stage 4 asymmetric competition terms",
            boundary_issue = any(abs(fit.params[index] - first(problem.bounds[index])) < 1e-5 || abs(fit.params[index] - last(problem.bounds[index])) < 1e-5 for index in eachindex(fit.params)),
            params = string((names = problem.names, values = fit.params, package_api = "GrowthParameterEstimation.run_joint_multistart")),
        ))
        push!(overlays, _reduced_overlay(fit, problem, environments, model_name))
    end

    route_a_path = joinpath(treated_dir, "coculture_treated_unlinked_model_ranking.csv")
    route_a_overlay_path = joinpath(treated_dir, "figures", "coculture_treated_joint_overlays.csv")
    isfile(route_a_path) || error("Stage 1 -> Stage 3 -> Stage 4 ranking is required")
    isfile(route_a_overlay_path) || error("Stage 1 -> Stage 3 -> Stage 4 overlays are required")
    route_a = CSV.read(route_a_path, DataFrame)
    route_a = route_a[Bool.(route_a.eligible_for_inheritance), :]
    route_a = sort!(route_a, :bic)
    for row in eachrow(route_a)
        push!(ranking_rows, (
            route = "Stage 1 -> Stage 3 -> Stage 4",
            model = String(row.model),
            pooling_mode = String(row.pooling_mode),
            bic = Float64(row.bic),
            raw_sse = Float64(row.ssr),
            scaled_sse = Float64(row.scaled_ssr),
            n_parameters = Int(row.n_parameters),
            n_points = Int(row.n_points),
            inherited_stage = "Stage 3 asymmetric competition with lineage-specific loss fixed from untreated coculture",
            fitted_stage = "Stage 4 treatment-response terms",
            boundary_issue = Bool(row.boundary_issue),
            params = String(row.params),
        ))
    end

    ranking = sort!(DataFrame(ranking_rows), :bic)
    ranking.delta_bic = ranking.bic .- minimum(ranking.bic)
    ranking.rank = collect(1:nrow(ranking))
    CSV.write(joinpath(comparison_dir, "reduced_stage_model_ranking.csv"), ranking)
    CSV.write(joinpath(comparison_dir, "treatment_first_multistart.csv"), vcat(multistart_parts...; cols = :union))
    CSV.write(joinpath(comparison_dir, "treatment_first_parameter_stability.csv"), vcat(identifiability_parts...; cols = :union))
    CSV.write(joinpath(comparison_dir, "treatment_first_overlays.csv"), vcat(overlays...; cols = :union))

    audit = DataFrame(
        route = ["Stage 1 -> Stage 3 -> Stage 4", "Stage 1 -> Stage 2 -> Stage 4"],
        omitted_stage = ["Stage 2 treated monoculture", "Stage 3 untreated coculture"],
        inherited_information = [
            "Stage 1 growth plus fixed Stage 3 asymmetric competition/loss",
            "Stage 1 growth plus fixed Stage 2 timing, Hill response, and cis state fraction",
        ],
        newly_fitted_information = [
            "Treatment terms from treated coculture",
            "Coculture context scale or asymmetric competition from treated coculture",
        ],
        comparison_observations = [168, 168],
        residual_scaling = ["component-trajectory peak", "component-trajectory peak"],
    )
    CSV.write(joinpath(comparison_dir, "inheritance_route_audit.csv"), audit)
    return (ranking = ranking, comparison_dir = comparison_dir)
end
