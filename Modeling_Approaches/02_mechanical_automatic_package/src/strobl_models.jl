# Literature benchmark models from Strobl et al., Cancer Research (2021).
#
# These models deliberately remain diagnostic. In the paper, S and R are
# drug-sensitive and fully resistant phenotypes. Here they are mapped to the
# separately measured A2780Naive and A2780cis lineages only to test whether the
# paper's low-parameter competition structures can explain the A2780 data.

const STROBL_REFERENCE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC8455086/"
const STROBL_SUPPLEMENT_URL = "https://people.maths.ox.ac.uk/maini/PKM%20publications/487-supp.pdf"

const STROBL_MODEL_VARIANTS = [
    (name = "strobl_birth_no_cost_no_turnover", cost_mode = :none, fit_cost = false, fit_turnover = false, density_dependent_death = false,
     paper_model = "2params_noCost_noTurnover"),
    (name = "strobl_birth_growth_cost_only", cost_mode = :growth, fit_cost = true, fit_turnover = false, density_dependent_death = false,
     paper_model = "3params_noTurnover"),
    (name = "strobl_birth_turnover_only", cost_mode = :none, fit_cost = false, fit_turnover = true, density_dependent_death = false,
     paper_model = "3params_noCost"),
    (name = "strobl_birth_growth_cost_turnover", cost_mode = :growth, fit_cost = true, fit_turnover = true, density_dependent_death = false,
     paper_model = "4params"),
    (name = "strobl_birth_capacity_cost_turnover", cost_mode = :capacity, fit_cost = true, fit_turnover = true, density_dependent_death = false,
     paper_model = "supplement_cost_in_KR"),
    (name = "strobl_birth_death_cost_turnover", cost_mode = :death, fit_cost = true, fit_turnover = true, density_dependent_death = false,
     paper_model = "supplement_cost_in_dR"),
    (name = "strobl_density_dependent_death", cost_mode = :growth, fit_cost = true, fit_turnover = true, density_dependent_death = true,
     paper_model = "supplement_density_dependent_death"),
]

_strobl_variant(name::AbstractString) = only(filter(variant -> variant.name == name, STROBL_MODEL_VARIANTS))

function _strobl_capacity_initialization(environments)
    by_density = Dict{String,Float64}()
    for density in ("20k", "30k")
        selected = filter(environment -> environment.density == density, environments)
        peaks = [maximum(environment.sensitive .+ environment.resistant) for environment in selected]
        by_density[density] = max(1.1 * maximum(peaks), 2.0 * maximum(environment.fixed_total_u0 for environment in selected))
    end
    p0 = [by_density["20k"], by_density["30k"]]
    bounds = [(0.5 * value, 4.0 * value) for value in p0]
    return p0, bounds
end

function _strobl_parameter_spec(environments, variant; treated::Bool)
    p0, bounds = _strobl_capacity_initialization(environments)
    names = Symbol[:K_20k, :K_30k]
    if variant.fit_cost
        push!(p0, 0.15); push!(bounds, (0.0, 0.8)); push!(names, :resistance_cost)
    end
    if variant.fit_turnover
        push!(p0, 0.10); push!(bounds, (0.0, 0.8)); push!(names, :turnover_fraction)
    end
    if treated
        push!(p0, 0.75); push!(bounds, (0.0, 2.0)); push!(names, :d_D)
    end
    return (p0 = p0, bounds = bounds, names = names)
end

function _strobl_parameter_values(p, names, variant, density, sensitive_baseline)
    index = Dict(name => position for (position, name) in enumerate(names))
    carrying_capacity = p[index[density == "20k" ? :K_20k : :K_30k]]
    cost = variant.fit_cost ? p[index[:resistance_cost]] : zero(eltype(p))
    turnover = variant.fit_turnover ? p[index[:turnover_fraction]] : zero(eltype(p))
    d_D = haskey(index, :d_D) ? p[index[:d_D]] : zero(eltype(p))
    r_sensitive = sensitive_baseline.r
    r_resistant = variant.cost_mode == :growth ? (one(cost) - cost) * r_sensitive : r_sensitive
    d_sensitive = turnover * r_sensitive
    d_resistant = variant.cost_mode == :death ? (one(cost) + cost) * d_sensitive : d_sensitive
    K_sensitive = carrying_capacity
    K_resistant = variant.cost_mode == :capacity ? (one(cost) - cost) * carrying_capacity : carrying_capacity
    return (
        r_sensitive = r_sensitive,
        r_resistant = r_resistant,
        d_sensitive = d_sensitive,
        d_resistant = d_resistant,
        K_sensitive = K_sensitive,
        K_resistant = K_resistant,
        d_D = d_D,
    )
end

function _strobl_rhs(sensitive, resistant, values, dose, density_dependent_death)
    total = sensitive + resistant
    # `dose` is D(t)/D_max, preserving Strobl's published 1 - 2*d_D*D/D_max term.
    sensitive_drug_factor = one(dose) - 2values.d_D * dose
    if density_dependent_death
        crowding = total / max(values.K_sensitive, 1e-8)
        return (
            values.r_sensitive * sensitive_drug_factor * sensitive - values.d_sensitive * crowding * sensitive,
            values.r_resistant * resistant - values.d_sensitive * crowding * resistant,
        )
    end
    return (
        values.r_sensitive * (one(total) - total / values.K_sensitive) * sensitive_drug_factor * sensitive - values.d_sensitive * sensitive,
        values.r_resistant * (one(total) - total / values.K_resistant) * resistant - values.d_resistant * resistant,
    )
end

function _strobl_multistarts(spec)
    low = [clamp(value * 0.8, first(bound), last(bound)) for (value, bound) in zip(spec.p0, spec.bounds)]
    high = [clamp(value * 1.2, first(bound), last(bound)) for (value, bound) in zip(spec.p0, spec.bounds)]
    return [spec.p0, low, high]
end

function _strobl_untreated_problem(environments, baselines, datasets, u0, variant)
    spec = _strobl_parameter_spec(environments, variant; treated = false)
    model! = function (du, u, p, t)
        for environment_index in eachindex(environments)
            sensitive_index = 2environment_index - 1
            resistant_index = sensitive_index + 1
            sensitive = max(u[sensitive_index], zero(u[sensitive_index]))
            resistant = max(u[resistant_index], zero(u[resistant_index]))
            values = _strobl_parameter_values(p, spec.names, variant, environments[environment_index].density, baselines[environment_index][1])
            du[sensitive_index], du[resistant_index] = _strobl_rhs(sensitive, resistant, values, zero(eltype(p)), variant.density_dependent_death)
        end
    end
    return merge(spec, (model = model!, datasets = datasets, u0 = u0))
end

function _fit_strobl_untreated_benchmarks(environments, baselines, datasets, u0, metadata; max_time_per_fit)
    rows = NamedTuple[]
    overlays = DataFrame[]
    parameter_rows = NamedTuple[]
    for variant in STROBL_MODEL_VARIANTS
        problem = _strobl_untreated_problem(environments, baselines, datasets, u0, variant)
        println("  Strobl untreated benchmark $(variant.name)...")
        multistart = GrowthParameterEstimation.run_joint_multistart(
            problem.model, problem.datasets, problem.u0, _strobl_multistarts(problem);
            bounds = problem.bounds,
            solver = Rodas5(),
            optimizer = :nelder_mead,
            maxiters = max(350, Int(round(max_time_per_fit * 45))),
            reltol = 1e-7,
            abstol = 1e-7,
            initial_time = 0.0,
        )
        fit = multistart.fit
        boundary_issue = any(
            abs(fit.params[index] - first(problem.bounds[index])) <= 0.02 * (last(problem.bounds[index]) - first(problem.bounds[index])) ||
            abs(last(problem.bounds[index]) - fit.params[index]) <= 0.02 * (last(problem.bounds[index]) - first(problem.bounds[index]))
            for index in eachindex(fit.params)
        )
        push!(rows, (
            cell_line = COCULTURE_SYSTEM_LABEL,
            model = variant.name,
            pooling_mode = "strobl_joint",
            eligible_for_inheritance = false,
            bic = Float64(fit.bic),
            ssr = Float64(fit.raw_sse),
            scaled_ssr = Float64(fit.scaled_sse),
            n_parameters = length(fit.params),
            n_points = sum(length(dataset.x) for dataset in datasets),
            n_densities = 2,
            n_mixes = length(unique(environment.mix for environment in environments)),
            parameter_scope = "Strobl literature benchmark; shared phenotype parameters with density-specific K",
            growth_inheritance = "A2780Naive Stage 1 r fixes the paper time scale; A2780cis is mapped to the paper resistant phenotype",
            initial_condition_strategy = "fixed experimental day-zero totals and nominal mixture fractions; paper n0/fR are therefore not refitted",
            residual_scaling = "component-trajectory peak",
            boundary_issue = boundary_issue,
            params = string((names = problem.names, values = fit.params, paper_model = variant.paper_model, reference = STROBL_REFERENCE_URL, package_api = "GrowthParameterEstimation.run_joint_multistart")),
        ))
        for (dataset_index, dataset) in enumerate(datasets)
            meta = metadata[dataset_index]
            fixed_u0 = meta.component == "sensitive" ? environments[meta.environment_index].u0_sensitive : environments[meta.environment_index].u0_resistant
            count = length(dataset.x) + 1
            push!(overlays, DataFrame(
                time = vcat(0.0, dataset.x), observed = vcat(fixed_u0, dataset.y), predicted = vcat(fixed_u0, fit.predictions[dataset_index]),
                fixed_day0_anchor = vcat(true, fill(false, length(dataset.x))), model = fill(variant.name, count),
                pooling_mode = fill("strobl_joint", count), density = fill(meta.density, count), mix = fill(meta.mix, count), component = fill(meta.component, count),
            ))
        end
        for (name, value, bound) in zip(problem.names, fit.params, problem.bounds)
            push!(parameter_rows, (
                stage = "coculture_untreated", model = variant.name, pooling_mode = "strobl_joint", density = "shared_or_named",
                parameter = String(name), effective_value = Float64(value), lower_bound = Float64(first(bound)), upper_bound = Float64(last(bound)),
                bound_position = (Float64(value) - first(bound)) / max(last(bound) - first(bound), eps()),
            ))
        end
    end
    return (rows = rows, overlays = overlays, parameter_rows = parameter_rows)
end

function _strobl_linked_problem(monoculture_environments, coculture_environments, start, variant)
    baselines = Dict(density => _coculture_monoculture_baselines(start, density) for density in ("20k", "30k"))
    spec = _strobl_parameter_spec(coculture_environments, variant; treated = true)
    naive_mono = [environment for environment in monoculture_environments if environment.cell_line == "A2780Naive"]
    cis_mono = [environment for environment in monoculture_environments if environment.cell_line == "A2780cis"]
    coculture_offset = length(naive_mono) + length(cis_mono)
    model! = function (du, u, p, t)
        for (index, environment) in enumerate(naive_mono)
            sensitive = max(u[index], zero(u[index]))
            values = _strobl_parameter_values(p, spec.names, variant, environment.density, baselines[environment.density][1])
            dose = environment.dose / 1.47
            du[index], _ = _strobl_rhs(sensitive, zero(sensitive), values, dose, variant.density_dependent_death)
        end
        for (index, environment) in enumerate(cis_mono)
            state_index = length(naive_mono) + index
            resistant = max(u[state_index], zero(u[state_index]))
            values = _strobl_parameter_values(p, spec.names, variant, environment.density, baselines[environment.density][1])
            _, du[state_index] = _strobl_rhs(zero(resistant), resistant, values, environment.dose / 1.47, variant.density_dependent_death)
        end
        for (index, environment) in enumerate(coculture_environments)
            sensitive_index = coculture_offset + 2index - 1
            resistant_index = sensitive_index + 1
            sensitive = max(u[sensitive_index], zero(u[sensitive_index]))
            resistant = max(u[resistant_index], zero(u[resistant_index]))
            values = _strobl_parameter_values(p, spec.names, variant, environment.density, baselines[environment.density][1])
            du[sensitive_index], du[resistant_index] = _strobl_rhs(sensitive, resistant, values, 1.0 / 1.47, variant.density_dependent_death)
        end
    end
    datasets = NamedTuple[]
    metadata = NamedTuple[]
    u0 = Float64[]
    for (index, environment) in enumerate(naive_mono)
        push!(u0, environment.fixed_u0)
        push!(datasets, (x = environment.times, y = environment.observed, state_index = index, residual_scale = environment.residual_scale))
        push!(metadata, (context = "monoculture", cell_line = environment.cell_line, density = environment.density, dose = environment.dose, mix = "", component = "total", fixed_u0 = environment.fixed_u0))
    end
    for (index, environment) in enumerate(cis_mono)
        push!(u0, environment.fixed_u0)
        state_index = length(naive_mono) + index
        push!(datasets, (x = environment.times, y = environment.observed, state_index = state_index, residual_scale = environment.residual_scale))
        push!(metadata, (context = "monoculture", cell_line = environment.cell_line, density = environment.density, dose = environment.dose, mix = "", component = "total", fixed_u0 = environment.fixed_u0))
    end
    for (index, environment) in enumerate(coculture_environments)
        append!(u0, [environment.u0_sensitive, environment.u0_resistant])
        sensitive_index = coculture_offset + 2index - 1
        resistant_index = sensitive_index + 1
        push!(datasets, (x = environment.times, y = environment.sensitive, state_index = sensitive_index, residual_scale = environment.sensitive_scale))
        push!(metadata, (context = "coculture", cell_line = COCULTURE_SYSTEM_LABEL, density = environment.density, dose = 1.0, mix = environment.mix, component = "sensitive", fixed_u0 = environment.u0_sensitive))
        push!(datasets, (x = environment.times, y = environment.resistant, state_index = resistant_index, residual_scale = environment.resistant_scale))
        push!(metadata, (context = "coculture", cell_line = COCULTURE_SYSTEM_LABEL, density = environment.density, dose = 1.0, mix = environment.mix, component = "resistant", fixed_u0 = environment.u0_resistant))
    end
    return merge(spec, (model = model!, datasets = datasets, metadata = metadata, u0 = u0))
end

function _fit_strobl_linked_benchmarks(monoculture_environments, coculture_environments, start; max_time_per_fit)
    rows = NamedTuple[]
    fits = Dict{String,Any}()
    overlays = DataFrame[]
    for variant in STROBL_MODEL_VARIANTS
        problem = _strobl_linked_problem(monoculture_environments, coculture_environments, start, variant)
        println("  Strobl linked-treatment benchmark $(variant.name)...")
        multistart = GrowthParameterEstimation.run_joint_multistart(
            problem.model, problem.datasets, problem.u0, _strobl_multistarts(problem);
            bounds = problem.bounds,
            solver = Rodas5(),
            optimizer = :nelder_mead,
            maxiters = max(350, Int(round(max_time_per_fit * 45))),
            reltol = 1e-7,
            abstol = 1e-7,
            initial_time = 0.0,
        )
        fit = multistart.fit
        boundary_issue = any(
            abs(fit.params[index] - first(problem.bounds[index])) <= 0.02 * (last(problem.bounds[index]) - first(problem.bounds[index])) ||
            abs(last(problem.bounds[index]) - fit.params[index]) <= 0.02 * (last(problem.bounds[index]) - first(problem.bounds[index]))
            for index in eachindex(fit.params)
        )
        fits[variant.name] = (fit = fit, problem = problem)
        push!(overlays, _linked_overlay(fit, problem, variant.name))
        push!(rows, (
            cell_line = COCULTURE_SYSTEM_LABEL,
            model = variant.name,
            pooling_mode = "strobl_joint",
            eligible_for_inheritance = false,
            diagnostic_model = true,
            bic = Float64(fit.bic),
            ssr = Float64(fit.raw_sse),
            scaled_ssr = Float64(fit.scaled_sse),
            n_parameters = length(fit.params),
            n_points = sum(length(dataset.x) for dataset in problem.datasets),
            n_densities = 2,
            n_mixes = 3,
            n_monoculture_trajectories = 12,
            n_coculture_component_trajectories = 12,
            inherited_timing_hypothesis = "none; immediate Strobl drug multiplier",
            parameter_scope = "Strobl literature benchmark jointly fitted to the same 24 trajectories",
            initial_condition_strategy = "fixed experimental day-zero totals and nominal mixture fractions; no latent A2780cis split",
            residual_scaling = "trajectory peak across combined 24-trajectory objective",
            boundary_issue = boundary_issue,
            params = string((names = problem.names, values = fit.params, paper_model = variant.paper_model, reference = STROBL_REFERENCE_URL, package_api = "GrowthParameterEstimation.run_joint_multistart")),
        ))
    end
    return (rows = rows, fits = fits, overlays = overlays)
end
