module A2780AdaptiveAdapter

using CSV
using DataFrames
using GrowthParameterEstimation
using JSON3

using ..AdaptiveSimulationEngine

export RankedCandidate,
       load_a2780_adaptive_config,
       write_a2780_adaptive_config,
       build_a2780_scenario,
       a2780_model_spec,
       a2780_segment_simulator

struct RankedCandidate
    stage::Int
    rank::Int
    model::String
    pooling::String
    bic::Float64
    parameters::Dict{String,Float64}
    boundary_issue::Bool
    identifiability::String
    evidence_level::String
    eligible::Bool
    exclusion_reason::String
end

function _untreated_coculture_library(csv_root, winner::RankedCandidate)
    overlay_path = joinpath(csv_root, "coculture_untreated", "figures", "coculture_untreated_joint_overlays.csv")
    overlay = CSV.read(overlay_path, DataFrame)
    selected = filter(row ->
        String(row.model) == winner.model && String(row.pooling_mode) == winner.pooling,
        overlay,
    )
    isempty(selected) && error("No untreated coculture overlays match the selected Stage 3 model")

    original = Dict{String,Any}()
    fitted = Dict{String,Any}()
    point_rows(rows, column) = [
        (t = Float64(row.time), y = Float64(row[column]))
        for row in eachrow(sort(rows, :time))
    ]
    for density in sort(unique(String.(selected.density)))
        density_rows = filter(row -> String(row.density) == density, selected)
        for mix in sort(unique(String.(density_rows.mix)))
            condition = filter(row -> String(row.mix) == mix, density_rows)
            naive = filter(row -> String(row.component) == "sensitive", condition)
            cis = filter(row -> String(row.component) == "resistant", condition)
            isempty(naive) && error("Missing untreated Naive overlay for $density, mix $mix")
            isempty(cis) && error("Missing untreated cis overlay for $density, mix $mix")
            key = "coculture_untreated__$(density)__$(mix)__0"
            original[key] = (
                key = key,
                label = "Untreated coculture | $density, mix $mix",
                category = "Untreated coculture",
                condition = "coculture_untreated",
                culture = "coculture",
                treated = false,
                density = density,
                mix = mix,
                cellLine = "A2780Naive+A2780cis",
                dose = 0.0,
                components = ["naive", "cis"],
                naive = point_rows(naive, :observed),
                cis = point_rows(cis, :observed),
            )
            fitted[key] = (
                model = "Asymmetric competition with loss and inherited monoculture growth",
                naive = point_rows(naive, :predicted),
                cis = point_rows(cis, :predicted),
            )
        end
    end
    return (original = original, fitted = fitted)
end

_bool(value) = lowercase(strip(string(value))) == "true"

function _parse_parameters(text)
    source = String(text)
    names_match = match(r"names\s*=\s*\[([^\]]*)\]", source)
    values_match = match(r"values\s*=\s*\[([^\]]*)\]", source)
    names_match === nothing && return Dict{String,Float64}()
    values_match === nothing && return Dict{String,Float64}()
    names = [replace(strip(name), ":" => "") for name in split(names_match.captures[1], ',')]
    values = tryparse.(Float64, strip.(split(values_match.captures[1], ',')))
    length(names) == length(values) || return Dict{String,Float64}()
    any(isnothing, values) && return Dict{String,Float64}()
    return Dict(name => value::Float64 for (name, value) in zip(names, values))
end

function _parse_parameter_vector(text)
    source = strip(String(text))
    startswith(source, "[") && endswith(source, "]") || return Float64[]
    values = tryparse.(Float64, strip.(split(source[2:(end - 1)], ',')))
    any(isnothing, values) && return Float64[]
    return Float64[value::Float64 for value in values]
end

function _stage1_baselines(model, pooling, values)
    nbase = model in ("logistic_growth", "gompertz_growth") ? 2 : 3
    model in ("lagged_theta_logistic_growth", "baranyi_theta_logistic_growth", "adaptation_theta_logistic_growth") && (nbase = 4)
    length(values) >= nbase || return NamedTuple[]
    base = values[1:nbase]
    function effective(density)
        current = copy(base)
        if pooling == "partial_5pct"
            length(values) >= nbase + 2 || return Float64[]
            sign = density == "20k" ? -1.0 : 1.0
            current[1] *= exp(sign * values[nbase + 1])
            current[2] *= exp(sign * values[nbase + 2])
        elseif pooling != "shared"
            return Float64[]
        end
        return current
    end
    rows = NamedTuple[]
    for density in ("20k", "30k")
        current = effective(density)
        isempty(current) && return NamedTuple[]
        theta = model in ("logistic_growth", "gompertz_growth") ? 1.0 : current[3]
        push!(rows, (
            density = density,
            model = model,
            pooling = pooling,
            r = current[1],
            K = current[2],
            theta = theta,
            lag_time = model == "lagged_theta_logistic_growth" ? current[4] : 0.0,
            baranyi_q0 = model == "baranyi_theta_logistic_growth" ? current[4] : 1.0,
            adaptation_rate = model == "adaptation_theta_logistic_growth" ? current[4] : 1.0,
        ))
    end
    return rows
end

function _stage1_candidates(csv_root)
    path = joinpath(csv_root, "monoculture_untreated", "monoculture_untreated_automatic_model_ranking.csv")
    isfile(path) || error("required Stage 1 ranking export is missing: $path")
    table = CSV.read(path, DataFrame)
    candidates = NamedTuple[]
    for cell_line in unique(String.(table.cell_line))
        rows = sort(filter(row -> String(row.cell_line) == cell_line, table), :bic)
        for (rank, row) in enumerate(eachrow(rows))
            values = _parse_parameter_vector(row.params)
            model, pooling = String(row.model), String(row.pooling_mode)
            baselines = _stage1_baselines(model, pooling, values)
            eligible = _bool(row.eligible_for_inheritance) && !_bool(row.diagnostic_model) &&
                !isempty(baselines) && isfinite(Float64(row.bic))
            push!(candidates, (
                stage = 1,
                rank = rank,
                cell_line = cell_line,
                model = model,
                pooling = pooling,
                bic = Float64(row.bic),
                baselines = baselines,
                boundary_issue = _bool(row.boundary_issue),
                identifiability = "see Stage 1 profile exports",
                evidence_level = "A2780 fitted growth",
                eligible = eligible,
                exclusion_reason = eligible ? "" : "not eligible for inherited simulation",
            ))
        end
    end
    return candidates
end

function _identifiability_summary(path, parameters)
    isfile(path) || return "not exported"
    table = CSV.read(path, DataFrame)
    :parameter in propertynames(table) || return "not exported"
    rows = filter(row -> String(row.parameter) in keys(parameters), table)
    isempty(rows) && return "not exported"
    labels = unique(String(row.identifiability) for row in eachrow(rows) if !ismissing(row.identifiability))
    isempty(labels) && return "not exported"
    any(contains("poorly_identified"), labels) && return "boundary-limited"
    return join(sort(labels), ", ")
end

function _ranked_candidates(stage, ranking_path; identifiability_path = "", status_allowed = true, evidence_level = "data-fitted")
    isfile(ranking_path) || error("required Stage $stage ranking export is missing: $ranking_path")
    table = CSV.read(ranking_path, DataFrame)
    candidates = RankedCandidate[]
    for (index, row) in enumerate(eachrow(table))
        parameters = :params in propertynames(table) ? _parse_parameters(row.params) : Dict{String,Float64}()
        finite_parameters = !isempty(parameters) && all(isfinite, values(parameters))
        finite_score = :bic in propertynames(table) && !ismissing(row.bic) && isfinite(Float64(row.bic))
        inherited = :eligible_for_inheritance in propertynames(table) ? _bool(row.eligible_for_inheritance) : true
        diagnostic = :diagnostic_model in propertynames(table) ? _bool(row.diagnostic_model) : false
        boundary = :boundary_issue in propertynames(table) ? _bool(row.boundary_issue) : false
        eligible = status_allowed && finite_parameters && finite_score && inherited && !diagnostic
        reasons = String[]
        status_allowed || push!(reasons, "stage status disallows inheritance")
        finite_score || push!(reasons, "non-finite ranking score")
        finite_parameters || push!(reasons, "missing or non-finite parameter export")
        inherited || push!(reasons, "not eligible for inheritance")
        diagnostic && push!(reasons, "diagnostic candidate")
        pooling_column = :pooling_mode in propertynames(table) ? :pooling_mode : (:pooling in propertynames(table) ? :pooling : nothing)
        model = String(row.model)
        push!(candidates, RankedCandidate(
            stage,
            :rank in propertynames(table) ? Int(row.rank) : index,
            model,
            pooling_column === nothing ? "" : String(row[pooling_column]),
            finite_score ? Float64(row.bic) : Inf,
            parameters,
            boundary,
            _identifiability_summary(identifiability_path, parameters),
            evidence_level,
            eligible,
            join(reasons, "; "),
        ))
    end
    return candidates
end

function _winner(candidates, stage)
    index = findfirst(candidate -> candidate.eligible, candidates)
    index === nothing && error("Stage $stage has no eligible exported model")
    return candidates[index], candidates[1:(index - 1)]
end

function _stage2_candidates(csv_root)
    path = joinpath(csv_root, "monoculture_treated", "monoculture_treated_joint_dose_model_ranking.csv")
    candidates = _ranked_candidates(2, path;
        identifiability_path = joinpath(csv_root, "monoculture_treated", "monoculture_treated_identifiability.csv"),
        evidence_level = "A2780 fitted; measured-dose range",
    )
    return candidates
end

function _baseline_rows(csv_root)
    path = joinpath(csv_root, "monoculture_untreated", "untreated_group_baselines.csv")
    isfile(path) || error("required Stage 1 baseline export is missing: $path")
    rows = CSV.read(path, DataFrame)
    all(row -> _bool(row.inheritance_allowed), eachrow(rows)) || error("Stage 1 baseline status disallows inheritance")
    for column in (:r, :K)
        all(row -> !ismissing(row[column]) && isfinite(Float64(row[column])), eachrow(rows)) ||
            error("Stage 1 contains non-finite $column values")
    end
    return rows
end

function _status_allowed(path, column)
    isfile(path) || return false
    table = CSV.read(path, DataFrame)
    nrow(table) > 0 || return false
    column in propertynames(table) || return false
    return _bool(table[1, column])
end

function _candidate_json(candidate)
    return (
        stage = candidate.stage,
        rank = candidate.rank,
        model = candidate.model,
        pooling = candidate.pooling,
        bic = candidate.bic,
        parameters = candidate.parameters,
        boundary_issue = candidate.boundary_issue,
        identifiability = candidate.identifiability,
        evidence_level = candidate.evidence_level,
        eligible = candidate.eligible,
        exclusion_reason = candidate.exclusion_reason,
    )
end

function load_a2780_adaptive_config(package_root)
    csv_root = joinpath(package_root, "outputs", "csv")
    baselines = _baseline_rows(csv_root)
    stage1 = _stage1_candidates(csv_root)

    stage2 = _stage2_candidates(csv_root)
    stage2_winner, stage2_skipped = _winner(stage2, 2)

    stage3 = _ranked_candidates(3,
        joinpath(csv_root, "coculture_untreated", "coculture_untreated_automatic_model_ranking.csv");
        identifiability_path = joinpath(csv_root, "coculture_untreated", "coculture_untreated_identifiability.csv"),
        evidence_level = "A2780 fitted competition",
    )
    stage3_winner, stage3_skipped = _winner(stage3, 3)

    stage4_status = joinpath(csv_root, "coculture_treated", "linked_treatment_status.csv")
    stage4 = _ranked_candidates(4,
        joinpath(csv_root, "coculture_treated", "linked_treatment_model_ranking.csv");
        identifiability_path = joinpath(csv_root, "coculture_treated", "linked_treatment_identifiability.csv"),
        status_allowed = _status_allowed(stage4_status, :inheritance_allowed),
        evidence_level = "A2780 linked fit",
    )
    stage4_winner, stage4_skipped = _winner(stage4, 4)

    baseline_json = [(
        cell_line = String(row.cell_line),
        density = String(row.density),
        model = String(row.best_model),
        pooling = String(row.pooling_mode),
        r = Float64(row.r),
        K = Float64(row.K),
        shape_parameter = ismissing(row.shape_parameter) ? nothing : String(row.shape_parameter),
        shape_value = ismissing(row.shape_value) || !isfinite(Float64(row.shape_value)) ? nothing : Float64(row.shape_value),
        theta = :theta in propertynames(baselines) && !ismissing(row.theta) && isfinite(Float64(row.theta)) ? Float64(row.theta) : nothing,
        lag_time = :lag_time in propertynames(baselines) && !ismissing(row.lag_time) && isfinite(Float64(row.lag_time)) ? Float64(row.lag_time) : nothing,
        baranyi_q0 = :baranyi_q0 in propertynames(baselines) && !ismissing(row.baranyi_q0) && isfinite(Float64(row.baranyi_q0)) ? Float64(row.baranyi_q0) : nothing,
        adaptation_rate = :adaptation_rate in propertynames(baselines) && !ismissing(row.adaptation_rate) && isfinite(Float64(row.adaptation_rate)) ? Float64(row.adaptation_rate) : nothing,
    ) for row in eachrow(baselines)]

    stage1_json = [candidate for candidate in stage1 if candidate.rank <= 12]
    advanced = Any[stage1_json...]
    append!(advanced, [_candidate_json(candidate) for candidate in vcat(stage2, stage3, stage4) if candidate.rank <= 12])
    skipped = [_candidate_json(candidate) for candidate in vcat(stage2_skipped, stage3_skipped, stage4_skipped)]
    return (
        schema_version = 2,
        generated_from = "Julia staged ranking/status artifacts",
        calibration_end_day = 14.0,
        observables = ["A2780Naive", "Total A2780cis", "Total population"],
        latent_states_hidden = ["cis sensitive-like", "cis tolerant-like"],
        experimental_library = _untreated_coculture_library(csv_root, stage3_winner),
        dose = (
            drug = "cisplatin",
            unit = "uM",
            presets = [0.67, 1.0, 1.47],
            effect_levels = [0.25, 0.50, 0.75],
            mapping = "piecewise linear concentration-to-IC-effect interpolation",
            measured_range = [0.67, 1.47],
            translational_unit = "relative exposure",
        ),
        winners = (
            stage1 = baseline_json,
            stage2 = _candidate_json(stage2_winner),
            stage3 = _candidate_json(stage3_winner),
            stage4 = _candidate_json(stage4_winner),
        ),
        skipped_nominal_winners = skipped,
        advanced_models = advanced,
        research = [
            (topic = "cisplatin PK/PD", citation = "El-Kareh and Secomb (2003)", url = "https://pubmed.ncbi.nlm.nih.gov/12659689/", role = "advanced intracellular uptake, binding, and delayed cytotoxicity candidate"),
            (topic = "adaptive competition", citation = "Strobl et al. (2021)", url = "https://pubmed.ncbi.nlm.nih.gov/33172930/", role = "competition, turnover, carrying-capacity proximity, resistant fraction, and fitness-cost sensitivity"),
            (topic = "ovarian heterogeneity", citation = "Bicaku et al. (2012)", url = "https://www.nature.com/articles/bjc2012207", role = "ovarian cell-line and cisplatin-response heterogeneity"),
            (topic = "ovarian adaptive therapy", citation = "Hockings et al. (2025)", url = "https://pubmed.ncbi.nlm.nih.gov/40299825/", role = "resistant-fitness and competition structure; not cisplatin PK"),
            (topic = "monitoring policies", citation = "Gallagher et al. (2025 preprint)", url = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11998818/", role = "discrete monitoring and time-varying thresholds"),
            (topic = "credibility", citation = "Braakman et al. (2022)", url = "https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12755/", role = "context of use, verification, validation, and uncertainty"),
        ],
    )
end

function write_a2780_adaptive_config(path, package_root)
    config = load_a2780_adaptive_config(package_root)
    open(path, "w") do io
        JSON3.pretty(io, config)
        println(io)
    end
    return config
end

_hill(dose, emax, ec50, hill) = emax * dose^hill / (max(ec50, eps())^hill + dose^hill)

function _dose_effect_level(dose)
    concentration = max(Float64(dose), 0.0)
    concentration <= 0.67 && return 0.25 * concentration / 0.67
    concentration <= 1.0 && return 0.25 + (concentration - 0.67) * (0.50 - 0.25) / (1.0 - 0.67)
    return 0.50 + (concentration - 1.0) * (0.75 - 0.50) / (1.47 - 1.0)
end

_shifted_fraction(fraction, shift) = inv(1 + exp(-(log(fraction / (1 - fraction)) + shift)))

const A2780_PARAMETER_NAMES = [
    :r_naive, :K_naive, :theta_naive,
    :r_cis, :K_cis, :theta_cis,
    :alpha_nc, :alpha_cn, :death_naive, :death_cis,
    :naive_emax, :naive_ec50, :naive_hill,
    :cis_emax_sensitive, :cis_emax_tolerant, :tolerant_growth_scale,
    :exposure_clearance, :damage_accumulation, :damage_repair,
    :beta_naive, :beta_cis, :naive_lambda, :naive_onset, :cis_onset,
    :naive_growth_adjustment_mode, :naive_growth_adjustment,
    :cis_growth_adjustment_mode, :cis_growth_adjustment,
]

function _growth_adjustment_encoding(row)
    if row.model == "lagged_theta_logistic_growth"
        return (1.0, something(row.lag_time, 0.0))
    elseif row.model == "baranyi_theta_logistic_growth"
        return (2.0, something(row.baranyi_q0, 1.0))
    elseif row.model == "adaptation_theta_logistic_growth"
        return (3.0, something(row.adaptation_rate, 1.0))
    end
    return (0.0, 1.0)
end

function _intrinsic_growth_adjustment(mode, parameter, r, t)
    mode < 0.5 && return 1.0
    mode < 1.5 && return t > max(parameter, 0.0) ? 1.0 : 0.0
    mode < 2.5 && return max(parameter, 1e-8) / (max(parameter, 1e-8) + exp(-max(r, 0.0) * max(t, 0.0)))
    return 1 - exp(-max(parameter, 1e-8) * max(t, 0.0))
end

function _parameter_vector(parameters)
    return [
        parameters.naive.r, parameters.naive.K, parameters.naive.theta,
        parameters.cis.r, parameters.cis.K, parameters.cis.theta,
        parameters.alpha_nc, parameters.alpha_cn, parameters.death_naive, parameters.death_cis,
        parameters.naive_emax, parameters.naive_ec50, parameters.naive_hill,
        parameters.cis_emax_sensitive, parameters.cis_emax_tolerant, parameters.tolerant_growth_scale,
        parameters.exposure_clearance, parameters.damage_accumulation, parameters.damage_repair,
        parameters.beta_naive, parameters.beta_cis, parameters.naive_lambda, parameters.naive_onset, parameters.cis_onset,
        parameters.naive.growth_adjustment_mode, parameters.naive.growth_adjustment,
        parameters.cis.growth_adjustment_mode, parameters.cis.growth_adjustment,
    ]
end

function a2780_model_spec()
    dynamics! = function (du, u, p, t, exposure)
        N = max(u[1], 0.0)
        S = max(u[2], 0.0)
        T = max(u[3], 0.0)
        effective_exposure = max(u[4], 0.0)
        damage = max(u[5], 0.0)
        treatment_clock = u[6]
        C = S + T
        commanded_dose = exposure(t)
        naive_load = max(N + p[7] * C, 0.0)
        cis_load = max(C + p[8] * N, 0.0)
        naive_adjustment = _intrinsic_growth_adjustment(p[25], p[26], p[1], t)
        cis_adjustment = _intrinsic_growth_adjustment(p[27], p[28], p[4], t)
        naive_growth = naive_adjustment * p[1] * N * max(0.0, 1 - (naive_load / p[2])^p[3]) - p[9] * N
        cis_growth = cis_adjustment * p[4] * C * max(0.0, 1 - (cis_load / p[5])^p[6]) - p[10] * C
        sensitive_share = C > 0 ? S / C : 1.0
        tolerant_share = C > 0 ? T / C : 0.0
        pharmacodynamic_active = commanded_dose > 0 || effective_exposure > 0.01
        active_concentration = commanded_dose > 0 ? commanded_dose : effective_exposure
        effect = _dose_effect_level(active_concentration)
        elapsed = max(treatment_clock, 0.0)
        naive_activation = pharmacodynamic_active && elapsed > p[23] ? 1 - exp(-p[22] * (elapsed - p[23])) : 0.0
        cis_activation = pharmacodynamic_active && elapsed > p[24] ? 1.0 : 0.0
        naive_scale = pharmacodynamic_active ? exp(clamp(p[20] * naive_load / p[2], -4.0, 4.0)) : 1.0
        cis_scale = pharmacodynamic_active ? exp(clamp(p[21] * cis_load / p[5], -4.0, 4.0)) : 1.0
        naive_kill = naive_scale * naive_activation * _hill(effect, p[11], p[12], p[13])
        sensitive_kill = cis_scale * cis_activation * _hill(effect, p[14], 0.5, 4.0)
        tolerant_kill = cis_scale * cis_activation * _hill(effect, p[15], 0.5, 4.0)
        du[1] = naive_growth - naive_kill * N
        du[2] = cis_growth * sensitive_share - sensitive_kill * S
        du[3] = (pharmacodynamic_active ? p[16] : 1.0) * cis_growth * tolerant_share - tolerant_kill * T
        du[4] = p[17] * (commanded_dose - effective_exposure)
        du[5] = p[18] * effective_exposure - p[19] * damage
        du[6] = commanded_dose > 0 ? 1.0 : 0.0
        return nothing
    end
    return GrowthParameterEstimation.ModelSpec(
        name = "a2780_data_selected_adaptive",
        ode! = dynamics!,
        param_names = A2780_PARAMETER_NAMES,
        bounds = fill((-Inf, Inf), length(A2780_PARAMETER_NAMES)),
        n_states = 6,
        observable = state -> state[1] + state[2] + state[3],
        base_growth_family = "A2780 data-selected linked model",
        state_names = [:naive, :cis_sensitive_like, :cis_tolerant_like, :effective_exposure, :damage_signal, :treatment_clock],
        metadata = Dict(:context_of_use => "research schedule simulation"),
    )
end

function a2780_segment_simulator(model, state, parameters, t0, t1, dose)
    current_state = length(state) == 5 ? vcat(Float64.(state), 0.0) : Float64.(state)
    if dose > 0 && current_state[4] <= 0.01
        current_state[6] = 0.0
    elseif dose == 0 && current_state[4] <= 0.01
        current_state[6] = 0.0
    end
    result = GrowthParameterEstimation.simulate(
        model,
        [t0, t1],
        _parameter_vector(parameters);
        u0 = current_state,
        exposure = GrowthParameterEstimation.ConstantExposure(Float64(dose)),
        reltol = 1e-8,
        abstol = 1e-8,
    )
    result.success || error("GrowthParameterEstimation segment simulation failed: $(result.reason)")
    current = result.states[:, end]
    return (state = current, effective_exposure = current[4], damage_signal = current[5])
end

function build_a2780_scenario(config, protocol; density = "20k", initial_state = nothing)
    baselines = config.winners.stage1
    naive_row = only(filter(row -> row.cell_line == "A2780Naive" && row.density == density, baselines))
    cis_row = only(filter(row -> row.cell_line == "A2780cis" && row.density == density, baselines))
    competition = config.winners.stage3.parameters
    treatment = config.winners.stage4.parameters
    contrast = get(competition, "log_contrast_density", 0.0)
    density_sign = density == "20k" ? -1.0 : 1.0
    cis_density_scale = exp(density_sign * get(treatment, "cis_log_contrast_density", 0.0))
    tolerant_fraction = _shifted_fraction(treatment["cis_f_tolerant0"], get(treatment, "cis_tolerant_logit_shift", 0.0))
    naive_growth_mode, naive_growth_adjustment = _growth_adjustment_encoding(naive_row)
    cis_growth_mode, cis_growth_adjustment = _growth_adjustment_encoding(cis_row)
    parameters = (
        naive = (
            model = naive_row.model,
            r = naive_row.r,
            K = naive_row.K,
            theta = something(naive_row.theta, something(naive_row.shape_value, 1.0)),
            growth_adjustment_mode = naive_growth_mode,
            growth_adjustment = naive_growth_adjustment,
        ),
        cis = (
            model = cis_row.model,
            r = cis_row.r,
            K = cis_row.K,
            theta = something(cis_row.theta, something(cis_row.shape_value, 1.0)),
            growth_adjustment_mode = cis_growth_mode,
            growth_adjustment = cis_growth_adjustment,
        ),
        alpha_nc = competition["alpha_sr"] * exp(density_sign * contrast),
        alpha_cn = competition["alpha_rs"] * exp(density_sign * contrast),
        death_naive = competition["death_sensitive"],
        death_cis = competition["death_resistant"],
        naive_emax = treatment["naive_emax"] * exp(density_sign * get(treatment, "naive_log_contrast_density", 0.0)),
        naive_ec50 = treatment["naive_ec50_effect"],
        naive_hill = treatment["naive_hill_n"],
        cis_emax_sensitive = treatment["cis_emax_sensitive"] * cis_density_scale,
        cis_emax_tolerant = treatment["cis_emax_tolerant"] * cis_density_scale * exp(clamp(get(treatment, "log_cis_tolerant_context", 0.0), -4.0, 4.0)),
        tolerant_growth_scale = exp(get(treatment, "log_cis_tolerant_growth_context", 0.0)),
        exposure_clearance = log(2) / 0.5,
        damage_accumulation = 1.0,
        damage_repair = log(2) / 1.0,
        beta_naive = get(treatment, "beta_naive", 0.0),
        beta_cis = get(treatment, "beta_cis", 0.0),
        naive_lambda = treatment["naive_lambda"],
        naive_onset = treatment["naive_t_onset"],
        cis_onset = treatment["cis_t_onset"],
        initial_tolerant_fraction = tolerant_fraction,
    )
    state = initial_state === nothing ? [50.0, 50.0 * (1 - tolerant_fraction), 50.0 * tolerant_fraction, 0.0, 0.0, 0.0] : initial_state
    observable = state -> [state[1], state[2] + state[3], state[1] + state[2] + state[3]]
    return SimulationScenario(a2780_model_spec(), parameters, state, protocol;
        observable = observable,
        segment_simulator = a2780_segment_simulator,
    )
end

end
