module Workflow

using CSV
using DataFrames
using Dates
using Random
using Statistics
using TOML
using Optimization
using OptimizationOptimJL

using ..DataLayer
using ..Exposure
using ..Registry
using ..Simulation

export FitCondition, PipelineConfig,
    PipelineStage,
       default_config, save_config, load_config,
    default_stages, default_population_stages, default_population_cellline_stages, summarize_datasets,
    validate_strict_schema, generate_qc_report, save_qc_report,
    save_run_manifest, load_run_manifest,
    bootstrap_stage_uncertainty,
    build_conditions, fit, rank_models, plot_topk, export_results, run_pipeline, run_staged_pipeline

struct FitCondition
    name::String
    time::Vector{Float64}
    count::Vector{Float64}
    error::Vector{Float64}
    u0::Vector{Float64}
    exposure::Exposure.AbstractExposure
    metadata::Dict{Symbol,Any}
end

struct PipelineConfig
    version::String
    model_names::Vector{String}
    n_starts::Int
    top_k::Int
    maxiters::Int
    reltol::Float64
    abstol::Float64
    weighted::Bool
    seed::Int
    output_dir::String
end

struct PipelineStage
    name::String
    description::String
    condition_filter::Function
    condition_cols::Vector{Symbol}
    model_names::Vector{String}
    shared_params::Vector{Symbol}
    fixed_params::Dict{Symbol,Float64}
    inherited_params::Dict{Symbol,Tuple{String,Symbol}}
end

function default_config(; output_dir::String = "results")
    return PipelineConfig(
        "1.0.0",
        Registry.list_models(),
        20,
        5,
        800,
        1e-8,
        1e-8,
        true,
        42,
        output_dir,
    )
end

function _stage_filter(; culture_type=nothing, treated=nothing, population_types::Vector{String}=String[], cell_lines::Vector{String}=String[])
    expected_culture = isnothing(culture_type) ? nothing : lowercase(String(culture_type))
    expected_populations = Set(lowercase.(population_types))
    expected_cell_lines = Set(lowercase.(cell_lines))

    return function (row)
        if expected_culture !== nothing
            if !haskey(parent(row), :culture_type)
                return false
            end
            row_culture = lowercase(string(row[:culture_type]))
            if row_culture != expected_culture
                return false
            end
        end

        if treated !== nothing
            dose_value = haskey(parent(row), :dose) ? Float64(row[:dose]) : 0.0
            if treated && dose_value <= 0.0
                return false
            elseif !treated && dose_value > 0.0
                return false
            end
        end

        if !isempty(expected_populations)
            if !haskey(parent(row), :population_type)
                return false
            end
            population_value = lowercase(string(row[:population_type]))
            if !(population_value in expected_populations)
                return false
            end
        end

        if !isempty(expected_cell_lines)
            if !haskey(parent(row), :cell_line)
                return false
            end
            row_cell_line = lowercase(string(row[:cell_line]))
            if !(row_cell_line in expected_cell_lines)
                return false
            end
        end

        return true
    end
end

function default_stages()
    return [
        PipelineStage(
            "untreated_monoculture",
            "Baseline monoculture fits on untreated conditions.",
            _stage_filter(culture_type="monoculture", treated=false),
            [:cell_line, :population_type, :density, :replicate],
            ["logistic_growth", "gompertz_growth"],
            [:r, :K],
            Dict{Symbol,Float64}(),
            Dict{Symbol,Tuple{String,Symbol}}(),
        ),
        PipelineStage(
            "treated_monoculture",
            "Treatment-response fits that can inherit untreated carrying capacity and growth parameters.",
            _stage_filter(culture_type="monoculture", treated=true),
            [:treatment_amount, :dose, :cell_line, :population_type, :density, :replicate],
            ["theta_logistic_hill_inhibition", "theta_logistic_hill_kill", "pkpd_inhibition", "transit_chain_erlang"],
            [:ic50, :hill],
            Dict{Symbol,Float64}(),
            Dict(
                :r => ("untreated_monoculture", :r),
                :K => ("untreated_monoculture", :K),
            ),
        ),
        PipelineStage(
            "untreated_coculture",
            "Competition model selection on untreated mixed-population conditions.",
            _stage_filter(culture_type="coculture", treated=false),
            [:cell_line, :density, :replicate],
            ["null_coculture", "lotka_volterra_competition"],
            Symbol[],
            Dict{Symbol,Float64}(),
            Dict{Symbol,Tuple{String,Symbol}}(),
        ),
        PipelineStage(
            "treated_coculture",
            "Full treated coculture model ranking with interaction and drug-response terms.",
            _stage_filter(culture_type="coculture", treated=true),
            [:treatment_amount, :dose, :cell_line, :density, :replicate],
            ["lotka_volterra_hill_competition", "sensitive_resistant"],
            [:ic50S, :ic50R, :hill],
            Dict{Symbol,Float64}(),
            Dict{Symbol,Tuple{String,Symbol}}(),
        ),
    ]
end

function default_population_stages(populations::Vector{String} = ["naive", "resistant"])
    stages = PipelineStage[]

    for population in populations
        pop_clean = lowercase(strip(population))
        untreated_name = "untreated_monoculture_" * pop_clean
        treated_name = "treated_monoculture_" * pop_clean

        push!(stages, PipelineStage(
            untreated_name,
            "Untreated monoculture fit for population " * pop_clean * " with global r and K.",
            _stage_filter(culture_type="monoculture", treated=false, population_types=[pop_clean]),
            [:cell_line, :density, :replicate],
            ["logistic_growth", "gompertz_growth"],
            [:r, :K],
            Dict{Symbol,Float64}(),
            Dict{Symbol,Tuple{String,Symbol}}(),
        ))

        push!(stages, PipelineStage(
            treated_name,
            "Treated monoculture fit for population " * pop_clean * "; r and K inherited, IC50 and Hill global within population.",
            _stage_filter(culture_type="monoculture", treated=true, population_types=[pop_clean]),
            [:treatment_amount, :dose, :cell_line, :density, :replicate],
            ["theta_logistic_hill_inhibition", "theta_logistic_hill_kill", "pkpd_inhibition", "transit_chain_erlang"],
            [:ic50, :hill],
            Dict{Symbol,Float64}(),
            Dict(
                :r => (untreated_name, :r),
                :K => (untreated_name, :K),
            ),
        ))
    end

    push!(stages, PipelineStage(
        "untreated_coculture",
        "Untreated coculture interaction model selection with global growth/competition terms.",
        _stage_filter(culture_type="coculture", treated=false),
        [:cell_line, :density, :replicate],
        ["null_coculture", "lotka_volterra_competition"],
        [:rS, :KS, :rR, :KR, :alpha_SR, :alpha_RS],
        Dict{Symbol,Float64}(),
        Dict{Symbol,Tuple{String,Symbol}}(),
    ))

    push!(stages, PipelineStage(
        "treated_coculture",
        "Treated coculture fit where dose/treatment amount varies but IC50 values remain global by subtype.",
        _stage_filter(culture_type="coculture", treated=true),
        [:treatment_amount, :dose, :cell_line, :density, :replicate],
        ["lotka_volterra_hill_competition", "sensitive_resistant"],
        [:ic50S, :ic50R, :hill],
        Dict{Symbol,Float64}(),
        Dict{Symbol,Tuple{String,Symbol}}(),
    ))

    return stages
end

function default_population_cellline_stages(
    df::DataFrame;
    populations::Vector{String} = ["naive", "resistant"],
)
    normalized = DataLayer.normalize_schema(df)
    available_lines = sort(unique(string.(normalized.cell_line)))

    stages = PipelineStage[]

    for population in populations
        pop_clean = lowercase(strip(population))
        untreated_name = "untreated_monoculture_" * pop_clean

        push!(stages, PipelineStage(
            untreated_name,
            "Untreated monoculture fit for population " * pop_clean * " with global r and K.",
            _stage_filter(culture_type="monoculture", treated=false, population_types=[pop_clean]),
            [:cell_line, :density, :replicate],
            ["logistic_growth", "gompertz_growth"],
            [:r, :K],
            Dict{Symbol,Float64}(),
            Dict{Symbol,Tuple{String,Symbol}}(),
        ))

        for cell_line in available_lines
            cell_tag = replace(lowercase(cell_line), r"[^a-z0-9]+" => "_")
            treated_name = "treated_monoculture_" * pop_clean * "_" * cell_tag

            push!(stages, PipelineStage(
                treated_name,
                "Treated monoculture fit for population " * pop_clean * " and cell line " * cell_line * "; IC50 remains fixed for this pair downstream.",
                _stage_filter(culture_type="monoculture", treated=true, population_types=[pop_clean], cell_lines=[cell_line]),
                [:treatment_amount, :dose, :density, :replicate],
                ["theta_logistic_hill_inhibition", "theta_logistic_hill_kill", "pkpd_inhibition", "transit_chain_erlang"],
                [:ic50, :hill],
                Dict{Symbol,Float64}(),
                Dict(
                    :r => (untreated_name, :r),
                    :K => (untreated_name, :K),
                ),
            ))
        end
    end

    for cell_line in available_lines
        cell_tag = replace(lowercase(cell_line), r"[^a-z0-9]+" => "_")
        naive_treated_name = "treated_monoculture_naive_" * cell_tag
        resistant_treated_name = "treated_monoculture_resistant_" * cell_tag

        push!(stages, PipelineStage(
            "untreated_coculture_" * cell_tag,
            "Untreated coculture interaction model selection for cell line " * cell_line * ".",
            _stage_filter(culture_type="coculture", treated=false, cell_lines=[cell_line]),
            [:density, :replicate],
            ["null_coculture", "lotka_volterra_competition"],
            [:rS, :KS, :rR, :KR, :alpha_SR, :alpha_RS],
            Dict{Symbol,Float64}(),
            Dict(
                :rS => ("untreated_monoculture_naive", :r),
                :KS => ("untreated_monoculture_naive", :K),
                :rR => ("untreated_monoculture_resistant", :r),
                :KR => ("untreated_monoculture_resistant", :K),
            ),
        ))

        push!(stages, PipelineStage(
            "treated_coculture_" * cell_tag,
            "Treated coculture for cell line " * cell_line * "; IC50S/IC50R inherited from treated monocultures for the same line.",
            _stage_filter(culture_type="coculture", treated=true, cell_lines=[cell_line]),
            [:treatment_amount, :dose, :density, :replicate],
            ["lotka_volterra_hill_competition", "sensitive_resistant"],
            [:ic50S, :ic50R, :hill],
            Dict{Symbol,Float64}(),
            Dict(
                :rS => ("untreated_monoculture_naive", :r),
                :rR => ("untreated_monoculture_resistant", :r),
                :KS => ("untreated_monoculture_naive", :K),
                :KR => ("untreated_monoculture_resistant", :K),
                :ic50S => (naive_treated_name, :ic50),
                :ic50R => (resistant_treated_name, :ic50),
            ),
        ))
    end

    return stages
end

function validate_strict_schema(
    df::DataFrame;
    required_metadata::Vector{Symbol} = copy(DataLayer.STRICT_REQUIRED_METADATA),
)
    normalized = DataLayer.normalize_schema(df)
    DataLayer.validate_timeseries(normalized)
    DataLayer.validate_required_metadata(normalized; required_metadata=required_metadata)
    return true
end

function _git_commit_hash()
    try
        return readchomp(`git rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _new_failure_log()
    return DataFrame(
        stage=String[],
        model=String[],
        condition=String[],
        reason=String[],
        detail=String[],
        hint=String[],
        timestamp=String[],
    )
end

function _record_failure!(df::DataFrame; stage::AbstractString, model::AbstractString, condition::AbstractString, reason::AbstractString, detail::AbstractString = "", hint::AbstractString = "")
    push!(df, (
        stage=String(stage),
        model=String(model),
        condition=String(condition),
        reason=String(reason),
        detail=String(detail),
        hint=String(hint),
        timestamp=string(now()),
    ))
    return df
end

function _config_to_dict(cfg::PipelineConfig)
    return Dict(
        "version" => cfg.version,
        "model_names" => cfg.model_names,
        "n_starts" => cfg.n_starts,
        "top_k" => cfg.top_k,
        "maxiters" => cfg.maxiters,
        "reltol" => cfg.reltol,
        "abstol" => cfg.abstol,
        "weighted" => cfg.weighted,
        "seed" => cfg.seed,
        "output_dir" => cfg.output_dir,
    )
end

function _serialize_param_bank(parameter_bank::Dict{String,Dict{Symbol,Float64}})
    out = Dict{String,Any}()
    for (stage_name, param_map) in parameter_bank
        out[stage_name] = Dict(String(k) => v for (k, v) in param_map)
    end
    return out
end

function _deserialize_param_bank(raw)
    out = Dict{String,Dict{Symbol,Float64}}()
    for (stage_name, param_map_any) in raw
        param_map = Dict{Symbol,Float64}()
        for (k, v) in param_map_any
            param_map[Symbol(k)] = Float64(v)
        end
        out[String(stage_name)] = param_map
    end
    return out
end

function _serialize_uncertainty_bank(uncertainty_bank::Dict{String,Dict{Symbol,Dict{String,Float64}}})
    out = Dict{String,Any}()
    for (stage_name, param_map) in uncertainty_bank
        out[stage_name] = Dict(String(param) => stats for (param, stats) in param_map)
    end
    return out
end

function _deserialize_uncertainty_bank(raw)
    out = Dict{String,Dict{Symbol,Dict{String,Float64}}}()
    for (stage_name, param_map_any) in raw
        param_map = Dict{Symbol,Dict{String,Float64}}()
        for (k, stats_any) in param_map_any
            stats_map = Dict{String,Float64}()
            for (skey, sval) in stats_any
                stats_map[String(skey)] = Float64(sval)
            end
            param_map[Symbol(k)] = stats_map
        end
        out[String(stage_name)] = param_map
    end
    return out
end

function _stage_result_summary(stage_result)
    return Dict(
        "name" => stage_result.name,
        "description" => stage_result.description,
        "status" => stage_result.status,
        "n_conditions" => stage_result.n_conditions,
        "candidate_models" => collect(stage_result.candidate_models),
        "selected_model" => isnothing(stage_result.selected_model) ? "" : String(stage_result.selected_model),
        "fixed_params" => Dict(String(k) => v for (k, v) in stage_result.fixed_params),
        "inherited_params" => Dict(String(k) => v for (k, v) in stage_result.inherited_params),
        "output_dir" => isnothing(stage_result.output_dir) ? "" : String(stage_result.output_dir),
    )
end

function save_run_manifest(
    path::AbstractString;
    config::PipelineConfig,
    stage_results,
    parameter_bank::Dict{String,Dict{Symbol,Float64}},
    uncertainty_bank::Dict{String,Dict{Symbol,Dict{String,Float64}}} = Dict{String,Dict{Symbol,Dict{String,Float64}}}(),
    failures::DataFrame = _new_failure_log(),
    resume_from_stage = nothing,
    completed::Bool = false,
    halted_stage = nothing,
)
    manifest = Dict(
        "timestamp" => string(now()),
        "git_commit_hash" => _git_commit_hash(),
        "completed" => completed,
        "resume_from_stage" => isnothing(resume_from_stage) ? "" : String(resume_from_stage),
        "halted_stage" => isnothing(halted_stage) ? "" : String(halted_stage),
        "config" => _config_to_dict(config),
        "parameter_bank" => _serialize_param_bank(parameter_bank),
        "uncertainty_bank" => _serialize_uncertainty_bank(uncertainty_bank),
        "stages" => [_stage_result_summary(sr) for sr in stage_results],
        "failure_count" => nrow(failures),
    )

    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, manifest)
    end
    return path
end

function load_run_manifest(path::AbstractString)
    raw = TOML.parsefile(path)
    parameter_bank = haskey(raw, "parameter_bank") ? _deserialize_param_bank(raw["parameter_bank"]) : Dict{String,Dict{Symbol,Float64}}()
    uncertainty_bank = haskey(raw, "uncertainty_bank") ? _deserialize_uncertainty_bank(raw["uncertainty_bank"]) : Dict{String,Dict{Symbol,Dict{String,Float64}}}()
    return (
        raw=raw,
        parameter_bank=parameter_bank,
        uncertainty_bank=uncertainty_bank,
        completed_stages=haskey(raw, "stages") ? [String(stage["name"]) for stage in raw["stages"] if get(stage, "status", "") == "completed"] : String[],
    )
end

function generate_qc_report(
    df::DataFrame;
    condition_cols::Vector{Symbol} = [:cell_line, :population_type, :culture_type, :treatment_amount, :density, :replicate],
)
    normalized = DataLayer.normalize_schema(df)
    missingness = DataFrame(
        column=String[],
        n_missing=Int[],
        frac_missing=Float64[],
    )

    for col in names(normalized)
        n_missing = count(ismissing, normalized[!, col])
        push!(missingness, (String(col), n_missing, n_missing / max(nrow(normalized), 1)))
    end

    cols = _resolve_condition_cols(normalized, condition_cols)
    grouped = groupby(normalized, cols)

    condition_summary = DataFrame(
        condition=String[],
        n_rows=Int[],
        n_replicates=Int[],
        min_time=Float64[],
        max_time=Float64[],
        monotone_time=Bool[],
    )
    issues = DataFrame(severity=String[], issue=String[], condition=String[], detail=String[])
    outliers = DataFrame(condition=String[], row_index=Int[], time=Float64[], count=Float64[], z_score=Float64[])

    for g in grouped
        condition = join(["$(c)=$(g[1, c])" for c in cols], " | ")
        time_values = Float64.(g.time)
        monotone = all(diff(time_values) .>= 0)
        push!(condition_summary, (
            condition,
            nrow(g),
            length(unique(g.replicate)),
            minimum(time_values),
            maximum(time_values),
            monotone,
        ))

        if !monotone
            push!(issues, ("error", "non_monotone_time", condition, "Time is not monotone within condition"))
        end

        counts = Float64.(g.count)
        if length(counts) >= 3
            count_std = std(counts)
            if isfinite(count_std) && count_std > 0
                zscores = abs.((counts .- mean(counts)) ./ count_std)
                for (idx, z) in enumerate(zscores)
                    if z > 3.0
                        push!(outliers, (condition, idx, time_values[idx], counts[idx], z))
                    end
                end
            end
        end
    end

    return (
        missingness=missingness,
        condition_summary=condition_summary,
        outliers=outliers,
        issues=issues,
    )
end

function save_qc_report(qc_report; output_dir::String)
    mkpath(output_dir)
    missingness_path = joinpath(output_dir, "qc_missingness.csv")
    summary_path = joinpath(output_dir, "qc_condition_summary.csv")
    outliers_path = joinpath(output_dir, "qc_outliers.csv")
    issues_path = joinpath(output_dir, "qc_issues.csv")

    CSV.write(missingness_path, qc_report.missingness)
    CSV.write(summary_path, qc_report.condition_summary)
    CSV.write(outliers_path, qc_report.outliers)
    CSV.write(issues_path, qc_report.issues)

    return (
        missingness=missingness_path,
        condition_summary=summary_path,
        outliers=outliers_path,
        issues=issues_path,
    )
end

function _bootstrap_resample(df::DataFrame, condition_cols::Vector{Symbol}, rng::AbstractRNG)
    grouped = groupby(df, condition_cols)
    parts = DataFrame[]
    for g in grouped
        idx = rand(rng, 1:nrow(g), nrow(g))
        sample = g[idx, :]
        sort!(sample, :time)
        push!(parts, sample)
    end
    return reduce(vcat, parts)
end

function bootstrap_stage_uncertainty(
    spec::Registry.ModelSpec,
    stage_df::DataFrame;
    condition_cols::Vector{Symbol},
    shared_params::Vector{Symbol} = copy(spec.param_names),
    fixed_params::Dict{Symbol,Float64} = Dict{Symbol,Float64}(),
    tie_constraints::Dict{Symbol,Symbol} = Dict{Symbol,Symbol}(),
    n_bootstrap::Int = 20,
    n_starts::Int = 5,
    maxiters::Int = 300,
    weighted::Bool = true,
    reltol::Float64 = 1e-8,
    abstol::Float64 = 1e-8,
    seed::Int = 2026,
)
    normalized = DataLayer.normalize_schema(stage_df)
    conditions = build_conditions(normalized; condition_cols=condition_cols)
    isempty(conditions) && return Dict{Symbol,Dict{String,Float64}}()

    rng = MersenneTwister(seed)
    samples = Dict{Symbol,Vector{Float64}}()

    for b in 1:n_bootstrap
        boot_df = _bootstrap_resample(normalized, condition_cols, rng)
        boot_conditions = build_conditions(boot_df; condition_cols=condition_cols)
        try
            fit_result = fit(
                spec,
                boot_conditions;
                shared_params=shared_params,
                fixed_params=fixed_params,
                tie_constraints=tie_constraints,
                n_starts=n_starts,
                maxiters=maxiters,
                weighted=weighted,
                reltol=reltol,
                abstol=abstol,
                seed=seed + b,
                top_k=1,
            )
            summary = _mean_parameter_summary(spec, fit_result)
            for (param, value) in summary
                push!(get!(samples, param, Float64[]), value)
            end
        catch
        end
    end

    stats = Dict{Symbol,Dict{String,Float64}}()
    for (param, values) in samples
        if isempty(values)
            continue
        end
        sorted_values = sort(values)
        stats[param] = Dict(
            "mean" => mean(values),
            "std" => (length(values) > 1 ? std(values) : 0.0),
            "ci_lower" => sorted_values[max(1, ceil(Int, 0.025 * length(sorted_values)))],
            "ci_upper" => sorted_values[min(length(sorted_values), ceil(Int, 0.975 * length(sorted_values)))],
            "n_success" => length(values),
        )
    end

    return stats
end

function summarize_datasets(
    df::DataFrame;
    group_cols::Vector{Symbol} = [:cell_line, :population_type, :culture_type, :treatment_amount, :dose, :density, :replicate],
    ic50_col::Symbol = :ic50_reference,
)
    DataLayer.validate_timeseries(df)

    cols = _resolve_condition_cols(df, group_cols)
    grouped = groupby(df, cols)
    rows = NamedTuple[]

    for g in grouped
        tmin = minimum(g.time)
        tmax = maximum(g.time)
        duration_days = tmax - tmin

        metadata = Dict{Symbol,Any}(c => g[1, c] for c in cols)
        dose_values = haskey(metadata, :treatment_amount) ? Float64(metadata[:treatment_amount]) : (haskey(metadata, :dose) ? Float64(metadata[:dose]) : 0.0)

        ic50_value = if ic50_col in Symbol.(names(g))
            vals = [v for v in g[!, ic50_col] if !ismissing(v)]
            isempty(vals) ? missing : mean(Float64.(vals))
        else
            missing
        end

        push!(rows, (
            n_timepoints = nrow(g),
            min_time = tmin,
            max_time = tmax,
            duration_days = duration_days,
            treatment_amount = dose_values,
            ic50_reference = ic50_value,
            condition = join(["$(c)=$(metadata[c])" for c in cols], " | "),
        ))
    end

    out = DataFrame(rows)
    if !isempty(out)
        sort!(out, [:duration_days, :treatment_amount], rev = [true, false])
    end
    return out
end

function save_config(path::AbstractString, cfg::PipelineConfig)
    dict = Dict(
        "version" => cfg.version,
        "model_names" => cfg.model_names,
        "n_starts" => cfg.n_starts,
        "top_k" => cfg.top_k,
        "maxiters" => cfg.maxiters,
        "reltol" => cfg.reltol,
        "abstol" => cfg.abstol,
        "weighted" => cfg.weighted,
        "seed" => cfg.seed,
        "output_dir" => cfg.output_dir,
    )
    open(path, "w") do io
        TOML.print(io, dict)
    end
    return path
end

function load_config(path::AbstractString)
    cfg = TOML.parsefile(path)
    return PipelineConfig(
        string(get(cfg, "version", "1.0.0")),
        String.(get(cfg, "model_names", Registry.list_models())),
        Int(get(cfg, "n_starts", 20)),
        Int(get(cfg, "top_k", 5)),
        Int(get(cfg, "maxiters", 800)),
        Float64(get(cfg, "reltol", 1e-8)),
        Float64(get(cfg, "abstol", 1e-8)),
        Bool(get(cfg, "weighted", true)),
        Int(get(cfg, "seed", 42)),
        string(get(cfg, "output_dir", "results")),
    )
end

function build_conditions(df::DataFrame; condition_cols::Vector{Symbol} = [:dose, :cell_line, :density, :replicate])
    normalized = DataLayer.normalize_schema(df)
    DataLayer.validate_timeseries(normalized)
    resolved_cols = _resolve_condition_cols(normalized, condition_cols)
    grouped = groupby(normalized, resolved_cols)

    conditions = FitCondition[]
    for g in grouped
        g_sorted = sort(g, :time)
        metadata = Dict{Symbol,Any}(c => g_sorted[1, c] for c in resolved_cols)
        cname = join(["$(c)=$(metadata[c])" for c in resolved_cols], " | ")

        dose_value = haskey(metadata, :dose) ? Float64(metadata[:dose]) : 0.0
        exposure = Exposure.ConstantExposure(dose_value)

        push!(conditions, FitCondition(
            cname,
            Float64.(g_sorted.time),
            Float64.(g_sorted.count),
            Float64.(g_sorted.error),
            [max(1e-12, Float64(g_sorted.count[1]))],
            exposure,
            metadata,
        ))
    end
    return conditions
end

function _resolve_condition_cols(df::DataFrame, requested::Vector{Symbol})
    available = Set(Symbol.(names(df)))
    cols = [col for col in requested if col in available]
    if isempty(cols)
        cols = [col for col in [:dose, :cell_line, :density, :replicate] if col in available]
    end
    isempty(cols) && error("No valid condition columns available for grouping")
    return cols
end

function _successful_model_names(rank_result)
    return [name for name in rank_result.ranking.model if haskey(rank_result.fits, name) && isfinite(rank_result.ranking[findfirst(==(name), rank_result.ranking.model), :bic])]
end

function _mean_parameter_summary(spec::Registry.ModelSpec, fit_result)
    values = Dict{Symbol,Vector{Float64}}()
    for per_condition in fit_result.per_condition
        per_condition.success || continue
        for (idx, name) in enumerate(spec.param_names)
            push!(get!(values, name, Float64[]), per_condition.params[idx])
        end
    end

    summary = Dict{Symbol,Float64}()
    for name in spec.param_names
        if haskey(values, name) && !isempty(values[name])
            summary[name] = mean(values[name])
        end
    end

    return summary
end

function _resolve_inherited_fixed_params(stage::PipelineStage, parameter_bank::Dict{String,Dict{Symbol,Float64}})
    merged = copy(stage.fixed_params)
    for (dest_param, (source_stage, source_param)) in stage.inherited_params
        if haskey(parameter_bank, source_stage) && haskey(parameter_bank[source_stage], source_param)
            merged[dest_param] = parameter_bank[source_stage][source_param]
        end
    end
    return merged
end

function _serialize_inherited_map(map::Dict{Symbol,Tuple{String,Symbol}})
    return Dict(String(k) => Dict("source_stage" => v[1], "source_param" => String(v[2])) for (k, v) in map)
end

function _select_stage_model(rank_result, stage::PipelineStage, selection_mode::Symbol, manual_choices::Dict{String,String})
    candidates = _successful_model_names(rank_result)
    isempty(candidates) && return nothing

    if selection_mode == :manual
        return get(manual_choices, stage.name, nothing)
    elseif selection_mode == :best_bic
        return candidates[1]
    else
        throw(ArgumentError("selection_mode must be :best_bic or :manual"))
    end
end

function run_staged_pipeline(
    data_input;
    stages::Vector{PipelineStage} = default_stages(),
    config::PipelineConfig = default_config(),
    selection_mode::Symbol = :best_bic,
    manual_choices::Dict{String,String} = Dict{String,String}(),
    export_stage_results::Bool = true,
    strict_schema::Bool = false,
    required_metadata::Vector{Symbol} = copy(DataLayer.STRICT_REQUIRED_METADATA),
    qc_before_fit::Bool = true,
    n_bootstrap::Int = 0,
    bootstrap_seed::Int = 2026,
    resume_from_stage = nothing,
    resume_manifest_path = nothing,
)
    df = if data_input isa DataFrame
        DataLayer.normalize_schema(data_input)
    elseif data_input isa AbstractString
        DataLayer.load_timeseries(data_input)
    else
        error("data_input must be a DataFrame or file path")
    end

    DataLayer.validate_timeseries(df)

    if strict_schema
        validate_strict_schema(df; required_metadata=required_metadata)
    end

    qc_report = qc_before_fit ? generate_qc_report(df) : nothing
    failure_log = _new_failure_log()
    qc_paths = nothing

    if qc_before_fit && export_stage_results
        qc_paths = save_qc_report(qc_report; output_dir=joinpath(config.output_dir, "diagnostics"))
    end

    parameter_bank = Dict{String,Dict{Symbol,Float64}}()
    uncertainty_bank = Dict{String,Dict{Symbol,Dict{String,Float64}}}()
    stage_results = NamedTuple[]
    halted_stage = nothing

    if !isnothing(resume_manifest_path)
        resume_state = load_run_manifest(resume_manifest_path)
        parameter_bank = merge(parameter_bank, resume_state.parameter_bank)
        uncertainty_bank = merge(uncertainty_bank, resume_state.uncertainty_bank)
        if isnothing(resume_from_stage) && !isempty(resume_state.completed_stages)
            last_completed = last(resume_state.completed_stages)
            idx = findfirst(s -> s.name == last_completed, stages)
            if !isnothing(idx) && idx < length(stages)
                resume_from_stage = stages[idx + 1].name
            end
        end
    end

    should_run = isnothing(resume_from_stage)

    for stage in stages
        if !should_run
            if stage.name == resume_from_stage
                should_run = true
            else
                push!(stage_results, (
                    name=stage.name,
                    description=stage.description,
                    status="resumed_prior",
                    n_conditions=0,
                    candidate_models=String[],
                    selected_model=nothing,
                    fixed_params=copy(stage.fixed_params),
                    inherited_params=_serialize_inherited_map(stage.inherited_params),
                    output_dir=nothing,
                    uncertainty=Dict{Symbol,Dict{String,Float64}}(),
                    result=nothing,
                ))
                continue
            end
        end

        mask = [stage.condition_filter(row) for row in eachrow(df)]
        stage_df = df[mask, :]

        if nrow(stage_df) == 0
            _record_failure!(failure_log; stage=stage.name, model="none", condition="all", reason="missing_stage_data", detail="No rows matched stage filter", hint="Check culture_type, population_type, cell_line, treatment_amount, and density metadata.")
            push!(stage_results, (
                name=stage.name,
                description=stage.description,
                status="skipped",
                n_conditions=0,
                candidate_models=String[],
                selected_model=nothing,
                fixed_params=copy(stage.fixed_params),
                inherited_params=_serialize_inherited_map(stage.inherited_params),
                output_dir=nothing,
                uncertainty=Dict{Symbol,Dict{String,Float64}}(),
                result=nothing,
            ))
            continue
        end

        condition_cols = _resolve_condition_cols(stage_df, stage.condition_cols)
        conditions = build_conditions(stage_df; condition_cols=condition_cols)
        fixed_params = _resolve_inherited_fixed_params(stage, parameter_bank)

        rank_result = rank_models(
            stage.model_names,
            conditions;
            top_k=config.top_k,
            n_starts=config.n_starts,
            maxiters=config.maxiters,
            weighted=config.weighted,
            reltol=config.reltol,
            abstol=config.abstol,
            seed=config.seed,
            shared_params=stage.shared_params,
            fixed_params=fixed_params,
        )

        candidate_models = _successful_model_names(rank_result)
        candidate_models = candidate_models[1:min(3, length(candidate_models))]
        selected_model = _select_stage_model(rank_result, stage, selection_mode, manual_choices)

        stage_output_dir = export_stage_results ? joinpath(config.output_dir, stage.name) : nothing
        if export_stage_results
            plot_topk(rank_result; conditions=conditions, top_k=min(config.top_k, max(length(candidate_models), 1)), output_dir=joinpath(stage_output_dir, "figures"))
            try
                export_results(rank_result; output_dir=stage_output_dir)
            catch err
                _record_failure!(failure_log; stage=stage.name, model="export", condition="all", reason="stage_export_failure", detail=sprint(showerror, err), hint="Inspect stage output directory permissions and rank_result contents.")
            end
        end

        append!(failure_log, rank_result.failures)

        if isnothing(selected_model)
            push!(stage_results, (
                name=stage.name,
                description=stage.description,
                status="awaiting_selection",
                n_conditions=length(conditions),
                candidate_models=candidate_models,
                selected_model=nothing,
                fixed_params=fixed_params,
                inherited_params=_serialize_inherited_map(stage.inherited_params),
                output_dir=stage_output_dir,
                uncertainty=Dict{Symbol,Dict{String,Float64}}(),
                result=rank_result,
            ))
            halted_stage = stage.name
            if export_stage_results
                save_run_manifest(
                    joinpath(config.output_dir, "run_manifest.toml");
                    config=config,
                    stage_results=stage_results,
                    parameter_bank=parameter_bank,
                    uncertainty_bank=uncertainty_bank,
                    failures=failure_log,
                    resume_from_stage=resume_from_stage,
                    completed=false,
                    halted_stage=halted_stage,
                )
            end
            break
        end

        if !(selected_model in candidate_models)
            error("Selected model $(selected_model) for stage $(stage.name) is not among successful candidates")
        end

        selected_fit = rank_result.fits[selected_model]
        selected_spec = Registry.get_model(selected_model)
        parameter_bank[stage.name] = _mean_parameter_summary(selected_spec, selected_fit)
        uncertainty = n_bootstrap > 0 ? bootstrap_stage_uncertainty(
            selected_spec,
            stage_df;
            condition_cols=condition_cols,
            shared_params=stage.shared_params,
            fixed_params=fixed_params,
            n_bootstrap=n_bootstrap,
            n_starts=min(config.n_starts, 5),
            maxiters=min(config.maxiters, 300),
            weighted=config.weighted,
            reltol=config.reltol,
            abstol=config.abstol,
            seed=bootstrap_seed + length(stage_results),
        ) : Dict{Symbol,Dict{String,Float64}}()
        uncertainty_bank[stage.name] = uncertainty

        push!(stage_results, (
            name=stage.name,
            description=stage.description,
            status="completed",
            n_conditions=length(conditions),
            candidate_models=candidate_models,
            selected_model=selected_model,
            fixed_params=fixed_params,
            inherited_params=_serialize_inherited_map(stage.inherited_params),
            output_dir=stage_output_dir,
            uncertainty=uncertainty,
            result=rank_result,
        ))

        if export_stage_results
            save_run_manifest(
                joinpath(config.output_dir, "run_manifest.toml");
                config=config,
                stage_results=stage_results,
                parameter_bank=parameter_bank,
                uncertainty_bank=uncertainty_bank,
                failures=failure_log,
                resume_from_stage=resume_from_stage,
                completed=false,
                halted_stage=nothing,
            )
        end
    end

    if export_stage_results
        mkpath(joinpath(config.output_dir, "diagnostics"))
        CSV.write(joinpath(config.output_dir, "diagnostics", "staged_failure_report.csv"), failure_log)
        save_run_manifest(
            joinpath(config.output_dir, "run_manifest.toml");
            config=config,
            stage_results=stage_results,
            parameter_bank=parameter_bank,
            uncertainty_bank=uncertainty_bank,
            failures=failure_log,
            resume_from_stage=resume_from_stage,
            completed=isnothing(halted_stage),
            halted_stage=halted_stage,
        )
    end

    return (
        config=config,
        stages=stage_results,
        parameter_bank=parameter_bank,
        uncertainty_bank=uncertainty_bank,
        halted_stage=halted_stage,
        completed=isnothing(halted_stage),
        qc_report=qc_report,
        qc_paths=qc_paths,
        failures=failure_log,
        manifest_path=export_stage_results ? joinpath(config.output_dir, "run_manifest.toml") : nothing,
    )
end

function _canonical_param(param::Symbol, ties::Dict{Symbol,Symbol})
    seen = Set{Symbol}()
    cur = param
    while haskey(ties, cur) && !(cur in seen)
        push!(seen, cur)
        cur = ties[cur]
    end
    return cur
end

function _build_layout(
    spec::Registry.ModelSpec,
    n_conditions::Int;
    shared_params::Vector{Symbol},
    fixed_params::Dict{Symbol,Float64},
    tie_constraints::Dict{Symbol,Symbol},
)
    pnames = spec.param_names

    shared_set = Set(shared_params)
    condition_specific = Set(setdiff(pnames, shared_params))

    bounds_map = Dict(name => spec.bounds[i] for (i, name) in enumerate(pnames))

    index_map = Dict{Tuple{Symbol,Int},Int}()
    theta_lb = Float64[]
    theta_ub = Float64[]
    theta_name = Tuple{Symbol,Int}[]

    next_idx = 1

    for name in pnames
        if haskey(fixed_params, name)
            continue
        end
        if name in shared_set
            canon = _canonical_param(name, tie_constraints)
            key = (canon, 0)
            if !haskey(index_map, key)
                index_map[key] = next_idx
                b = bounds_map[canon]
                push!(theta_lb, b[1]); push!(theta_ub, b[2])
                push!(theta_name, key)
                next_idx += 1
            end
        end
    end

    for cond_idx in 1:n_conditions
        for name in pnames
            if haskey(fixed_params, name)
                continue
            end
            if name in condition_specific
                canon = _canonical_param(name, tie_constraints)
                key = (canon, cond_idx)
                if !haskey(index_map, key)
                    index_map[key] = next_idx
                    b = bounds_map[canon]
                    push!(theta_lb, b[1]); push!(theta_ub, b[2])
                    push!(theta_name, key)
                    next_idx += 1
                end
            end
        end
    end

    return index_map, theta_lb, theta_ub, shared_set, condition_specific
end

function _expand_params(
    theta::AbstractVector{T},
    spec::Registry.ModelSpec,
    cond_idx::Int,
    index_map::Dict{Tuple{Symbol,Int},Int},
    shared_set::Set{Symbol},
    fixed_params::Dict{Symbol,Float64},
    tie_constraints::Dict{Symbol,Symbol},
) where {T}
    p = Vector{T}(undef, length(spec.param_names))
    for (i, name) in enumerate(spec.param_names)
        if haskey(fixed_params, name)
            p[i] = T(fixed_params[name])
            continue
        end
        canon = _canonical_param(name, tie_constraints)
        if name in shared_set
            p[i] = theta[index_map[(canon, 0)]]
        else
            p[i] = theta[index_map[(canon, cond_idx)]]
        end
    end
    return p
end

function _initial_theta(index_map, lb::Vector{Float64}, ub::Vector{Float64}, rng::AbstractRNG)
    theta = similar(lb)
    for i in eachindex(lb)
        width = max(ub[i] - lb[i], 1e-12)
        theta[i] = lb[i] + rand(rng) * width
    end
    return theta
end

function fit(
    spec::Registry.ModelSpec,
    conditions::Vector{FitCondition};
    shared_params::Vector{Symbol} = copy(spec.param_names),
    fixed_params::Dict{Symbol,Float64} = Dict{Symbol,Float64}(),
    tie_constraints::Dict{Symbol,Symbol} = Dict{Symbol,Symbol}(),
    n_starts::Int = 20,
    maxiters::Int = 800,
    weighted::Bool = true,
    reltol::Float64 = 1e-8,
    abstol::Float64 = 1e-8,
    seed::Int = 42,
    top_k::Int = 5,
)
    rng = MersenneTwister(seed)
    n_conditions = length(conditions)
    n_obs = sum(length(c.time) for c in conditions)

    index_map, lb, ub, shared_set, condition_specific = _build_layout(
        spec,
        n_conditions;
        shared_params=shared_params,
        fixed_params=fixed_params,
        tie_constraints=tie_constraints,
    )

    if isempty(lb)
        error("No free parameters remain after fixed/tie constraints")
    end

    failure_log = _new_failure_log()

    function objective(theta, _)
        sse = 0.0
        penalty = 0.0

        for (ci, cond) in enumerate(conditions)
            p = _expand_params(theta, spec, ci, index_map, shared_set, fixed_params, tie_constraints)
            sim = Simulation.simulate(
                spec,
                cond.time,
                p;
                u0=cond.u0,
                exposure=cond.exposure,
                reltol=reltol,
                abstol=abstol,
            )

            if !sim.success
                penalty += 1e12
                continue
            end

            resid = cond.count .- sim.observed
            if weighted
                w = 1.0 ./ max.(cond.error .^ 2, 1e-12)
                sse += sum(w .* (resid .^ 2))
            else
                sse += sum(resid .^ 2)
            end
        end

        return sse + penalty
    end

    loss = OptimizationFunction((x, p) -> objective(x, p), Optimization.AutoForwardDiff())

    rows = NamedTuple[]
    best = nothing

    for start_id in 1:n_starts
        theta0 = if start_id == 1
            (lb .+ ub) ./ 2
        else
            _initial_theta(index_map, lb, ub, rng)
        end

        try
            prob = OptimizationProblem(loss, theta0; lb=lb, ub=ub)
            result = solve(prob, OptimizationOptimJL.Fminbox(OptimizationOptimJL.BFGS()); maxiters=maxiters)

            theta_hat = collect(result.u)
            obj = objective(theta_hat, nothing)

            k = length(theta_hat)
            aic = n_obs * log(max(obj, 1e-12) / n_obs) + 2k
            bic = n_obs * log(max(obj, 1e-12) / n_obs) + k * log(n_obs)

            row = (
                start_id=start_id,
                objective=obj,
                converged=result.retcode == ReturnCode.Success,
                retcode=string(result.retcode),
                aic=aic,
                bic=bic,
                params=theta_hat,
            )
            push!(rows, row)

            if isnothing(best) || obj < best.objective
                best = row
            end
        catch err
            _record_failure!(failure_log; stage="optimization", model=spec.name, condition="all", reason=Simulation.classify_failure(err), detail=sprint(showerror, err), hint="Check initial bounds, parameter identifiability, and solver stability.")
        end
    end

    if isempty(rows)
        error("All optimization starts failed for model $(spec.name)")
    end

    sorted_rows = sort(rows, by=r -> r.objective)
    keep = sorted_rows[1:min(top_k, length(sorted_rows))]

    per_condition = NamedTuple[]
    for (ci, cond) in enumerate(conditions)
        pbest = _expand_params(best.params, spec, ci, index_map, shared_set, fixed_params, tie_constraints)
        sim = Simulation.simulate(
            spec,
            cond.time,
            pbest;
            u0=cond.u0,
            exposure=cond.exposure,
            reltol=reltol,
            abstol=abstol,
        )
        if !sim.success
            _record_failure!(failure_log; stage="solver", model=spec.name, condition=cond.name, reason=sim.reason, detail="Simulation failed during best-fit evaluation", hint="Inspect condition-specific time range, bounds, and observation scaling.")
        end

        push!(per_condition, (
            condition=cond.name,
            params=pbest,
            observed=sim.observed,
            residuals=sim.success ? (cond.count .- sim.observed) : fill(NaN, length(cond.time)),
            success=sim.success,
            reason=sim.reason,
        ))
    end

    return (
        model=spec.name,
        best=best,
        top_fits=keep,
        n_obs=n_obs,
        n_params=length(best.params),
        per_condition=per_condition,
        failures=failure_log,
    )
end

function rank_models(
    model_names::Vector{String},
    conditions::Vector{FitCondition};
    top_k::Int = 5,
    kwargs...
)
    rows = NamedTuple[]
    fit_map = Dict{String,Any}()
    all_failures = _new_failure_log()

    for name in model_names
        spec = Registry.get_model(name)
        try
            fres = fit(spec, conditions; top_k=top_k, kwargs...)
            fit_map[name] = fres
            push!(rows, (
                model=name,
                sse=fres.best.objective,
                weighted_sse=fres.best.objective,
                aic=fres.best.aic,
                bic=fres.best.bic,
                n_params=fres.n_params,
            ))
            append!(all_failures, fres.failures)
        catch err
            push!(rows, (
                model=name,
                sse=Inf,
                weighted_sse=Inf,
                aic=Inf,
                bic=Inf,
                n_params=0,
            ))
            _record_failure!(all_failures; stage="model_fit", model=name, condition="all", reason=Simulation.classify_failure(err), detail=sprint(showerror, err), hint="Check model-to-data mismatch, strict schema metadata, and inherited parameter availability.")
        end
    end

    ranking = DataFrame(rows)
    sort!(ranking, :bic)
    best_bic = minimum(ranking.bic)
    ranking.delta_bic = ranking.bic .- best_bic

    return (ranking=ranking, fits=fit_map, failures=all_failures)
end

function _load_plots_or_nothing()
    return isdefined(Main, :Plots) ? getfield(Main, :Plots) : nothing
end

function plot_topk(
    rank_result;
    conditions::Vector{FitCondition},
    top_k::Int = 5,
    output_dir::String = "results/figures",
)
    mkpath(output_dir)
    ranking = rank_result.ranking
    ranked_models = ranking.model[1:min(top_k, nrow(ranking))]
    top_models = [m for m in ranked_models if haskey(rank_result.fits, m)]

    plot_module = _load_plots_or_nothing()
    has_plots = !isnothing(plot_module)
    generated = String[]

    if isempty(top_models)
        return generated
    end

    for cond in conditions
        overlay = DataFrame(time=cond.time, observed=cond.count)
        for m in top_models
            fit_info = rank_result.fits[m]
            cond_hit = findfirst(pc -> pc.condition == cond.name, fit_info.per_condition)
            if !isnothing(cond_hit)
                overlay[!, Symbol("pred_" * m)] = fit_info.per_condition[cond_hit].observed
            end
        end

        csv_path = joinpath(output_dir, replace(cond.name, '|' => '_', ' ' => '_') * "_overlay.csv")
        CSV.write(csv_path, overlay)
        push!(generated, csv_path)

        if has_plots
            p = Base.invokelatest(
                plot_module.plot,
                cond.time,
                cond.count;
                seriestype=:scatter,
                label="Observed",
                title=cond.name,
                xlabel="Time",
                ylabel="Count",
                legend=:topleft,
                linewidth=2,
                markersize=4,
                size=(1200, 700),
                foreground_color_legend=nothing,
            )
            for m in top_models
                col = Symbol("pred_" * m)
                if col in names(overlay)
                    Base.invokelatest(plot_module.plot!, p, overlay.time, overlay[!, col]; label=m, linewidth=2.5)
                end
            end
            png_path = replace(csv_path, ".csv" => ".png")
            Base.invokelatest(plot_module.savefig, p, png_path)
            push!(generated, png_path)
        end
    end

    bic_rows = [(model=m, bic=ranking[findfirst(==(m), ranking.model), :bic]) for m in top_models]
    bic_df = DataFrame(bic_rows)
    bic_csv = joinpath(output_dir, "top_models_bic.csv")
    CSV.write(bic_csv, bic_df)
    push!(generated, bic_csv)

    if has_plots
        p_bic = Base.invokelatest(
            plot_module.bar,
            bic_df.model,
            bic_df.bic;
            legend=false,
            xlabel="Model",
            ylabel="BIC",
            title="Top Model BIC Comparison",
            size=(1200, 700),
            linecolor=:black,
            linewidth=0.8,
            bar_width=0.6,
        )
        bic_png = joinpath(output_dir, "top_models_bic.png")
        Base.invokelatest(plot_module.savefig, p_bic, bic_png)
        push!(generated, bic_png)
    end

    return generated
end

function export_results(
    rank_result;
    output_dir::String = "results",
)
    tables_dir = joinpath(output_dir, "tables")
    params_dir = joinpath(output_dir, "params")
    diag_dir = joinpath(output_dir, "diagnostics")
    fig_dir = joinpath(output_dir, "figures")

    for d in (tables_dir, params_dir, diag_dir, fig_dir)
        mkpath(d)
    end

    ranking_path = joinpath(tables_dir, "model_ranking.csv")
    CSV.write(ranking_path, rank_result.ranking)

    ranked_models = collect(rank_result.ranking.model)
    successful_models = [m for m in ranked_models if haskey(rank_result.fits, m)]
    isempty(successful_models) && error("No successful model fits available to export.")

    best_model = successful_models[1]
    best_fit = rank_result.fits[best_model]
    best_bic = rank_result.ranking[findfirst(==(best_model), rank_result.ranking.model), :bic]

    starts_df = DataFrame(best_fit.top_fits)
    starts_path = joinpath(tables_dir, "top_fit_starts.csv")
    CSV.write(starts_path, starts_df)

    p_df = DataFrame(param_index=collect(1:length(best_fit.best.params)), value=best_fit.best.params)
    params_path = joinpath(params_dir, "best_params.csv")
    CSV.write(params_path, p_df)

    fail_path = joinpath(diag_dir, "failure_report.csv")
    CSV.write(fail_path, rank_result.failures)

    summary = DataFrame(
        key=["best_model", "best_bic", "generated_at"],
        value=[best_model, string(best_bic), string(now())],
    )
    summary_path = joinpath(output_dir, "best_model_summary.csv")
    CSV.write(summary_path, summary)

    return (
        ranking=ranking_path,
        starts=starts_path,
        params=params_path,
        failures=fail_path,
        summary=summary_path,
        figures=fig_dir,
    )
end

function run_pipeline(
    data_input;
    config::PipelineConfig = default_config(),
    include_models::Vector{String} = String[],
    exclude_models::Vector{String} = String[],
    strict_schema::Bool = false,
    required_metadata::Vector{Symbol} = copy(DataLayer.STRICT_REQUIRED_METADATA),
    qc_before_fit::Bool = true,
)
    df = if data_input isa DataFrame
        DataLayer.normalize_schema(data_input)
    elseif data_input isa AbstractString
        DataLayer.load_timeseries(data_input)
    else
        error("data_input must be a DataFrame or file path")
    end

    DataLayer.validate_timeseries(df)
    if strict_schema
        validate_strict_schema(df; required_metadata=required_metadata)
    end

    qc_report = qc_before_fit ? generate_qc_report(df) : nothing
    qc_paths = qc_before_fit ? save_qc_report(qc_report; output_dir=joinpath(config.output_dir, "diagnostics")) : nothing
    conditions = build_conditions(df)

    model_names = isempty(include_models) ? copy(config.model_names) : include_models
    if !isempty(exclude_models)
        model_names = setdiff(model_names, exclude_models)
    end

    rank_result = rank_models(
        model_names,
        conditions;
        top_k=config.top_k,
        n_starts=config.n_starts,
        maxiters=config.maxiters,
        weighted=config.weighted,
        reltol=config.reltol,
        abstol=config.abstol,
        seed=config.seed,
    )

    successful_models = [m for m in rank_result.ranking.model if haskey(rank_result.fits, m)]
    if isempty(successful_models)
        return (
            config=config,
            conditions=conditions,
            ranking=rank_result.ranking,
            failures=rank_result.failures,
            qc_report=qc_report,
            qc_paths=qc_paths,
            plots=String[],
            exports=nothing,
        )
    end

    plot_paths = plot_topk(rank_result; conditions=conditions, top_k=config.top_k, output_dir=joinpath(config.output_dir, "figures"))
    export_paths = export_results(rank_result; output_dir=config.output_dir)

    return (
        config=config,
        conditions=conditions,
        ranking=rank_result.ranking,
        failures=rank_result.failures,
        qc_report=qc_report,
        qc_paths=qc_paths,
        plots=plot_paths,
        exports=export_paths,
    )
end

end # module Workflow
