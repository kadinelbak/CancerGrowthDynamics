module JointSummary

using DataFrames
using Statistics

export summarize_joint_bic, summarize_joint_bic_by_group, summarize_joint_parameter_stability,
    summarize_pooling_bic, symmetric_relative_pair

"""
    symmetric_relative_pair(center, log_contrast)

Return symmetric low/high environment values around a positive geometric center.
For `abs(log_contrast) <= log(1.05)`, each value remains within five percent of
the center.
"""
symmetric_relative_pair(center::Real, log_contrast::Real) = (
    low = center * exp(-log_contrast),
    high = center * exp(log_contrast),
)

"""
    summarize_pooling_bic(df; group_col=:cell_line, top_n=5,
        independent_mode="independent_diagnostic", inadequacy_delta=10.0)

Rank eligible shared/partial-pooling candidates by BIC within each biological
group. Independent fits are retained as diagnostics. A group is marked
`inadequate_pooling` when its best independent fit improves BIC by at least
`inadequacy_delta` over the best eligible fit.
"""
function summarize_pooling_bic(
    df::AbstractDataFrame;
    group_col::Symbol = :cell_line,
    model_col::Symbol = :model,
    pooling_col::Symbol = :pooling_mode,
    bic_col::Symbol = :bic,
    eligible_col::Symbol = :eligible_for_inheritance,
    top_n::Int = 5,
    independent_mode::AbstractString = "independent_diagnostic",
    inadequacy_delta::Real = 10.0,
)
    top_n > 0 || throw(ArgumentError("top_n must be positive"))
    _require_columns(df, [group_col, model_col, pooling_col, bic_col, eligible_col])
    ranking_parts = DataFrame[]
    status_rows = NamedTuple[]

    for grp in groupby(DataFrame(df), group_col)
        group_value = first(grp[!, group_col])
        valid = grp[[value isa Real && isfinite(Float64(value)) && abs(Float64(value)) < 9.99e11 for value in grp[!, bic_col]], :]
        eligible = valid[Bool.(valid[!, eligible_col]), :]
        isempty(eligible) && continue
        sort!(eligible, bic_col)
        best_eligible_bic = Float64(first(eligible[!, bic_col]))
        eligible.delta_bic = Float64.(eligible[!, bic_col]) .- best_eligible_bic
        eligible.rank_within_cell_line = collect(1:nrow(eligible))
        shown = first(eligible, min(top_n, nrow(eligible)))
        push!(ranking_parts, shown)

        independent = valid[String.(valid[!, pooling_col]) .== String(independent_mode), :]
        best_independent_bic = isempty(independent) ? NaN : minimum(Float64.(independent[!, bic_col]))
        improvement = isfinite(best_independent_bic) ? best_eligible_bic - best_independent_bic : NaN
        inadequate = isfinite(improvement) && improvement >= inadequacy_delta
        winner = first(eligible)
        push!(status_rows, (
            cell_line = string(group_value),
            winning_model = string(winner[model_col]),
            winning_pooling_mode = string(winner[pooling_col]),
            winning_bic = best_eligible_bic,
            best_independent_bic = best_independent_bic,
            independent_bic_improvement = improvement,
            inadequacy_delta = Float64(inadequacy_delta),
            inadequate_pooling = inadequate,
            inheritance_allowed = !inadequate,
        ))
    end

    ranking = isempty(ranking_parts) ? DataFrame() : vcat(ranking_parts...; cols = :union)
    !isempty(ranking) && sort!(ranking, [group_col, :rank_within_cell_line])
    return (ranking = ranking, status = DataFrame(status_rows))
end

function _require_columns(df::AbstractDataFrame, columns)
    available = Set(propertynames(df))
    missing_columns = setdiff(Set(Symbol.(columns)), available)
    isempty(missing_columns) || throw(ArgumentError("Missing required columns: $(sort!(collect(missing_columns)))"))
end

"""
    summarize_joint_bic_by_group(df; group_col=:cell_line, top_n=5,
                                 model_col=:model, bic_col=:bic,
                                 environment_cols=[:density])

Rank models separately within a biological grouping such as cell line while
aggregating BIC across that group's independently fitted environments. The
returned table contains at most `top_n` complete model summaries per group.
"""
function summarize_joint_bic_by_group(
    df::AbstractDataFrame;
    group_col::Symbol = :cell_line,
    top_n::Int = 5,
    model_col::Symbol = :model,
    bic_col::Symbol = :bic,
    environment_cols::Vector{Symbol} = [:density],
)
    top_n > 0 || throw(ArgumentError("top_n must be positive"))
    _require_columns(df, [group_col, model_col, bic_col, environment_cols...])
    summaries = DataFrame[]
    rank_col = Symbol("rank_within_", String(group_col))

    for grp in groupby(DataFrame(df), group_col)
        group_value = first(grp[!, group_col])
        summary = summarize_joint_bic(
            grp;
            model_col = model_col,
            bic_col = bic_col,
            environment_cols = environment_cols,
        ).aggregate
        isempty(summary) && continue
        complete = summary[summary.complete_environment_coverage, :]
        candidates = isempty(complete) ? summary : complete
        shown = first(candidates, min(top_n, nrow(candidates)))
        insertcols!(shown, 1, group_col => fill(group_value, nrow(shown)))
        insertcols!(shown, 2, rank_col => collect(1:nrow(shown)))
        push!(summaries, shown)
    end

    isempty(summaries) && return DataFrame()
    result = vcat(summaries...; cols = :union)
    sort!(result, [group_col, rank_col])
    return result
end

"""
    summarize_joint_bic(df; model_col=:model, bic_col=:bic,
                        environment_cols=[:cell_line, :density])

Summarize model BIC values across independently fitted environments. Returns
`(long, matrix, aggregate, winners)`. `aggregate.total_bic` is the sum across
complete environments, and `delta_total_bic` is measured from the best total.
"""
function summarize_joint_bic(
    df::AbstractDataFrame;
    model_col::Symbol = :model,
    bic_col::Symbol = :bic,
    environment_cols::Vector{Symbol} = [:cell_line, :density],
)
    _require_columns(df, [model_col, bic_col, environment_cols...])
    rows = NamedTuple[]
    winner_rows = NamedTuple[]

    for grp in groupby(DataFrame(df), environment_cols)
        environment = join([string(first(grp[!, col])) for col in environment_cols], " | ")
        valid = grp[[value isa Real && isfinite(Float64(value)) && abs(Float64(value)) < 9.99e11 for value in grp[!, bic_col]], :]
        isempty(valid) && continue
        ordered = sort(DataFrame(valid), bic_col)
        best_bic = Float64(first(ordered[!, bic_col]))
        for (rank, row) in enumerate(eachrow(ordered))
            push!(rows, (
                model = string(row[model_col]),
                environment = environment,
                bic = Float64(row[bic_col]),
                delta_bic = Float64(row[bic_col]) - best_bic,
                environment_rank = rank,
            ))
        end
        winner = first(ordered)
        push!(winner_rows, (
            environment = environment,
            winning_model = string(winner[model_col]),
            winning_bic = Float64(winner[bic_col]),
        ))
    end

    long = DataFrame(rows)
    isempty(long) && return (long = long, matrix = DataFrame(), aggregate = DataFrame(), winners = DataFrame(winner_rows))
    n_environments = length(unique(long.environment))
    aggregate_rows = NamedTuple[]
    for grp in groupby(long, :model)
        bics = Float64.(grp.bic)
        deltas = Float64.(grp.delta_bic)
        ranks = Float64.(grp.environment_rank)
        push!(aggregate_rows, (
            model = first(grp.model),
            environments_fit = nrow(grp),
            complete_environment_coverage = nrow(grp) == n_environments,
            total_bic = sum(bics),
            mean_bic = mean(bics),
            mean_environment_rank = mean(ranks),
            environment_wins = count(==(1.0), ranks),
            sum_environment_delta_bic = sum(deltas),
            worst_environment_delta_bic = maximum(deltas),
        ))
    end

    aggregate = DataFrame(aggregate_rows)
    complete_totals = aggregate.total_bic[aggregate.complete_environment_coverage]
    best_total = isempty(complete_totals) ? minimum(aggregate.total_bic) : minimum(complete_totals)
    aggregate.delta_total_bic = aggregate.total_bic .- best_total
    aggregate.relative_support = [
        delta <= 2 ? "substantial" :
        delta <= 6 ? "moderate" :
        delta <= 10 ? "weak" : "little"
        for delta in aggregate.delta_total_bic
    ]
    sort!(aggregate, [:delta_total_bic, :mean_environment_rank])

    matrix = unstack(select(long, :model, :environment, :bic), :model, :environment, :bic)
    matrix = leftjoin(select(aggregate, :model, :total_bic, :delta_total_bic, :mean_environment_rank, :environment_wins, :relative_support), matrix; on = :model)
    sort!(matrix, :delta_total_bic)
    winners = sort(DataFrame(winner_rows), :environment)
    return (long = long, matrix = matrix, aggregate = aggregate, winners = winners)
end

"""
    summarize_joint_parameter_stability(parameter_df)

Summarize independently fitted parameter values across environments. The input
must contain `model`, `parameter`, `value`, `lower_bound`, and `upper_bound`.
`bound_range_fraction` reports the observed cross-environment span as a
fraction of the allowed optimization interval.
"""
function summarize_joint_parameter_stability(
    parameter_df::AbstractDataFrame;
    model_col::Symbol = :model,
    parameter_col::Symbol = :parameter,
    value_col::Symbol = :value,
    lower_col::Symbol = :lower_bound,
    upper_col::Symbol = :upper_bound,
)
    _require_columns(parameter_df, [model_col, parameter_col, value_col, lower_col, upper_col])
    rows = NamedTuple[]

    for grp in groupby(DataFrame(parameter_df), [model_col, parameter_col])
        values = Float64.(grp[!, value_col])
        lower = minimum(Float64.(grp[!, lower_col]))
        upper = maximum(Float64.(grp[!, upper_col]))
        allowed_span = upper - lower
        value_min = minimum(values)
        value_max = maximum(values)
        value_mean = mean(values)
        value_sd = length(values) > 1 ? std(values) : 0.0
        observed_span = value_max - value_min
        bound_fraction = allowed_span > 0 ? observed_span / allowed_span : NaN
        boundary_margin = allowed_span > 0 ? min(value_min - lower, upper - value_max) / allowed_span : NaN
        coefficient_of_variation = abs(value_mean) > 1e-12 ? value_sd / abs(value_mean) : NaN
        stability_class = !isfinite(bound_fraction) ? "not_assessed" :
            bound_fraction <= 0.10 ? "tight" :
            bound_fraction <= 0.25 ? "moderate" : "broad"
        push!(rows, (
            model = string(first(grp[!, model_col])),
            parameter = string(first(grp[!, parameter_col])),
            environments_fit = length(values),
            value_mean = value_mean,
            value_sd = value_sd,
            value_min = value_min,
            value_max = value_max,
            observed_span = observed_span,
            coefficient_of_variation = coefficient_of_variation,
            lower_bound = lower,
            upper_bound = upper,
            bound_range_fraction = bound_fraction,
            minimum_boundary_margin_fraction = boundary_margin,
            near_optimization_bound = isfinite(boundary_margin) && boundary_margin < 0.05,
            stability_class = stability_class,
        ))
    end

    result = DataFrame(rows)
    !isempty(result) && sort!(result, [:model, :parameter])
    return result
end

end
