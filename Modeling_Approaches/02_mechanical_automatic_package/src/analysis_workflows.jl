module AnalysisWorkflows

using CSV
using DataFrames
using GrowthParameterEstimation
using Plots
using Statistics

using ..IOUtils

export run_condition_analysis!

_safe_label(v) = (v === missing || v === nothing || isempty(strip(string(v)))) ? "missing" : string(v)

function _param_numbers(params_text)
    return [parse(Float64, m.match) for m in eachmatch(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", string(params_text))]
end

function _predict_untreated(model_name::AbstractString, params, x::AbstractVector, y0::Real)
    length(x) == 0 && return Float64[]
    r = length(params) >= 1 ? Float64(params[1]) : 0.0
    K = length(params) >= 2 ? max(Float64(params[2]), 1e-8) : max(Float64(y0), 1.0)
    y = Float64[max(Float64(y0), 1e-8)]
    for i in 2:length(x)
        dt = Float64(x[i] - x[i - 1])
        prev = max(y[end], 1e-8)
        dydt = if model_name == "logistic_growth"
            r * prev * max(0.0, 1 - prev / K)
        elseif model_name == "gompertz_growth"
            r * prev * log(K / prev)
        elseif model_name == "logistic_simple_death"
            d = length(params) >= 3 ? Float64(params[3]) : 0.0
            r * prev * max(0.0, 1 - prev / K) - d * prev
        elseif model_name == "allee_growth"
            A = length(params) >= 3 ? max(Float64(params[3]), 1e-8) : max(0.2 * prev, 1e-8)
            r * prev * max(0.0, 1 - prev / K) * (prev / A - 1)
        elseif model_name in ("theta_logistic_growth", "generalized_logistic_growth")
            theta = length(params) >= 3 ? max(Float64(params[3]), 1e-8) : 1.0
            r * prev * max(0.0, 1 - (prev / K)^theta)
        else
            return Float64[]
        end
        push!(y, max(prev + dt * dydt, 0.0))
    end
    return y
end

function _write_untreated_multi_model_overlays(decoded::DataFrame, ranking::DataFrame, condition::AbstractString, out)
    condition in ("monoculture_untreated", "coculture_untreated") || return (path = nothing, csv_path = nothing)
    (:time in propertynames(decoded) && :count in propertynames(decoded)) || return (path = nothing, csv_path = nothing)

    group_cols = Symbol[]
    :cell_line in propertynames(decoded) && push!(group_cols, :cell_line)
    :density in propertynames(decoded) && push!(group_cols, :density)
    :mix in propertynames(decoded) && any(!isempty(_safe_label(v)) && _safe_label(v) != "missing" for v in decoded.mix) && push!(group_cols, :mix)

    groups = isempty(group_cols) ? [decoded] : collect(groupby(decoded, group_cols))
    fig_dir = joinpath(out.images, "figures")
    csv_dir = joinpath(out.csv, "figures")
    mkpath(fig_dir)
    mkpath(csv_dir)

    overlay_rows = DataFrame[]
    plot_paths = String[]
    combined = plot(
        xlabel = "Time (day)",
        ylabel = "Cell count",
        title = "$(condition) observed means with fitted model overlays",
        legend = :outertopright,
    )

    for (idx, grp) in enumerate(groups)
        group_label_parts = ["$(c)=$(_safe_label(first(grp[!, c])))" for c in group_cols]
        group_label = isempty(group_label_parts) ? "pooled" : join(group_label_parts, " | ")
        selector = trues(nrow(ranking))
        for c in group_cols
            c in propertynames(ranking) || continue
            selector .&= [string(v) == _safe_label(first(grp[!, c])) for v in ranking[!, c]]
        end
        rank_sub = ranking[selector, :]
        nrow(rank_sub) == 0 && continue
        sort!(rank_sub, :bic)
        top_models = first(rank_sub, min(3, nrow(rank_sub)))
        required_models = Set(["theta_logistic_growth", "generalized_logistic_growth"])
        missing_required = filter(row -> String(row.model) in required_models && !(String(row.model) in Set(String.(top_models.model))), rank_sub)
        if nrow(missing_required) > 0
            top_models = vcat(top_models, missing_required; cols = :union)
            sort!(top_models, :bic)
        end

        mean_df = combine(groupby(grp, :time), :count => mean => :mean_count)
        sort!(mean_df, :time)
        x = Float64.(mean_df.time)
        y = Float64.(mean_df.mean_count)
        subgroup_plot = plot(
            x,
            y;
            seriestype = :scatter,
            markerstrokewidth = 0,
            ms = 4,
            label = "observed mean",
            xlabel = "Time (day)",
            ylabel = "Cell count",
            title = group_label,
            legend = :outertopright,
        )
        scatter!(combined, x, y; ms = 3, markerstrokewidth = 0, label = "$(group_label) observed")

        for row in eachrow(top_models)
            params = _param_numbers(row.params)
            yhat = _predict_untreated(String(row.model), params, x, first(y))
            length(yhat) == length(x) || continue
            label = "$(row.model) BIC=$(round(Float64(row.bic); digits = 2))"
            plot!(subgroup_plot, x, yhat; lw = 2, label = label)
            plot!(combined, x, yhat; lw = 1.8, label = "$(group_label) $(row.model)")
            push!(overlay_rows, DataFrame(
                condition = fill(condition, length(x)),
                group = fill(group_label, length(x)),
                model = fill(String(row.model), length(x)),
                bic = fill(Float64(row.bic), length(x)),
                time = x,
                observed_mean = y,
                predicted = yhat,
            ))
        end

        safe_group = "subgroup_$(idx)"
        sub_path = joinpath(fig_dir, "$(condition)_$(safe_group)_multi_model_overlay.png")
        savefig(subgroup_plot, sub_path)
        push!(plot_paths, sub_path)
    end

    combined_path = joinpath(fig_dir, "$(condition)_multi_model_overlay_all_groups.png")
    savefig(combined, combined_path)
    push!(plot_paths, combined_path)

    overlay_csv_path = joinpath(csv_dir, "$(condition)_multi_model_overlay_predictions.csv")
    if !isempty(overlay_rows)
        CSV.write(overlay_csv_path, vcat(overlay_rows...))
    end
    return (path = combined_path, csv_path = isempty(overlay_rows) ? nothing : overlay_csv_path, plot_paths = plot_paths)
end

function run_condition_analysis!(decoded::DataFrame, fit_result, condition::AbstractString; start::AbstractString = pwd())
    out = IOUtils.condition_output_dirs(condition; start)

    ranking = fit_result.ranking
    nrow(ranking) == 0 && error("No ranking rows for condition: $(condition)")

    top = first(sort(ranking, :bic), min(3, nrow(ranking)))
    top_path = joinpath(out.metrics, "$(condition)_analysis_top3_by_bic.csv")
    CSV.write(top_path, top)

    # Sensitivity on best candidate model (available in ranking table)
    best_row = first(sort(ranking, :bic), 1)
    best_model = String(best_row.model[1])

    # Package-level sensitivity utility call; fall back gracefully if signature differs.
    sens_df = DataFrame()
    try
        sens = GrowthParameterEstimation.parameter_sensitivity_analysis(decoded; model_name = best_model)
        sens_df = DataFrame(sens)
    catch
        sens_df = DataFrame(model = [best_model], note = ["Sensitivity call not available with current signature; update wrapper."])
    end

    sens_path = joinpath(out.metrics, "$(condition)_parameter_sensitivity.csv")
    CSV.write(sens_path, sens_df)

    bic_plot = bar(
        string.(top.model),
        top.bic;
        legend = false,
        xlabel = "Model",
        ylabel = "BIC",
        title = "$(condition) automatic top models by BIC",
        xrotation = 20,
    )
    bic_plot_path = joinpath(out.images, "$(condition)_automatic_top_bic_bar.png")
    savefig(bic_plot, bic_plot_path)

    # Coverage diagnostics: ensure all experimental groups are represented and visualized.
    group_cols = Symbol[]
    :cell_line in names(decoded) && push!(group_cols, :cell_line)
    :density in names(decoded) && push!(group_cols, :density)
    :dose in names(decoded) && push!(group_cols, :dose)
    :mix in names(decoded) && any(!isempty(_safe_label(v)) && _safe_label(v) != "missing" for v in decoded.mix) && push!(group_cols, :mix)

    has_replicate = :replicate in propertynames(decoded)
    coverage_df = if isempty(group_cols)
        DataFrame(
            group = ["pooled"],
            n_rows = [nrow(decoded)],
            n_replicates = [has_replicate ? length(unique(decoded.replicate)) : missing],
        )
    elseif has_replicate
        combine(groupby(decoded, group_cols), nrow => :n_rows, :replicate => (x -> length(unique(x))) => :n_replicates)
    else
        combine(groupby(decoded, group_cols), nrow => :n_rows)
    end
    coverage_path = joinpath(out.metrics, "$(condition)_group_coverage_counts.csv")
    CSV.write(coverage_path, coverage_df)

    # Mean trajectory plot by available groups.
    p_cov = plot(
        xlabel = "Time (day)",
        ylabel = "Cell count",
        title = "$(condition) mean trajectories by group",
        legend = :outertopright,
    )
    if isempty(group_cols)
        mean_df = combine(groupby(decoded, :time), :count => mean => :mean_count)
        sort!(mean_df, :time)
        plot!(p_cov, mean_df.time, mean_df.mean_count; lw = 2.5, label = "pooled")
    else
        for grp in groupby(decoded, group_cols)
            mean_df = combine(groupby(grp, :time), :count => mean => :mean_count)
            sort!(mean_df, :time)
            label = join(["$(c)=$(_safe_label(first(grp[!, c])))" for c in group_cols], " | ")
            plot!(p_cov, mean_df.time, mean_df.mean_count; lw = 2, label = label)
        end
    end
    coverage_plot_path = joinpath(out.images, "$(condition)_group_mean_coverage_plot.png")
    savefig(p_cov, coverage_plot_path)

    untreated_overlay = _write_untreated_multi_model_overlays(decoded, ranking, condition, out)
    extra_outputs = String[]
    untreated_overlay.path !== nothing && push!(extra_outputs, untreated_overlay.path)
    untreated_overlay.csv_path !== nothing && push!(extra_outputs, untreated_overlay.csv_path)

    IOUtils.write_manifest_row(
        condition = condition,
        step = "analysis",
        outputs = vcat([top_path, sens_path, bic_plot_path, coverage_path, coverage_plot_path], extra_outputs),
        start = start,
    )

    return (
        top_path = top_path,
        sensitivity_path = sens_path,
        image_path = bic_plot_path,
        coverage_path = coverage_path,
        coverage_plot_path = coverage_plot_path,
        untreated_overlay_path = untreated_overlay.path,
        untreated_overlay_csv_path = untreated_overlay.csv_path,
        sensitivity = sens_df,
    )
end

end
