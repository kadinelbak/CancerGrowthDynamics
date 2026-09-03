module SampleAwareReport

using CSV
using DataFrames
using Plots
using Statistics

using ..IOUtils
using ..StagedA2780Workflow

export render_a2780_sample_report_html

const SAMPLE_REPORT_NAME = "a2780_sample_aware_staged_model_comparison.html"
const SAMPLE_IMAGE_DIR = "sample_aware"
const NAIVE_COLOR = :crimson
const CIS_COLOR = :steelblue
const FIT_COLOR = :black

_as_string(value) = value === missing ? "" : String(value)
_as_float(value) = Float64(value)

function _sample_summary(samples::DataFrame, keys::Vector{Symbol})
    grouped = groupby(samples, vcat(keys, [:time]))
    return combine(
        grouped,
        :count => (values -> mean(Float64.(values))) => :sample_mean,
        :count => (values -> length(values) > 1 ? std(Float64.(values); corrected = true) : 0.0) => :sample_sd,
        :count => length => :n_wells,
    )
end

function _ols_slope(times, values)
    x = Float64.(times)
    y = Float64.(values)
    length(x) >= 2 || return NaN
    x_center = x .- mean(x)
    denominator = sum(abs2, x_center)
    denominator > 0 || return NaN
    return sum(x_center .* (y .- mean(y))) / denominator
end

function _last_window(rows::AbstractDataFrame; points::Int = 5)
    ordered = sort(rows, :time)
    return last(ordered, min(points, nrow(ordered)))
end

function _trajectory_metrics(
    stage::AbstractString,
    samples::DataFrame,
    summary::DataFrame,
    predictions::DataFrame,
    keys::Vector{Symbol},
)
    joined = innerjoin(summary, predictions; on = vcat(keys, [:time]), makeunique = true)
    rows = NamedTuple[]
    for group in groupby(joined, keys)
        cv_values = [
            Float64(row.sample_sd) / abs(Float64(row.sample_mean))
            for row in eachrow(group)
            if row.n_wells > 1 && abs(Float64(row.sample_mean)) > 1e-9
        ]
        coverage = mean(abs.(Float64.(group.predicted) .- Float64.(group.sample_mean)) .<= Float64.(group.sample_sd))
        sample_group = samples
        for key in keys
            sample_group = sample_group[string.(sample_group[!, key]) .== string(first(group[!, key])), :]
        end
        late_slopes = Float64[]
        peak = max(maximum(Float64.(sample_group.count)), 1.0)
        threshold = 0.01 * peak
        for well in groupby(sample_group, :replicate)
            late = _last_window(well)
            push!(late_slopes, _ols_slope(late.time, late.count))
        end
        positive = count(slope -> isfinite(slope) && slope > threshold, late_slopes)
        negative = count(slope -> isfinite(slope) && slope < -threshold, late_slopes)
        flat = length(late_slopes) - positive - negative
        prediction_late = _last_window(select(group, vcat(keys, [:time, :predicted])))
        prediction_slope = _ols_slope(prediction_late.time, prediction_late.predicted)
        key_values = NamedTuple{Tuple(keys)}(Tuple(first(group[!, key]) for key in keys))
        push!(rows, merge(key_values, (
            stage = String(stage),
            n_wells = isempty(late_slopes) ? 0 : length(late_slopes),
            median_cv = isempty(cv_values) ? NaN : median(cv_values),
            fit_within_one_sd = coverage,
            positive_wells = positive,
            flat_wells = flat,
            negative_wells = negative,
            mean_late_slope = isempty(late_slopes) ? NaN : mean(late_slopes),
            predicted_late_slope = prediction_slope,
            direction_consensus = isempty(late_slopes) ? NaN : maximum((positive, flat, negative)) / length(late_slopes),
        )))
    end
    return DataFrame(rows)
end

function _load_stage_data(root::AbstractString)
    csv_root = joinpath(root, "outputs", "csv")

    stage1_samples = CSV.read(joinpath(csv_root, "monoculture_untreated", "monoculture_untreated_a2780_decoded.csv"), DataFrame)
    stage1_predictions = CSV.read(joinpath(csv_root, "monoculture_untreated", "figures", "monoculture_untreated_all_groups_best_overlays.csv"), DataFrame)
    rename!(stage1_predictions, :prediction => :predicted)

    stage2_samples = CSV.read(joinpath(csv_root, "monoculture_treated", "monoculture_treated_a2780_decoded.csv"), DataFrame)
    stage2_predictions = CSV.read(joinpath(csv_root, "monoculture_treated", "figures", "monoculture_treated_best_joint_model_by_environment_overlays.csv"), DataFrame)
    stage2_predictions = stage2_predictions[.!Bool.(stage2_predictions.fixed_day0_anchor), :]
    stage2_timing_predictions = CSV.read(joinpath(csv_root, "monoculture_treated", "figures", "monoculture_treated_timing_hypothesis_overlays.csv"), DataFrame)
    stage2_timing_ranking = sort!(CSV.read(joinpath(csv_root, "monoculture_treated", "monoculture_treated_timing_hypothesis_ranking.csv"), DataFrame), :bic)
    stage2_timing_winner = String(first(stage2_timing_ranking.model))

    stage3_samples = CSV.read(joinpath(csv_root, "coculture_untreated", "coculture_untreated_recovered_design.csv"), DataFrame)
    stage3_predictions = CSV.read(joinpath(csv_root, "coculture_untreated", "figures", "coculture_untreated_joint_overlays.csv"), DataFrame)
    stage3_ranking = sort!(CSV.read(joinpath(csv_root, "coculture_untreated", "coculture_untreated_pooling_top5.csv"), DataFrame), :bic)
    stage3_winner = first(stage3_ranking)
    stage3_predictions = stage3_predictions[
        (String.(stage3_predictions.model) .== String(stage3_winner.model)) .&
        (String.(stage3_predictions.pooling_mode) .== String(stage3_winner.pooling_mode)) .&
        .!Bool.(stage3_predictions.fixed_day0_anchor),
        :,
    ]

    stage4_samples = CSV.read(joinpath(csv_root, "coculture_treated", "coculture_treated_recovered_design.csv"), DataFrame)
    stage4_predictions = CSV.read(joinpath(csv_root, "coculture_treated", "figures", "linked_treatment_combined_overlays.csv"), DataFrame)
    stage4_status = CSV.read(joinpath(csv_root, "coculture_treated", "linked_treatment_status.csv"), DataFrame)
    stage4_winner = String(first(stage4_status.winning_model))
    stage4_predictions = stage4_predictions[
        (String.(stage4_predictions.model) .== stage4_winner) .&
        (String.(stage4_predictions.context) .== "coculture") .&
        .!Bool.(stage4_predictions.fixed_day0_anchor),
        :,
    ]

    return (
        stage1 = (samples = stage1_samples, predictions = stage1_predictions, keys = [:cell_line, :density]),
        stage2 = (
            samples = stage2_samples,
            predictions = stage2_predictions,
            timing_predictions = stage2_timing_predictions,
            timing_winner = stage2_timing_winner,
            keys = [:cell_line, :density, :dose],
        ),
        stage3 = (samples = stage3_samples, predictions = stage3_predictions, keys = [:density, :mix, :component]),
        stage4 = (samples = stage4_samples, predictions = stage4_predictions, keys = [:density, :mix, :component]),
    )
end

function _legend_panel(; coculture::Bool = false)
    panel = plot(
        axis = nothing,
        framestyle = :none,
        legend = :bottom,
        legend_column = 3,
        background_color = :white,
        foreground_color_legend = nothing,
        background_color_legend = nothing,
        margin = 0Plots.mm,
    )
    if coculture
        plot!(panel, [NaN], [NaN]; color = NAIVE_COLOR, alpha = 0.25, lw = 1, label = "Naive well")
        plot!(panel, [NaN], [NaN]; color = NAIVE_COLOR, fillalpha = 0.18, ribbon = [1.0], lw = 1.5, label = "Naive mean +/- SD")
        plot!(panel, [NaN], [NaN]; color = NAIVE_COLOR, lw = 2.7, label = "Naive fit")
        plot!(panel, [NaN], [NaN]; color = CIS_COLOR, alpha = 0.25, lw = 1, label = "cis well")
        plot!(panel, [NaN], [NaN]; color = CIS_COLOR, fillalpha = 0.18, ribbon = [1.0], lw = 1.5, label = "cis mean +/- SD")
        plot!(panel, [NaN], [NaN]; color = CIS_COLOR, lw = 2.7, label = "cis fit")
    else
        plot!(panel, [NaN], [NaN]; color = :gray45, alpha = 0.35, lw = 1, label = "Well-level sample")
        plot!(panel, [NaN], [NaN]; color = CIS_COLOR, fillalpha = 0.18, ribbon = [1.0], lw = 1.5, label = "Mean +/- 1 SD")
        plot!(panel, [NaN], [NaN]; color = FIT_COLOR, lw = 2.7, label = "Selected staged fit")
    end
    return panel
end

function _draw_sample_layer!(panel, sample_rows, summary_rows, prediction_rows; color, fit_color = FIT_COLOR)
    for well in groupby(sample_rows, :replicate)
        ordered = sort(well, :time)
        plot!(panel, ordered.time, ordered.count; color = color, alpha = 0.22, lw = 0.9, label = false)
    end
    summary_rows = sort(summary_rows, :time)
    lower = min.(Float64.(summary_rows.sample_sd), Float64.(summary_rows.sample_mean))
    upper = Float64.(summary_rows.sample_sd)
    plot!(
        panel,
        summary_rows.time,
        summary_rows.sample_mean;
        ribbon = (lower, upper),
        color = color,
        fillalpha = 0.16,
        lw = 1.5,
        label = false,
    )
    scatter!(panel, summary_rows.time, summary_rows.sample_mean; color = color, ms = 2.6, markerstrokewidth = 0, label = false)
    prediction_rows = sort(prediction_rows, :time)
    plot!(panel, prediction_rows.time, prediction_rows.predicted; color = fit_color, lw = 2.5, label = false)
    return panel
end

function _mono_figure(stage, summary, output_path; treated::Bool)
    panels = Any[]
    cell_lines = ("A2780Naive", "A2780cis")
    densities = ("20k", "30k")
    doses = treated ? (0.67, 1.0, 1.47) : (0.0,)
    for cell_line in cell_lines, density in densities, dose in doses
        samples = stage.samples[
            (String.(stage.samples.cell_line) .== cell_line) .&
            (String.(stage.samples.density) .== density) .&
            isapprox.(Float64.(stage.samples.dose), dose; atol = 1e-8),
            :,
        ]
        summary_rows = summary[
            (String.(summary.cell_line) .== cell_line) .&
            (String.(summary.density) .== density) .&
            (treated ? isapprox.(Float64.(summary.dose), dose; atol = 1e-8) : trues(nrow(summary))),
            :,
        ]
        predictions = stage.predictions[
            (String.(stage.predictions.cell_line) .== cell_line) .&
            (String.(stage.predictions.density) .== density) .&
            (treated ? isapprox.(Float64.(stage.predictions.dose), dose; atol = 1e-8) : trues(nrow(stage.predictions))),
            :,
        ]
        title = treated ? "$(cell_line), $(density), $(dose) uM" : "$(cell_line), $(density)"
        panel = plot(title = title, titlefontsize = 8, xlabel = "Time (day)", ylabel = "Measured population", legend = false)
        _draw_sample_layer!(panel, samples, summary_rows, predictions; color = CIS_COLOR)
        push!(panels, panel)
    end
    legend_panel = _legend_panel()
    rows = treated ? 4 : 2
    columns = treated ? 3 : 2
    layout = @layout [grid(rows, columns); legend{0.09h}]
    figure = plot(
        panels...,
        legend_panel;
        layout = layout,
        size = treated ? (1500, 1600) : (1200, 920),
        margin = 4Plots.mm,
        plot_title = treated ? "Treated monoculture: well-level samples and selected fits" : "Untreated monoculture: well-level samples and selected fits",
    )
    mkpath(dirname(output_path))
    savefig(figure, output_path)
    return output_path
end

function _timing_figure(stage, summary, output_path)
    hypothesis_order = (
        "cis_onset_only",
        "independent_onset_gradual",
        "shared_onset_gradual",
        "cis_gradual_only",
        "partial_onset_0_5day",
    )
    hypothesis_labels = Dict(
        "cis_onset_only" => "cis onset only",
        "independent_onset_gradual" => "independent onset + ramp",
        "shared_onset_gradual" => "shared onset + ramp",
        "cis_gradual_only" => "cis gradual only",
        "partial_onset_0_5day" => "onsets within 0.5 day",
    )
    alternative_colors = (:darkorange2, :seagreen3, :orchid3, :slateblue3)
    panels = Any[]
    for cell_line in ("A2780Naive", "A2780cis"), density in ("20k", "30k"), dose in (0.67, 1.0, 1.47)
        selector = (String.(stage.samples.cell_line) .== cell_line) .&
            (String.(stage.samples.density) .== density) .&
            isapprox.(Float64.(stage.samples.dose), dose; atol = 1e-8)
        summary_selector = (String.(summary.cell_line) .== cell_line) .&
            (String.(summary.density) .== density) .&
            isapprox.(Float64.(summary.dose), dose; atol = 1e-8)
        samples = stage.samples[selector, :]
        summary_rows = sort(summary[summary_selector, :], :time)
        panel = plot(
            title = "$(cell_line), $(density), $(dose) uM",
            titlefontsize = 8,
            xlabel = "Time (day)",
            ylabel = "Measured population",
            legend = false,
        )
        for well in groupby(samples, :replicate)
            ordered = sort(well, :time)
            plot!(panel, ordered.time, ordered.count; color = :gray45, alpha = 0.2, lw = 0.8, label = false)
        end
        lower = min.(Float64.(summary_rows.sample_sd), Float64.(summary_rows.sample_mean))
        plot!(panel, summary_rows.time, summary_rows.sample_mean; ribbon = (lower, Float64.(summary_rows.sample_sd)), color = CIS_COLOR, fillalpha = 0.15, lw = 1.4, label = false)
        scatter!(panel, summary_rows.time, summary_rows.sample_mean; color = CIS_COLOR, ms = 2.4, markerstrokewidth = 0, label = false)
        alternative_index = 0
        for hypothesis in hypothesis_order
            rows = stage.timing_predictions[
                (String.(stage.timing_predictions.cell_line) .== cell_line) .&
                (String.(stage.timing_predictions.density) .== density) .&
                isapprox.(Float64.(stage.timing_predictions.dose), dose; atol = 1e-8) .&
                (String.(stage.timing_predictions.timing_hypothesis) .== hypothesis),
                :,
            ]
            rows = sort(rows, :time)
            if hypothesis == stage.timing_winner
                plot!(panel, rows.time, rows.predicted; color = FIT_COLOR, lw = 2.6, label = false)
            else
                alternative_index += 1
                plot!(panel, rows.time, rows.predicted; color = alternative_colors[alternative_index], lw = 1.15, linestyle = :dash, alpha = 0.9, label = false)
            end
        end
        push!(panels, panel)
    end
    legend_panel = plot(axis = nothing, framestyle = :none, legend = :bottom, legend_column = 4, margin = 0Plots.mm)
    plot!(legend_panel, [NaN], [NaN]; color = :gray45, alpha = 0.3, lw = 1, label = "Well-level sample")
    plot!(legend_panel, [NaN], [NaN]; color = CIS_COLOR, fillalpha = 0.18, ribbon = [1.0], lw = 1.5, label = "Mean +/- 1 SD")
    plot!(legend_panel, [NaN], [NaN]; color = FIT_COLOR, lw = 2.6, label = "Winner: $(hypothesis_labels[stage.timing_winner])")
    alternative_index = 0
    for hypothesis in hypothesis_order
        hypothesis == stage.timing_winner && continue
        alternative_index += 1
        plot!(legend_panel, [NaN], [NaN]; color = alternative_colors[alternative_index], lw = 1.15, linestyle = :dash, label = hypothesis_labels[hypothesis])
    end
    layout = @layout [grid(4, 3); legend{0.1h}]
    figure = plot(
        panels...,
        legend_panel;
        layout = layout,
        size = (1500, 1640),
        margin = 4Plots.mm,
        plot_title = "Treated monoculture timing audit: wells, variability, and five hypotheses",
    )
    mkpath(dirname(output_path))
    savefig(figure, output_path)
    return output_path
end

function _coculture_figure(stage, summary, output_path; treated::Bool)
    panels = Any[]
    for density in ("20k", "30k"), mix in ("25-75", "50-50", "75-25")
        panel = plot(
            title = "$(density), mix $(mix)",
            titlefontsize = 9,
            xlabel = "Time (day)",
            ylabel = mix == "25-75" ? "Measured population" : "",
            legend = false,
        )
        for (component, color) in (("sensitive", NAIVE_COLOR), ("resistant", CIS_COLOR))
            selector = (String.(stage.samples.density) .== density) .& (String.(stage.samples.mix) .== mix) .& (String.(stage.samples.component) .== component)
            summary_selector = (String.(summary.density) .== density) .& (String.(summary.mix) .== mix) .& (String.(summary.component) .== component)
            prediction_selector = (String.(stage.predictions.density) .== density) .& (String.(stage.predictions.mix) .== mix) .& (String.(stage.predictions.component) .== component)
            _draw_sample_layer!(panel, stage.samples[selector, :], summary[summary_selector, :], stage.predictions[prediction_selector, :]; color = color, fit_color = color)
        end
        push!(panels, panel)
    end
    legend_panel = _legend_panel(; coculture = true)
    layout = @layout [grid(2, 3); legend{0.11h}]
    figure = plot(
        panels...,
        legend_panel;
        layout = layout,
        size = (1500, 940),
        margin = 4Plots.mm,
        plot_title = treated ? "Treated coculture: well-level samples, variability, and selected fit" : "Untreated coculture: well-level samples, variability, and selected fit",
    )
    mkpath(dirname(output_path))
    savefig(figure, output_path)
    return output_path
end

function _percent(value)
    return isfinite(value) ? string(round(100value; digits = 1), "%") : "not available"
end

function _stage_summary_html(metrics::DataFrame)
    stages = [
        ("Stage 1", "Untreated monoculture"),
        ("Stage 2", "Treated monoculture"),
        ("Stage 3", "Untreated coculture"),
        ("Stage 4", "Treated coculture"),
    ]
    rows = String[]
    for (stage, label) in stages
        subset = metrics[String.(metrics.stage) .== stage, :]
        well_range = "$(minimum(subset.n_wells))-$(maximum(subset.n_wells))"
        push!(rows, "<tr><td>$(stage)</td><td>$(label)</td><td>$(well_range)</td><td>$(_percent(median(subset.median_cv)))</td><td>$(_percent(mean(subset.fit_within_one_sd)))</td><td>$(_percent(median(subset.direction_consensus)))</td></tr>")
    end
    return """
<div class="table-wrap"><table>
<thead><tr><th>Stage</th><th>Experiment</th><th>Wells per trajectory</th><th>Median between-well CV</th><th>Fit inside mean +/- SD</th><th>Median late-direction consensus</th></tr></thead>
<tbody>$(join(rows, "\n"))</tbody>
</table></div>
"""
end

function _direction_label(row)
    counts = Dict("increasing" => row.positive_wells, "approximately flat" => row.flat_wells, "decreasing" => row.negative_wells)
    winner = first(sort(collect(keys(counts)); by = key -> (-counts[key], key)))
    tied = count(==(maximum(values(counts))), values(counts)) > 1
    return tied ? "mixed" : winner
end

function _stage4_direction_html(metrics::DataFrame)
    rows = String[]
    for row in eachrow(sort(metrics, [:density, :mix, :component]))
        component = String(row.component) == "sensitive" ? "A2780Naive" : "A2780cis"
        sample_direction = _direction_label(row)
        fit_direction = row.predicted_late_slope > 0 ? "increasing" : row.predicted_late_slope < 0 ? "decreasing" : "flat"
        agreement = sample_direction == fit_direction ? "direction agrees" : sample_direction == "mixed" ? "samples are mixed" : "direction differs"
        push!(rows, "<tr><td>$(row.density)</td><td>$(row.mix)</td><td>$(component)</td><td>$(row.positive_wells) / $(row.flat_wells) / $(row.negative_wells)</td><td>$(round(row.mean_late_slope; digits = 1))</td><td>$(round(row.predicted_late_slope; digits = 1))</td><td>$(agreement)</td></tr>")
    end
    return """
<div class="table-wrap"><table>
<thead><tr><th>Density</th><th>Mix</th><th>Lineage</th><th>Increasing / flat / decreasing wells</th><th>Mean well late slope</th><th>Fit late slope</th><th>Interpretation</th></tr></thead>
<tbody>$(join(rows, "\n"))</tbody>
</table></div>
"""
end

function _findings_html(metrics::DataFrame)
    stage4 = metrics[String.(metrics.stage) .== "Stage 4", :]
    directional_matches = count(row -> _direction_label(row) == (row.predicted_late_slope > 0 ? "increasing" : row.predicted_late_slope < 0 ? "decreasing" : "flat"), eachrow(stage4))
    mixed = count(row -> _direction_label(row) == "mixed", eachrow(stage4))
    mixed_label = mixed == 1 ? "trajectory has" : "trajectories have"
    return """
<section>
<div class="stage-heading"><span>Experimental findings</span><h2>What the well-level samples add</h2></div>
<p>The bands quantify reproducibility across wells after tiles have already been averaged within each well. They are <strong>plus or minus one between-well standard deviation</strong>, not confidence intervals and not uncertainty in the fitted ODE parameters.</p>
$(_stage_summary_html(metrics))
<h3>Direction of the treated-coculture samples</h3>
<p>Late direction is estimated independently in each well from its final five measurements. A slope is classified as approximately flat when its magnitude is no more than 1% of that trajectory's peak count per day. The selected Stage 4 curve agrees with the majority well direction in <strong>$(directional_matches) of $(nrow(stage4))</strong> lineage-environment trajectories; <strong>$(mixed)</strong> $(mixed_label) no unique well-level direction. This is why a smooth mean curve can look convincing while the biological direction remains uncertain.</p>
$(_stage4_direction_html(stage4))
<h3>Interpretation</h3>
<p><strong>Stages 1 and 2:</strong> the well trajectories show how much untreated growth and monoculture treatment response vary before coculture parameters are introduced. Wide bands identify environments where one shared curve should be interpreted as a population-average response rather than a typical well.</p>
<p><strong>Stage 3:</strong> agreement across wells supports the broad asymmetric-competition pattern, while disagreements near the final days limit claims about lineage-specific loss. The error bands make clear whether an apparent downturn is replicated or driven by a subset of wells.</p>
<p><strong>Stage 4:</strong> the selected linked-treatment mechanism remains the canonical BIC winner, but the sample-level directions are not uniformly reproduced. Longer observation and independent experimental repeats are needed before late decline, rebound, or tolerant-state growth is treated as established biology.</p>
<p><strong>Model-selection boundary:</strong> rankings are intentionally inherited from the validated staged analysis. Recomputing ordinary BIC after treating repeated measurements from the same well as independent would inflate the observation count and overstate evidence. A future ranking change should use a hierarchical likelihood or a well-level bootstrap.</p>
</section>
"""
end

function _replace_stage_figure(html::String, old_name::String, new_name::String, caption::String)
    relative_path = "../images/$(SAMPLE_IMAGE_DIR)/$(basename(new_name))"
    old_image_pattern = Regex("src=\"[^\"]*" * replace(old_name, "." => "\\.") * "\"")
    html = replace(html, old_image_pattern => "src=\"$(relative_path)\"")
    figure_pattern = Regex("(<img src=\"" * replace(relative_path, "." => "\\.") * "\"[^>]*><figcaption>)(.*?)(</figcaption>)")
    return replace(html, figure_pattern => SubstitutionString("\\1" * caption * "\\3"))
end

function render_a2780_sample_report_html(; start::AbstractString = pwd(), mirror_directory::Union{Nothing,AbstractString} = nothing)
    root = IOUtils.package_root(start)
    report_dir = joinpath(root, "outputs", "reports")
    image_dir = joinpath(root, "outputs", "images", SAMPLE_IMAGE_DIR)
    csv_dir = joinpath(root, "outputs", "csv", SAMPLE_IMAGE_DIR)
    mkpath.((report_dir, image_dir, csv_dir))

    stages = _load_stage_data(root)
    summaries = Dict{Symbol,DataFrame}()
    metrics_parts = DataFrame[]
    for (number, key) in enumerate((:stage1, :stage2, :stage3, :stage4))
        stage = getproperty(stages, key)
        summary = _sample_summary(stage.samples, stage.keys)
        summaries[key] = summary
        push!(metrics_parts, _trajectory_metrics("Stage $(number)", stage.samples, summary, stage.predictions, stage.keys))
    end
    metrics = vcat(metrics_parts...; cols = :union)
    CSV.write(joinpath(csv_dir, "sample_timepoint_summaries.csv"), vcat([
        hcat(DataFrame(stage = fill("Stage $(index)", nrow(summaries[key]))), summaries[key]; makeunique = true)
        for (index, key) in enumerate((:stage1, :stage2, :stage3, :stage4))
    ]...; cols = :union))
    CSV.write(joinpath(csv_dir, "sample_trajectory_findings.csv"), metrics)

    stage1_image = _mono_figure(stages.stage1, summaries[:stage1], joinpath(image_dir, "stage1_untreated_monoculture_samples.png"); treated = false)
    stage2_image = _mono_figure(stages.stage2, summaries[:stage2], joinpath(image_dir, "stage2_treated_monoculture_samples.png"); treated = true)
    stage2_timing_image = _timing_figure(stages.stage2, summaries[:stage2], joinpath(image_dir, "stage2_timing_hypotheses_samples.png"))
    stage3_image = _coculture_figure(stages.stage3, summaries[:stage3], joinpath(image_dir, "stage3_untreated_coculture_samples.png"); treated = false)
    stage4_image = _coculture_figure(stages.stage4, summaries[:stage4], joinpath(image_dir, "stage4_treated_coculture_samples.png"); treated = true)

    canonical_path = StagedA2780Workflow.render_a2780_report_html(; start)
    html = read(canonical_path, String)
    html = replace(html,
        "<title>A2780 staged model comparison</title>" => "<title>A2780 sample-aware staged model comparison</title>",
        "<h1>A2780 staged model comparison</h1>" => "<h1>A2780 sample-aware staged model comparison</h1>",
        "<p>Top-five mechanistic equations, BIC rankings, and fitted graph grids for the four-stage analysis.</p>" => "<p>The same four-stage equations and validated rankings, shown against tile-averaged well samples with between-well variability bands.</p>",
        "<p class=\"artifact-note\">Detailed parameters, diagnostics, provenance, and inheritance audits are retained in <code>outputs/csv</code>.</p>" => "<p class=\"artifact-note\"><a class=\"back-home\" href=\"a2780_staged_model_comparison.html\">Open the canonical mean-trajectory report</a><br>Detailed sample summaries and direction audits are retained in <code>outputs/csv/sample_aware</code>.</p>",
    )
    html = _replace_stage_figure(html, "monoculture_untreated_pooling_model_grid.png", "stage1_untreated_monoculture_samples.png", "Every faint line is one well after tile averaging. The band is the across-well mean plus or minus one standard deviation; the black line is the selected Stage 1 fit.")
    html = _replace_stage_figure(html, "monoculture_treated_best_joint_model_by_environment.png", "stage2_treated_monoculture_samples.png", "Treated monoculture well trajectories and between-well variability for every cell-line, density, and dose environment. The black line is the selected Stage 2 fit.")
    html = _replace_stage_figure(html, "monoculture_treated_timing_hypothesis_grid.png", "stage2_timing_hypotheses_samples.png", "The same treated-monoculture wells and variability bands, with all five timing hypotheses overlaid. The black line is the BIC winner; thinner dashed lines show the alternative onset and gradual-activation rules.")
    html = _replace_stage_figure(html, "coculture_untreated_best_mechanistic_fit_grid.png", "stage3_untreated_coculture_samples.png", "Untreated coculture well trajectories with lineage-specific mean plus or minus one standard deviation bands and the selected Stage 3 fits.")
    html = _replace_stage_figure(html, "linked_treatment_coculture_grid.png", "stage4_treated_coculture_samples.png", "Treated coculture well trajectories with lineage-specific variability bands and the selected linked Stage 4 fit.")
    html = replace(html, "<section class=\"equation-summary\">" => _findings_html(metrics) * "\n<section class=\"equation-summary\">")

    report_path = joinpath(report_dir, SAMPLE_REPORT_NAME)
    write(report_path, html)

    if mirror_directory !== nothing
        mirror_report_dir = joinpath(mirror_directory, "outputs", "reports")
        mirror_image_dir = joinpath(mirror_directory, "outputs", "images", SAMPLE_IMAGE_DIR)
        mirror_csv_dir = joinpath(mirror_directory, "outputs", "csv", SAMPLE_IMAGE_DIR)
        mkpath.((mirror_report_dir, mirror_image_dir, mirror_csv_dir))
        cp(report_path, joinpath(mirror_report_dir, SAMPLE_REPORT_NAME); force = true)
        for path in (stage1_image, stage2_image, stage2_timing_image, stage3_image, stage4_image)
            cp(path, joinpath(mirror_image_dir, basename(path)); force = true)
        end
        cp(joinpath(csv_dir, "sample_timepoint_summaries.csv"), joinpath(mirror_csv_dir, "sample_timepoint_summaries.csv"); force = true)
        cp(joinpath(csv_dir, "sample_trajectory_findings.csv"), joinpath(mirror_csv_dir, "sample_trajectory_findings.csv"); force = true)
    end
    return report_path
end

end
