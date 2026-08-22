module StagedA2780Workflow

using CSV
using BlackBoxOptim
using DataFrames
using Dates
using JSON3
using OrdinaryDiffEq
using Plots
using Random
using Statistics

using ..IOUtils
using ..FitWorkflows
using ..AnalysisWorkflows
using ..ModelRegistry
using GrowthParameterEstimation

export STAGED_A2780_CONDITIONS,
       decode_a2780_condition,
       run_a2780_staged_goal!,
       refresh_a2780_output_summary!,
       render_a2780_report_html

const STAGED_A2780_CONDITIONS = [
    "monoculture_untreated",
    "monoculture_treated",
    "coculture_untreated",
    "coculture_treated",
]

const FAILURE_SENTINEL = 1e12
const BEST_ENVIRONMENT_FIGURE = "monoculture_treated_best_joint_model_by_environment.png"
const JOINT_MODEL_LABELS = Dict(
    "joint_anchored_linear_kill" => "Anchored linear kill",
    "joint_anchored_hill_kill" => "Anchored Hill kill",
    "joint_time_decay_dose_scaled" => "Time-decay kill",
    "joint_ic_effect_logistic_hill" => "IC-effect logistic Hill",
    "joint_ic_effect_hill_ramp" => "IC-effect Hill ramp",
    "joint_ic_effect_hill_ramp_onset" => "Delayed Hill ramp",
    "joint_ic_effect_transit_death" => "Transit damage/death",
    "joint_ic_effect_two_population" => "Sensitive/tolerant",
    "joint_intracellular_platinum_pkpd" => "Intracellular platinum PK/PD",
)

function _is_valid_metric(x)
    x === missing && return false
    xf = try
        Float64(x)
    catch
        return false
    end
    return isfinite(xf) && abs(xf) < FAILURE_SENTINEL * 0.999
end

function _ranking_metric_col(df::DataFrame)
    :bic in propertynames(df) && return :bic
    :aic in propertynames(df) && return :aic
    :ssr in propertynames(df) && return :ssr
    :sse in propertynames(df) && return :sse
    return nothing
end

function _valid_ranking(df::DataFrame)
    metric_col = _ranking_metric_col(df)
    metric_col === nothing && return DataFrame()
    return filter(row -> _is_valid_metric(row[metric_col]), df)
end

function _cell_text(v)
    v === missing && return ""
    return strip(String(v))
end

function _looks_like_well_label(v)
    txt = _cell_text(v)
    return match(r"^[A-Ha-h][0-9]{1,2}$", txt) !== nothing
end

function _condition_is_coculture(v)
    return occursin("coculture", lowercase(_cell_text(v)))
end

function _is_a2780_row(row)
    function val(c)
        c in propertynames(row) || return ""
        return _cell_text(row[c])
    end
    haystack = lowercase(join([
        val(:cell_line),
        val(:source_file),
        val(:condition),
        val(:mix),
    ], " "))
    return occursin("a2780", haystack)
end

function _infer_cached_metadata_from_source!(decoded::DataFrame)
    :source_file in propertynames(decoded) || return decoded
    n = nrow(decoded)
    :cell_line in propertynames(decoded) ? (decoded.cell_line = Any[decoded.cell_line[i] for i in 1:n]) : (decoded.cell_line = Any[missing for _ in 1:n])
    :density in propertynames(decoded) ? (decoded.density = Any[decoded.density[i] for i in 1:n]) : (decoded.density = Any[missing for _ in 1:n])
    :dose in propertynames(decoded) ? (decoded.dose = Any[decoded.dose[i] for i in 1:n]) : (decoded.dose = Any[missing for _ in 1:n])
    :mix in propertynames(decoded) ? (decoded.mix = Any[decoded.mix[i] for i in 1:n]) : (decoded.mix = Any[missing for _ in 1:n])
    for i in 1:n
        src = _cell_text(decoded.source_file[i])
        lower_src = lowercase(src)
        if isempty(_cell_text(decoded.cell_line[i]))
            if occursin("a2780cis", lower_src)
                decoded.cell_line[i] = "A2780cis"
            elseif occursin("a2780", lower_src)
                decoded.cell_line[i] = "A2780Naive"
            end
        end
        m = match(r"(\d{2})-(\d{2})", lower_src)
        if m !== nothing && (
            isempty(_cell_text(decoded.mix[i])) ||
            (_condition_is_coculture(decoded.condition[i]) && _looks_like_well_label(decoded.mix[i]))
        )
            decoded.mix[i] = "$(m.captures[1])-$(m.captures[2])"
        end
        if !_condition_is_coculture(decoded.condition[i]) && _looks_like_well_label(decoded.mix[i])
            decoded.mix[i] = ""
        end
        if isempty(_cell_text(decoded.density[i])) && occursin("coculture", _cell_text(decoded.condition[i]))
            decoded.density[i] = "unknown_density"
        end
        if (decoded.dose[i] === missing || Float64(decoded.dose[i]) == 0.0) && (occursin("ic50", lower_src) || occursin("1um", lower_src))
            decoded.dose[i] = 1.0
        end
    end

    # Cached coculture outputs contain pooled combined files plus the specific
    # condition files. Keep specific files for separated fitting when present.
    if any(startswith.(lowercase.(String.(decoded.source_file)), "measure_"))
        filter!(row -> startswith(lowercase(_cell_text(row.source_file)), "measure_"), decoded)
    end
    return decoded
end

function decode_a2780_condition(condition::AbstractString; start::AbstractString = pwd())
    decoded = try
        IOUtils.decode_condition_dataframe(condition; start = start)
    catch e
        out = IOUtils.condition_output_dirs(condition; start = start)
        cached_path = joinpath(out.csv, "$(condition)_automatic_decoded.csv")
        if isfile(cached_path)
            @warn "Processed_Datasets decode failed; using cached decoded output" condition cached_path exception=e
            CSV.read(cached_path, DataFrame)
        else
            rethrow(e)
        end
    end
    decoded = _infer_cached_metadata_from_source!(decoded)
    filtered = filter(_is_a2780_row, decoded)
    if nrow(filtered) == 0 && :cell_line in propertynames(decoded)
        nonblank = filter(row -> !isempty(_cell_text(row.cell_line)), decoded)
        if nrow(nonblank) > 0 && all(row -> occursin("a2780", lowercase(_cell_text(row.cell_line))), nonblank)
            filtered = nonblank
        end
    end
    if nrow(filtered) == 0 && condition in ("coculture_untreated", "coculture_treated")
        @warn "No explicit A2780 labels found in cached coculture data; retaining decoded coculture rows for diagnostics" condition
        filtered = decoded
    end
    nrow(filtered) == 0 && error("No A2780 rows decoded for condition: $(condition)")
    return filtered
end

function _coverage(decoded::DataFrame)
    cols = Symbol[]
    for c in (:condition, :cell_line, :density, :dose, :mix, :replicate)
        c in propertynames(decoded) && push!(cols, c)
    end
    isempty(cols) && return DataFrame(group = ["pooled"], n_rows = [nrow(decoded)])
    return combine(groupby(decoded, cols), nrow => :n_rows)
end

function _write_failure(out, condition::AbstractString, message::AbstractString)
    diag_dir = joinpath(out.csv, "diagnostics")
    mkpath(diag_dir)
    path = joinpath(diag_dir, "failure_report.csv")
    CSV.write(path, DataFrame(condition = [condition], message = [message], timestamp_utc = [Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS")]))
    return path
end

function _param_numbers(params_text)
    txt = String(params_text)
    return [parse(Float64, m.match) for m in eachmatch(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", txt)]
end

function _best_growth_baseline(ranking::DataFrame)
    valid = _valid_ranking(ranking)
    nrow(valid) == 0 && return nothing
    sort!(valid, _ranking_metric_col(valid))
    best = first(valid, 1)
    nums = _param_numbers(best.params[1])
    length(nums) < 2 && return nothing
    return (model = String(best.model[1]), params = nums[1:2])
end

function _measurement_columns(decoded::DataFrame)
    time_col = (:time in propertynames(decoded)) ? :time : :day
    meta_cols = Set([:dose, :time, :day, :replicate, :cell_line, :density, :sample_id, :file, :source_file, :condition, :mix])
    numeric_cols = [c for c in propertynames(decoded) if !(c in meta_cols) && eltype(decoded[!, c]) <: Real]
    area_col = :count in propertynames(decoded) ? :count : first(numeric_cols)
    return time_col, area_col
end

function _solve_observable(ode!, p, u0, x; obs = u -> u[1])
    prob = ODEProblem(ode!, u0, (x[1], x[end]), p)
    sol = solve(prob, Rodas5(); reltol = 1e-7, abstol = 1e-7, saveat = x, maxiters = 50_000)
    sol.retcode == ReturnCode.Success || error("ODE solve failed")
    return [obs(u) for u in sol.u]
end

function _bic(y, yhat, k)
    ssr = sum((y .- yhat) .^ 2)
    n = length(y)
    return n * log(max(ssr, 1e-20) / n) + k * log(n), ssr
end

function _fit_custom_treated_model(model_name, x, y, dose, r_anchor, K_anchor; max_time = 10.0)
    if model_name == "simeoni_transit_compartment"
        ode! = function (du, u, p, t)
            k_kill, k_transit = p
            live = max(u[1], 0.0)
            z2 = max(u[2], 0.0)
            z3 = max(u[3], 0.0)
            growth = r_anchor * live * max(0.0, 1 - live / max(K_anchor, 1e-8))
            drug_flux = k_kill * max(dose, 0.0) * live
            du[1] = growth - drug_flux
            du[2] = drug_flux - k_transit * z2
            du[3] = k_transit * (z2 - z3)
        end
        bounds = [(0.0, 8.0), (1e-3, 5.0)]
        p0 = [0.5, 0.5]
        obs = u -> u[1] + u[2]
        u0 = [max(y[1], 1.0), 0.0, 0.0]
    elseif model_name == "time_decay_kill_fixed"
        ode! = function (du, u, p, t)
            k_kill, lambda = p
            N = max(u[1], 0.0)
            growth = r_anchor * N * max(0.0, 1 - N / max(K_anchor, 1e-8))
            kill = k_kill * exp(-lambda * max(t - x[1], 0.0))
            du[1] = growth - kill * N
        end
        bounds = [(0.0, 4.0), (0.0, 2.0)]
        p0 = [0.2, 0.1]
        obs = u -> u[1]
        u0 = [max(y[1], 1.0)]
    elseif model_name == "time_decay_kill_dose_scaled"
        ode! = function (du, u, p, t)
            k_kill, lambda = p
            N = max(u[1], 0.0)
            growth = r_anchor * N * max(0.0, 1 - N / max(K_anchor, 1e-8))
            kill = k_kill * max(dose, 0.0) * exp(-lambda * max(t - x[1], 0.0))
            du[1] = growth - kill * N
        end
        bounds = [(0.0, 4.0), (0.0, 2.0)]
        p0 = [0.2, 0.1]
        obs = u -> u[1]
        u0 = [max(y[1], 1.0)]
    else
        error("Unknown custom treated model: $(model_name)")
    end

    function loss(p_vec)
        pv = Vector{Float64}(p_vec)
        try
            yhat = _solve_observable(ode!, pv, u0, x; obs = obs)
            return sum((y .- yhat) .^ 2)
        catch
            return FAILURE_SENTINEL
        end
    end

    result = BlackBoxOptim.bboptimize(
        loss;
        SearchRange = bounds,
        NumDimensions = length(p0),
        Method = :de_rand_1_bin,
        MaxTime = max_time,
        TraceMode = :silent,
    )
    p_opt = Vector{Float64}(result.archive_output.best_candidate)
    yhat = _solve_observable(ode!, p_opt, u0, x; obs = obs)
    bic, ssr = _bic(y, yhat, length(p_opt))
    return (params = p_opt, bic = bic, ssr = ssr, predicted = yhat)
end

function _append_treated_monoculture_literature_models!(ranking::DataFrame, decoded::DataFrame; start::AbstractString, max_time_per_fit::Float64)
    baseline_map = FitWorkflows._load_untreated_monoculture_cellline_baselines(; start = start)
    baseline_map === nothing && return ranking
    time_col, area_col = _measurement_columns(decoded)
    colset = Set(propertynames(decoded))
    cell_col = :cell_line in colset ? :cell_line : nothing
    density_col = :density in colset ? :density : nothing
    cell_col === nothing && return ranking

    gcols = Symbol[cell_col]
    density_col !== nothing && push!(gcols, density_col)
    :dose in colset && push!(gcols, :dose)

    rows = NamedTuple[]
    overlay_rows = DataFrame[]
    for grp in groupby(decoded, gcols)
        dose_val = :dose in propertynames(grp) ? Float64(first(grp.dose)) : 0.0
        means = combine(groupby(grp, time_col), area_col => mean => :observed)
        sort!(means, time_col)
        x = Float64.(means[!, time_col])
        y = Float64.(means.observed)
        cell_label = _cell_text(first(grp[!, cell_col]))
        isempty(cell_label) && (cell_label = "pooled")
        density_label = density_col === nothing ? "" : _cell_text(first(grp[!, density_col]))
        r_anchor, K_anchor = FitWorkflows._baseline_rk_for_group(baseline_map, cell_label, density_label, 0.3, max(maximum(y), 1.0))

        for model_name in ("simeoni_transit_compartment", "time_decay_kill_fixed", "time_decay_kill_dose_scaled")
            fit = _fit_custom_treated_model(model_name, x, y, dose_val, r_anchor, K_anchor; max_time = max(2.0, min(max_time_per_fit / 3, 8.0)))
            push!(rows, (
                model = model_name,
                dose = dose_val,
                bic = fit.bic,
                ssr = fit.ssr,
                params = string((fit = fit.params, r_anchor = r_anchor, K_anchor = K_anchor)),
                cell_line = cell_label,
                density = density_label,
            ))
            push!(overlay_rows, DataFrame(
                time = x,
                observed = y,
                predicted = fit.predicted,
                model = fill(model_name, length(x)),
                cell_line = fill(cell_label, length(x)),
                density = fill(density_label, length(x)),
                dose = fill(dose_val, length(x)),
            ))
        end
    end

    if !isempty(rows)
        ranking = vcat(ranking, DataFrame(rows); cols = :union)
        sort!(ranking, :bic)
    end
    if !isempty(overlay_rows)
        out = IOUtils.condition_output_dirs("monoculture_treated"; start = start)
        fig_dir = joinpath(out.csv, "figures")
        mkpath(fig_dir)
        CSV.write(joinpath(fig_dir, "monoculture_treated_literature_model_overlays.csv"), vcat(overlay_rows...))
    end
    return ranking
end

function _write_stage_manifest(rows; start::AbstractString = pwd())
    outdir = joinpath(IOUtils.package_root(start), "outputs", "reports")
    mkpath(outdir)
    path = joinpath(outdir, "a2780_staged_manifest.csv")
    CSV.write(path, DataFrame(rows))
    return path
end

function _write_coverage(condition, coverage; start)
    out = IOUtils.condition_output_dirs(condition; start = start)
    path = joinpath(out.metrics, "$(condition)_a2780_coverage.csv")
    CSV.write(path, coverage)
    return path
end

function _write_notebook_summary(summary_df; start)
    outdir = joinpath(IOUtils.package_root(start), "outputs", "reports")
    mkpath(outdir)
    path = joinpath(outdir, "a2780_stage_overview.csv")
    CSV.write(path, summary_df)
    return path
end

function run_a2780_staged_goal!(; start::AbstractString = pwd(), max_time_per_fit::Float64 = 12.0)
    manifest_rows = NamedTuple[]
    summaries = NamedTuple[]
    baseline = nothing
    prior_ok = true

    for condition in STAGED_A2780_CONDITIONS
        out = IOUtils.condition_output_dirs(condition; start = start)
        status = "started"
        message = ""
        decoded_rows = 0
        fit_rows = 0
        best_model = ""
        best_bic = NaN
        best_ssr = NaN
        outputs = String[]

        if !prior_ok
            status = "skipped"
            message = "Prior stage did not produce a finite valid fit."
            failure_path = _write_failure(out, condition, message)
            push!(outputs, failure_path)
        else
            try
                decoded = decode_a2780_condition(condition; start = start)
                decoded_rows = nrow(decoded)
                decoded_path = joinpath(out.csv, "$(condition)_a2780_decoded.csv")
                CSV.write(decoded_path, decoded)
                push!(outputs, decoded_path)

                cov = _coverage(decoded)
                push!(outputs, _write_coverage(condition, cov; start = start))

                fit = FitWorkflows.run_condition_fit!(decoded, condition; start = start, untreated_baseline = baseline, max_time_per_fit = max_time_per_fit)
                ranking = fit.ranking
                condition == "monoculture_untreated" && hasproperty(fit, :inheritance_allowed) &&
                    !fit.inheritance_allowed && error("Independent untreated fits indicate inadequate density pooling; treated inheritance stopped.")

                valid = _valid_ranking(ranking)
                nrow(valid) == 0 && error("No finite non-placeholder fits for $(condition).")
                metric_col = _ranking_metric_col(valid)
                sort!(valid, metric_col)
                fit_rows = nrow(ranking)
                best_model = String(valid.model[1])
                best_bic = :bic in propertynames(valid) ? Float64(valid.bic[1]) : NaN
                best_ssr = :ssr in propertynames(valid) ? Float64(valid.ssr[1]) : (:sse in propertynames(valid) ? Float64(valid.sse[1]) : NaN)
                push!(outputs, fit.ranking_path, fit.best_path)

                analysis = AnalysisWorkflows.run_condition_analysis!(decoded, merge(fit, (; ranking = ranking)), condition; start = start)
                push!(outputs, analysis.top_path, analysis.sensitivity_path, analysis.image_path, analysis.coverage_path, analysis.coverage_plot_path)

                if condition == "monoculture_untreated"
                    baseline = _best_growth_baseline(ranking)
                    baseline === nothing && error("Monoculture untreated did not produce r/K baseline parameters.")
                end

                status = "completed"
            catch e
                status = "failed"
                message = sprint(showerror, e)
                push!(outputs, _write_failure(out, condition, message))
                prior_ok = false
            end
        end

        push!(manifest_rows, (
            timestamp_utc = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
            condition = condition,
            status = status,
            message = message,
            output_count = length(outputs),
            outputs = join(outputs, ";"),
        ))
        push!(summaries, (
            condition = condition,
            status = status,
            decoded_rows = decoded_rows,
            fit_rows = fit_rows,
            best_model = best_model,
            best_bic = best_bic,
            best_ssr = best_ssr,
            message = message,
        ))
    end

    manifest_path = _write_stage_manifest(manifest_rows; start = start)
    overview = DataFrame(summaries)
    overview_path = _write_notebook_summary(overview; start = start)
    return (overview = overview, overview_path = overview_path, manifest_path = manifest_path)
end

function _html_escape(s)
    return replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
end

function _table_html(df::DataFrame; limit::Int = 50)
    shown = first(df, min(limit, nrow(df)))
    headers = join(["<th>$(_html_escape(c))</th>" for c in names(shown)], "")
    rows = String[]
    for r in eachrow(shown)
        push!(rows, "<tr>" * join(["<td>$(_html_escape(r[c]))</td>" for c in names(shown)], "") * "</tr>")
    end
    return "<table><thead><tr>$(headers)</tr></thead><tbody>$(join(rows, "\n"))</tbody></table>"
end

function _copy_plot_refs(condition; start)
    out = IOUtils.condition_output_dirs(condition; start = start)
    imgs = isdir(out.images) ? sort(filter(f -> endswith(lowercase(f), ".png"), readdir(out.images; join = true))) : String[]
    fig_dir = joinpath(out.images, "figures")
    if isdir(fig_dir)
        append!(imgs, sort(filter(f -> endswith(lowercase(f), ".png"), readdir(fig_dir; join = true))))
    end
    return filter(img -> basename(img) != BEST_ENVIRONMENT_FIGURE, imgs)
end

function _render_a2780_report_html_detailed(; start::AbstractString = pwd())
    report_dir = joinpath(IOUtils.package_root(start), "outputs", "reports")
    mkpath(report_dir)
    overview_path = joinpath(report_dir, "a2780_stage_overview.csv")
    overview = isfile(overview_path) ? CSV.read(overview_path, DataFrame) : DataFrame()

    sections = String[]
    push!(sections, "<h1>A2780 staged growth estimation</h1>")
    push!(sections, "<p>Stages: monoculture untreated, monoculture treated, coculture untreated, coculture treated.</p>")
    push!(sections, "<p>Custom treated monoculture models include Simeoni-style transit compartments, time-decay treatment effects, IC-effect Hill ramps, delayed extinction, and separate sensitive/tolerant populations. Joint fits use the local GrowthParameterEstimation package.</p>")
    push!(sections, "<p><strong>Dose encoding:</strong> the corrected experiment mapping is 0.67 uM to IC25, 1.0 uM to IC50, and 1.47 uM to IC75. IC-aware models fit effect levels 0.25, 0.50, and 0.75 while retaining concentration for display.</p>")
    !isempty(overview) && push!(sections, "<h2>Stage overview</h2>" * _table_html(overview))

    for condition in STAGED_A2780_CONDITIONS
        out = IOUtils.condition_output_dirs(condition; start = start)
        push!(sections, "<h2>$(_html_escape(condition))</h2>")
        cov_path = joinpath(out.metrics, "$(condition)_a2780_coverage.csv")
        isfile(cov_path) && condition != "monoculture_treated" && push!(sections, "<h3>Coverage</h3>" * _table_html(CSV.read(cov_path, DataFrame)))
        rank_path = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
        if isfile(rank_path) && !(condition in STAGED_A2780_CONDITIONS)
            rank = CSV.read(rank_path, DataFrame)
            metric_col = _ranking_metric_col(rank)
            metric_col !== nothing && sort!(rank, metric_col)
            push!(sections, "<h3>Model ranking</h3>" * _table_html(rank; limit = 15))
        end
        if condition == "monoculture_untreated"
            top5_path = joinpath(out.csv, "monoculture_untreated_pooling_top5.csv")
            status_path = joinpath(out.csv, "monoculture_untreated_pooling_status.csv")
            parameters_path = joinpath(out.csv, "monoculture_untreated_pooling_parameter_estimates.csv")
            initial_path = joinpath(out.csv, "monoculture_untreated_initial_condition_diagnostics.csv")
            isfile(top5_path) && push!(sections, "<h3>Top five model and pooling combinations per cell line</h3><p>Each growth family is fitted jointly to the 20k and 30k trajectories. Shared and symmetric plus-or-minus 5% partial pooling are eligible; independent fits are diagnostic only.</p>" * _table_html(CSV.read(top5_path, DataFrame); limit = 10))
            isfile(status_path) && push!(sections, "<h3>Pooling adequacy and inheritance</h3>" * _table_html(CSV.read(status_path, DataFrame); limit = 10))
            if isfile(parameters_path)
                parameters = CSV.read(parameters_path, DataFrame)
                shown_columns = filter(col -> col in propertynames(parameters), [:cell_line, :model, :pooling_mode, :density, :parameter, :center_value, :effective_value, :deviation_percent, :parameter_scope])
                push!(sections, "<h3>Center and effective density parameters</h3>" * _table_html(select(parameters, shown_columns); limit = 40))
            end
            if isfile(initial_path)
                initial = CSV.read(initial_path, DataFrame)
                shown_columns = filter(col -> col in propertynames(initial), [:cell_line, :density, :nominal_density, :fixed_u0, :u0_time_day, :first_observed_value, :first_observed_to_fixed_u0, :K, :u0_over_K, :residual_scale, :u0_strategy])
                push!(sections, "<h3>Starting-density diagnostics</h3><p>Each trajectory starts at day 0 with the experiment-design count: 67 for 20k and 100 for 30k. The first measured point remains an observation, not an estimated or substituted initial state. Density changes state evolution through N/K without requiring a different intrinsic growth law.</p>" * _table_html(select(initial, shown_columns); limit = 20))
            end
            pooling_figure = joinpath(out.images, "figures", "monoculture_untreated_pooling_model_grid.png")
            if isfile(pooling_figure)
                rel = replace(relpath(pooling_figure, report_dir), "\\" => "/")
                push!(sections, "<h3>Untreated model overlay graph table</h3><figure><img src=\"$rel\" alt=\"Untreated model and pooling overlays\"><figcaption>Observed means with the five leading eligible model and pooling combinations in each cell-line and density environment.</figcaption></figure>")
            end
        end
        if condition == "monoculture_treated"
            joint_path = joinpath(out.csv, "monoculture_treated_joint_dose_model_ranking.csv")
            if isfile(joint_path)
                joint_rank = CSV.read(joint_path, DataFrame)
                top5_path = joinpath(out.csv, "monoculture_treated_joint_cell_line_top5.csv")
                top5 = isfile(top5_path) ?
                    CSV.read(top5_path, DataFrame) :
                    GrowthParameterEstimation.summarize_joint_bic_by_group(
                        joint_rank;
                        group_col = :cell_line,
                        top_n = 5,
                        environment_cols = [:density],
                    )
                display_columns = [
                    :cell_line,
                    :rank_within_cell_line,
                    :model,
                    :pooling_mode,
                    :bic,
                    :delta_bic,
                    :scaled_ssr,
                    :ssr,
                    :n_parameters,
                    :growth_family,
                    :effect_parameter_scope,
                    :boundary_issue,
                ]
                shown_columns = filter(col -> col in propertynames(top5), display_columns)
                shown = select(top5, shown_columns)
                :model in propertynames(shown) && (shown.model = [get(JOINT_MODEL_LABELS, string(model), string(model)) for model in shown.model])
                rename!(shown, Dict(
                    :rank_within_cell_line => :rank,
                    :bic => :BIC,
                    :delta_bic => :delta_BIC,
                ))
                push!(sections, "<h3>Top five density-coupled models by cell line</h3><p>Each candidate is fitted simultaneously to 20k and 30k at IC25, IC50, and IC75. The exact winning untreated growth family and density-specific r/K are fixed; eligible treatment fits compare shared effects against plus-or-minus 5% amplitude pooling. Trajectory-peak residual scaling prevents larger-count curves from dominating BIC.</p>" * _table_html(shown; limit = 10))
                status_path = joinpath(out.csv, "monoculture_treated_pooling_status.csv")
                isfile(status_path) && push!(sections, "<h3>Treated pooling adequacy</h3>" * _table_html(CSV.read(status_path, DataFrame); limit = 10))
                parameter_path = joinpath(out.csv, "monoculture_treated_joint_parameter_estimates.csv")
                if isfile(parameter_path)
                    parameters = CSV.read(parameter_path, DataFrame)
                    shown_parameter_columns = filter(col -> col in propertynames(parameters), [:cell_line, :model, :pooling_mode, :density, :parameter, :center_value, :effective_value, :deviation_pct, :bound_position, :parameter_scope])
                    push!(sections, "<h3>Treatment center and effective parameters</h3>" * _table_html(select(parameters, shown_parameter_columns); limit = 50))
                end
                u0_path = joinpath(out.csv, "monoculture_treated_joint_initial_condition_diagnostics.csv")
                if isfile(u0_path)
                    u0_diagnostics = CSV.read(u0_path, DataFrame)
                    shown_u0_columns = filter(col -> col in propertynames(u0_diagnostics), [:cell_line, :density, :dose, :ic_label, :nominal_density, :fixed_u0, :u0_time_day, :first_observed_value, :u0_over_K, :residual_scale, :nominal_density_ratio, :fixed_u0_ratio, :first_observed_ratio, :initial_condition_strategy])
                    push!(sections, "<h3>Initial-condition controls</h3><p>Every treated trajectory starts at day 0 with 67 cells for 20k or 100 cells for 30k; u0 is fixed and is not counted as a kinetic parameter. The exact untreated growth law accounts for density through N/K, while treatment amplitudes are shared or constrained within five percent.</p>" * _table_html(select(u0_diagnostics, shown_u0_columns); limit = 20))
                end
                inheritance_path = joinpath(out.csv, "monoculture_treated_inheritance_audit.csv")
                isfile(inheritance_path) && push!(sections, "<h3>Untreated-to-treated inheritance audit</h3><p>This table records the exact untreated winning family and density-specific r, K, and shape value fixed in every treated fit.</p>" * _table_html(CSV.read(inheritance_path, DataFrame); limit = 20))
                best_figure_path = joinpath(out.images, "figures", BEST_ENVIRONMENT_FIGURE)
                if isfile(best_figure_path)
                    rel = replace(relpath(best_figure_path, report_dir), "\\" => "/")
                    push!(sections, "<h3>Mechanistic best-fit graph table</h3><p>Rows are cell-line and starting-density trajectories; columns are IC25, IC50, and IC75. Both density rows for a cell line use the same density-coupled BIC winner. Anchored Hill cannot reproduce a rise followed by decline because its drug effect is time-invariant; delayed Hill-ramp and transit models can.</p><figure><img src=\"$rel\" alt=\"Mechanistic best-fit graph table\"><figcaption>Cell-line-level mechanistic BIC winners fitted across both starting densities and all three doses.</figcaption></figure>")
                end
            end
        end
        if condition in ("coculture_untreated", "coculture_treated")
            top5_path = joinpath(out.csv, "$(condition)_pooling_top5.csv")
            status_path = joinpath(out.csv, "$(condition)_pooling_status.csv")
            parameter_path = joinpath(out.csv, "$(condition)_joint_parameter_estimates.csv")
            initial_path = joinpath(out.csv, "$(condition)_initial_mix_diagnostics.csv")
            provenance_path = joinpath(out.csv, "$(condition)_density_provenance_validation.csv")
            if isfile(top5_path)
                top5 = CSV.read(top5_path, DataFrame)
                shown_columns = filter(column -> column in propertynames(top5), [:cell_line, :rank_within_cell_line, :model, :pooling_mode, :bic, :delta_bic, :scaled_ssr, :ssr, :n_parameters, :boundary_issue, :parameter_scope, :untreated_interaction_model])
                push!(sections, "<h3>Top five whole-system mechanistic models</h3><p>Each candidate is one coupled sensitive/resistant model fitted simultaneously across 20k and 30k and all three mix ratios. Each lineage uses its exact untreated-monoculture growth family, r, and K; coculture does not rescale carrying capacity. Shared and symmetric plus-or-minus 5% pooling apply only to interaction or treatment parameters, while independent density fits remain diagnostic.</p>" * _table_html(select(top5, shown_columns); limit = 5))
            end
            isfile(status_path) && push!(sections, "<h3>Cross-density pooling adequacy</h3>" * _table_html(CSV.read(status_path, DataFrame); limit = 5))
            if isfile(initial_path)
                initial = CSV.read(initial_path, DataFrame)
                shown_columns = filter(column -> column in propertynames(initial), [:density, :mix, :nominal_sensitive_fraction, :fixed_total_u0, :fixed_sensitive_u0, :fixed_resistant_u0, :first_observed_sensitive, :first_observed_resistant, :first_observed_total, :observed_sensitive_fraction, :fraction_difference, :initial_condition_strategy])
                push!(sections, "<h3>Starting population and mixture diagnostics</h3><p>The two component populations are separate ODE states. At day 0 their total is fixed to 67 for 20k or 100 for 30k and split by the nominal mixture proportion. First measured values remain observations and are shown only as a consistency diagnostic.</p>" * _table_html(select(initial, shown_columns); limit = 6))
            end
            inheritance_name = condition == "coculture_untreated" ? "coculture_untreated_monoculture_inheritance_audit.csv" : "coculture_treated_inheritance_audit.csv"
            inheritance_path = joinpath(out.csv, inheritance_name)
            if isfile(inheritance_path)
                inheritance_description = condition == "coculture_untreated" ?
                    "Exact density-specific monoculture family, r, K, and shape values fixed for both sensitive and resistant states, with no carrying-capacity rescaling." :
                    "Exact untreated-coculture interaction family, pooling mode, and interaction parameters fixed before adding treatment effects."
                push!(sections, "<h3>Baseline inheritance audit</h3><p>$(inheritance_description)</p>" * _table_html(CSV.read(inheritance_path, DataFrame); limit = 30))
            end
            if isfile(provenance_path)
                provenance = CSV.read(provenance_path, DataFrame)
                max_error = :absolute_check_error in propertynames(provenance) ? maximum(filter(isfinite, Float64.(coalesce.(provenance.absolute_check_error, NaN))); init = 0.0) : NaN
                push!(sections, "<p><strong>Density recovery check:</strong> cached basenames had lost their parent 20k/30k folders. Replicate blocks recover those labels; for treated data the maximum difference from the independent day aggregate is $(round(max_error; digits = 4)).</p>")
            end
            if isfile(parameter_path) && isfile(top5_path)
                parameters = CSV.read(parameter_path, DataFrame)
                top5 = CSV.read(top5_path, DataFrame)
                eligible = :eligible_for_inheritance in propertynames(top5) ? top5[Bool.(top5.eligible_for_inheritance), :] : top5
                winner = eligible[argmin(eligible.bic), :]
                winner_parameters = parameters[
                    (String.(parameters.model) .== String(winner.model)) .&
                    (String.(parameters.pooling_mode) .== String(winner.pooling_mode)),
                    :,
                ]
                push!(sections, "<h3>Winning effective parameters</h3>" * _table_html(winner_parameters; limit = 30))
            end
            if condition == "coculture_treated"
                linked_ranking_path = joinpath(out.csv, "linked_treatment_top5.csv")
                linked_status_path = joinpath(out.csv, "linked_treatment_status.csv")
                linked_shared_path = joinpath(out.csv, "linked_treatment_shared_parameter_audit.csv")
                linked_effective_path = joinpath(out.csv, "linked_treatment_effective_parameter_inheritance.csv")
                linked_identifiability_path = joinpath(out.csv, "linked_treatment_identifiability.csv")
                linked_provenance_path = joinpath(out.csv, "linked_treatment_data_provenance.csv")
                if isfile(linked_ranking_path)
                    push!(sections, "<h3>Primary linked monoculture-coculture treatment analysis</h3><p>This is the primary stage-four analysis. One combined objective fits all 12 treated-monoculture trajectories and all 12 treated-coculture component trajectories. A2780Naive retains its delayed Hill-ramp treatment module and A2780cis retains its sensitive/tolerant module. Their intrinsic drug parameters use the same parameter indices in both contexts. Coculture candidates may add only explicit competition/load modifiers; the fully free context model is diagnostic and cannot win inheritance.</p>" * _table_html(CSV.read(linked_ranking_path, DataFrame); limit = 5))
                end
                if isfile(linked_provenance_path)
                    provenance = CSV.read(linked_provenance_path, DataFrame)
                    shown_columns = filter(
                        column -> column in propertynames(provenance),
                        [:artifact, :exists, :bytes, :sha256, :source_mode, :canonical_processed_data_present, :canonical_source_commit, :treated_monoculture_label_correction],
                    )
                    push!(sections, "<h3>Canonical data provenance</h3><p>The source tree was restored from Git history. SHA-256 hashes identify the decoded artifacts and fitted stage-two inputs used by the linked objective.</p>" * _table_html(select(provenance, shown_columns); limit = 10))
                end
                isfile(linked_status_path) && push!(sections, "<h3>Drug-inheritance adequacy</h3><p>If the fully free diagnostic improves BIC by at least 10, inherited drug response is flagged as inadequate rather than silently replaced.</p>" * _table_html(CSV.read(linked_status_path, DataFrame); limit = 2))
                isfile(linked_shared_path) && push!(sections, "<h3>Shared intrinsic drug parameters</h3><p>Sequential stage-two estimates initialize the combined fit. The joint estimate is informed by all treated data while remaining one shared intrinsic parameter across monoculture and coculture.</p>" * _table_html(CSV.read(linked_shared_path, DataFrame); limit = 30))
                isfile(linked_effective_path) && push!(sections, "<h3>Exact parameter-index inheritance audit</h3><p>A ratio of one confirms that monoculture and coculture call the same intrinsic parameter after density adjustment. Context effects are represented separately by mechanistic modifiers.</p>" * _table_html(CSV.read(linked_effective_path, DataFrame); limit = 30))
                isfile(linked_identifiability_path) && push!(sections, "<h3>Linked-fit boundary assessment</h3><p>Positive upper limits and explicitly allowed negative modifier limits are profiled in both directions. Expanded bounds are accepted only for a BIC improvement of at least two and an estimate at least five percent inside the new interval.</p>" * _table_html(CSV.read(linked_identifiability_path, DataFrame); limit = 30))
                linked_figures = [
                    ("linked_treatment_monoculture_grid.png", "Linked treated-monoculture graph table", "Rows are cell line and density; columns are IC25, IC50, and IC75."),
                    ("linked_treatment_coculture_grid.png", "Linked treated-coculture graph table", "Rows are density and columns are mixture; both lineage observations and predictions are shown."),
                    ("linked_treatment_hypothesis_comparison_grid.png", "Strict inheritance versus selected modifier", "Dashed curves use strict transferred drug response; solid curves use the best eligible coculture modifier."),
                ]
                for (filename, title, caption) in linked_figures
                    path = joinpath(out.images, "figures", filename)
                    if isfile(path)
                        rel = replace(relpath(path, report_dir), "\\" => "/")
                        push!(sections, "<h3>$(title)</h3><figure><img src=\"$rel\" alt=\"$(title)\"><figcaption>$(caption)</figcaption></figure>")
                    end
                end
                comparison_path = joinpath(out.csv, "coculture_treated_nonadditive_model_comparison.csv")
                if isfile(comparison_path)
                    comparison = CSV.read(comparison_path, DataFrame)
                    push!(sections, "<h3>Exploratory coculture-only treatment analysis</h3><p>This earlier comparison re-estimates treatment parameters from the single-dose coculture data alone. It is retained for sensitivity analysis but is no longer the primary mechanistic result.</p>" * _table_html(comparison; limit = 3))
                end
                comparison_figure_path = joinpath(out.images, "figures", "coculture_treated_nonadditive_simulation_grid.png")
                if isfile(comparison_figure_path)
                    rel = replace(relpath(comparison_figure_path, report_dir), "\\" => "/")
                    push!(sections, "<h3>Non-additive simulation graph table</h3><figure><img src=\"$rel\" alt=\"Additive and interaction-scaled treated coculture simulations\"><figcaption>Each line is generated from a whole-cell-line joint fit across both seeding densities and all mix ratios. Colors identify cell populations; line styles identify the additive, competitor-scaled, and load-scaled treatment hypotheses.</figcaption></figure>")
                end
            end
            identifiability_path = joinpath(out.csv, "$(condition)_identifiability.csv")
            if isfile(identifiability_path) && isfile(top5_path)
                identifiability = CSV.read(identifiability_path, DataFrame)
                top5 = CSV.read(top5_path, DataFrame)
                eligible = :eligible_for_inheritance in propertynames(top5) ? top5[Bool.(top5.eligible_for_inheritance), :] : top5
                winner = eligible[argmin(eligible.bic), :]
                winner_identifiability = identifiability[
                    (String.(identifiability.model) .== String(winner.model)) .&
                    (String.(identifiability.pooling_mode) .== String(winner.pooling_mode)),
                    :,
                ]
                push!(sections, "<h3>Boundary and identifiability assessment</h3><p>Upper bounds were profiled at 1.5x and 2x when estimates landed within two percent of a bound. Expansion was accepted only for a BIC improvement of at least two and a move at least five percent inside the new bound.</p>" * _table_html(winner_identifiability; limit = 20))
            end
            figure_path = joinpath(out.images, "figures", "$(condition)_best_mechanistic_fit_grid.png")
            if isfile(figure_path)
                rel = replace(relpath(figure_path, report_dir), "\\" => "/")
                push!(sections, "<h3>Best mechanistic fit graph table</h3><figure><img src=\"$rel\" alt=\"$(_html_escape(condition)) best mechanistic fit grid\"><figcaption>Rows are seeding densities and columns are nominal mix ratios. Sensitive and cis-resistant observations and model predictions are shown together in every environment.</figcaption></figure>")
            end
        end
    end

    html = """
<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>A2780 staged growth estimation</title>
<style>
body { font-family: Segoe UI, Arial, sans-serif; margin: 32px; color: #1d2428; }
table { border-collapse: collapse; margin: 12px 0 28px; width: 100%; font-size: 13px; }
th, td { border: 1px solid #cfd8dc; padding: 6px 8px; text-align: left; vertical-align: top; }
th { background: #eef4f5; }
img { display: block; max-width: 100%; height: auto; border: 1px solid #ccd4d6; }
figure { margin: 18px 0 28px; }
figcaption { color: #536267; font-size: 12px; margin-top: 6px; overflow-wrap: anywhere; }
.graph-grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 20px; align-items: start; margin: 12px 0 32px; }
.graph-grid figure { min-width: 0; margin: 0; }
@media (max-width: 640px) {
  body { margin: 16px; }
  .graph-grid { grid-template-columns: minmax(0, 1fr); }
}
    .back-home { display: inline-block; margin: 0 0 18px; color: var(--accent, #2563eb); font-weight: 700; text-decoration: none; }
    .back-home:hover { text-decoration: underline; }
</style>
</head>
<body>
$(join(sections, "\n"))
</body>
</html>
"""
    html_path = joinpath(report_dir, "a2780_staged_model_comparison.html")
    write(html_path, html)
    return html_path
end

const REPORT_MODEL_LABELS = Dict(
    "logistic_growth" => "Logistic growth",
    "theta_logistic_growth" => "Theta-logistic growth",
    "lagged_theta_logistic_growth" => "Hard-lag theta-logistic growth",
    "baranyi_theta_logistic_growth" => "Baranyi-adjusted theta-logistic growth",
    "adaptation_theta_logistic_growth" => "Smooth adaptation theta-logistic growth",
    "gompertz_growth" => "Gompertz growth",
    "joint_ic_effect_hill_ramp_onset" => "Delayed Hill-ramp kill",
    "joint_ic_effect_transit_death" => "Transit damage/death",
    "joint_ic_effect_two_population" => "Sensitive/tolerant populations",
    "joint_intracellular_platinum_pkpd" => "Intracellular platinum PK/PD",
    "independent_onset_gradual" => "Independent onset plus gradual activation",
    "shared_onset_gradual" => "Shared onset plus gradual activation",
    "partial_onset_0_5day" => "Onsets within 0.5 day plus gradual activation",
    "cis_gradual_only" => "Resistant gradual activation only",
    "cis_onset_only" => "Resistant onset only",
    "lv_symmetric_competition" => "Symmetric competition",
    "lv_asymmetric_competition" => "Asymmetric competition",
    "lv_asymmetric_competition_death" => "Asymmetric competition with death",
    "strict_inheritance" => "Strict drug-effect inheritance",
    "competitor_scaled" => "Competitor-scaled drug effect",
    "load_scaled" => "Total-load-scaled drug effect",
    "load_plus_context_amplitude" => "Load-scaled plus context amplitude",
    "tolerant_context_shift" => "Coculture tolerant-state shift",
    "subpopulation_load_scaled" => "Subpopulation-specific load scaling",
    "load_plus_tolerant_context" => "Load scaling plus tolerant-state shift",
    "load_plus_tolerant_growth_context" => "Load plus tolerant-state and growth plasticity",
    "fully_free_context_diagnostic" => "Fully free context diagnostic",
)

const REPORT_MODEL_EQUATIONS = Dict(
    "logistic_growth" => raw"\(\displaystyle \frac{dX}{dt}=rX\left(1-\frac{X}{K}\right)\)",
    "theta_logistic_growth" => raw"\(\displaystyle \frac{dX}{dt}=rX\left[1-\left(\frac{X}{K}\right)^{\theta}\right]\)",
    "lagged_theta_logistic_growth" => raw"\(\displaystyle \frac{dX}{dt}=\mathbf{1}_{t>\tau}\,rX\left[1-\left(\frac{X}{K}\right)^{\theta}\right]\)",
    "baranyi_theta_logistic_growth" => raw"\(\displaystyle \frac{dX}{dt}=\alpha_B(t)rX\left[1-\left(\frac{X}{K}\right)^{\theta}\right],\quad \alpha_B(t)=\frac{q_0}{q_0+e^{-rt}}\)",
    "adaptation_theta_logistic_growth" => raw"\(\displaystyle \frac{dX}{dt}=A(t)rX\left[1-\left(\frac{X}{K}\right)^{\theta}\right],\quad A(t)=1-e^{-\lambda_A t}\)",
    "gompertz_growth" => raw"\(\displaystyle \frac{dX}{dt}=rX\ln\!\left(\frac{K}{X}\right)\)",
    "joint_ic_effect_hill_ramp_onset" => raw"\(\displaystyle \frac{dX}{dt}=G_i(X)-A_i(t)H_i(z)X,\quad A_i(t)=\mathbf{1}_{t>t_{\mathrm{on},i}}\left[1-e^{-\lambda_i(t-t_{\mathrm{on},i})}\right]\)",
    "joint_ic_effect_transit_death" => raw"\(\displaystyle \frac{dP}{dt}=G_i(P)-A_i(t)H_i(z)P,\quad \frac{dD}{dt}=A_i(t)H_i(z)P-k_{\mathrm{clear}}D,\quad \widehat y=P+\tfrac12D\)",
    "joint_ic_effect_two_population" => raw"\(\displaystyle \frac{dS}{dt}=G(S+T)\frac{S}{S+T}-A_C(t)H_{CS}(z)S,\quad \frac{dT}{dt}=G(S+T)\frac{T}{S+T}-A_C(t)H_{CT}(z)T\)",
    "joint_intracellular_platinum_pkpd" => raw"\(\displaystyle \frac{dc_i}{dt}=D-k_{\mathrm{efflux}}c_i,\quad \frac{dc_k}{dt}=c_i-k_{\mathrm{repair}}c_k,\quad \frac{dX}{dt}=G_i(X)-H(c_k)X\)",
    "independent_onset_gradual" => raw"\(A_N(t)=R(t;\lambda_N,t_{\mathrm{on},N}),\quad A_C(t)=R(t;\lambda_C,t_{\mathrm{on},C})\)",
    "shared_onset_gradual" => raw"\(A_N(t)=R(t;\lambda_N,t_{\mathrm{on}}),\quad A_C(t)=R(t;\lambda_C,t_{\mathrm{on}})\)",
    "partial_onset_0_5day" => raw"\(t_{\mathrm{on},N}=\bar t_{\mathrm{on}}-\delta_t,\quad t_{\mathrm{on},C}=\bar t_{\mathrm{on}}+\delta_t,\quad |\delta_t|\le 0.5\,\mathrm{d}\)",
    "cis_gradual_only" => raw"\(A_N(t)=R(t;\lambda_N,t_{\mathrm{on},N}),\quad A_C(t)=1-e^{-\lambda_Ct}\)",
    "cis_onset_only" => raw"\(A_N(t)=R(t;\lambda_N,t_{\mathrm{on},N}),\quad A_C(t)=\mathbf{1}_{t>t_{\mathrm{on},C}}\)",
    "lv_symmetric_competition" => raw"\(\displaystyle \frac{dN}{dt}=G_N(N,N+\alpha C),\quad \frac{dC}{dt}=G_C(C,C+\alpha N)\)",
    "lv_asymmetric_competition" => raw"\(\displaystyle \frac{dN}{dt}=G_N(N,N+\alpha_{NC}C),\quad \frac{dC}{dt}=G_C(C,C+\alpha_{CN}N)\)",
    "lv_asymmetric_competition_death" => raw"\(\displaystyle \frac{dN}{dt}=G_N(N,N+\alpha_{NC}C)-d_NN,\quad \frac{dC}{dt}=G_C(C,C+\alpha_{CN}N)-d_CC\)",
    "strict_inheritance" => raw"\(\displaystyle \frac{dX_i}{dt}=G_i^{\mathrm{co}}-E_i^{\mathrm{mono}}(t,z)X_i\)",
    "competitor_scaled" => raw"\(\displaystyle \frac{dX_i}{dt}=G_i^{\mathrm{co}}-e^{\beta_i\alpha_{ij}X_j/K_i}E_i^{\mathrm{mono}}(t,z)X_i\)",
    "load_scaled" => raw"\(\displaystyle \frac{dX_i}{dt}=G_i^{\mathrm{co}}-e^{\beta_iL_i/K_i}E_i^{\mathrm{mono}}(t,z)X_i\)",
    "load_plus_context_amplitude" => raw"\(\displaystyle \frac{dX_i}{dt}=G_i^{\mathrm{co}}-e^{\beta_iL_i/K_i+\gamma_i}E_i^{\mathrm{mono}}(t,z)X_i\)",
    "tolerant_context_shift" => raw"\(f_{T0}^{\mathrm{co}}=\operatorname{logit}^{-1}[\operatorname{logit}(f_{T0})+\delta_f],\quad H_{CT}^{\mathrm{co}}=e^{\gamma_T}H_{CT}\)",
    "subpopulation_load_scaled" => raw"\(M_N=e^{\beta_NL_N/K_N},\quad M_{CS}=e^{\beta_{CS}L_C/K_C},\quad M_{CT}=e^{\beta_{CT}L_C/K_C}\)",
    "load_plus_tolerant_context" => raw"\(M_i=e^{\beta_iL_i/K_i},\quad f_{T0}^{\mathrm{co}}=\operatorname{logit}^{-1}[\operatorname{logit}(f_{T0})+\delta_f],\quad H_{CT}^{\mathrm{co}}=e^{\gamma_T}H_{CT}\)",
    "load_plus_tolerant_growth_context" => raw"\(\displaystyle M_C=e^{\beta_C\,L_C/K_C},\quad \frac{dT}{dt}=\rho_TG_C^{\mathrm{co}}\frac{T}{C}-M_CA_C(t)H_{CT}^{\mathrm{co}}(z)T,\quad \rho_T=e^{\gamma_{r,T}}\)",
    "fully_free_context_diagnostic" => raw"\(\displaystyle \frac{dX_i}{dt}=G_i^{\mathrm{co}}-E_i^{\mathrm{co}}(t,z)X_i,\quad \mathbf q_{\mathrm{drug}}^{\mathrm{co}}\ne\mathbf q_{\mathrm{drug}}^{\mathrm{mono}}\)",
)

function _report_parameter_count(row)
    for column in (:n_parameters, :parameter_count)
        column in propertynames(row) || continue
        value = row[column]
        value === missing && continue
        parsed = tryparse(Int, string(value))
        parsed === nothing || return parsed
    end
    return typemax(Int)
end

function _report_ranked_rows(df::DataFrame)
    valid = [_is_valid_metric(value) && abs(Float64(value)) < 1e11 for value in df.bic]
    ranked = copy(df[valid, :])
    sort!(ranked, :bic)
    ranked[!, :report_rank] = collect(1:nrow(ranked))
    ranked[!, :report_delta_bic] = Float64.(ranked.bic) .- minimum(Float64.(ranked.bic))
    return ranked
end

function _report_display_rows(df::DataFrame; limit::Int = 5)
    ranked = _report_ranked_rows(df)
    isempty(ranked) && return ranked
    shown_indices = collect(1:min(limit, nrow(ranked)))
    counts = [_report_parameter_count(row) for row in eachrow(ranked)]
    simplest_count = minimum(counts)
    simplest_index = findfirst(==(simplest_count), counts)
    simplest_index in shown_indices || push!(shown_indices, simplest_index)
    return ranked[shown_indices, :]
end

function _report_ranking_table(df::DataFrame; display_subset::Bool = true)
    ranked = display_subset ? _report_display_rows(df) : _report_ranked_rows(df)
    model_names = String.(ranked.model)
    ranks = Int.(ranked.report_rank)
    pooling = :pooling_mode in propertynames(ranked) ? String.(ranked.pooling_mode) : fill("", nrow(ranked))
    counts = [_report_parameter_count(row) for row in eachrow(ranked)]
    all_counts = [_report_parameter_count(row) for row in eachrow(_report_ranked_rows(df))]
    simplest_count = minimum(all_counts)
    roles = [rank == 1 && count == simplest_count ? "Winner; simplest" :
        (rank == 1 ? "Winner" : (count == simplest_count ? "Simplest candidate" : "Leading candidate"))
        for (rank, count) in zip(ranks, counts)]
    return DataFrame(
        ID = ["M$(rank)" for rank in ranks],
        Model = [get(REPORT_MODEL_LABELS, model, model) for model in model_names],
        Pooling = pooling,
        Equation = [get(REPORT_MODEL_EQUATIONS, model, "See the model source and CSV parameter artifact") for model in model_names],
        delta_BIC = round.(Float64.(ranked.report_delta_bic); digits = 3),
        Role = roles,
    )
end

function _report_bic_plot_html(df::DataFrame; title::AbstractString = "Model comparison")
    shown = _report_display_rows(df)
    isempty(shown) && return ""
    max_delta = maximum(Float64.(shown.report_delta_bic))
    rows = String[]
    for row in eachrow(shown)
        delta = Float64(row.report_delta_bic)
        width = max_delta > 0 ? max(1.5, 100 * delta / max_delta) : 1.5
        label = get(REPORT_MODEL_LABELS, String(row.model), String(row.model))
        push!(rows, """<div class="bic-row"><div class="bic-label"><strong>M$(row.report_rank)</strong><span>$(_html_escape(label))</span></div><div class="bic-track"><span class="bic-bar" style="width:$(round(width; digits=2))%"></span></div><output>$(round(delta; digits=2))</output></div>""")
    end
    return """<figure class="bic-figure"><h4>$(_html_escape(title)): Delta BIC</h4><div class="bic-axis" aria-label="$(_html_escape(title)) Delta BIC bar plot">$(join(rows))</div><figcaption>Horizontal bars show Delta BIC relative to M1 within this table. Shorter is better; M1 is zero by definition.</figcaption></figure>"""
end

function _report_top_five_html(path::AbstractString; by_cell_line::Bool = false, label::AbstractString = "Model comparison")
    isfile(path) || return "<p class=\"missing\">Ranking CSV not found: $(_html_escape(path))</p>"
    ranking = CSV.read(path, DataFrame)
    if by_cell_line && :cell_line in propertynames(ranking)
        blocks = String[]
        for group in groupby(ranking, :cell_line; sort = true)
            cell_line = String(first(group.cell_line))
            group_df = DataFrame(group)
            push!(blocks, "<h3>$(_html_escape(cell_line))</h3><div class=\"table-wrap\">" *
                _table_html(_report_ranking_table(group_df); limit = 6) * "</div>" *
                _report_bic_plot_html(group_df; title = "$(label), $(cell_line)"))
        end
        return join(blocks, "\n")
    end
    return "<div class=\"table-wrap\">" * _table_html(_report_ranking_table(ranking); limit = 6) * "</div>" *
        _report_bic_plot_html(ranking; title = label)
end

function _report_stage_figure(path::AbstractString, report_dir::AbstractString, alt::AbstractString, caption::AbstractString)
    isfile(path) || return "<p class=\"missing\">Graph not found: $(_html_escape(path))</p>"
    relative_path = replace(relpath(path, report_dir), "\\" => "/")
    return "<figure><img src=\"$relative_path\" alt=\"$(_html_escape(alt))\"><figcaption>$(_html_escape(caption))</figcaption></figure>"
end

function _percentile(values::Vector{Float64}, probability::Float64)
    isempty(values) && return NaN
    sorted = sort(values)
    index = clamp(ceil(Int, probability * length(sorted)), 1, length(sorted))
    return sorted[index]
end

function _treated_coculture_endpoint_bootstrap(csv_root::AbstractString, winner::AbstractString; n_bootstrap::Int = 5000, seed::Int = 4040)
    decoded_path = joinpath(csv_root, "coculture_treated", "coculture_treated_a2780_decoded.csv")
    overlay_path = joinpath(csv_root, "coculture_treated", "figures", "linked_treatment_combined_overlays.csv")
    isfile(decoded_path) && isfile(overlay_path) || return DataFrame()
    decoded = CSV.read(decoded_path, DataFrame)
    wells = filter(row -> occursin("_well_day_averages", String(row.source_file)), decoded)
    isempty(wells) && return DataFrame()
    endpoint_day = maximum(Float64.(wells.time))
    endpoint = wells[Float64.(wells.time) .== endpoint_day, :]
    overlay = CSV.read(overlay_path, DataFrame)
    selected = overlay[(String.(overlay.model) .== winner) .&
        (String.(overlay.context) .== "coculture") .&
        (Float64.(overlay.time) .== endpoint_day), :]
    rng = MersenneTwister(seed)
    rows = NamedTuple[]
    for group in groupby(endpoint, [:density, :mix, :cell_line]; sort = true)
        values = Float64.(group.count)
        samples = [mean(rand(rng, values, length(values))) for _ in 1:n_bootstrap]
        lineage = String(first(group.cell_line))
        component = lineage == "A2780Naive" ? "sensitive" : "resistant"
        prediction_rows = selected[(String.(selected.density) .== String(first(group.density))) .&
            (String.(selected.mix) .== String(first(group.mix))) .&
            (String.(selected.component) .== component), :]
        prediction = isempty(prediction_rows) ? NaN : Float64(first(prediction_rows.predicted))
        lower = _percentile(samples, 0.025)
        upper = _percentile(samples, 0.975)
        push!(rows, (
            density = String(first(group.density)),
            mix = String(first(group.mix)),
            lineage = lineage,
            endpoint_day = endpoint_day,
            n_wells = length(values),
            observed_mean = mean(values),
            ci95_lower = lower,
            ci95_upper = upper,
            model_prediction = prediction,
            model_within_observed_ci95 = isfinite(prediction) && lower <= prediction <= upper,
            bootstrap_method = "nonparametric well-resampling within condition",
            n_bootstrap = n_bootstrap,
            seed = seed,
        ))
    end
    result = DataFrame(rows)
    CSV.write(joinpath(csv_root, "coculture_treated", "linked_treatment_endpoint_bootstrap.csv"), result)
    return result
end

function _endpoint_bootstrap_html(endpoint::DataFrame, winner::AbstractString)
    isempty(endpoint) && return "<p class=\"missing\">Endpoint bootstrap could not be generated.</p>"
    shown = select(endpoint,
        :density => :Density,
        :mix => :Mix,
        :lineage => :Lineage,
        :n_wells => :Wells,
        :observed_mean => ByRow(x -> round(x; digits = 1)) => :Observed_mean,
        :ci95_lower => ByRow(x -> round(x; digits = 1)) => :CI95_lower,
        :ci95_upper => ByRow(x -> round(x; digits = 1)) => :CI95_upper,
        :model_prediction => ByRow(x -> round(x; digits = 1)) => :Model_prediction,
        :model_within_observed_ci95 => :Prediction_inside_CI)
    return """
<div class="uncertainty-audit">
<h3>Day-$(Int(first(endpoint.endpoint_day))) endpoint bootstrap: 95% confidence intervals</h3>
<p>For each density, mixture, and lineage, the six well-level endpoint measurements were resampled with replacement 5,000 times. The interval is the 2.5th to 97.5th percentile of the bootstrapped well mean. The separately exported aggregate row was excluded, preventing double-counting. The selected <code>$(_html_escape(winner))</code> prediction is shown for comparison.</p>
<div class="table-wrap">$(_table_html(shown; limit = nrow(shown)))</div>
<p class="audit-warning"><strong>Interpretation limit:</strong> these are confidence intervals for the observed endpoint mean, not confidence intervals for every fitted parameter and not proof that the latent sensitive/tolerant mechanism is identifiable.</p>
</div>
"""
end

function _linked_sensitivity_html(csv_root::AbstractString)
    path = joinpath(csv_root, "coculture_treated", "linked_treatment_identifiability.csv")
    isfile(path) || return ""
    sensitivity = CSV.read(path, DataFrame)
    shown = select(sensitivity,
        :parameter => :Parameter,
        :estimate => ByRow(x -> round(Float64(x); sigdigits = 5)) => :Estimate,
        :lower_bound => ByRow(x -> round(Float64(x); sigdigits = 5)) => :Lower_bound,
        :upper_bound => ByRow(x -> round(Float64(x); sigdigits = 5)) => :Upper_bound,
        :bound_position => ByRow(x -> round(Float64(x); digits = 3)) => :Bound_position,
        :identifiability => :Profile_status)
    return """
<div class="uncertainty-audit">
<h3>GrowthParameterEstimation two-sided sensitivity and bound profile</h3>
<p>The winning linked fit was passed through <code>profile_joint_fit_bounds_two_sided</code>. Each fitted dimension is challenged toward both bounds and classified by its position in the accepted interval. A value near 0 or 1 indicates that the optimum remains close to a bound and is practically weakly identified on that side.</p>
<div class="table-wrap compact-parameters">$(_table_html(shown; limit = nrow(shown)))</div>
</div>
"""
end

function _structured_parameter_summary(parameter_df::DataFrame, row)
    isempty(parameter_df) && return nothing
    mask = trues(nrow(parameter_df))
    if :model in propertynames(parameter_df) && :model in propertynames(row)
        mask .&= String.(parameter_df.model) .== String(row.model)
    elseif :timing_hypothesis in propertynames(parameter_df) && :model in propertynames(row)
        mask .&= String.(parameter_df.timing_hypothesis) .== String(row.model)
    end
    for column in (:cell_line, :pooling_mode)
        column in propertynames(parameter_df) && column in propertynames(row) || continue
        mask .&= String.(parameter_df[!, column]) .== String(row[column])
    end
    selected = parameter_df[mask, :]
    isempty(selected) && return nothing
    value_column = :effective_value in propertynames(selected) ? :effective_value :
        (:estimate in propertynames(selected) ? :estimate : nothing)
    value_column === nothing && return nothing
    entries = String[]
    for parameter_row in eachrow(selected)
        density = :density in propertynames(parameter_row) ? "$(parameter_row.density): " : ""
        value = round(Float64(parameter_row[value_column]); sigdigits = 7)
        push!(entries, "$(density)$(parameter_row.parameter)=$(value)")
    end
    return join(entries, "; ")
end

function _appendix_ranking_html(path::AbstractString; title::AbstractString, by_cell_line::Bool = false, parameter_path::Union{Nothing,AbstractString} = nothing)
    isfile(path) || return ""
    ranking = CSV.read(path, DataFrame)
    parameter_df = parameter_path !== nothing && isfile(parameter_path) ? CSV.read(parameter_path, DataFrame) : DataFrame()
    groups = by_cell_line && :cell_line in propertynames(ranking) ? groupby(ranking, :cell_line; sort = true) : [ranking]
    blocks = String[]
    for group in groups
        group_df = DataFrame(group)
        ranked = _report_ranked_rows(group_df)
        pooling = :pooling_mode in propertynames(ranked) ? String.(ranked.pooling_mode) : fill("", nrow(ranked))
        boundary = :boundary_issue in propertynames(ranked) ? string.(ranked.boundary_issue) : fill("", nrow(ranked))
        fallback_params = :params in propertynames(ranked) ? String.(ranked.params) : fill("Not exported", nrow(ranked))
        params = [something(_structured_parameter_summary(parameter_df, row), fallback)
            for (row, fallback) in zip(eachrow(ranked), fallback_params)]
        appendix = DataFrame(
            ID = ["M$(rank)" for rank in ranked.report_rank],
            Model = [get(REPORT_MODEL_LABELS, String(model), String(model)) for model in ranked.model],
            Code = String.(ranked.model),
            Pooling = pooling,
            BIC = round.(Float64.(ranked.bic); digits = 3),
            delta_BIC = round.(Float64.(ranked.report_delta_bic); digits = 3),
            Free_parameters = [_report_parameter_count(row) for row in eachrow(ranked)],
            Boundary_issue = boundary,
            Parameters = params,
        )
        suffix = by_cell_line ? ": $(String(first(group_df.cell_line)))" : ""
        push!(blocks, "<h3>$(_html_escape(title * suffix))</h3><div class=\"table-wrap appendix-table\">$(_table_html(appendix; limit = nrow(appendix)))</div>")
    end
    return join(blocks)
end

function _report_stage4_expanded_equations_html()
    return raw"""
<div class="notation-key equation-detail"><h3>Stage 4 expanded treated-coculture equations</h3>
<p>The BIC table below names the candidate-specific change. The full linked system inherits Stage 1 growth, Stage 2 treatment timing/Hill response, and Stage 3 coculture competition/death.</p>
<div class="math">\[C=S+T,\qquad L_N=N+\alpha_{NC,d}C,\qquad L_C=C+\alpha_{CN,d}N\]</div>
<div class="math">\[G_N^{\mathrm{co}}=G_N(N,L_N,t)-d_NN,\qquad G_C^{\mathrm{co}}=G_C(C,L_C,t)-d_CC\]</div>
<div class="math">\[M_N=e^{\beta_NL_N/K_{N,d}},\qquad M_C=e^{\beta_CL_C/K_C}\]</div>
<div class="math">\[f_{T0}^{\mathrm{co}}=\operatorname{logit}^{-1}\!\left(\operatorname{logit}(f_{T0})+\delta_f\right),\qquad H_{CT}^{\mathrm{co}}(z)=e^{\gamma_T}H_{CT}(z),\qquad \rho_T=e^{\gamma_{r,T}}\]</div>
<div class="math">\[\frac{dN}{dt}=G_N^{\mathrm{co}}-M_NA_N(t)H_N(z)N\]</div>
<div class="math">\[\frac{dS}{dt}=G_C^{\mathrm{co}}\frac{S}{C}-M_CA_C(t)H_{CS}(z)S\]</div>
<div class="math">\[\frac{dT}{dt}=\rho_TG_C^{\mathrm{co}}\frac{T}{C}-M_CA_C(t)H_{CT}^{\mathrm{co}}(z)T\]</div>
<p><strong>How to read the candidate rows:</strong> strict inheritance sets the extra coculture modifiers to their neutral values; load-scaled candidates estimate one or more \(\beta\) terms; tolerant-state candidates also estimate \(\delta_f\), \(\gamma_T\), and/or \(\gamma_{r,T}\).</p>
</div>
"""
end

function render_a2780_report_html(; start::AbstractString = pwd())
    root = IOUtils.package_root(start)
    report_dir = joinpath(root, "outputs", "reports")
    mkpath(report_dir)
    csv_root = joinpath(root, "outputs", "csv")
    image_root = joinpath(root, "outputs", "images")
    timing_ranking_path = joinpath(csv_root, "monoculture_treated", "monoculture_treated_timing_hypothesis_ranking.csv")
    timing_figure_path = joinpath(image_root, "monoculture_treated", "figures", "monoculture_treated_timing_hypothesis_grid.png")
    timing_winner = if isfile(timing_ranking_path)
        timing_ranking = sort!(CSV.read(timing_ranking_path, DataFrame), :bic)
        nrow(timing_ranking) > 0 ? String(first(timing_ranking.model)) : "independent_onset_gradual"
    else
        "independent_onset_gradual"
    end
    cis_activation_equation = timing_winner == "cis_onset_only" ?
        raw"A_C(t)=\mathbf{1}_{t>t_{\mathrm{on},C}}" :
        (timing_winner == "cis_gradual_only" ?
            raw"A_C(t)=1-e^{-\lambda_Ct}" :
            raw"A_C(t)=\mathbf{1}_{t>t_{\mathrm{on},C}}\left[1-e^{-\lambda_C(t-t_{\mathrm{on},C})}\right]")
    linked_status_path = joinpath(csv_root, "coculture_treated", "linked_treatment_status.csv")
    linked_winner = if isfile(linked_status_path)
        linked_status = CSV.read(linked_status_path, DataFrame)
        nrow(linked_status) > 0 ? String(first(linked_status.winning_model)) : "load_scaled"
    else
        "load_scaled"
    end
    linked_context_description = get(Dict(
        "strict_inheritance" => "no additional coculture treatment modifier",
        "competitor_scaled" => "competitor-load drug scaling",
        "load_scaled" => "total-load drug scaling",
        "load_plus_context_amplitude" => "total-load scaling and context amplitude",
        "tolerant_context_shift" => "initial tolerant-fraction and tolerant-kill shifts",
        "subpopulation_load_scaled" => "subpopulation-specific load scaling",
        "load_plus_tolerant_context" => "total-load scaling, an initial tolerant-fraction shift, and tolerant-kill attenuation",
        "load_plus_tolerant_growth_context" => "total-load scaling, tolerant-state shifts, and tolerant-growth plasticity",
        "fully_free_context_diagnostic" => "a fully separate diagnostic treatment vector",
    ), linked_winner, "the candidate-specific context terms listed in its ranking row")
    seed_audit_path = joinpath(csv_root, "coculture_treated", "linked_treatment_stage2_seed_audit.csv")
    stage4_inheritance_note = ""
    if isfile(seed_audit_path)
        seed_audit = CSV.read(seed_audit_path, DataFrame)
        cis_seed = seed_audit[String.(seed_audit.cell_line) .== "A2780cis", :]
        if nrow(cis_seed) > 0 && :nominal_stage2_winner in propertynames(cis_seed)
            row = first(cis_seed)
            if String(row.nominal_stage2_winner) != String(row.treatment_family)
                stage4_inheritance_note = " Stage 2 nominally selected $(row.nominal_stage2_winner), but Stage 4 uses the best population-balance-compatible $(row.treatment_family) candidate (Delta BIC=$(round(Float64(row.compatible_delta_bic); digits = 2))); the mechanisms are not separated by Delta BIC >= 2."
            end
        end
    end
    untreated_baseline_path = joinpath(csv_root, "monoculture_untreated", "untreated_group_baselines.csv")
    untreated_baselines = isfile(untreated_baseline_path) ? CSV.read(untreated_baseline_path, DataFrame) : DataFrame()
    function selected_growth_model(cell_line, fallback)
        isempty(untreated_baselines) && return fallback
        rows = untreated_baselines[String.(untreated_baselines.cell_line) .== cell_line, :]
        isempty(rows) ? fallback : String(first(rows.best_model))
    end
    naive_growth_model = selected_growth_model("A2780Naive", "logistic_growth")
    cis_growth_model = selected_growth_model("A2780cis", "theta_logistic_growth")
    naive_growth_label = get(REPORT_MODEL_LABELS, naive_growth_model, naive_growth_model)
    cis_growth_label = get(REPORT_MODEL_LABELS, cis_growth_model, cis_growth_model)
    naive_growth_equation = get(REPORT_MODEL_EQUATIONS, naive_growth_model, "")
    cis_growth_equation = get(REPORT_MODEL_EQUATIONS, cis_growth_model, "")

    stages = [
        (
            number = 1,
            title = "Untreated monoculture",
            note = "Top five model and pooling combinations per cell line. Both seeding densities are fitted jointly with fixed day-zero populations.",
            notation = raw"""
<div class="notation-key"><h3>Stage 1 notation key</h3><dl>
<dt>\(X_i(t)\)</dt><dd>Observed total population for cell line \(i\); \(i=N\) means A2780Naive and \(i=C\) means A2780cis.</dd>
<dt>\(\dot X_i=dX_i/dt\)</dt><dd>Rate of change of that total population, in measured cells per day.</dd>
<dt>\(G_i(X_i,t)\)</dt><dd>Candidate intrinsic growth law: logistic, theta-logistic, Gompertz, hard-lag theta-logistic, Baranyi-adjusted theta-logistic, or smooth-adaptation theta-logistic.</dd>
<dt>\(r_i,K_i,\theta_i\)</dt><dd>Intrinsic growth rate, carrying capacity, and optional theta-logistic shape.</dd>
<dt>\(\tau_i,q_{0,i},\lambda_{A,i}\)</dt><dd>Optional hard-lag duration, Baranyi initial physiological-state parameter, or smooth adaptation rate. Only the parameter belonging to the candidate being fitted is used.</dd>
</dl></div>
""",
            ranking = joinpath(csv_root, "monoculture_untreated", "monoculture_untreated_pooling_model_ranking.csv"),
            by_cell_line = true,
            figure = joinpath(image_root, "monoculture_untreated", "figures", "monoculture_untreated_pooling_model_grid.png"),
            caption = "Models shown: A2780Naive $(naive_growth_label) and A2780cis $(cis_growth_label), each using its selected density-pooling mode. These are joint 20k/30k untreated fits.",
        ),
        (
            number = 2,
            title = "Treated monoculture",
            note = "Top five whole-cell-line fits per cell line. Each model fits 20k and 30k across IC25, IC50, and IC75 while inheriting untreated growth.",
            notation = raw"""
<div class="notation-key"><h3>Stage 2 notation key</h3><dl>
<dt>\(X_i(t)\)</dt><dd>Total population in a one-state treatment model for cell line \(i\).</dd>
<dt>\(P_i(t),D_i(t)\)</dt><dd>Proliferating/live and damaged-but-still-visible compartments used only by the transit model.</dd>
<dt>\(S(t),T(t)\)</dt><dd>Drug-sensitive and drug-tolerant latent A2780cis subpopulations used only by the population-balance model; \(C=S+T\).</dd>
<dt>\(\widehat y_i(t)\)</dt><dd>Predicted observable: \(X_i\), \(P_i+\tfrac12D_i\), or \(S+T\), depending on the candidate.</dd>
<dt>\(G_i\)</dt><dd>The exact growth family and effective density-specific parameters inherited from Stage 1.</dd>
<dt>\(A_i(t)\)</dt><dd>Dimensionless time-activation term. It is 0 before onset and then turns drug action on either abruptly or gradually, depending on the timing candidate.</dd>
<dt>\(R(t;\lambda_i,t_{\mathrm{on},i})\)</dt><dd>Delayed ramp shorthand: \(\mathbf{1}_{t>t_{\mathrm{on},i}}[1-e^{-\lambda_i(t-t_{\mathrm{on},i})}]\).</dd>
<dt>\(H_i(z)\)</dt><dd>Hill dose-response term: \(H_i(z)=E_{\max,i}z^{h_i}/(EC_{50,i}^{h_i}+z^{h_i}+\varepsilon)\). It sets dose strength, not timing.</dd>
<dt>\(A_i(t)H_i(z)\)</dt><dd>The active per-capita treatment effect applied to a live/proliferating state. Transit models route this effect into a damaged-visible compartment before clearance.</dd>
</dl><p><strong>Why the derivatives differ:</strong> \(dX/dt\), \(dP/dt\), and \(dS/dt\) are rates for different biological state variables, not interchangeable names for the same quantity.</p></div>
""",
            ranking = joinpath(csv_root, "monoculture_treated", "monoculture_treated_joint_dose_model_ranking.csv"),
            by_cell_line = true,
            figure = joinpath(image_root, "monoculture_treated", "figures", "monoculture_treated_best_joint_model_by_environment.png"),
            caption = "Models shown: A2780Naive delayed Hill-ramp treatment and A2780cis sensitive/tolerant treatment, with the selected pooling mode for each lineage. One joint winner is used across both densities and all three doses for each cell line.",
        ),
        (
            number = 3,
            title = "Untreated coculture",
            note = "Top five coupled competition models across both densities and all three mix ratios. Lineage growth is inherited from untreated monoculture.",
            notation = raw"""
<div class="notation-key"><h3>Stage 3 notation key</h3><dl>
<dt>\(N(t),C(t)\)</dt><dd>A2780Naive and total A2780cis populations, respectively.</dd>
<dt>\(\dot N,\dot C\)</dt><dd>Lineage-specific rates of change. These are coupled because each lineage contributes to the other's crowding.</dd>
<dt>\(L_N=N+\alpha_{NC}C\)</dt><dd>Effective competitive load perceived by A2780Naive; \(L_N\) is a load, not another population compartment.</dd>
<dt>\(L_C=C+\alpha_{CN}N\)</dt><dd>Effective competitive load perceived by A2780cis.</dd>
<dt>\(G_N,G_C\)</dt><dd>Stage-1 lineage growth laws evaluated under the coculture competitive loads.</dd>
<dt>\(d_N,d_C\)</dt><dd>Optional first-order lineage loss rates in the competition-plus-death candidate.</dd>
</dl></div>
""",
            ranking = joinpath(csv_root, "coculture_untreated", "coculture_untreated_pooling_model_ranking.csv"),
            by_cell_line = false,
            figure = joinpath(image_root, "coculture_untreated", "figures", "coculture_untreated_best_mechanistic_fit_grid.png"),
            caption = "Model shown: asymmetric competition with lineage-specific loss using partial_5pct pooling, fitted jointly across 20k/30k and all three starting mixtures.",
        ),
        (
            number = 4,
            title = "Treated coculture",
            note = "Top five linked treatment hypotheses across all treated monoculture and coculture environments. Coculture treatment is IC50 at 1.0 uM (Hill effect signal z = 0.50). Stage-2 drug parameters are inherited within a +/-5% validation margin; every additional context effect is named explicitly.$(stage4_inheritance_note)",
            notation = raw"""
<div class="notation-key"><h3>Stage 4 notation key</h3><dl>
<dt>\(N(t),C(t)\)</dt><dd>Total A2780Naive and A2780cis populations. For the resistant population-balance model, \(C=S+T\).</dd>
<dt>\(X_i(t)\)</dt><dd>Compact lineage placeholder used in the BIC table: \(X_i=N\) for \(i=N\), and \(X_i=C\) for \(i=C\).</dd>
<dt>\(M_N,M_C\)</dt><dd>Dimensionless coculture modifiers multiplying the inherited Stage-2 drug effects.</dd>
<dt>\(A_N(t),A_C(t)\)</dt><dd>Cell-line-specific treatment timing rules selected in Stage 2. A2780Naive uses delayed ramp activation; A2780cis uses the BIC-selected timing rule from the Stage 2 audit.</dd>
<dt>\(H_N(z),H_{CS}(z),H_{CT}(z)\)</dt><dd>Inherited Stage-2 Hill dose-response terms for Naive, cis sensitive-like, and cis tolerant-like treatment effects. Their products with \(A_i(t)\) are the actual time- and dose-dependent kill/damage rates.</dd>
<dt>\(f_{T0}^{\mathrm{co}},H_{CT}^{\mathrm{co}}\)</dt><dd>Coculture-adjusted initial tolerant fraction and tolerant-state Hill kill term.</dd>
<dt>\(\rho_T\)</dt><dd>Multiplier on tolerant-state growth in treated coculture; \(\rho_T=1\) means no growth-context change.</dd>
</dl></div>
""",
            ranking = joinpath(csv_root, "coculture_treated", "linked_treatment_model_ranking.csv"),
            by_cell_line = false,
            figure = joinpath(image_root, "coculture_treated", "figures", "linked_treatment_coculture_grid.png"),
            caption = "Model shown: $(get(REPORT_MODEL_LABELS, linked_winner, linked_winner)) with linked_global pooling, fitted across both seeding densities and all mixture environments.",
        ),
    ]

    sections = String[
        "<header><a class=\"back-home\" href=\"../../../../index.html\">Back to reports home</a><h1>A2780 staged model comparison</h1><p>Mechanistic equations, Delta-BIC rankings, and explicitly labeled fitted graph grids for the four-stage analysis.</p><p class=\"artifact-note\">Each teaching table shows the five leading candidates plus the simplest tested candidate when it is not already present. Absolute BIC values and complete parameter vectors are retained in the appendix and <code>outputs/csv</code>.</p></header>",
    ]
    for stage in stages
        table = _report_top_five_html(stage.ranking; by_cell_line = stage.by_cell_line, label = "Stage $(stage.number) $(stage.title)")
        graph = _report_stage_figure(stage.figure, report_dir, stage.title, stage.caption)
        stage_details = stage.number == 4 ? _report_stage4_expanded_equations_html() : ""
        uncertainty = if stage.number == 4
            endpoint = _treated_coculture_endpoint_bootstrap(csv_root, linked_winner)
            _linked_sensitivity_html(csv_root) * _endpoint_bootstrap_html(endpoint, linked_winner)
        else
            ""
        end
        push!(sections, "<section><div class=\"stage-heading\"><span>Stage $(stage.number)</span><h2>$(_html_escape(stage.title))</h2></div><p>$(_html_escape(stage.note))</p>$(stage.notation)$(stage_details)$(table)$(graph)$(uncertainty)</section>")
        if stage.number == 2 && isfile(timing_ranking_path)
            timing_table = _report_top_five_html(timing_ranking_path; label = "Stage 2 timing audit")
            timing_graph = _report_stage_figure(
                timing_figure_path,
                report_dir,
                "Treated monoculture timing hypotheses",
                "The top three timing architectures are overlaid in every cell-line, density, and dose panel. The BIC-selected resistant timing is $(get(REPORT_MODEL_LABELS, timing_winner, timing_winner)).",
            )
            timing_notation = raw"""
<div class="notation-key"><h3>Timing-audit notation key</h3><dl>
<dt>\(R(t;\lambda_i,t_{\mathrm{on},i})\)</dt><dd>Delayed gradual activation: \(R(t;\lambda_i,t_{\mathrm{on},i})=\mathbf{1}_{t>t_{\mathrm{on},i}}\left[1-e^{-\lambda_i(t-t_{\mathrm{on},i})}\right]\). It equals 0 before onset and approaches 1 after onset.</dd>
<dt>\(\mathbf{1}_{t>t_{\mathrm{on},i}}\)</dt><dd>Pure onset switch. It is 0 before the fitted onset time and 1 after onset, with no gradual ramp.</dd>
<dt>\(\lambda_i\)</dt><dd>Activation speed in day\(^{-1}\). Larger \(\lambda_i\) means the delayed effect reaches full strength faster.</dd>
</dl></div>
"""
            push!(sections, "<section><div class=\"stage-heading\"><span>Stage 2 timing audit</span><h2>Onset and gradual-effect combinations</h2></div><p>All five timing hypotheses use one joint objective over the same 12 treated-monoculture trajectories. The A2780Naive delayed Hill-ramp and A2780cis sensitive/tolerant population model are otherwise unchanged. Boundary flags and multistart results remain in <code>outputs/csv/monoculture_treated</code>.</p>$(timing_notation)$(timing_table)$(timing_graph)</section>")
        end
    end
    push!(sections, """
<section class="equation-summary">
<div class="stage-heading"><span>Summary</span><h2>Equation progression</h2></div>
<p>The selected models build sequentially. Every growth, density, dose, delay, population-balance, competition, death, context-scaling, initial-condition, observation, and model-selection term used by the winning pipeline is defined below.</p>
<div class="progression">
  <article>
    <h3>Design, density pooling, and initial conditions</h3>
    <div class="math">\\[s_d=\\begin{cases}-1,&d=20\\mathrm{k}\\\\+1,&d=30\\mathrm{k}\\end{cases}\\]</div>
    <div class="math">\\[p_d=p_{\\mathrm{center}}e^{s_d\\delta_p},\\qquad |\\delta_p|\\le \\ln(1.05)\\]</div>
    <div class="math">\\[U_d=\\begin{cases}67,&d=20\\mathrm{k}\\\\100,&d=30\\mathrm{k}\\end{cases}\\]</div>
    <div class="math">\\[X_i(0)=U_d\\qquad\\text{(monoculture cell line }i\\text{)}\\]</div>
    <div class="math">\\[N(0)=f_{\\mathrm{mix}}U_d,\\qquad C(0)=(1-f_{\\mathrm{mix}})U_d\\qquad\\text{(coculture)}\\]</div>
    <div class="math">\\[S(0)=(1-f_{T0})C(0),\\qquad T(0)=f_{T0}C(0)\\]</div>
    <div class="math">\\[y_N(t)=N(t),\\qquad y_C(t)=C(t)=S(t)+T(t)\\]</div>
    <div class="math">\\[[x]_+=\\max(x,0)\\]</div>
    <p>The density contrast is symmetric and bounded at plus-or-minus five percent. Untreated partial pooling applies it to <code>r</code> and <code>K</code>; treated partial pooling applies it to treatment amplitudes; untreated coculture partial pooling applies it to both competition coefficients together. Starting populations are fixed and are never kinetic parameters.</p>
  </article>
  <article>
    <h3>1. Intrinsic untreated monoculture growth</h3>
    <p><strong>A2780Naive: $(_html_escape(naive_growth_label)).</strong></p>
    <div class="math">$(naive_growth_equation)</div>
    <p><strong>A2780cis: $(_html_escape(cis_growth_label)).</strong></p>
    <div class="math">$(cis_growth_equation)</div>
    <p>The displayed equations are read from the selected Stage 1 baseline artifact when this report is built. The exact winning family, including any fitted lag or adaptation term, and its density-specific effective parameters are inherited by all later stages. The Baranyi candidate is used only as a phenomenological lag description here; its original biological interpretation was developed for microbial growth and is not evidence of an ovarian-cancer mechanism.</p>
  </article>
  <article>
    <h3>2. Add treated-monoculture effects</h3>
    <div class="math">\\[z=\\begin{cases}0.25,&D=0.67\\,\\mu\\mathrm{M}\\;(\\mathrm{IC}_{25})\\\\0.50,&D=1.00\\,\\mu\\mathrm{M}\\;(\\mathrm{IC}_{50})\\\\0.75,&D=1.47\\,\\mu\\mathrm{M}\\;(\\mathrm{IC}_{75})\\end{cases}\\]</div>
    <div class="math">\\[A_N(t)=\\mathbf{1}_{t>t_{\\mathrm{on},N}}\\left[1-e^{-\\lambda_N(t-t_{\\mathrm{on},N})}\\right]\\]</div>
    <div class="math">\\[$(cis_activation_equation)\\]</div>
    <div class="math">\\[\\varepsilon=10^{-12}\\]</div>
    <div class="math">\\[H_{N,d}(z)=\\frac{E_{\\max,N,d}z^{h_N}}{EC_{50,N}^{h_N}+z^{h_N}+\\varepsilon}\\]</div>
    <div class="math">\\[H_{CS,d}(z)=\\frac{E_{\\max,CS,d}z^4}{0.5^4+z^4+\\varepsilon},\\qquad H_{CT,d}(z)=\\frac{E_{\\max,CT,d}z^4}{0.5^4+z^4+\\varepsilon}\\]</div>
    <div class="math">\\[E_{\\max,i,d}=E_{\\max,i,\\mathrm{center}}e^{s_d\\delta_{E,i}},\\qquad |\\delta_{E,i}|\\le\\ln(1.05)\\]</div>
    <div class="math">\\[\\frac{dN}{dt}=G_N(N,t)-A_N(t)H_{N,d}(z)N\\]</div>
    <div class="math">\\[C=S+T\\]</div>
    <div class="math">\\[\\frac{dS}{dt}=G_C(C,t)\\frac{S}{C}-A_C(t)H_{CS,d}(z)S\\]</div>
    <div class="math">\\[\\frac{dT}{dt}=G_C(C,t)\\frac{T}{C}-A_C(t)H_{CT,d}(z)T\\]</div>
    <p>The A2780Naive winner is delayed Hill-ramp kill. The A2780cis winner has separate sensitive and tolerant kill amplitudes and an estimated initial tolerant fraction. A five-model timing audit selected <code>$(timing_winner)</code> by BIC; onset, gradual activation, shared-onset, and bounded-onset alternatives remain in the ranking rather than being silently discarded.</p>
  </article>
  <article>
    <h3>3. Add untreated-coculture competition and death</h3>
    <div class="math">\\[C=S+T\\]</div>
    <div class="math">\\[\\alpha_{NC,d}=\\alpha_{NC,\\mathrm{center}}e^{s_d\\delta_\\alpha},\\qquad \\alpha_{CN,d}=\\alpha_{CN,\\mathrm{center}}e^{s_d\\delta_\\alpha}\\]</div>
    <div class="math">\\[L_N=N+\\alpha_{NC,d}C,\\qquad L_C=C+\\alpha_{CN,d}N\\]</div>
    <div class="math">\\[G_N^{\\mathrm{co}}=G_N(N,L_N,t)-d_NN\\]</div>
    <div class="math">\\[G_C^{\\mathrm{co}}=G_C(C,L_C,t)-d_CC\\]</div>
    <div class="math">\\[\\frac{dN}{dt}=G_N^{\\mathrm{co}},\\qquad \\frac{dC}{dt}=G_C^{\\mathrm{co}}\\qquad\\text{(untreated coculture)}\\]</div>
    <p>The winner is asymmetric competition with lineage-specific death. The two competition coefficients may differ by direction; their common density contrast is bounded at plus-or-minus five percent. Death rates are shared across density. Every selected monoculture growth parameter, including a lag/adaptation parameter when present, remains fixed.</p>
  </article>
  <article>
    <h3>4. Fully expanded treated-coculture winner</h3>
    <div class="math">\\[D=1.0\\,\\mu\\mathrm{M}=\\mathrm{IC}_{50},\\qquad z=0.50\\]</div>
    <div class="math">\\[M_N=\\exp\\!\\left(\\operatorname{clamp}\\!\\left[\\beta_N\\frac{L_N}{K_{N,d}},-4,4\\right]\\right),\\qquad M_C=\\exp\\!\\left(\\operatorname{clamp}\\!\\left[\\beta_C\\frac{L_C}{K_C},-4,4\\right]\\right)\\]</div>
    <div class="math">\\[f_{T0}^{\\mathrm{co}}=\\operatorname{logit}^{-1}\\!\\left[\\operatorname{logit}(f_{T0})+\\delta_f\\right]\\]</div>
    <div class="math">\\[H_{CT,d}^{\\mathrm{co}}(0.50)=e^{\\gamma_T}H_{CT,d}(0.50),\\qquad \\rho_T=e^{\\gamma_{r,T}}\\]</div>
    <div class="math">\\[\\frac{dN}{dt}=G_N^{\\mathrm{co}}-M_NA_N(t)H_{N,d}(0.50)N\\]</div>
    <div class="math">\\[\\frac{dS}{dt}=G_C^{\\mathrm{co}}\\frac{S}{C}-M_CA_C(t)H_{CS,d}(0.50)S\\]</div>
    <div class="math">\\[\\frac{dT}{dt}=\\rho_TG_C^{\\mathrm{co}}\\frac{T}{C}-M_CA_C(t)H_{CT,d}^{\\mathrm{co}}(0.50)T\\]</div>
    <div class="math">\\[C=S+T,\\qquad L_N=N+\\alpha_{NC,d}C,\\qquad L_C=C+\\alpha_{CN,d}N\\]</div>
    <p>This is the selected <code>$(linked_winner)</code> model. Intrinsic growth and untreated competition/death parameters are fixed from Stages 1 and 3. Stage 2 drug parameters and the BIC-selected timing architecture are inherited within a plus-or-minus five-percent validation window. Its fitted context structure is $(_html_escape(linked_context_description)); modifiers not named by this winner are fixed at their neutral values.</p>
  </article>
  <article>
    <h3>Fitting, scale normalization, and BIC selection</h3>
    <div class="math">\\[s_j=\\max_t y_j(t)\\]</div>
    <div class="math">\\[\\mathrm{SSE}_{\\mathrm{scaled}}=\\sum_j\\sum_t\\left[\\frac{y_j(t)-\\widehat y_j(t)}{s_j}\\right]^2\\]</div>
    <div class="math">\\[\\mathrm{SSE}_{\\mathrm{raw}}=\\sum_j\\sum_t\\left[y_j(t)-\\widehat y_j(t)\\right]^2\\]</div>
    <div class="math">\\[\\mathrm{BIC}=n\\ln\\!\\left(\\frac{\\mathrm{SSE}_{\\mathrm{scaled}}}{n}\\right)+k\\ln n\\]</div>
    <p><strong>Small-to-large trajectory normalization.</strong> Each trajectory \\(j\\) is divided by its own observed peak \\(s_j\\). A 500-cell trajectory and a 4,000-cell trajectory therefore contribute comparable relative errors instead of the larger curve dominating merely because its residuals have larger units. The scaled SSE is dimensionless and drives optimization and BIC; raw SSE in squared cell-count units is also exported for scale-aware interpretation.</p>
    <p><strong>What is jointly fitted.</strong> Every candidate is fitted across all trajectories stated in its stage with one objective. In the canonical report, image tiles are averaged within wells and then biological samples are averaged at each time point. The sample-aware report preserves sample-level trajectories after within-well averaging and displays their between-sample uncertainty bands.</p>
    <p><strong>BIC counting.</strong> Here \\(n\\) is the total number of fitted time-point observations across all joint trajectories, not the number of wells, tiles, or environments. The count \\(k\\) includes every freely optimized center, shape, treatment, context, and pooling-contrast parameter. Fixed inherited parameters and fixed day-zero populations are not counted again. Lower BIC is better only relative to candidates fitted to the same observations and objective; it is not proof that the winning biological mechanism is true.</p>
    <p><strong>How to read Delta BIC.</strong> For candidate <code>m</code>, <code>Delta BIC(m) = BIC(m) - minimum BIC</code>, so the table winner is always zero. As descriptive evidence bands, 0-2 means little separation from the winner, 2-6 indicates positive separation, 6-10 strong separation, and values above 10 very strong separation. These are relative model-selection heuristics, not probabilities, confidence intervals, or proof of mechanism. Comparisons are valid only within a table whose candidates use the same observations and residual objective.</p>
    <p><strong>Initial conditions and time origin.</strong> Day zero is fixed at 67 measured cells for 20k seeding and 100 for 30k seeding. The first observed point is day 1. These anchors are passed through the package's fixed-initial-time and initial-state builder interface and are not estimated kinetic parameters.</p>
    <p><strong>Pooling.</strong> Shared models use one parameter value across 20k and 30k. Partial pooling permits only \\(r\\) and \\(K\\), or the stage-specific named effects, to differ through symmetric log contrasts bounded at plus-or-minus five percent: \\(p_{20k}=p_c e^{-\\delta_p}\\), \\(p_{30k}=p_c e^{\\delta_p}\\), \\(|\\delta_p|\\le\\ln(1.05)\\). Fully independent density fits are retained as diagnostics but cannot silently become inherited defaults.</p>
    <p><strong>Fit validity.</strong> A fit is rejected if any parameter, prediction, SSE, or BIC is non-finite, or if its objective retains the failure sentinel. Boundary profiles expand requested bounds and accept an expansion only when the improvement is scientifically material under the configured BIC and margin rules. Multistart results and parameter-stability summaries are retained where that stage uses multistart optimization.</p>
    <p><strong>Uncertainty and weighting.</strong> Error ribbons show empirical variation; they are not inverse-variance weights in the canonical scaled-SSE objective. Later-stage fits inherit only eligible finite winners and record the exact source family and parameters in inheritance audit CSVs. This prevents a downstream stage from quietly replacing a selected lag, growth, dose-response, or competition law.</p>
  </article>
</div>
<div class="model-guide">
  <h2>Symbol and function guide</h2>
  <p>This glossary states what every symbol does biologically and mathematically. A symbol with subscript <code>d</code> can differ between 20k and 30k only through the explicitly bounded density-pooling rule.</p>
  <div class="guide-grid">
    <article>
      <h3>Populations and indices</h3>
      <dl>
        <dt>\\(t\\)</dt><dd>Time in days from the fixed day-zero initial condition.</dd>
        <dt>\\(X_i(t)\\)</dt><dd>Generic total population for a one-state monoculture model of cell line \\(i\\). This symbol avoids assigning a biological compartment that the model does not contain.</dd>
        <dt>\\(N(t)\\)</dt><dd>A2780Naive population when both lineages are modeled explicitly.</dd>
        <dt>\\(P_i(t),D_i(t)\\)</dt><dd>Proliferating/live and damaged-visible states used only in the transit-compartment model.</dd>
        <dt>\\(S(t),T(t)\\)</dt><dd>Latent drug-sensitive and drug-tolerant A2780cis subpopulations used only in the population-balance model.</dd>
        <dt>\\(C(t)=S(t)+T(t)\\)</dt><dd>Total observed A2780cis population. The experiment does not separately observe \\(S\\) and \\(T\\).</dd>
        <dt>\\(\\widehat y_i(t)\\)</dt><dd>Model prediction on the observed scale. It equals \\(X_i\\), \\(P_i+\\tfrac12D_i\\), or \\(S+T\\), according to model structure.</dd>
        <dt>\\(\\dot X=dX/dt\\)</dt><dd>Derivative notation for a state's rate of change. Different letters mean distinct biological states, not inconsistent names for one quantity.</dd>
        <dt>\\(i\\)</dt><dd>Lineage index: \\(N\\) for A2780Naive or \\(C\\) for A2780cis.</dd>
        <dt>\\(d\\)</dt><dd>Starting-density environment: 20k or 30k.</dd>
      </dl>
    </article>
    <article>
      <h3>Intrinsic growth</h3>
      <dl>
        <dt><code>r_i</code></dt><dd>Intrinsic per-capita growth rate, with units day^-1.</dd>
        <dt><code>K_i</code></dt><dd>Carrying capacity on the measured cell-count scale.</dd>
        <dt><code>theta_i</code></dt><dd>Dimensionless theta-logistic shape. \\(\\theta=1\\) recovers ordinary logistic growth; other values shift how sharply crowding suppresses growth.</dd>
        <dt><code>tau_i</code></dt><dd>Hard-lag duration in days. The hard-lag candidate sets intrinsic growth to zero until \\(t>\\tau_i\\), then switches the fitted theta-logistic law on.</dd>
        <dt><code>q0_i</code></dt><dd>Positive Baranyi initial-state parameter. It defines \\(\\alpha_B(t)=q_0/(q_0+e^{-r_it})\\), which smoothly raises the fraction of the intrinsic growth rate expressed over time.</dd>
        <dt><code>lambda_A,i</code></dt><dd>Smooth adaptation rate in day\\(^{-1}\\). It defines \\(A_i(t)=1-e^{-\\lambda_{A,i}t}\\), equivalent to \\(dA_i/dt=\\lambda_{A,i}(1-A_i)\\) with \\(A_i(0)=0\\).</dd>
        <dt><code>[x]_+</code></dt><dd>Positive part \\([x]_+=\\max(x,0)\\), used in monoculture fits to prevent the inherited growth term from becoming artificial negative crowding death above <code>K</code>.</dd>
      </dl>
    </article>
    <article>
      <h3>Starting density and mixture</h3>
      <dl>
        <dt><code>U_d</code></dt><dd>Fixed total at day zero: 67 for 20k and 100 for 30k.</dd>
        <dt><code>f_mix</code></dt><dd>Nominal A2780Naive fraction from the 25:75, 50:50, or 75:25 mixture label.</dd>
        <dt><code>f_T0</code></dt><dd>Estimated initial tolerant fraction within the A2780cis population.</dd>
        <dt><code>s_d</code></dt><dd>Density sign, -1 for 20k and +1 for 30k, used to create symmetric deviations around a cell-line center.</dd>
        <dt><code>delta_p</code></dt><dd>Log-scale density contrast bounded by \\(|\\delta_p|\\le\\ln(1.05)\\), so an eligible effective parameter can deviate by at most five percent in either direction.</dd>
      </dl>
    </article>
    <article>
      <h3>Drug dose and Hill response</h3>
      <dl>
        <dt><code>D</code></dt><dd>Physical drug concentration in uM.</dd>
        <dt><code>z</code></dt><dd>Normalized IC effect level: 0.25, 0.50, or 0.75 for IC25, IC50, or IC75.</dd>
        <dt><code>Emax_i</code></dt><dd>Maximum per-capita drug-effect amplitude, with units day^-1.</dd>
        <dt><code>EC50_i</code></dt><dd>Effect level or concentration at which the Hill term reaches half of <code>Emax</code>.</dd>
        <dt><code>h_i</code></dt><dd>Dimensionless Hill coefficient controlling dose-response steepness.</dd>
        <dt><code>H_i(z)</code></dt><dd>Hill dose-response function. It converts dose/effect level into a potential per-capita kill or damage rate; it does not encode treatment delay.</dd>
      </dl>
    </article>
    <article>
      <h3>Treatment activation A(t)</h3>
      <div class="math">\\[A_N(t)=\\mathbf{1}_{t>t_{\\mathrm{on},N}}\\left[1-e^{-\\lambda_N(t-t_{\\mathrm{on},N})}\\right]\\]</div>
      <div class="math">\\[$(cis_activation_equation)\\]</div>
      <dl>
        <dt><code>A_i(t)</code></dt><dd>Dimensionless fraction of the eventual treatment effect active at time <code>t</code>. The Naive lineage uses a delayed ramp. The resistant lineage uses the BIC-selected <code>$(timing_winner)</code> timing rule.</dd>
        <dt><code>t_onset,i</code></dt><dd>Delay in days before drug action starts.</dd>
        <dt><code>lambda_i</code></dt><dd>Activation rate in day^-1. Larger values switch treatment on faster. After onset, the activation half-time is \\(\\ln(2)/\\lambda_i\\).</dd>
        <dt>\\(k_i(t,z)=A_i(t)H_i(z)\\)</dt><dd>Actual time- and dose-dependent per-capita treatment rate applied to the population.</dd>
      </dl>
    </article>
    <article>
      <h3>Coculture interaction</h3>
      <dl>
        <dt><code>alpha_NC</code></dt><dd>Effect of one cis cell on the crowding experienced by A2780Naive.</dd>
        <dt><code>alpha_CN</code></dt><dd>Effect of one A2780Naive cell on the crowding experienced by A2780cis.</dd>
        <dt>\\(L_N=N+\\alpha_{NC}C\\)</dt><dd>Effective competitive load perceived by A2780Naive.</dd>
        <dt>\\(L_C=C+\\alpha_{CN}N\\)</dt><dd>Effective competitive load perceived by A2780cis.</dd>
        <dt><code>d_N,d_C</code></dt><dd>Additional lineage-specific first-order loss rates in day^-1 in the competition-with-death model.</dd>
      </dl>
    </article>
    <article>
      <h3>Coculture treatment modifier</h3>
      <dl>
        <dt><code>beta_N,beta_C</code></dt><dd>Dimensionless coefficients linking normalized competitive load to drug-effect scaling.</dd>
        <dt><code>M_i</code></dt><dd>Dimensionless multiplier on inherited drug action. \\(M_i=1\\) means no coculture modification, above one strengthens drug action, and below one weakens it.</dd>
        <dt><code>delta_f</code></dt><dd>Log-odds shift in the latent tolerant fraction specifically in treated coculture: \\(f_{T0}^{\\mathrm{co}}=\\operatorname{logit}^{-1}[\\operatorname{logit}(f_{T0})+\\delta_f]\\).</dd>
        <dt><code>gamma_T</code></dt><dd>Log multiplier on tolerant-compartment drug kill in coculture. Its exponential directly scales <code>H_CT</code>.</dd>
        <dt><code>rho_T</code></dt><dd>Positive multiplier on tolerant-compartment growth in treated coculture. It equals <code>exp(gamma_r,T)</code>; the selected estimate is approximately 1.25.</dd>
        <dt><code>clamp(x,-4,4)</code></dt><dd>Numerical/biological guard limiting the log multiplier before exponentiation.</dd>
        <dt>Current sign</dt><dd>The fitted <code>beta_N</code> and <code>beta_C</code> are negative, so larger competitive load attenuates effective drug action in the selected model.</dd>
      </dl>
    </article>
    <article>
      <h3>Fit and selection quantities</h3>
      <dl>
        <dt><code>y_j(t)</code></dt><dd>Observed value for trajectory <code>j</code>.</dd>
        <dt><code>yhat_j(t)</code></dt><dd>ODE prediction for the same observation.</dd>
        <dt><code>scale_j</code></dt><dd>Peak observed value of trajectory <code>j</code>, used so high-count trajectories do not dominate optimization.</dd>
        <dt><code>n</code></dt><dd>Total observations in the joint objective.</dd>
        <dt><code>k</code></dt><dd>Number of parameters estimated in that candidate, including pooling contrasts and context modifiers.</dd>
        <dt><code>BIC</code></dt><dd>Bayesian information criterion; lower is preferred because it rewards fit while penalizing additional fitted parameters.</dd>
      </dl>
    </article>
  </div>

  <h2>How every candidate model works</h2>
  <article class="model-family">
    <h3>Stage 1: untreated monoculture candidates</h3>
    <p class="equation-label"><strong>Logistic</strong></p><div class="math">\\[\\frac{dX}{dt}=rX\\left[1-\\frac{X}{K}\\right]_+\\]</div>
    <p>Growth is nearly exponential at low population and slows linearly with occupancy \\(X/K\\). It has the fewest mechanistic parameters.</p>
    <p class="equation-label"><strong>Theta-logistic</strong></p><div class="math">\\[\\frac{dX}{dt}=rX\\left[1-\\left(\\frac{X}{K}\\right)^\\theta\\right]_+\\]</div>
    <p>Adds a shape parameter. This can move most growth suppression closer to or farther from carrying capacity without changing the low-density growth rate.</p>
    <p class="equation-label"><strong>Gompertz</strong></p><div class="math">\\[\\frac{dX}{dt}=rX\\ln\\!\\left(\\frac{K}{X}\\right)\\]</div>
    <p>Per-capita growth decreases logarithmically with population, producing asymmetric sigmoidal trajectories.</p>
    <p class="equation-label"><strong>Simple-death diagnostic</strong></p><div class="math">\\[\\frac{dX}{dt}=rX\\left[1-\\frac{X}{K}\\right]_+-dX\\]</div>
    <p>Tests whether an untreated first-order loss term is needed. It is diagnostic and is inherited only with strong identifiable BIC support.</p>
    <p class="equation-label"><strong>Allee diagnostic</strong></p><div class="math">\\[\\frac{dX}{dt}=rX\\left[1-\\frac{X}{K}\\right]_+\\left(\\frac{X}{A_{\\mathrm{Allee}}}-1\\right)\\]</div>
    <p>Tests for suppressed or negative growth below a critical population threshold. It is not treated as a default baseline.</p>
    <p><strong>Pooling:</strong> shared and partial-pooling versions use the same biological ODE. Partial pooling adds bounded density contrasts to <code>r</code> and <code>K</code>; independent fits are diagnostic only.</p>
  </article>
  <article class="model-family">
    <h3>Stage 2: treated monoculture candidates</h3>
    <p class="equation-label"><strong>Anchored linear kill</strong></p><div class="math">\\[\\frac{dX}{dt}=G(X)-k_{\\mathrm{kill}}DX\\]</div>
    <p>Immediate kill proportional to physical concentration. It is simple but cannot create a delayed rise-then-decline trajectory.</p>
    <p class="equation-label"><strong>Anchored concentration-Hill kill</strong></p><div class="math">\\[\\frac{dX}{dt}=G(X)-\\frac{E_{\\max}D^h}{EC_{50}^h+D^h}X\\]</div>
    <p>Immediate saturating concentration response. It remains monotonic in time at fixed dose.</p>
    <p class="equation-label"><strong>Time-decay kill</strong></p><div class="math">\\[\\frac{dX}{dt}=G(X)-k_{\\mathrm{kill}}De^{-\\lambda t}X\\]</div>
    <p>Drug action is strongest initially and wanes over time, representing clearance, adaptation, or loss of exposure and allowing later regrowth.</p>
    <p class="equation-label"><strong>Immediate IC-Hill</strong></p><div class="math">\\[\\frac{dX}{dt}=G(X)-H(z)X\\]</div>
    <p>Uses normalized IC effect levels rather than physical concentration, with no onset delay.</p>
    <p class="equation-label"><strong>Hill ramp</strong></p><div class="math">\\[\\frac{dX}{dt}=G(X)-\\left(1-e^{-\\lambda t}\\right)H(z)X\\]</div>
    <p>Drug action accumulates from day zero.</p>
    <p class="equation-label"><strong>Delayed Hill ramp</strong></p><div class="math">\\[\\frac{dX}{dt}=G(X)-A(t)H(z)X\\]</div>
    <p>Adds an explicit onset time before the accumulating effect. This is the selected A2780Naive mechanism.</p>
    <p class="equation-label"><strong>Transit damage</strong></p><div class="math">\\[\\frac{dP}{dt}=G(P)-A(t)H(z)P,\\qquad \\frac{dD}{dt}=A(t)H(z)P-k_{\\mathrm{clear}}D\\]</div>
    <div class="math">\\[\\widehat y=P+\\frac12D\\]</div>
    <p>Cells leave the proliferating live compartment, remain partly visible while damaged, and are cleared later. The candidate fixes Hill <code>EC50=0.5</code> and <code>h=4</code>.</p>
    <p class="equation-label"><strong>Sensitive/tolerant populations</strong></p><div class="math">\\[\\frac{dS}{dt}=G(S+T)\\frac{S}{S+T}-A_C(t)H_{CS}(z)S,\\qquad \\frac{dT}{dt}=G(S+T)\\frac{T}{S+T}-A_C(t)H_{CT}(z)T\\]</div>
    <div class="math">\\[\\widehat y=S+T\\]</div>
    <p>Both latent states share one inherited total-growth law but have different kill amplitudes. Their initial split is controlled by <code>f_T0</code>. This is the selected A2780cis mechanism.</p>
    <p class="equation-label"><strong>Intracellular platinum PK/PD</strong></p><div class="math">\\[\\frac{dc_i}{dt}=D-k_{\\mathrm{efflux}}c_i,\\qquad \\frac{dc_k}{dt}=c_i-k_{\\mathrm{repair}}c_k,\\qquad \\frac{dX}{dt}=G(X)-H(c_k)X\\]</div>
    <p>This literature-grounded diagnostic candidate adds lumped intracellular platinum and DNA-bound platinum states. Uptake and binding scales are fixed because cell-count-only trajectories cannot separately identify them from kill amplitude; intracellular platinum or DNA-adduct measurements would be needed for a fully mechanistic fit.</p>
    <p><strong>Timing audit:</strong> independent onset plus ramp; shared onset plus ramp; onset differences bounded within 0.5 day; resistant gradual-only; and resistant onset-only.</p>
    <p>These timing candidates refit the same 12-trajectory objective through <code>GrowthParameterEstimation.run_joint_multistart</code>. BIC therefore compares the timing architectures with their different parameter counts on an identical data set.</p>
  </article>
  <article class="model-family">
    <h3>Stage 3: untreated coculture candidates</h3>
    <p class="equation-label"><strong>Symmetric competition</strong></p><div class="math">\\[\\alpha_{NC}=\\alpha_{CN}=\\alpha,\\qquad \\frac{dN}{dt}=G_N(N+\\alpha C),\\qquad \\frac{dC}{dt}=G_C(C+\\alpha N)\\]</div>
    <p>Both lineages exert the same cross-crowding strength.</p>
    <p class="equation-label"><strong>Asymmetric competition</strong></p><div class="math">\\[\\alpha_{NC}\\ne\\alpha_{CN},\\qquad \\frac{dN}{dt}=G_N(N+\\alpha_{NC}C),\\qquad \\frac{dC}{dt}=G_C(C+\\alpha_{CN}N)\\]</div>
    <p>Each lineage can affect the other differently.</p>
    <p class="equation-label"><strong>Asymmetric competition plus death</strong></p><div class="math">\\[\\frac{dN}{dt}=G_N(L_N)-d_NN,\\qquad \\frac{dC}{dt}=G_C(L_C)-d_CC\\]</div>
    <p>Adds lineage-specific loss beyond competitive growth suppression. This is the selected Stage 3 model. All six density/mix environments are fitted together.</p>
  </article>
  <article class="model-family">
    <h3>Stage 4: linked treated-coculture hypotheses</h3>
    <p class="equation-label"><strong>Strict inheritance</strong></p><div class="math">\\[M_N=M_C=1\\]</div>
    <p>The same intrinsic drug vector acts in monoculture and coculture with no context modification.</p>
    <p class="equation-label"><strong>Competitor-scaled</strong></p><div class="math">\\[M_N=e^{\\beta_N\\alpha_{NC}C/K_N},\\qquad M_C=e^{\\beta_C\\alpha_{CN}N/K_C}\\]</div>
    <p>Only the opposing lineage's contribution modifies drug action.</p>
    <p class="equation-label"><strong>Load-scaled</strong></p><div class="math">\\[M_N=e^{\\beta_NL_N/K_N},\\qquad M_C=e^{\\beta_CL_C/K_C}\\]</div>
    <p>Total effective load, including self and competitor, modifies drug action. This is the BIC winner.</p>
    <p class="equation-label"><strong>Load plus context amplitude</strong></p><div class="math">\\[M_i=e^{\\beta_iL_i/K_i+\\gamma_i}\\]</div>
    <p>Adds a density-independent coculture shift <code>gamma_i</code> on top of load scaling.</p>
    <p class="equation-label"><strong>Tolerant context shift</strong></p><div class="math">\\[f_{T0}^{\\mathrm{co}}=\\operatorname{logit}^{-1}[\\operatorname{logit}(f_{T0})+\\delta_f],\\qquad H_{CT}^{\\mathrm{co}}=e^{\\gamma_T}H_{CT}\\]</div>
    <p>Allows coculture to enrich the latent tolerant state and attenuate its kill rate while preserving the stage-2 sensitive/tolerant model.</p>
    <p class="equation-label"><strong>Subpopulation load scaling</strong></p><div class="math">\\[M_{CS}=e^{\\beta_{CS}L_C/K_C},\\qquad M_{CT}=e^{\\beta_{CT}L_C/K_C}\\]</div>
    <p>Tests whether competitive load scales sensitive and tolerant kill differently without changing their initial composition.</p>
    <p class="equation-label"><strong>Load plus tolerant context</strong></p><div class="math">\\[M_N,M_C,\\delta_f,\\gamma_T\\quad\\text{jointly define load scaling and tolerant-state modification}\\]</div>
    <p>Combines total-load drug scaling with the resistant-state shift.</p>
    <p class="equation-label"><strong>Load plus tolerant growth context</strong></p><div class="math">\\[M_C=e^{\\beta_C\\,L_C/K_C},\\qquad \\left.\\frac{dT}{dt}\\right|_{\\mathrm{growth}}=\\rho_TG_C^{\\mathrm{co}}\\frac{T}{C},\\qquad \\rho_T=e^{\\gamma_{r,T}}\\]</div>
    <p>Adds one bounded phenotype-specific growth-plasticity multiplier to the preceding model. This is the selected Stage 4 mechanism; BIC pays for the additional fitted parameter.</p>
    <p class="equation-label"><strong>Fully free context diagnostic</strong></p><div class="math">\\[\\mathbf q_{\\mathrm{drug}}^{\\mathrm{co}}\\ne\\mathbf q_{\\mathrm{drug}}^{\\mathrm{mono}}\\]</div>
    <p>Fits separate complete drug vectors by context. It tests whether inheritance fails, but is deliberately diagnostic because it gives up mechanistic sharing.</p>
    <p><strong>Important:</strong> Stage 4 does not attach an arbitrary coefficient to every term. Growth and untreated interaction parameters are fixed from earlier stages, and the Stage 2 drug vector is restricted to a plus-or-minus five-percent inheritance window. Candidate-specific context terms are counted in BIC, and the fully free context vector remains diagnostic only.</p>
  </article>

  <h2>How the four-stage handoff works</h2>
  <ol class="handoff">
    <li><strong>Stage 1 learns intrinsic growth.</strong> Each cell line is fitted jointly across 20k and 30k. Its winning family and effective <code>r/K/theta</code> values become the biological growth baseline.</li>
    <li><strong>Stage 2 learns treatment response.</strong> The Stage 1 growth law is fixed while all three monoculture doses and both densities determine delayed, Hill, transit, or population-balance parameters. A second joint BIC audit chooses whether resistant treatment needs onset, gradual activation, both, or timing shared with Naive.</li>
    <li><strong>Stage 3 learns untreated interaction.</strong> Stage 1 growth remains fixed while both densities and all mixture ratios determine competition and death parameters.</li>
    <li><strong>Stage 4 combines the mechanisms.</strong> Stage 1 growth and Stage 3 interaction are fixed. Stage 2 drug parameters are inherited within plus-or-minus five percent. BIC then tests strict inheritance, load effects, resistant-state shifts, tolerant growth plasticity, and a fully separate diagnostic vector across the combined treated data.</li>
  </ol>
</div>
</section>
""")

    appendix_blocks = [
        _appendix_ranking_html(stages[1].ranking; title = "Stage 1 untreated monoculture", by_cell_line = true,
            parameter_path = joinpath(csv_root, "monoculture_untreated", "monoculture_untreated_pooling_parameter_estimates.csv")),
        _appendix_ranking_html(stages[2].ranking; title = "Stage 2 treated monoculture", by_cell_line = true,
            parameter_path = joinpath(csv_root, "monoculture_treated", "monoculture_treated_joint_parameter_estimates.csv")),
        _appendix_ranking_html(timing_ranking_path; title = "Stage 2 timing audit",
            parameter_path = joinpath(csv_root, "monoculture_treated", "monoculture_treated_timing_hypothesis_parameters.csv")),
        _appendix_ranking_html(stages[3].ranking; title = "Stage 3 untreated coculture",
            parameter_path = joinpath(csv_root, "coculture_untreated", "coculture_untreated_joint_parameter_estimates.csv")),
        _appendix_ranking_html(stages[4].ranking; title = "Stage 4 treated coculture"),
    ]
    push!(sections, """
<section class="report-appendix" id="model-appendix">
<div class="stage-heading"><span>Appendix</span><h2>Complete model, BIC, and parameter audit</h2></div>
<p>This appendix is the archival view. Unlike the teaching tables above, it lists every finite tested candidate, its absolute BIC, Delta BIC, free-parameter count, boundary flag, and complete exported fitted parameter vector. Model IDs restart within each table and match the corresponding Delta-BIC chart.</p>
$(join(appendix_blocks))
</section>
""")

    html = """
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>A2780 staged model comparison</title>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>
<style>
:root { color-scheme: light; --ink: #1d2529; --muted: #56656a; --line: #cbd5d8; --soft: #eef3f4; --accent: #a21f3d; }
* { box-sizing: border-box; }
body { font-family: "Segoe UI", Arial, sans-serif; margin: 0; color: var(--ink); background: #fff; line-height: 1.45; }
header, section { width: min(1440px, calc(100% - 48px)); margin: 0 auto; }
header { padding: 40px 0 24px; border-bottom: 3px solid var(--ink); }
h1 { margin: 0 0 8px; font-size: 32px; letter-spacing: 0; }
h2 { margin: 2px 0 0; font-size: 24px; letter-spacing: 0; }
h3 { margin: 24px 0 8px; font-size: 17px; letter-spacing: 0; }
p { max-width: 980px; color: var(--muted); }
section { padding: 34px 0 48px; border-bottom: 1px solid var(--line); }
.stage-heading { display: flex; gap: 14px; align-items: baseline; }
.stage-heading span { color: var(--accent); font-size: 13px; font-weight: 700; text-transform: uppercase; }
.artifact-note { margin-bottom: 0; }
.table-wrap { width: 100%; overflow-x: auto; }
table { border-collapse: collapse; margin: 10px 0 24px; width: 100%; min-width: 860px; font-size: 13px; }
th, td { border: 1px solid var(--line); padding: 8px 10px; text-align: left; vertical-align: top; }
th { background: var(--soft); font-weight: 650; }
td:nth-child(4) { min-width: 360px; font-family: Consolas, "Courier New", monospace; font-size: 12px; overflow-wrap: anywhere; }
figure { margin: 28px 0 0; }
img { display: block; width: 100%; height: auto; border: 1px solid var(--line); }
figcaption { color: var(--muted); font-size: 12px; margin-top: 8px; }
.bic-figure { margin: 8px 0 34px; padding: 16px 18px; border: 1px solid var(--line); }
.bic-figure h4 { margin: 0 0 14px; font-size: 14px; }
.bic-axis { display: grid; gap: 9px; }
.bic-row { display: grid; grid-template-columns: minmax(220px, 0.36fr) minmax(180px, 1fr) 64px; gap: 10px; align-items: center; }
.bic-label { display: grid; grid-template-columns: 34px minmax(0, 1fr); gap: 5px; font-size: 12px; }
.bic-label span { overflow-wrap: anywhere; }
.bic-track { height: 14px; background: var(--soft); border-left: 2px solid var(--ink); }
.bic-bar { display: block; height: 100%; min-width: 3px; background: var(--accent); }
.bic-row output { font-variant-numeric: tabular-nums; font-size: 12px; }
.uncertainty-audit { margin-top: 28px; padding-top: 20px; border-top: 2px solid var(--line); }
.uncertainty-audit h3 { margin-top: 0; }
.audit-warning { padding: 10px 12px; border-left: 3px solid var(--accent); background: var(--soft); }
.report-appendix { width: min(1720px, calc(100% - 48px)); }
.appendix-table table { min-width: 1260px; }
.appendix-table td:last-child { min-width: 620px; font-family: Consolas, "Courier New", monospace; font-size: 11px; overflow-wrap: anywhere; }
.compact-parameters table { min-width: 760px; }
.notation-key { margin: 20px 0 24px; padding: 12px 0 14px; border-top: 2px solid var(--line); border-bottom: 1px solid var(--line); }
.notation-key h3 { margin: 0 0 10px; font-size: 15px; }
.notation-key dl { grid-template-columns: minmax(160px, 0.25fr) minmax(0, 1fr); }
.notation-key p { margin: 12px 0 0; }
.notation-key mjx-container { margin: 0 !important; }
.progression { border-top: 1px solid var(--line); margin-top: 24px; }
.progression article { padding: 18px 0; border-bottom: 1px solid var(--line); }
.progression article h3 { margin-top: 0; }
.progression article code { display: block; width: fit-content; max-width: 100%; margin: 7px 0; padding: 7px 9px; background: var(--soft); overflow-wrap: anywhere; }
.progression article p { margin-bottom: 0; }
.math { overflow-x: auto; overflow-y: hidden; margin: 10px 0; padding: 8px 4px; color: var(--ink); }
.math mjx-container[display="true"] { margin: 0.45em 0 !important; text-align: left !important; }
.equation-label { margin: 18px 0 2px; color: var(--ink); }
td mjx-container { font-size: 104% !important; }
.model-guide { margin-top: 42px; padding-top: 30px; border-top: 3px solid var(--ink); }
.model-guide > h2 { margin-top: 34px; }
.guide-grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 18px 32px; margin: 20px 0 38px; }
.guide-grid article { border-top: 2px solid var(--line); padding-top: 10px; }
.guide-grid h3 { margin-top: 0; }
dl { display: grid; grid-template-columns: minmax(120px, 0.34fr) minmax(0, 1fr); gap: 8px 14px; margin: 0; }
dt { font-weight: 650; }
dd { margin: 0; color: var(--muted); }
.model-family { padding: 18px 0 26px; border-bottom: 1px solid var(--line); }
.model-family code { display: block; width: fit-content; max-width: 100%; margin: 9px 0 4px; padding: 7px 9px; background: var(--soft); overflow-wrap: anywhere; }
.model-family p { margin: 4px 0 14px; }
.handoff { padding-left: 22px; max-width: 1100px; }
.handoff li { margin: 12px 0; padding-left: 5px; }
.missing { color: #8b1d2c; font-weight: 600; }
code { font-family: Consolas, "Courier New", monospace; }
@media (max-width: 700px) {
  header, section { width: min(100% - 28px, 1440px); }
  header { padding-top: 24px; }
  h1 { font-size: 26px; }
  h2 { font-size: 21px; }
  .stage-heading { display: block; }
  .guide-grid { grid-template-columns: minmax(0, 1fr); }
  .bic-row { grid-template-columns: minmax(0, 1fr) 50px; }
  .bic-label { grid-column: 1 / -1; }
  dl { grid-template-columns: minmax(0, 1fr); gap: 3px; }
  .notation-key dl { grid-template-columns: minmax(0, 1fr); }
  dd { max-width: 100%; margin-bottom: 8px; overflow-x: auto; overflow-y: hidden; }
}
    .back-home { display: inline-block; margin: 0 0 18px; color: var(--accent, #2563eb); font-weight: 700; text-decoration: none; }
    .back-home:hover { text-decoration: underline; }
</style>
</head>
<body>
$(join(sections, "\n"))
</body>
</html>
"""
    html_path = joinpath(report_dir, "a2780_staged_model_comparison.html")
    write(html_path, html)
    return html_path
end

"""Rebuild the staged overview and manifest from already validated fit artifacts."""
function refresh_a2780_output_summary!(; start::AbstractString = pwd())
    summaries = NamedTuple[]
    manifest_rows = NamedTuple[]
    for condition in STAGED_A2780_CONDITIONS
        out = IOUtils.condition_output_dirs(condition; start = start)
        decoded_path = joinpath(out.csv, "$(condition)_a2780_decoded.csv")
        pooling_ranking_path = joinpath(out.csv, "$(condition)_pooling_model_ranking.csv")
        ranking_path = isfile(pooling_ranking_path) ? pooling_ranking_path : joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
        outputs = filter(isfile, [decoded_path, ranking_path])
        decoded_rows = isfile(decoded_path) ? nrow(CSV.read(decoded_path, DataFrame)) : 0
        ranking = isfile(ranking_path) ? _valid_ranking(CSV.read(ranking_path, DataFrame)) : DataFrame()
        status = isempty(ranking) ? "failed" : "completed"
        message = isempty(ranking) ? "No finite validated ranking artifact." : "Reused validated executed outputs."
        best_model = ""
        best_bic = NaN
        best_ssr = NaN
        if !isempty(ranking)
            metric = _ranking_metric_col(ranking)
            sort!(ranking, metric)
            best_model = String(first(ranking.model))
            best_bic = :bic in propertynames(ranking) ? Float64(first(ranking.bic)) : NaN
            best_ssr = :ssr in propertynames(ranking) ? Float64(first(ranking.ssr)) : (:sse in propertynames(ranking) ? Float64(first(ranking.sse)) : NaN)
        end
        push!(summaries, (
            condition = condition,
            status = status,
            decoded_rows = decoded_rows,
            fit_rows = nrow(ranking),
            best_model = best_model,
            best_bic = best_bic,
            best_ssr = best_ssr,
            message = message,
        ))
        push!(manifest_rows, (
            timestamp_utc = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
            condition = condition,
            status = status,
            message = message,
            output_count = length(outputs),
            outputs = join(outputs, ";"),
        ))
    end
    overview = DataFrame(summaries)
    overview_path = _write_notebook_summary(overview; start = start)
    manifest_path = _write_stage_manifest(manifest_rows; start = start)
    return (overview = overview, overview_path = overview_path, manifest_path = manifest_path)
end

end
