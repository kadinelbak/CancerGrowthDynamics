module FitWorkflows

using CSV
using DataFrames
using Statistics
using OrdinaryDiffEq
using BlackBoxOptim
using GrowthParameterEstimation

using ..IOUtils
using ..ModelRegistry

export load_untreated_baseline, run_condition_fit!

# ---------------------------------------------------------------------------
# Untreated baseline loading
# ---------------------------------------------------------------------------

function load_untreated_baseline(; start::AbstractString = pwd())
    root = IOUtils.package_root(start)
    base_dir = joinpath(root, "outputs", "csv", "monoculture_untreated")
    summary_path = joinpath(base_dir, "best_model_summary.csv")
    params_path = joinpath(base_dir, "params", "best_params.csv")

    if !isfile(summary_path) || !isfile(params_path)
        return nothing
    end

    summary_df = CSV.read(summary_path, DataFrame)
    params_df  = CSV.read(params_path,  DataFrame)
    summary_map = Dict(string(row.key) => row.value for row in eachrow(summary_df))

    model_name = get(summary_map, "best_model", missing)
    params     = Float64.(params_df.value)

    (ismissing(model_name) || isempty(params)) && return nothing

    return (
        model       = String(model_name),
        params      = params,
        summary_path = summary_path,
        params_path  = params_path,
    )
end

# ---------------------------------------------------------------------------
# Write staged baseline CSV (for traceability)
# ---------------------------------------------------------------------------

function _write_staged_baseline(out, untreated_baseline)
    untreated_baseline === nothing && return nothing

    rows = Tuple{String,String}[("baseline_model", untreated_baseline.model)]
    append!(rows, [("param_$(i)", string(untreated_baseline.params[i])) for i in eachindex(untreated_baseline.params)])

    staged_path = joinpath(out.csv, "staged_untreated_baseline.csv")
    CSV.write(staged_path, DataFrame(key = first.(rows), value = last.(rows)))
    return staged_path
end

_nonmissing_strings(vals) = String[v for v in vals if v !== nothing]

function _load_untreated_monoculture_cellline_baselines(; start::AbstractString = pwd())
    root = IOUtils.find_repo_root(start)
    auto_path = joinpath(IOUtils.package_root(start), "outputs", "csv", "monoculture_untreated", "untreated_group_baselines.csv")
    if isfile(auto_path)
        df_auto = CSV.read(auto_path, DataFrame)
        required_auto = Set([:cell_line, :density, :r, :K])
        if all(c -> c in names(df_auto), required_auto)
            by_cell_density = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
            by_cell = Dict{String,Tuple{Float64,Float64}}()

            for row in eachrow(df_auto)
                cell = String(row.cell_line)
                density = String(row.density)
                r = Float64(row.r)
                K = Float64(row.K)
                by_cell_density[(cell, density)] = (r, K)
            end

            cell_groups = groupby(df_auto, :cell_line)
            for g in cell_groups
                cell = String(first(g.cell_line))
                by_cell[cell] = (mean(Float64.(g.r)), mean(Float64.(g.K)))
            end

            return (path = auto_path, by_cell_density = by_cell_density, by_cell = by_cell)
        end
    end

    candidates = [
        joinpath(root, "Modeling_Approaches", "01_mechanical_manual", "Untreated MonoCulture", "outputs", "csv", "untreated_monoculture_logistic_parameters.csv"),
        joinpath(root, "Modeling_Approaches", "01_mechanical_manual", "Untreated MonoCulture", "untreated_monoculture_logistic_parameters.csv"),
    ]
    existing = filter(isfile, candidates)
    isempty(existing) && return nothing
    path = first(existing)

    df = CSV.read(path, DataFrame)
    required = Set([:CellLine, :Density, :r, :K])
    all(c -> c in names(df), required) || return nothing

    by_cell_density = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
    by_cell = Dict{String,Tuple{Float64,Float64}}()

    for row in eachrow(df)
        cell = String(row.CellLine)
        density = String(row.Density)
        r = Float64(row.r)
        K = Float64(row.K)
        by_cell_density[(cell, density)] = (r, K)
    end

    cell_groups = groupby(df, :CellLine)
    for g in cell_groups
        cell = String(first(g.CellLine))
        by_cell[cell] = (mean(Float64.(g.r)), mean(Float64.(g.K)))
    end

    return (path = path, by_cell_density = by_cell_density, by_cell = by_cell)
end

function _baseline_rk_for_group(
    baseline_map,
    cell_line::Union{Missing,AbstractString},
    density::Union{Missing,AbstractString},
    fallback_r::Float64,
    fallback_K::Float64,
)
    baseline_map === nothing && return (fallback_r, fallback_K)
    ismissing(cell_line) && return (fallback_r, fallback_K)

    cell = String(cell_line)
    den = ismissing(density) ? "" : String(density)

    if haskey(baseline_map.by_cell_density, (cell, den))
        return baseline_map.by_cell_density[(cell, den)]
    end
    if haskey(baseline_map.by_cell, cell)
        return baseline_map.by_cell[cell]
    end
    return (fallback_r, fallback_K)
end

function _fit_anchored_linear_kill(
    x::Vector{Float64},
    y::Vector{Float64},
    dose::Float64,
    r0::Float64,
    K0::Float64;
    solver = Rodas5(),
    max_time::Float64 = 12.0,
)
    function ode!(du, u, p, t)
        k_kill = p[1]
        N = max(u[1], 0.0)
        growth = r0 * N * max(0.0, 1 - N / max(K0, 1e-8))
        kill = k_kill * max(dose, 0.0) * N
        du[1] = growth - kill
    end

    u0 = [max(y[1], 1.0)]
    prob = ODEProblem(ode!, u0, (x[1], x[end]), [0.05])

    function loss(p_vec)
        pv = Vector{Float64}(p_vec)
        try
            sol = solve(remake(prob; p = pv), solver;
                        reltol = 1e-6, abstol = 1e-6,
                        saveat = x, maxiters = 50_000)
            sol.retcode == ReturnCode.Success || return 1e12
            yhat = first.(sol.u)
            return sum((y .- yhat) .^ 2)
        catch
            return 1e12
        end
    end

    result = bboptimize(
        loss;
        SearchRange = [(0.0, 8.0)],
        NumDimensions = 1,
        Method = :de_rand_1_bin,
        MaxTime = max_time,
        TraceMode = :silent,
    )
    p_opt = Vector{Float64}(result.archive_output.best_candidate)

    bic, ssr = try
        sol2 = solve(remake(prob; p = p_opt), solver;
                     reltol = 1e-10, abstol = 1e-10, saveat = x)
        yhat2 = first.(sol2.u)
        ssr2 = sum((y .- yhat2) .^ 2)
        n, k = length(x), length(p_opt)
        n * log(max(ssr2, 1e-20) / n) + k * log(n), ssr2
    catch
        1e12, 1e12
    end

    return (params = p_opt, bic = bic, ssr = ssr)
end

function _fit_anchored_hill_kill(
    x::Vector{Float64},
    y::Vector{Float64},
    dose::Float64,
    r0::Float64,
    K0::Float64;
    solver = Rodas5(),
    max_time::Float64 = 12.0,
)
    function ode!(du, u, p, t)
        emax, ec50, hill_n = p
        N = max(u[1], 0.0)
        C = max(dose, 0.0)
        kill = emax * (C^hill_n / (ec50^hill_n + C^hill_n + 1e-12))
        growth = r0 * N * max(0.0, 1 - N / max(K0, 1e-8))
        du[1] = growth - kill * N
    end

    u0 = [max(y[1], 1.0)]
    prob = ODEProblem(ode!, u0, (x[1], x[end]), [0.3, max(dose, 0.1), 1.0])

    function loss(p_vec)
        pv = Vector{Float64}(p_vec)
        try
            sol = solve(remake(prob; p = pv), solver;
                        reltol = 1e-6, abstol = 1e-6,
                        saveat = x, maxiters = 50_000)
            sol.retcode == ReturnCode.Success || return 1e12
            yhat = first.(sol.u)
            return sum((y .- yhat) .^ 2)
        catch
            return 1e12
        end
    end

    result = bboptimize(
        loss;
        SearchRange = [(0.0, 1.5), (1e-3, 10.0), (0.1, 6.0)],
        NumDimensions = 3,
        Method = :de_rand_1_bin,
        MaxTime = max_time,
        TraceMode = :silent,
    )
    p_opt = Vector{Float64}(result.archive_output.best_candidate)

    bic, ssr = try
        sol2 = solve(remake(prob; p = p_opt), solver;
                     reltol = 1e-10, abstol = 1e-10, saveat = x)
        yhat2 = first.(sol2.u)
        ssr2 = sum((y .- yhat2) .^ 2)
        n, k = length(x), length(p_opt)
        n * log(max(ssr2, 1e-20) / n) + k * log(n), ssr2
    catch
        1e12, 1e12
    end

    return (params = p_opt, bic = bic, ssr = ssr)
end

# ---------------------------------------------------------------------------
# Per-dose model fitting using local model specs
# ---------------------------------------------------------------------------

"""
Fit a single LocalModelSpec to one dose group (x, y time-series).
Returns (params, bic, ssr).
"""
function _fit_spec_to_dose(
    spec::ModelRegistry.LocalModelSpec,
    x::Vector{Float64},
    y::Vector{Float64},
    dose::Float64,
    r0::Float64,
    K0::Float64;
    solver   = Rodas5(),
    max_time::Float64 = 45.0,
)
    p0  = spec.p0_factory(r0, K0, dose)
    obs = spec.observable          # (u::Vector) -> scalar
    nst = spec.n_states

    # Constant-exposure closure, baked-in for this dose group
    exposure_fn = Returns(dose)
    ode4! = (du, u, p, t) -> spec.ode!(du, u, p, t, exposure_fn)

    u0    = zeros(nst)
    u0[1] = max(y[1], 1.0)
    tspan = (x[1], x[end])
    prob  = ODEProblem(ode4!, u0, tspan, p0)

    n_dims = length(p0)

    function loss(p_vec)
        pv = Vector{Float64}(p_vec)
        try
            sol = solve(remake(prob; p = pv), solver;
                        reltol = 1e-6, abstol = 1e-6,
                        saveat = x, maxiters = 50_000)
            sol.retcode == ReturnCode.Success || return 1e12
            yhat = [obs(u) for u in sol.u]
            length(yhat) == length(y) || return 1e12
            return sum((y .- yhat) .^ 2)
        catch
            return 1e12
        end
    end

    result = bboptimize(
        loss;
        SearchRange  = collect(zip(first.(spec.bounds), last.(spec.bounds))),
        NumDimensions = n_dims,
        Method       = :de_rand_1_bin,
        MaxTime      = max_time,
        TraceMode    = :silent,
    )
    p_opt = Vector{Float64}(result.archive_output.best_candidate)

    # Compute final BIC/SSR with tighter tolerances
    bic, ssr = try
        sol2 = solve(remake(prob; p = p_opt), solver;
                     reltol = 1e-10, abstol = 1e-10, saveat = x)
        yhat2 = [obs(u) for u in sol2.u]
        ssr2 = sum((y .- yhat2) .^ 2)
        n, k = length(x), length(p_opt)
        bic2 = n * log(max(ssr2, 1e-20) / n) + k * log(n)
        bic2, ssr2
    catch
        1e12, 1e12
    end

    return (params = p_opt, bic = bic, ssr = ssr)
end

"""
Fit all local model specs to each dose group in `decoded` and return a
ranked DataFrame with columns: model, dose, bic, ssr, params.
"""
function _run_treated_fitting(
    decoded::DataFrame,
    out;
    start::AbstractString = pwd(),
    untreated_baseline = nothing,
    max_time_per_fit::Float64 = 45.0,
)
    r0 = isnothing(untreated_baseline) ? 0.3 : untreated_baseline.params[1]
    K0 = isnothing(untreated_baseline) ? 3000.0 : untreated_baseline.params[2]
    baseline_map = _load_untreated_monoculture_cellline_baselines(; start = start)

    specs = ModelRegistry.local_model_specs()

    # Detect the measurement column (area or first numeric non-metadata col)
    meta_cols = Set([:dose, :time, :day, :replicate, :cell_line, :density, :sample_id, :file])
    numeric_cols = [c for c in propertynames(decoded)
                    if !(c in meta_cols) && eltype(decoded[!, c]) <: Real]
    time_col = :time
    :day in propertynames(decoded) && !(:time in propertynames(decoded)) && (time_col = :day)
    area_col = isempty(numeric_cols) ? error("No numeric measurement column found in decoded data") :
                                       numeric_cols[1]

    rows = NamedTuple{(:model, :dose, :bic, :ssr, :params, :cell_line, :density),
                      Tuple{String, Float64, Float64, Float64, String, String, String}}[]

    colset = Set(Symbol.(names(decoded)))
    cell_col = :cell_line in colset ? :cell_line : (:cellline in colset ? :cellline : nothing)
    density_col = :density in colset ? :density : nothing
    has_cell = cell_col !== nothing
    has_density = density_col !== nothing

    gcols = Any[]
    has_cell && push!(gcols, cell_col)
    has_density && push!(gcols, density_col)
    push!(gcols, :dose)

    for grp in groupby(decoded, gcols)
        dose_val = Float64(first(grp.dose))
        grouped = combine(groupby(grp, time_col), area_col => mean => :y_mean)
        sort!(grouped, time_col)
        x = Float64.(grouped[!, time_col])
        y = Float64.(grouped[!, :y_mean])

        cell_label = has_cell ? String(first(grp[!, cell_col])) : "pooled"
        den_label = has_density ? String(first(grp[!, density_col])) : "pooled"
        r_group, K_group = _baseline_rk_for_group(baseline_map, cell_label, den_label, r0, K0)

        for spec in specs
            println("  Fitting $(spec.name) at dose=$(dose_val), cell=$(cell_label), density=$(den_label)...")
            fit = _fit_spec_to_dose(spec, x, y, dose_val, r_group, K_group; max_time = max_time_per_fit)
            push!(rows, (
                model  = spec.name,
                dose   = dose_val,
                bic    = fit.bic,
                ssr    = fit.ssr,
                params = string((fit = fit.params, r_anchor = r_group, K_anchor = K_group)),
                cell_line = cell_label,
                density = den_label,
            ))
        end
    end

    # Additional fast anchored models using untreated monoculture r/K per cell line and density.
    # These are fitted by (cell_line, density, dose) groups and only estimate drug-effect terms.
    if has_cell
        gcols = has_density ? [cell_col, density_col, :dose] : [cell_col, :dose]
        for grp in groupby(decoded, gcols)
            dose_val = Float64(first(grp.dose))
            grouped = combine(groupby(grp, time_col), area_col => mean => :y_mean)
            sort!(grouped, time_col)
            x = Float64.(grouped[!, time_col])
            y = Float64.(grouped[!, :y_mean])

            cell_label = String(first(grp[!, cell_col]))
            den_label = has_density ? String(first(grp[!, density_col])) : ""
            r_anchor, K_anchor = _baseline_rk_for_group(baseline_map, cell_label, den_label, r0, K0)
            anchored_budget = max(1.0, min(max_time_per_fit / 3, 4.0))

            println("  Fitting anchored_linear_kill at dose=$(dose_val), cell=$(cell_label), density=$(den_label)...")
            lin_fit = _fit_anchored_linear_kill(x, y, dose_val, r_anchor, K_anchor; max_time = anchored_budget)
            push!(rows, (
                model = "anchored_linear_kill",
                dose = dose_val,
                bic = lin_fit.bic,
                ssr = lin_fit.ssr,
                params = string((k_kill = lin_fit.params[1], r_anchor = r_anchor, K_anchor = K_anchor)),
                cell_line = cell_label,
                density = den_label,
            ))

            println("  Fitting anchored_hill_kill at dose=$(dose_val), cell=$(cell_label), density=$(den_label)...")
            hill_fit = _fit_anchored_hill_kill(x, y, dose_val, r_anchor, K_anchor; max_time = anchored_budget)
            push!(rows, (
                model = "anchored_hill_kill",
                dose = dose_val,
                bic = hill_fit.bic,
                ssr = hill_fit.ssr,
                params = string((emax = hill_fit.params[1], ec50 = hill_fit.params[2], hill_n = hill_fit.params[3], r_anchor = r_anchor, K_anchor = K_anchor)),
                cell_line = cell_label,
                density = den_label,
            ))
        end
    end

    return sort!(DataFrame(rows), :bic)
end

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

function run_condition_fit!(
    decoded::DataFrame,
    condition::AbstractString;
    start::AbstractString = pwd(),
    untreated_baseline = nothing,
    max_time_per_fit::Float64 = 45.0,
)
    ModelRegistry.ensure_model_registry!()
    out = IOUtils.condition_output_dirs(condition; start)

    if condition in ("monoculture_treated", "coculture_treated") && untreated_baseline === nothing
        untreated_baseline = load_untreated_baseline(; start = start)
    end

    staged_baseline_path = _write_staged_baseline(out, untreated_baseline)

    # Use local fitting loop for treated conditions
    if condition in ("monoculture_treated", "coculture_treated")
        ranking_df = _run_treated_fitting(decoded, out;
                                          start = start,
                                          untreated_baseline = untreated_baseline,
                                          max_time_per_fit = max_time_per_fit)
        ranking_path = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
        best_path    = joinpath(out.csv, "$(condition)_automatic_best_models_top10.csv")
        CSV.write(ranking_path, ranking_df)
        CSV.write(best_path, first(ranking_df, min(10, nrow(ranking_df))))

        IOUtils.write_manifest_row(
            condition = condition,
            step      = "fit",
            outputs   = _nonmissing_strings([ranking_path, best_path, staged_baseline_path]),
            start     = start,
        )
        return (
            result              = nothing,
            ranking             = ranking_df,
            ranking_path        = ranking_path,
            best_path           = best_path,
            untreated_baseline  = untreated_baseline,
            staged_baseline_path = staged_baseline_path,
            used_cached_results = false,
            fit_api_available   = true,
        )
    end

    # Untreated conditions: fit per experimental subgroup so r/K can carry over
    # to treated fits by cell line and density.
    model_map = Dict(
        "monoculture_untreated" => ["logistic_growth"],
        "coculture_untreated"   => ["logistic_growth"],
    )
    include_models = get(model_map, condition, ["logistic_growth"])

    # Detect columns
    time_col = :time
    :day in propertynames(decoded) && !(:time in propertynames(decoded)) && (time_col = :day)
    meta_cols = Set([:dose, :time, :day, :replicate, :cell_line, :density, :sample_id, :file])
    numeric_cols = [c for c in propertynames(decoded) if !(c in meta_cols) && eltype(decoded[!, c]) <: Real]
    area_col = isempty(numeric_cols) ? error("No numeric measurement column") : numeric_cols[1]

    builtin_map = Dict(
        "logistic_growth" => GrowthParameterEstimation.logistic_growth!,
        "gompertz_growth" => GrowthParameterEstimation.gompertz_growth!,
    )

    colset = Set(Symbol.(names(decoded)))
    cell_col = :cell_line in colset ? :cell_line : (:cellline in colset ? :cellline : nothing)
    density_col = :density in colset ? :density : nothing

    subgroup_cols = Symbol[]
    cell_col !== nothing && push!(subgroup_cols, cell_col)
    density_col !== nothing && push!(subgroup_cols, density_col)

    subgroups = isempty(subgroup_cols) ? [decoded] : collect(groupby(decoded, subgroup_cols))

    all_candidates = NamedTuple[]
    for (idx, subdf) in enumerate(subgroups)
        grouped = combine(groupby(subdf, time_col), area_col => mean => :y_mean)
        sort!(grouped, time_col)
        x = Float64.(grouped[!, time_col])
        y = Float64.(grouped[!, :y_mean])

        specs_dict = Dict{String,NamedTuple}()
        for mname in include_models
            haskey(builtin_map, mname) || continue
            specs_dict[mname] = (
                model  = builtin_map[mname],
                p0     = [0.3, max(y[end], 100.0)],
                bounds = [(1e-6, 5.0), (1.0, 1e7)],
            )
        end

        subgroup_tmp = joinpath(out.csv, "tables", "subgroup_$(condition)_$(idx)_ranking.csv")
        mkpath(dirname(subgroup_tmp))
        fits = GrowthParameterEstimation.compare_models_dict(
            x, y, specs_dict;
            show_stats = false,
            output_csv = subgroup_tmp,
        )

        cell_label = cell_col === nothing ? "pooled" : String(first(subdf[!, cell_col]))
        den_label = density_col === nothing ? "pooled" : String(first(subdf[!, density_col]))

        for (mname, v) in pairs(fits)
            push!(all_candidates, (
                model = String(mname),
                bic = Float64(v.bic),
                ssr = Float64(v.ssr),
                params_vec = Vector{Float64}(v.params),
                params = string(v.params),
                cell_line = cell_label,
                density = den_label,
            ))
        end
    end

    ranking_df = sort!(DataFrame(
        model = [c.model for c in all_candidates],
        bic = [c.bic for c in all_candidates],
        ssr = [c.ssr for c in all_candidates],
        params = [c.params for c in all_candidates],
        cell_line = [c.cell_line for c in all_candidates],
        density = [c.density for c in all_candidates],
    ), :bic)

    # Build per-subgroup untreated baseline table used by treated models.
    subgroup_baselines = NamedTuple[]
    grouped_candidates = groupby(ranking_df, [:cell_line, :density])
    for g in grouped_candidates
        g_sorted = sort(g, :bic)
        best = first(g_sorted, 1)
        params_txt = String(best.params[1])
        # Parse out first two numeric entries from the parameter vector string.
        nums = [parse(Float64, m.match) for m in eachmatch(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", params_txt)]
        if length(nums) >= 2
            push!(subgroup_baselines, (
                cell_line = String(best.cell_line[1]),
                density = String(best.density[1]),
                best_model = String(best.model[1]),
                r = nums[1],
                K = nums[2],
                bic = Float64(best.bic[1]),
                ssr = Float64(best.ssr[1]),
            ))
        end
    end

    subgroup_baseline_df = DataFrame(subgroup_baselines)
    subgroup_baseline_path = joinpath(out.csv, "untreated_group_baselines.csv")
    CSV.write(subgroup_baseline_path, subgroup_baseline_df)

    if !isempty(all_candidates)
        best_idx = argmin([c.bic for c in all_candidates])
        best = all_candidates[best_idx]
        summary_path = joinpath(out.csv, "best_model_summary.csv")
        params_dir = joinpath(out.csv, "params")
        mkpath(params_dir)
        params_path = joinpath(params_dir, "best_params.csv")
        CSV.write(summary_path, DataFrame(key = ["best_model", "best_cell_line", "best_density"], value = [best.model, best.cell_line, best.density]))
        CSV.write(params_path, DataFrame(param = ["p$(i)" for i in eachindex(best.params_vec)], value = best.params_vec))
    end

    ranking_path = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
    best_path    = joinpath(out.csv, "$(condition)_automatic_best_models_top10.csv")
    CSV.write(ranking_path, ranking_df)
    CSV.write(best_path, first(ranking_df, min(10, nrow(ranking_df))))

    IOUtils.write_manifest_row(
        condition = condition,
        step      = "fit",
        outputs   = _nonmissing_strings([ranking_path, best_path, subgroup_baseline_path, staged_baseline_path]),
        start     = start,
    )
    return (
        result              = nothing,
        ranking             = ranking_df,
        ranking_path        = ranking_path,
        best_path           = best_path,
        subgroup_baseline_path = subgroup_baseline_path,
        untreated_baseline  = untreated_baseline,
        staged_baseline_path = staged_baseline_path,
        used_cached_results = false,
        fit_api_available   = true,
    )
end

end
