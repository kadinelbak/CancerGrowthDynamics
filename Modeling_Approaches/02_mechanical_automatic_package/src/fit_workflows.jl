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
    untreated_baseline = nothing,
    max_time_per_fit::Float64 = 45.0,
)
    r0 = isnothing(untreated_baseline) ? 0.3 : untreated_baseline.params[1]
    K0 = isnothing(untreated_baseline) ? 3000.0 : untreated_baseline.params[2]

    specs = ModelRegistry.local_model_specs()

    # Detect the measurement column (area or first numeric non-metadata col)
    meta_cols = Set([:dose, :time, :day, :replicate, :cell_line, :density, :sample_id, :file])
    numeric_cols = [c for c in propertynames(decoded)
                    if !(c in meta_cols) && eltype(decoded[!, c]) <: Real]
    time_col = :time
    :day in propertynames(decoded) && !(:time in propertynames(decoded)) && (time_col = :day)
    area_col = isempty(numeric_cols) ? error("No numeric measurement column found in decoded data") :
                                       numeric_cols[1]

    rows = NamedTuple{(:model, :dose, :bic, :ssr, :params),
                      Tuple{String, Float64, Float64, Float64, String}}[]

    for dose_val in sort(unique(decoded[!, :dose]))
        sub = filter(r -> r[:dose] == dose_val, decoded)
        grouped = combine(groupby(sub, time_col), area_col => mean => :y_mean)
        sort!(grouped, time_col)
        x = Float64.(grouped[!, time_col])
        y = Float64.(grouped[!, :y_mean])

        for spec in specs
            println("  Fitting $(spec.name) at dose=$(dose_val)...")
            fit = _fit_spec_to_dose(spec, x, y, dose_val, r0, K0; max_time = max_time_per_fit)
            push!(rows, (
                model  = spec.name,
                dose   = dose_val,
                bic    = fit.bic,
                ssr    = fit.ssr,
                params = string(fit.params),
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
                                          untreated_baseline = untreated_baseline,
                                          max_time_per_fit = max_time_per_fit)
        ranking_path = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
        best_path    = joinpath(out.csv, "$(condition)_automatic_best_models_top10.csv")
        CSV.write(ranking_path, ranking_df)
        CSV.write(best_path, first(ranking_df, min(10, nrow(ranking_df))))

        IOUtils.write_manifest_row(
            condition = condition,
            step      = "fit",
            outputs   = filter(!isnothing, [ranking_path, best_path, staged_baseline_path]),
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

    # Untreated conditions: use the package's built-in compare_models_dict
    model_map = Dict(
        "monoculture_untreated" => ["logistic_growth", "gompertz_growth"],
        "coculture_untreated"   => ["logistic_growth", "gompertz_growth"],
    )
    include_models = get(model_map, condition, ["logistic_growth", "gompertz_growth"])

    # Detect columns
    time_col = :time
    :day in propertynames(decoded) && !(:time in propertynames(decoded)) && (time_col = :day)
    meta_cols = Set([:dose, :time, :day, :replicate, :cell_line, :density, :sample_id, :file])
    numeric_cols = [c for c in propertynames(decoded) if !(c in meta_cols) && eltype(decoded[!, c]) <: Real]
    area_col = isempty(numeric_cols) ? error("No numeric measurement column") : numeric_cols[1]

    grouped = combine(groupby(decoded, time_col), area_col => mean => :y_mean)
    sort!(grouped, time_col)
    x = Float64.(grouped[!, time_col])
    y = Float64.(grouped[!, :y_mean])

    builtin_map = Dict(
        "logistic_growth" => GrowthParameterEstimation.logistic_growth!,
        "gompertz_growth" => GrowthParameterEstimation.gompertz_growth!,
    )

    specs_dict = Dict{String,NamedTuple}()
    for mname in include_models
        haskey(builtin_map, mname) || continue
        specs_dict[mname] = (
            model  = builtin_map[mname],
            p0     = [0.3, max(y[end], 100.0)],
            bounds = [(1e-6, 5.0), (1.0, 1e7)],
        )
    end

    fits = GrowthParameterEstimation.compare_models_dict(
        x, y, specs_dict;
        show_stats = false,
        output_csv = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv"),
    )

    ranking_rows = [(model = n, bic = v.bic, ssr = v.ssr, params = string(v.params))
                    for (n, v) in pairs(fits)]
    ranking_df   = sort!(DataFrame(ranking_rows), :bic)
    ranking_path = joinpath(out.csv, "$(condition)_automatic_model_ranking.csv")
    best_path    = joinpath(out.csv, "$(condition)_automatic_best_models_top10.csv")
    CSV.write(ranking_path, ranking_df)
    CSV.write(best_path, first(ranking_df, min(10, nrow(ranking_df))))

    IOUtils.write_manifest_row(
        condition = condition,
        step      = "fit",
        outputs   = filter(!isnothing, [ranking_path, best_path, staged_baseline_path]),
        start     = start,
    )
    return (
        result              = fits,
        ranking             = ranking_df,
        ranking_path        = ranking_path,
        best_path           = best_path,
        untreated_baseline  = untreated_baseline,
        staged_baseline_path = staged_baseline_path,
        used_cached_results = false,
        fit_api_available   = true,
    )
end

end
