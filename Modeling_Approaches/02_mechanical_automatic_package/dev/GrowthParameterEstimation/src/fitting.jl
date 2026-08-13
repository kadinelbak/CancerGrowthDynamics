# Fitting module - Contains all ODE fitting and comparison functions
module Fitting

using StatsBase
using CSV
using DataFrames
using DifferentialEquations
using LsqFit
using RecursiveArrayTools
using DiffEqParamEstim
using Optimization
using ForwardDiff
using OptimizationOptimJL

using ..Models

export setUpProblem, calculate_bic, pQuickStat, run_single_fit,
    compare_models, compare_datasets, compare_models_dict, fit_three_datasets,
    run_joint_fit, run_joint_multistart, profile_joint_fit_bounds, profile_joint_fit_bounds_two_sided,
    compare_joint_models_dict

"""
    setUpProblem(model, x, y, solver, u0, p0, tspan, bounds; max_time=100.0, maxiters=10_000)

Set up and solve an ODE fitting problem using Optimization.jl (BFGS/Fminbox).
Returns the optimized parameters, the dense solution, and the optimized problem.
"""
function setUpProblem(model, x, y, solver, u0, p0, tspan, bounds; max_time::Real = 100.0, maxiters::Integer = 10_000)
    p0_vec = collect(p0)
    prob = ODEProblem(model, u0, tspan, p0_vec)
    ndims = length(p0_vec)

    loss = build_loss_objective(
        prob, solver,
        L2Loss(x, y),
        Optimization.AutoForwardDiff();
        maxiters = maxiters,
        verbose  = false,
    )

    if bounds === nothing
        optprob = Optimization.OptimizationProblem(loss, p0_vec)
        result  = Optimization.solve(optprob, OptimizationOptimJL.BFGS(); maxiters = maxiters)
    else
        lb = Float64[first(b) for b in bounds]
        ub = Float64[last(b) for b in bounds]
        optprob = Optimization.OptimizationProblem(loss, p0_vec; lb = lb, ub = ub)
        result  = Optimization.solve(optprob, OptimizationOptimJL.Fminbox(OptimizationOptimJL.BFGS()); maxiters = maxiters)
    end

    p_opt = collect(result.u)
    if length(p_opt) != ndims
        error("Optimizer returned $(length(p_opt)) parameters but $(ndims) were expected. Check bounds/initial guess configuration.")
    end
    prob_opt = remake(prob; p = p_opt, u0 = [y[1]], tspan = tspan)
    x_dense  = range(x[1], x[end], length = 1000)
    sol_opt  = solve(prob_opt, solver; reltol = 1e-12, abstol = 1e-12, saveat = x_dense)

    return p_opt, sol_opt, prob_opt
end

"""
    calculate_bic(prob, x, y, solver, p)

Compute BIC and SSR for a solved ODE model with parameters `p`.
"""
function calculate_bic(prob, x, y, solver, p)
    sol   = solve(remake(prob; p = p), solver; reltol = 1e-15, abstol = 1e-15, saveat = x)
    resid = y .- getindex.(sol.u, 1)
    ssr   = sum(resid .^ 2)
    k     = length(p)
    n     = length(x)
    bic   = n * log(ssr / n) + k * log(n)
    return bic, ssr
end

"""
    pQuickStat(x, y, p, sol, prob, bic, ssr)

Print a short summary of optimized parameters and fit statistics.
"""
function pQuickStat(x, y, p, sol, prob, bic, ssr)
    println("Optimized params: ", p)
    println("SSR: ", ssr)
    println("BIC: ", bic)
end

"""
    run_single_fit(x, y, p0; model=Models.logistic_growth!, fixed_params=nothing,
                   solver=Rodas5(), bounds=nothing, max_time=100.0, show_stats=true)

Fit a single model to `x`, `y` data with optional fixed parameters and bounds.
"""
function run_single_fit(
    x::Vector{<:Real},
    y::Vector{<:Real},
    p0::Vector{<:Real};
    model            = Models.logistic_growth!,
    fixed_params     = nothing,
    solver           = Tsit5(),
    bounds           = nothing,
    max_time::Real   = 100.0,
    show_stats::Bool = true,
)
    # Handle fixed parameters by wrapping the model
    if fixed_params !== nothing
        original_model = model
        n_total_params = length(p0) + length(fixed_params)

        model = function(du, u, p_free, t)
            p_full  = zeros(n_total_params)
            free_ix = 1
            for i in 1:n_total_params
                if haskey(fixed_params, i)
                    p_full[i] = fixed_params[i]
                else
                    p_full[i] = p_free[free_ix]
                    free_ix  += 1
                end
            end
            original_model(du, u, p_full, t)
        end

        free_indices = [i for i in 1:length(p0) if !haskey(fixed_params, i)]
        p0 = p0[free_indices]
        if bounds !== nothing
            bounds = bounds[free_indices]
        end
    end

    nparams = length(p0)

    _default_upper(p::Real) = max(10.0, abs(Float64(p)) * 100)

    if bounds === nothing
        bounds = [(0.0, _default_upper(p0[i])) for i in 1:nparams]
    else
        length(bounds) == nparams || throw(ArgumentError("bounds must have length $nparams"))
        bounds = [
            (
                isfinite(b[1]) ? Float64(b[1]) : 0.0,
                isfinite(b[2]) ? Float64(b[2]) : _default_upper(p0[i]),
            ) for (i, b) in enumerate(bounds)
        ]
    end

    x      = Float64.(x)
    y      = Float64.(y)
    tspan  = (x[1], x[end])
    u0     = [y[1]]

    p_opt, sol_opt, prob_opt = setUpProblem(model, x, y, solver, u0, p0, tspan, bounds; max_time = max_time)
    bic, ssr = calculate_bic(prob_opt, x, y, solver, p_opt)
    show_stats && pQuickStat(x, y, p_opt, sol_opt, prob_opt, bic, ssr)

    return (params = p_opt, bic = bic, ssr = ssr, solution = sol_opt)
end

"""
    compare_models(x, y, name1, model1, p0_1, name2, model2, p0_2; ...)

Fit two models to the same dataset and return a summary plus the best model.
"""
function compare_models(
    x::Vector{<:Real},
    y::Vector{<:Real},
    name1::String, model1::Function, p0_1::Vector{<:Real},
    name2::String, model2::Function, p0_2::Vector{<:Real};
    solver             = Rodas5(),
    bounds1            = nothing,
    bounds2            = nothing,
    fixed_params1      = nothing,
    fixed_params2      = nothing,
    show_stats::Bool   = false,
    output_csv::String = "model_comparison.csv",
)
    fit1 = run_single_fit(x, y, p0_1; model = model1, fixed_params = fixed_params1,
                          solver = solver, bounds = bounds1, show_stats = show_stats)
    fit2 = run_single_fit(x, y, p0_2; model = model2, fixed_params = fixed_params2,
                          solver = solver, bounds = bounds2, show_stats = show_stats)

    println("=== $name1 ===")
    println("Params: $(fit1.params), BIC: $(fit1.bic), SSR: $(fit1.ssr)")
    println("=== $name2 ===")
    println("Params: $(fit2.params), BIC: $(fit2.bic), SSR: $(fit2.ssr)")

    df_out = DataFrame(
        Model  = [name1, name2],
        Params = [string(fit1.params), string(fit2.params)],
        BIC    = [fit1.bic, fit2.bic],
        SSR    = [fit1.ssr, fit2.ssr],
    )
    CSV.write(output_csv, df_out)
    println("Results saved to $output_csv")

    best_model = fit1.bic <= fit2.bic ?
        (name = name1, params = fit1.params, bic = fit1.bic, ssr = fit1.ssr, solution = fit1.solution) :
        (name = name2, params = fit2.params, bic = fit2.bic, ssr = fit2.ssr, solution = fit2.solution)

    return (
        model1     = (name = name1, params = fit1.params, bic = fit1.bic, ssr = fit1.ssr, solution = fit1.solution),
        model2     = (name = name2, params = fit2.params, bic = fit2.bic, ssr = fit2.ssr, solution = fit2.solution),
        best_model = best_model,
    )
end

"""
    compare_datasets(x1, y1, name1, model1, p0_1, x2, y2, name2, model2, p0_2; ...)

Fit a model to two datasets and write a CSV summary.
"""
function compare_datasets(
    x1::Vector{<:Real}, y1::Vector{<:Real}, name1::String, model1::Function, p0_1::Vector{<:Real},
    x2::Vector{<:Real}, y2::Vector{<:Real}, name2::String, model2::Function, p0_2::Vector{<:Real};
    solver             = Rodas5(),
    bounds1            = nothing,
    bounds2            = nothing,
    fixed_params1      = nothing,
    fixed_params2      = nothing,
    show_stats::Bool   = false,
    output_csv::String = "dataset_comparison.csv",
)
    fit1 = run_single_fit(x1, y1, p0_1; model = model1, fixed_params = fixed_params1,
                          solver = solver, bounds = bounds1, show_stats = show_stats)
    fit2 = run_single_fit(x2, y2, p0_2; model = model2, fixed_params = fixed_params2,
                          solver = solver, bounds = bounds2, show_stats = show_stats)

    println("=== $name1 ===")
    println("Params: $(fit1.params), BIC: $(fit1.bic), SSR: $(fit1.ssr)")
    println("=== $name2 ===")
    println("Params: $(fit2.params), BIC: $(fit2.bic), SSR: $(fit2.ssr)")

    df_out = DataFrame(
        Dataset = [name1, name2],
        Params  = [string(fit1.params), string(fit2.params)],
        BIC     = [fit1.bic, fit2.bic],
        SSR     = [fit1.ssr, fit2.ssr],
    )
    CSV.write(output_csv, df_out)
    println("Results saved to $output_csv")
end

"""
    compare_models_dict(x, y, specs; default_solver=Rodas5(), show_stats=false, output_csv=\"all_models_comparison.csv\")

Fit each model in `specs` (a Dict of NamedTuples) to `x`, `y` and write summary tables.
"""
function compare_models_dict(
    x::Vector{<:Real},
    y::Vector{<:Real},
    specs::Dict{String,<:NamedTuple};
    default_solver        = Rodas5(),
    show_stats::Bool      = false,
    output_csv::String    = "all_models_comparison.csv",
)
    fits = Dict{String,Any}()
    results = NamedTuple[]

    for (name, spec) in specs
        solver_i = haskey(spec, :solver) ? spec.solver : default_solver
        fit = run_single_fit(
            x, y, spec.p0;
            model        = spec.model,
            fixed_params = get(spec, :fixed_params, nothing),
            solver       = solver_i,
            bounds       = get(spec, :bounds, nothing),
            show_stats   = show_stats,
        )
        fits[name] = fit
        push!(results, (Model = name, Params = fit.params, BIC = fit.bic, SSR = fit.ssr))
    end

    df_summary = DataFrame(
        Model  = [r.Model for r in results],
        Params = [string(r.Params) for r in results],
        BIC    = [r.BIC for r in results],
        SSR    = [r.SSR for r in results],
    )
    println("\nBIC Summary:")
    display(df_summary[:, [:Model, :BIC]])

    CSV.write(output_csv, df_summary)
    println("Summary saved to $output_csv")

    pred_rows = NamedTuple[]
    for (name, fit) in pairs(fits)
        for (t, u) in zip(fit.solution.t, fit.solution.u)
            push!(pred_rows, (Model = name, Time = t, Prediction = u[1]))
        end
    end
    df_preds = DataFrame(pred_rows)
    preds_csv = replace(output_csv, r"\\.csv$" => "_predictions.csv")
    CSV.write(preds_csv, df_preds)
    println("Predictions saved to $preds_csv")

    return fits
end

"""
    fit_three_datasets(x1, y1, name1, x2, y2, name2, x3, y3, name3, p0; ...)

Fit the same model to three datasets and return all results.
"""
function fit_three_datasets(
    x1::Vector{<:Real}, y1::Vector{<:Real}, name1::String,
    x2::Vector{<:Real}, y2::Vector{<:Real}, name2::String,
    x3::Vector{<:Real}, y3::Vector{<:Real}, name3::String,
    p0::Vector{<:Real};
    model                = Models.logistic_growth!,
    fixed_params         = nothing,
    solver               = Rodas5(),
    bounds               = nothing,
    show_stats::Bool     = false,
    output_csv::String   = "three_datasets_comparison.csv",
)
    fit1 = run_single_fit(x1, y1, p0; model = model, fixed_params = fixed_params,
                          solver = solver, bounds = bounds, show_stats = show_stats)
    fit2 = run_single_fit(x2, y2, p0; model = model, fixed_params = fixed_params,
                          solver = solver, bounds = bounds, show_stats = show_stats)
    fit3 = run_single_fit(x3, y3, p0; model = model, fixed_params = fixed_params,
                          solver = solver, bounds = bounds, show_stats = show_stats)

    println("=== $name1 ===")
    println("Params: $(fit1.params), BIC: $(fit1.bic), SSR: $(fit1.ssr)")
    println("=== $name2 ===")
    println("Params: $(fit2.params), BIC: $(fit2.bic), SSR: $(fit2.ssr)")
    println("=== $name3 ===")
    println("Params: $(fit3.params), BIC: $(fit3.bic), SSR: $(fit3.ssr)")

    df_out = DataFrame(
        Dataset = [name1, name2, name3],
        Params  = [string(fit1.params), string(fit2.params), string(fit3.params)],
        BIC     = [fit1.bic, fit2.bic, fit3.bic],
        SSR     = [fit1.ssr, fit2.ssr, fit3.ssr],
    )
    CSV.write(output_csv, df_out)
    println("Results saved to $output_csv")

    return (fit1 = fit1, fit2 = fit2, fit3 = fit3)
end

"""
    fit_three_datasets(x_datasets, y_datasets; p0=[0.1, 100.0], model=Models.logistic_growth!, fixed_params=nothing, solver=Rodas5(), bounds=nothing)

Convenience wrapper for fitting the same model to many datasets at once.
"""
function fit_three_datasets(
    x_datasets::Vector{<:Vector{<:Real}},
    y_datasets::Vector{<:Vector{<:Real}};
    p0::Vector{<:Real}     = [0.1, 100.0],
    model                  = Models.logistic_growth!,
    fixed_params           = nothing,
    solver                 = Rodas5(),
    bounds                 = nothing,
)
    n_datasets = length(x_datasets)
    @assert length(y_datasets) == n_datasets "Number of x and y datasets must match"

    individual_fits = Vector{Any}()
    for i in 1:n_datasets
        try
            fit_result = run_single_fit(
                x_datasets[i], y_datasets[i], p0;
                model        = model,
                fixed_params = fixed_params,
                solver       = solver,
                bounds       = bounds,
                show_stats   = false,
            )
            push!(individual_fits, (dataset = i, fit_result = fit_result))
        catch e
            println("Warning: Failed to fit dataset $i: $e")
            push!(individual_fits, (dataset = i, fit_result = nothing))
        end
    end

    successful_fits = filter(f -> f.fit_result !== nothing, individual_fits)

    if !isempty(successful_fits)
        all_params = [f.fit_result.params for f in successful_fits]
        mean_params = [mean([p[i] for p in all_params]) for i in 1:length(all_params[1])]
        std_params  = [std([p[i] for p in all_params]) for i in 1:length(all_params[1])]
        mean_ssr    = mean([f.fit_result.ssr for f in successful_fits])

        summary = (
            mean_params   = mean_params,
            std_params    = std_params,
            mean_ssr      = mean_ssr,
            n_successful  = length(successful_fits),
            n_total       = n_datasets,
        )
    else
        summary = (
            mean_params   = Float64[],
            std_params    = Float64[],
            mean_ssr      = Inf,
            n_successful  = 0,
            n_total       = n_datasets,
        )
    end

    return (individual_fits = individual_fits, summary = summary)
end

"""
    run_joint_fit(model, dataset_specs, u0, p0; solver=Tsit5(), bounds=nothing,
                  u0_builder=nothing, initial_time=nothing, show_stats=false,
                  maxiters=10_000)

Fit a single parameter vector `p` for a multi-state ODE model against multiple datasets.

`dataset_specs` is a vector of NamedTuples with fields:
- `x`: observation times
- `y`: observations
- either `state_index`: 1-based state index in the solution vector
- or `observable`: a function accepting `(u, p, t)`, `(u, p)`, or `u`
- optional `residual_scale`: a positive fixed scale applied to residuals for that
  dataset. The default `1.0` preserves ordinary SSE.

`u0_builder`, when supplied, maps the parameter vector to the initial state. This supports
joint mixture models whose fitted parameters include an initial population fraction.
`initial_time` sets the time associated with `u0` when it precedes the first
observation. This permits day-zero seeding states without adding synthetic
observations to the likelihood or BIC sample size.
`reltol` and `abstol` control ODE solve accuracy during optimization and prediction.
Set `optimizer=:nelder_mead` for bounded derivative-free screening of nonsmooth or
stiff models; the default remains `:bfgs` for backward compatibility.

Returns `(params, bic, sse, solution, predictions, save_times)`.
"""
function run_joint_fit(
    model::Function,
    dataset_specs::Vector{<:NamedTuple},
    u0::Vector{<:Real},
    p0::Vector{<:Real};
    solver            = Tsit5(),
    bounds            = nothing,
    u0_builder        = nothing,
    initial_time      = nothing,
    show_stats::Bool  = false,
    maxiters::Integer = 10_000,
    reltol::Real      = 1e-10,
    abstol::Real      = 1e-10,
    optimizer::Symbol = :bfgs,
)
    isempty(dataset_specs) && throw(ArgumentError("dataset_specs cannot be empty"))

    processed = [
        (
            x = Float64.(collect(ds.x)),
            y = Float64.(collect(ds.y)),
            state_index = haskey(ds, :state_index) ? Int(ds.state_index) : 0,
            observable = haskey(ds, :observable) ? ds.observable : nothing,
            residual_scale = haskey(ds, :residual_scale) ? max(Float64(ds.residual_scale), eps(Float64)) : 1.0,
        ) for ds in dataset_specs
    ]

    n_states = length(u0)
    for ds in processed
        if ds.observable === nothing
            ds.state_index >= 1 || throw(ArgumentError("each dataset requires state_index >= 1 or an observable"))
            ds.state_index <= n_states || throw(ArgumentError("state_index $(ds.state_index) exceeds u0 length $n_states"))
        end
        length(ds.x) == length(ds.y) || throw(ArgumentError("x and y length mismatch in dataset_specs"))
    end

    initial_state(p) = begin
        values = u0_builder === nothing ? u0 : u0_builder(p)
        length(values) == n_states || throw(ArgumentError("u0_builder must return $n_states states"))
        collect(values)
    end

    joint_observable(ds, state, p, t) = begin
        ds.observable === nothing && return state[ds.state_index]
        if applicable(ds.observable, state, p, t)
            return ds.observable(state, p, t)
        elseif applicable(ds.observable, state, p)
            return ds.observable(state, p)
        elseif applicable(ds.observable, state)
            return ds.observable(state)
        end
        throw(ArgumentError("observable must accept (u, p, t), (u, p), or u"))
    end

    save_times = sort(unique(vcat([ds.x for ds in processed]...)))
    fit_initial_time = initial_time === nothing ? save_times[1] : Float64(initial_time)
    fit_initial_time <= save_times[1] ||
        throw(ArgumentError("initial_time must be no later than the first observation"))
    tspan = (fit_initial_time, save_times[end])
    prob = ODEProblem(model, Float64.(initial_state(p0)), tspan, Float64.(p0))

    function objective(p)
        sol = try
            solve(remake(prob; p = p, u0 = initial_state(p)), solver; saveat = save_times, reltol = reltol, abstol = abstol)
        catch
            return 1e12
        end
        sol.retcode == ReturnCode.Success || return 1e12

        sse = 0.0
        for ds in processed
            for (xi, yi) in zip(ds.x, ds.y)
                idx = findfirst(t -> isapprox(t, xi; atol = 1e-10, rtol = 1e-10), save_times)
                idx === nothing && return 1e12
                yhat = joint_observable(ds, sol.u[idx], p, xi)
                isfinite(yhat) && yhat >= 0 || return 1e12
                sse += ((yi - yhat) / ds.residual_scale)^2
            end
        end
        return isfinite(sse) ? sse : 1e12
    end

    nparams = length(p0)
    _default_upper(p::Real) = max(10.0, abs(Float64(p)) * 100)

    if bounds === nothing
        bounds = [(0.0, _default_upper(p0[i])) for i in 1:nparams]
    else
        length(bounds) == nparams || throw(ArgumentError("bounds must have length $nparams"))
        bounds = [
            (
                isfinite(b[1]) ? Float64(b[1]) : 0.0,
                isfinite(b[2]) ? Float64(b[2]) : _default_upper(p0[i]),
            ) for (i, b) in enumerate(bounds)
        ]
    end

    lb = Float64[first(b) for b in bounds]
    ub = Float64[last(b) for b in bounds]
    p_opt = if optimizer == :bfgs
        loss = Optimization.OptimizationFunction((x, _) -> objective(x), Optimization.AutoForwardDiff())
        optprob = Optimization.OptimizationProblem(loss, Float64.(p0); lb = lb, ub = ub)
        result = Optimization.solve(optprob, OptimizationOptimJL.Fminbox(OptimizationOptimJL.BFGS()); maxiters = maxiters)
        collect(result.u)
    elseif optimizer == :nelder_mead
        bounded(z) = lb .+ (ub .- lb) ./ (1 .+ exp.(-clamp.(z, -30.0, 30.0)))
        fractions = clamp.((Float64.(p0) .- lb) ./ (ub .- lb), 1e-6, 1 - 1e-6)
        z0 = log.(fractions ./ (1 .- fractions))
        loss = Optimization.OptimizationFunction((z, _) -> objective(bounded(z)))
        optprob = Optimization.OptimizationProblem(loss, z0)
        result = Optimization.solve(optprob, OptimizationOptimJL.NelderMead(); maxiters = maxiters)
        collect(bounded(result.u))
    else
        throw(ArgumentError("optimizer must be :bfgs or :nelder_mead"))
    end
    sse = objective(p_opt)
    n = sum(length(ds.x) for ds in processed)
    k = length(p_opt)
    bic = n * log(max(sse, 1e-12) / n) + k * log(n)

    sol_opt = solve(remake(prob; p = p_opt, u0 = initial_state(p_opt)), solver; saveat = save_times, reltol = reltol, abstol = abstol)
    predictions = [
        [joint_observable(ds, sol_opt.u[findfirst(t -> isapprox(t, xi; atol = 1e-10, rtol = 1e-10), save_times)], p_opt, xi)
         for xi in ds.x]
        for ds in processed
    ]
    raw_sse = sum(
        sum((yi - yhat)^2 for (yi, yhat) in zip(ds.y, prediction))
        for (ds, prediction) in zip(processed, predictions)
    )

    if show_stats
        println("Optimized params: ", p_opt)
        println("Joint SSE: ", sse)
        println("Joint BIC: ", bic)
    end

    return (
        params = p_opt,
        bic = bic,
        sse = sse,
        scaled_sse = sse,
        raw_sse = raw_sse,
        solution = sol_opt,
        predictions = predictions,
        save_times = save_times,
        initial_time = fit_initial_time,
    )
end

"""
    run_joint_multistart(model, dataset_specs, u0, starts; kwargs...)

Run `run_joint_fit` from several initial parameter vectors and retain the
finite fit with the lowest BIC. All starts must have the same parameter count.
The returned `summary` records successful and failed starts for auditability.
"""
function run_joint_multistart(
    model::Function,
    dataset_specs::Vector{<:NamedTuple},
    u0::Vector{<:Real},
    starts::Vector{<:AbstractVector};
    kwargs...,
)
    isempty(starts) && throw(ArgumentError("starts cannot be empty"))
    parameter_count = length(first(starts))
    all(length(start) == parameter_count for start in starts) ||
        throw(ArgumentError("all multistart vectors must have equal length"))

    fits = Any[]
    successful_indices = Int[]
    rows = NamedTuple[]
    for (start_index, start) in enumerate(starts)
        fit = try
            run_joint_fit(model, dataset_specs, u0, Float64.(start); kwargs...)
        catch error
            push!(fits, nothing)
            push!(rows, (
                start_index = start_index,
                status = "failed",
                bic = NaN,
                scaled_sse = NaN,
                raw_sse = NaN,
                params = "",
                message = sprint(showerror, error),
            ))
            continue
        end
        valid = isfinite(fit.bic) && isfinite(fit.scaled_sse) &&
            isfinite(fit.raw_sse) && fit.scaled_sse < 9.99e11
        push!(fits, fit)
        valid && push!(successful_indices, start_index)
        push!(rows, (
            start_index = start_index,
            status = valid ? "completed" : "invalid",
            bic = valid ? Float64(fit.bic) : NaN,
            scaled_sse = valid ? Float64(fit.scaled_sse) : NaN,
            raw_sse = valid ? Float64(fit.raw_sse) : NaN,
            params = valid ? string(Float64.(fit.params)) : "",
            message = valid ? "" : "Non-finite or failure-sentinel joint fit.",
        ))
    end
    isempty(successful_indices) && error("No finite multistart fit completed")
    best_start_index = successful_indices[
        argmin([Float64(fits[index].bic) for index in successful_indices])
    ]
    return (
        fit = fits[best_start_index],
        fits = fits,
        summary = DataFrame(rows),
        best_start_index = best_start_index,
    )
end

"""
    profile_joint_fit_bounds(model, dataset_specs, u0, p0;
        bounds, parameter_names, explicit_upper_profiles=Dict(),
        trigger_fraction=0.02, expansion_factors=[1.5, 2.0],
        min_bic_improvement=2.0, min_new_bound_margin=0.05, kwargs...)

Run a joint fit and sequentially profile upper bounds for parameters whose
estimate lies within `trigger_fraction` of the active bound. An expanded fit is
accepted only when it improves BIC by at least `min_bic_improvement` and moves
the estimate at least `min_new_bound_margin` inside the expanded interval.

Returns `(fit, bounds, profile, identifiability)`. `profile` records every
attempt, and `identifiability` marks parameters that remain near a bound.
"""
function profile_joint_fit_bounds(
    model::Function,
    dataset_specs::Vector{<:NamedTuple},
    u0::Vector{<:Real},
    p0::Vector{<:Real};
    bounds::Vector{<:Tuple},
    parameter_names::Vector{Symbol},
    explicit_upper_profiles::AbstractDict = Dict{Symbol,Vector{Float64}}(),
    profile_parameters::Union{Nothing,Vector{Symbol}} = nothing,
    trigger_fraction::Real = 0.02,
    expansion_factors::Vector{Float64} = [1.5, 2.0],
    min_bic_improvement::Real = 2.0,
    min_new_bound_margin::Real = 0.05,
    kwargs...,
)
    length(bounds) == length(p0) == length(parameter_names) ||
        throw(ArgumentError("bounds, p0, and parameter_names must have equal length"))
    0 <= trigger_fraction < 1 || throw(ArgumentError("trigger_fraction must be in [0, 1)"))
    0 <= min_new_bound_margin < 1 || throw(ArgumentError("min_new_bound_margin must be in [0, 1)"))

    active_bounds = [(Float64(lo), Float64(hi)) for (lo, hi) in bounds]
    best = run_joint_fit(model, dataset_specs, u0, Float64.(p0); bounds = active_bounds, kwargs...)
    rows = NamedTuple[]

    for index in eachindex(parameter_names)
        name = parameter_names[index]
        profile_parameters !== nothing && !(name in profile_parameters) && continue
        lower, upper = active_bounds[index]
        span = upper - lower
        span > 0 || continue
        base_margin = (upper - Float64(best.params[index])) / span
        base_margin <= trigger_fraction || continue

        candidates = haskey(explicit_upper_profiles, name) ?
            sort(unique(Float64.(explicit_upper_profiles[name]))) :
            sort(unique(upper .* expansion_factors))
        candidates = filter(>(upper), candidates)

        accepted = false
        for candidate_upper in candidates
            candidate_bounds = copy(active_bounds)
            candidate_bounds[index] = (lower, candidate_upper)
            candidate_p0 = clamp.(Float64.(best.params), first.(candidate_bounds), last.(candidate_bounds))
            candidate = run_joint_fit(
                model,
                dataset_specs,
                u0,
                candidate_p0;
                bounds = candidate_bounds,
                kwargs...,
            )
            candidate_span = candidate_upper - lower
            candidate_margin = (candidate_upper - Float64(candidate.params[index])) / candidate_span
            bic_improvement = Float64(best.bic) - Float64(candidate.bic)
            accept = bic_improvement >= min_bic_improvement && candidate_margin >= min_new_bound_margin
            push!(rows, (
                parameter = String(name),
                original_upper = upper,
                candidate_upper = candidate_upper,
                original_bic = Float64(best.bic),
                candidate_bic = Float64(candidate.bic),
                bic_improvement = bic_improvement,
                estimate = Float64(candidate.params[index]),
                upper_margin_fraction = candidate_margin,
                accepted = accept,
            ))
            if accept
                best = candidate
                active_bounds = candidate_bounds
                accepted = true
                break
            end
        end

        if !accepted && isempty(candidates)
            push!(rows, (
                parameter = String(name),
                original_upper = upper,
                candidate_upper = upper,
                original_bic = Float64(best.bic),
                candidate_bic = Float64(best.bic),
                bic_improvement = 0.0,
                estimate = Float64(best.params[index]),
                upper_margin_fraction = base_margin,
                accepted = false,
            ))
        end
    end

    identifiability = DataFrame(
        parameter = String.(parameter_names),
        estimate = Float64.(best.params),
        lower_bound = first.(active_bounds),
        upper_bound = last.(active_bounds),
    )
    identifiability.bound_position = [
        (identifiability.estimate[i] - identifiability.lower_bound[i]) /
        max(identifiability.upper_bound[i] - identifiability.lower_bound[i], eps(Float64))
        for i in 1:nrow(identifiability)
    ]
    identifiability.identifiability = [
        position <= trigger_fraction || position >= 1 - trigger_fraction ?
            "poorly_identified_at_bound" : "interior"
        for position in identifiability.bound_position
    ]

    return (
        fit = best,
        bounds = active_bounds,
        profile = DataFrame(rows),
        identifiability = identifiability,
    )
end

"""
    profile_joint_fit_bounds_two_sided(model, dataset_specs, u0, p0;
        bounds, parameter_names, explicit_lower_profiles=Dict(), kwargs...)

Run the standard upper-bound profile and then test explicitly supplied lower
bounds. Lower expansion is deliberately explicit because many biological
parameters have a physical lower limit of zero. A candidate is accepted only
when it improves BIC and moves the estimate away from the expanded boundary.
"""
function profile_joint_fit_bounds_two_sided(
    model::Function,
    dataset_specs::Vector{<:NamedTuple},
    u0::Vector{<:Real},
    p0::Vector{<:Real};
    bounds::Vector{<:Tuple},
    parameter_names::Vector{Symbol},
    explicit_upper_profiles::AbstractDict = Dict{Symbol,Vector{Float64}}(),
    explicit_lower_profiles::AbstractDict = Dict{Symbol,Vector{Float64}}(),
    physical_lower_limits::AbstractDict = Dict{Symbol,Float64}(),
    profile_parameters::Union{Nothing,Vector{Symbol}} = nothing,
    trigger_fraction::Real = 0.02,
    min_bic_improvement::Real = 2.0,
    min_new_bound_margin::Real = 0.05,
    kwargs...,
)
    upper = profile_joint_fit_bounds(
        model,
        dataset_specs,
        u0,
        p0;
        bounds = bounds,
        parameter_names = parameter_names,
        explicit_upper_profiles = explicit_upper_profiles,
        profile_parameters = profile_parameters,
        trigger_fraction = trigger_fraction,
        min_bic_improvement = min_bic_improvement,
        min_new_bound_margin = min_new_bound_margin,
        kwargs...,
    )
    best = upper.fit
    active_bounds = copy(upper.bounds)
    rows = NamedTuple[]
    for row in eachrow(upper.profile)
        push!(rows, (
            parameter = String(row.parameter),
            direction = "upper",
            original_bound = Float64(row.original_upper),
            candidate_bound = Float64(row.candidate_upper),
            original_bic = Float64(row.original_bic),
            candidate_bic = Float64(row.candidate_bic),
            bic_improvement = Float64(row.bic_improvement),
            estimate = Float64(row.estimate),
            bound_margin_fraction = Float64(row.upper_margin_fraction),
            accepted = Bool(row.accepted),
        ))
    end

    for index in eachindex(parameter_names)
        name = parameter_names[index]
        profile_parameters !== nothing && !(name in profile_parameters) && continue
        haskey(explicit_lower_profiles, name) || continue
        lower, upper_bound = active_bounds[index]
        span = upper_bound - lower
        span > 0 || continue
        lower_margin = (Float64(best.params[index]) - lower) / span
        lower_margin <= trigger_fraction || continue
        physical_lower = Float64(get(physical_lower_limits, name, -Inf))
        candidates = sort(unique(Float64.(explicit_lower_profiles[name])); rev = true)
        candidates = filter(candidate -> physical_lower <= candidate < lower, candidates)
        for candidate_lower in candidates
            candidate_bounds = copy(active_bounds)
            candidate_bounds[index] = (candidate_lower, upper_bound)
            candidate_p0 = clamp.(Float64.(best.params), first.(candidate_bounds), last.(candidate_bounds))
            candidate = run_joint_fit(
                model,
                dataset_specs,
                u0,
                candidate_p0;
                bounds = candidate_bounds,
                kwargs...,
            )
            candidate_span = upper_bound - candidate_lower
            candidate_margin = (Float64(candidate.params[index]) - candidate_lower) / candidate_span
            bic_improvement = Float64(best.bic) - Float64(candidate.bic)
            accept = bic_improvement >= min_bic_improvement && candidate_margin >= min_new_bound_margin
            push!(rows, (
                parameter = String(name),
                direction = "lower",
                original_bound = lower,
                candidate_bound = candidate_lower,
                original_bic = Float64(best.bic),
                candidate_bic = Float64(candidate.bic),
                bic_improvement = bic_improvement,
                estimate = Float64(candidate.params[index]),
                bound_margin_fraction = candidate_margin,
                accepted = accept,
            ))
            if accept
                best = candidate
                active_bounds = candidate_bounds
                break
            end
        end
    end

    identifiability = DataFrame(
        parameter = String.(parameter_names),
        estimate = Float64.(best.params),
        lower_bound = first.(active_bounds),
        upper_bound = last.(active_bounds),
    )
    identifiability.bound_position = [
        (identifiability.estimate[i] - identifiability.lower_bound[i]) /
        max(identifiability.upper_bound[i] - identifiability.lower_bound[i], eps(Float64))
        for i in 1:nrow(identifiability)
    ]
    identifiability.identifiability = [
        position <= trigger_fraction || position >= 1 - trigger_fraction ?
            "poorly_identified_at_bound" : "interior"
        for position in identifiability.bound_position
    ]
    return (
        fit = best,
        bounds = active_bounds,
        profile = DataFrame(rows),
        identifiability = identifiability,
    )
end

"""
    compare_joint_models_dict(dataset_specs, u0, specs; default_solver=Tsit5(), show_stats=false, output_csv="joint_model_comparison.csv")

Fit multiple joint models and write a BIC/SSE summary CSV.

`specs` maps model names to NamedTuples like `(model=<f!>, p0=<vector>, bounds=<vector>)`.
A spec may also override `u0` or provide `u0_builder` for parameterized initial states.
"""
function compare_joint_models_dict(
    dataset_specs::Vector{<:NamedTuple},
    u0::Vector{<:Real},
    specs::Dict{String,<:NamedTuple};
    default_solver      = Tsit5(),
    show_stats::Bool    = false,
    output_csv::Union{Nothing,String} = "joint_model_comparison.csv",
    maxiters::Integer = 10_000,
)
    fits = Dict{String,Any}()
    rows = NamedTuple[]

    for (name, spec) in specs
        solver_i = haskey(spec, :solver) ? spec.solver : default_solver
        u0_i = haskey(spec, :u0) ? spec.u0 : u0
        u0_builder_i = haskey(spec, :u0_builder) ? spec.u0_builder : nothing
        fit = run_joint_fit(
            spec.model,
            dataset_specs,
            u0_i,
            spec.p0;
            solver = solver_i,
            bounds = get(spec, :bounds, nothing),
            u0_builder = u0_builder_i,
            show_stats = show_stats,
            maxiters = maxiters,
        )
        fits[name] = fit
        push!(rows, (Model = name, Params = fit.params, BIC = fit.bic, SSE = fit.sse))
    end

    df_summary = DataFrame(
        Model = [r.Model for r in rows],
        Params = [string(r.Params) for r in rows],
        BIC = [r.BIC for r in rows],
        SSE = [r.SSE for r in rows],
    )

    sort!(df_summary, :BIC)
    if output_csv !== nothing
        CSV.write(output_csv, df_summary)
        println("Joint comparison summary saved to $output_csv")
    end

    return fits
end

end # module Fitting
