# Analysis module - Contains all statistical analysis and validation functions
module Analysis

using StatsBase
using DifferentialEquations
using Distributions
using Random

# Import models and fitting from other modules
using ..Models
using ..Fitting

export leave_one_out_validation, k_fold_cross_validation, parameter_sensitivity_analysis,
       residual_analysis, enhanced_bic_analysis

"""
leave_one_out_validation(
    x::Vector{<:Real},
    y::Vector{<:Real},
    p0::Vector{<:Real};
    model                = Models.logistic_growth!,
    fixed_params         = nothing,
    solver               = Rodas5(),
    bounds               = nothing,
    show_stats::Bool     = false
)

Performs leave-one-out cross-validation by fitting the model with each data point 
removed and predicting that point. Returns validation metrics including RMSE and R².
"""
function leave_one_out_validation(
    x::Vector{<:Real},
    y::Vector{<:Real},
    p0::Vector{<:Real};
    model                = Models.logistic_growth!,
    fixed_params         = nothing,
    solver               = Rodas5(),
    bounds               = nothing,
    show_stats::Bool     = false
)
    n = length(x)
    predictions = zeros(n)
    fit_params = Vector{Vector{Float64}}(undef, n)
    
    println("Performing leave-one-out cross-validation...")
    
    for i in 1:n
        # Create training set (all points except i)
        x_train = [x[j] for j in 1:n if j != i]
        y_train = [y[j] for j in 1:n if j != i]
        
        try
            # Fit model on training data
            fit_result = Fitting.run_single_fit(
                x_train, y_train, p0;
                model        = model,
                fixed_params = fixed_params,
                solver       = solver,
                bounds       = bounds,
                show_stats   = false
            )
            
            fit_params[i] = fit_result.params
            
            # Predict the left-out point
            tspan = (minimum(x_train), x[i])
            if x[i] < minimum(x_train)
                tspan = (x[i], maximum(x_train))
            elseif x[i] > maximum(x_train)
                tspan = (minimum(x_train), x[i])
            end
            
            prob = ODEProblem(model, [y_train[1]], tspan, fit_result.params)
            sol = solve(prob, solver, saveat=[x[i]])
            
            if length(sol.u) > 0
                predictions[i] = sol.u[end][1]
            else
                predictions[i] = NaN
            end
            
        catch e
            println("Warning: Failed to fit model excluding point $i: $e")
            predictions[i] = NaN
            fit_params[i] = fill(NaN, length(p0))
        end
        
        if show_stats
            println("LOO iteration $i/$n completed")
        end
    end
    
    # Calculate validation metrics
    valid_indices = .!isnan.(predictions)
    y_valid = y[valid_indices]
    pred_valid = predictions[valid_indices]
    
    rmse = sqrt(mean((y_valid .- pred_valid).^2))
    mae = mean(abs.(y_valid .- pred_valid))
    
    # R-squared
    ss_res = sum((y_valid .- pred_valid).^2)
    ss_tot = sum((y_valid .- mean(y_valid)).^2)
    r_squared = 1 - ss_res/ss_tot
    
    # Parameter variability
    param_std = [std([fit_params[i][j] for i in 1:n if !any(isnan, fit_params[i])]) 
                 for j in 1:length(p0)]
    
    results = (
        predictions = predictions,
        actual = y,
        rmse = rmse,
        mae = mae,
        r_squared = r_squared,
        fit_params = fit_params,
        param_std = param_std,
        n_valid = sum(valid_indices)
    )
    
    if show_stats
        println("\n=== Leave-One-Out Cross-Validation Results ===")
        println("RMSE: $(round(rmse, digits=4))")
        println("MAE: $(round(mae, digits=4))")
        println("R²: $(round(r_squared, digits=4))")
        println("Valid predictions: $(sum(valid_indices))/$n")
        println("Parameter standard deviations: $(round.(param_std, digits=4))")
    end
    
    return results
end

"""
k_fold_cross_validation(
    x::Vector{<:Real},
    y::Vector{<:Real},
    p0::Vector{<:Real};
    k_folds::Int         = 5,
    model                = Models.logistic_growth!,
    fixed_params         = nothing,
    solver               = Rodas5(),
    bounds               = nothing,
    show_stats::Bool     = false
)

Performs k-fold cross-validation by splitting data into k folds and using each 
fold as validation set while training on the remaining folds.
"""
function k_fold_cross_validation(
    x::Vector{<:Real},
    y::Vector{<:Real},
    p0::Vector{<:Real};
    k_folds::Int         = 5,
    model                = Models.logistic_growth!,
    fixed_params         = nothing,
    solver               = Rodas5(),
    bounds               = nothing,
    show_stats::Bool     = false
)
    n = length(x)
    fold_size = div(n, k_folds)
    
    # Randomly shuffle indices
    indices = Random.randperm(n)
    
    all_predictions = Float64[]
    all_actual = Float64[]
    fold_metrics = NamedTuple[]
    
    println("Performing $k_folds-fold cross-validation...")
    
    for k in 1:k_folds
        # Define validation fold
        start_idx = (k-1) * fold_size + 1
        end_idx = k == k_folds ? n : k * fold_size
        val_indices = indices[start_idx:end_idx]
        train_indices = setdiff(indices, val_indices)
        
        # Split data
        x_train, y_train = x[train_indices], y[train_indices]
        x_val, y_val = x[val_indices], y[val_indices]

        # Ensure monotonic time order for ODE fitting/prediction
        train_order = sortperm(x_train)
        x_train = x_train[train_order]
        y_train = y_train[train_order]

        val_order = sortperm(x_val)
        x_val = x_val[val_order]
        y_val = y_val[val_order]
        
        try
            # Fit on training data
            fit_result = Fitting.run_single_fit(
                x_train, y_train, p0;
                model        = model,
                fixed_params = fixed_params,
                solver       = solver,
                bounds       = bounds,
                show_stats   = false
            )
            
            # Predict validation set
            tspan = (minimum(x_train), maximum([maximum(x_train), maximum(x_val)]))
            prob = ODEProblem(model, [y_train[1]], tspan, fit_result.params)
            
            predictions = Float64[]
            for xi in x_val
                sol = solve(prob, solver, saveat=[xi])
                if length(sol.u) > 0
                    push!(predictions, sol.u[end][1])
                else
                    push!(predictions, NaN)
                end
            end
            
            # Calculate fold metrics
            valid_mask = .!isnan.(predictions)
            if sum(valid_mask) > 0
                pred_valid = predictions[valid_mask]
                val_valid = y_val[valid_mask]
                
                rmse = sqrt(mean((val_valid .- pred_valid).^2))
                mae = mean(abs.(val_valid .- pred_valid))
                
                push!(fold_metrics, (fold=k, rmse=rmse, mae=mae, n_valid=sum(valid_mask)))
                append!(all_predictions, pred_valid)
                append!(all_actual, val_valid)
            end
            
        catch e
            println("Warning: Failed to fit fold $k: $e")
            push!(fold_metrics, (fold=k, rmse=NaN, mae=NaN, n_valid=0))
        end
        
        if show_stats
            println("Fold $k/$k_folds completed")
        end
    end
    
    # Overall metrics
    overall_rmse = sqrt(mean((all_actual .- all_predictions).^2))
    overall_mae = mean(abs.(all_actual .- all_predictions))
    
    ss_res = sum((all_actual .- all_predictions).^2)
    ss_tot = sum((all_actual .- mean(all_actual)).^2)
    r_squared = 1 - ss_res/ss_tot
    
    results = (
        fold_metrics = fold_metrics,
        overall_rmse = overall_rmse,
        overall_mae = overall_mae,
        r_squared = r_squared,
        predictions = all_predictions,
        actual = all_actual
    )
    
    if show_stats
        println("\n=== $k_folds-Fold Cross-Validation Results ===")
        println("Overall RMSE: $(round(overall_rmse, digits=4))")
        println("Overall MAE: $(round(overall_mae, digits=4))")
        println("Overall R²: $(round(r_squared, digits=4))")
        
        fold_rmse = [m.rmse for m in fold_metrics if !isnan(m.rmse)]
        if length(fold_rmse) > 0
            println("Fold RMSE: mean=$(round(mean(fold_rmse), digits=4)), std=$(round(std(fold_rmse), digits=4))")
        end
    end
    
    return results
end

"""
parameter_sensitivity_analysis(
    x::Vector{<:Real},
    y::Vector{<:Real},
    fit_result::NamedTuple;
    perturbation::Float64    = 0.1,
    model                   = Models.logistic_growth!,
    solver                  = Rodas5()
)

Analyzes how sensitive model predictions are to changes in fitted parameters.
Perturbs each parameter by ±perturbation fraction and measures the effect on predictions.
"""
function parameter_sensitivity_analysis(
    x::Vector{<:Real},
    y::Vector{<:Real},
    fit_result::NamedTuple;
    perturbation::Float64    = 0.1,
    model                   = Models.logistic_growth!,
    solver                  = Rodas5()
)
    params = fit_result.params
    n_params = length(params)
    
    # Get baseline predictions
    tspan = (x[1], x[end])
    x_dense = range(x[1], x[end], length=100)
    prob_baseline = ODEProblem(model, [y[1]], tspan, params)
    sol_baseline = solve(prob_baseline, solver, saveat=x_dense)
    baseline_pred = getindex.(sol_baseline.u, 1)
    
    # Storage for sensitivity results
    sensitivity_metrics = Dict{Int, NamedTuple}()
    
    println("Performing parameter sensitivity analysis...")
    println("Base parameters: $(round.(params, digits=4))")
    println("Perturbation: ±$(perturbation*100)%")
    
    for i in 1:n_params
        # Perturb parameter up and down
        params_up = copy(params)
        params_down = copy(params)
        
        params_up[i] *= (1 + perturbation)
        params_down[i] *= (1 - perturbation)
        
        try
            # Solve with perturbed parameters
            prob_up = ODEProblem(model, [y[1]], tspan, params_up)
            prob_down = ODEProblem(model, [y[1]], tspan, params_down)
            
            sol_up = solve(prob_up, solver, saveat=x_dense)
            sol_down = solve(prob_down, solver, saveat=x_dense)
            
            pred_up = getindex.(sol_up.u, 1)
            pred_down = getindex.(sol_down.u, 1)
            
            # Calculate sensitivity metrics
            # Absolute change
            abs_change_up = abs.(pred_up .- baseline_pred)
            abs_change_down = abs.(pred_down .- baseline_pred)
            
            # Relative change
            rel_change_up = abs_change_up ./ abs.(baseline_pred .+ 1e-10)
            rel_change_down = abs_change_down ./ abs.(baseline_pred .+ 1e-10)
            
            # Summary metrics
            max_abs_change = max(maximum(abs_change_up), maximum(abs_change_down))
            mean_abs_change = mean([abs_change_up; abs_change_down])
            max_rel_change = max(maximum(rel_change_up), maximum(rel_change_down))
            mean_rel_change = mean([rel_change_up; rel_change_down])
            
            # Sensitivity index (normalized by parameter change)
            param_rel_change = perturbation
            sensitivity_index = max_rel_change / param_rel_change
            
            sensitivity_metrics[i] = (
                param_index = i,
                param_value = params[i],
                max_abs_change = max_abs_change,
                mean_abs_change = mean_abs_change,
                max_rel_change = max_rel_change,
                mean_rel_change = mean_rel_change,
                sensitivity_index = sensitivity_index,
                pred_up = pred_up,
                pred_down = pred_down
            )
            
            println("Parameter $i: SI = $(round(sensitivity_index, digits=3)), Max rel change = $(round(max_rel_change*100, digits=2))%")
            
        catch e
            println("Warning: Failed sensitivity analysis for parameter $i: $e")
            sensitivity_metrics[i] = (
                param_index = i,
                param_value = params[i],
                max_abs_change = NaN,
                mean_abs_change = NaN,
                max_rel_change = NaN,
                mean_rel_change = NaN,
                sensitivity_index = NaN,
                pred_up = fill(NaN, length(x_dense)),
                pred_down = fill(NaN, length(x_dense))
            )
        end
    end
    
    # Rank parameters by sensitivity
    valid_metrics = [m for m in values(sensitivity_metrics) if !isnan(m.sensitivity_index)]
    sorted_metrics = sort(valid_metrics, by = m -> m.sensitivity_index, rev=true)
    
    println("\n=== Parameter Sensitivity Ranking ===")
    for (rank, metric) in enumerate(sorted_metrics)
        println("$rank. Parameter $(metric.param_index) (value=$(round(metric.param_value, digits=4))): SI=$(round(metric.sensitivity_index, digits=3))")
    end

    return (
        sensitivity_metrics = sensitivity_metrics,
        ranking = sorted_metrics,
        baseline_predictions = baseline_pred,
        x_dense = x_dense
    )
end

"""
residual_analysis(
    x::Vector{<:Real},
    y::Vector{<:Real},
    fit_result::NamedTuple;
    model              = Models.logistic_growth!,
    solver             = Rodas5(),
    outlier_threshold::Float64 = 2.0
)

Performs comprehensive residual analysis for model diagnostics.
Calculates residuals and identifies outliers.
"""
function residual_analysis(
    x::Vector{<:Real},
    y::Vector{<:Real},
    fit_result::NamedTuple;
    model              = Models.logistic_growth!,
    solver             = Rodas5(),
    outlier_threshold::Float64 = 2.0
)
    params = fit_result.params
    
    # Get model predictions at data points
    tspan = (x[1], x[end])
    prob = ODEProblem(model, [y[1]], tspan, params)
    sol = solve(prob, solver, saveat=x)
    y_pred = getindex.(sol.u, 1)
    
    # Calculate residuals
    residuals = y .- y_pred
    
    # Calculate standardized residuals
    residual_std = std(residuals)
    standardized_residuals = residuals ./ residual_std
    
    # Identify outliers
    outlier_indices = findall(abs.(standardized_residuals) .> outlier_threshold)
    n_outliers = length(outlier_indices)
    
    # Calculate residual statistics
    residual_stats = (
        mean_residual = mean(residuals),
        std_residual = residual_std,
        rmse = sqrt(mean(residuals.^2)),
        mae = mean(abs.(residuals)),
        max_abs_residual = maximum(abs.(residuals)),
        min_residual = minimum(residuals),
        max_residual = maximum(residuals),
        n_outliers = n_outliers,
        outlier_threshold = outlier_threshold
    )
    
    # Test for normality (Shapiro-Wilk approximation)
    n = length(residuals)
    if n >= 3
        sorted_residuals = sort(standardized_residuals)
        # Simple normality test: check if residuals follow normal distribution
        expected_quantiles = [quantile(Normal(0,1), (i-0.5)/n) for i in 1:n]
        correlation_coef = cor(sorted_residuals, expected_quantiles)
        normality_pvalue = 1 - abs(1 - correlation_coef)  # Approximate p-value
    else
        correlation_coef = NaN
        normality_pvalue = NaN
    end
    
    # Test for autocorrelation (Durbin-Watson approximation)
    if n > 1
        diff_residuals = diff(residuals)
        durbin_watson = sum(diff_residuals.^2) / sum(residuals.^2)
        # DW ≈ 2 indicates no autocorrelation, < 2 positive autocorr, > 2 negative
        autocorr_concern = (durbin_watson < 1.5 || durbin_watson > 2.5)
    else
        durbin_watson = NaN
        autocorr_concern = false
    end
    
    # Create diagnostic summary
    println("=== Residual Analysis Summary ===")
    println("RMSE: $(round(residual_stats.rmse, digits=4))")
    println("MAE: $(round(residual_stats.mae, digits=4))")
    println("Mean residual: $(round(residual_stats.mean_residual, digits=4))")
    println("Std residual: $(round(residual_stats.std_residual, digits=4))")
    println("Max |residual|: $(round(residual_stats.max_abs_residual, digits=4))")
    println("Outliers (|std_resid| > $outlier_threshold): $n_outliers")
    
    if !isnan(correlation_coef)
        println("Normality correlation: $(round(correlation_coef, digits=3))")
        println("Normality concern: $(correlation_coef < 0.95 ? "Yes" : "No")")
    end
    
    if !isnan(durbin_watson)
        println("Durbin-Watson statistic: $(round(durbin_watson, digits=3))")
        println("Autocorrelation concern: $(autocorr_concern ? "Yes" : "No")")
    end
    
    if n_outliers > 0
        println("\nOutlier details:")
        for idx in outlier_indices
            println("  Point $idx: x=$(x[idx]), y=$(y[idx]), residual=$(round(residuals[idx], digits=3)), std_resid=$(round(standardized_residuals[idx], digits=3))")
        end
    end
    
    return (
        residuals = residuals,
        standardized_residuals = standardized_residuals,
        predicted_values = y_pred,
        outlier_indices = outlier_indices,
        statistics = residual_stats,
        normality_correlation = correlation_coef,
        durbin_watson = durbin_watson,
        autocorrelation_concern = autocorr_concern
    )
end

"""
enhanced_bic_analysis(
    x::Vector{<:Real},
    y::Vector{<:Real};
    models = [Models.logistic_growth!, Models.gompertz_growth!, Models.exponential_growth_with_delay!],
    model_names = ["Logistic", "Gompertz", "Exponential+Delay"],
    p0_values = [[0.1, 100.0], [0.1, 100.0], [0.1, 1.0, 1.0]],
    solver = Rodas5(),
    population_size::Int = 150,
    max_time::Float64 = 60.0
)

Performs comprehensive BIC-based model comparison and selection.
Fits multiple models, calculates information criteria, and provides model recommendations.
"""
function enhanced_bic_analysis(
    x::Vector{<:Real},
    y::Vector{<:Real};
    models = [Models.logistic_growth!, Models.gompertz_growth!, Models.exponential_growth_with_delay!],
    model_names = ["Logistic", "Gompertz", "Exponential+Delay"],
    p0_values = [[0.1, 100.0], [0.1, 100.0], [0.1, 1.0, 1.0]],
    solver = Rodas5(),
    population_size::Int = 150,
    max_time::Float64 = 60.0
)
    n_models = length(models)
    results = []
    
    println("=== Enhanced BIC Analysis ===")
    println("Comparing $n_models models...")
    
    for i in 1:n_models
        model = models[i]
        model_name = i <= length(model_names) ? model_names[i] : "Model_$i"
        p0 = i <= length(p0_values) ? p0_values[i] : [0.1, 100.0]
        
        println("\nFitting $model_name...")
        
        try
            # Fit the model
            fit_result = Fitting.run_single_fit(x, y, p0; 
                                        model=model, 
                                        solver=solver,
                                        show_stats=false)
            
            # Calculate additional information criteria
            n = length(y)
            k = length(fit_result.params)
            
            # AIC (Akaike Information Criterion)
            aic = n * log(fit_result.ssr / n) + 2 * k
            
            # AICc (corrected AIC for small samples)
            if n > k + 1
                aicc = aic + (2 * k * (k + 1)) / (n - k - 1)
            else
                aicc = Inf  # AICc undefined for small samples
            end
            
            # BIC (already calculated in run_single_fit)
            bic = fit_result.bic
            
            # Calculate R-squared
            ss_res = fit_result.ssr
            ss_tot = sum((y .- mean(y)).^2)
            r_squared = 1 - ss_res / ss_tot
            
            # Calculate adjusted R-squared
            adj_r_squared = 1 - (ss_res / (n - k)) / (ss_tot / (n - 1))
            
            # Model complexity penalty (normalized)
            complexity_penalty = k / n
            
            # Calculate prediction accuracy at data points
            tspan = (x[1], x[end])
            prob = ODEProblem(model, [y[1]], tspan, fit_result.params)
            sol = solve(prob, solver, saveat=x)
            y_pred = getindex.(sol.u, 1)
            
            rmse = sqrt(mean((y .- y_pred).^2))
            mae = mean(abs.(y .- y_pred))
            
            model_result = (
                model_name = model_name,
                model = model,
                params = fit_result.params,
                ssr = fit_result.ssr,
                bic = bic,
                aic = aic,
                aicc = aicc,
                r_squared = r_squared,
                adj_r_squared = adj_r_squared,
                rmse = rmse,
                mae = mae,
                n_params = k,
                complexity_penalty = complexity_penalty,
                y_pred = y_pred,
                fit_success = true
            )
            
            push!(results, model_result)
            
            println("  Success: BIC=$(round(bic, digits=2)), R²=$(round(r_squared, digits=4)), RMSE=$(round(rmse, digits=4))")
            
        catch e
            println("  Failed: $e")
            # Add failed model result
            model_result = (
                model_name = model_name,
                model = model,
                params = fill(NaN, length(p0)),
                ssr = Inf,
                bic = Inf,
                aic = Inf,
                aicc = Inf,
                r_squared = -Inf,
                adj_r_squared = -Inf,
                rmse = Inf,
                mae = Inf,
                n_params = length(p0),
                complexity_penalty = length(p0) / length(y),
                y_pred = fill(NaN, length(y)),
                fit_success = false
            )
            push!(results, model_result)
        end
    end
    
    # Filter successful fits
    successful_results = filter(r -> r.fit_success, results)
    
    if isempty(successful_results)
        println("\nError: No models fitted successfully!")
        return (results = results, ranking = [], recommendations = "No successful fits")
    end
    
    # Rank models by different criteria
    bic_ranking = sort(successful_results, by = r -> r.bic)
    aic_ranking = sort(successful_results, by = r -> r.aic)
    aicc_ranking = sort(successful_results, by = r -> r.aicc)
    r2_ranking = sort(successful_results, by = r -> r.r_squared, rev=true)
    
    # Calculate BIC weights (evidence ratios)
    min_bic = minimum([r.bic for r in successful_results])
    bic_weights = []
    for result in successful_results
        delta_bic = result.bic - min_bic
        weight = exp(-0.5 * delta_bic)
        push!(bic_weights, weight)
    end
    total_weight = sum(bic_weights)
    bic_weights = bic_weights ./ total_weight
    
    # Print comprehensive results
    println("\n=== Model Comparison Results ===")
    println("Model Rankings by BIC (lower is better):")
    for (i, result) in enumerate(bic_ranking)
        weight_idx = findfirst(r -> r.model_name == result.model_name, successful_results)
        weight = bic_weights[weight_idx]
        println("$i. $(result.model_name): BIC=$(round(result.bic, digits=2)), Weight=$(round(weight, digits=3)), R²=$(round(result.r_squared, digits=4))")
    end
    
    # Model selection recommendations
    best_bic = bic_ranking[1]
    println("\n=== Model Selection Recommendations ===")
    println("Best model by BIC: $(best_bic.model_name)")
    
    # Check for substantial evidence
    if length(bic_ranking) > 1
        delta_bic = bic_ranking[2].bic - bic_ranking[1].bic
        if delta_bic > 10
            evidence = "Very strong"
        elseif delta_bic > 6
            evidence = "Strong"
        elseif delta_bic > 2
            evidence = "Positive"
        else
            evidence = "Weak"
        end
        println("Evidence strength: $evidence (ΔBIC = $(round(delta_bic, digits=2)))")
    end
    
    # Additional recommendations
    top_weight = maximum(bic_weights)
    if top_weight > 0.9
        recommendation = "Clear best model: $(best_bic.model_name)"
    elseif top_weight > 0.7
        recommendation = "Probable best model: $(best_bic.model_name), but consider alternatives"
    else
        recommendation = "Model uncertainty present. Consider model averaging or additional data."
    end
    
    println("Recommendation: $recommendation")
    
    return (
        results = results,
        successful_results = successful_results,
        bic_ranking = bic_ranking,
        aic_ranking = aic_ranking,
        r2_ranking = r2_ranking,
        bic_weights = bic_weights,
        best_model = best_bic,
        recommendation = recommendation
    )
end

end # module Analysis
