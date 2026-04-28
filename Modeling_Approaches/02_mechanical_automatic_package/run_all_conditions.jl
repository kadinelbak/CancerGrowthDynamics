using MechanicalAutomaticModeling

conditions = ["monoculture_untreated", "monoculture_treated", "coculture_untreated", "coculture_treated"]

for cond in conditions
    println("\n" * "="^60)
    println("CONDITION: $cond")
    println("="^60)

    state = Dict{String,Any}()
    try
        decode_condition_dataframe!(state, cond)
        println("  decode_ok  : true  (rows=$(nrow(state["data"])))")
    catch e
        println("  decode_ok  : false  ($(e))")
        continue
    end

    try
        run_condition_fit!(state, cond)
        n = nrow(state["fit_results"])
        best = state["fit_results"][argmin(state["fit_results"].bic), :]
        println("  fit_ok     : true  (fit_rows=$n)")
        println("  best_model : $(best.model_name)  BIC=$(round(best.bic, digits=2))")
    catch e
        println("  fit_ok     : false  ($(e))")
        continue
    end

    try
        run_condition_analysis!(state, cond)
        n = nrow(state["sensitivity_results"])
        println("  analysis_ok: true  (sensitivity_rows=$n)")
    catch e
        println("  analysis_ok: false  ($(e))")
    end

    println("  DONE")
end

println("\n" * "="^60)
println("All conditions complete.")
