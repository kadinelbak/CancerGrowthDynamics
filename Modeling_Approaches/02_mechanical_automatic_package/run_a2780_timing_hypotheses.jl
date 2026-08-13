using Pkg
Pkg.activate(@__DIR__)

using CSV
using DataFrames

include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

const FW = MechanicalAutomaticModeling.FitWorkflows

start = joinpath(@__DIR__, "notebooks")
max_time = parse(Float64, get(ENV, "A2780_MAX_TIME_PER_FIT", "1.5"))
result = FW._fit_treated_timing_hypotheses(
    start = start,
    max_time_per_fit = max_time,
)
ranking = result.ranking
nrow(ranking) == length(FW.TREATED_TIMING_HYPOTHESES) ||
    error("Not all timing hypotheses completed")
all(isfinite, Float64.(ranking.bic)) ||
    error("Timing ranking contains non-finite BIC")
all(Float64.(ranking.scaled_ssr) .< 9.99e11) ||
    error("Timing ranking contains failure sentinels")

println("timing_rows=", nrow(ranking))
println("winner=", first(ranking.model))
println("winner_bic=", first(ranking.bic))
println("figure_path=", result.figure_path)
