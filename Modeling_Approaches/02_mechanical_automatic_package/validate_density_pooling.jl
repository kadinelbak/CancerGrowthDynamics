using Pkg
Pkg.activate(@__DIR__)

using Test
using DataFrames
using OrdinaryDiffEq
using GrowthParameterEstimation

include(joinpath(@__DIR__, "dev", "GrowthParameterEstimation", "test", "pooling_profile_test.jl"))
include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

@test isdefined(GrowthParameterEstimation, :run_joint_fit)
@test isdefined(GrowthParameterEstimation, :profile_joint_fit_bounds)
@test isdefined(GrowthParameterEstimation, :summarize_pooling_bic)
@test isdefined(MechanicalAutomaticModeling.FitWorkflows, :_run_density_aware_untreated_fitting)
@test isdefined(MechanicalAutomaticModeling.FitWorkflows, :_run_density_aware_treated_monoculture_fitting)

println("density_pooling_validation=passed")
