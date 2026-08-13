include(joinpath(@__DIR__, "src", "MechanicalAutomaticModeling.jl"))
using .MechanicalAutomaticModeling

required_models = Set([
    "simeoni_transit_compartment",
    "time_decay_kill_fixed",
    "time_decay_kill_dose_scaled",
])
available_models = Set([spec.name for spec in MechanicalAutomaticModeling.ModelRegistry.local_model_specs()])
missing = setdiff(required_models, available_models)
isempty(missing) || error("Missing treated monoculture literature models: $(collect(missing))")

println("module loaded ok")
println("staged conditions: ", join(MechanicalAutomaticModeling.StagedA2780Workflow.STAGED_A2780_CONDITIONS, ", "))
