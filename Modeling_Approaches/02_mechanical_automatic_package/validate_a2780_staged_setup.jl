using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()
get(ENV, "A2780_PRECOMPILE", "false") == "true" && Pkg.precompile()

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

for (dose, expected_effect, expected_label) in ((0.67, 0.25, "IC25"), (1.0, 0.50, "IC50"), (1.47, 0.75, "IC75"))
    metadata = MechanicalAutomaticModeling.FitWorkflows._treated_dose_metadata(dose)
    metadata.effect_level == expected_effect || error("Incorrect effect level for $(dose) uM")
    metadata.ic_label == expected_label || error("Incorrect IC label for $(dose) uM")
end

println("loaded ok")
println("staged conditions: ", join(MechanicalAutomaticModeling.StagedA2780Workflow.STAGED_A2780_CONDITIONS, ", "))
println("literature models: ", join(sort(collect(required_models)), ", "))
