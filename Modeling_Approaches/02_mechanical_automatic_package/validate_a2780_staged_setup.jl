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

baseline = (model = "logistic_growth", r = 0.5, K = 4000.0, shape_value = NaN)
joint_specs = MechanicalAutomaticModeling.FitWorkflows._joint_treated_model_specs(
    [0.67, 1.0, 1.47, 0.67, 1.0, 1.47],
    [0.25, 0.50, 0.75, 0.25, 0.50, 0.75],
    [1, 1, 1, 2, 2, 2],
    NamedTuple[baseline for _ in 1:6],
    0.0,
    "shared",
)
@assert haskey(joint_specs, "joint_intracellular_platinum_pkpd")
pkpd = joint_specs["joint_intracellular_platinum_pkpd"]
@assert pkpd.layout == :platinum_pkpd
@assert pkpd.param_names == [:emax, :ec50_bound, :hill_n, :k_efflux, :k_repair]

datasets = [(x = [1.0, 2.0, 3.0], y = [100.0, 120.0, 130.0], residual_scale = 130.0) for _ in 1:6]
mapped, pkpd_u0, pkpd_u0_builder = MechanicalAutomaticModeling.FitWorkflows._joint_model_inputs(
    pkpd,
    datasets,
    fill(100.0, 6),
    [1, 1, 1, 2, 2, 2],
)
@assert length(mapped) == 6
@assert length(pkpd_u0) == 18
@assert pkpd_u0_builder === nothing
du = similar(pkpd_u0)
pkpd.model(du, pkpd_u0, pkpd.p0, 0.0)
@assert all(isfinite, du)
@assert du[2] > 0 && du[3] == 0

for (dose, expected_effect, expected_label) in ((0.67, 0.25, "IC25"), (1.0, 0.50, "IC50"), (1.47, 0.75, "IC75"))
    metadata = MechanicalAutomaticModeling.FitWorkflows._treated_dose_metadata(dose)
    metadata.effect_level == expected_effect || error("Incorrect effect level for $(dose) uM")
    metadata.ic_label == expected_label || error("Incorrect IC label for $(dose) uM")
end

println("loaded ok")
println("staged conditions: ", join(MechanicalAutomaticModeling.StagedA2780Workflow.STAGED_A2780_CONDITIONS, ", "))
println("literature models: ", join(sort(collect(required_models)), ", "))
