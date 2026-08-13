using CSV
using DataFrames
using MechanicalAutomaticModeling

const SW = MechanicalAutomaticModeling.StagedA2780Workflow
const FW = MechanicalAutomaticModeling.FitWorkflows
const IO = MechanicalAutomaticModeling.IOUtils

@assert IO._infer_metadata_from_path(joinpath("Processed_Datasets", "Treated MonoCulture", "20k", "IC25", "A2780Naive.csv")).dose == 1.47
@assert IO._infer_metadata_from_path(joinpath("Processed_Datasets", "Treated MonoCulture", "20k", "IC50", "A2780Naive.csv")).dose == 1.0
@assert IO._infer_metadata_from_path(joinpath("Processed_Datasets", "Treated MonoCulture", "20k", "IC75", "A2780Naive.csv")).dose == 0.67

expected_total(density) = density == "20k" ? 67.0 : density == "30k" ? 100.0 : error("Unexpected density $(density)")
same_number(left, right) = (isfinite(left) && isfinite(right)) ? isapprox(left, right; rtol = 1e-10, atol = 1e-10) : (!isfinite(left) && !isfinite(right))
safe_string(value) = ismissing(value) ? "" : String(value)

const LINKED_CONTEXT_MODIFIED_MODELS = Set([
    "tolerant_context_shift",
    "load_plus_tolerant_context",
    "load_plus_tolerant_growth_context",
])

const LINKED_CONTEXT_MODIFIED_PARAMETERS = Set([
    "cis_emax_tolerant",
    "cis_f_tolerant0",
])

const LINKED_PROVENANCE_ARTIFACTS = Set([
    "monoculture_treated_a2780_decoded.csv",
    "coculture_treated_a2780_decoded.csv",
    "monoculture_treated_joint_parameter_estimates.csv",
    "monoculture_treated_pooling_status.csv",
    "monoculture_treated_timing_hypothesis_ranking.csv",
    "monoculture_treated_timing_hypothesis_parameters.csv",
])

function assert_linked_ranking_policy(linked_ranking, linked_top5)
    @assert nrow(linked_ranking) == length(FW.LINKED_TREATMENT_HYPOTHESES)
    @assert Set(String.(linked_ranking.model)) == Set(FW.LINKED_TREATMENT_HYPOTHESES)
    @assert nrow(linked_top5) == FW.LINKED_TREATMENT_VISIBLE_LIMIT
    @assert String.(linked_top5.model) == String.(first(linked_ranking, FW.LINKED_TREATMENT_VISIBLE_LIMIT).model)
end

function assert_linked_effective_inheritance(linked_effective, winning_model)
    context_modified_model = String(winning_model) in LINKED_CONTEXT_MODIFIED_MODELS
    for row in eachrow(linked_effective)
        ratio = Float64(row.ratio)
        if !isfinite(ratio)
            @assert !isfinite(Float64(row.treated_monoculture_value))
            @assert !isfinite(Float64(row.treated_coculture_intrinsic_value))
        elseif context_modified_model && String(row.parameter) in LINKED_CONTEXT_MODIFIED_PARAMETERS
            @assert ratio > 0.0
        else
            @assert isapprox(ratio, 1.0; rtol = 1e-12, atol = 1e-12)
        end
    end
end

untreated_mono_out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs("monoculture_untreated"; start = pwd())
treated_mono_out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs("monoculture_treated"; start = pwd())
mono_baselines = CSV.read(joinpath(untreated_mono_out.csv, "untreated_group_baselines.csv"), DataFrame)

for (out, filename) in (
    (untreated_mono_out, "monoculture_untreated_initial_condition_diagnostics.csv"),
    (treated_mono_out, "monoculture_treated_joint_initial_condition_diagnostics.csv"),
)
    initial = CSV.read(joinpath(out.csv, filename), DataFrame)
    @assert all(row -> same_number(Float64(row.fixed_u0), expected_total(String(row.density))), eachrow(initial))
    @assert all(Float64.(initial.u0_time_day) .== 0.0)
end

treated_mono_audit = CSV.read(joinpath(treated_mono_out.csv, "monoculture_treated_inheritance_audit.csv"), DataFrame)
@assert nrow(treated_mono_audit) == nrow(mono_baselines)
for audit in eachrow(treated_mono_audit)
    baseline = mono_baselines[
        (String.(mono_baselines.cell_line) .== String(audit.treated_cell_line)) .&
        (String.(mono_baselines.density) .== String(audit.density)),
        :,
    ]
    @assert nrow(baseline) == 1
    baseline = first(baseline)
    @assert String(audit.untreated_growth_family) == String(baseline.best_model)
    @assert same_number(Float64(audit.untreated_r), Float64(baseline.r))
    @assert same_number(Float64(audit.untreated_K), Float64(baseline.K))
    @assert safe_string(audit.untreated_shape_parameter) == safe_string(baseline.shape_parameter)
    @assert same_number(Float64(audit.untreated_shape_value), Float64(baseline.shape_value))
end

for condition in ("coculture_untreated", "coculture_treated")
    decoded = SW.decode_a2780_condition(condition; start = pwd())
    recovered, provenance = FW._recover_coculture_design(decoded, condition)
    environments, initial = FW._coculture_environments(recovered)
    @assert length(environments) == 6
    @assert Set(environment.density for environment in environments) == Set(["20k", "30k"])
    @assert Set(environment.mix for environment in environments) == Set(["25-75", "50-50", "75-25"])
    @assert all(environment -> same_number(environment.fixed_total_u0, expected_total(environment.density)), environments)
    @assert all(environment -> same_number(environment.u0_sensitive + environment.u0_resistant, expected_total(environment.density)), environments)
    @assert all(environment -> same_number(environment.u0_sensitive / environment.fixed_total_u0, FW._mix_fraction(environment.mix)), environments)
    @assert all(isfinite, initial.observed_sensitive_fraction)
    if condition == "coculture_treated"
        finite_errors = filter(isfinite, provenance.absolute_check_error)
        @assert !isempty(finite_errors)
        @assert maximum(finite_errors) < 0.02
    end

    out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs(condition; start = pwd())
    ranking = CSV.read(joinpath(out.csv, "$(condition)_pooling_model_ranking.csv"), DataFrame)
    @assert nrow(ranking) >= (condition == "coculture_untreated" ? 9 : 21)
    @assert all(isfinite, ranking.bic)
    @assert all(abs.(ranking.bic) .< 9.99e11)
    @assert all(isfinite, ranking.ssr)
    @assert all(ranking.ssr .< 9.99e11)
    status = first(CSV.read(joinpath(out.csv, "$(condition)_pooling_status.csv"), DataFrame))
    @assert !Bool(status.inadequate_pooling)
    @assert String(status.winning_pooling_mode) in ("shared", "partial_5pct")

    overlay = CSV.read(joinpath(out.csv, "figures", "$(condition)_joint_overlays.csv"), DataFrame)
    winner = overlay[
        (String.(overlay.model) .== String(status.winning_model)) .&
        (String.(overlay.pooling_mode) .== String(status.winning_pooling_mode)),
        :,
    ]
    coverage = combine(groupby(winner, [:density, :mix, :component]), :time => (values -> length(unique(values))) => :n_times)
    @assert nrow(coverage) == 12
    @assert all(coverage.n_times .>= 14)
    for trajectory in groupby(winner, [:density, :mix, :component])
        anchor = trajectory[Bool.(trajectory.fixed_day0_anchor), :]
        @assert nrow(anchor) == 1
        @assert Float64(first(anchor.time)) == 0.0
        fraction = FW._mix_fraction(String(first(anchor.mix)))
        expected = expected_total(String(first(anchor.density))) * (String(first(anchor.component)) == "sensitive" ? fraction : 1 - fraction)
        @assert same_number(Float64(first(anchor.observed)), expected)
        @assert same_number(Float64(first(anchor.predicted)), expected)
    end
    @assert isfile(joinpath(out.images, "figures", "$(condition)_best_mechanistic_fit_grid.png"))
end

coculture_untreated_out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs("coculture_untreated"; start = pwd())
coculture_treated_out = MechanicalAutomaticModeling.IOUtils.condition_output_dirs("coculture_treated"; start = pwd())
coculture_mono_audit = CSV.read(joinpath(coculture_untreated_out.csv, "coculture_untreated_monoculture_inheritance_audit.csv"), DataFrame)
@assert nrow(coculture_mono_audit) == 4
@assert all(String.(coculture_mono_audit.carrying_capacity_mode) .== "exact_fixed_untreated_K_no_rescaling")
for audit in eachrow(coculture_mono_audit)
    baseline = mono_baselines[
        (String.(mono_baselines.cell_line) .== String(audit.monoculture_cell_line)) .&
        (String.(mono_baselines.density) .== String(audit.density)),
        :,
    ]
    @assert nrow(baseline) == 1
    baseline = first(baseline)
    @assert String(audit.growth_family) == String(baseline.best_model)
    @assert same_number(Float64(audit.r), Float64(baseline.r))
    @assert same_number(Float64(audit.K), Float64(baseline.K))
    @assert safe_string(audit.shape_parameter) == safe_string(baseline.shape_parameter)
    @assert same_number(Float64(audit.shape_value), Float64(baseline.shape_value))
end

coculture_baseline = first(CSV.read(joinpath(coculture_untreated_out.csv, "coculture_untreated_joint_baseline.csv"), DataFrame))
@assert !occursin("capacity_scale", String(coculture_baseline.parameter_names))
@assert occursin("no carrying-capacity rescaling", String(coculture_baseline.monoculture_growth_inheritance))
coculture_treated_audit = CSV.read(joinpath(coculture_treated_out.csv, "coculture_treated_inheritance_audit.csv"), DataFrame)
baseline_names = split(String(coculture_baseline.parameter_names), ";")
baseline_values = parse.(Float64, split(String(coculture_baseline.parameter_values), ";"))
@assert String.(coculture_treated_audit.parameter) == baseline_names
@assert all(same_number.(Float64.(coculture_treated_audit.inherited_value), baseline_values))
@assert all(String.(coculture_treated_audit.inherited_model) .== String(coculture_baseline.model))
@assert all(String.(coculture_treated_audit.inherited_pooling_mode) .== String(coculture_baseline.pooling_mode))

treated_ranking = CSV.read(joinpath("outputs", "csv", "coculture_treated", "coculture_treated_pooling_model_ranking.csv"), DataFrame)
@assert "dual_transit_damage" in String.(treated_ranking.model)
@assert "dual_transit_competitor_scaled" in String.(treated_ranking.model)
@assert "dual_transit_load_scaled" in String.(treated_ranking.model)
@assert "dual_time_decay_kill" in String.(treated_ranking.model)
@assert "sensitive_tolerant_transition" in String.(treated_ranking.model)

scaled_models = treated_ranking[in.(String.(treated_ranking.model), Ref(["dual_transit_competitor_scaled", "dual_transit_load_scaled"])), :]
@assert nrow(scaled_models) == 6
@assert Set(String.(scaled_models.pooling_mode)) == Set(["shared", "partial_5pct", "independent_diagnostic"])
@assert all(Int.(scaled_models.n_densities) .== 2)
@assert all(Int.(scaled_models.n_mixes) .== 3)
@assert length(unique(Int.(scaled_models.n_points))) == 1
@assert first(Int.(scaled_models.n_points)) >= 144
@assert all(occursin.("multiplicative", String.(scaled_models.treatment_interaction_mode)))

nonadditive = CSV.read(joinpath(coculture_treated_out.csv, "coculture_treated_nonadditive_model_comparison.csv"), DataFrame)
@assert nrow(nonadditive) == 3
@assert Set(String.(nonadditive.model)) == Set(["dual_transit_damage", "dual_transit_competitor_scaled", "dual_transit_load_scaled"])
@assert all(isfinite, Float64.(nonadditive.bic))
@assert all(Float64.(nonadditive.delta_bic) .>= 0.0)
@assert isfile(joinpath(coculture_treated_out.images, "figures", "coculture_treated_nonadditive_simulation_grid.png"))

treated_parameters = CSV.read(joinpath(coculture_treated_out.csv, "coculture_treated_joint_parameter_estimates.csv"), DataFrame)
for model_name in ("dual_transit_competitor_scaled", "dual_transit_load_scaled")
    model_parameters = treated_parameters[String.(treated_parameters.model) .== model_name, :]
    @assert Set(["beta_sensitive", "beta_resistant"]) <= Set(String.(model_parameters.parameter))
end

linked_ranking = CSV.read(joinpath(coculture_treated_out.csv, "linked_treatment_model_ranking.csv"), DataFrame)
linked_top5 = CSV.read(joinpath(coculture_treated_out.csv, "linked_treatment_top5.csv"), DataFrame)
assert_linked_ranking_policy(linked_ranking, linked_top5)
@assert all(isfinite, Float64.(linked_ranking.bic))
@assert all(isfinite, Float64.(linked_ranking.ssr))
@assert all(Int.(linked_ranking.n_points) .== 336)
@assert all(Int.(linked_ranking.n_monoculture_trajectories) .== 12)
@assert all(Int.(linked_ranking.n_coculture_component_trajectories) .== 12)
@assert all(Bool.(linked_ranking.eligible_for_inheritance[String.(linked_ranking.model) .!= "fully_free_context_diagnostic"]))
@assert !Bool(only(linked_ranking.eligible_for_inheritance[String.(linked_ranking.model) .== "fully_free_context_diagnostic"]))

linked_status = first(CSV.read(joinpath(coculture_treated_out.csv, "linked_treatment_status.csv"), DataFrame))
@assert String(linked_status.winning_model) != "fully_free_context_diagnostic"
@assert String(linked_status.winning_pooling_mode) == "linked_global"
@assert isfinite(Float64(linked_status.winning_bic))

linked_winner = first(linked_ranking[String.(linked_ranking.model) .== String(linked_status.winning_model), :])
serialized_values = match(r"values = \[([^\]]+)\]", String(linked_winner.params))
@assert serialized_values !== nothing
ranking_parameters = parse.(Float64, strip.(split(serialized_values.captures[1], ",")))
linked_parameter_audit = CSV.read(joinpath(coculture_treated_out.csv, "linked_treatment_shared_parameter_audit.csv"), DataFrame)
@assert length(ranking_parameters) == nrow(linked_parameter_audit)
@assert all(isapprox.(ranking_parameters, Float64.(linked_parameter_audit.estimate); rtol = 1e-12, atol = 1e-12))
f_tolerant = linked_parameter_audit[String.(linked_parameter_audit.parameter) .== "cis_f_tolerant0", :]
@assert nrow(f_tolerant) == 1
@assert 0.0 <= Float64(only(f_tolerant.estimate)) <= 0.95

linked_effective = CSV.read(joinpath(coculture_treated_out.csv, "linked_treatment_effective_parameter_inheritance.csv"), DataFrame)
@assert nrow(linked_effective) == 20
@assert Set(String.(linked_effective.lineage)) == Set(["A2780Naive", "A2780cis"])
@assert Set(String.(linked_effective.density)) == Set(["20k", "30k"])
assert_linked_effective_inheritance(linked_effective, linked_status.winning_model)
@assert all(String.(linked_effective.inheritance_mode) .== "same parameter index in one combined global objective")

linked_overlay = CSV.read(joinpath(coculture_treated_out.csv, "figures", "linked_treatment_combined_overlays.csv"), DataFrame)
winner_overlay = linked_overlay[String.(linked_overlay.model) .== String(linked_status.winning_model), :]
mono_coverage = combine(groupby(winner_overlay[String.(winner_overlay.context) .== "monoculture", :], [:cell_line, :density, :dose]), nrow => :n_rows)
co_coverage = combine(groupby(winner_overlay[String.(winner_overlay.context) .== "coculture", :], [:density, :mix, :component]), nrow => :n_rows)
@assert nrow(mono_coverage) == 12
@assert nrow(co_coverage) == 12
for trajectory in groupby(winner_overlay, [:context, :cell_line, :density, :dose, :mix, :component])
    anchor = trajectory[Bool.(trajectory.fixed_day0_anchor), :]
    @assert nrow(anchor) == 1
    @assert Float64(first(anchor.time)) == 0.0
    @assert same_number(Float64(first(anchor.observed)), Float64(first(anchor.predicted)))
    if String(first(anchor.context)) == "monoculture"
        @assert same_number(Float64(first(anchor.observed)), expected_total(String(first(anchor.density))))
    else
        fraction = FW._mix_fraction(String(first(anchor.mix)))
        expected = expected_total(String(first(anchor.density))) * (String(first(anchor.component)) == "sensitive" ? fraction : 1 - fraction)
        @assert same_number(Float64(first(anchor.observed)), expected)
    end
end

linked_provenance = CSV.read(joinpath(coculture_treated_out.csv, "linked_treatment_data_provenance.csv"), DataFrame)
@assert Set(String.(linked_provenance.artifact)) == LINKED_PROVENANCE_ARTIFACTS
@assert all(Bool.(linked_provenance.exists))
@assert all(length.(String.(linked_provenance.sha256)) .== 64)
@assert all(Bool.(linked_provenance.canonical_processed_data_present))
@assert all(String.(linked_provenance.canonical_source_commit) .== FW.LINKED_TREATMENT_CANONICAL_DATA_COMMIT)
@assert all(occursin.("restored from Git history", String.(linked_provenance.source_mode)))
@assert all(String.(linked_provenance.treated_monoculture_label_correction) .== "legacy IC25 folder -> 1.47 uM/IC75; legacy IC75 folder -> 0.67 uM/IC25")
for image_name in (
    "linked_treatment_monoculture_grid.png",
    "linked_treatment_coculture_grid.png",
    "linked_treatment_hypothesis_comparison_grid.png",
)
    @assert isfile(joinpath(coculture_treated_out.images, "figures", image_name))
end

println("A2780 coculture joint validation passed")
