using GrowthParameterEstimation
using DataFrames
using CSV

# -----------------------------------------------------------------------------
# 1) Load your dataset
# -----------------------------------------------------------------------------
# Required core columns are normalized by normalize_schema().
# In strict_schema mode, include metadata such as:
# - culture_type
# - population_type
# - cell_line
# - treatment_amount
# - density
raw = CSV.read("path/to/your_dataset.csv", DataFrame)
df = normalize_schema(raw)

# -----------------------------------------------------------------------------
# 2) Optional quick checks before fitting
# -----------------------------------------------------------------------------
validate_timeseries(df)
validate_required_metadata(df)

qc = generate_qc_report(df)
qc_paths = save_qc_report(qc; output_dir = "results/qc")
println("QC saved to: ", qc_paths)

# -----------------------------------------------------------------------------
# 3) Configure and run one-shot pipeline
# -----------------------------------------------------------------------------
cfg = default_config(output_dir = "results/one_shot_run")

# Option A: include all registered models in config.model_names
result = run_pipeline(
    df;
    config = cfg,
    strict_schema = true,
    qc_before_fit = true,
)

println("Top ranking rows:")
println(first(result.ranking, min(5, nrow(result.ranking))))
println("Failures logged: ", nrow(result.failures))
println("Exports: ", result.exports)

# -----------------------------------------------------------------------------
# 4) Optional: staged run with provenance + resume support
# -----------------------------------------------------------------------------
stages = default_population_cellline_stages(df; populations = ["naive", "resistant"])

staged = run_staged_pipeline(
    df;
    stages = stages,
    config = default_config(output_dir = "results/staged_run"),
    selection_mode = :best_bic,
    strict_schema = true,
    qc_before_fit = true,
    n_bootstrap = 20,
    export_stage_results = true,
)

println("Staged completed: ", staged.completed)
println("Manifest path: ", staged.manifest_path)

# Resume pattern (if previous run halted):
# resumed = run_staged_pipeline(
#     df;
#     stages = stages,
#     config = default_config(output_dir = "results/staged_run"),
#     selection_mode = :manual,
#     manual_choices = Dict("treated_monoculture_naive_a549" => "theta_logistic_hill_inhibition"),
#     strict_schema = true,
#     resume_manifest_path = staged.manifest_path,
#     resume_from_stage = staged.halted_stage,
# )
