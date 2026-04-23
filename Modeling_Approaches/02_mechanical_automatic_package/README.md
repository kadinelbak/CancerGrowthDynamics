# 02 Mechanical Automatic Modelling (Package)

This track implements package-driven mechanistic modeling with GrowthParameterEstimation.jl.

## Implemented Structure
- `src/`
- `notebooks/`
- `configs/`
- `outputs/csv/`
- `outputs/images/`
- `outputs/metrics/`
- `outputs/manifests/`

## Condition Notebooks
- `notebooks/monoculture_untreated_auto.ipynb`
- `notebooks/monoculture_treated_auto.ipynb`
- `notebooks/coculture_untreated_auto.ipynb`
- `notebooks/coculture_treated_auto.ipynb`

Each notebook follows the same flow:
1. Activate package environment.
2. Decode condition-specific processed data into a standardized DataFrame.
3. Fit a condition-appropriate model set via GrowthParameterEstimation.jl.
4. Run analysis and sensitivity export.
5. Write summary and run manifest entries.

## Shared Modules
- `src/io_utils.jl`: root discovery, decode helper, output directory management, run manifest writer.
- `src/model_registry.jl`: custom model registration hook for full model-zoo extension.
- `src/fit_workflows.jl`: condition-level fit wrappers and ranking export.
- `src/analysis_workflows.jl`: post-fit metrics, sensitivity export, and analysis plotting.
- `src/MechanicalAutomaticModeling.jl`: entry module including all shared modules.

## Outputs
All condition outputs are written under condition-specific subfolders:
- `outputs/csv/<condition>/`
- `outputs/images/<condition>/`
- `outputs/metrics/<condition>/`

Cross-condition run tracking is appended to:
- `outputs/manifests/run_manifest.csv`

## Environment
Dependencies are defined in `Project.toml` and pinned to GrowthParameterEstimation `0.3`.

## Next Implementation Tasks
- Register complete custom treated/coculture model zoo in `src/model_registry.jl`.
- Extend fit wrappers to staged pipelines where untreated parameters should be inherited.
- Add bootstrap uncertainty wrappers once final sensitivity signature is confirmed against installed package version.
