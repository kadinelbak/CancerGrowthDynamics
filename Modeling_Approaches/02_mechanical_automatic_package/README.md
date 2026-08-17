# 02 Mechanical Automatic Modelling (Package)

This track implements package-driven mechanistic modeling with
GrowthParameterEstimation.jl.

GrowthParameterEstimation provides the reusable fitting machinery. This
repository keeps the A2780-specific scientific workflow: data decoding, staged
model sequencing, inheritance policy, treatment and coculture hypotheses,
validation scripts, plots, and report rendering.

## Implemented Structure
- `src/`
- `notebooks/`
- `configs/`
- `outputs/csv/`
- `outputs/images/`
- `outputs/metrics/`
- `outputs/manifests/`

## Condition Notebooks
- `notebooks/a2780_staged_model_comparison.ipynb`
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
- `src/coculture_joint.jl`: coupled sensitive/resistant fitting across seeding densities and mix ratios.
- `src/linked_treatment_joint.jl`: linked treated-monoculture and treated-coculture hypothesis tests for A2780 drug-inheritance analysis.
- `src/analysis_workflows.jl`: post-fit metrics, sensitivity export, and analysis plotting.
- `src/adaptive_simulation_engine.jl`: model-agnostic treatment windows, feedback decisions, ensembles, reverse lookup, and outcome metrics.
- `src/a2780_adaptive_adapter.jl`: strict staged-artifact selection and A2780 observable adapter for the simulator.
- `src/MechanicalAutomaticModeling.jl`: entry module including all shared modules.

## Adaptive Simulator

`AdaptiveSimulationEngine` is deliberately independent of A2780 names, files,
HTML, and ODE libraries. A `SimulationScenario` receives a segment-simulation
callback so treatment, refresh, and monitoring events can be orchestrated
without coupling the engine to a particular solver.

The A2780 adapter resolves each default from ranking and status exports. It
rejects candidates with missing/non-finite parameters, diagnostic-only status,
or invalid inheritance and records why a higher-ranked candidate was skipped.
The default browser output contains only A2780Naive, total A2780cis, and total
population. Any internal cis substructure is hidden unless the user enables an
advanced diagnostic view.

Rebuild the browser configuration and Pages mirror with:

```julia
julia --project=. scripts/build_adaptive_simulator_report.jl
```

The JSON and JavaScript configuration files beside the report are generated
from Julia artifacts. Results after day 14 are extrapolation and use an
ensemble sensitivity band. Experimental dose values outside 0.67-1.47 uM are
flagged as extrapolation; translational mode reports normalized exposure only.

## Outputs
All condition outputs are written under condition-specific subfolders:
- `outputs/csv/<condition>/`
- `outputs/images/<condition>/`
- `outputs/metrics/<condition>/`

Cross-condition run tracking is appended to:
- `outputs/manifests/run_manifest.csv`

The staged A2780 report compares shared, symmetric plus-or-minus 5% density
pooling, and independent diagnostic fits. Coculture models use separate measured
sensitive and cis-resistant states, fix both first observations as initial
conditions, and fit all 20k/30k by 25-75/50-50/75-25 environments jointly.
Untreated coculture compares shared-resource and Lotka-Volterra competition
families. Treated coculture inherits the winning interaction law and compares
constant kill, time-decay, delayed-ramp, transit-damage, and
sensitive/tolerant-transition mechanisms.

## Package Boundary
Keep reusable numerical mechanics in GrowthParameterEstimation:
- joint fitting, multistart ranking, and profiling
- fixed initial-time and `u0_builder` support
- trajectory residual scaling and raw/scaled SSE reporting
- finite-fit rejection, BIC summaries, and parameter-stability summaries

Keep A2780 scientific decisions here:
- metadata decoding and treatment-dose correction
- stage ordering and inheritance rules
- A2780 model hypotheses and context modifiers
- plots, reports, validation scripts, and generated outputs

## Linked Treatment Output Policy
The linked treatment analysis fits all nine hypotheses and writes the complete
ranking to `linked_treatment_model_ranking.csv`. Reports and visible summary
tables show only the top five rows via `linked_treatment_top5.csv` and the
legacy `coculture_treated_automatic_best_models_top10.csv` artifact.

## Environment
Dependencies are defined in `Project.toml` and pinned to registered
GrowthParameterEstimation `0.4.1`.

## Validation
Run `validate_density_pooling.jl` for monoculture pooling tests and
`validate_a2780_coculture_joint.jl` for coculture provenance, finite-ranking,
pooling, model-coverage, and graph-table checks.
