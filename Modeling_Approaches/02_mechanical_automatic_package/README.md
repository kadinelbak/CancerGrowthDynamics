# 02 Mechanical Automatic Modelling (Package)

Autonomous compartment models of A2780 (cisplatin-sensitive) and A2780cis (resistant) ovarian
cancer cells, fitted stage by stage to the New Datasets experiments, then used for optimal control
of the cisplatin schedule (and, over 30 days, of fasting). No model uses time directly: delays and
accumulation come from states driven by the inputs, so every model can be driven by any dose
schedule D(t). Fitting uses GrowthParameterEstimation.jl, pinned in `Manifest.toml`.

## Run everything

```
julia run_pipeline.jl                  # raw data → fits → analyses → control → tests → Final Report
julia run_pipeline.jl --list           # the steps
julia run_pipeline.jl --from analyses  # start later, reusing earlier outputs
julia run_pipeline.jl --resume         # reuse unchanged fits, bootstraps and control checkpoints
```

Each step runs in its own process with its log in `outputs/logs/`. `outputs/run_manifest.toml`
records the git commit, Julia and package versions, thread count, SHA-256 of every data file and the
time and exit code of each step. A full run takes about 5–6 hours on 20 threads.

**Reproducibility.** Data selection, stages and hold-outs are configuration; bootstrap replicates and
multistarts are seeded; optimisers stop on iteration caps. Each optimiser call also has a wall-clock
safety limit (300 s for the stage fits, 120 s per bootstrap fit), because without one the bounded BFGS
can run for a very long time on a few candidates. A fit that reaches 95% of its limit may give a
slightly different answer on another machine, so every such fit is listed in
`outputs/zoo_fits/time_limit_hits.csv` and counted in the run manifest; a run with 0 hits is
reproducible, which `test/test_regression.jl` checks against `test/reference_results.toml`.

## Configuration (no data, stage or hold-out is hard-coded)

| File | What it defines |
|---|---|
| `data/sources.toml` | every raw data file (and sheet), tagged with its experimental conditions: condition, line, mono/co culture, seeding ratio, run, seeding, serum (FBS %), dose schedule, time origin, QC rules |
| `data/schedules.toml` | dose schedules D(t) as piecewise-constant steps |
| `data/qc_rules.toml` | exclusion rules (applied only to the sources that list them) |
| `data/measurements.csv` | the tidy table built from the above by `data/build_tidy.jl`: one row per measurement, excluded points kept with their reason |
| `stages.toml` | the lines (symbol, measured IC50), the selection rule, and every fitting stage: its data (a tag selector), what it inherits, its candidates and free parameters, constraints |
| `validation.toml` | forecasts, leave-one-out blocks, the held-out treated cocultures and the intermittent-schedule prediction |
| `control_config.toml` | optimal control: decision interval, doses, toxicity proxy, weights, fasting and host parameters, each labelled DATA, LITERATURE or ASSUMPTION |

### Adding data
1. Add a `[[source]]` entry to `data/sources.toml` with the file, the sheet for `.xlsx`, and its tags;
   add a schedule to `data/schedules.toml` if it is new.
2. If it belongs to an existing stage, the stage's tag selector picks it up. A new ratio, run or
   condition becomes a new environment of that stage. Otherwise add a stage to `stages.toml`.
3. `julia run_pipeline.jl`.

## Stages (pages show the thesis names)

| id | page name | data | model |
|---|---|---|---|
| s1 | 1 | untreated monocultures | growth law per line |
| s1b | 2a | low-resource monocultures | Stage 1 growth + every low-resource mechanism (A/B shared-structure test) |
| s2 | 2b | treated monocultures, 3 doses jointly | Stage 1 growth + every cisplatin mechanism; one structure shared by both lines |
| s3 | 3a | untreated (Ce0) + 1 µM (Ce1) cocultures | both lines, shared K; the resistant kill is never switched off (floor 50% of Stage 2b) |
| s3b | 3b | low-resource cocultures | both lines with their Stage 2a mechanisms |
| s4 | 4 | low resource + 1 µM coculture | both mechanisms per line, inherited; a few refitted parameters |

At every stage the carried model is the candidate with the fewest free parameters within 3 BIC of
the lowest (lower BIC breaks ties). Every candidate is refitted to 200 residual block-bootstrap
replicates. The Hill kill uses each line's measured IC50 as EC50 (A2780 1 µM, A2780cis 8.34 µM).

## Main scripts

| Script | Role |
|---|---|
| `zoo_data.jl` | selects data by tag from the tidy table |
| `zoo_models.jl` | the model engine: growth laws, delay kernels, actions, inputs D(t), equations as TeX |
| `run_zoo_fits.jl` | runs the stages of `stages.toml`, bootstraps, writes `outputs/zoo_fits/` and the site bundle |
| `zoo_analyses.jl` | package bootstrap, sensitivity (local, profile, Sobol), validation (`validation.toml`) |
| `run_control.jl`, `control_report.jl` | 14-day dose control (sequential MPC with a toxicity proxy) |
| `run_control_fasting.jl`, `control_fasting_report.jl` | 30-day dose + fasting control |
| `../../build_final_report.jl` | the Final Report from every chapter |

## Tests

| Test | Checks |
|---|---|
| `test/test_data.jl` | the tidy table reproduces the earlier hard-coded loaders exactly |
| `test/test_stages.jl` | `stages.toml` rebuilds every candidate of the saved run (`DRYRUN=1`) |
| `test/test_validation_config.jl` | `validation.toml` selects the intended data |
| `test/check_intermittent.jl` | the engine reproduces the earlier intermittent-schedule predictions |
| `test/test_control.jl`, `test/test_control_fasting.jl` | plants, linearisation, gradients and controllers |
| `test/test_regression.jl` | a rerun reproduces the reference results (`test/reference_results.toml`) |

---

# Legacy (time-based and earlier package workflows)

The material below describes the earlier notebook and staged-report workflows. They are kept for
reference; the time-based models are archived in `archive_time_models/`.

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
- `src/sample_aware_report.jl`: well-preserving staged summaries, between-well variability bands, direction audits, and report rendering.
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
julia +1.10.4 --compiled-modules=no --startup-file=no --project=. scripts/build_adaptive_simulator_report.jl
```

The JSON and JavaScript configuration files beside the report are generated
from Julia artifacts. Results after day 14 are extrapolation and use an
ensemble sensitivity band. Experimental dose values outside 0.67-1.47 uM are
flagged as extrapolation; translational mode reports normalized exposure only.

## Sample-Aware Staged Report

The sample-aware companion report keeps image tiles averaged within each well,
but does not average the replicate wells into a single observed trajectory. It
shows each well, the across-well mean plus or minus one standard deviation, and
the selected canonical staged fit. Its rankings and equations are inherited
from the validated staged report; the well-level view audits reproducibility
and does not recompute ordinary BIC with correlated measurements treated as
independent.

Rebuild the report, CSV findings, figures, and Pages mirror with:

```julia
julia +1.10.4 --compiled-modules=no --startup-file=no --project=. scripts/build_sample_aware_report.jl
```

The generated report is
`outputs/reports/a2780_sample_aware_staged_model_comparison.html`.

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

The checked-in manifest was resolved with Julia `1.10.4`. Use that runtime for
reproducible local validation.

### Windows Application Control

On Windows systems where Application Control blocks Julia-generated package
image DLLs under `.julia/compiled`, run this project with compiled modules
disabled. This bypasses Julia's cache files without changing Windows security
policy or the scientific environment:

```powershell
./scripts/run_tests_windows.ps1
```

For report or validation scripts, use the same launch options:

```powershell
julia +1.10.4 --compiled-modules=no --startup-file=no --project=. path/to/script.jl
```

The first cache-free load is slower because code is compiled in memory. Do not
use an unqualified Julia `1.12` invocation with this manifest: it re-resolves
dependencies and can create native package images that the local Windows policy
then rejects as `Bad Image` files.

## Validation
Run `validate_density_pooling.jl` for monoculture pooling tests and
`validate_a2780_coculture_joint.jl` for coculture provenance, finite-ranking,
pooling, model-coverage, and graph-table checks.
