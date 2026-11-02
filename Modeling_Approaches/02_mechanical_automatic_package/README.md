# 02 Mechanical Automatic Modelling (Package)

Autonomous compartment models of A2780 (cisplatin-sensitive) and A2780cis (resistant) ovarian
cancer cells, fitted stage by stage to the New Datasets experiments, then used for optimal control
of the cisplatin schedule (and, over 120 days, of fasting). No model uses time directly: delays and
accumulation come from states driven by the inputs, so every model can be driven by any dose
schedule D(t). Fitting uses GrowthParameterEstimation.jl, pinned in `Manifest.toml`.

## Setup (clean clone)

Julia 1.12.7 (juliaup). The raw data live in `New Datasets/` at the repository root (see `data/sources.toml`).

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # installs the pinned packages, incl. GrowthParameterEstimation
```

Checkpoints (`*.jls`) are not in git, so a clean clone starts with `julia run_pipeline.jl` (no `--from`).
Long runs: keep the machine awake and the lid open (a closed lid sleeps it); `sensitivity/run_awake.ps1` runs a
script detached with a Windows keep-awake request.

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

**Reused data are flagged.** A control step whose inputs are unchanged skips its work and writes
`outputs/reuse/<step>.toml`; `run_pipeline.jl` copies it into the step's `reused_from` (time of the run whose outputs
are reused) in the manifest. When the Final Report is built, `report_history.jl` reads the manifest and records a
`[provenance]` section in `version.toml` (each of fits, analyses, control, long control, 30/60-day control and sweeps:
new in this version, reused with unchanged inputs, reused from an earlier run via `--from`/`--only`, or resumed with
`--resume`). It is shown at the top of the Final Report and on the history card, which carries a "partly reused" tag if
anything was reused. Versions built without a manifest show "provenance not recorded".

**Headline scenario.** `control_config.toml` `[scenarios] headline` (user, 2026-11-02: `"B"`) names the main
drug-response scenario. Both scenarios are still computed; the control, sweep and fasting reports list the headline first
as "main result" and the other as "sensitivity case", and the Final Report history adds a `q+:robust_main` row.

**Reproducibility.** Data selection, stages and hold-outs are configuration; bootstrap replicates and
multistarts are seeded; optimisers stop on iteration caps only (BFGS 1500, Nelder–Mead 2500) — there is no
wall-clock limit (since 2026-10-07), so results do not depend on machine load, threads or sleep.
`test/test_regression.jl` fails unless `test/reference_results.toml` was made (by `test/make_reference.jl`)
with the same configuration and data (SHA-256 fingerprint); it refits Stages 1 and 2b and the first
bootstrap replicates from scratch and compares the control tables. After an intended change, accept the
new run with `make_reference.jl`, then rerun the test.

**Starting values (user, 2026-10-07).** Day 0 is the nominal seeding — 100 cells per well at 30k, 66.7 at 20k,
split by the seeding ratio in cocultures — and is the initial condition, not a data point. A series whose first
measured count is more than 3-fold above or below it gets its starting value estimated (one extra parameter
per series, counted in the BIC; `zoo_data.jl` `x0_info`, `bootstrap_tools.jl` `fit_cand`).

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

### Sensitivity runs (without touching the main configuration)

Copy a configuration file, change it, and point the scripts at the copy:

| Variable | Replaces |
|---|---|
| `SOURCES_TOML`, `QC_TOML` | `data/sources.toml`, `data/qc_rules.toml` (read by `data/build_tidy.jl`; `SOURCES_TOML` also by the fits) |
| `DATA_TABLE` | `data/measurements.csv` (written by `build_tidy.jl`, read by the fits) |
| `STAGES_TOML` | `stages.toml` |
| `ZOO_OUT` | `outputs/zoo_fits/`; a run with `ZOO_OUT` publishes nothing (no docs copy, site bundle or Final Report) |

Options in `stages.toml [selection]` (defaults reproduce the reference run):

- `residual_weighting = "none" | "sd" | "power"` (with `weight_power`, default 0.55): per-point residual
  scales, the replicate SD floored at the series' median SD, or max(y, 50)^p. Weighted BICs compare
  models within one scheme only.
- `bootstrap_generator = "most_complex" | "selected"`: which candidate generates the bootstrap replicates.

`[starting_values] estimate = "mono"` estimates the day-0 count of every monoculture series as an extra
parameter (counted in the BIC) instead of fixing it at the nominal seeding; `"auto"` estimates it only when
the first count is far from the seeding, `"all"` also for cocultures (user 2026-11-02; `outputs/zoo_fits/starting_values.csv`).
`[selection] extra_starts = 20` adds seeded Latin-hypercube starts to the fixed starts of every main fit
(`EXTRA_STARTS` overrides, 0 = off), drawn within `extra_starts_spread` (10×) of the fixed starts and each limited
to `extra_starts_work_factor` (20×) the work of the slowest fixed-start fit, so a start in a stiff corner is dropped
rather than run for hours; `outputs/zoo_fits/fit_starts.csv` lists per candidate how many starts reached the best BIC.

A coculture candidate with `kill_growth_coupled = true` multiplies each line's drug kill by its growth
factor (Norton–Simon: 1 − x for logistic, 1 − x^θ for theta-logistic growth, x = crowding/K).

Example: β_C floor at 25% without the outlier exclusions,
```
SOURCES_TOML=sens/sources.toml DATA_TABLE=sens/measurements.csv julia --project=. data/build_tidy.jl
STAGES_TOML=sens/stages.toml SOURCES_TOML=sens/sources.toml DATA_TABLE=sens/measurements.csv ZOO_OUT=sens/fits julia -t 20 --project=. run_zoo_fits.jl
```

## Stages (pages show the thesis names)

| id | page name | data | model |
|---|---|---|---|
| s1 | 1 | untreated monocultures | growth law per line |
| s1b | 2a | low-resource monocultures | Stage 1 growth + every low-resource mechanism (A/B shared-structure test) |
| s2 | 2b | treated monocultures, 4 doses jointly (A2780cis 1 µM day 14 left out, `sources.toml exclude_days`) | Stage 1 growth + every cisplatin mechanism; one structure shared by both lines |
| s3 | 3a | untreated (Ce0) + 1 µM (Ce1) cocultures | both lines, shared K; the resistant kill is refitted but kept at or above `betaC_floor` × its Stage 2b value (user 2026-11-02; no "_noC" candidates) |
| s3b | 3b | low-resource cocultures | both lines with their Stage 2a mechanisms; K and the crowding coefficients are taken from Stage 3a, only the resource sensitivities ε may be refitted (user 2026-11-02) |
| s4 | 4 | low resource + 1 µM coculture | both mechanisms per line, inherited; a few refitted parameters |

At every stage the carried model is the candidate with the fewest free parameters within 2 BIC of
the lowest (lower BIC breaks ties). Every candidate is refitted to 200 residual block-bootstrap
replicates. The Hill kill uses each line's measured IC50 as EC50 (A2780 1 µM, A2780cis 8.34 µM).

## Main scripts

| Script | Role |
|---|---|
| `zoo_data.jl` | selects data by tag from the tidy table |
| `zoo_models.jl` | the model engine: growth laws, delay kernels, actions, inputs D(t), equations as TeX |
| `run_zoo_fits.jl` | runs the stages of `stages.toml`, bootstraps, writes `outputs/zoo_fits/` and the site bundle |
| `zoo_analyses.jl` | package bootstrap, sensitivity (local, profile, Sobol), validation (`validation.toml`) |
| `run_control.jl`, `control_report.jl` | 14-day dose control (sequential MPC with a toxicity proxy); `CONTROL_DAYS=30` / `60` for the longer runs |
| `run_control_sweep.jl`, `control_sweep_report.jl` | control sweep over burden-barrier shape, threshold, run length, look-ahead mode and objective (`SWEEP_DATA` + `ZOO_STATE` for the all-doses fit) |
| `run_control_fasting.jl`, `control_fasting_report.jl` | long (120-day, `long_days`) dose + fasting control, read off at days 14/30/60/120 |
| `../../build_final_report.jl` | the Final Report from every chapter, opening with the modelling decisions; every version is archived in `docs/outputs/final_report/history/` (`../../report_history.jl`) |
| `decisions.toml` | dated log of modelling decisions, shown in the Final Report and its history |

## Tests

| Test | Checks |
|---|---|
| `test/test_data.jl` | the tidy table reproduces the source values of the earlier loaders (QC since 2026-10-07 may only add exclusions; day 0 is the seeding) |
| `test/test_stages.jl` | `stages.toml` rebuilds every candidate of the saved run (`DRYRUN=1`) |
| `test/test_validation_config.jl` | `validation.toml` selects the intended data |
| `test/check_intermittent.jl` | retired (2026-10-07): compared against the earlier scale-fitted intermittent predictions; the intermittent data are now predicted with the scale fixed at 1 (Validation page) |
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
