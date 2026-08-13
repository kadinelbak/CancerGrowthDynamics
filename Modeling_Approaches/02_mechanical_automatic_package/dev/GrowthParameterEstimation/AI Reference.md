# GrowthParameterEstimation.jl — AI Reference

This document explains the **main exported functions** in the package in:
- plain English (what each function does), and
- math form (the main equations/metrics used internally).

---

## 1) Data layer

### `load_timeseries(path; kwargs...)`
**English:** Reads a CSV and normalizes it to the package schema.

**Math:** No model equation; this is an I/O + preprocessing step:
\[
\text{DataFrame}_{raw} \xrightarrow{\text{normalize\_schema}} \text{DataFrame}_{standardized}
\]

### `normalize_schema(df; column_map, defaults)`
**English:** Renames mapped columns, fills missing required columns with defaults, enforces canonical column order/types.

**Math:** Applies deterministic transforms elementwise, e.g.
\[
t_i \leftarrow \mathrm{Float64}(t_i), \quad y_i \leftarrow \mathrm{Float64}(y_i)
\]
and for missing errors:
\[
\sigma_i \leftarrow 1.0 \quad \text{if missing.}
\]

### `validate_timeseries(df)`
**English:** Checks schema validity (required columns, finite numeric values, nonnegative counts, positive errors, monotone time per condition).

**Math:** Key constraints:
\[
y_i \ge 0, \quad \sigma_i > 0, \quad t_{j+1}-t_j \ge 0 \text{ (within condition)}.
\]

---

## 2) Exposure layer

### `build_exposure(kind; kwargs...)`
**English:** Factory that returns an exposure function object (`constant`, `pulse`, `stepped`, `decay`).

**Math:** Constructs a function \(E(t)\) based on `kind`.

### `evaluate_exposure(exposure, times)`
**English:** Evaluates the exposure profile at all requested times.

**Math:**
\[
\mathbf{E} = [E(t_1), E(t_2), \dots, E(t_n)].
\]

#### Exposure function families

- **Constant:**
  \[
  E(t)=c
  \]

- **Pulse:**
  \[
  E(t)=
  \begin{cases}
  A, & t_{start}\le t\le t_{end}\\
  0, & \text{otherwise}
  \end{cases}
  \]

- **Stepped:** piecewise constant by change points \(\tau_k\):
  \[
  E(t)=v_k \text{ for } \tau_k \le t < \tau_{k+1}.
  \]

- **Decaying:**
  \[
  E(t)=
  \begin{cases}
  0, & t<t_0\\
  c_0 e^{-\lambda (t-t_0)}, & t\ge t_0
  \end{cases}
  \]

---

## 3) Registry + simulation layer

### `register_model(spec; overwrite=false)`
**English:** Adds a `ModelSpec` to the in-memory model registry.

**Math:** No direct equation; stores model components:
\[
\{f(u,p,t,E),\; \text{bounds},\; h(u,p,t)\}
\]

### `get_model(name)`
**English:** Returns the `ModelSpec` for a registered model name.

### `list_models()`
**English:** Returns sorted registered model names.

### `simulate(spec, times, params; u0, exposure, reltol, abstol, enforce_nonnegative)`
**English:** Solves the model ODE and maps states to observed outputs.

**Math:** Integrates
\[
\frac{du}{dt}=f(u,p,t,E(t)), \quad u(t_0)=u_0
\]
then computes observations
\[
\hat y_i = h(u(t_i),p,t_i).
\]

If `enforce_nonnegative=true`, negative states are clamped post-solve.

---

## 4) Observation helpers

### `observed_signal(spec, state, p, t)`
**English:** Converts latent state to measured signal via a map and scale.

**Math:**
\[
\hat y = s\cdot g(u,p,t)
\]
where `s = spec.scale`, `g = spec.map_fn`.

### `viable_total(state, p, t)`
**English:** Uses first state component as observed value.

**Math:**
\[
\hat y = u_1.
\]

### `sum_states(indices)`
**English:** Returns a mapping function that sums selected state components.

**Math:**
\[
\hat y = \sum_{i\in \mathcal I} u_i.
\]

---

## 5) Core model ODE functions

Let \(N=u_1\).

### `logistic_growth!(du,u,p,t)`
\[
\frac{dN}{dt}=rN\left(1-\frac{N}{K}\right)
\]

### `logistic_growth_with_death!(du,u,p,t)`
\[
\frac{dN}{dt}=rN\left(1-\frac{N}{K}\right)-dN
\]

### `gompertz_growth!(du,u,p,t)`
\[
\frac{dN}{dt}=aN\log\left(\frac{K}{N}\right)
\]
(implemented with safety guards near invalid domains).

### `gompertz_growth_with_death!(du,u,p,t)`
\[
\frac{dN}{dt}=aN\log\left(\frac{K}{N}\right)-dN
\]

### `exponential_growth!(du,u,p,t)`
\[
\frac{dN}{dt}=rN
\]

### `exponential_growth_with_delay!(du,u,p,t)`
**English:** Logistic-style growth switched on after lag.

\[
\frac{dN}{dt}=\mathbf{1}_{t\ge t_{lag}}\;rN\left(1-\frac{N}{K}\right)
\]

### `logistic_growth_with_delay!(du,u,p,t)`
\[
\frac{dN}{dt}=\mathbf{1}_{t\ge t_{lag}}\;rN\left(1-\frac{N}{K}\right)
\]

### `exponential_growth_with_death_and_delay!(du,u,p,t)`
\[
\frac{dN}{dt}=\mathbf{1}_{t\ge t_{lag}}\;rN\left(1-\frac{N}{K}\right)-dN
\]

---

## 6) Fitting helpers

### `setUpProblem(model, x, y, solver, u0, p0, tspan, bounds; max_time, maxiters)`
**English:** Builds and solves a bounded optimization problem for parameters.

**Math:** Minimizes L2 loss between observations and model predictions:
\[
\min_p \sum_{i=1}^n \left(y_i-\hat y_i(p)\right)^2.
\]

### `calculate_bic(prob, x, y, solver, p)`
**English:** Computes SSR and BIC for fixed parameters.

**Math:**
\[
\mathrm{SSR}=\sum_{i=1}^n (y_i-\hat y_i)^2,
\]
\[
\mathrm{BIC}=n\log\left(\frac{\mathrm{SSR}}{n}\right)+k\log n,
\]
where \(k\) is number of fitted parameters.

### `pQuickStat(x, y, p, sol, prob, bic, ssr)`
**English:** Prints fitted parameters and fit metrics.

### `run_single_fit(x, y, p0; model, fixed_params, solver, bounds, max_time, show_stats)`
**English:** High-level single-model fit wrapper; returns `(params, bic, ssr, solution)`.

**Math:** Solves the same objective as `setUpProblem` and reports BIC/SSR above.

### `compare_models(...)`
**English:** Fits two models on same dataset and returns both fits + best by BIC.

**Math:**
\[
\text{best} = \arg\min_{m\in\{1,2\}} \mathrm{BIC}_m.
\]

### `compare_datasets(...)`
**English:** Fits (possibly same) model to two datasets and writes side-by-side metrics.

### `compare_models_dict(x, y, specs; ...)`
**English:** Batch-fits a dictionary of model specs, writes BIC/SSR summary and predictions.

### `fit_three_datasets(...)`
**English:** Fits one model to three datasets (named overload) or many datasets (vector overload), then summarizes fit quality/parameter variability.

**Math (multi-dataset summary overload):**
\[
\bar p_j = \frac{1}{M}\sum_{m=1}^M p_{m,j},
\quad
s_j = \sqrt{\frac{1}{M-1}\sum_{m=1}^M (p_{m,j}-\bar p_j)^2}.
\]

### `run_joint_fit(model, dataset_specs, u0, p0; solver, bounds, show_stats, maxiters)`
**English:** Fits one shared parameter vector across multiple observed datasets/states in a multi-state ODE.

Each dataset spec contains `(x, y, state_index)`.

**Math:** Uses a joint objective:
\[
\min_p \sum_{d=1}^{D}\sum_{i=1}^{n_d}\left(y_{d,i}-\hat y_{d,i}(p)\right)^2,
\]
where \(\hat y_{d,i}(p)\) comes from the state component indexed by `state_index`.

Reports joint BIC:
\[
\mathrm{BIC}_{joint}=n_{tot}\log\left(\frac{\mathrm{SSE}_{joint}}{n_{tot}}\right)+k\log(n_{tot}).
\]

### `compare_joint_models_dict(dataset_specs, u0, specs; default_solver, show_stats, output_csv)`
**English:** Batch-fits several joint models and writes a BIC/SSE summary CSV.

**Math:** Applies the same joint objective for each model and ranks by lower joint BIC.

---

## 7) Analysis helpers

### `leave_one_out_validation(x, y, p0; ...)`
**English:** Refit model \(n\) times, each time leaving one point out, then predict that left-out point.

**Math:** For each \(i\), fit on \(\mathcal D\setminus\{i\}\), predict \(\hat y_i\), then compute
\[
\mathrm{RMSE}=\sqrt{\frac{1}{n_v}\sum (y_i-\hat y_i)^2},
\quad
\mathrm{MAE}=\frac{1}{n_v}\sum |y_i-\hat y_i|,
\]
\[
R^2=1-\frac{\sum (y_i-\hat y_i)^2}{\sum (y_i-\bar y)^2}.
\]

### `k_fold_cross_validation(x, y, p0; k_folds, ...)`
**English:** Splits data into \(k\) folds, trains on \(k-1\), validates on held-out fold, aggregates errors.

**Math:** Same RMSE/MAE/\(R^2\) formulas over pooled validation predictions.

### `parameter_sensitivity_analysis(x, y, fit_result; perturbation, ...)`
**English:** Perturbs each parameter by \(\pm\delta\) and measures prediction change.

**Math:** For parameter \(p_j\):
\[
p_j^{\pm}=p_j(1\pm\delta),
\]
relative change in prediction:
\[
\Delta_{rel}(t)=\frac{|\hat y^{\pm}(t)-\hat y_{base}(t)|}{|\hat y_{base}(t)|+\varepsilon},
\]
sensitivity index:
\[
SI_j = \frac{\max_t \Delta_{rel}(t)}{\delta}.
\]

### `residual_analysis(x, y, fit_result; outlier_threshold, ...)`
**English:** Computes residual diagnostics, standardized residuals, outlier flags, and rough normality/autocorrelation checks.

**Math:**
\[
r_i=y_i-\hat y_i,
\quad
z_i=\frac{r_i}{s_r},
\]
outlier if \(|z_i|>\tau\).

---

## 8) Staged workflow + hardening (v0.3.0)

### `run_pipeline(df; config, models_to_compare, selection_mode, manual_choices, strict_schema, qc_before_fit, ...)`
**English:** End-to-end convenience workflow for condition construction, fitting, and ranking across all conditions.

**Math:** Per condition/model objective remains weighted least squares over trajectory mismatch:
\[
\min_p \sum_{i=1}^{n} w_i\left(y_i-\hat y_i(p)\right)^2,
\]
with model ranking typically based on BIC.

### `run_staged_pipeline(df; stages, selection_mode, manual_choices, strict_schema, qc_before_fit, n_bootstrap, resume_manifest_path, resume_from_stage, ...)`
**English:** Multi-stage orchestration where selected parameters/models are inherited by downstream stages. Supports fully automatic ranking or checkpoint/manual decisions at each stage.

**Math:** For stage index \(s\), downstream parameter vector can be decomposed as:
\[
p^{(s)} = \left[p^{(s)}_{free},\; p^{(s)}_{inherited}\right],
\]
where inherited components are fixed to estimates from prior stages.

### `default_population_stages(...)` and `default_population_cellline_stages(df; ...)`
**English:** Build biologically structured stage templates (for example, population-global growth terms and per-cell-line treatment potency terms).

### `validate_required_metadata(df)` and `validate_strict_schema(df)`
**English:** Fail-fast preflight checks for required metadata and strict schema consistency.

**Math:** Implements constraint checks such as domain-validity predicates:
\[
\forall i:\; t_i\in\mathbb{R}_{\ge 0},\; y_i\in\mathbb{R}_{\ge 0},\; \sigma_i\in\mathbb{R}_{>0}
\]
plus required metadata presence by row/condition.

### `generate_qc_report(df)` and `save_qc_report(report; output_dir=...)`
**English:** Summarizes data quality and condition coverage before model fitting.

### `save_run_manifest(...)` and `load_run_manifest(path)`
**English:** Persist and reload provenance metadata (config, selected models, inherited maps, timestamps, stage status) to make staged runs reproducible and resumable.

### `bootstrap_stage_uncertainty(stage_result; n_bootstrap, rng, ...)`
**English:** Estimates uncertainty around stage-level parameter estimates via bootstrap resampling.

**Math:** If \(\hat p^{(b)}\) is the estimate from bootstrap replicate \(b\), then summary moments are:
\[
\bar p = \frac{1}{B}\sum_{b=1}^{B}\hat p^{(b)},
\quad
\widehat{\mathrm{SE}}(p) = \sqrt{\frac{1}{B-1}\sum_{b=1}^{B}(\hat p^{(b)}-\bar p)^2}.
\]

Durbin–Watson:
\[
DW=\frac{\sum_{i=2}^{n}(r_i-r_{i-1})^2}{\sum_{i=1}^{n}r_i^2}.
\]

### `enhanced_bic_analysis(x, y; models, model_names, p0_values, solver, ...)`
**English:** Fits multiple candidate models and compares by BIC/AIC/AICc/\(R^2\), with model recommendation.

**Math:**
\[
\mathrm{AIC}=n\log(\mathrm{SSR}/n)+2k,
\]
\[
\mathrm{AICc}=\mathrm{AIC}+\frac{2k(k+1)}{n-k-1},
\]
\[
w_m \propto \exp\left(-\frac12\Delta\mathrm{BIC}_m\right), \quad \sum_m w_m=1.
\]

---

## 8) Workflow APIs

### `default_config(; output_dir="results")`
**English:** Creates default pipeline settings.

### `save_config(path, cfg)` / `load_config(path)`
**English:** Save/load pipeline config to/from TOML.

### `build_conditions(df; condition_cols=[:dose,:cell_line,:density,:replicate])`
**English:** Groups normalized data into `FitCondition`s (one per experimental condition).

**Math:** Partitions dataset into groups \(\mathcal D_c\) and sorts each by time.

### `fit(spec, conditions; shared_params, fixed_params, tie_constraints, ...)`
**English:** Jointly fits one model across multiple conditions with shared/fixed/tied parameter constraints.

**Math:** Minimizes global objective
\[
\min_\theta \sum_{c=1}^{C} \sum_{i=1}^{n_c} w_{c,i}\big(y_{c,i}-\hat y_{c,i}(\theta)\big)^2,
\]
with
\[
w_{c,i}=\frac{1}{\sigma_{c,i}^2}
\]
when `weighted=true`.

### `rank_models(model_names, conditions; top_k, kwargs...)`
**English:** Fits each candidate model and returns BIC-ranked summary + failure log.

**Math:** ranks by ascending BIC.

### `plot_topk(rank_result; conditions, top_k, output_dir)`
**English:** Writes top-model overlay CSVs/plots per condition + BIC bar chart outputs.

### `export_results(rank_result; output_dir)`
**English:** Exports ranking, best-fit starts, parameter table, failure report, and summary files.

### `run_pipeline(data_input; config, include_models, exclude_models)`
**English:** End-to-end execution: normalize/validate data → build conditions → rank models → plot/export outputs.

---

## 9) Practical notes

- Time vectors should be monotone for stable ODE fitting and diagnostics.
- Lower BIC/AIC is better; BIC weights approximate relative model support.
- If no model fits succeed in a workflow run, plotting/export may be skipped and failures should be inspected.
