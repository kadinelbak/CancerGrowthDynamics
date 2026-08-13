# Architecture Plan and Gap Closure Roadmap

## Core Goal
Build a reusable modeling package that ingests treated growth data, simulates multiple model families, fits consistently, compares fairly, and exports reproducible publication-ready outputs.

## Current Implemented Architecture

### 1) Data layer
- `normalize_schema`, `validate_timeseries`, `load_timeseries`
- Enforces canonical columns: `time,count,error,dose,cell_line,density,replicate,unit_time,unit_count`
- Handles defaults and finite/nonnegative checks

### 2) Drug exposure layer
- Unified callable exposure objects:
  - `ConstantExposure`
  - `PulseExposure`
  - `SteppedExposure`
  - `DecayingExposure`

### 3) Model layer
- `ModelSpec` registry and `register_model` interface
- Unified dynamics signature: `dynamics!(du, u, p, t, exposure)`
- Built-in families include baseline Hill, S/R, damage-repair, adaptive IC50, PK-PD, transit-chain, bi-exponential

### 4) Simulation layer
- Shared tolerances and aligned sampling at observation times
- Failure categorization (`solver_failure`, `biological_domain_failure`, `unknown_failure`)
- Nonnegative floor support

### 5) Fitting layer
- Multi-start optimization (`Fminbox(BFGS)`)
- Shared vs condition-specific parameter fitting across grouped conditions
- Fixed parameters and tie constraints
- Weighted objective support via error columns

### 6) Selection layer
- `rank_models` computes and ranks by SSE/weighted SSE/AIC/BIC and ΔBIC
- Keeps top-k fits per model

### 7) Visualization layer
- `plot_topk` creates overlay CSVs for each condition
- Auto-generates PNG overlays when `Plots.jl` is available

### 8) Export layer
- Stable folder layout via `export_results`:
  - `tables/`
  - `params/`
  - `diagnostics/`
  - `figures/`

### 9) One-command pipeline
- `run_pipeline`: load/normalize/validate -> build conditions -> fit/rank -> plot -> export

### 10) Reproducibility
- `PipelineConfig` with explicit version
- Deterministic seeds in fitting/ranking
- Config read/write via TOML

## Public API Surface
- `register_model`
- `simulate`
- `fit`
- `rank_models`
- `plot_topk`
- `export_results`
- `run_pipeline`

## Remaining Deepening Work (Production-Grade)
1. Add true DDE solver-backed models (current delay coverage includes transit-chain surrogate).
2. Add practical identifiability suite: profile likelihood, bootstrap, parameter-correlation heatmaps.
3. Add richer observation/noise models (state-space/likelihood variants).
4. Add richer diagnostics plots and theme-consistent publication styling defaults.
5. Add benchmark datasets and larger integration tests.
