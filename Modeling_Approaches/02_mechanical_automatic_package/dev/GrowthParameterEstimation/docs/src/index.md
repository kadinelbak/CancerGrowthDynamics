# GrowthParameterEstimation.jl

Current release: `v0.3.0`.

```@meta
CurrentModule = GrowthParameterEstimation
```

Welcome to the in-progress documentation for `GrowthParameterEstimation.jl`. The package provides ODE models and utilities for fitting, comparing, and analyzing growth dynamics. The documentation site is generated with [Documenter.jl](https://juliadocs.github.io/Documenter.jl/) directly from the package sources so it can track the public API as it evolves.

## What's New In v0.3.0

- Hardened staged workflows with strict schema checks, QC generation, and persisted manifests.
- Resume support for long-running staged pipelines.
- Bootstrap uncertainty estimation at stage level.
- Population and cell-line stage templates with parameter inheritance helpers.
- Coculture competition model support and simulation sweep utilities.

## Breaking Changes

- `v0.3.0` is a breaking pre-1.0 minor release relative to `v0.2.x`.
- Staged and workflow-oriented entry points now expect stricter canonical metadata and schema handling.
- Workflow exports now use the structured output layout documented in `CHANGELOG.md`.

See template assets in the repository for practical entry points:

- `examples/pipeline_one_shot_template.jl`
- `tests/pipeline_step_by_step_template.ipynb`

## Getting Started

```julia
using GrowthParameterEstimation

# choose a model and fit it to data
x = 0.0:1.0:5.0
y = [5.0, 8.0, 12.0, 18.0, 22.0, 28.0]
result = run_single_fit(collect(x), collect(y), [0.2, 50.0]; show_stats=false)
```

## API

```@autodocs
Modules = [GrowthParameterEstimation]
```
