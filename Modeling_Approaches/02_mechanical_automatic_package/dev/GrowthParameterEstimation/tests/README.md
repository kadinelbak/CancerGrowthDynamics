# Practice Notebook

This folder now contains maintained practice notebooks:

- `function_tour.ipynb`
- `pipeline_step_by_step_template.ipynb`

They provide a guided API walkthrough, a synthetic joint-fitting example, and a
step-by-step staged pipeline template.

For a script-based end-to-end template, see:

- `../examples/pipeline_one_shot_template.jl`

## Canonical automated tests

Automated package tests are maintained in:

- `test/runtests.jl`

Run from repository root:

```julia
julia --project=. test/runtests.jl
```
