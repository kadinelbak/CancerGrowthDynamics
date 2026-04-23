# 01 Mechanical Manual Modelling

This folder contains the full manual/mechanistic modelling workflow migrated from the legacy `Modelling Data Notebooks` location.

## Scope
- Manual Julia notebooks for untreated/treated mono- and co-culture modelling.
- Associated CSV parameter tables and generated plots.
- Local Julia environments used by notebook workflows.

## Current Structure
- `BaseCode/`: shared setup scripts.
- `SimpleModels/`: baseline coculture mechanistic models.
- `Untreated MonoCulture/`: untreated monoculture fitting workflow.
- `Treated MonoCulture/`: treated monoculture fitting workflow.
- `Untreated CoCulture/`: untreated coculture fitting workflow.
- `Treated Coculture/`: treated coculture transfer/fitting workflow.

## Week 1 Migration Notes
- Source files were moved from the old path into this folder.
- Hard-coded source references were refactored from `Modelling Data Notebooks/...` to `Modeling_Approaches/01_mechanical_manual/...` where applicable.
- A compatibility symlink named `Modelling Data Notebooks` is kept at repository root to avoid breaking legacy references during transition.

## Run Guidance
1. Start with notebook-specific environment setup (where `Project.toml` exists).
2. Execute notebook cells top-to-bottom to regenerate outputs.
3. Use each notebook's `outputs/csv` and `outputs/images` folders for exported artifacts.

## Handoff to Next Approaches
Use outputs from this folder as the fixed baseline for:
- `02_mechanical_automatic_package`
- `03_automatic_neural_network`
