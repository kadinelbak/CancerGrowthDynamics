# 03 Automatic Neural Network Modelling

This is a data-driven modelling track that does **not** use GrowthParameterEstimation.jl. It has three goals:
1. Predict next-day growth.
2. Recover an analytical form and feature importance.
3. Provide a differentiable, dose-controllable model for optimal control.

## Data (`src/data.jl`)
The loader rebuilds every 20k and 30k A2780/A2780cis series from the rawest source in the repo, then applies quality control in a fixed order. The full audit is in `outputs/data/report.html`.

**Sources:**
- **From raw per-image files:** untreated and treated monoculture (0, 0.62, 1 and 1.47 uM), intermittent (on→off and off→on at 1 uM), and untreated (ce0) and treated (ce1) co-culture.
  - Co-culture wells are assigned by plate column, per `Well Layouts`: columns 1–3 are 20k and 4–6 are 30k, holding 75/25, 50/50 and 25/75 in order.
  - The legacy co-culture `Averages` files split wells by imaging tile, which mixes 20k and 30k wells. They are not used.
- **From the LR workbooks:** run 1 (20k and 30k) and run 2 (30k). No raw images exist for these.
- **Not loaded:**
  - The New Datasets monoculture CSVs and the ce0/ce1 workbooks, because they duplicate the raw data.
  - 40k data, TykNu, the `New Datasets/Validation` files, and day 0 (a seeding value, not a measurement).

**Quality control, in order:**
1. **Well spikes:** a well-day is dropped if it jumps more than 2-fold against its own neighbouring days while the other wells don't.
2. **Plate-wide artefact days:** a day is dropped from a whole plate when the median dip across its series is below −0.25 (log).
3. **Single points:** a last point more than 2.5-fold off trend, or an interior point more than 2-fold off both neighbours.
4. **Day 8–9 drops:** day 8 or 9 is dropped from a series only when it dips more than 0.25 (log) below the line from day 7 to day 10.

**Density:** seeding density enters through each well's starting cell number, plus one growth-rate multiplier for 20k wells (`r20_over_r30`). Stage 2 reports its interval and a hold-out test that fits without any 20k wells, then predicts them.

**Splits:**

| Split | Contents |
|---|---|
| train: drug | monocultures, intermittent, LR monoculture |
| train: competition | untreated co-culture (ce0, LR co-culture), used only to fit the competition coefficients |
| **validation** | ce1 treated co-culture (held out) |

## Model
UDE equations (`src/ude.jl`):

```
dS/dt = r_S S (1 − (S + a_SR R)/K_S) − δ_S(t) S − kill_S(v, ρ) S
dR/dt = r_R R (1 − (R + a_RS S)/K_R) − δ_R(t) R − kill_R(v, ρ) R
dv/dt = γ (u(t) − v)                               u(t) = dose schedule
kill  = v · softplus(NN(v, ρ, line))               3-8-8-1 tanh network
```

- δ(t) is a delayed death term, active only in low-resource wells.
- log(N + floor) accounts for the image-analysis detection floor.
- The competition coefficients `a` are separate for normal and low resource.
- Comparators: a Hill kill (mechanistic), no drug, and a pure neural ODE.

## Pipeline
Run the scripts in order, from this folder. Use threads for 03 to 05.

| Script | Stage | Output |
|---|---|---|
| `julia --project=. scripts/run_01_data.jl` | 0: harmonised data and plots | `outputs/data/` |
| `julia -t 12 --project=. scripts/run_02_nextday.jl` | 1: next-day benchmark (persistence, ridge, boosted trees, MLP) | `outputs/nextday/` |
| `julia -t 14 --project=. scripts/run_03_ude.jl` | 2: UDE ensembles, hold-out folds, validation | `outputs/ude/` |
| `julia -t 12 --project=. scripts/run_04_symbolic.jl` | 3: sparse regression to a closed form, then refit as an ODE | `outputs/symbolic/` |
| `julia -t 4 --project=. scripts/run_05_importance.jl` | 4: permutation, TreeSHAP, UDE sensitivity and ablation | `outputs/importance/` |
| `julia -t 8 --project=. scripts/run_06_control.jl` | 5: optimal daily dose schedule | `outputs/control/` |
| `julia --project=. scripts/run_07_verdict.jl` | Verdict: neural vs mechanistic scorecard, with guidance on where each should be used | `outputs/verdict/` |

Each stage writes a `report.html`.

**Environment variables:**
- `NN_FAST=1` runs a quick smoke version of every stage.
- `NN_OUT=<dir>` redirects the outputs.

## Neural vs mechanistic
The neural models (MLP, a UDE with a neural kill term, and a pure neural ODE) are a **competing method** to the mechanistic models, not a step towards them. `outputs/verdict/report.html` scores the best of each on every test and explains where each should be used.

## Notes
- The fits use fixed-step RK4 with ForwardDiff gradients. That avoids the heavy SciML stack on Julia 1.12, and it makes the same model directly differentiable for the control problem.
- The sparse regression is an exhaustive best-subset search over a 14-term library, with a STLSQ (SINDy) sweep shown alongside.
- The optimal-control set-up mirrors `02_mechanical_automatic_package/src/optimal_control_one.jl`:
  - 13 daily doses
  - mean dose fixed at 1 uM, using the same projection
  - the objective is total cells at the final time
