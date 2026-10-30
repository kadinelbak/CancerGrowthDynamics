# 03 Automatic Neural Network Modelling

This is a data-driven modelling track that does **not** use GrowthParameterEstimation.jl. It has four goals:
1. Predict next-day growth.
2. Recover an analytical form and feature importance.
3. Provide a differentiable, dose-controllable model for optimal control.
4. Compare neural and mechanistic models fairly, on the same held-out data (Stages 6–7).

## Rules every model follows
- **Autonomous: no clock time.** No model has a time-dependent term, and no model or feature takes time as an input. The models drive optimal control, so they must be restartable from any state. Delays come only from states driven by the dose (exposure v) or by the cells (resource Q). Each line has at most 3 compartments: cells, v and Q.
  - Next-day features: `day`, `since_switch` (time since the last dose change) and `cumdose` (∫dose dt, which equals dose × day under a constant dose and never saturates, so it acts as a treatment clock) were removed on 2026-10-07. The dose history enters only through the filtered exposure v, which saturates and is a state of the input.
- **Starting values (user, 2026-10-07).** Day 0 is not a measurement and is never in the loss.
  - Every series starts at t = 0 from its nominal seeding: 100 cells per 30k well and 66.7 per 20k well, times the line's share in a coculture (25:75 gives 25 and 75).
  - The starting value is estimated (a bounded free parameter, counted in the BIC) only when the first kept count (day d, normally 1) is more than 3^d times above or below nominal. This is the zoo's rule (`zoo_data.jl` `x0_info`, bounds [n0/10, 3·max(n0, first count)]).
  - In the current data this applies to 5 high-dose or untreated A2780 monocultures, every low-resource series and every treated A2780 coculture series.
  - **Held-out series always start from the nominal seeding**, as in the zoo's `validation.toml`. Held-out data never set a starting value.
- **Intermittent values are cell counts** (already normalised to cell size), on the same scale as the other data. There is no per-trajectory offset.
- **Model choice:** among the candidates within 2 BIC of the lowest, the one with the fewest free parameters is kept.
- **Seeds:** every neural fit is repeated from `NN_SEEDS` (default 5) random initialisations, and so is every hold-out fold. There are two runtime reductions, both made on 2026-10-07 because the full rerun shared the machine with the mechanistic refit: the input-ablation retrains of Stage 4 use 3 seeds (`NN_ABL_SEEDS`), and the learned-resource UDE is not refitted in the hold-out folds. Scores are reported as mean ± SD over seeds and for the seed-ensemble prediction. Stochastic next-day models (boosted trees, MLP ensemble) are also refitted with 5 seeds. Hill, the closed form and no-drug are deterministic single fits; an 8-member Bayesian bootstrap of Hill and of the closed form gives parameter intervals.
- **Info tips:** every table and figure on the pages has an "i" tip that says how to read it and which direction is better. `html_table` and `figure` in `src/report.jl` refuse to build without one.

## Data (`src/data.jl`)
The loader rebuilds every 20k and 30k A2780/A2780cis series from the rawest source in the repo, then applies quality control in a fixed order. The full audit is in `outputs/data/report.html`.

**Sources:**
- **From raw per-image files:** untreated and treated monoculture (0, 0.62, 1 and 1.47 uM), intermittent (on→off and off→on at 1 uM), and untreated (ce0) and treated (ce1) co-culture. Co-culture wells are assigned by plate column (columns 1–3 = 20k, 4–6 = 30k; 75/25, 50/50, 25/75).
- **From the LR workbooks:** run 1 (20k and 30k) and run 2 (30k).
- **Stage 6 only:** the mechanistic zoo's tidy table `02_mechanical_automatic_package/data/measurements.csv` (`load_tidy`). It supplies the 30k monocultures, Ce0 and Ce1 with exactly the zoo's points and QC flags.
- **Not loaded:** the duplicate New Datasets files, 40k data, TykNu and `New Datasets/Validation`.
- **Known difference:** in the 30k Ce0 75S25R resistant series, the workbook (and so the tidy table) reads about 0.8× the raw-image rebuild. Stage 6 uses the workbook values, so its training data match the zoo's.

**Quality control, in order:** well spikes, plate-wide artefact days, single-point outliers, then day 8–9 drops. Details and counts are in the Stage 0 report.

**Splits:**

| Split | Contents |
|---|---|
| train: drug | monocultures, intermittent, LR monoculture |
| train: competition | untreated co-culture (ce0, LR co-culture), used only for the competition coefficients |
| **validation** | ce1 treated co-culture (held out) |

## Model (`src/ude.jl`)
```
dS/dt = r_S S (1 − (S + a_SR R)/K_S) − d_S(Q) S − kill_S(v, ρ) S
dR/dt = r_R R (1 − (R + a_RS S)/K_R) − d_R(Q) R − kill_R(v, ρ) R
dv/dt = γ (u − v)                      drug exposure, driven by the dose u
dQ/dt = −c Q (S + R) / 3500            resource left (low-resource wells only; Q = 1 in normal medium)
d_L(Q) = δ_L (1 − Q)                   low-resource death grows as the resource is used up
kill  = v · softplus(NN(v, ρ, line))   3-8-8-1 tanh network (UDE)
```
- **Competitors:**
  - **Neural:**
    - the UDE above;
    - a UDE that also learns the resource death;
    - a pure neural ODE with per-capita growth = NN(ρ, v, line, low, Q). It has no competition term, so it is fitted to the drug and untreated co-culture units together.
  - **Mechanistic:** a Hill kill, the closed form recovered in Stage 3, and no drug.
- The fits use fixed-step RK4 with ForwardDiff gradients (Adam warm-up, then L-BFGS), on log counts with σ = √(0.15² + (SEM/mean)²).

## Pipeline
Run from this folder. Long stages take `-t 6` (6 threads), so a concurrent mechanistic refit keeps its cores. `scripts/run_all.ps1 [-Fast] [-From n] [-Threads n] [-Out dir]` runs everything in order and logs to `outputs/pipeline_run.log`.

| Script | Stage | Output | Runtime (6 threads, 2026-10-07, sharing the machine with a 16-thread mechanistic refit) |
|---|---|---|---|
| `julia --project=. scripts/run_01_data.jl` | 0: harmonised data, QC, starting values | `outputs/data/` | ~1 min |
| `julia -t 6 --project=. scripts/run_02_nextday.jl` | 1: next-day benchmark (persistence, ridge, boosted trees, MLP; 5 seeds) | `outputs/nextday/` | 5 min |
| `julia -t 6 --project=. scripts/run_03_ude.jl` | 2: UDE / UDE + NN resource / pure neural ODE (5 seeds each), Hill (+ 8-member bootstrap), 4 hold-out folds (UDE and pure NODE × 5 seeds, Hill, no drug), validation. Finished fits are checkpointed in `outputs/ude/checkpoints/`; `NN_RESUME=1` reuses them | `outputs/ude/` | ~2.5 h (fits; ~5 min when resumed from checkpoints) |
| `julia -t 6 --project=. scripts/run_04_symbolic.jl` | 3: sparse regression to a closed form, refit as an ODE | `outputs/symbolic/` | 40 min |
| `julia -t 6 --project=. scripts/run_05_importance.jl` | 4: permutation, TreeSHAP, UDE sensitivity, input ablation (3 seeds per variant) | `outputs/importance/` | 1 h 50 min |
| `julia -t 6 --project=. scripts/run_06_control.jl` | 5: optimal daily dose schedule | `outputs/control/` | 4 min |
| `julia -t 4 --project=../02_mechanical_automatic_package scripts/export_mech_heldout.jl` | exports the zoo's held-out Ce1 prediction point by point (reads the zoo only) | `outputs/comparison/mech_heldout_*.csv` | ~1 min |
| `julia -t 6 --project=. scripts/run_07_comparison.jl` | 6: three-way comparison on the same held-out points | `outputs/comparison/` | 12 min |
| `julia -t 6 --project=. scripts/run_08_learning_curve.jl` | 7: learning curve, UDE vs Hill / closed form at 25/50/75% (4 random subsets each) and 100% (5 UDE seeds, reused from Stage 2) of the training wells | `outputs/learning_curve/` | 37 min |
| `julia -t 6 --project=. scripts/run_03b_folds.jl` | 2b: hold-out folds of Stage 2 (Hill, no drug, UDE and pure NODE × 5 seeds; run after Stage 7) | `outputs/ude/folds.html`, `holdout_folds.csv` | ~1.5 h (62 min in the final run plus 18 fits carried over from an interrupted run) |
| `julia --project=. scripts/run_09_verdict.jl` | 8: held-out scorecard, neural vs mechanistic | `outputs/verdict/` | <1 min |
| `julia scripts/publish_docs.jl` | copies the pages and figures to `docs/outputs/nn_track/` | | seconds |

Each stage writes a `report.html`.

**Environment variables:**
- `NN_FAST=1` runs a quick smoke version of every stage.
- `NN_OUT=<dir>` redirects the outputs.
- `NN_SEEDS=<n>` sets the number of seeds (default 5).
- `NN_LC_SUBSETS=<n>` sets the random subsets per learning-curve fraction (default 4).
- `NN_ABL_SEEDS=<n>` sets the seeds per ablation variant (default 3).
- `NN_RESUME=1` reuses the Stage 2 checkpoints.

## Stage 6: the fair comparison
- **Held out:** every model predicts the same 84 points, the 30k Ce1 treated cocultures (6 trajectories, days 1–14) of the zoo's tidy table, from the nominal seeding.
- **Training:** every model trains on the zoo's own training data for that prediction: the 30k monocultures and Ce0.
- **Models:** (a) the zoo's carried prediction, per point from `export_mech_heldout.jl`; (b) this track's Hill and closed form; (c) the UDE; (d) the autonomous pure neural ODE.
- **Metrics:** RMSE in cells, nRMSE, log-RMSE, and paired-bootstrap 95% CIs over trajectories for the differences against the best model. The page also lists the remaining differences between the arms: the loss scale, the use of Ce1 in the zoo's model choice, what is refitted on Ce0, and the detection floor.

## Neural vs mechanistic
The neural models are a **competing method**, not a step towards the mechanistic model. `outputs/verdict/report.html` scores only held-out tests. In-sample fits and BIC are shown but never counted. A difference within the noise (below 0.02 log units or the seed SD, or a bootstrap CI that includes 0) is a tie.
