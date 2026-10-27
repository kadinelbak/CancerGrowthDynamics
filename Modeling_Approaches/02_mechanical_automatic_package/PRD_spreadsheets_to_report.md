# PRD: Spreadsheets in → fits, qualification and control graphs out

**Owner:** Kadin El Bakkouri · **Implementer:** Fable (orchestrating `coder`/`scout` subagents) · **Date:** 2026-10-10 (rev. 2)
**Repo:** `C:\Research\CancerGrowthDynamics` (GitHub `kadinelbak/CancerGrowthDynamics`), work from `main`
**Package dir (PKG):** `Modeling_Approaches/02_mechanical_automatic_package`

> **Rev. 2 decisions (owner, 2026-10-10)** — these override anything else in this document:
> 1. **Any number of cell lines and any number of treatments** in v1 (not just a sensitive/resistant pair and one drug).
> 2. **Onboarding must be AI-free by default** and as easy as possible: drop files in, fill in a pre-filled tag
>    sheet in Excel, run. An AI assist may exist but is optional and must produce the same tag sheet.
> 3. **Coculture is optional.** Many studies will be monoculture only.
> 4. **Users supply their own EC50 values** (per line, per treatment).
> 5. The output folder (`runs/`, §5.6) is **git-ignored**; a user never needs to know about it beyond "your
>    results are in this folder".
> 6. **The A2780 study keeps its own configuration and publish step** (§9). The generic pipeline must never
>    change its numbers or its website.
> 7. **The generic pipeline is a personal project for other people to use, made outside the research lab.**
>    It must be separable from the lab's study:
>    - **Separate code.** Put the generic code (`cgd.jl`, the tagging, reader, validator, stage-plan, report and
>      control code it needs) in a self-contained folder, `study_pipeline/`, that could become its own repository
>      without edits.
>    - **One-way dependency.** The A2780 study in this repo uses that folder. The folder never depends on the
>      A2780 study.
>    - **No lab data shipped.** The demo study and all test fixtures are **synthetic**. Lab spreadsheets (A2780,
>      TykNu) are used only in tests that run inside this repo, and are skipped if the files are absent.
>    - **No lab specifics.** No lab names, branding or A2780-specific text in the generic code, templates or
>      report.
>    - **Add a LICENSE file** to `study_pipeline/` (owner to choose; default MIT) and a README written for a
>      first-time user.

---

## 1. Problem

The pipeline already runs from one command (`julia run_pipeline.jl`), is configuration-driven *for one dataset*
(data selection, stages, hold-outs and control settings are TOML files), is seeded, and `test/test_regression.jl`
reproduces the 2026-10-03 reference (31/31).

**It is not usable by someone with a different set of spreadsheets.** Audit of the current code:

| Barrier | Where (examples) |
|---|---|
| Exactly two lines named `A2780`/`A2780cis` hard-coded; "line 2 is the resistant one" assumed | `zoo_data.jl:35,76,93-94,177,190`, `run_zoo_fits.jl:177-178,202,242`, `zoo_analyses.jl:392,577,891`, `run_control_fasting.jl:70`, `control_model.jl:9,39` |
| One drug (cisplatin, µM) and one second input (low resource = FBS %) hard-wired | `zoo_models.jl`, `zoo_data.jl:159` `ur_of(fbs)`, `zoo_common.jl:214-216`, `control_fasting.jl:393-399` |
| Condition names (`untreated`, `low_resource`, `treated`, `Ce0`, `Ce1`, `combination`, `intermittent`) are fixed data structures | `zoo_data.jl:80-121` |
| Stage ids `s1…s4` and their inheritance re-encoded in code; Stage 4 is a bespoke path | `zoo_analyses.jl:72-73,365,463,729`, `zoo_selection.jl:82`, `zoo_stage4.jl`, `export_zoo_bundle.jl`, `control_model.jl:179,197`, `validation.toml` group ids |
| Ratio notation `25S75R` and sheet suffix map (`1-3` = 25:75) hard-coded | `zoo_data.jl:35,97-103` |
| Seeding → starting count (100 at 30k) is a constant | `zoo_data.jl:28-37`, `validation.toml:191`, `control_config.toml:101` |
| 14-day horizon, dose levels, axis label "Mean Cells", report text | `zoo_common.jl:225`, `control_config.toml:30-32`, `zoo_report.jl:23`, `control_report.jl:51` |
| Reader: one series per file/sheet, header in row 1, column literally `Day` | `data/build_tidy.jl:23-34,88-91` |
| Metadata (line, dose, ratio, run, seeding) exists only in file/sheet names and is hand-typed into `sources.toml` | `New Datasets/*` |
| No input validation; the "checks" step compares to legacy A2780 numbers and fails on any new data | `test/test_data.jl` |
| Outputs written into the repo's `outputs/` and the website `docs/`; Final Report has 9 hard-coded chapters | `run_zoo_fits.jl:36`, `zoo_analyses.jl:33`, `control_report.jl:7`, `build_final_report.jl`, `stage_names.jl` |

A real second study already sits in the repo and cannot be run today: `New Datasets/ALL Tyknu-cpr LR data.xlsx`
(TykNu / TykNu-cpr, same layout as the A2780 low-resource workbook).

## 2. Goal

> A user who has never seen this repo installs Julia, drops a folder of spreadsheets into `studies/<name>/raw/`,
> runs one command, fills in the tag sheet it generates (most cells pre-filled), runs the command again, and gets a
> folder with fits, a qualification scorecard and optimal-control graphs, plus one HTML report. Re-running on the
> same inputs gives identical numbers. **No AI, no account, no API key is required.**

### Acceptance tests (all must pass)

1. **A2780 parity.** The current A2780 dataset run through the new study path reproduces
   `test/reference_results.toml` (31/31), also on 4 threads. The general N-line / N-treatment code with N = 2 lines,
   1 drug + low resource must give the *same numbers* as today.
2. **TykNu, zero-code.** `ALL Tyknu-cpr LR data.xlsx` → `init` → tag sheet (user only confirms guesses / types
   EC50 if any treated data) → `run` → report. Only the stages its data support run; the rest are listed as
   "skipped: no data" with the reason.
3. **Monoculture-only, many lines, many treatments.** A synthetic study generated by a test script: 4 lines,
   2 drugs (one with a combination arm), 3 doses each, no coculture, in a *deliberately different* spreadsheet
   layout (wide sheet, header on row 3, column `Time (h)`). After tagging, fits recover the generating parameters
   inside their bootstrap 95% intervals and the carried model is the generating structure.
4. **Coculture with 3 lines.** A synthetic study with a 3-line coculture recovers the generating interaction
   structure (or a model within the 2-BIC rule of it).
5. **Bad input rejected before fitting**, with file / sheet / row messages: missing time column, text in numeric
   cells, non-monotone or duplicated days, duplicate series, coculture well with a line missing, tag sheet row
   with an unknown line or treatment, treated data with no EC50 and too few doses to estimate one.
6. **Determinism.** Two runs of the same study (same commit, same Manifest) produce byte-identical result CSVs
   (timestamps and absolute paths only in the manifest).
7. **Clean clone, Windows and macOS/Linux.** `git clone` → install Julia (pinned version, via juliaup) →
   `cgd init studies/demo` → `cgd run studies/demo` works with no other manual step. A bundled tiny demo study
   (`studies/demo/`) runs in under 10 minutes.

## 3. User flow

```
1. put spreadsheets in studies/<name>/raw/         (any mix of .xlsx / .csv, any layout from §5.2)
2. cgd init studies/<name>                         deterministic scan, no AI
      → studies/<name>/tags.xlsx                   one row per data series, guesses pre-filled (highlighted),
                                                    dropdowns for every categorical column, + Lines / Treatments /
                                                    Schedules / Study sheets
3. open tags.xlsx in Excel / LibreOffice, fix highlighted cells, fill EC50s, save
4. cgd run studies/<name>                          validates the tags first; stops with plain-language errors
      → studies/<name>/results/<run_id>/           fits, qualification, control, report.html
```

- `cgd` is a thin wrapper (`cgd.bat` on Windows, `cgd` shell script elsewhere) around
  `julia --project=<PKG> cgd.jl ...`. Double-clicking `cgd.bat` with no arguments asks for a study folder
  (simple console prompt) and runs `init` or `run` depending on whether `tags.xlsx` is complete.
- Re-running `init` on a study whose `tags.xlsx` exists **never overwrites user edits**: it only appends rows for
  new files and marks rows whose file disappeared.
- Users who already have tidy data can skip tagging by using the long-format template (§5.3).

## 4. Scope

**In scope (v1)**
- Any number of lines (populations), each with a name and optional role labels (e.g. "resistant") used only for
  reporting and control weights — never for model structure.
- Any number of treatments. A treatment has a *kind*:
  - `drug` — concentration input, Hill-type kill/inhibition with per-line EC50 (the current cisplatin modules,
    generalised);
  - `resource` — a stress / deprivation input scaled 0–1 between a "normal" and a "minimum" level (the current
    low-resource modules, generalised from FBS %).
- Any doses, schedules (constant or piecewise steps), seedings, runs, horizons, time units.
- Monoculture-only studies; cocultures of any subset of lines; combinations of any treatments.
- AI-free tagging (§5.1); optional AI assist (§5.4).
- Self-contained results folders; one generic study report; the A2780 website/Final Report kept as a separate
  `publish` step.

**Out of scope (v1)**
- Per-well raw rows as model input (onboarding may aggregate them to mean ± SD per day; no per-well model).
- New model *families* beyond generalising the existing zoo; the NN/UDE track; a graphical app.
- Pharmacokinetic drug–drug interaction models beyond what the combination candidates in §5.5 express.

## 5. Requirements

### 5.1 Tagging system (AI-free, primary onboarding path)

**Why a spreadsheet:** every user already has Excel or LibreOffice; dropdowns make it hard to type something
invalid; it needs no install, account or key; it is reviewable and committable; it is the single source of truth
for what each file means.

`cgd init` writes `tags.xlsx` with these sheets (column order fixed, documented in a `README` sheet inside the
workbook):

**`Series`** — one row per data series found in `raw/`:

| column | meaning | pre-fill from |
|---|---|---|
| `include` | yes/no | yes |
| `file`, `sheet`, `column` | where the series is (`column` only for wide sheets) | scan |
| `time_column`, `value_column`, `sd_column`, `sem_column`, `n_column` | which columns hold what | header matching (`day`, `time`, `hour`, `mean`, `count`, `cells`, `sd`, `sem`, `n`…) |
| `header_row`, `time_unit` | layout | scan; `h` if column says hour |
| `line` | which line this series counts (dropdown from `Lines`) | name tokens matched against `Lines` |
| `well_group` | series counted in the same wells share a value (coculture; blank = monoculture) | same file + same ratio/run tokens |
| `seeding_total` | cells seeded per well (all lines) | `30k`, `20k` tokens |
| `seeding_fraction` | this line's share at seeding (1 for mono) | ratio tokens (`25S75R`, `1-3`, `25Sensitive`) |
| `treatment_1 … treatment_k` | dose of each treatment (one column per row of `Treatments`; 0/blank = none) | `1um`, `IC50`, `0.5% FBS`… |
| `schedule` | blank = constant dose from t=0; else a name from `Schedules` | — |
| `run` | experiment repeat | leading `1 -`, `2 -`, `run2` tokens |
| `day0_is_seeding` | is the first row the nominal seeding, not a measurement | yes if first SD = 0 or value ≈ nominal |
| `notes` | free text, ignored | — |

**`Lines`** — `name`, `role_label` (optional), and one `EC50_<treatment>` column per drug-kind treatment
(units of that treatment). Pre-filled with line names guessed from tokens; EC50 cells left blank.

**`Treatments`** — `name`, `kind` (`drug` | `resource`), `unit`, for `resource`: `normal_level`, `min_level`.

**`Schedules`** — `name`, then `t_start`/`dose` pairs (piecewise-constant). Pre-filled for any schedule-like tokens
(`on0_off7`).

**`Study`** — `title`, `time_unit`, `count_unit`, `reference_seeding`, `reference_count` (count at day 0 for a
monoculture at `reference_seeding`; default from the data), `x0_estimate_fold` (default 3), `data_days`.

Rules:
- Every guessed cell is highlighted yellow with a comment saying *what it was guessed from*
  ("dose 1.47 from filename `A2780_IC25_1.47um_30k.csv`"). Unguessable required cells are red. `run` refuses
  while any red cell remains and lists them.
- All categorical columns use Excel data validation (dropdowns) referencing the `Lines`, `Treatments`,
  `Schedules` sheets.
- `tags.xlsx` is the only file a normal user edits. Internally it compiles to today's TOML configs
  (`sources.toml`, `schedules.toml`, a generated stage plan) under `results/<run_id>/inputs/`, so the existing
  fitting code keeps its configuration model. A CSV export (`tags/*.csv`) is written alongside for diffing in git.
- **Recommended file naming** is documented (e.g. `<line>_<treatment><dose><unit>_<seeding>k_run<k>.csv`); files
  named this way pre-fill completely, but no naming convention is ever *required*.
- The token parser must handle the known quirks of the current data: BOM and trailing space in `Day `, extra
  blank-header columns, sheet names like `2 - A2780cis 30k 3-1` (run 2, line, seeding, ratio code), coculture lines
  in separate sheets of one workbook (`25Sensitive_Untreated` / `75Resistant_Untreated`), day-0 rows that are the
  nominal seeding (`99.8`, `24.95`), intermittent values already normalised to cell size (do not divide by 144).

### 5.2 Accepted spreadsheet layouts (reader)

Extend `data/build_tidy.jl` (or a new reader it calls) to accept, via the tag sheet's columns:
- one series per file or per sheet (today);
- **wide**: one time column and several series columns in one sheet (one `Series` row per column);
- **long**: one sheet with tag columns + time + value (one `Series` row per distinct tag combination, or the
  canonical template, §5.3);
- header on any row; any time-column name; time in days or hours (converted to days); extra/blank columns ignored.

### 5.3 Canonical long template (skip tagging)

Ship `templates/study_template.xlsx` with the `Lines`, `Treatments`, `Schedules`, `Study` sheets of §5.1 and a
`Data` sheet with columns
`line, well_group, run, seeding_total, seeding_fraction, <one column per treatment>, schedule, time, mean, sd, sem, n_wells`.
A study containing only this workbook needs no `init`.

### 5.4 Optional AI assist (not required, never in the numerical path)

`cgd init --assist`: if the `claude` CLI is on PATH, after the deterministic scan it is asked (headless,
`claude -p`, prompt file in the repo, read-only access to `raw/`) to fill **only the cells the scanner left red or
yellow**, writing its reasoning into the cell comments. Output is the same `tags.xlsx`; cells it filled are
highlighted blue so the user reviews them. It must never fill EC50s, doses or units it cannot cite from the files.
If `claude` is absent or fails, `init` behaves exactly as without `--assist`. Provenance (model id, prompt SHA) is
recorded in `tags.xlsx`'s README sheet and the run manifest. This is a convenience; acceptance tests 2–5 must pass
without it.

### 5.5 Generalising the models and the stage plan (the core engineering work)

**Principle:** model structure depends on *which data exist*, never on names. After this work the strings
`A2780`, `A2780cis`, `cisplatin`, `Ce0`, `Ce1`, `25S75R`, `fbs`, and stage ids appear only in `studies/a2780/`
and in A2780 publish-only code. Add a test that greps the fitting/analysis/control code and fails otherwise.

**Auto-generated stage plan** (from the compiled tags; a study may supply an explicit `stages.toml` override —
the A2780 study does, to stay identical to the reference):

| Stage (generic name) | Runs when the data contain | Per | Inherits |
|---|---|---|---|
| Growth | untreated monoculture of a line | line | — |
| Single-treatment response | monoculture of a line under exactly one treatment (any doses/schedules, fitted jointly) | line × treatment | Growth |
| Untreated coculture | wells with ≥ 2 lines and no treatment | well-group line set | Growth of each line |
| Treated coculture | wells with ≥ 2 lines under one treatment | line set × treatment | Growth + that Single-treatment response of each line + Untreated coculture (if present, else fit interactions here) |
| Combination | any series with ≥ 2 treatments active | line (or line set) × treatment set | the Single-treatment responses involved (+ coculture stages if a coculture) |
| Schedule validation (not fitted) | series whose schedule is held out in `validation` settings | — | the carried models |

- A stage whose required inheritance is missing is **skipped** with a reason (e.g. Combination for line X skipped:
  no single-treatment data for treatment B). A line with no untreated data can still be fitted if the user marks a
  growth law in the `Lines` sheet; otherwise its stages are skipped.
- Candidate catalogues stay those of `zoo_models.jl` (`DRUG_ZOO` for `drug`, `LR_ZOO` for `resource`), keyed by
  treatment *kind*. Combination candidates: independent action (additive kills), Bliss-type multiplicative
  survival, and the current Stage 4 interaction (χ) term; selected by the 2-BIC rule.
- **Coculture with N lines:** Lotka–Volterra competition with shared K and an N×N crowding matrix. Candidates are
  generated in increasing complexity — one-for-one crowding; one shared α; per-line α (row-wise); full pairwise α
  (only when ≥ N² informative environments, else not generated); plus the existing coculture-death and refit-kill
  variants. The 2-BIC parsimony rule selects. Generation is capped (config `max_coculture_candidates`) and the cap
  is reported.
- `joint_rule` ("shared" / "ab") becomes a per-stage setting in the generated plan, default "ab" (shared structure
  tested, kept only if it wins), instead of being keyed by stage id in `zoo_selection.jl:82`.
- Stage 4's bespoke path (`zoo_stage4.jl`, `run_zoo_stage4.jl`) is folded into the generic Combination stage, or
  kept as the implementation of that stage but parameterised by the treatments/lines it combines.
- Validation (`validation.toml`) is generated per stage: forecast on every fitted group; leave-one-out over every
  tag with ≥ 3 values in that stage (dose, ratio, run, line set); held-out condition and schedule-prediction checks
  only when the user marks series `holdout = yes` in the tag sheet (add that column). Entries with no data are
  skipped with a note.
- **Constants → config:** `XMAX`/plot horizon → `Study.data_days`; nominal starting count →
  `reference_count × seeding_total / reference_seeding × seeding_fraction`; `ur_of(fbs)` →
  `(normal_level − level)/(normal_level − min_level)` per resource treatment; dose levels, `u_max_data` and policy
  names → derived from the doses present; "Mean Cells" → `Study.count_unit`; colours by line index from a
  colour-blind-safe palette.

### 5.6 EC50s

- User-supplied in the `Lines` sheet, one per (line, drug-kind treatment), in that treatment's unit. Used as fixed
  values exactly as today's per-line EC50.
- If blank: estimate it (one extra free parameter, counted in the BIC, flagged "estimated" in the report) **only if**
  that line has ≥ 3 distinct non-zero doses of the treatment; otherwise `run` stops with a message naming the line
  and treatment.
- A sensitivity option `ec50_mode = fixed | estimate` per treatment in the `Treatments` sheet lets a user estimate
  even when a value was given, and the report then shows given vs estimated.
- Never default two lines to the same EC50 silently.

### 5.7 Optimal control (generalised)

- Runs when at least one `drug` treatment has a fitted single-treatment response. Controls are the drug
  treatments (and `resource` treatments if the user marks them `controllable = yes` in `Treatments`); the plant is
  the carried model of the richest stage available: coculture if present, else **each line's monoculture model
  summed** (no interaction; stated in the report).
- **Objective = a literature-based preset, chosen by the user** in the `Study` sheet (`control_objective`
  dropdown). Each preset's weights are labelled LITERATURE or ASSUMPTION and can be overridden in
  `studies/<name>/control.toml`. Each preset is documented in the report with its citation:
  - `burden` (default when no line is labelled resistant): minimise integrated + terminal total burden with a
    quadratic dose penalty and the host-toxicity constraint (Swan 1990, *Math Biosci* 101:237; de Pillis &
    Radunskaya 2001, the toxicity model already in `control_config.toml`).
  - `time_to_progression` (default when any line is labelled resistant): the current A2780 objective (total
    burden + time-to-progression + terminal resistant cost), following adaptive-therapy reasoning (Gatenby et al.
    2009, *Cancer Res* 69:4894; Zhang et al. 2017, *Nat Commun* 8:1816).
  - `custom`: per-line weights (`control_weight` in `Lines`) and term weights set by the user.

  Presets must not change the A2780 study, which keeps its own `control_config.toml`. Before relying on a
  citation, the implementer checks it and fixes the reference if it is wrong.
- Dose levels per drug default to the measured doses + 0; `u_max` defaults to the highest measured dose; anything
  above is "demonstration / extrapolated" as today.
- Search cost grows with the number of controlled treatments: blocked search stays exact for one or two controls;
  above that, use the existing continuous SQP warm-started from the measured-dose plan, and say so in the report.
- Settings labelled `ASSUMPTION` in `control_config.toml` become study-overridable defaults
  (`configs/defaults.toml` ← `studies/<name>/control.toml`, optional).

### 5.8 Outputs: the results folder

A plain-language explanation goes into the README: *"Each time you run a study, its results go into a new folder
`studies/<name>/results/<date-time>/`. Old results are never overwritten. Open `report.html` in that folder."*
(This replaces the `runs/` naming of rev. 1.)

```
studies/<name>/results/<run_id>/     run_id = <UTC yyyymmdd-hhmm>-<git short sha>; git-ignored
  report.html        start here: data overview, fits, qualification scorecard, control, skipped stages
  inputs/            copies of raw files, tags.xlsx, the compiled TOML configs
  data/              measurements.csv (tidy; exclusions kept with reasons), validation_report.md
  fits/              per stage: candidates, BIC, carried model, parameters + bootstrap CIs; figures (PNG + SVG)
  qualification/     scorecard.csv + figures
  control/           schedules, trajectories, comparisons, figures (if run)
  manifest.toml      commit, dirty flag, Julia + Manifest pin, threads, seeds, SHA-256 of every input, step
                     timings, time-limit hits, skipped stages + reasons, AI-assist provenance if used
  logs/
```

- `run` writes nothing to `PKG/outputs/` or `docs/`. `cgd publish a2780 <run_id>` copies into `docs/`, runs
  `stage_names.jl`, `export_zoo_bundle.jl`, `build_final_report.jl` (A2780-only, as today).
- `--from`, `--only`, `--resume`, hold files carry over, scoped to a results folder.
- The A2780 reference run's `manifest.toml` and scorecard are committed under `studies/a2780/reference/`.
- Plots follow `CancerGrowthDynamics/plotting_guidelines.md`: data with error bars, fitted band, excluded points
  marked not hidden, axes in the study's units.

### 5.9 Validation of inputs (deterministic, before any fit)

`run` first compiles and validates the tags and writes `data/validation_report.md`.
**Errors (stop):** red cells left; file/sheet/column missing; text in numeric cells; negative, non-monotone or
duplicated times; duplicate series; unknown line/treatment/schedule; a `well_group` with one line, or its lines
under different treatments/schedules; seeding fractions of a well group not summing to 1 (±0.01); drug data with
no EC50 and < 3 doses.
**Warnings (continue, shown in report):** < 5 time points; SD = 0 after day 0; first measured point > 3-fold off
nominal (x0 will be estimated); QC exclusions > 20% of a series; treated series not below untreated at the same
seeding by the last day.
QC thresholds (`qc_rules.toml` `floor` = 50 cells) scale with `reference_count`.
`test/test_data.jl` is split: generic validation for every study; the legacy A2780 comparison moves into the
A2780 study's regression test.

### 5.10 Qualification scorecard

One table, one row per carried model, PASS / WARN / FAIL against configurable thresholds:

| Check | Default pass rule |
|---|---|
| Fit converged within iteration caps | no cap hit on the primary fit |
| Selection stable | carried model selected in ≥ 50% of bootstrap replicates |
| Parameters identifiable | 95% CI upper/lower < 100; parameters on bounds/floors flagged |
| Residuals unbiased | no runs-test failure at α = 0.01 |
| Forecast (last 30% of days) | RMSE ≤ 2× in-sample |
| Leave-one-block-out | per-block RMSE ≤ 2× in-sample |
| Held-out condition / schedule (if marked) | beats the no-treatment baseline |
| Bootstrap time-limit hits | ≤ 2% |

Every column and "i" tip states whether lower / higher / closer-to-target is better (existing `tip(text, better)`
rule in `info_tips.jl`; unknown columns stop the build).

### 5.11 Reproducibility and setup

- Pinned Julia version (checked by `setup.jl` with a clear message), committed `Manifest.toml`,
  `GrowthParameterEstimation.jl` pinned by git-tree SHA.
- `cgd setup` (first run): instantiate + precompile; optionally build a sysimage so later runs start fast.
- Seeds derived from (stage, line set, candidate, replicate), never thread id; results independent of `--threads`.
- Optimisers stop on iteration caps only; no wall-clock limit on primary fits (keep current behaviour).
- `cgd run --detached` starts long runs with a keep-awake request (reuse `sensitivity/run_awake.ps1` on Windows,
  `caffeinate`/`systemd-inhibit` elsewhere) — background tasks are killed after ~3 h and sleep has paused runs.
- A `--quick` mode (fewer bootstrap replicates and multistarts, clearly labelled "preview") so new users see
  results in minutes before committing to a multi-hour run.

## 6. Domain rules to preserve (non-negotiable; encode in code and tests)

1. **Model selection:** carry forward the model with the **fewest free parameters within 2 BIC of the lowest**
   (ties: lower BIC). No separate "simple vs complex" comparisons.
2. **Autonomous models:** no explicit time terms (so any schedule D(t) can drive them for control);
   **≤ 3 compartments per line**.
3. Treatment mechanisms are shared across stages through inheritance; washout/reset compartments not supported.
4. **Day 0 = nominal seeding**, the initial condition, not a data point; x0 estimated only when the first measured
   point is > 3-fold off nominal, and that parameter counts in the BIC.
5. **EC50 per line per drug**, user-supplied (§5.6).
6. Values already normalised to cell size are ≈ counts; never rescale them.
7. Defaults: residual weighting `power`, p = 0.55; bootstrap generated from the carried model.
8. Thesis stage names (1, 2a, 2b, 3a, 3b, 4) only in the A2780 publish step; generic reports use the generic stage
   names of §5.5.

## 7. Implementation plan (each milestone ends green)

Branch `study-pipeline` from `main`.

| # | Milestone | Done when |
|---|---|---|
| M0 | Baseline | A2780 regression passes on `main`; grep test (§5.5) in report mode lists every hard-coded string. |
| M1 | Results folder + study folder | `studies/a2780/` holds today's configs (explicit `stages.toml`); `cgd run studies/a2780` writes to `results/`; `publish` reproduces today's `docs/`. **Regression 31/31.** |
| M2 | Tagging: `init`, `tags.xlsx`, compiler to TOML, reader layouts (§5.1–5.3), validator (§5.9) | A2780 raw files → `init` → tags that compile to a `sources.toml` semantically equal to today's (test diffs them); acceptance test 5 passes. |
| M3 | De-hard-code names and constants (N = 2 path still) | Grep test enforcing; **regression 31/31**. |
| M4 | Auto stage plan + skipping; EC50 handling (§5.6) | TykNu runs end to end (acceptance 2). For information only, report how the auto plan for A2780 differs from its explicit `stages.toml`; A2780 always uses its own. |
| M5 | N lines / N treatments: monoculture + combination | Acceptance test 3. **Regression 31/31.** |
| M6 | N-line coculture | Acceptance test 4. **Regression 31/31.** |
| M7 | Generalised control (§5.7) + scorecard (§5.10) + generic report | Control and scorecard render for A2780, TykNu and synthetic studies. |
| M8 | Setup, demo study, `--quick`, `--assist`, README quick start, determinism + clean-clone tests | Acceptance 6 and 7. |

- M1, M3, M5, M6 touch numerical code: the regression test is the safety net after every subagent change (fast
  regression during development; one full detached run per milestone).
- M2 (I/O + tagging) is independent of M3 and can run in parallel.

## 8. Risks

| Risk | Mitigation |
|---|---|
| Refactor silently changes A2780 numbers | Regression after every change; byte-identical determinism test; optimiser settings untouched. |
| N-line coculture / many-treatment candidate explosion, long runs | Candidates generated by increasing complexity with a reported cap; 2-BIC rule; `--quick` preview. |
| Wrong tags (dose or line) fit silently | Guessed cells highlighted with their evidence; red cells block; validator plausibility warnings (§5.9). |
| Unidentifiable parameters with sparse new data | Scorecard flags them; EC50 estimation only with ≥ 3 doses. |
| QC thresholds tuned to A2780 count scale | Scale `floor` by `reference_count`; report exclusion rate per series. |
| Long runs killed by sleep or tool limits | `--detached` + keep-awake; `--resume`. |

## 9. Resolved questions (owner, 2026-10-10)

1. **Combination candidates:** independent action, Bliss and the current χ interaction. Confirmed.
2. **Control objective:** literature-based presets as defaults, which the user can override (§5.7).
3. **The A2780 study keeps its own configuration permanently:** its explicit `stages.toml`,
   `validation.toml`, `control_config.toml` and publish step. The auto stage plan is never applied to it.
   M4's "auto plan matches" check is informational only.
