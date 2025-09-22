# Treated MonoCulture Splitter

This folder contains `split_treated_monoculture.py`, which organizes Treated MonoCulture CSVs into nested directories by:

- Seeding density: 20k vs 30k (A/C → 20k, B/D → 30k)
- IC level: IC50 vs IC75 (A/B → IC50, C/D → IC75)
- Line: A2780Naive (columns 1–3) vs A2780cis (columns 4–6)

IC25 data are intentionally excluded (to be handled later).

## Inputs
- Prefers processed CSVs in `Processed_Datasets/Treated MonoCulture/` if present; otherwise falls back to raw `Datasets/Treated MonoCulture/`.
- Parses well from image names like `..._A1.tif` and tolerates occasional typos in IC tokens (e.g., `IC75reated`).

## Outputs
Files are written under:
```
Processed_Datasets/Treated MonoCulture/
  20k/
    IC50/ A2780Naive.csv, A2780cis.csv
    IC75/ A2780Naive.csv, A2780cis.csv
  30k/
    IC50/ A2780Naive.csv, A2780cis.csv
    IC75/ A2780Naive.csv, A2780cis.csv
    IC25/ A2780Naive.csv, A2780cis.csv  (when running the IC25 script)
```

## Run
From the repo root, activate your environment and run:

```powershell
# Windows PowerShell
. .\.venv\Scripts\Activate.ps1
python "Processing Files/Treated MonoCulture/split_treated_monoculture.py"
```

If using the VS Code task added by the assistant, run the task named:
- "Split Treated MonoCulture (venv)"
- "Split Treated MonoCulture IC25 (venv)" (for the IC25 mapping below)

## Notes
- IC25 rows (detected via `IC25treated` token) are excluded.
- If both processed and raw datasets exist, the processed ones are used.
- Headers are preserved; rows are sorted by image name for determinism.

---

## IC25 (1.47 uM) special splitter

Use `split_treated_monoculture_ic25.py` to extract and split the two 1.47 uM CSVs into `IC25` folders.

Mapping provided for IC25:
- A4–A6 → 20k Naive
- B4–B6 → 30k Naive
- A1–A3 → 20k Cis
- B1–B3 → 30k Cis

Run with PowerShell:

```powershell
. .\.venv\Scripts\Activate.ps1
python "Processing Files/Treated MonoCulture/split_treated_monoculture_ic25.py"
```
