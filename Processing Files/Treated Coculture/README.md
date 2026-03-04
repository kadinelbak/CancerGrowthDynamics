# Treated CoCulture Processing Pipeline

This folder contains a 4-step processing pipeline that **keeps original flat CSVs intact** and writes organized outputs under:

- `Processed_Datasets/Treated Coculture/20k/...`
- `Processed_Datasets/Treated Coculture/30k/...`

## Step order

1. `split_treated_coculture_by_seeding.py`
   - Splits by tile index into `20k` and `30k`
   - Supports both tile numbering styles:
     - zero-indexed `Tile-0..5` => 20k=`0..2`, 30k=`3..5`
     - one-indexed `Tile-1..6` => 20k=`1..3`, 30k=`4..6`
2. `split_treated_coculture_cis_status.py`
   - Copies density files into `cis` and `non_cis`
3. `compute_treated_coculture_well_day_averages.py`
   - Creates per-file `*_well_day_averages.csv`
4. `compute_treated_coculture_day_averages.py`
   - Creates per-file `*_day_averages.csv`

## Reports produced

- `Processed_Datasets/Treated Coculture/split_treated_coculture_report.csv`
- `Processed_Datasets/Treated Coculture/cis_split_summary.csv`
- `Processed_Datasets/Treated Coculture/averages_generation_summary.csv`
- `Processed_Datasets/Treated Coculture/day_averages_generation_summary.csv`
