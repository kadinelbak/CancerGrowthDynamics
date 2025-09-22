#!/usr/bin/env python3
"""
Batch-fix headers in Processed_Datasets to standardize on Cells.

Actions:
- If a CSV contains 'Area µm^2' (or legacy Mean/SD/SEM columns with 'Area µm^2'),
  rename to 'Cells' equivalents without changing numeric values (we assume
  inputs in Processed_Datasets are already converted to Cells by process_area_data.py).
- Known patterns handled:
  * Raw-like files: 'Area µm^2' -> 'Cells'
  * Sample averages: 'Mean Area µm^2' -> 'Mean Cells', 'SD Area µm^2' -> 'SD Cells', 'SEM Area µm^2' -> 'SEM Cells'
  * Day averages: 'Mean Area µm^2' -> 'Mean Cells', 'SD Area µm^2' -> 'SD Cells', 'SEM Area µm^2' -> 'SEM Cells'
  * Intermittent variants: 'Average_Area_um2' -> 'Average_Cells', etc.

This script edits files in place under Processed_Datasets/.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

RENAMES = {
    # raw
    'Area µm^2': 'Cells',
    # sample/day averages (untreated monoculture)
    'Mean Area µm^2': 'Mean Cells',
    'SD Area µm^2': 'SD Cells',
    'SEM Area µm^2': 'SEM Cells',
    # intermittent scripts
    'Average_Area_um2': 'Average_Cells',
    'StdDev_Area_um2': 'StdDev_Cells',
    'Min_Area_um2': 'Min_Cells',
    'Max_Area_um2': 'Max_Cells',
    'Mean_Area_um2': 'Mean_Cells',
    'Std_Dev_um2': 'Std_Dev_Cells',
    'Std_Error_um2': 'Std_Error_Cells',
    'CI_95_Margin_um2': 'CI_95_Margin_Cells',
    'CI_95_Lower_um2': 'CI_95_Lower_Cells',
    'CI_95_Upper_um2': 'CI_95_Upper_Cells',
}


def rewrite_csv_headers(path: Path) -> bool:
    try:
        with path.open('r', encoding='utf-8-sig', newline='') as f:
            reader = csv.reader(f)
            rows = list(reader)
        if not rows:
            return False
        header = rows[0]
        new_header = [RENAMES.get(h, h) for h in header]
        if new_header == header:
            return False
        with path.open('w', encoding='utf-8-sig', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(new_header)
            writer.writerows(rows[1:])
        print(f"Fixed headers in: {path}")
        return True
    except Exception as e:
        print(f"Failed to fix {path}: {e}")
        return False


def main(argv):
    root = Path(__file__).resolve().parent
    # find Processed_Datasets up the tree
    base = root
    for _ in range(6):
        if (base / 'Processed_Datasets').exists():
            break
        if base.parent == base:
            break
        base = base.parent
    target = base / 'Processed_Datasets'
    if not target.exists():
        print(f"Processed_Datasets not found at {target}")
        return 1
    paths = list(target.rglob('*.csv'))
    print(f"Scanning {len(paths)} CSVs under {target}...")
    changed = 0
    for p in paths:
        if rewrite_csv_headers(p):
            changed += 1
    print(f"Done. Updated headers in {changed} files.")
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
