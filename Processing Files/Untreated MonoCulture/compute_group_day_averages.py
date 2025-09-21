#!/usr/bin/env python3
"""
Compute per-file day averages with error bars from per-well sample averages.

For each Averages folder (e.g., .../20k/Averages/), this script:
- Finds each *_sample_averages.csv (per input dataset)
- Aggregates per day across wells (three per day) to produce Mean/SD/SEM
- Writes <sample>_day_averages.csv next to the input
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from statistics import mean, stdev
from typing import Dict, List


def read_per_well_means(sample_csv: Path) -> Dict[int, List[float]]:
    day_to_means: Dict[int, List[float]] = {}
    with sample_csv.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                day = int((row.get("Day") or "").strip())
                m = float((row.get("Mean Area µm^2") or row.get("Mean") or "").strip())
            except Exception:
                continue
            day_to_means.setdefault(day, []).append(m)
    return day_to_means


def write_per_file_day_averages(avg_dir: Path) -> int:
    sample_files = sorted(avg_dir.glob("*_sample_averages.csv"))
    if not sample_files:
        return 0

    wrote = 0
    for s in sample_files:
        day_to_means = read_per_well_means(s)
        days = sorted(day_to_means.keys())
        out_path = s.with_name(s.stem.replace("_sample_averages", "_day_averages") + ".csv")
        with out_path.open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Day", "N Samples", "Mean Area µm^2", "SD Area µm^2", "SEM Area µm^2"])
            for day in days:
                vals = day_to_means[day]
                n = len(vals)
                if n == 0:
                    continue
                mval = mean(vals)
                sd = stdev(vals) if n >= 2 else 0.0
                sem = sd / math.sqrt(n) if n > 0 else 0.0
                writer.writerow([day, n, f"{mval:.6f}", f"{sd:.6f}", f"{sem:.6f}"])
        wrote += len(days)
        print(f"Wrote {len(days):4d} days -> {out_path}")
    return wrote


def main() -> int:
    here = Path(__file__).resolve().parent
    root = here
    # climb up to find Processed_Datasets
    for _ in range(5):
        if (root / "Processed_Datasets").exists():
            break
        if root.parent == root:
            break
        root = root.parent

    targets = [
        root / "Processed_Datasets" / "Untreated MonoCulture" / "20k" / "Averages",
        root / "Processed_Datasets" / "Untreated MonoCulture" / "30k" / "Averages",
    ]

    wrote = 0
    for d in targets:
        if not d.exists():
            continue
        count = write_per_file_day_averages(d)
        wrote += count
    if wrote == 0:
        print("No averages written (no Averages folders or no sample averages files found).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
