#!/usr/bin/env python3
"""
Compute per-file day averages with error bars from per-well sample averages
for Treated MonoCulture (supports IC25/IC50/IC75 under 20k/30k).

For each Averages folder (e.g., .../20k/IC50/Averages/), this script:
- Finds each *_sample_averages.csv (per input dataset)
- Aggregates per day across wells (generally three wells per day) to produce Mean/SD/SEM
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
                m_str = (row.get("Mean Cells") or row.get("Mean Area µm^2") or row.get("Mean") or "").strip()
                m = float(m_str)
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
            writer.writerow(["Day", "N Samples", "Mean Cells", "SD Cells", "SEM Cells"])
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
    for _ in range(6):
        if (root / "Processed_Datasets").exists():
            break
        if root.parent == root:
            break
        root = root.parent

    base = root / "Processed_Datasets" / "Treated MonoCulture"
    wrote = 0
    for seeding in ("20k", "30k"):
        for ic in ("IC25", "IC50", "IC75"):
            avg_dir = base / seeding / ic / "Averages"
            if not avg_dir.exists():
                continue
            wrote += write_per_file_day_averages(avg_dir)

    if wrote == 0:
        print("No averages written (no Averages folders or no sample averages files found).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
