#!/usr/bin/env python3
"""
Compute per-day, per-well sample averages (6-tile means) for Treated MonoCulture.

For each leaf CSV (e.g., 20k/IC50/A2780Naive.csv), this script:
1) Parses Day and Well from the Image column (e.g., ..._Day5_..._A2.tif)
2) Aggregates tiles to a single mean per (Day, Well)

Writes outputs under a sibling Averages folder:
  <dir>/Averages/<basename>_sample_averages.csv

Columns: Day, Well, N Tiles, Mean Cells, SD Cells, SEM Cells
"""
from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path
from statistics import mean, stdev
from typing import Dict, List, Tuple

DAY_RE = re.compile(r"_Day(\d+)_", re.IGNORECASE)
WELL_RE = re.compile(r"_([ABCD][1-6])\.tif$", re.IGNORECASE)


def parse_day_and_well(image: str) -> Tuple[int, str] | None:
    d = DAY_RE.search(image)
    w = WELL_RE.search(image)
    if not d or not w:
        return None
    return int(d.group(1)), w.group(1).upper()


def per_day_sample_stats(input_csv: Path, output_csv: Path) -> int:
    tiles: Dict[Tuple[int, str], List[float]] = {}

    try:
        fh = input_csv.open("r", encoding="utf-8-sig", newline="")
    except Exception:
        fh = input_csv.open("r", encoding="latin-1", newline="")
    with fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header:
            return 0
        for row in reader:
            if not row or len(row) < 2:
                continue
            image, cells_str = row[0], row[1]
            parsed = parse_day_and_well(image)
            if not parsed:
                continue
            try:
                val = float(cells_str)
            except ValueError:
                continue
            tiles.setdefault(parsed, []).append(val)

    def well_sort_key(w: str):
        return (w[0], int(w[1]))

    rows: List[Tuple[int, str, int, float, float, float]] = []
    for (day, well), vals in tiles.items():
        if not vals:
            continue
        ntiles = len(vals)
        m = mean(vals)
        sd = stdev(vals) if ntiles >= 2 else 0.0
        sem = sd / math.sqrt(ntiles) if ntiles > 0 else 0.0
        rows.append((day, well, ntiles, m, sd, sem))

    rows.sort(key=lambda r: (r[0], well_sort_key(r[1])))

    with output_csv.open("w", encoding="utf-8-sig", newline="") as fo:
        writer = csv.writer(fo)
        writer.writerow(["Day", "Well", "N Tiles", "Mean Cells", "SD Cells", "SEM Cells"])
        for day, well, ntiles, m, sd, sem in rows:
            writer.writerow([day, well, ntiles, f"{m:.6f}", f"{sd:.6f}", f"{sem:.6f}"])

    return len(rows)


def find_treated_inputs(root: Path) -> List[Path]:
    base = root / "Processed_Datasets" / "Treated MonoCulture"
    inputs: List[Path] = []
    if not base.exists():
        return inputs
    for seeding in ("20k", "30k"):
        sdir = base / seeding
        if not sdir.exists():
            continue
        for ic in ("IC25", "IC50", "IC75"):
            icdir = sdir / ic
            if not icdir.exists():
                continue
            for csv_path in icdir.glob("*.csv"):
                name = csv_path.name.lower()
                if name.endswith("_sample_averages.csv") or name.endswith("_day_averages.csv") or name == "day_averages.csv":
                    continue
                inputs.append(csv_path)
    return inputs


def main(argv: List[str]) -> int:
    here = Path(__file__).resolve().parent
    root = here
    for _ in range(6):
        if (root / "Processed_Datasets").exists():
            break
        if root.parent == root:
            break
        root = root.parent

    inputs: List[Path] = [Path(a) for a in argv] if argv else find_treated_inputs(root)
    if not inputs:
        print("No Treated MonoCulture input CSVs found to process.")
        return 1

    print(f"Found {len(inputs)} Treated leaf CSVs to process:")
    for i, p in enumerate(inputs, 1):
        print(f"  {i}. {p}")

    total_rows = 0
    for inp in inputs:
        try:
            out_dir = inp.parent / "Averages"
            out_dir.mkdir(parents=True, exist_ok=True)
            out = out_dir / (inp.stem + "_sample_averages.csv")
            print(f"Processing: {inp} -> {out}")
            n = per_day_sample_stats(inp, out)
            total_rows += n
            print(f"Wrote {n:4d} (Day,Well) rows -> {out}")
        except Exception as e:
            print(f"ERROR processing {inp}: {e}")

    print(f"Done. Total (Day,Well) rows written: {total_rows}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
