#!/usr/bin/env python3
"""
Compute per-day, per-well sample averages and error bars (SD, SEM) from processed
Untreated MonoCulture CSVs. Each sample corresponds to a single well's tiles
on a given day (typically 6 tiles), producing three sample averages per day
for each input file.

For each input CSV (e.g., 20k/A2780Naive.csv), the script:
1) Parses Day and Well from the Image column (e.g., ..._Day5_..._A2.tif)
2) Aggregates the tiles to a single mean per (Day, Well)

Outputs a file in an Averages subfolder named <basename>_sample_averages.csv with columns:
Day, Well, N Tiles, Mean Area µm^2, SD Area µm^2, SEM Area µm^2

No original files are modified or removed.
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
WELL_RE = re.compile(r"_([AB][1-6])\.tif$", re.IGNORECASE)


def parse_day_and_well(image: str) -> Tuple[int, str] | None:
    d = DAY_RE.search(image)
    w = WELL_RE.search(image)
    if not d or not w:
        return None
    return int(d.group(1)), w.group(1).upper()


def per_day_sample_stats(input_csv: Path, output_csv: Path) -> int:
    """Compute per (Day, Well) statistics over tiles and write CSV.
    Returns number of (Day, Well) rows written.
    """
    # Accumulate tiles -> areas per (Day, Well)
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
            image, area_str = row[0], row[1]
            parsed = parse_day_and_well(image)
            if not parsed:
                continue
            try:
                area = float(area_str)
            except ValueError:
                continue
            tiles.setdefault(parsed, []).append(area)

    # Prepare per (Day, Well) stats
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
        writer.writerow(["Day", "Well", "N Tiles", "Mean Area µm^2", "SD Area µm^2", "SEM Area µm^2"])
        for day, well, ntiles, m, sd, sem in rows:
            writer.writerow([day, well, ntiles, f"{m:.6f}", f"{sd:.6f}", f"{sem:.6f}"])

    return len(rows)


def find_default_inputs(root: Path) -> List[Path]:
    base = root / "Processed_Datasets" / "Untreated MonoCulture"
    inputs: List[Path] = []
    for group in ("20k", "30k"):
        d = base / group
        if not d.exists():
            continue
        for p in d.glob("*.csv"):
            name = p.name.lower()
            if name.endswith("_sample_averages.csv") or name.endswith("_day_averages.csv") or name == "day_averages.csv":
                continue
            inputs.append(p)
    return inputs


def main(argv: List[str]) -> int:
    here = Path(__file__).resolve().parent
    root = here
    for _ in range(5):
        if (root / "Processed_Datasets").exists():
            break
        if root.parent == root:
            break
        root = root.parent

    inputs: List[Path] = [Path(a) for a in argv] if argv else find_default_inputs(root)
    if not inputs:
        print("No input CSVs found to process.")
        return 1

    print(f"Found {len(inputs)} base CSVs to process:")
    for i, p in enumerate(inputs, 1):
        print(f"  {i}. {p}")

    total_rows = 0
    for inp in inputs:
        try:
            # Place outputs in an Averages subfolder under the same directory (20k or 30k)
            out_dir = inp.parent / "Averages"
            out_dir.mkdir(parents=True, exist_ok=True)
            out = out_dir / (inp.stem + "_sample_averages.csv")
            print(f"Processing: {inp} -> {out}")
            n = per_day_sample_stats(inp, out)
            total_rows += n
            print(f"Wrote {n:4d} (Day,Well) rows -> {out}")
        except Exception as e:
            print(f"ERROR processing {inp}: {e}")

    print(f"Done. Total (Day,Well) rows written across files: {total_rows}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
