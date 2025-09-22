#!/usr/bin/env python3
"""
Split Treated MonoCulture datasets into nested folders by seeding density (20k/30k)
and IC level (IC50/IC75), then into A2780Naive vs A2780cis based on well columns.

Rules provided:
- Wells: A and C -> 20k; B and D -> 30k
- Columns 1-3 are A2780Naive; 4-6 are A2780cis
- Rows A/B belong to IC50; rows C/D belong to IC75
- Exclude IC25 for now (to be organized later)

Input detection:
- Prefer processed files at Processed_Datasets/Treated MonoCulture/*.csv
  (e.g., normalized with 'Cells' column), else fall back to raw Datasets/...
- IC level also parsed from the Image name token like "IC50treated" or the typo
  variant "IC75reated". However, we trust row->IC mapping as the source of truth.

Output:
Processed_Datasets/Treated MonoCulture/
  ├─ 20k/
  │   ├─ IC50/
  │   │    ├─ A2780Naive.csv
  │   │    └─ A2780cis.csv
  │   └─ IC75/
  │        ├─ A2780Naive.csv
  │        └─ A2780cis.csv
  └─ 30k/
      ├─ IC50/
      │    ├─ A2780Naive.csv
      │    └─ A2780cis.csv
      └─ IC75/
           ├─ A2780Naive.csv
           └─ A2780cis.csv

Original CSVs are not modified.

You can run this from anywhere inside the repo.
"""
from __future__ import annotations

import csv
import glob
import re
from pathlib import Path
from typing import Dict, List, Tuple

# Match well coordinate suffix like ..._A1.tif or ..._D6.tif
WELL_RE = re.compile(r"_([ABCD])([1-6])\.tif$", re.IGNORECASE)

# Match IC token, tolerate missing 't' in 'treated' (seen as 'IC75reated') and optional spaces
IC_RE = re.compile(r"_IC\s*(25|50|75)\s*t?reated[ _]", re.IGNORECASE)

Seeding = str  # '20k' | '30k'
ICLevel = str  # 'IC50' | 'IC75'
LineName = str  # 'A2780Naive' | 'A2780cis'


def find_project_root(start: Path) -> Path:
    """Walk upwards to the directory that contains 'Datasets'."""
    current = start
    for _ in range(6):
        if (current / "Datasets").exists():
            return current
        if current.parent == current:
            break
        current = current.parent
    return start


def well_to_seeding(row_letter: str) -> Seeding:
    row_letter = row_letter.upper()
    return "20k" if row_letter in ("A", "C") else "30k"


def row_to_ic(row_letter: str) -> ICLevel:
    row_letter = row_letter.upper()
    return "IC50" if row_letter in ("A", "B") else "IC75"


def col_to_line(col_num: int) -> LineName:
    return "A2780Naive" if col_num in (1, 2, 3) else "A2780cis"


def parse_ic_from_name(image_name: str) -> str | None:
    m = IC_RE.search(image_name)
    if not m:
        return None
    return f"IC{m.group(1)}"


def classify(image_name: str) -> Tuple[Seeding, ICLevel, LineName] | None:
    """Classify a single image name into (seeding, icLevel, line).

    Excludes IC25 based on either IC token in name or row-derived mapping when
    token is missing. If the well coordinate cannot be parsed, returns None.
    """
    m = WELL_RE.search(image_name)
    if not m:
        return None

    row_letter = m.group(1).upper()
    col_num = int(m.group(2))

    seeding = well_to_seeding(row_letter)
    ic_from_row = row_to_ic(row_letter)
    line = col_to_line(col_num)

    # If IC token is present in the filename, use it for exclusion of IC25.
    # Otherwise trust row mapping (A/B -> IC50, C/D -> IC75).
    ic_token = parse_ic_from_name(image_name)
    if ic_token == "IC25":
        return None  # Exclude IC25 rows entirely for now

    # Derive final IC level (prefer row mapping; token is advisory only).
    ic_level: ICLevel = ic_from_row
    return seeding, ic_level, line


def split_csvs(dataset_dir: Path, out_root: Path) -> None:
    # Ensure output base exists
    (out_root / "20k" / "IC50").mkdir(parents=True, exist_ok=True)
    (out_root / "20k" / "IC75").mkdir(parents=True, exist_ok=True)
    (out_root / "30k" / "IC50").mkdir(parents=True, exist_ok=True)
    (out_root / "30k" / "IC75").mkdir(parents=True, exist_ok=True)

    # Accumulate rows per bucket
    buckets: Dict[Tuple[Seeding, ICLevel, LineName], List[List[str]]] = {
        ("20k", "IC50", "A2780Naive"): [],
        ("20k", "IC50", "A2780cis"): [],
        ("20k", "IC75", "A2780Naive"): [],
        ("20k", "IC75", "A2780cis"): [],
        ("30k", "IC50", "A2780Naive"): [],
        ("30k", "IC50", "A2780cis"): [],
        ("30k", "IC75", "A2780Naive"): [],
        ("30k", "IC75", "A2780cis"): [],
    }

    header: List[str] | None = None

    csv_files = sorted(dataset_dir.glob("*.csv"))
    if not csv_files:
        print(f"No CSV files found in {dataset_dir}")
        return

    skipped_no_well = 0
    skipped_ic25 = 0
    total_rows = 0

    for csv_path in csv_files:
        # Read with BOM tolerance, fallback to latin-1 if necessary
        try:
            f = csv_path.open("r", encoding="utf-8-sig", newline="")
        except Exception:
            f = csv_path.open("r", encoding="latin-1", newline="")
        with f:
            reader = csv.reader(f)
            try:
                local_header = next(reader)
            except StopIteration:
                continue

            if header is None:
                header = local_header
            expected_cols = len(header)

            for row in reader:
                if not row or len(row) != expected_cols:
                    continue
                total_rows += 1
                image_name = row[0]
                cls = classify(image_name)
                if cls is None:
                    # Determine if IC25 or malformed well
                    ic_token = parse_ic_from_name(image_name) or ""
                    if ic_token == "IC25":
                        skipped_ic25 += 1
                    else:
                        skipped_no_well += 1
                    continue
                buckets[cls].append(row)

    if header is None:
        print("No data found; nothing to write.")
        return

    # Write outputs
    outputs = []
    for seeding in ("20k", "30k"):
        for ic in ("IC50", "IC75"):
            for line in ("A2780Naive", "A2780cis"):
                out_dir = out_root / seeding / ic
                out_file = out_dir / f"{line}.csv"
                rows = buckets[(seeding, ic, line)]
                rows.sort(key=lambda r: r[0])
                with out_file.open("w", encoding="utf-8-sig", newline="") as wf:
                    writer = csv.writer(wf)
                    writer.writerow(header)
                    writer.writerows(rows)
                outputs.append((out_file, len(rows)))

    # Summary
    print("Split Treated MonoCulture complete:\n-------------------------------")
    for path, count in outputs:
        print(f"{count:5d} rows -> {path}")
    print("-------------------------------")
    print(f"Total input rows: {total_rows}")
    print(f"Excluded IC25   : {skipped_ic25}")
    print(f"Skipped (no well): {skipped_no_well}")


def main() -> None:
    here = Path(__file__).resolve().parent
    project_root = find_project_root(here)

    processed_dir = project_root / "Processed_Datasets" / "Treated MonoCulture"
    raw_dir = project_root / "Datasets" / "Treated MonoCulture"

    # Prefer processed (normalized) if present, else use raw
    dataset_dir = processed_dir if processed_dir.exists() else raw_dir
    if not dataset_dir.exists():
        raise SystemExit(
            f"Treated MonoCulture directory not found in either:\n  {processed_dir}\n  {raw_dir}"
        )

    out_root = project_root / "Processed_Datasets" / "Treated MonoCulture"
    split_csvs(dataset_dir, out_root)


if __name__ == "__main__":
    main()
