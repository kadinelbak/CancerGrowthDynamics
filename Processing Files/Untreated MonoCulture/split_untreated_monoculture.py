#!/usr/bin/env python3
"""
Split Untreated MonoCulture datasets into 20k and 30k CSVs by well label (A1-6, B1-6).

Rules:
- A1-A3  -> 20k/A2780Naive.csv
- A4-A6  -> 20k/A2780cis.csv
- B1-B3  -> 30k/A2780Naive.csv
- B4-B6  -> 30k/A2780cis.csv

Original CSVs are not modified; new CSVs are written under:
Processed_Datasets/Untreated MonoCulture/{20k,30k}/

This script is self-rooting: it searches upwards to find the repository root that contains
'Datasets' and 'Processed_Datasets'. You can place and run it from anywhere in this repo.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, List, Tuple

# Patterns like ..._A1.tif or ..._B6.tif at the end of the filename
WELL_RE = re.compile(r"_([AB])([1-6])\.tif$", re.IGNORECASE)

# Output mapping keys
Seeding = str  # '20k' | '30k'
LineName = str  # 'A2780Naive' | 'A2780cis'


def find_project_root(start: Path) -> Path:
    """Walk upwards from start to find a directory that has 'Datasets'."""
    current = start
    for _ in range(5):  # search up to 5 levels just in case
        if (current / "Datasets").exists():
            return current
        if current.parent == current:
            break
        current = current.parent
    return start  # fallback


def classify_image(image_name: str) -> Tuple[Seeding, LineName] | None:
    """Return (seeding, line) for an image filename or None if no well pattern found."""
    m = WELL_RE.search(image_name)
    if not m:
        return None
    row = m.group(1).upper()
    col = int(m.group(2))

    seeding: Seeding = "20k" if row == "A" else "30k"
    line: LineName = "A2780Naive" if col in (1, 2, 3) else "A2780cis"
    return seeding, line


def split_csvs(dataset_dir: Path, out_root: Path) -> None:
    # Ensure output directories exist
    (out_root / "20k").mkdir(parents=True, exist_ok=True)
    (out_root / "30k").mkdir(parents=True, exist_ok=True)

    # Prepare accumulators
    buckets: Dict[Tuple[Seeding, LineName], List[List[str]]] = {
        ("20k", "A2780Naive"): [],
        ("20k", "A2780cis"): [],
        ("30k", "A2780Naive"): [],
        ("30k", "A2780cis"): [],
    }

    header: List[str] | None = None

    # Process all CSVs in the Untreated MonoCulture dataset directory
    for csv_path in sorted(dataset_dir.glob("*.csv")):
        # Read with utf-8-sig to tolerate BOM; fallback to latin-1 if needed
        try:
            content_iter = csv_path.open("r", encoding="utf-8-sig", newline="")
        except Exception:
            content_iter = csv_path.open("r", encoding="latin-1", newline="")
        with content_iter as f:
            reader = csv.reader(f)
            try:
                local_header = next(reader)
            except StopIteration:
                continue  # empty file

            # First non-empty header encountered becomes the canonical header
            if header is None:
                header = local_header

            # Validate column count matches header; if not, skip row defensively
            expected_cols = len(header)
            for row in reader:
                if not row or len(row) != expected_cols:
                    continue
                image_name = row[0]
                cls = classify_image(image_name)
                if cls is None:
                    continue
                buckets[cls].append(row)

    if header is None:
        print("No CSV data found; nothing to write.")
        return

    # Write out the four CSVs
    outputs = [
        ("20k", "A2780Naive"),
        ("20k", "A2780cis"),
        ("30k", "A2780Naive"),
        ("30k", "A2780cis"),
    ]

    for seeding, line in outputs:
        out_dir = out_root / seeding
        out_file = out_dir / f"{line}.csv"
        rows = buckets[(seeding, line)]
        # Sort rows for determinism (by image name)
        rows.sort(key=lambda r: r[0])
        with out_file.open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)
        print(f"Wrote {len(rows):5d} rows -> {out_file}")


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    project_root = find_project_root(here)

    # Prefer to read input CSVs from Processed_Datasets if present; fallback to raw Datasets.
    processed_dir = project_root / "Processed_Datasets" / "Untreated MonoCulture"
    raw_dir = project_root / "Datasets" / "Untreated MonoCulture"

    if processed_dir.exists():
        dataset_dir = processed_dir
    elif raw_dir.exists():
        dataset_dir = raw_dir
    else:
        raise SystemExit(
            f"Neither processed nor raw dataset directories found:\n  {processed_dir}\n  {raw_dir}"
        )

    out_root = project_root / "Processed_Datasets" / "Untreated MonoCulture"

    split_csvs(dataset_dir, out_root)
