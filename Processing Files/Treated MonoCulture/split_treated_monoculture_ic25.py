#!/usr/bin/env python3
"""
Split Treated MonoCulture IC25 (1.47 uM) CSVs into 20k/30k and A2780Naive/A2780cis.

User-specified mapping for IC25:
  - tiles A4-A6 -> 20k Naive
  - tiles B4-B6 -> 30k Naive
  - tiles A1-A3 -> 20k Cis
  - tiles B1-B3 -> 30k Cis

Assumptions:
  - There are two CSVs with '1.47uMTreated' in their filenames:
      * one "naive" (no 'cis' in filename)
      * one "cis" (contains 'cis' in filename)
  - Prefer processed files in Processed_Datasets/Treated MonoCulture if present,
    else read raw Datasets/Treated MonoCulture.
  - Outputs are written to Processed_Datasets/Treated MonoCulture/{20k,30k}/IC25/.
"""
from __future__ import annotations

import csv
from pathlib import Path
import re
from typing import Dict, List, Tuple

WELL_RE = re.compile(r"_([AB])([1-6])\.tif$", re.IGNORECASE)


def find_project_root(start: Path) -> Path:
    current = start
    for _ in range(6):
        if (current / "Datasets").exists():
            return current
        if current.parent == current:
            break
        current = current.parent
    return start


def load_csv_rows(csv_path: Path) -> Tuple[List[str], List[List[str]]]:
    # Read with BOM tolerance, fallback to latin-1 if necessary
    try:
        f = csv_path.open("r", encoding="utf-8-sig", newline="")
    except Exception:
        f = csv_path.open("r", encoding="latin-1", newline="")
    with f:
        reader = csv.reader(f)
        header = next(reader, None)
        if header is None:
            return [], []
        rows = [r for r in reader if r and len(r) == len(header)]
        return header, rows


def classify_ic25(image_name: str, is_naive_file: bool) -> Tuple[str, str] | None:
    """Return (seeding, line) for IC25 well based on custom mapping.

    Mapping for IC25:
      - A4-A6 -> 20k Naive         (use only when is_naive_file=True)
      - B4-B6 -> 30k Naive         (use only when is_naive_file=True)
      - A1-A3 -> 20k Cis           (use only when is_naive_file=False)
      - B1-B3 -> 30k Cis           (use only when is_naive_file=False)
    """
    m = WELL_RE.search(image_name)
    if not m:
        return None
    row = m.group(1).upper()
    col = int(m.group(2))

    if is_naive_file:
        if row == "A" and col in (4, 5, 6):
            return "20k", "A2780Naive"
        if row == "B" and col in (4, 5, 6):
            return "30k", "A2780Naive"
        return None
    else:
        if row == "A" and col in (1, 2, 3):
            return "20k", "A2780cis"
        if row == "B" and col in (1, 2, 3):
            return "30k", "A2780cis"
        return None


def split_ic25(dataset_dir: Path, out_root: Path) -> None:
    # Find the two 1.47uM files
    candidates = sorted(dataset_dir.glob("*1.47uMTreated*.csv"))
    if not candidates:
        print(f"No IC25 (1.47uM) CSVs found in {dataset_dir}")
        return

    naive_csv: Path | None = None
    cis_csv: Path | None = None
    for p in candidates:
        name_lower = p.name.lower()
        if "cis" in name_lower:
            cis_csv = p
        else:
            naive_csv = p

    if naive_csv is None or cis_csv is None:
        print("Warning: Expected two files (naive and cis) with '1.47uMTreated' in the name.")
        print("Found:")
        for p in candidates:
            print(f" - {p}")
        # Continue with any that are present

    # Prepare output dirs
    (out_root / "20k" / "IC25").mkdir(parents=True, exist_ok=True)
    (out_root / "30k" / "IC25").mkdir(parents=True, exist_ok=True)

    # Buckets: (seeding, line) -> rows
    buckets: Dict[Tuple[str, str], List[List[str]]] = {
        ("20k", "A2780Naive"): [],
        ("30k", "A2780Naive"): [],
        ("20k", "A2780cis"): [],
        ("30k", "A2780cis"): [],
    }

    header: List[str] | None = None

    def process_one(csv_path: Path, is_naive_file: bool) -> Tuple[int, int]:
        nonlocal header
        if csv_path is None:
            return (0, 0)
        h, rows = load_csv_rows(csv_path)
        if not h:
            return (0, 0)
        if header is None:
            header = h
        kept = 0
        for r in rows:
            cls = classify_ic25(r[0], is_naive_file)
            if cls is None:
                continue
            buckets[cls].append(r)
            kept += 1
        return (len(rows), kept)

    total_naive, kept_naive = process_one(naive_csv, True)
    total_cis, kept_cis = process_one(cis_csv, False)

    if header is None:
        print("No IC25 data extracted.")
        return

    # Write outputs
    outputs: List[Tuple[Path, int]] = []
    for seeding in ("20k", "30k"):
        for line in ("A2780Naive", "A2780cis"):
            out_dir = out_root / seeding / "IC25"
            out_file = out_dir / f"{line}.csv"
            rows = buckets[(seeding, line)]
            rows.sort(key=lambda r: r[0])
            with out_file.open("w", encoding="utf-8-sig", newline="") as f:
                w = csv.writer(f)
                w.writerow(header)
                w.writerows(rows)
            outputs.append((out_file, len(rows)))

    print("IC25 split complete:\n-----------------")
    for path, count in outputs:
        print(f"{count:5d} rows -> {path}")
    print("-----------------")
    print(f"From naive file: kept {kept_naive}/{total_naive} rows")
    print(f"From cis file  : kept {kept_cis}/{total_cis} rows")


def main() -> None:
    here = Path(__file__).resolve().parent
    root = find_project_root(here)

    processed_dir = root / "Processed_Datasets" / "Treated MonoCulture"
    raw_dir = root / "Datasets" / "Treated MonoCulture"

    # Use raw if processed doesn't contain the 1.47 files
    dataset_dir = processed_dir if any(processed_dir.glob("*1.47uMTreated*.csv")) else raw_dir
    if not dataset_dir.exists():
        raise SystemExit(
            f"Treated MonoCulture directory not found in either:\n  {processed_dir}\n  {raw_dir}"
        )

    out_root = root / "Processed_Datasets" / "Treated MonoCulture"
    split_ic25(dataset_dir, out_root)


if __name__ == "__main__":
    main()
