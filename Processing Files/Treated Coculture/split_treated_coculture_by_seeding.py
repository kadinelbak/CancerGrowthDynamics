"""Split treated co-culture CSVs into 20k and 30k by tile index while preserving originals.

Rule from experiment layout:
- 20k uses tiles 1-3
- 30k uses tiles 4-6

Many files are exported as zero-indexed tiles (Tile-0..Tile-5). This script detects
indexing style per file:
- If tiles include 0, mapping is 20k: 0-2 and 30k: 3-5
- Otherwise mapping is 20k: 1-3 and 30k: 4-6

Input files:
  Processed_Datasets/Treated Coculture/measure_*.csv

Outputs:
  Processed_Datasets/Treated Coculture/20k/<same_filename>
  Processed_Datasets/Treated Coculture/30k/<same_filename>
  Processed_Datasets/Treated Coculture/split_treated_coculture_report.csv
"""
from __future__ import annotations

import re
from pathlib import Path
import pandas as pd

IMG_TILE_RE = re.compile(r"_Tile-(\d+)_", re.IGNORECASE)


def find_repo_root(start: Path) -> Path:
    cur = start
    for _ in range(12):
        if (cur / "Processed_Datasets").exists() or (cur / "Datasets").exists():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return start


def extract_tile(image_name: str) -> int | None:
    m = IMG_TILE_RE.search(str(image_name))
    return int(m.group(1)) if m else None


def tile_groups_from_values(tile_values: pd.Series) -> tuple[set[int], set[int], str]:
    valid = sorted({int(v) for v in tile_values.dropna().astype(int).tolist()})
    if not valid:
        return set(), set(), "no-tiles"
    if 0 in valid:
        return {0, 1, 2}, {3, 4, 5}, "zero-indexed"
    return {1, 2, 3}, {4, 5, 6}, "one-indexed"


def process_file(src_csv: Path, out_20k: Path, out_30k: Path) -> dict:
    df = pd.read_csv(src_csv)
    if "Image" not in df.columns:
        return {
            "file": src_csv.name,
            "status": "missing-image-column",
            "rows_total": len(df),
            "rows_20k": 0,
            "rows_30k": 0,
        }

    tiles = df["Image"].map(extract_tile)
    g20, g30, indexing = tile_groups_from_values(tiles)
    if not g20 and not g30:
        return {
            "file": src_csv.name,
            "status": "no-tiles-parsed",
            "rows_total": len(df),
            "rows_20k": 0,
            "rows_30k": 0,
        }

    df = df.assign(__Tile=tiles)
    df20 = df[df["__Tile"].isin(g20)].drop(columns=["__Tile"])
    df30 = df[df["__Tile"].isin(g30)].drop(columns=["__Tile"])

    out_20k.mkdir(parents=True, exist_ok=True)
    out_30k.mkdir(parents=True, exist_ok=True)

    wrote20 = False
    wrote30 = False
    if not df20.empty:
        df20.to_csv(out_20k / src_csv.name, index=False)
        wrote20 = True
    if not df30.empty:
        df30.to_csv(out_30k / src_csv.name, index=False)
        wrote30 = True

    return {
        "file": src_csv.name,
        "status": "ok",
        "indexing": indexing,
        "rows_total": len(df),
        "rows_20k": len(df20),
        "rows_30k": len(df30),
        "wrote_20k": wrote20,
        "wrote_30k": wrote30,
    }


def main() -> None:
    repo_root = find_repo_root(Path(__file__).resolve().parent)
    base = repo_root / "Processed_Datasets" / "Treated Coculture"
    src_files = sorted(p for p in base.glob("measure_*.csv") if p.is_file())

    out_20k = base / "20k"
    out_30k = base / "30k"

    records = [process_file(p, out_20k, out_30k) for p in src_files]
    report = pd.DataFrame(records)
    report_path = base / "split_treated_coculture_report.csv"
    report.to_csv(report_path, index=False)

    print(f"Processed files: {len(src_files)}")
    print(f"Report written: {report_path}")
    if not report.empty:
        print(report[["file", "indexing", "rows_20k", "rows_30k", "status"]])


if __name__ == "__main__":
    main()
