"""Compute per-well per-day tile averages for treated co-culture datasets.

Reads from:
  Processed_Datasets/Treated Coculture/<density>/<cis_status>/measure_*.csv

Parses from Image name:
  Day<d>_<CHANNEL>_Tile-<t>_<Well>.tif

Produces per-file averages in:
  Processed_Datasets/Treated Coculture/<density>/<cis_status>/Averages/
      <basename>_well_day_averages.csv

Also writes combined files per density/cis status:
  combined_well_day_averages.csv

Value column handling:
- Prefers 'Cells µm^2'
- falls back to 'Cells'
- falls back to 'cell_area_um2'

Output metric name is standardized as 'mean_value'.
"""
from __future__ import annotations

import re
from pathlib import Path
import pandas as pd

DENSITIES = ["20k", "30k"]
CIS_STATES = ["cis", "non_cis"]
IMG_RE = re.compile(r"Day(\d+)_([A-Za-z ]+)_Tile-(\d+)_([A-D][1-6])\.tif$", re.IGNORECASE)
VALUE_COL_CANDIDATES = ["Cells µm^2", "Cells", "cell_area_um2"]


def find_repo_root(start: Path) -> Path:
    cur = start
    for _ in range(12):
        if (cur / "Processed_Datasets").exists() or (cur / "Datasets").exists():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return start


def parse_image(image_name: str):
    m = IMG_RE.search(str(image_name))
    if not m:
        return None, None, None, None
    return int(m.group(1)), m.group(2).strip(), int(m.group(3)), m.group(4).upper()


def pick_value_column(columns: list[str]) -> str | None:
    for candidate in VALUE_COL_CANDIDATES:
        if candidate in columns:
            return candidate
    return None


def well_sort_key(well: str) -> tuple[int, int]:
    if not isinstance(well, str) or len(well) < 2:
        return (99, 99)
    row_order = {"A": 0, "B": 1, "C": 2, "D": 3}
    row = row_order.get(well[0].upper(), 99)
    try:
        col = int(well[1:])
    except ValueError:
        col = 99
    return (row, col)


def process_file(csv_path: Path, out_dir: Path) -> dict | None:
    df = pd.read_csv(csv_path)
    if "Image" not in df.columns:
        return None

    value_col = pick_value_column(list(df.columns))
    if value_col is None:
        return None

    parsed = df["Image"].map(parse_image)
    parsed_df = pd.DataFrame(parsed.tolist(), columns=["Day", "Channel", "Tile", "Well"], index=df.index)
    df = pd.concat([df, parsed_df], axis=1)

    df = df.dropna(subset=["Day", "Well", value_col]).copy()
    if df.empty:
        return None

    grp = (
        df.groupby(["Day", "Well"], as_index=False)
        .agg(
            mean_value=(value_col, "mean"),
            n_tiles=("Tile", "nunique"),
            n_channels=("Channel", "nunique"),
        )
    )
    grp["mean_value"] = grp["mean_value"].round(2)
    grp = grp.sort_values(by=["Day", "Well"], key=lambda c: [well_sort_key(v) if c.name == "Well" else v for v in c])

    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{csv_path.stem}_well_day_averages.csv"
    grp.to_csv(out_file, index=False)

    return {
        "file": csv_path.name,
        "value_column": value_col,
        "rows_in": len(df),
        "days": grp["Day"].nunique(),
        "wells": grp["Well"].nunique(),
        "output": str(out_file),
    }


def main() -> None:
    repo_root = find_repo_root(Path(__file__).resolve().parent)
    base = repo_root / "Processed_Datasets" / "Treated Coculture"

    summary_records = []
    for density in DENSITIES:
        for status in CIS_STATES:
            folder = base / density / status
            if not folder.exists():
                continue

            avg_dir = folder / "Averages"
            combined_parts = []

            for csv_file in sorted(folder.glob("measure_*.csv")):
                rec = process_file(csv_file, avg_dir)
                if rec is None:
                    continue
                summary_records.append({**rec, "density": density, "cis_status": status})
                part = pd.read_csv(rec["output"])
                part.insert(0, "source_file", csv_file.name)
                combined_parts.append(part)

            if combined_parts:
                combined = pd.concat(combined_parts, ignore_index=True)
                combined.to_csv(avg_dir / "combined_well_day_averages.csv", index=False)

    summary_df = pd.DataFrame(summary_records)
    out = base / "averages_generation_summary.csv"
    summary_df.to_csv(out, index=False)
    print(f"Summary written: {out}")
    if not summary_df.empty:
        print(summary_df[["density", "cis_status", "file", "days", "wells", "value_column"]])


if __name__ == "__main__":
    main()
