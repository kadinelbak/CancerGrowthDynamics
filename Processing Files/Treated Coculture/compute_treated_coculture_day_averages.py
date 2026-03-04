"""Compute per-day averages across wells for treated co-culture datasets.

Uses per-file inputs:
  Processed_Datasets/Treated Coculture/<density>/<cis_status>/Averages/*_well_day_averages.csv

For each per-file table:
  - day_mean_value = mean(mean_value over wells) per Day
  - n_wells = count unique Well per Day

Outputs per-file:
  <basename>_day_averages.csv

Outputs per-folder:
  combined_day_averages.csv

Global summary:
  Processed_Datasets/Treated Coculture/day_averages_generation_summary.csv
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd

DENSITIES = ["20k", "30k"]
CIS_STATES = ["cis", "non_cis"]


def find_repo_root(start: Path) -> Path:
    cur = start
    for _ in range(12):
        if (cur / "Processed_Datasets").exists() or (cur / "Datasets").exists():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return start


def process_file(path: Path):
    df = pd.read_csv(path)
    required = {"Day", "Well", "mean_value"}
    if not required.issubset(df.columns):
        return None

    df["Day"] = pd.to_numeric(df["Day"], errors="coerce")
    grp = (
        df.dropna(subset=["Day"])
        .groupby("Day", as_index=False)
        .agg(day_mean_value=("mean_value", "mean"), n_wells=("Well", "nunique"))
        .sort_values("Day")
    )
    if grp.empty:
        return None

    grp["day_mean_value"] = grp["day_mean_value"].round(2)
    out_path = path.parent / (path.stem.replace("_well_day_averages", "") + "_day_averages.csv")
    grp.to_csv(out_path, index=False)
    return out_path, grp


def main() -> None:
    repo_root = find_repo_root(Path(__file__).resolve().parent)
    base = repo_root / "Processed_Datasets" / "Treated Coculture"

    rows = []
    for density in DENSITIES:
        for status in CIS_STATES:
            avg_dir = base / density / status / "Averages"
            if not avg_dir.exists():
                continue

            combined_parts = []
            for input_file in sorted(avg_dir.glob("*_well_day_averages.csv")):
                result = process_file(input_file)
                if not result:
                    continue
                out_path, grp = result
                grp2 = grp.copy()
                grp2.insert(0, "source_file", input_file.name)
                combined_parts.append(grp2)

                rows.append(
                    {
                        "density": density,
                        "cis_status": status,
                        "source_file": input_file.name,
                        "output_file": out_path.name,
                        "days": int(grp["Day"].nunique()),
                        "mean_first_day": float(grp.loc[grp["Day"].idxmin(), "day_mean_value"]),
                    }
                )

            if combined_parts:
                combined = pd.concat(combined_parts, ignore_index=True)
                combined.to_csv(avg_dir / "combined_day_averages.csv", index=False)

    summary = pd.DataFrame(rows)
    out_summary = base / "day_averages_generation_summary.csv"
    summary.to_csv(out_summary, index=False)
    print(f"Summary written: {out_summary}")
    if not summary.empty:
        print(summary.head())


if __name__ == "__main__":
    main()
