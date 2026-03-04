"""Split treated co-culture 20k/30k files into cis and non_cis subfolders.

Classification:
- filename containing 'A2780cis' (case-insensitive) -> cis
- otherwise -> non_cis

Copies files (does not move) to preserve existing files and avoid information loss.

Input:
  Processed_Datasets/Treated Coculture/{20k,30k}/measure_*.csv

Outputs:
  Processed_Datasets/Treated Coculture/{20k,30k}/cis/<same_filename>
  Processed_Datasets/Treated Coculture/{20k,30k}/non_cis/<same_filename>
  Processed_Datasets/Treated Coculture/cis_split_summary.csv
"""
from __future__ import annotations

from pathlib import Path
import shutil
import pandas as pd

DENSITIES = ["20k", "30k"]


def find_repo_root(start: Path) -> Path:
    cur = start
    for _ in range(12):
        if (cur / "Processed_Datasets").exists() or (cur / "Datasets").exists():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return start


def main() -> None:
    repo_root = find_repo_root(Path(__file__).resolve().parent)
    base = repo_root / "Processed_Datasets" / "Treated Coculture"

    records = []
    for density in DENSITIES:
        folder = base / density
        if not folder.exists():
            continue

        cis_dir = folder / "cis"
        non_dir = folder / "non_cis"
        cis_dir.mkdir(parents=True, exist_ok=True)
        non_dir.mkdir(parents=True, exist_ok=True)

        for csv in sorted(folder.glob("measure_*.csv")):
            if not csv.is_file():
                continue
            is_cis = "a2780cis" in csv.name.lower()
            target_dir = cis_dir if is_cis else non_dir
            dst = target_dir / csv.name
            shutil.copy2(csv, dst)
            records.append(
                {
                    "density": density,
                    "file": csv.name,
                    "cis": is_cis,
                    "destination": str(dst),
                    "action": "copied",
                }
            )

    summary = pd.DataFrame(records)
    out = base / "cis_split_summary.csv"
    summary.to_csv(out, index=False)
    print(f"Summary written: {out}")
    if not summary.empty:
        print(summary[["density", "file", "cis", "action"]])


if __name__ == "__main__":
    main()
