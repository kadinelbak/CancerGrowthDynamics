"""Organize 20k and 30k untreated co-culture CSVs into cis / non_cis.

Files containing 'A2780cis' (case-insensitive) are classified as cis, otherwise non_cis.
Existing files are copied (not removed) to preserve flat listing; adjust MOVE_FILES=True to move instead.
Produces a summary CSV in each density folder.
"""
from __future__ import annotations
from pathlib import Path
import shutil
import pandas as pd

ROOT = Path(r"c:/Users/MainFrameTower/Desktop/CancerGrowthDynamics")
BASE = ROOT / "Processed_Datasets" / "Untreated CoCulture"
DENSITIES = ["20k", "30k"]
MOVE_FILES = False  # set True to move instead of copy


def classify_and_place(density: str):
    folder = BASE / density
    cis_dir = folder / 'cis'
    non_dir = folder / 'non_cis'
    cis_dir.mkdir(exist_ok=True)
    non_dir.mkdir(exist_ok=True)

    records = []
    for csv in folder.glob('measure_*.csv'):
        if csv.is_dir():
            continue
        name_lower = csv.name.lower()
        is_cis = 'a2780cis' in name_lower
        target_dir = cis_dir if is_cis else non_dir
        target_path = target_dir / csv.name
        action = 'copied'
        if MOVE_FILES:
            shutil.move(str(csv), target_path)
            action = 'moved'
        else:
            shutil.copy2(csv, target_path)
        records.append({
            'file': csv.name,
            'density': density,
            'cis': is_cis,
            'destination': str(target_path),
            'action': action,
        })
    return records


def main():
    all_records = []
    for d in DENSITIES:
        all_records.extend(classify_and_place(d))
    df = pd.DataFrame(all_records)
    out = BASE / 'cis_split_summary.csv'
    df.to_csv(out, index=False)
    print(f"Summary written: {out}")
    print(df[['density','file','cis','action']])

if __name__ == '__main__':
    main()
