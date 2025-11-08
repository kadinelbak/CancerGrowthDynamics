"""Compute per-day averages across wells for untreated co-culture datasets.

Uses the existing per-file well_day averages produced in each:
  Processed_Datasets/Untreated CoCulture/<density>/<cis_status>/Averages/*_well_day_averages.csv

For each such file:
  - Read Day, Well, mean_cells
  - Compute day_mean_cells = mean(mean_cells over wells for that day)
  - Count wells contributing (n_wells)
  - Round day_mean_cells to 2 decimals
  - Output: <basename>_day_averages.csv in same Averages folder

Also creates a combined file per (density, cis_status):
  combined_day_averages.csv with columns:
    source_file, Day, day_mean_cells, n_wells

Skips files without at least one Day.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

def find_repo_root(start: Path) -> Path:
    cur = start
    for _ in range(10):
        if (cur / "Processed_Datasets").exists() or (cur / "Datasets").exists():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return start

DENSITIES = ["20k","30k"]
CIS_STATES = ["cis","non_cis"]


def process_well_day_file(path: Path):
    df = pd.read_csv(path)
    required = {'Day','Well','mean_cells'}
    if not required.issubset(df.columns):
        return None
    # Ensure Day numeric
    df['Day'] = pd.to_numeric(df['Day'], errors='coerce')
    grp = df.dropna(subset=['Day']).groupby('Day', as_index=False).agg(
        day_mean_cells=('mean_cells','mean'),
        n_wells=('Well','nunique')
    )
    if grp.empty:
        return None
    grp['day_mean_cells'] = grp['day_mean_cells'].round(2)
    grp = grp.sort_values('Day')
    out_path = path.parent / (path.stem.replace('_well_day_averages','') + '_day_averages.csv')
    grp.to_csv(out_path, index=False)
    return out_path, grp


def main():
    repo_root = find_repo_root(Path(__file__).resolve().parent)
    base = repo_root / "Processed_Datasets" / "Untreated CoCulture"
    summary_rows = []
    for density in DENSITIES:
        for state in CIS_STATES:
            avg_dir = base / density / state / 'Averages'
            if not avg_dir.exists():
                continue
            combined_parts = []
            for f in sorted(avg_dir.glob('*_well_day_averages.csv')):
                res = process_well_day_file(f)
                if not res:
                    continue
                out_path, grp = res
                grp2 = grp.copy()
                grp2.insert(0,'source_file', f.name)
                combined_parts.append(grp2)
                summary_rows.append({
                    'density': density,
                    'cis_status': state,
                    'source_file': f.name,
                    'output_file': out_path.name,
                    'days': grp['Day'].nunique(),
                    'mean_first_day': grp.loc[grp['Day'].idxmin(),'day_mean_cells']
                })
            if combined_parts:
                combined = pd.concat(combined_parts, ignore_index=True)
                combined.to_csv(avg_dir / 'combined_day_averages.csv', index=False)
    summary_df = pd.DataFrame(summary_rows)
    summary_path = base / 'day_averages_generation_summary.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f"Summary written: {summary_path}")
    if not summary_df.empty:
        print(summary_df.head())

if __name__ == '__main__':
    main()
