"""Reprocess raw Untreated CoCulture datasets:
- Read from Datasets/Untreated CoCulture
- Divide 'Area µm^2' by 144
- Rename column to 'Cells'
- Write to Processed_Datasets/Untreated CoCulture (overwrite existing)
- Idempotent: if 'Cells' already present and 'Area µm^2' absent, skip unless --force
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

ROOT = Path(r"c:/Users/MainFrameTower/Desktop/CancerGrowthDynamics")
RAW_DIR = ROOT / "Datasets" / "Untreated CoCulture"
PROC_DIR = ROOT / "Processed_Datasets" / "Untreated CoCulture"


def reprocess(force: bool=False):
    if not RAW_DIR.exists():
        print(f"Raw directory missing: {RAW_DIR}")
        return
    PROC_DIR.mkdir(parents=True, exist_ok=True)
    csvs = sorted(RAW_DIR.glob('measure_*.csv'))
    if not csvs:
        print("No raw untreated co-culture CSVs found.")
        return
    rows = []
    for src in csvs:
        try:
            df = pd.read_csv(src)
        except Exception as e:
            print(f"[ERROR] Failed reading {src.name}: {e}")
            continue
        status = ''
        if 'Area µm^2' in df.columns:
            df['Cells'] = df['Area µm^2'] / 144.0
            # Keep original Area? We can drop to avoid confusion
            df = df.drop(columns=['Area µm^2'])
            status = 'converted'
        elif 'Cells' in df.columns and not force:
            status = 'skipped-already-cells'
        elif 'Cells' in df.columns and force:
            # Assume needs reconversion only if an original backup existed (not tracked here)
            status = 'forced-no-area'
        else:
            status = 'no-target-column'
        out_path = PROC_DIR / src.name
        if status.startswith('converted') or (status.startswith('skipped') and not out_path.exists()):
            df.to_csv(out_path, index=False)
        elif status.startswith('converted'):
            df.to_csv(out_path, index=False)
        rows.append({'file': src.name, 'rows': len(df), 'status': status})
        print(f"{src.name}: {status}")
    rep = pd.DataFrame(rows)
    rep_path = PROC_DIR / 'reprocess_report.csv'
    rep.to_csv(rep_path, index=False)
    print(f"Report written: {rep_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--force', action='store_true', help='Force reprocess if only Cells column present')
    args = ap.parse_args()
    reprocess(force=args.force)

if __name__ == '__main__':
    main()
