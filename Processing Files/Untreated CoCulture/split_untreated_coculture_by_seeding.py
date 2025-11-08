"""Split untreated co-culture processed CSVs into 20k and 30k subsets based on well IDs.

Rule:
  Wells A1-A3, B1-B3, C1-C3  -> 20k
  Wells A4-A6, B4-B6, C4-C6  -> 30k

We identify wells from the image filename embedded in the 'Image' column, expecting a pattern ending with _A3.tif etc.
Outputs are written to:
  Processed_Datasets/Untreated CoCulture/20k/<original_basename> (filtered)
  Processed_Datasets/Untreated CoCulture/30k/<original_basename>

If a file contains rows for both groups both output files are produced; if only one group is present only that file is written.
"""
from __future__ import annotations
import re
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

WELLS_20K = {f"{row}{col}" for row in "ABC" for col in (1,2,3)}  # A1..C3
WELLS_30K = {f"{row}{col}" for row in "ABC" for col in (4,5,6)}  # A4..C6
WELL_PATTERN = re.compile(r"_([ABC][1-6])\.tif$", re.IGNORECASE)

def extract_well(image_name: str) -> str | None:
    m = WELL_PATTERN.search(image_name)
    return m.group(1).upper() if m else None


def process_file(path: Path) -> dict:
    df = pd.read_csv(path)
    if 'Image' not in df.columns:
        print(f"[WARN] Skipping {path.name}: missing 'Image' column")
        return {"file": str(path), "status": "missing-image-col"}
    wells = df['Image'].map(extract_well)
    df = df.assign(__Well=wells)
    # Partition
    df_20 = df[df['__Well'].isin(WELLS_20K)].drop(columns=['__Well'])
    df_30 = df[df['__Well'].isin(WELLS_30K)].drop(columns=['__Well'])

    written = []
    if not df_20.empty:
        out_path = OUT_20K / path.name
        df_20.to_csv(out_path, index=False)
        written.append(str(out_path))
    if not df_30.empty:
        out_path = OUT_30K / path.name
        df_30.to_csv(out_path, index=False)
        written.append(str(out_path))

    summary = {
        "file": str(path),
        "rows_total": len(df),
        "rows_20k": len(df_20),
        "rows_30k": len(df_30),
        "written": written or None,
    }
    print(summary)
    return summary


def main():
    repo_root = find_repo_root(Path(__file__).resolve().parent)
    src_dir = repo_root / "Processed_Datasets" / "Untreated CoCulture"
    out_20k = src_dir / "20k"
    out_30k = src_dir / "30k"
    out_20k.mkdir(parents=True, exist_ok=True)
    out_30k.mkdir(parents=True, exist_ok=True)

    # rebind module-level outputs for process_file
    global OUT_20K, OUT_30K
    OUT_20K, OUT_30K = out_20k, out_30k

    csv_files = [p for p in src_dir.glob('measure_*.csv') if p.is_file() and p.name not in ('20k', '30k')]
    results = [process_file(p) for p in csv_files]
    # Simple report
    report_df = pd.DataFrame(results)
    report_path = src_dir / 'split_coculture_report.csv'
    report_df.to_csv(report_path, index=False)
    print(f"Report written to {report_path}")

if __name__ == '__main__':
    main()
