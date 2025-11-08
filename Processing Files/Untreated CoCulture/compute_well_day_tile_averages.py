"""Compute per-well per-day tile averages for untreated co-culture processed datasets.

For each density (20k, 30k) and cis status (cis, non_cis), read each measure_*.csv file, parse:
  Image pattern: <...>_Day<day>_<CHANNEL>_Tile-<tile>_<Well>.tif
Where <CHANNEL> can contain spaces (e.g., 'FITC', 'TX RED').
Extract Day (int), Channel (str), Tile (int), Well (e.g., A3) and compute mean Cells over tiles for each (Day, Well).
NOTE: If multiple channels exist for the same Day/Well/Tile, their values are both included in the averaging; n_tiles reflects the count of unique tile indices (not multiplied by channels). This effectively averages across channels. Adjust if channel-specific summaries are desired.
Outputs:
  <density>/<cis|non_cis>/Averages/<basename>_well_day_averages.csv
Also writes a combined file per density/cis status consolidating all input files:
  <density>/<cis|non_cis>/Averages/combined_well_day_averages.csv with columns:
    source_file, Day, Well, mean_cells, n_tiles
Rounding: mean_cells rounded to 2 decimal places.
Ordering: Rows sorted by Day ascending then Well alphanumerically (A1..A6,B1..,C6).
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

DENSITIES = ["20k", "30k"]
CIS_STATES = ["cis", "non_cis"]
# Capture day, channel (any letters/spaces), tile, well
IMG_RE = re.compile(r"Day(\d+)_([A-Za-z ]+)_Tile-(\d+)_([ABC][1-6])\.tif$", re.IGNORECASE)


def parse_row(image: str):
    m = IMG_RE.search(image)
    if not m:
        return None, None, None, None
    day = int(m.group(1))
    channel = m.group(2).strip()
    tile = int(m.group(3))
    well = m.group(4).upper()
    return day, channel, tile, well


def sort_well_key(well: str):
    if not isinstance(well, str) or len(well) < 2:
        return (99, 99)
    row = well[0]
    col = well[1:]
    row_order = {'A':0,'B':1,'C':2}.get(row, 99)
    try:
        col_num = int(col)
    except ValueError:
        col_num = 99
    return (row_order, col_num)


def process_file(csv_path: Path, out_dir: Path):
    df = pd.read_csv(csv_path)
    if 'Image' not in df.columns or 'Cells' not in df.columns:
        return None
    parsed = df['Image'].map(parse_row)
    df[['Day','Channel','Tile','Well']] = pd.DataFrame(parsed.tolist(), index=df.index)
    missing = df['Day'].isna().sum()
    # Group ignoring channel (averaging across available channels per well/day)
    grp = (df.dropna(subset=['Day','Well'])
             .groupby(['Day','Well'], as_index=False)
             .agg(mean_cells=('Cells','mean'),
                  n_tiles=('Tile','nunique'),
                  channel_count=('Channel','nunique')))
    grp['mean_cells'] = grp['mean_cells'].round(2)
    grp = grp.sort_values(by=['Day','Well'], key=lambda col: [sort_well_key(v) if col.name=='Well' else v for v in col])
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{csv_path.stem}_well_day_averages.csv"
    grp.to_csv(out_file, index=False)
    return {
        'file': csv_path.name,
        'rows': len(df),
        'parsed_rows': len(df) - missing,
        'unique_days': grp['Day'].nunique(),
        'unique_wells': grp['Well'].nunique(),
        'output': str(out_file)
    }


def main():
    repo_root = find_repo_root(Path(__file__).resolve().parent)
    base = repo_root / "Processed_Datasets" / "Untreated CoCulture"
    summary_records = []
    for density in DENSITIES:
        for state in CIS_STATES:
            folder = base / density / state
            if not folder.exists():
                continue
            averages_dir = folder / 'Averages'
            combined_rows = []
            for csv in sorted(folder.glob('measure_*.csv')):
                rec = process_file(csv, averages_dir)
                if rec:
                    summary_records.append({**rec, 'density': density, 'cis_status': state})
                    per_file_avg = pd.read_csv(rec['output'])
                    per_file_avg.insert(0, 'source_file', csv.name)
                    combined_rows.append(per_file_avg)
            if combined_rows:
                combined = pd.concat(combined_rows, ignore_index=True)
                combined = combined.sort_values(by=['source_file','Day','Well'], key=lambda col: [sort_well_key(v) if col.name=='Well' else v for v in col])
                combined.to_csv(averages_dir / 'combined_well_day_averages.csv', index=False)
    summary_df = pd.DataFrame(summary_records)
    summary_path = base / 'averages_generation_summary.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f"Summary written: {summary_path}")
    if not summary_df.empty:
        print(summary_df[['density','cis_status','file','unique_days','unique_wells']])

if __name__ == '__main__':
    main()
