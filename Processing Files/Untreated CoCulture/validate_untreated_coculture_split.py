"""Validate that the 20k and 30k split CSVs contain only the expected wells.

20k expected wells: A1-3, B1-3, C1-3
30k expected wells: A4-6, B4-6, C4-6

Outputs a CSV summary and prints anomalies.
"""
from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

ROOT = Path(r"c:/Users/MainFrameTower/Desktop/CancerGrowthDynamics")
BASE = ROOT / "Processed_Datasets" / "Untreated CoCulture"
WELL_PATTERN = re.compile(r"_([ABC][1-6])\.tif$", re.IGNORECASE)
EXPECTED_20K = {f"{r}{c}" for r in "ABC" for c in (1,2,3)}
EXPECTED_30K = {f"{r}{c}" for r in "ABC" for c in (4,5,6)}

records = []
for density, expected in [("20k", EXPECTED_20K), ("30k", EXPECTED_30K)]:
    folder = BASE / density
    for csv in sorted(folder.glob('measure_*.csv')):
        df = pd.read_csv(csv)
        if 'Image' not in df.columns:
            records.append({
                'density': density,
                'file': csv.name,
                'rows': len(df),
                'missing_image_col': True,
                'unique_wells': None,
                'invalid_wells': None,
            })
            continue
        wells = df['Image'].str.extract(WELL_PATTERN, expand=False).str.upper()
        unique = sorted(w for w in wells.dropna().unique())
        invalid = sorted(set(unique) - expected)
        counts = wells.value_counts().to_dict()
        records.append({
            'density': density,
            'file': csv.name,
            'rows': len(df),
            'unique_wells': ' '.join(unique),
            'invalid_wells': ' '.join(invalid) if invalid else '',
            'per_well_counts': counts,
        })

summary = pd.DataFrame(records)
summary_path = BASE / 'validation_split_summary.csv'
summary.to_csv(summary_path, index=False)
print(f"Validation summary written to {summary_path}")

if 'invalid_wells' in summary.columns:
    anomalies = summary[summary['invalid_wells'].astype(str).str.len() > 0]
    if not anomalies.empty:
        print("Anomalies detected:")
        print(anomalies[['density','file','invalid_wells']])
    else:
        print("No invalid wells detected.")
else:
    print("Column 'invalid_wells' missing in summary (no records processed?)")
