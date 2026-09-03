from __future__ import annotations

import csv
import html
import json
import math
import re
import statistics
import zipfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
from xml.etree import ElementTree as ET


BASE_DIR = Path(r"C:\Users\elbak\OneDrive\Desktop\Research\Results\CancerGrowthDynamics\New Datasets")
OUT_DIR = Path(r"C:\Users\elbak\OneDrive\Desktop\Research\Results\CancerGrowthDynamics\outputs\new_data_report")
OUT_DIR.mkdir(parents=True, exist_ok=True)

NS = {
    "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "p": "http://schemas.openxmlformats.org/package/2006/relationships",
}


def to_number(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, str) and value.startswith("#"):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if text == "":
        return None
    try:
        if re.fullmatch(r"[-+]?\d+", text):
            return int(text)
        return float(text)
    except Exception:
        return value


def format_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isfinite(value):
            if abs(value) >= 1000:
                return f"{value:,.0f}"
            if abs(value) >= 100:
                return f"{value:,.1f}"
            if abs(value) >= 10:
                return f"{value:,.2f}"
            return f"{value:.3f}".rstrip("0").rstrip(".")
        return str(value)
    return str(value)


def cell_ref_to_col(ref: str) -> int:
    col = 0
    for ch in ref:
        if ch.isalpha():
            col = col * 26 + (ord(ch.upper()) - 64)
        else:
            break
    return col


@dataclass
class SheetData:
    source_file: str
    sheet_name: str
    rows: list[dict[str, Any]]
    headers: list[str]
    kind: str
    notes: str

    @property
    def preview_rows(self) -> list[dict[str, Any]]:
        return self.rows[:3]


def parse_csv(path: Path) -> SheetData:
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        rows = list(reader)
    headers = [h.strip() for h in (rows[0] if rows else [])]
    data: list[dict[str, Any]] = []
    for raw in rows[1:]:
        row = {}
        for i, header in enumerate(headers):
            value = raw[i] if i < len(raw) else ""
            row[header] = to_number(value)
        data.append(row)
    return SheetData(
        source_file=path.name,
        sheet_name="CSV",
        rows=data,
        headers=headers,
        kind="Mono-culture time course",
        notes="Single-sheet CSV time course at 30,000 cells/mL.",
    )


def parse_xlsx(path: Path) -> list[SheetData]:
    out: list[SheetData] = []
    with zipfile.ZipFile(path) as zf:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("a:si", NS):
                text = "".join(t.text or "" for t in si.iterfind(".//a:t", NS))
                shared_strings.append(text)

        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_targets: dict[str, str] = {}
        for rel in rels.findall("p:Relationship", NS):
            rel_targets[rel.attrib["Id"]] = rel.attrib["Target"]

        for sheet in wb.find("a:sheets", NS):
            sheet_name = sheet.attrib["name"]
            rel_id = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
            target = rel_targets[rel_id]
            if not target.startswith("xl/"):
                target = "xl/" + target
            ws_root = ET.fromstring(zf.read(target))

            header_map: OrderedDict[int, str] = OrderedDict()
            rows_data: list[dict[str, Any]] = []
            for row in ws_root.findall(".//a:sheetData/a:row", NS):
                row_map: dict[str, Any] = {}
                for cell in row.findall("a:c", NS):
                    ref = cell.attrib.get("r", "")
                    col_idx = cell_ref_to_col(ref)
                    cell_type = cell.attrib.get("t")
                    v = cell.find("a:v", NS)
                    value: Any = None
                    if cell_type == "s" and v is not None:
                        idx = int(v.text)
                        value = shared_strings[idx] if idx < len(shared_strings) else ""
                    elif cell_type == "inlineStr":
                        isel = cell.find("a:is", NS)
                        value = "".join(t.text or "" for t in isel.iterfind(".//a:t", NS)) if isel is not None else ""
                    else:
                        value = v.text if v is not None else None

                    if row.attrib.get("r") == "1":
                        header_map[col_idx] = str(value).strip() if value is not None else ""
                    else:
                        header = header_map.get(col_idx, f"Column {col_idx}")
                        row_map[header] = to_number(value)
                if row.attrib.get("r") != "1":
                    rows_data.append(row_map)

            headers = list(header_map.values())
            out.append(
                SheetData(
                    source_file=path.name,
                    sheet_name=sheet_name,
                    rows=rows_data,
                    headers=headers,
                    kind="Low-resource co-culture workbook" if "LR" in path.name else "Workbook",
                    notes="Sheet-level time course with day, sample count, mean cells, SD, and SEM.",
                )
            )
    return out


def classify_file(name: str) -> str:
    lower = name.lower()
    if lower.endswith(".csv"):
        if "untreated" in lower:
            return "Mono-culture untreated"
        if "ic25" in lower:
            return "Mono-culture treated, IC25"
        if "ic50" in lower:
            return "Mono-culture treated, IC50"
        if "ic75" in lower:
            return "Mono-culture treated, IC75"
        return "Mono-culture"
    if lower.startswith("ce0_"):
        return "Co-culture untreated (ce0)"
    if lower.startswith("ce1_"):
        return "Co-culture treated at IC50, 1 uM (ce1)"
    if "lr data" in lower:
        return "Low-resource workbook"
    return "Workbook"


def load_all_data() -> list[SheetData]:
    sheets: list[SheetData] = []
    for path in sorted(BASE_DIR.iterdir()):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        if suffix == ".csv":
            sheets.append(parse_csv(path))
        elif suffix == ".xlsx":
            sheets.extend(parse_xlsx(path))
    return sheets


def numeric_rows(sheet: SheetData) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for row in sheet.rows:
        day = to_number(row.get("Day")) or to_number(row.get("Day "))
        mean = to_number(row.get("Mean Cells"))
        sem = to_number(row.get("SEM Cells"))
        sd = to_number(row.get("SD Cells"))
        n = to_number(row.get("N Samples"))
        if day is None or mean is None:
            continue
        try:
            day_f = float(day)
            mean_f = float(mean)
            sem_f = float(sem) if sem is not None else None
            sd_f = float(sd) if sd is not None else None
            n_f = float(n) if n is not None else None
        except (TypeError, ValueError):
            continue
        rows.append(
            {
                "Day": day_f,
                "Mean Cells": mean_f,
                "SEM Cells": sem_f,
                "SD Cells": sd_f,
                "N Samples": n_f,
            }
        )
    return rows


def svg_escape(text: Any) -> str:
    return html.escape(format_value(text), quote=True)


def make_svg_chart(sheet: SheetData) -> str:
    rows = numeric_rows(sheet)
    if not rows:
        return "<div class='empty-chart'>No numeric data found.</div>"

    days = [r["Day"] for r in rows]
    means = [r["Mean Cells"] for r in rows]
    errs = [r["SEM Cells"] for r in rows]

    w, h = 920, 420
    m = {"l": 78, "r": 28, "t": 68, "b": 72}
    cw = w - m["l"] - m["r"]
    ch = h - m["t"] - m["b"]
    xmin, xmax = min(days), max(days)
    if xmin == xmax:
        xmin -= 1
        xmax += 1
    ymin = 0.0
    ymax = max(m for m in means if m is not None)
    if any(e is not None for e in errs):
        ymax = max(ymax, max((m + e) for m, e in zip(means, errs) if e is not None))
    ymax = ymax * 1.15 if ymax > 0 else 1.0
    if ymax == 0:
        ymax = 1.0

    def x_map(x: float) -> float:
        return m["l"] + (x - xmin) / (xmax - xmin) * cw

    def y_map(y: float) -> float:
        return m["t"] + (1 - (y - ymin) / (ymax - ymin)) * ch

    grid_lines = []
    y_ticks = 5
    x_ticks = min(8, len(days))
    for i in range(y_ticks + 1):
        frac = i / y_ticks
        y = ymin + (ymax - ymin) * frac
        yy = y_map(y)
        grid_lines.append(
            f"<line x1='{m['l']}' y1='{yy:.1f}' x2='{w - m['r']}' y2='{yy:.1f}' class='grid'/>"
            f"<text x='{m['l'] - 12}' y='{yy + 4:.1f}' class='axis-label axis-left'>{format_value(y)}</text>"
        )

    x_positions = sorted(set(days))
    for x in x_positions:
        xx = x_map(x)
        grid_lines.append(
            f"<line x1='{xx:.1f}' y1='{m['t']}' x2='{xx:.1f}' y2='{h - m['b']}' class='grid-v'/>"
            f"<text x='{xx:.1f}' y='{h - m['b'] + 22}' text-anchor='middle' class='axis-label'>{format_value(x)}</text>"
        )

    points = []
    path_pts = []
    for d, mean, err in zip(days, means, errs):
        xx = x_map(d)
        yy = y_map(mean)
        path_pts.append(f"{xx:.1f},{yy:.1f}")
        points.append(f"<circle cx='{xx:.1f}' cy='{yy:.1f}' r='4.2' class='point'/>")
        if err is not None:
            lo = max(ymin, mean - err)
            hi = min(ymax, mean + err)
            points.append(f"<line x1='{xx:.1f}' y1='{y_map(lo):.1f}' x2='{xx:.1f}' y2='{y_map(hi):.1f}' class='error-bar'/>")
            points.append(f"<line x1='{xx - 6:.1f}' y1='{y_map(lo):.1f}' x2='{xx + 6:.1f}' y2='{y_map(lo):.1f}' class='error-cap'/>")
            points.append(f"<line x1='{xx - 6:.1f}' y1='{y_map(hi):.1f}' x2='{xx + 6:.1f}' y2='{y_map(hi):.1f}' class='error-cap'/>")

    line = "<polyline points='" + " ".join(path_pts) + "' class='series'/>"
    title = f"{sheet.source_file} | {sheet.sheet_name}"
    subtitle = sheet.notes

    return f"""
    <svg viewBox="0 0 {w} {h}" class="chart" role="img" aria-label="{html.escape(title)}">
      <rect x="0" y="0" width="{w}" height="{h}" rx="18" class="chart-bg"/>
      <text x="30" y="34" class="chart-title">{html.escape(title)}</text>
      <text x="30" y="56" class="chart-subtitle">{html.escape(subtitle)}</text>
      <line x1="{m['l']}" y1="{h - m['b']}" x2="{w - m['r']}" y2="{h - m['b']}" class="axis"/>
      <line x1="{m['l']}" y1="{m['t']}" x2="{m['l']}" y2="{h - m['b']}" class="axis"/>
      {''.join(grid_lines)}
      {line}
      {''.join(points)}
      <text x="{w/2:.1f}" y="{h - 24}" text-anchor="middle" class="axis-title">Day</text>
      <text transform="translate(24 {h/2:.1f}) rotate(-90)" text-anchor="middle" class="axis-title">Mean Cells</text>
    </svg>
    """


def make_table(sheet: SheetData) -> str:
    headers = sheet.headers[:]
    preview = sheet.preview_rows
    if not headers and preview:
        headers = list(preview[0].keys())
    head_html = "".join(f"<th>{html.escape(str(h))}</th>" for h in headers)
    body_rows = []
    for row in preview:
        cells = []
        for h in headers:
            cells.append(f"<td>{html.escape(format_value(row.get(h)))}</td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    if not body_rows:
        body_rows.append(f"<tr><td colspan='{max(1, len(headers))}'>No preview rows available.</td></tr>")
    return f"""
    <table class="preview-table">
      <thead><tr>{head_html}</tr></thead>
      <tbody>{''.join(body_rows)}</tbody>
    </table>
    """


def make_manifest(sheets: list[SheetData]) -> list[dict[str, Any]]:
    files: OrderedDict[str, list[SheetData]] = OrderedDict()
    for sheet in sheets:
        files.setdefault(sheet.source_file, []).append(sheet)
    manifest = []
    for file_name, file_sheets in files.items():
        manifest.append(
            {
                "file": file_name,
                "classification": classify_file(file_name),
                "sheet_count": len(file_sheets),
                "sheets": [s.sheet_name for s in file_sheets],
            }
        )
    return manifest


def build_html(sheets: list[SheetData]) -> str:
    manifest = make_manifest(sheets)
    total_files = len(manifest)
    total_sheets = len(sheets)
    total_points = sum(len(numeric_rows(s)) for s in sheets)

    grouped: OrderedDict[str, list[SheetData]] = OrderedDict()
    for sheet in sheets:
        grouped.setdefault(sheet.source_file, []).append(sheet)

    sections = []
    for file_name, file_sheets in grouped.items():
        classification = classify_file(file_name)
        cards = []
        for sheet in file_sheets:
            cards.append(
                f"""
                <section class="card">
                  <div class="card-meta">
                    <div class="card-kicker">{html.escape(classification)}</div>
                    <div class="card-caption">{html.escape(file_name)} | {html.escape(sheet.sheet_name)}</div>
                    <div class="card-note">{html.escape(sheet.notes)}</div>
                  </div>
                  {make_svg_chart(sheet)}
                  <div class="preview-wrap">
                    <div class="preview-label">Sheet preview</div>
                    {make_table(sheet)}
                  </div>
                </section>
                """
            )
        sections.append(
            f"""
            <details class="workbook" open>
              <summary>
                <span class="summary-title">{html.escape(file_name)}</span>
                <span class="summary-subtitle">{html.escape(classification)} · {len(file_sheets)} sheet(s)</span>
              </summary>
              <div class="card-grid">{''.join(cards)}</div>
            </details>
            """
        )

    manifest_rows = []
    for item in manifest:
        manifest_rows.append(
            "<tr>"
            f"<td>{html.escape(item['file'])}</td>"
            f"<td>{html.escape(item['classification'])}</td>"
            f"<td>{item['sheet_count']}</td>"
            f"<td>{html.escape(', '.join(item['sheets']))}</td>"
            "</tr>"
        )

    html_doc = f"""
    <!doctype html>
    <html lang="en">
    <head>
      <meta charset="utf-8"/>
      <meta name="viewport" content="width=device-width, initial-scale=1"/>
      <title>New Data Report</title>
      <style>
        :root {{
          --bg: #f4f1ea;
          --panel: #fffdf8;
          --ink: #1f2937;
          --muted: #5b6472;
          --line: #d8d1c3;
          --accent: #1f6f78;
          --accent-2: #8a5a44;
          --grid: #e7dfd3;
          --shadow: 0 10px 30px rgba(20, 27, 38, 0.08);
        }}
        * {{ box-sizing: border-box; }}
        body {{
          margin: 0;
          font-family: "Segoe UI", "Trebuchet MS", Arial, sans-serif;
          background:
            radial-gradient(circle at top left, rgba(31, 111, 120, 0.10), transparent 26%),
            radial-gradient(circle at top right, rgba(138, 90, 68, 0.10), transparent 22%),
            var(--bg);
          color: var(--ink);
        }}
        .wrap {{
          max-width: 1480px;
          margin: 0 auto;
          padding: 28px 22px 56px;
        }}
        header {{
          display: grid;
          gap: 12px;
          margin-bottom: 22px;
        }}
        h1 {{
          margin: 0;
          font-size: 2.2rem;
          letter-spacing: -0.02em;
        }}
        .lede {{
          max-width: 1000px;
          line-height: 1.5;
          color: var(--muted);
          margin: 0;
          font-size: 1.03rem;
        }}
        .stats {{
          display: flex;
          flex-wrap: wrap;
          gap: 10px;
          margin-top: 4px;
        }}
        .pill {{
          padding: 8px 12px;
          border-radius: 999px;
          background: rgba(255,255,255,0.78);
          border: 1px solid var(--line);
          box-shadow: var(--shadow);
          font-size: 0.92rem;
        }}
        .section-note {{
          background: rgba(255,255,255,0.75);
          border: 1px solid var(--line);
          border-radius: 18px;
          padding: 16px 18px;
          box-shadow: var(--shadow);
          margin-bottom: 18px;
          line-height: 1.55;
        }}
        .section-note strong {{ color: var(--accent); }}
        table.meta, table.preview-table {{
          width: 100%;
          border-collapse: collapse;
          background: rgba(255,255,255,0.92);
          border: 1px solid var(--line);
          border-radius: 14px;
          overflow: hidden;
          box-shadow: var(--shadow);
        }}
        table.meta th, table.meta td, table.preview-table th, table.preview-table td {{
          padding: 8px 10px;
          border-bottom: 1px solid #ece6da;
          vertical-align: top;
          font-size: 0.88rem;
        }}
        table.meta th, table.preview-table th {{
          background: #f8f4eb;
          text-align: left;
          position: sticky;
          top: 0;
        }}
        details.workbook {{
          margin: 18px 0 26px;
          border: 1px solid var(--line);
          border-radius: 20px;
          background: rgba(255,255,255,0.55);
          overflow: hidden;
          box-shadow: var(--shadow);
        }}
        details.workbook > summary {{
          list-style: none;
          cursor: pointer;
          padding: 16px 18px;
          display: flex;
          flex-direction: column;
          gap: 4px;
          background: linear-gradient(90deg, rgba(31,111,120,0.06), rgba(138,90,68,0.06));
        }}
        details.workbook > summary::-webkit-details-marker {{ display: none; }}
        .summary-title {{
          font-size: 1.1rem;
          font-weight: 700;
        }}
        .summary-subtitle {{
          color: var(--muted);
          font-size: 0.92rem;
        }}
        .card-grid {{
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(460px, 1fr));
          gap: 18px;
          padding: 18px;
        }}
        .card {{
          background: var(--panel);
          border: 1px solid var(--line);
          border-radius: 18px;
          padding: 16px;
          box-shadow: var(--shadow);
          display: grid;
          gap: 14px;
          align-content: start;
        }}
        .card-meta {{
          display: grid;
          gap: 4px;
        }}
        .card-kicker {{
          color: var(--accent-2);
          font-weight: 700;
          font-size: 0.82rem;
          text-transform: uppercase;
          letter-spacing: 0.08em;
        }}
        .card-caption {{
          font-size: 1.02rem;
          font-weight: 700;
        }}
        .card-note {{
          color: var(--muted);
          font-size: 0.92rem;
          line-height: 1.45;
        }}
        svg.chart {{
          width: 100%;
          height: auto;
          display: block;
        }}
        .chart-bg {{
          fill: #fffefa;
          stroke: var(--line);
          stroke-width: 1;
        }}
        .chart-title {{
          font-size: 17px;
          font-weight: 700;
          fill: var(--ink);
        }}
        .chart-subtitle {{
          font-size: 12px;
          fill: var(--muted);
        }}
        .axis {{
          stroke: #2b3340;
          stroke-width: 1.2;
        }}
        .grid {{
          stroke: var(--grid);
          stroke-width: 1;
          stroke-dasharray: 3 4;
        }}
        .grid-v {{
          stroke: #eee7db;
          stroke-width: 1;
          stroke-dasharray: 3 5;
        }}
        .axis-label {{
          fill: #667085;
          font-size: 11px;
        }}
        .axis-left {{
          text-anchor: end;
        }}
        .axis-title {{
          fill: #2b3340;
          font-size: 13px;
          font-weight: 700;
        }}
        .series {{
          fill: none;
          stroke: var(--accent);
          stroke-width: 2.6;
        }}
        .point {{
          fill: var(--accent-2);
          stroke: white;
          stroke-width: 1.4;
        }}
        .error-bar {{
          stroke: rgba(138,90,68,0.85);
          stroke-width: 1.6;
        }}
        .error-cap {{
          stroke: rgba(138,90,68,0.85);
          stroke-width: 1.6;
        }}
        .preview-wrap {{
          display: grid;
          gap: 8px;
        }}
        .preview-label {{
          font-weight: 700;
          color: var(--accent);
          font-size: 0.92rem;
        }}
        .empty-chart {{
          padding: 20px;
          border: 1px dashed var(--line);
          border-radius: 16px;
          color: var(--muted);
        }}
        .footnote {{
          color: var(--muted);
          font-size: 0.88rem;
          line-height: 1.5;
          margin-top: 10px;
        }}
        @media (max-width: 720px) {{
          .card-grid {{
            grid-template-columns: 1fr;
            padding: 12px;
          }}
          .wrap {{ padding: 18px 12px 40px; }}
        }}
      </style>
    </head>
    <body>
      <div class="wrap">
        <header>
          <h1>New Data Report</h1>
          <p class="lede">
            This report reads every CSV and XLSX file in <code>New Datasets</code>, plots each source with a title that
            includes the file and sheet name, and shows a small sheet preview underneath. The labeling follows the
            experiment naming in your notes: <strong>LR</strong> means low resource, <strong>ce0</strong> is untreated,
            <strong>ce1</strong> is treated at IC50 (1 uM), and <strong>1-1</strong>, <strong>1-3</strong>, and
            <strong>3-1</strong> are co-culture sensitive:resistant ratios.
          </p>
          <div class="stats">
            <div class="pill">{total_files} files</div>
            <div class="pill">{total_sheets} plotted sheets / tables</div>
            <div class="pill">{total_points} day-level observations</div>
          </div>
        </header>

        <div class="section-note">
          <strong>How to read the labels:</strong> the mono-culture files are the standalone CSV time courses, while the
          co-culture workbooks either split into untreated/treated paired sheets or contain many LR sheets for different
          ratios and cell types. Every chart is labeled with the exact file and sheet name it came from, and each sheet
          preview shows the first three rows when available.
        </div>

        <section>
          <h2>File Manifest</h2>
          <table class="meta">
            <thead>
              <tr><th>File</th><th>Classification</th><th>Sheets</th><th>Sheet names</th></tr>
            </thead>
            <tbody>
              {''.join(manifest_rows)}
            </tbody>
          </table>
          <div class="footnote">
            Preview rows are extracted directly from the spreadsheet content. For XLSX files, the sheet names and the
            first row headers are read from the workbook XML so the report reflects the source structure rather than a
            hand-entered summary.
          </div>
        </section>

        {''.join(sections)}
      </div>
    </body>
    </html>
    """
    return html_doc


def main() -> None:
    sheets = load_all_data()
    html_doc = build_html(sheets)
    report_path = OUT_DIR / "report.html"
    report_path.write_text(html_doc, encoding="utf-8")

    manifest_path = OUT_DIR / "manifest.json"
    manifest_path.write_text(json.dumps(make_manifest(sheets), indent=2), encoding="utf-8")

    print(report_path)
    print(manifest_path)


if __name__ == "__main__":
    main()
