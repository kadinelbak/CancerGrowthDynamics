"""Export the numeric trajectories in the low-resource Excel workbooks.

This intentionally uses only the Python standard library: it makes the source
workbooks auditable without imposing an additional Python dependency on users.
"""

from __future__ import annotations

import csv
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
REL = "http://schemas.openxmlformats.org/package/2006/relationships"


def col_index(ref: str) -> int:
    letters = re.match(r"[A-Z]+", ref).group(0)
    answer = 0
    for letter in letters:
        answer = answer * 26 + ord(letter) - 64
    return answer - 1


def workbook_rows(path: Path):
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared = ["".join(node.itertext()) for node in root.findall("m:si", NS)]
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        targets = {item.attrib["Id"]: item.attrib["Target"] for item in rels.findall(f"{{{REL}}}Relationship")}
        book = ET.fromstring(z.read("xl/workbook.xml"))
        for sheet in book.findall("m:sheets/m:sheet", NS):
            target = targets[sheet.attrib[f"{{{NS['r']}}}id"]]
            xml_path = "xl/" + target.lstrip("/")
            root = ET.fromstring(z.read(xml_path))
            rows = []
            for row in root.findall("m:sheetData/m:row", NS):
                values = {}
                for cell in row.findall("m:c", NS):
                    raw = cell.find("m:v", NS)
                    value = "" if raw is None else raw.text
                    if cell.attrib.get("t") == "s" and value:
                        value = shared[int(value)]
                    elif cell.attrib.get("t") == "inlineStr":
                        value = "".join(cell.itertext())
                    values[col_index(cell.attrib["r"])] = value
                if values:
                    rows.append([values.get(i, "") for i in range(max(values) + 1)])
            yield sheet.attrib["name"], rows


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def export_workbook(path: Path, out):
    for sheet, rows in workbook_rows(path):
        # Data sheets use a time column followed by one or more numeric series.
        header_row = next((i for i, row in enumerate(rows) if any(isinstance(cell, str) and "day" in cell.lower() for cell in row)), None)
        if header_row is None:
            continue
        headers = rows[header_row]
        for row in rows[header_row + 1:]:
            if not row or number(row[0] if row else None) is None:
                continue
            day = number(row[0])
            for col, heading in enumerate(headers[1:], start=1):
                if not heading or col >= len(row):
                    continue
                value = number(row[col])
                if value is not None:
                    out.writerow([path.name, sheet, heading, day, value])


def main():
    source = Path(sys.argv[1])
    destination = Path(sys.argv[2])
    files = sorted(source.glob("*.xlsx"))
    with destination.open("w", newline="", encoding="utf-8") as handle:
        out = csv.writer(handle)
        out.writerow(["workbook", "sheet", "series", "day", "value"])
        for workbook in files:
            export_workbook(workbook, out)


if __name__ == "__main__":
    main()
