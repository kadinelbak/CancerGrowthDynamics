#!/usr/bin/env python3
"""Build the interactive reader page for a manuscript draft (.docx).

Usage:
    python3 docs/papers/build_paper.py "Working Draft - equations 2.docx"

The text is converted with pandoc exactly as written; this script only adds
structure for the reader (equation blocks with info toggles, section ids,
status badges, citation links, figure zoom). Equation explanations live in
docs/papers/from-culture-to-control/equation-notes.js and are keyed by the
equation number printed in the draft, e.g. "2.14", so renumbering equations in
the draft means renumbering the keys there too.

Needs pandoc and BeautifulSoup (bs4).
"""

import datetime
import html
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

from bs4 import BeautifulSoup, NavigableString

HERE = pathlib.Path(__file__).resolve().parent
OUT_DIR = HERE / "from-culture-to-control"
TITLE = "From Culture to Control"
SUBTITLE = ("Staged Mechanistic Modeling and Optimal Dosing of Cisplatin-Sensitive and "
            "-Resistant Ovarian Cancer Cells")

STATUS_RE = re.compile(r"\s*\[(Drafted|Now|Update|Pending|Write last|Bullets|To write)(?:;[^\]]*)?\]\s*$")
EQ_NUM_RE = re.compile(r"^\((\d+\.\d+)\)$")
PLACEHOLDER_RE = re.compile(r"\((?:ADD [^)]*|LINK IT|INSERT SOURCE)\)|\{[^{}]{1,12}\}")
CITE_RE = re.compile(
    r"\b(de Pillis|de Groot|[A-Z][A-Za-zü′'-]+)"
    r"(?: et al\.| (?:&|and) [A-Z][A-Za-z-]+)?,?\s\(?((?:19|20)\d{2})([ab])?"
)


def slug(text):
    text = re.sub(r"[^\w\s.-]", "", text.lower())
    return re.sub(r"[\s.]+", "-", text).strip("-")


def pandoc(docx, media_dir):
    result = subprocess.run(
        ["pandoc", str(docx), "-t", "html5", "--mathjax", "--wrap=none",
         f"--extract-media={media_dir}"],
        check=True, capture_output=True, text=True)
    return result.stdout


def inline_math(span):
    """Turn a display-math span that sits inside running text into inline math."""
    tex = span.get_text()
    tex = tex.removeprefix("\\[").removesuffix("\\]")
    tex = re.sub(r"(\\[ ,;!]|\s)+$", "", tex)
    span.string = f"\\({tex}\\)"
    span["class"] = ["math", "inline"]


def equation_block(soup, tex, number=None, key=None):
    key = key or number
    block = soup.new_tag("div", attrs={"class": "eq", "data-eq": key, "id": f"eq-{slug(key)}"})
    body = soup.new_tag("div", attrs={"class": "eq-body"})
    body.append(NavigableString(tex))
    block.append(body)
    side = soup.new_tag("div", attrs={"class": "eq-side"})
    if number:
        num = soup.new_tag("a", attrs={"class": "eq-num", "href": f"#eq-{slug(key)}",
                                       "title": "Link to this equation"})
        num.string = f"({number})"
        side.append(num)
    button = soup.new_tag("button", attrs={
        "class": "eq-info", "type": "button", "aria-expanded": "false",
        "aria-label": f"Explain equation {number or ''}".strip(),
        "title": "Explain this equation"})
    button.string = "i"
    side.append(button)
    block.append(side)
    return block


COPIES = []  # (extracted media path, published file name), filled by transform


def transform(fragment, image_name):
    soup = BeautifulSoup(fragment, "html.parser")

    # Numbered equations: pandoc gives a one-row table of [display math | (n.n)].
    for table in soup.find_all("table"):
        rows = table.find_all("tr")
        if len(rows) != 1:
            continue
        cells = rows[0].find_all("td")
        if len(cells) != 2:
            continue
        match = EQ_NUM_RE.match(cells[1].get_text(strip=True))
        span = cells[0].find("span", class_="math")
        if match and span:
            table.replace_with(equation_block(soup, span.get_text(), match.group(1)))

    # Unnumbered display equations that stand alone in a paragraph.
    unnumbered = 0
    for p in soup.find_all("p"):
        spans = p.find_all("span", class_="display")
        if len(spans) == 1 and p.get_text(strip=True) == spans[0].get_text(strip=True):
            unnumbered += 1
            p.replace_with(equation_block(soup, spans[0].get_text(), key=f"u{unnumbered}"))

    # Display math left inside running text, tables or lists reads as inline math.
    for span in soup.find_all("span", class_="display"):
        inline_math(span)

    # Equations written inside the model-zoo tables get their own info toggles.
    for table in soup.find_all("table"):
        caption = table.find_previous_sibling("p")
        label = caption.get_text() if caption else ""
        m = re.match(r"Table (2\.2[ab])", label)
        if not m:
            continue
        table["data-table"] = m.group(1)
        for i, row in enumerate(table.find("tbody").find_all("tr"), start=1):
            cell = row.find_all("td")[1]
            key = f"{m.group(1)}-{i}"
            button = soup.new_tag("button", attrs={
                "class": "eq-info eq-info-inline", "type": "button", "data-eq": key,
                "aria-expanded": "false", "aria-label": f"Explain {row.find('td').get_text(strip=True)}",
                "title": "Explain this term"})
            button.string = "i"
            cell.append(button)
            row["id"] = f"eq-{slug(key)}"
            row["data-eq-row"] = key

    # Headings: drop the bold wrapper, give ids, turn [Status] into a badge.
    used = set()
    for h in soup.find_all(["h1", "h2", "h3", "h4"]):
        for strong in h.find_all("strong"):
            strong.unwrap()
        text = h.get_text().strip()
        if not text:
            h.decompose()
            continue
        status = None
        m = STATUS_RE.search(text)
        if m:
            status = m.group(1)
            text = text[: m.start()].strip()
        base = "sec-" + slug(text)[:60]
        hid, n = base, 2
        while hid in used:
            hid, n = f"{base}-{n}", n + 1
        used.add(hid)
        h.clear()
        h["id"] = hid
        h.append(NavigableString(text))
        if status:
            h["data-status"] = status
            badge = soup.new_tag("span", attrs={"class": f"status status-{slug(status)}"})
            badge.string = status
            h.append(NavigableString(" "))
            h.append(badge)

    # References: one id per entry, keyed by first author surname and year.
    refs = {}
    ref_head = soup.find(lambda t: t.name == "h2" and t.get_text(strip=True) == "References")
    if ref_head:
        for p in ref_head.find_next_siblings():
            if p.name != "p":
                break
            text = p.get_text()
            m = re.match(r"\s*((?:de )?[^,]+),.*?\((\d{4})([ab]?)\)", text)
            if not m:
                continue
            surname = m.group(1).strip()
            rid = f"ref-{slug(surname)}-{m.group(2)}{m.group(3)}"
            p["id"] = rid
            p["class"] = ["reference"]
            refs[(surname, m.group(2) + m.group(3))] = rid
            refs.setdefault((surname, m.group(2)), rid)

    # Citations, placeholders and inline status notes in running text.
    def rewrite_text(node):
        text = str(node)
        pieces, last = [], 0
        for m in re.finditer(f"{CITE_RE.pattern}|{PLACEHOLDER_RE.pattern}|"
                             r"\[(?:Pending|Write last|Now|Update|Drafted)[^\]]*\]", text):
            pieces.append(html.escape(text[last:m.start()], quote=False))
            whole = m.group(0)
            if m.group(1):
                rid = refs.get((m.group(1), (m.group(2) or "") + (m.group(3) or ""))) \
                    or refs.get((m.group(1), m.group(2) or ""))
                pieces.append(f'<a class="cite" href="#{rid}">{html.escape(whole)}</a>'
                              if rid else html.escape(whole))
            elif whole.startswith("["):
                pieces.append(f'<span class="draft-note">{html.escape(whole)}</span>')
            else:
                pieces.append(f'<span class="placeholder" title="Placeholder in the draft">'
                              f'{html.escape(whole)}</span>')
            last = m.end()
        if last == 0:
            return
        pieces.append(html.escape(text[last:], quote=False))
        node.replace_with(BeautifulSoup("".join(pieces), "html.parser"))

    for node in list(soup.find_all(string=True)):
        parent = node.parent
        if parent is None or parent.name in ("a", "script", "style", "h1", "h2", "h3", "h4"):
            continue
        if parent.find_parent(class_=["math", "eq-body"]) or \
                {"math", "eq-body"} & set(parent.get("class") or []):
            continue
        if ref_head and parent.find_parent("p", class_="reference"):
            continue
        rewrite_text(node)

    # Wide tables scroll sideways instead of stretching the page.
    for table in soup.find_all("table"):
        for col in table.find_all("colgroup"):
            col.decompose()
        wrap = soup.new_tag("div", attrs={"class": "table-wrap"})
        table.wrap(wrap)

    # Figures: each image is copied next to the page under a name taken from its caption
    # ("Figure 3.4." -> figure-3-4.png; Figure 2.3 keeps its pipeline name), zoomable, caption attached.
    for img in soup.find_all("img"):
        holder = img.find_parent("p") or img
        caption = holder.find_next_sibling("p")
        cap_text = caption.get_text() if caption else ""
        m = re.match(r"Figure (\d+)\.(\d+)", cap_text)
        if m is None or (m.group(1), m.group(2)) == ("2", "3"):
            name, fid, alt = image_name, "fig-2-3", "Pipeline from raw images to dosing schedules in nine steps"
        else:
            name, fid = f"figure-{m.group(1)}-{m.group(2)}.png", f"fig-{m.group(1)}-{m.group(2)}"
            alt = cap_text.split(".", 2)[-1].strip()
        COPIES.append((img["src"], name))
        img["src"] = name
        img["alt"] = alt
        img["class"] = ["zoomable"]
        del img["style"]
        figure = soup.new_tag("figure", attrs={"class": "figure", "id": fid})
        holder.replace_with(figure)
        figure.append(img)
        if caption and m is not None:
            figcap = soup.new_tag("figcaption")
            for child in list(caption.contents):
                figcap.append(child)
            caption.decompose()
            figure.append(figcap)

    for a in soup.find_all("a", href=re.compile(r"^https?://")):
        a["target"] = "_blank"
        a["rel"] = "noopener"
        a["href"] = re.sub(r"[?&]utm_source=chatgpt\.com", "", a["href"])

    return str(soup)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    docx = pathlib.Path(sys.argv[1])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    image_name = "figure-2-3-pipeline.png"
    with tempfile.TemporaryDirectory() as tmp:
        fragment = pandoc(docx, tmp)
        body = transform(fragment, image_name)
        for src, name in COPIES:
            path = pathlib.Path(src)
            if not path.is_absolute():
                path = pathlib.Path(tmp) / src
            if not path.exists():
                path = next(pathlib.Path(tmp).rglob(pathlib.Path(src).name))
            shutil.copy(path, OUT_DIR / name)
    template = (HERE / "paper_template.html").read_text()
    page = (template
            .replace("{{TITLE}}", TITLE)
            .replace("{{SUBTITLE}}", SUBTITLE)
            .replace("{{BUILT}}", datetime.date.today().isoformat())
            .replace("{{SOURCE}}", html.escape(docx.name))
            .replace("{{BODY}}", body))
    (OUT_DIR / "index.html").write_text(page)
    print(f"wrote {OUT_DIR / 'index.html'}")


if __name__ == "__main__":
    main()
