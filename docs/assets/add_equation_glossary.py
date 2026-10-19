"""Add the click-to-see-variables equation script to every report page under docs/.

Run after regenerating reports:  python3 docs/assets/add_equation_glossary.py
Safe to run again: pages that already load the script are left alone. Pages are found by
walking docs/, and the script tag uses a path relative to each page so it works on GitHub
Pages and anywhere else the folder is hosted.
"""
import os
import re

DOCS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MARK = "equation-glossary.js"


def tag_for(page):
    rel = os.path.relpath(os.path.join(DOCS, "assets", MARK), os.path.dirname(page)).replace(os.sep, "/")
    return f'<script defer src="{rel}"></script>'


def add(page):
    with open(page, encoding="utf-8") as f:
        html = f.read()
    if MARK in html:
        return False
    tag = tag_for(page)
    # the first </head> is the page's own; later ones can sit inside embedded report strings
    m = re.search(r"</head\s*>", html, re.I)
    if m:
        html = html[: m.start()] + tag + "\n" + html[m.start():]
    else:
        m = list(re.finditer(r"</body\s*>", html, re.I))
        html = html[: m[-1].start()] + tag + "\n" + html[m[-1].start():] if m else html + "\n" + tag + "\n"
    with open(page, "w", encoding="utf-8") as f:
        f.write(html)
    return True


def main():
    changed = 0
    for d, _, files in os.walk(DOCS):
        for name in files:
            if name.endswith(".html") and add(os.path.join(d, name)):
                changed += 1
    print(f"added the equation glossary script to {changed} page(s)")


if __name__ == "__main__":
    main()
