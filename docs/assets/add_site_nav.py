"""Give the main pages of the site the same top bar (docs/assets/site-nav.js).

Run after regenerating reports:  python3 docs/assets/add_site_nav.py
Safe to run again: a page that already loads the script is
left alone. Only the main pages get the bar; archives, the neural-network track and the old-data reports keep
their own navigation.
"""
import os
import re

DOCS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MARK = "site-nav.js"
PKG = "Modeling_Approaches/02_mechanical_automatic_package/outputs/"
PAGES = ["index.html", "papers/index.html", "papers/from-culture-to-control/index.html",
         "education/index.html", "education/pipeline-map.html", "education/pipeline-equations.html",
         "outputs/results/summary.html", "outputs/results/presentation.html", "outputs/results/staged.html", "outputs/results/explorer.html",
         "outputs/final_report/report.html", "outputs/final_report/history/index.html", "outputs/time_keeping/index.html",
         PKG + "control/report.html", PKG + "control_30d/report.html", PKG + "control_60d/report.html",
         PKG + "control_120d/report.html", PKG + "control_sweep/report.html", PKG + "findings/report.html",
         PKG + "zoo_fits/report.html", PKG + "zoo_analyses/bootstrap/report.html",
         PKG + "zoo_analyses/sensitivity/report.html", PKG + "zoo_analyses/validation/report.html",
         PKG + "zoo_analyses/resistant_kill/report.html"]


def add(page):
    with open(page, encoding="utf-8") as f:
        html = f.read()
    if MARK in html[:20000]:
        return False
    rel = os.path.relpath(os.path.join(DOCS, "assets", MARK), os.path.dirname(page)).replace(os.sep, "/")
    tag = f'<script defer src="{rel}"></script>'
    m = re.search(r"</head\s*>", html, re.I)   # the first </head> is the page's own
    html = html[: m.start()] + tag + "\n" + html[m.start():] if m else tag + "\n" + html
    with open(page, "w", encoding="utf-8") as f:
        f.write(html)
    return True


def main():
    n = sum(add(os.path.join(DOCS, p)) for p in PAGES if os.path.isfile(os.path.join(DOCS, p)))
    print(f"added the site navigation bar to {n} page(s)")


if __name__ == "__main__":
    main()
