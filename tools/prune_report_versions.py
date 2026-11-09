"""Remove duplicate or incomplete Final Report versions from docs/outputs/final_report/history.

    python3 tools/prune_report_versions.py [--remove ID ...] [--note ID=TEXT ...] [--dry-run]

Besides the versions named with --remove, every "*_unrecorded" version is removed when
  * its folder is missing (the version list would link to nothing), or
  * its page is the same as a recorded version's page once the site-nav script tag is ignored
    (build_final_report.jl used to archive such copies after add_site_nav changed the live page).

For each removed version it deletes history/<id>/, drops it from the version list in history/index.html
and, in the live report.html and the archived pages, drops its entry from the "Report history" list and
marks the remaining mentions of it "since removed". --note sets a version's note (shown on the version list).
"""
import argparse, hashlib, json, os, re, shutil, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FR = os.path.join(ROOT, "docs", "outputs", "final_report")
HIST = os.path.join(FR, "history")
INDEX = os.path.join(HIST, "index.html")
NAV = re.compile(r"""<script defer src=["'][^"']*site-nav\.js["']></script>\n?""")


def read(p):
    with open(p, encoding="utf-8") as f:
        return f.read()


def write(p, s):
    with open(p, "w", encoding="utf-8") as f:
        f.write(s)


def sha(s):
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def load_index():
    s = read(INDEX)
    m = re.search(r"const H = (\{.*?\});\n", s, re.S)
    return s, m, json.loads(m.group(1))


def page(id):
    p = os.path.join(HIST, id, "report.html")
    return read(p) if os.path.isfile(p) else None


def find_duplicates(versions):
    recorded = {}
    for v in versions:
        if v["id"].endswith("_unrecorded"):
            continue
        pg = page(v["id"])
        if pg is not None:
            recorded[sha(NAV.sub("", pg))] = v["id"]
    out = {}
    for v in versions:
        id = v["id"]
        if not id.endswith("_unrecorded"):
            continue
        pg = page(id)
        if pg is None:
            out[id] = "folder missing"
        else:
            twin = recorded.get(sha(NAV.sub("", pg)))
            if twin:
                out[id] = f"copy of {twin}"
    return out


def drop_blocks(s, id):
    """Drop <details class='ver'> entries and history-table rows of version id."""
    key = f"<details class='ver'><summary><b>{id}</b>"
    while (i := s.find(key)) >= 0:
        depth, j = 0, i
        while True:
            a, b = s.find("<details", j), s.find("</details>", j)
            if a >= 0 and a < b:
                depth, j = depth + 1, a + 8
            else:
                depth, j = depth - 1, b + len("</details>")
                if depth == 0:
                    break
        s = s[:i].rstrip(" ") + s[j:]
    for href in (f"'../{id}/report.html'", f"'history/{id}/report.html'"):
        while (p := s.find(href)) >= 0:
            a, b = s.rfind("<tr", 0, p), s.find("</tr>", p)
            if a < 0 or b < 0:
                break
            s = s[:a] + s[b + 5:]
    key = f"<option value='{id}'"
    while (p := s.find(key)) >= 0:
        s = s[:p] + s[s.find("</option>", p) + len("</option>"):]
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--remove", nargs="*", default=[], metavar="ID=REASON",
                    help="version to remove, optionally with a reason (ID=REASON)")
    ap.add_argument("--note", nargs="*", default=[], metavar="ID=TEXT")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    s, m, H = load_index()
    ids = [v["id"] for v in H["versions"]]
    remove = find_duplicates(H["versions"])
    for r in a.remove:
        id, _, why = r.partition("=")
        if id not in ids:
            sys.exit(f"unknown version {id}")
        remove[id] = why or "removed"
    notes = dict(n.split("=", 1) for n in a.note)
    for id in notes:
        if id not in ids or id in remove:
            sys.exit(f"cannot set a note on {id}")

    for id, why in remove.items():
        print(f"remove {id}: {why}")
    for id, t in notes.items():
        print(f"note   {id}: {t}")
    if a.dry_run:
        return

    for id in remove:
        shutil.rmtree(os.path.join(HIST, id), ignore_errors=True)
    H["versions"] = [v for v in H["versions"] if v["id"] not in remove]
    for v in H["versions"]:
        if v["id"] in notes:
            v["note"] = notes[v["id"]]
            f = os.path.join(HIST, v["id"], "version.toml")
            if os.path.isfile(f):
                t = read(f)
                q = json.dumps(notes[v["id"]], ensure_ascii=False)
                t, n = re.subn(r'^note = ".*"$', lambda _: f"note = {q}", t, count=1, flags=re.M)
                write(f, t if n else f"note = {q}\n" + t)
    js = json.dumps(H, ensure_ascii=False, separators=(",", ":"))
    write(INDEX, s[:m.start(1)] + js + s[m.end(1):])

    pages = [os.path.join(FR, "report.html")] + [os.path.join(HIST, v["id"], "report.html") for v in H["versions"]]
    for p in pages:
        if not os.path.isfile(p):
            continue
        t0 = t = read(p)
        hit = [id for id in remove if id in t]
        for id in hit:
            t = drop_blocks(t, id)
        if hit:   # one pass, so a reason naming another removed version is not rewritten again
            t = re.sub("|".join(map(re.escape, hit)), lambda m: f"{m.group(0)}, since removed", t)
        if t != t0:
            write(p, t)
            print(f"updated {os.path.relpath(p, ROOT)}")


if __name__ == "__main__":
    main()
