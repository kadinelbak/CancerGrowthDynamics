"""Add the research-quality rows to Final Report versions that were archived before they existed.

One-off (2026-11-01). report_history.jl now records research-quality metrics with every new version
(research_quality: held-out error, parameters pinned, selection stability, control robustness, scenario
agreement). This computes the same rows, with the same keys and value formats, from the output CSVs at
each archived version's commit, adds the ones a version.toml lacks (nothing else in the file changes), and
rewrites history/index.html.   python tools/backfill_quality.py

Which commit's outputs belong to a version: a version written by build_final_report.jl records the HEAD it
was built on, but its outputs were committed afterwards, together with the archived version, so the commit
that added the version's folder is tried first, then the recorded commit. A commit is used only if its
control outputs reproduce the version's recorded 14-day MPC totals; a version no commit matches (e.g. one
built between a refit and the control rerun that was committed with it) gets no quality rows. A commit that
lacks an output gives no row for that metric.
"""
import csv, io, json, math, statistics, subprocess, sys, tomllib
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PKG = "Modeling_Approaches/02_mechanical_automatic_package"
HIST = ROOT / "docs/outputs/final_report/history"


def git_csv(commit, rel):
    r = subprocess.run(["git", "-C", str(ROOT), "show", f"{commit}:{PKG}/outputs/{rel}"], capture_output=True)
    if r.returncode != 0:
        return []
    return list(csv.DictReader(io.StringIO(r.stdout.decode("utf-8"))))


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return math.nan


def pct(a, b):
    return "n/a" if b == 0 else f"{round(100 * a / b)}% ({a} of {b})"


def jround(x, d):
    """Julia's round(x; digits=d) printed the way Julia prints a Float64."""
    v = round(x, d)
    return repr(float(v))


def mpc_beats(m, b, run_days):
    def ttp(r):
        v = num(r.get("time_to_progression_day", ""))
        return v if math.isfinite(v) else run_days + 1
    tm, tb = ttp(m), ttp(b)
    return tm > tb + 1e-6 or (abs(tm - tb) <= 1e-6 and num(m["total_final"]) < num(b["total_final"]))


def mpc_wins(rows, keys, run_days_of):
    groups = {}
    for r in rows:
        groups.setdefault(tuple(r.get(k, "") for k in keys), []).append(r)
    won = n = 0
    for g in groups.values():
        pol = lambda name: next((r for r in g if r.get("policy") == name), None)
        m, a, c = pol("MPC, measured doses only"), pol("No treatment"), pol("Constant 1 µM")
        if m is None or a is None or c is None:
            continue
        n += 1
        rd = run_days_of(m)
        won += mpc_beats(m, a, rd) and mpc_beats(m, c, rd)
    return won, n


def quality(commit):
    out = []
    v = [num(r["nrmse_heldout"]) for r in git_csv(commit, "zoo_analyses/validation_blocked.csv") if r.get("status") == "completed"]
    v = [x for x in v if math.isfinite(x)]
    if v:
        out.append(("q-:heldout", "Quality · held-out error, median NRMSE over blocked checks (lower is better)",
                    f"{jround(statistics.median(v), 3)} (n = {len(v)}, worst {jround(max(v), 3)})"))
    ce = git_csv(commit, "zoo_fits/heldout_ce1.csv")
    r = next((r for r in ce if r.get("row", "").startswith("prediction")), None)
    if r:
        out.append(("q-:heldout_ce1", "Quality · held-out coculture run (Ce1) NRMSE (lower is better)", jround(num(r["nrmse"]), 3)))
    pr = git_csv(commit, "zoo_analyses/sensitivity_profile.csv")
    if pr:
        ok = sum(r.get("status") == "identifiable" for r in pr)
        miss = [f"{r['param']} ({r['group']})" for r in pr if r.get("status") != "identifiable"]
        out.append(("q+:pinned", "Quality · parameters pinned by their profile (higher is better)",
                    pct(ok, len(pr)) + ("; not pinned: " + ", ".join(miss) if miss else "")))
    bs = [r for r in git_csv(commit, "zoo_analyses/bootstrap_selection.csv") if r.get("optimal") == "true"]
    if bs:
        p = [num(r["pct_chosen"]) for r in bs]
        j = p.index(min(p))
        out.append(("q+:stability", "Quality · selection stability, mean % of bootstrap refits that pick the carried model (higher is better)",
                    f"{round(sum(p) / len(p))}% (lowest {round(p[j])}%, {bs[j]['group']})"))
    cc = [r for r in git_csv(commit, "control/control_outcomes.csv") if r.get("environment") == "normal"]
    if cc:
        w, n = mpc_wins(cc, ["scenario", "ratio"], lambda _: 14.0)
        if n:
            out.append(("q+:robust", "Quality · 14-day MPC beats no treatment and constant 1 µM, scenarios × ratios (higher is better)", pct(w, n)))
    sw = [r for r in git_csv(commit, "control_sweep/sweep_outcomes.csv") if r.get("environment", "normal") == "normal"]
    if sw:
        def rd(r):
            d = num(r.get("run_days", "14"))
            return d if math.isfinite(d) else 14.0
        w, n = mpc_wins(sw, ["data", "variant", "scenario", "ratio"], rd)
        if n:
            out.append(("q+:robust_sweep", "Quality · MPC beats both baselines across the control sweep (higher is better)", pct(w, n)))
    if cc:
        same = tot = 0
        for ratio in dict.fromkeys(r.get("ratio", "") for r in cc):
            def f(sc):
                r = next((r for r in cc if r["scenario"] == sc and r["ratio"] == ratio and r["policy"] == "MPC, measured doses only"), None)
                return None if r is None else r["schedule_uM"].split()
            a, b = f("A"), f("B")
            if a is None or b is None or len(a) != len(b):
                continue
            same += sum(num(x) == num(y) for x, y in zip(a, b))
            tot += len(a)
        if tot:
            out.append(("q+:agree", "Quality · scenario agreement, 14-day MPC dose blocks the same in A and B (higher is better)", pct(same, tot)))
    return out


def git_out(*a):
    r = subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True)
    return r.stdout.decode("utf-8").strip() if r.returncode == 0 else ""


def matching_commit(v):
    """the commit whose outputs are this version's (see the module docstring), or None"""
    added = git_out("log", "--diff-filter=A", "--format=%h", "--", f"docs/outputs/final_report/history/{v['id']}/version.toml").splitlines()
    cands = [c for c in ([added[-1]] if added else []) + [v.get("commit", "")] if c]
    want = {r["key"].split(":")[1]: r["value"] for r in v.get("results", []) if r["key"].startswith("mpc_total:")}
    for c in cands:
        if not want:
            return c
        cc = git_csv(c, "control/control_outcomes.csv")
        got = {r["scenario"]: str(round(num(r["total_final"]))) for r in cc
               if r.get("environment") == "normal" and r.get("ratio") == "50:50" and r.get("policy") == "MPC, measured doses only"}
        if all(got.get(k) == x for k, x in want.items()):
            return c
    return None


def tq(s):
    return json.dumps(s, ensure_ascii=False)


def main():
    for f in sorted(HIST.glob("*/version.toml")):
        text = f.read_text(encoding="utf-8")
        v = tomllib.loads(text)
        have = {r["key"] for r in v.get("results", [])}
        c = matching_commit(v)
        rows = [r for r in quality(c) if r[0] not in have] if c else []
        if rows:
            block = "".join(f"\n[[results]]\nkey = {tq(k)}\nlabel = {tq(l)}\nvalue = {tq(x)}\n" for k, l, x in rows)
            i = text.find("\n[[chapters]]")
            text = text + block if i < 0 else text[:i] + "\n" + block.lstrip("\n").rstrip("\n") + "\n" + text[i:]
            f.write_text(text, encoding="utf-8")
            tomllib.loads(text)  # still valid
        print(v["id"], "outputs from", c or "no matching commit", f"+{len(rows)} quality rows")
    versions = sorted((tomllib.loads(f.read_text(encoding="utf-8")) for f in HIST.glob("*/version.toml")), key=lambda v: v["id"])
    dec = tomllib.loads((ROOT / PKG / "decisions.toml").read_text(encoding="utf-8")).get("decision", [])
    data = dict(versions=versions, decisions=dec, generated=datetime.now().strftime("%Y-%m-%d %H:%M"))
    tpl = (ROOT / "report_history_template.html").read_text(encoding="utf-8")
    js = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    (HIST / "index.html").write_text(tpl.replace("/*HISTORY_JSON*/null", js), encoding="utf-8")
    print("wrote", HIST / "index.html", len(versions), "versions")


if __name__ == "__main__":
    sys.exit(main())
