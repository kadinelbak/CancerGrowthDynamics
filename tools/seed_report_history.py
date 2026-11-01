"""Seed docs/outputs/final_report/history/ with the Final Report versions that predate report_history.jl.

One-off (2026-10-27): every committed version of docs/outputs/final_report/report.html is taken from git,
with the stages.toml / control_config.toml settings and the output CSVs of the same commit, and archived
exactly as build_final_report.jl archives a new build (same keys, value formats and link rebasing as
report_history.jl). Then history/index.html is written.   python tools/seed_report_history.py
"""
import csv, hashlib, io, json, re, subprocess, sys, tomllib
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PKG = "Modeling_Approaches/02_mechanical_automatic_package"
HIST = ROOT / "docs/outputs/final_report/history"
REPORT = "docs/outputs/final_report/report.html"


def git(*a, binary=False):
    r = subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True)
    if r.returncode != 0:
        return None
    return r.stdout if binary else r.stdout.decode("utf-8")


def fmt(x):
    if x is None:
        return "not set"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, (list, tuple)):
        return ", ".join(fmt(v) for v in x)
    if isinstance(x, float):
        return str(int(x)) if x.is_integer() else repr(x)
    return str(x)


def settings(c):
    st = git("show", f"{c}:{PKG}/stages.toml"); cc = git("show", f"{c}:{PKG}/control_config.toml")
    S = tomllib.loads(st) if st else {}; C = tomllib.loads(cc) if cc else {}
    stg = lambda i: next((s for s in S.get("stage", []) if s.get("id") == i), {})
    sel = S.get("selection", {}); sch = C.get("schedule", {}); can = C.get("cancer", {})
    prog = C.get("progression", {}); w = C.get("weights", {})
    rule = sel.get("rule_bic")
    if rule is None:   # before 2026-10-27 the rule lived in zoo_models.jl
        m = re.search(r"const RULE_BIC = ([0-9.]+)", git("show", f"{c}:{PKG}/zoo_models.jl") or "")
        rule = float(m.group(1)) if m else None
    wt = sel.get("residual_weighting", "none")
    ex = stg("s2").get("exclude", [])
    lines = S.get("lines", {})
    d = float(sch.get("decision_interval_days", 2.0))
    look = (f"{fmt(sum(sch['lookahead_blocks']) * d)} days (blocks {fmt(sch['lookahead_blocks'])})" if "lookahead_blocks" in sch
            else f"{fmt(sch['horizon_days'])} days" if "horizon_days" in sch else None)
    days = []
    for k in ("closed_loop_days", "run_days", "long_days"):
        v = sch.get(k)
        days += [float(x) for x in (v if isinstance(v, list) else [v] if v is not None else [])]
    days = sorted(set(days))
    bc = can.get("B_cancer_fraction_of_K")
    obj = ", ".join(f"{k} = {fmt(v)}" for k, v in (("qP", prog.get("qP", 0.0)), ("qRT", prog.get("qRT", 0.0)), ("qB", w.get("qB")), ("qR", w.get("qR")),
                                                   ("qT", w.get("qT")), ("ru", w.get("ru")), ("terminal factor", w.get("terminal_factor"))) if v is not None)
    return [
        ("rule_bic", "Model selection rule", "not recorded" if rule is None else f"fewest parameters within {fmt(rule)} BIC"),
        ("weighting", "Residual weighting", f"power: max(y, 50)^{fmt(sel.get('weight_power', 0.55))}" if wt == "power" else fmt(wt)),
        ("s2_exclude", "Stage 2b series left out", "none" if not ex else "; ".join(f"{e.get('line', '?')} {fmt(e.get('dose', '?'))} µM" for e in ex)),
        ("kill_coupled", "Kill coupled to growth (Norton–Simon)", fmt(stg("s2").get("kill_growth_coupled", False))),
        ("betaC_floor", "Resistant kill floor in coculture (× Stage 2b)", fmt(stg("s3").get("betaC_floor"))),
        ("ec50", "EC50, A2780 / A2780cis (µM)", "not recorded" if not lines else f"{fmt(lines.get('A2780', {}).get('ec50'))} / {fmt(lines.get('A2780cis', {}).get('ec50'))}"),
        ("B_cancer", "B_cancer (burden penalty level)", "not set" if bc is None else f"{fmt(bc)} K" if bc > 0 else "none"),
        ("B_prog", "B_prog (progression threshold)", f"{fmt(prog['threshold_fraction_of_K'])} K" if "threshold_fraction_of_K" in prog else "no progression goal"),
        ("barrier", "Burden barrier shape", fmt(can.get("barrier", "relaxed_log"))),
        ("lookahead", "Look-ahead", fmt(look)),
        ("to_end", "Look-ahead ends with the run", fmt(sch.get("lookahead_to_end", False))),
        ("run_days", "Control run lengths (days)", "not set" if not days else fmt(days)),
        ("objective", "Objective weights", obj),
        ("B_tox", "Toxicity limit B_tox", fmt(C.get("toxicity", {}).get("B_tox"))),
        ("doses", "Measured doses (µM)", fmt(C.get("dose", {}).get("levels_data"))),
    ]


def rows(c, path):
    t = git("show", f"{c}:{path}")
    return list(csv.DictReader(io.StringIO(t))) if t else []


def results(c):
    out = []
    for r in rows(c, f"{PKG}/outputs/zoo_fits/zoo_model_selection.csv"):
        if r.get("optimal") == "true":
            out.append(("carried:" + r["group"], "Carried model · " + r["group"], r["candidate"]))
    cc = rows(c, f"{PKG}/outputs/control/control_outcomes.csv")
    for scen in ("A", "B"):
        r = next((r for r in cc if r.get("environment") == "normal" and r.get("scenario") == scen and r.get("ratio") == "50:50"
                  and r.get("policy") == "MPC, measured doses only"), None)
        if r is None:
            continue
        out.append((f"mpc:{scen}", f"14-day MPC schedule · scenario {scen}, 50:50 (µM)", r["schedule_uM"]))
        out.append((f"mpc_total:{scen}", f"14-day MPC N + C at the end · scenario {scen}, 50:50", str(round(float(r["total_final"])))))
        if "time_to_progression_day" in r:
            try:
                v = float(r["time_to_progression_day"])
            except ValueError:
                v = float("nan")
            out.append((f"mpc_ttp:{scen}", f"14-day MPC time to progression · scenario {scen}, 50:50",
                        "not progressed" if v != v else f"day {fmt(round(v, 1))}"))
    return out


def rebase(page):   # as report_history.jl rebase_archived
    i = page.find("const PARTS")
    head, rest = (page, "") if i < 0 else (page[:i], page[i:])
    for a, b in (("href='../../index.html'", "href='../../../../index.html'"), ("href='../archive", "href='../../../archive"),
                 ("href='history/", "href='../"), ("href='report.html'", "href='../../report.html'")):
        head = head.replace(a, b)
    rest = rest.replace('<base href=\\"../../', '<base href=\\"../../../../').replace('src:"../../', 'src:"../../../../')
    return head + rest


def tq(s):
    return json.dumps(s, ensure_ascii=False)


def write_toml(rec, f):
    o = [f"{k} = {tq(rec[k])}" for k in ("id", "date", "commit", "source", "note", "content_sha", "report_sha")]
    for grp in ("settings", "results"):
        for k, l, v in rec[grp]:
            o += ["", f"[[{grp}]]", f"key = {tq(k)}", f"label = {tq(l)}", f"value = {tq(v)}"]
    for t, s in rec["chapters"]:
        o += ["", "[[chapters]]", f"title = {tq(t)}", f"stamp = {tq(s)}"]
    f.write_text("\n".join(o) + "\n", encoding="utf-8")


def main():
    log = git("log", "--format=%h|%ad|%s", "--date=format:%Y-%m-%d %H:%M", "--", REPORT).strip().splitlines()
    for line in reversed(log):
        c, date, subj = line.split("|", 2)
        raw = git("show", f"{c}:{REPORT}", binary=True)
        page = raw.decode("utf-8")
        vid = datetime.strptime(date, "%Y-%m-%d %H:%M").strftime("%Y-%m-%d_%H%M")
        chapters = [(json.loads(t), json.loads(s)) for t, s in re.findall(r'title:("(?:[^"\\]|\\.)*"),blurb:.*?,stamp:("[^"]*")', page)]
        d = HIST / vid; d.mkdir(parents=True, exist_ok=True)
        (d / "report.html").write_text(rebase(page), encoding="utf-8")
        rec = dict(id=vid, date=date, commit=c, source="git history (seeded 2026-10-27 by tools/seed_report_history.py)", note=subj,
                   content_sha="", report_sha=hashlib.sha256(raw).hexdigest(), settings=settings(c), results=results(c), chapters=chapters)
        write_toml(rec, d / "version.toml")
        print(vid, c, len(chapters), "chapters", len(rec["settings"]), "settings", len(rec["results"]), "results")
    versions = []
    for f in sorted(HIST.glob("*/version.toml")):
        v = tomllib.loads(f.read_text(encoding="utf-8"))
        versions.append(v)
    versions.sort(key=lambda v: v["id"])
    dec = tomllib.loads((ROOT / PKG / "decisions.toml").read_text(encoding="utf-8")).get("decision", [])
    data = dict(versions=versions, decisions=dec, generated=datetime.now().strftime("%Y-%m-%d %H:%M"))
    tpl = (ROOT / "report_history_template.html").read_text(encoding="utf-8")
    js = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    (HIST / "index.html").write_text(tpl.replace("/*HISTORY_JSON*/null", js), encoding="utf-8")
    print("wrote", HIST / "index.html", len(versions), "versions")


if __name__ == "__main__":
    sys.exit(main())
