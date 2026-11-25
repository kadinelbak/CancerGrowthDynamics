#!/usr/bin/env python3
"""Build staged.js for the staged presentation (staged.html).

Run after a pipeline run, from the repository root:
    python3 docs/outputs/results/build_staged.py
It gathers what staged.html cannot read from bundle.js: the profile-likelihood, local and
global (Sobol) sensitivity tables, the bootstrap selection shares, the controller's equations
and plant parameters (from the 120-day control report) and the optimal-control figures of the
headline scenario (SVGs copied out of the 14-day and 120-day control reports) plus their outcome
rows. Everything else (equations, fitted values, bootstrap ranges, data, fitted curves) comes
from bundle.js at page load, so the deck always matches the latest run.
"""
import csv, json, re, html, pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = ROOT / "Modeling_Approaches/02_mechanical_automatic_package/outputs"
HERE = pathlib.Path(__file__).resolve().parent
import tomllib
with open(ROOT / "Modeling_Approaches/02_mechanical_automatic_package/control_config.toml", "rb") as _f:
    HEADLINE = tomllib.load(_f).get("scenarios", {}).get("headline", "A")   # control_config.toml [scenarios] headline
RATIO = "50:50"

def rows(p):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def num(x):
    try:
        v = float(x)
        return None if v != v else v
    except (TypeError, ValueError):
        return None

# ── validation tables ──────────────────────────────────────────────────────────
prof = [{k: (num(r[k]) if k in ("estimate", "ci_lo", "ci_hi", "threshold", "grid_lo", "grid_hi", "bound_position") else r[k])
         for k in ("group", "param", "estimate", "ci_lo", "ci_hi", "status", "threshold", "grid_lo", "grid_hi", "at_bound", "near_bound", "bound_position")}
        for r in rows(OUT / "zoo_analyses/sensitivity_profile.csv")]
local = [{"group": r["group"], "param": r["param"], "value": num(r["value"]), "si": num(r["si"]),
          "max_rel_change": num(r["max_rel_change"])} for r in rows(OUT / "zoo_analyses/sensitivity_local.csv")]
glob = [{"scope": r["scope"], "output": r["output"], "param": r["param"], "S1": num(r["S1"]), "ST": num(r["ST"]),
         "S1_lo": num(r["S1_lo"]), "S1_hi": num(r["S1_hi"]), "ST_lo": num(r["ST_lo"]), "ST_hi": num(r["ST_hi"]), "n": int(r["n"])}
        for r in rows(OUT / "zoo_analyses/sensitivity_global.csv")]
boot = [{"group": r["group"], "candidate": r["candidate"], "k": int(r["k"]), "pct_chosen": num(r["pct_chosen"]),
         "pct_lowest_bic": num(r["pct_lowest_bic"]), "optimal": r["optimal"] == "true"}
        for r in rows(OUT / "zoo_analyses/bootstrap_selection.csv")]
bootp = [{"group": r["group"], "model": r["model"], "param": r["param"], "estimate": num(r["estimate"]),
          "lo": num(r["block_lo"]), "hi": num(r["block_hi"])} for r in rows(OUT / "zoo_analyses/bootstrap_parameters.csv")]

# ── control reports: equations, plant parameters, figures ──────────────────────
def section(s, start, end):
    a = s.index(start); b = s.index(end, a)
    return s[a:b]

def svgs_after(s, heading, n):
    """The first n SVG figures after the heading text (the report's figures for one scenario)."""
    a = s.index(heading)
    out = []
    for m in re.finditer(r"<svg.*?</svg>", s[a:], re.S):
        svg = m.group(0)
        if len(svg) < 1000:
            continue
        out.append(svg)
        if len(out) == n:
            break
    return out

def svg_title(svg):
    m = re.search(r"<text[^>]*font-weight='bold'[^>]*>([^<]*)</text>", svg)
    return html.unescape(m.group(1)) if m else ""

r120 = (OUT / "control_120d/report.html").read_text(encoding="utf-8")
r14 = (OUT / "control/report.html").read_text(encoding="utf-8")

eqsec = section(r120, "<h2>The model the controller uses", "<h2>Checks against")
intro = re.sub(r"<[^>]+>", "", section(eqsec, "<p>", "<h3>")).strip()
blocks = []   # [heading, [tex...]]
for part in re.split(r"<h3>", eqsec)[1:]:
    head = re.sub(r"<[^>]+>", "", part.split("</h3>")[0]).strip()
    tex = [html.unescape(t) for t in re.findall(r"\\\[(.*?)\\\]", part, re.S)]
    blocks.append([head, tex])
plant = []
for tr in re.findall(r"<tr>(.*?)</tr>", section(eqsec, "<table>", "</table>"), re.S):
    cells = [html.unescape(re.sub(r"<[^>]+>", "", c)).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)]
    if len(cells) == 5:
        plant.append(cells)

# objective, penalties and barrier text of the 14-day report (the dose-only controller)
obj14 = section(r14, "<h2>Objective, soft penalties and barriers", "<h2>Results")
obj14_tex = [html.unescape(t) for t in re.findall(r"\\\[(.*?)\\\]", obj14, re.S)]

scen_head14 = f"<h3>Scenario {HEADLINE}:"
scen_head120 = f"<h3>Scenario {HEADLINE}:"
f14 = svgs_after(r14, scen_head14, 5)      # total, final-vs-dose, N and C, dose schedules, toxicity
f120 = svgs_after(r120, scen_head120, 7)   # total, final-vs-dose, N and C, schedules, toxicity, weight, protection
OTHER = "A" if HEADLINE == "B" else "B"
f120o = svgs_after(r120, f"<h3>Scenario {OTHER}:", 7)
figs14 = {svg_title(v): v for v in f14}
figs120 = {svg_title(v): v for v in f120}
figs120o = {svg_title(v): v for v in f120o}

# ── outcome rows of the headline scenario, 50:50 ───────────────────────────────
keep14 = ("policy", "schedule_uM", "N_final", "C_final", "total_final", "resistant_fraction", "max_toxicity",
          "cumulative_dose_uM_days", "first_progression_day", "time_to_progression_day", "control_stopped_day")
out14 = [{k: (num(r[k]) if k not in ("policy", "schedule_uM") else r[k]) for k in keep14}
         for r in rows(OUT / "control/control_outcomes.csv")
         if r["environment"] == "normal" and r["scenario"] == HEADLINE and r["ratio"] == RATIO]
keep120 = ("policy", "N_final", "C_final", "total_final", "resistant_fraction", "max_toxicity", "min_weight",
           "cumulative_dose_uM_days", "fasting_days", "first_progression_day", "time_to_progression_day", "control_stopped_day")
out120 = [{k: (num(r[k]) if k != "policy" else r[k]) for k in keep120}
          for r in rows(OUT / "control_120d/fasting_outcomes.csv")
          if r["scenario"] == HEADLINE and r["ratio"] == RATIO]
out120o = [{k: (num(r[k]) if k != "policy" else r[k]) for k in keep120}
           for r in rows(OUT / "control_120d/fasting_outcomes.csv")
           if r["scenario"] == OTHER and r["ratio"] == RATIO]

# ── control settings worth stating on the slides (control_config.toml) ─────────
cfg = (ROOT / "Modeling_Approaches/02_mechanical_automatic_package/control_config.toml").read_text(encoding="utf-8")
def cfgv(key):
    m = re.search(rf"^{key}\s*=\s*([^#\n]+)", cfg, re.M)
    return m.group(1).strip() if m else None
settings = {k: cfgv(k) for k in ("decision_interval_days", "lookahead_blocks", "lookahead_to_end", "u_max_data", "u_max_demo",
                                 "threshold_fraction_of_K", "step_width_fraction_of_K", "qP", "qRT", "substeps", "tiebreak_scale",
                                 "confirm_decisions", "settle_margin_days", "sustain_days", "qB", "qR", "qT", "ru", "terminal_factor",
                                 "mu_cancer", "mu_tox", "rv", "qW", "mu_host", "headline")}

# the 120-day report's barrier line is cut at a raw "<" in its HTML; take the complete definition from the 14-day report
bmu = next((t[t.index("b_\\mu(z)="):t.index("\\end{cases}") + len("\\end{cases}")] for t in obj14_tex if "\\begin{cases}" in t and "\\end{cases}" in t), None)
for b in blocks:
    for i, t in enumerate(b[1]):
        if "\\begin{cases}" in t and "\\end{cases}" not in t and bmu:
            b[1][i] = t[:t.index("b_\\mu(z)=")] + bmu

generated = re.search(r'"generated":"([^"]+)"', (HERE / "bundle.js").read_text(encoding="utf-8")).group(1)
data = {"generated": generated, "headline": HEADLINE, "ratio": RATIO, "profile": prof, "local": local, "global": glob,
        "bootSelection": boot, "bootParams": bootp, "ctrlIntro": intro, "ctrlBlocks": blocks, "plant": plant,
        "obj14": obj14_tex, "figs14": figs14, "figs120": figs120, "figs120o": figs120o, "other": OTHER, "out14": out14, "out120": out120, "out120o": out120o, "settings": settings}
js = "// written by build_staged.py — data for the staged presentation (staged.html)\nwindow.STAGED = " + json.dumps(data, ensure_ascii=False) + ";\n"
(HERE / "staged.js").write_text(js, encoding="utf-8")
print(f"staged.js: {len(js)/1e6:.2f} MB; figures 14 d: {list(figs14)}; 120 d: {list(figs120)}; outcomes {len(out14)}/{len(out120)}; blocks {[b[0][:30] for b in blocks]}; plant rows {len(plant)}")
