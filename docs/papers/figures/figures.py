"""Results figures for the thesis (Section 3), from the committed data and run outputs.

Usage: python3 docs/papers/figures/figures.py OUTDIR   (needs numpy and matplotlib)
The models are Python copies of the selected zoo models at their point estimates; the 14-day
control curves reproduce control_outcomes.csv to the cell.
"""
import csv
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from models import *

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
PKGO = f"{PKG}/outputs"

plt.rcParams.update({
    "font.family": "Inter", "font.size": 9, "axes.titlesize": 9.5, "axes.titleweight": "bold",
    "axes.labelsize": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#8a8a85", "axes.linewidth": 0.8, "xtick.color": "#52606d", "ytick.color": "#52606d",
    "axes.labelcolor": "#1f2933", "axes.grid": True, "grid.color": "#e6e6e3", "grid.linewidth": 0.6,
    "legend.frameon": False, "legend.fontsize": 8, "lines.linewidth": 1.8, "savefig.dpi": 300,
})
SENS, RES = "#2a78d6", "#eb6834"          # A2780 (sensitive), A2780cis (resistant)
AQUA, YELLOW, MAG, GRAY = "#1baf7a", "#eda100", "#e87ba4", "#8a8a85"
BLUES = ["#9ec5f4", "#5598e7", "#256abf", "#0d366b"]   # dose 0, 0.62, 1, 1.47
ORANGES = ["#f6c2a8", "#ef8a5e", "#c9501f", "#7f3010"]
DOSES = [("const_0", "0"), ("const_0_62", "0.62 µM"), ("const_1", "1 µM"), ("const_1_47", "1.47 µM")]

rows = load_measurements()
st = starts()


def pts(ax, d, v, sem, ex, color, label=None, hollow=False):
    m = d > 0
    ax.errorbar(d[m & ~ex], v[m & ~ex], yerr=sem[m & ~ex], fmt="o", ms=3.6, color=color, mfc="white" if hollow else color,
                mew=1.1, elinewidth=0.8, capsize=0, label=label, zorder=3)
    if (m & ex).any():
        ax.plot(d[m & ex], v[m & ex], "x", color=color, ms=4, zorder=3)


def label_panel(ax, s):
    ax.set_title(s, loc="left")


def finish(fig, name):
    fig.savefig(f"{OUT}/{name}.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def cells_axis(ax, ymax=None, xlabel=True):
    if xlabel:
        ax.set_xlabel("Day")
    ax.set_xlim(0, 14.5)
    ax.set_xticks([0, 2, 4, 6, 8, 10, 12, 14])
    if ymax:
        ax.set_ylim(0, ymax)


# ---------------------------------------------------------------- Fig 3.1 data
fig, axs = plt.subplots(2, 4, figsize=(10, 5.2), sharex=True)
ax = axs[0, 0]
for line, col in (("A2780", SENS), ("A2780cis", RES)):
    d, v, s, ex = series(rows, condition="untreated", line=line, culture="mono")
    pts(ax, d, v, s, ex, col, line)
    ax.plot(d[d > 0], v[d > 0], color=col, lw=1, alpha=.5)
label_panel(ax, "a  Untreated monoculture"); ax.legend(loc="upper left")
for ax, line, ramp, title in ((axs[0, 1], "A2780", BLUES, "b  A2780 + cisplatin"),
                              (axs[0, 2], "A2780cis", ORANGES, "c  A2780cis + cisplatin")):
    for (sch, lab), col in zip(DOSES, ramp):
        cond = "untreated" if sch == "const_0" else "treated"
        d, v, s, ex = series(rows, condition=cond, line=line, culture="mono", schedule=sch)
        pts(ax, d, v, s, ex, col, lab); ax.plot(d[d > 0], v[d > 0], color=col, lw=1, alpha=.5)
    label_panel(ax, title); ax.legend(loc="upper left", title="Dose", title_fontsize=8)
ax = axs[0, 3]
for line, col in (("A2780", SENS), ("A2780cis", RES)):
    d, v, s, ex = series(rows, condition="low_resource", line=line, culture="mono")
    pts(ax, d, v, s, ex, col, line); ax.plot(d[d > 0], v[d > 0], color=col, lw=1, alpha=.5)
label_panel(ax, "d  Low resource (0.5% FBS)"); ax.legend(loc="upper left")
RATIOS = [("75S25R", "75:25"), ("50S50R", "50:50"), ("25S75R", "25:75")]
STYLES = ["-", "--", ":"]
for ax, cond, title, extra in ((axs[1, 0], "Ce0", "e  Coculture, untreated", {}),
                               (axs[1, 1], "Ce1", "f  Coculture, 1 µM", {}),
                               (axs[1, 2], "low_resource", "g  Coculture, low resource (run 1)", {"run": "1"})):
    for (rk, rl), ls in zip(RATIOS, STYLES):
        for line, col in (("A2780", SENS), ("A2780cis", RES)):
            d, v, s, ex = series(rows, condition=cond, line=line, culture="co", ratio=rk, **extra)
            m = (d > 0) & ~ex
            ax.plot(d[m], v[m], ls, color=col, lw=1.4, marker="o", ms=2.6)
    label_panel(ax, title)
ax = axs[1, 3]
for line, col in (("A2780", SENS), ("A2780cis", RES)):
    d, v, s, ex = series(rows, condition="combination", line=line, culture="co")
    pts(ax, d, v, s, ex, col, line); ax.plot(d[d > 0], v[d > 0], color=col, lw=1, alpha=.5)
label_panel(ax, "h  Low resource + 1 µM, 50:50"); ax.legend(loc="upper left")
from matplotlib.lines import Line2D
axs[1, 0].legend(handles=[Line2D([], [], color=SENS, label="A2780"), Line2D([], [], color=RES, label="A2780cis")] +
                 [Line2D([], [], color=GRAY, ls=ls, label=f"S:R {rl}") for (rk, rl), ls in zip(RATIOS, STYLES)],
                 loc="upper left", fontsize=7)
for ax in axs.flat:
    cells_axis(ax, 4800, xlabel=ax in axs[1])
for ax in axs[:, 0]:
    ax.set_ylabel("Cells per image")
fig.tight_layout()
finish(fig, "fig3_1_data")

# ---------------------------------------------------------------- Fig 3.2 Stage 1 + 2a
fig, axs = plt.subplots(1, 4, figsize=(10, 2.8), sharey=True)
specs = [("a  A2780 · logistic", dict(condition="untreated", line="A2780", culture="mono"), logistic, [100], SENS),
         ("b  A2780cis · theta-logistic", dict(condition="untreated", line="A2780cis", culture="mono"), theta, [100], RES),
         ("c  A2780 · low resource", dict(condition="low_resource", line="A2780", culture="mono"), lr_N, [st["mono_lowres_A2780"], 0, 0], SENS),
         ("d  A2780cis · low resource", dict(condition="low_resource", line="A2780cis", culture="mono"), lr_C, [st["mono_lowres_A2780cis"], 0], RES)]
for ax, (title, sel, f, y0, col) in zip(axs, specs):
    d, v, s, ex = series(rows, **sel)
    pts(ax, d, v, s, ex, col, "Data (mean ± SEM)")
    t, y = rk4(f, y0, 14, const(0))
    ax.plot(t, y[:, 0], color=col, label="Selected model")
    label_panel(ax, title); cells_axis(ax, 4800)
axs[0].set_ylabel("Cells per image"); axs[0].legend(loc="lower right", fontsize=7)
fig.tight_layout(); finish(fig, "fig3_2_stage1_2a")

# ---------------------------------------------------------------- Fig 3.3 Stage 2b
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for ax, line, ramp, f, title in ((axs[0], "A2780", BLUES, drug_N, "a  Stage 2b · A2780"),
                                 (axs[1], "A2780cis", ORANGES, drug_C, "b  Stage 2b · A2780cis")):
    for (sch, lab), col in list(zip(DOSES, ramp))[1:]:
        dose = {"const_0_62": 0.62, "const_1": 1.0, "const_1_47": 1.47}[sch]
        d, v, s, ex = series(rows, condition="treated", line=line, culture="mono", schedule=sch)
        key = f"mono_treated_{line}_{dose}"
        pts(ax, d, v, s, ex, col, lab)
        t, y = rk4(f, [st.get(key, 100.0), 0, 0, 0], 14, const(dose))
        ax.plot(t, y[:, 0], color=col)
    label_panel(ax, title); cells_axis(ax, 4800)
    ax.legend(loc="upper left", title="Dose (lines = model)", title_fontsize=8)
axs[0].set_ylabel("Cells per image")
fig.tight_layout(); finish(fig, "fig3_3_stage2b")

# ---------------------------------------------------------------- Fig 3.4 Stage 3a
fig, axs = plt.subplots(2, 3, figsize=(9, 5.2), sharex=True, sharey=True)
for i, (cond, D) in enumerate((("Ce0", 0.0), ("Ce1", 1.0))):
    for j, (rk, rl) in enumerate(RATIOS):
        ax = axs[i, j]
        y0 = [st[f"co_{cond}_{rk}_A2780"], st[f"co_{cond}_{rk}_A2780cis"], 0, 0, 0, 0, 0, 0]
        t, y = rk4(cocult, y0, 14, const(D))
        for k, (line, col) in enumerate((("A2780", SENS), ("A2780cis", RES))):
            d, v, s, ex = series(rows, condition=cond, line=line, culture="co", ratio=rk)
            pts(ax, d, v, s, ex, col, line)
            ax.plot(t, y[:, k], color=col)
        label_panel(ax, f"{'abc'[j] if i == 0 else 'def'[j]}  {'Untreated' if D == 0 else '1 µM'}, S:R {rl}")
        cells_axis(ax, 4800, xlabel=i == 1)
for ax in axs[:, 0]:
    ax.set_ylabel("Cells per image")
axs[0, 0].legend(loc="upper left")
fig.tight_layout(); finish(fig, "fig3_4_stage3a")

# ---------------------------------------------------------------- Fig 3.5 Stage 3b + 4
fig, axs = plt.subplots(1, 4, figsize=(10, 2.8), sharey=True)
for j, (rk, rl) in enumerate(RATIOS):
    ax = axs[j]
    y0 = [st[f"co_lowres_run1_{rk}_A2780"], st[f"co_lowres_run1_{rk}_A2780cis"], 0, 0, 0]
    t, y = rk4(cocult_lr, y0, 14, const(0))
    for k, (line, col) in enumerate((("A2780", SENS), ("A2780cis", RES))):
        d, v, s, ex = series(rows, condition="low_resource", line=line, culture="co", ratio=rk, run="1")
        pts(ax, d, v, s, ex, col, line); ax.plot(t, y[:, k], color=col)
    label_panel(ax, f"{'abc'[j]}  Stage 3b · low resource, {rl}"); cells_axis(ax, 4800)
ax = axs[3]
y0 = [st["co_combination_50S50R_A2780"], st["co_combination_50S50R_A2780cis"], 0, 0, 0, 0, 0, 0]
t, y = rk4(cocult_s4, y0, 14, const(1.0))
for k, (line, col) in enumerate((("A2780", SENS), ("A2780cis", RES))):
    d, v, s, ex = series(rows, condition="combination", line=line, culture="co")
    pts(ax, d, v, s, ex, col, line); ax.plot(t, y[:, k], color=col)
label_panel(ax, "d  Stage 4 · low resource + 1 µM"); cells_axis(ax, 4800)
axs[0].set_ylabel("Cells per image"); axs[3].legend(loc="upper right")
fig.tight_layout(); finish(fig, "fig3_5_stage3b_4")

# ---------------------------------------------------------------- Fig 3.6 stability + sensitivity
sel = [r for r in csv.DictReader(open(f"{PKGO}/zoo_analyses/bootstrap_selection.csv")) if r["optimal"] == "true"]
names = {"s1-a2780": "1 · A2780", "s1-a2780cis": "1 · A2780cis", "s1b-a2780": "2a · A2780", "s1b-a2780cis": "2a · A2780cis",
         "s2-a2780": "2b · A2780", "s2-a2780cis": "2b · A2780cis", "s3": "3a · coculture", "s3b": "3b · low-res. coculture",
         "s4": "4 · low res. + cisplatin"}
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.4), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]
yy = np.arange(len(sel))[::-1]
ch = [float(r["pct_chosen"]) for r in sel]; lo = [float(r["pct_lowest_bic"]) for r in sel]
ax.barh(yy + 0.19, ch, 0.36, color=SENS, label="Selected by the rule")
ax.barh(yy - 0.19, lo, 0.36, color=AQUA, label="Lowest BIC")
for y_, c in zip(yy, ch):
    ax.text(c + 1.5, y_ + 0.19, f"{c:.0f}%", va="center", fontsize=7.5, color="#1f2933")
ax.set_yticks(yy); ax.set_yticklabels([names[r["group"]] for r in sel])
ax.set_xlim(0, 110); ax.set_xlabel("% of 200 bootstrap replicates"); ax.grid(axis="y", visible=False)
ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2, fontsize=7.5); label_panel(ax, "a  How often the selected model is picked again")
loc = list(csv.DictReader(open(f"{PKGO}/zoo_analyses/sensitivity_local.csv")))
loc.sort(key=lambda r: -float(r["si"])); loc = loc[:10]
stage = {"s1-a2780": "1", "s1-a2780cis": "1", "s1b-a2780": "2a", "s1b-a2780cis": "2a", "s2-a2780": "2b",
         "s2-a2780cis": "2b", "s3": "3a", "s3b": "3b", "s4": "4"}
ax = axs[1]
yy = np.arange(len(loc))[::-1]
ax.barh(yy, [float(r["si"]) for r in loc], 0.6, color=SENS)
for y_, r in zip(yy, loc):
    ax.text(float(r["si"]) + 0.6, y_, f"{float(r['si']):.1f}", va="center", fontsize=7.5)
ax.set_yticks(yy); ax.set_yticklabels([f"{r['param']}  ({stage[r['group']]})" for r in loc])
ax.set_xlabel("Local sensitivity index"); ax.grid(axis="y", visible=False); ax.set_xlim(0, 44)
label_panel(ax, "b  Ten most sensitive parameters (stage)")
fig.tight_layout(); finish(fig, "fig3_6_stability")

# ---------------------------------------------------------------- Fig 3.7 validation
fig, axs = plt.subplots(1, 4, figsize=(10.5, 2.9), gridspec_kw={"width_ratios": [1, 1, 0.9, 1.1]})
for ax, line, f, title in ((axs[0], "A2780", drug_N, "a  Intermittent dosing · A2780"),
                           (axs[1], "A2780cis", drug_C, "b  Intermittent dosing · A2780cis")):
    for sch, col, lab, dose_fn in (("on0_off7", SENS if line == "A2780" else RES, "1 µM days 0–7, then off", lambda t: 1.0 if t < 7 else 0.0),
                                   ("off0_on7", AQUA, "Off, then 1 µM from day 7", lambda t: 0.0 if t < 7 else 1.0)):
        d, v, s, ex = series(rows, condition="intermittent", line=line, schedule=sch)
        pts(ax, d, v, s, ex, col, lab)
        t, y = rk4(f, [100, 0, 0, 0], 14, dose_fn)
        ax.plot(t, y[:, 0], color=col)
    ax.axvline(7, color=GRAY, lw=0.8, ls="--")
    label_panel(ax, title); cells_axis(ax, 4800)
axs[0].set_ylabel("Cells per image"); axs[0].legend(loc="upper left", fontsize=7)
ax = axs[2]
vals = [("Fitted", 366, SENS), ("Predicted", 1132, RES), ("No-drug\nref.", 975, GRAY)]
ax.bar(range(3), [v for _, v, _ in vals], 0.6, color=[c for *_, c in vals])
for i, (_, v, _) in enumerate(vals):
    ax.text(i, v + 25, f"{v:,}", ha="center", fontsize=7.5)
ax.set_xticks(range(3)); ax.set_xticklabels([l for l, *_ in vals], fontsize=7.5); ax.grid(axis="x", visible=False)
ax.set_ylabel("RMSE (cells)"); label_panel(ax, "c  Held-out treated coculture")
ax = axs[3]
fr = [("1 · A2780", 2.05), ("1 · A2780cis", 1.32), ("2a · A2780", 2.18), ("2a · A2780cis", 1.11), ("2b · A2780", 1.41),
      ("2b · A2780cis", 2.53), ("3a", 1.29), ("3b", 1.65), ("4", 1.32)]
yy = np.arange(len(fr))[::-1]
ax.hlines(yy, 1, [v for _, v in fr], color="#cbd2d9", lw=1.2)
ax.plot([v for _, v in fr], yy, "o", color=SENS, ms=6)
ax.axvline(1, color=GRAY, lw=0.8)
ax.set_yticks(yy); ax.set_yticklabels([l for l, _ in fr], fontsize=7.5); ax.grid(axis="y", visible=False)
ax.set_xlabel("Forecast RMSE ÷ in-sample RMSE"); ax.set_xlim(0.8, 2.8); label_panel(ax, "d  Forecasting the last days")
fig.tight_layout(); finish(fig, "fig3_7_validation")

# ---------------------------------------------------------------- Fig 3.8 control 14 d
POL = [("No treatment", [0] * 7, GRAY, "-"), ("Constant 1 µM", [1] * 7, RES, "-"),
       ("Periodic 2 on / 2 off", [1, 0, 1, 0, 1, 0, 1], YELLOW, "-"),
       ("Threshold adaptive", None, AQUA, "-"), ("MPC", None, SENS, "-")]
outc = list(csv.DictReader(open(f"{PKGO}/control/control_outcomes.csv")))
def sched(scen, pol):
    for r in outc:
        if r["environment"] == "normal" and r["scenario"] == scen and r["ratio"] == "50:50" and r["policy"].startswith(pol):
            return [float(x) for x in r["schedule_uM"].split()]
fig, axs = plt.subplots(3, 2, figsize=(8.5, 6.6), sharex=True, gridspec_kw={"height_ratios": [0.75, 1, 1]})
for j, scen in enumerate("AB"):
    for name, lv, col, ls in POL:
        if name == "Threshold adaptive":
            lv = sched(scen, "Threshold")
        if name == "MPC":
            lv = sched(scen, "MPC, measured")
        PB = dict(P); PB.update(K3=4552.468382416471, aNC3=0.7448137076008907, aCN3=1.5177611370618196,
                                dN=0.19751551419039187, dC=0.07718422648569812)   # Ce0-only crowding fit
        kw = {} if scen == "A" else {"p": PB, "betaC": P["bC2"], "betaN": P["bN2"]}
        t, y = rk4(lambda yy, D: cocult(yy, D, **kw), [50, 50, 0, 0, 0, 0, 0, 0], 14, schedule(lv))
        lw = 2.6 if name == "MPC" else 1.5
        axs[0, j].step(np.arange(0, 16, 2), lv + [lv[-1]], where="post", color=col, lw=lw, label=name)
        axs[1, j].plot(t, y[:, 0] + y[:, 1], color=col, lw=lw, label=name)
        axs[2, j].plot(t, y[:, 1], color=col, lw=lw, label=name)
    axs[1, j].axhline(0.8 * P["K3"], color="#1f2933", lw=0.8, ls="--")
    axs[1, j].text(0.3, 0.8 * P["K3"] + 60, "Progression (0.8 K)", fontsize=7, color="#1f2933")
    axs[0, j].set_title(f"{'ab'[j]}  Scenario {scen}: {'no resistant kill (Stage 3a)' if scen == 'A' else 'monoculture kill for both lines'}", loc="left")
    axs[0, j].set_ylim(-0.1, 1.6); axs[1, j].set_ylim(0, 4600); axs[2, j].set_ylim(0, 4200)
    axs[2, j].set_xlabel("Day"); axs[2, j].set_xticks(range(0, 15, 2))
axs[0, 0].set_ylabel("Dose (µM)"); axs[1, 0].set_ylabel("Total cells (N + C)"); axs[2, 0].set_ylabel("Resistant cells (C)")
axs[2, 1].legend(loc="upper left", fontsize=7.5)
fig.tight_layout(); finish(fig, "fig3_8_control14")

# ---------------------------------------------------------------- Fig 3.9 control 120 d
pol120 = [("No treatment", 3476, 0), ("Constant 1 µM", 4303, 4245), ("Maximum 1.47 µM", 4331, 4331),
          ("Periodic 2 on / 2 off", 3189, 0), ("Threshold adaptive", 4256, 4103), ("Fasting only", 2330, 0),
          ("Periodic + fasting", 2302, 0), ("Constant + fasting", 1109, 0), ("MPC, dose only", 2944, 13),
          ("MPC, dose + fasting", 2814, 0)]
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6), gridspec_kw={"width_ratios": [1.3, 1]})
ax = axs[0]
yy = np.arange(len(pol120))[::-1]
ax.barh(yy, [a for _, a, _ in pol120], 0.62, color=[SENS if n.startswith("MPC") else "#9ec5f4" for n, *_ in pol120], label="Total cells")
ax.barh(yy, [c for *_, c in pol120], 0.62, color=RES, label="of which resistant")
for y_, (n, a, c) in zip(yy, pol120):
    ax.text(a + 60, y_, f"{a:,}", va="center", fontsize=7.5)
ax.axvline(0.8 * P["K3"], color="#1f2933", lw=0.8, ls="--")
ax.set_yticks(yy); ax.set_yticklabels([n for n, *_ in pol120]); ax.grid(axis="y", visible=False)
ax.set_xlim(0, 5000); ax.set_xlabel("Cells on day 120 (scenario A, 50:50)")
ax.legend(loc="lower right", fontsize=7.5); label_panel(ax, "a  120 days with cisplatin and fasting")
ax = axs[1]
ko = ["0 (fitted)", "0.5", "2"]
x = np.arange(3)
ax.bar(x - 0.2, [0, 4239, 4288], 0.38, color=RES, label="Fasting plan replayed")
ax.bar(x + 0.2, [0, 14, 9], 0.38, color=SENS, label="Re-optimised")
for xi, (a, b, fd) in zip(x, [(0, 0, 8), (4239, 14, 4), (4288, 9, 0)]):
    ax.text(xi - 0.2, a + 60, f"{a:,}", ha="center", fontsize=7.5)
    ax.text(xi + 0.2, b + 120, f"{b} cells, {fd} fast days", ha="center", va="bottom", rotation=90, fontsize=7)
ax.set_xticks(x); ax.set_xticklabels(ko); ax.set_xlabel("Recovery rate after refeeding, k_off (/day)")
ax.set_ylabel("Resistant cells on day 120"); ax.set_ylim(0, 5000); ax.grid(axis="x", visible=False)
ax.legend(loc="upper left", fontsize=7.5); label_panel(ax, "b  The fasting benefit needs slow recovery")
fig.tight_layout(); finish(fig, "fig3_9_control120")

# ---------------------------------------------------------------- Fig 3.10 NN vs mechanistic
tests = [("Next-day, treated coculture", 0.2245, 0.2501, "Neural"), ("Held-out Ce1 (20k + 30k)", 0.3897, 0.3944, "Tie"),
         ("Unseen schedule on→off", 0.1968, 0.2348, "Neural"), ("Unseen schedule off→on", 0.3555, 0.3853, "Tie"),
         ("Unseen dose (1 µM)", 0.217, 0.2034, "Tie"), ("Unseen density (20k)", 0.4369, 0.4387, "Tie"),
         ("Three-way, held-out 30k Ce1", 0.3479, 0.6372, "Neural")]
fig, ax = plt.subplots(figsize=(6.8, 3.0))
yy = np.arange(len(tests))[::-1]
for y_, (n, nn, me, w) in zip(yy, tests):
    ax.plot([nn, me], [y_, y_], color="#cbd2d9", lw=1.5, zorder=1)
    ax.text(0.68, y_, w, va="center", fontsize=7.5, color="#1f2933" if w != "Tie" else "#52606d")
ax.scatter([t[1] for t in tests], yy, s=36, color=SENS, label="Best neural", zorder=2)
ax.scatter([t[2] for t in tests], yy, s=36, color=RES, marker="s", label="Best mechanistic", zorder=2)
ax.set_yticks(yy); ax.set_yticklabels([t[0] for t in tests]); ax.grid(axis="y", visible=False)
ax.set_xlim(0.15, 0.74); ax.set_xlabel("Held-out log RMSE (lower is better)")
ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=2, fontsize=7.5); ax.set_title("Neural vs mechanistic on held-out tests (winner at right)", loc="left")
fig.tight_layout(); finish(fig, "fig3_10_nn")
print("done", sorted(os.listdir(OUT)))
