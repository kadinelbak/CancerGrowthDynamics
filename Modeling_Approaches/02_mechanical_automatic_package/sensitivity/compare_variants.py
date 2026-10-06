"""Compare the model-assumption variants: selected models, BIC, key parameters and the held-out Ce1 test.
   python sensitivity/compare_variants.py  ->  outputs/sensitivity/comparison.md (and printed)
V0 is the main run (outputs/zoo_fits). BICs compare only within one loss: unweighted (V0-V4) or SD-weighted (V5, V6).
"""
import csv, os

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(PKG, 'outputs', 'sensitivity')
VARS = [('V0_current', os.path.join(PKG, 'outputs', 'zoo_fits'), 'raw', 'floor 50% (current)'),
        ('V1_floor0', os.path.join(OUT, 'V1_floor0'), 'raw', 'no floor'),
        ('V2_floor25', os.path.join(OUT, 'V2_floor25'), 'raw', 'floor 25%'),
        ('V3_floor75', os.path.join(OUT, 'V3_floor75'), 'raw', 'floor 75%'),
        ('V4_ns', os.path.join(OUT, 'V4_ns'), 'raw', 'Norton-Simon kill (2b, 3a, 4)'),
        ('V6_sd', os.path.join(OUT, 'V6_sd'), 'sd', 'current, SD-weighted loss'),
        ('V5_ns_sd', os.path.join(OUT, 'V5_ns_sd'), 'sd', 'Norton-Simon, SD-weighted loss'),
        ('V8_sd_floor0', os.path.join(OUT, 'V8_sd_floor0'), 'sd', 'current, SD-weighted, no floor'),
        ('V7_ns_sd_floor0', os.path.join(OUT, 'V7_ns_sd_floor0'), 'sd', 'Norton-Simon, SD-weighted, no floor')]
GROUPS = ['s2-a2780', 's2-a2780cis', 's3', 's3b', 's4']
MAIN_HELDOUT = {'prediction': 711.1, 'in-sample': 351.6, 'reference': 870.3}   # run_heldout.jl on the main run


def rows(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f))


def summarize(d):
    sel = {r['group']: r for r in rows(os.path.join(d, 'zoo_model_selection.csv')) if r['optimal'] == 'true'}
    par = rows(os.path.join(d, 'zoo_parameters.csv'))
    pv = lambda g, name: next((float(r['estimate']) for r in par if r['group'] == g and r['candidate'] == sel[g]['candidate'] and r['param'] == name), None)
    ho = os.path.join(d, 'heldout_ce1.csv')
    held = {}
    if os.path.isfile(ho):
        for r in rows(ho):
            held[r.get('fit') or r.get('kind') or r.get('row') or list(r.values())[0]] = r
    return sel, pv, held


lines = ['| Variant | Loss | 2b structure | 3a model (k) | β_C (3a) | ε_N (4) | ΣBIC 2b+3a+4 | Held-out Ce1 RMSE |', '|---|---|---|---|---|---|---|---|']
for name, d, loss, desc in VARS:
    if not os.path.isfile(os.path.join(d, 'zoo_model_selection.csv')):
        lines.append(f'| {name} ({desc}) | {loss} | not run yet | | | | | |'); continue
    sel, pv, held = summarize(d)
    tot = sum(float(sel[g]['bic']) for g in ('s2-a2780', 's2-a2780cis', 's3', 's4') if g in sel)
    hp = None
    for k, r in held.items():
        if 'pred' in k.lower():
            hp = r
    hrmse = f"{float(hp['rmse']):.0f}" if hp and 'rmse' in hp else (f"{MAIN_HELDOUT['prediction']:.0f}" if name == 'V0_current' else '—')
    bc = pv('s3', 'β_C'); en = pv('s4', 'ε_N^R') or pv('s4', 'ε_N')
    lines.append(f"| {name} ({desc}) | {loss} | {sel['s2-a2780cis']['label']} | {sel['s3']['candidate']} ({sel['s3']['k']}) | "
                 f"{'' if bc is None else f'{bc:.3g}'} | {'' if en is None else f'{en:.3g}'} | {tot:.1f} | {hrmse} |")
txt = '\n'.join(lines)
open(os.path.join(OUT, 'comparison.md'), 'w', encoding='utf-8').write(txt + '\n')
print(txt)
