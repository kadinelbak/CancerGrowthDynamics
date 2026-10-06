"""Write the stages.toml copies for the model-assumption comparison (outputs/sensitivity/<variant>/stages.toml).

Only configuration changes, never code: each variant is the main stages.toml with one or two edits.
  V1_floor0     resistant kill in coculture refitted with no floor (betaC_floor = 0)
  V2_floor25    floor at 25% of the Stage 2b value
  V3_floor75    floor at 75%
  V4_ns         Norton-Simon: drug kill x growth factor (1 - crowding/K) in Stages 2b, 3a and 4
  V5_ns_sd      V4 with residuals weighted by the replicate SD
  V6_sd         the current model with residuals weighted by the replicate SD
  V7_ns_sd_floor0  V5 with no floor on the coculture resistant kill
  V8_sd_floor0     V6 with no floor
The current model (V0) is the main run itself (outputs/zoo_fits), so it is not refitted.
"""
import os, re

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = open(os.path.join(PKG, 'stages.toml'), encoding='utf-8').read()


def floor(s, f):
    s2, n = re.subn(r'(?m)^betaC_floor = [0-9.]+', f'betaC_floor = {f}', s)
    assert n == 1
    return s2


def ns(s):
    for sid in ('s2', 's3', 's4'):
        s, n = re.subn(rf'(?m)^(id = "{sid}"\n)', rf'\1kill_growth_coupled = true   # Norton-Simon variant\n', s)
        assert n == 1, sid
    return s


def sd(s):
    s2, n = re.subn(r'(?m)^(rule_bic = [0-9.]+.*\n)', r'\1residual_weighting = "sd"   # variant: per-point weights from the replicate SD\n', s)
    assert n == 1
    return s2


VARIANTS = {
    'V4_ns': ns(BASE),
    'V1_floor0': floor(BASE, 0.0),
    'V5_ns_sd': sd(ns(BASE)),
    'V6_sd': sd(BASE),
    'V2_floor25': floor(BASE, 0.25),
    'V3_floor75': floor(BASE, 0.75),
    'V7_ns_sd_floor0': floor(sd(ns(BASE)), 0.0),
    'V8_sd_floor0': floor(sd(BASE), 0.0),
}

if __name__ == '__main__':
    for name, text in VARIANTS.items():
        d = os.path.join(PKG, 'outputs', 'sensitivity', name)
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, 'stages.toml'), 'w', encoding='utf-8', newline='\n').write(text)   # LF, as the original
        print(name)
