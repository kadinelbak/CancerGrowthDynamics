"""Python re-implementation of the selected zoo models (point estimates), for drawing figures."""
import csv
import numpy as np

import pathlib
PKG = str(pathlib.Path(__file__).resolve().parents[3] / "Modeling_Approaches/02_mechanical_automatic_package")

P = dict(r_N=0.817091919028476, K_N=3941.2697475302666,
         r_C=0.44803069922918215, K_C=3797.7784680941418, th=5.486412073489721,
         c_N=0.5255520504496792, w=1.031484275095739, z=0.050000000000282,
         dLR=0.1613623397361551, uR=0.9548,
         bN2=0.8292943884778216, kN=0.7199177321323673,
         bC2=0.4848022223460061, kC=0.2591676734335748,
         K3=4378.533340475707, aNC3=0.679246019225933, aCN3=1.4705609356683877,
         dN=0.1683317918911323, dC=0.02593089159604885, bN3=0.38552720699401943,
         K3b=5289.523677129859, aNC3b=4.311741666876674, aCN3b=0.08249107184819118,
         chi=1.1519145721966415)


def rk4(f, y0, t_end, dose_fn, h=0.01):
    """Integrate dy/dt = f(y, D) with a piecewise-constant dose D = dose_fn(t)."""
    n = int(round(t_end / h))
    ys = np.empty((n + 1, len(y0)))
    y = np.array(y0, float)
    ys[0] = y
    for i in range(n):
        t = i * h
        D = dose_fn(t + 1e-9)
        k1 = f(y, D); k2 = f(y + h / 2 * k1, D); k3 = f(y + h / 2 * k2, D); k4 = f(y + h * k3, D)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ys[i + 1] = y
    return np.linspace(0, t_end, n + 1), ys


def const(d):
    return lambda t: d


def schedule(levels, step=2.0):
    return lambda t: levels[min(int(t // step), len(levels) - 1)]


# ---- monoculture ----------------------------------------------------------
def logistic(y, D, p=P):
    N = y[0]
    return np.array([p["r_N"] * N * (1 - N / p["K_N"])])


def theta(y, D, p=P):
    C = y[0]
    return np.array([p["r_C"] * C * (1 - (max(C, 0) / p["K_C"]) ** p["th"])])


def lr_N(y, D, p=P):  # second-order effect lowers K
    N, S, dS = y
    Keff = p["K_N"] / (1 + p["c_N"] * S)
    return np.array([p["r_N"] * N * (1 - N / Keff), dS,
                     p["w"] ** 2 * (p["uR"] - S) - 2 * p["z"] * p["w"] * dS])


def lr_C(y, D, p=P):  # cumulative exposure -> loss
    C, S = y
    g = p["r_C"] * C * (1 - (max(C, 0) / p["K_C"]) ** p["th"])
    return np.array([g - p["dLR"] * S * C, p["uR"]])


def drug_N(y, D, p=P):
    N, a1, a2, a3 = y
    free = 1 - N / p["K_N"]
    k = p["kN"]
    return np.array([p["r_N"] * N * free - free * p["bN2"] * a3 * N,
                     k * (D - a1), k * (a1 - a2), k * (a2 - a3)])


def drug_C(y, D, p=P, beta=None):
    C, a1, a2, a3 = y
    free = 1 - (max(C, 0) / p["K_C"]) ** p["th"]
    k = p["kC"]
    b = p["bC2"] if beta is None else beta
    return np.array([p["r_C"] * C * free - free * b * a3 * C,
                     k * (D - a1), k * (a1 - a2), k * (a2 - a3)])


# ---- coculture (Stage 3a; scenario A plant) --------------------------------
def cocult(y, D, p=P, betaC=0.0, betaN=None):
    N, C, n1, n2, n3, c1, c2, c3 = y
    bN = p["bN3"] if betaN is None else betaN
    fN = 1 - (N + p["aNC3"] * C) / p["K3"]
    fC = 1 - (max(C + p["aCN3"] * N, 0) / p["K3"]) ** p["th"]
    kN, kC = p["kN"], p["kC"]
    return np.array([
        p["r_N"] * N * fN - p["dN"] * N - fN * bN * n3 * N,
        p["r_C"] * C * fC - p["dC"] * C - fC * betaC * c3 * C,
        kN * (D - n1), kN * (n1 - n2), kN * (n2 - n3),
        kC * (D - c1), kC * (c1 - c2), kC * (c2 - c3)])


def cocult_lr(y, D, p=P):  # Stage 3b
    N, C, S, dS, SC = y
    Keff = p["K3b"] / (1 + p["c_N"] * S)
    fN = 1 - (N + p["aNC3b"] * C) / Keff
    fC = 1 - (max(C + p["aCN3b"] * N, 0) / p["K3b"]) ** p["th"]
    return np.array([p["r_N"] * N * fN,
                     p["r_C"] * C * fC - p["dLR"] * SC * C,
                     dS, p["w"] ** 2 * (p["uR"] - S) - 2 * p["z"] * p["w"] * dS, p["uR"]])


def cocult_s4(y, D, p=P):  # Stage 4: LR + cisplatin, chi on the sensitive kill
    N, C, S, dS, SC, n1, n2, n3 = y
    Keff = p["K3"] / (1 + p["c_N"] * S)
    fN = 1 - (N + p["aNC3"] * C) / Keff
    fC = 1 - (max(C + p["aCN3"] * N, 0) / p["K3"]) ** p["th"]
    kN = p["kN"]
    uhat = p["uR"] / 0.955
    return np.array([
        p["r_N"] * N * fN - p["dN"] * N - fN * (1 + p["chi"] * uhat) * p["bN3"] * n3 * N,
        p["r_C"] * C * fC - p["dC"] * C - p["dLR"] * SC * C,
        dS, p["w"] ** 2 * (p["uR"] - S) - 2 * p["z"] * p["w"] * dS, p["uR"],
        kN * (D - n1), kN * (n1 - n2), kN * (n2 - n3)])


# ---- data ------------------------------------------------------------------
def load_measurements():
    rows = list(csv.DictReader(open(f"{PKG}/data/measurements.csv")))
    for r in rows:
        for k in ("day", "value", "sd", "sem"):
            r[k] = float(r[k]) if r[k] not in ("", "NaN") else np.nan
    return rows


def series(rows, **kw):
    sel = [r for r in rows if all(r[k] == v for k, v in kw.items())]
    sel.sort(key=lambda r: r["day"])
    return (np.array([r["day"] for r in sel]), np.array([r["value"] for r in sel]),
            np.array([r["sem"] for r in sel]), np.array([r["excluded"] == "true" for r in sel]))


def starts():
    return {r["series"]: float(r["start"]) for r in csv.DictReader(open(f"{PKG}/outputs/zoo_fits/starting_values.csv"))}
