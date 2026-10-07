#!/usr/bin/env python3
"""
Compare the uniform-pressure run with a corrected test run of the same case
(hydrostatic boundary, option A, or gravity off, option B).

Reads full_mesh.uge, pflotran_co2.h5 and each run's own pflotran_co2.in, so the
result does not depend on the analysis scripts. The flux uses each run's own
gravity vector g from its deck (PFLOTRAN's default -9.8068 along z if absent):

    q(i->j) = area/d * [(p_i - p_j) + rho_ij * g . (x_j - x_i)]

    python3 src/compare_hydrotest.py --old <uniform-pressure run> --new <test run>

    Put this file in revision/src/ (or revision/), then from anywhere:
    python3 ~/co2-basalt/revision/src/compare_hydrotest.py
"""
import argparse, os, sys
import numpy as np
import h5py
from scipy.stats import spearmanr, rankdata

# revision/ is the parent of src/, whichever of the two holds this file
_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
SRC = os.path.join(REV, "src")

G, SEED = 9.80665, 1e-6
CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")


def read_uge(path):
    with open(path) as f:
        n = int(f.readline().split()[1])
        cells = np.loadtxt(f, max_rows=n)
        m = int(f.readline().split()[1])
        conn = np.loadtxt(f, max_rows=m)
    return (cells[:, 1:4], conn[:, 0].astype(int) - 1,
            conn[:, 1].astype(int) - 1, conn[:, 5])


def last_group(h):
    t = [k for k in h if "Time" in k]
    if not t:
        sys.exit("no time groups in the output")
    return h[max(t, key=lambda k: float(k.split("Time")[1].split()[0]))]


def field(g, prefix, n):
    k = next((k for k in g if k.startswith(prefix) and isinstance(g[k], h5py.Dataset)), None)
    return None if k is None else np.asarray(g[k][:], float).ravel()[:n]


def gravity(run):
    import re
    deck = open(os.path.join(run, "pflotran_co2.in")).read()
    m = re.search(r"^[ \t]*GRAVITY[ \t]+(\S+)[ \t]+(\S+)[ \t]+(\S+)", deck, re.M | re.I)
    if not m:
        return np.array([0.0, 0.0, -9.8068])
    return np.array([float(v.replace("d", "e").replace("D", "E")) for v in m.groups()])


def load(run):
    gv = gravity(run)
    xyz, i, j, area = read_uge(os.path.join(run, "full_mesh.uge"))
    n = len(xyz)
    with h5py.File(os.path.join(run, "pflotran_co2.h5")) as h:
        g = last_group(h)
        p = field(g, "Liquid Pressure", n)
        rho = field(g, "Liquid Density", n)
        rho = np.full(n, 998.0) if rho is None else rho
        carb = {}
        for m in CARB:
            k = next((k for k in g if k.lower().startswith(m.lower()) and "vf" in k.lower()
                      and isinstance(g[k], h5py.Dataset)), None)
            if k is not None:
                carb[m] = np.asarray(g[k][:], float).ravel()[:n] - SEED
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1); d[d <= 0] = np.finfo(float).eps
    q = area / d * ((p[i] - p[j]) + 0.5 * (rho[i] + rho[j]) * ((xyz[j] - xyz[i]) @ gv))
    through = np.zeros(n); np.add.at(through, i, np.abs(q)); np.add.at(through, j, np.abs(q))
    total = sum(carb.values()) if carb else np.zeros(n)
    return dict(xyz=xyz, i=i, j=j, q=q, p=p, through=0.5 * through, carb=carb, total=total,
                gv=gv, rho=float(np.mean(rho)))


def report(name, r):
    xyz, x, q, i, j = r["xyz"], r["xyz"][:, 0], r["q"], r["i"], r["j"]
    A = np.c_[xyz, np.ones(len(x))]
    c = np.linalg.lstsq(A, r["p"], rcond=None)[0]
    sel = r["total"] > 0
    print(f"\n=== {name}   (gravity {r['gv'][0]:g} {r['gv'][1]:g} {r['gv'][2]:g})")
    print(f"  fitted dP/dz {c[2]:9.0f} Pa/m   (equilibrium would be {r['rho'] * r['gv'][2]:,.0f})")
    ratios = []
    for X in (-6, -2, 2, 6):
        lo, hi = np.minimum(x[i], x[j]), np.maximum(x[i], x[j])
        cr = (lo < X) & (hi >= X)        # counts connections touching the plane once
        net = (q[cr] * np.sign(x[j][cr] - x[i][cr])).sum()
        ratios.append(np.abs(q[cr]).sum() / abs(net) if net else np.inf)
    print("  gross/net across x = -6, -2, 2, 6 m:", "  ".join(f"{v:8.1f}" for v in ratios),
          "   (1 = no circulation)")
    print(f"  precipitating cells {int(sel.sum()):6d}   mean net carbonate per cell {r['total'].mean():.3e}")
    if r["carb"]:
        tot = sum(v[v > 0].sum() for v in r["carb"].values())
        print("  shares:", "  ".join(f"{m} {100 * v[v > 0].sum() / tot:5.1f}%" for m, v in r["carb"].items()) if tot else "none")
    e = np.linspace(x.min(), x.max(), 6)
    print("  precipitating cells per x band (west to east):",
          [int((sel & (x >= a) & (x <= b)).sum()) for a, b in zip(e[:-1], e[1:])])
    if sel.any():
        rk = rankdata(r["through"]) / len(x) * 100
        print(f"  median flux rank of precipitating cells {np.median(rk[sel]):5.1f}%")
    return sel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", default=os.path.join(REV, "runs", "A_p32_100_s1181"))
    ap.add_argument("--new", default=os.path.join(REV, "runs_hydrotest", "A_p32_100_s1181"))
    ap.add_argument("--label", default="corrected (test)")
    a = ap.parse_args()
    for r in (a.old, a.new):
        for f in ("full_mesh.uge", "pflotran_co2.h5", "pflotran_co2.in"):
            if not os.path.isfile(os.path.join(r, f)):
                sys.exit(f"missing {os.path.join(r, f)}")
    old, new = load(a.old), load(a.new)
    s_old = report("uniform pressure, gravity on (current)", old)
    s_new = report(a.label, new)

    both = s_old & s_new
    print("\n=== change")
    print(f"  mean net carbonate per cell: new/old = {new['total'].mean() / old['total'].mean():.3g}"
          if old["total"].mean() else "  old run has no carbonate")
    print(f"  precipitating cells shared by both runs: {int(both.sum())} "
          f"of {int(s_old.sum())} old and {int(s_new.sum())} new")
    try:
        sys.path.insert(0, SRC); import flowfield
        sol = flowfield.solve(a.new, verbose=False)
        qf = np.abs(np.asarray(sol["flux"], float))
        up, dn = np.asarray(sol["up"], int), np.asarray(sol["dn"], int)
        ff = np.zeros(len(new["through"])); np.add.at(ff, up, qf); np.add.at(ff, dn, qf)
        print(f"  Spearman(new simulated flux, flowfield) = {spearmanr(new['through'], ff).correlation:.3f}"
              "   (near 1: the index and age describe the simulated flow)")
    except Exception as ex:
        print(f"  flowfield comparison skipped: {ex}")


if __name__ == "__main__":
    main()
