#!/usr/bin/env python3
"""Steady Darcy flow + mean groundwater age from a PFLOTRAN .uge mesh.

Gives the residence-time field for an a priori trapping index -- computable
BEFORE any reactive run, from the archived meshes alone.

Mean age (Goode 1996, WRR 32:289): integrating div(q*a) = phi over cell i,
    sum_out(q) * a_i - sum_in(q) * a_j = V_i * phi,   a = 0 in injected water.

Limitations: inlet is a distributed mass SOURCE (as run_pflotran.py does), not
fixed head; first-order upwind so a_i is the age at the cell's OUTFLOW face,
~half a cell above the cell average; mean age is a first moment, so a bimodal
transit-time spread is invisible to it.
"""
import os, sys, argparse
import numpy as np

PERM, POROSITY = 8.33e-8, 0.50      # run_pflotran.py:248, 526
MU, RHO = 5.47e-4, 988.0            # water at 50 C
INJ_FRACTION = 0.20                 # left 20% of domain


def read_uge(path):
    with open(path) as f:
        h = f.readline().split()
        if h[0].upper() != "CELLS":
            raise ValueError(f"{path}: expected 'CELLS n', got {h[:2]}")
        n = int(h[1])
        cells = np.empty((n, 4))
        for i in range(n):
            p = f.readline().split()
            cells[i] = (float(p[1]), float(p[2]), float(p[3]), float(p[4]))
        h = f.readline().split()
        if h[0].upper() != "CONNECTIONS":
            raise ValueError(f"{path}: expected 'CONNECTIONS m', got {h[:2]}")
        m = int(h[1])
        ids = np.empty((m, 2), dtype=np.int64); area = np.empty(m)
        for j in range(m):
            p = f.readline().split()
            ids[j] = (int(p[0]) - 1, int(p[1]) - 1)     # .uge is 1-based
            area[j] = float(p[5])
    return cells, ids, area


def read_ex(path):
    if not os.path.isfile(path):
        return np.empty(0, dtype=np.int64), np.empty(0)
    rows = []
    with open(path) as f:
        first = f.readline().split()
        if first and first[0].upper() != "CONNECTIONS":
            f.seek(0)
        for line in f:
            p = line.split()
            if len(p) >= 5:
                rows.append((int(p[0]) - 1, float(p[4])))
    if not rows:
        return np.empty(0, dtype=np.int64), np.empty(0)
    return (np.array([a for a, _ in rows], dtype=np.int64),
            np.array([b for _, b in rows]))


def solve(dfn_dir, verbose=True, inj_fraction=INJ_FRACTION, water_rate_kg_s=0.01):
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve

    cells, ids, area = read_uge(os.path.join(dfn_dir, "full_mesh.uge"))
    n = len(cells); xyz, vol = cells[:, :3], cells[:, 3]
    d = np.linalg.norm(xyz[ids[:, 0]] - xyz[ids[:, 1]], axis=1)
    d[d <= 0] = np.finfo(float).eps
    T = PERM * area / (MU * d)

    out_idx, _ = read_ex(os.path.join(dfn_dir, "boundary_right_e.ex"))
    if len(out_idx) == 0:
        raise RuntimeError(f"{dfn_dir}: no usable boundary_right_e.ex -- "
                           "CLOSED-SYSTEM case, residence time undefined")
    xmin, xmax = xyz[:, 0].min(), xyz[:, 0].max()
    inlet = np.where(xyz[:, 0] <= xmin + inj_fraction * (xmax - xmin))[0]
    if len(inlet) == 0:
        raise RuntimeError(f"{dfn_dir}: injection region contains no cells")

    Q = water_rate_kg_s / RHO
    src = np.zeros(n); src[inlet] = Q * vol[inlet] / vol[inlet].sum()

    is_out = np.zeros(n, bool); is_out[out_idx] = True
    r, c, v = [], [], []; b = src.copy()
    for (i, j), t in zip(ids, T):
        for a_, b_ in ((i, j), (j, i)):
            if not is_out[a_]:
                r += [a_, a_]; c += [a_, b_]; v += [t, -t]
    for i in out_idx:
        r.append(i); c.append(i); v.append(1.0); b[i] = 0.0
    p = spsolve(coo_matrix((v, (r, c)), shape=(n, n)).tocsr(), b)
    if not np.all(np.isfinite(p)):
        raise RuntimeError("pressure solve produced non-finite values")

    q = T * (p[ids[:, 0]] - p[ids[:, 1]])
    up = np.where(q >= 0, ids[:, 0], ids[:, 1])
    dn = np.where(q >= 0, ids[:, 1], ids[:, 0])
    aq = np.abs(q)
    i_out = np.zeros(n); np.add.at(i_out, up, aq)
    i_in = np.zeros(n); np.add.at(i_in, dn, aq)
    bnd = np.zeros(n)
    bnd[out_idx] = np.maximum(i_in[out_idx] + src[out_idx] - i_out[out_idx], 0.0)
    denom = i_out + bnd                      # TOTAL outflow; src must NOT be added

    r2, c2, v2 = [], [], []; rhs = vol * POROSITY
    stag = denom <= 0
    for k in np.where(~stag)[0]:
        r2.append(k); c2.append(k); v2.append(denom[k])
    for u, dn_, f in zip(up, dn, aq):
        if not stag[dn_]:
            r2.append(dn_); c2.append(u); v2.append(-f)
    for i in np.where(stag)[0]:
        r2.append(i); c2.append(i); v2.append(1.0); rhs[i] = 0.0
    age = spsolve(coo_matrix((v2, (r2, c2)), shape=(n, n)).tocsr(), rhs)
    age = np.where(np.isfinite(age) & (age >= 0), age, np.nan)

    qo = float(bnd.sum()); YR = 3.156e7
    res = dict(n_cells=n, n_conn=len(ids), n_inlet=len(inlet), n_outlet=len(out_idx),
               pore_volume_m3=float((vol * POROSITY).sum()),
               q_inlet_m3s=Q, q_outlet_m3s=qo,
               mass_balance_rel=abs(Q - qo) / Q,
               stagnant_frac=float(stag.mean()),
               pv_over_q_yr=float((vol * POROSITY).sum() / Q / YR),
               age_median_yr=float(np.nanmedian(age) / YR),
               age_p90_yr=float(np.nanpercentile(age, 90) / YR),
               age=age, flux=aq, vol=vol, up=up, dn=dn,
               inlet=inlet, outlet=out_idx, bnd_out=bnd, stagnant=stag)
    if verbose:
        print(f"  cells {n:,}  conn {len(ids):,}  inlet {len(inlet):,}  outlet {len(out_idx):,}")
        print(f"  pore volume    {res['pore_volume_m3']:.4e} m^3")
        print(f"  Q in / out     {Q:.4e} / {qo:.4e} m^3/s  (imbalance {res['mass_balance_rel']:.2e})")
        print(f"  PV/Q           {res['pv_over_q_yr']:.4f} yr  (bulk turnover)")
        print(f"  stagnant cells {100*res['stagnant_frac']:.1f}%")
        print(f"  mean age       median {res['age_median_yr']:.4f} yr, p90 {res['age_p90_yr']:.4f} yr")
        fw = float(np.nansum(bnd[out_idx] * age[out_idx]) / max(bnd[out_idx].sum(), 1e-30) / YR)
        print(f"  flux-wtd outlet age {fw:.4f} yr")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dfn_dir"); ap.add_argument("--npz"); ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    r = solve(a.dfn_dir)
    if r["mass_balance_rel"] > 1e-6:
        print(f"  WARNING: mass imbalance {100*r['mass_balance_rel']:.3f}%")
    if a.npz:
        np.savez_compressed(a.npz, **{k: v for k, v in r.items() if isinstance(v, np.ndarray)})
        print(f"  wrote {a.npz}")


if __name__ == "__main__":
    main()
