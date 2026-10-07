#!/usr/bin/env python3
"""
Validation of the steady-state flow field and of the flow-field associations.

    python3 src/flow_field_validation.py --root runs_gravityoff --old-root runs_dirichlet_20260924_setonix --old-prefix A_

(1) Does flowfield.solve() describe the flow each simulation actually had?
    The simulated flux of every connection comes from the run's own pressure field
    at the first output >= 0.05 y:  q = T [(p_i + rho g z_i) - (p_j + rho g z_j)],
    with rho g = 0 when the deck has gravity off. T is flowfield's transmissibility.
    It is compared with solve()'s flux: Spearman of |q| over connections, the
    flux-weighted share of connections flowing in the same direction, and the ratio
    of total |q| (circulation adds flux that injection alone does not produce).
    Done for the corrected runs and for the circulating runs (--old-root).

(2) Is "old water" only "far water"? Cells are split into 10 bands of equal volume by
    distance from the inlet face (x). Within each band, each cell's age percentile
    is computed within that band. The carbonate-weighted mean of the within-band
    percentile (50 = no preference) measures the age preference at fixed distance.
    The carbonate-weighted distance percentile is reported alongside.

(3) Does shape99 predict efficiency beyond network size? Partial Spearman correlation
    of shape99 and efficiency, controlling for pore volume and number of cells.

Writes results/flow_field_validation.csv.
"""
import argparse, csv, glob, os, re, sys
import numpy as np, h5py
from scipy.stats import spearmanr, rankdata, t as tdist

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb
import flowfield
from inlet_sensitivity import shape99

G = 9.8068


def gravity_on(d):
    s = open(os.path.join(d, "pflotran_co2.in")).read()
    m = re.search(r"^\s*GRAVITY\s+([-0-9.eEdD+]+)\s+([-0-9.eEdD+]+)\s+([-0-9.eEdD+]+)", s, re.M | re.I)
    if not m:
        return True                                   # PFLOTRAN default: gravity on
    return any(abs(float(x.replace("d", "e").replace("D", "e"))) > 0 for x in m.groups())


def simulated_flux(d, xyz, ids, T):
    snaps = cb.snapshots(d)
    t, p, k = next((s for s in snaps if s[0] >= 0.05), snaps[-1])
    with h5py.File(p, "r") as f:
        g = f[k]
        pn = next(n for n in g if n.startswith("Liquid Pressure"))
        P = np.asarray(g[pn][:], float).ravel()
        dn = next((n for n in g if n.startswith("Liquid Density")), None)
        rho = np.asarray(g[dn][:], float).ravel() if dn else np.full(len(P), 997.0)
    rg = G * (rho[ids[:, 0]] + rho[ids[:, 1]]) / 2 if gravity_on(d) else 0.0
    return T * ((P[ids[:, 0]] - P[ids[:, 1]]) + rg * (xyz[ids[:, 0], 2] - xyz[ids[:, 1], 2])), t


def check1(d):
    cells, ids, area = flowfield.read_uge(os.path.join(d, "full_mesh.uge"))
    xyz = cells[:, :3]
    dist = np.linalg.norm(xyz[ids[:, 0]] - xyz[ids[:, 1]], axis=1); dist[dist <= 0] = np.finfo(float).eps
    T = flowfield.PERM * area / (flowfield.MU * dist)
    sol = flowfield.solve(d, verbose=False)
    q_sim, t = simulated_flux(d, xyz, ids, T)
    # solve() flux along ids, signed i -> j
    up = np.asarray(sol["up"]); q_ab = np.asarray(sol["flux"], float)
    q_sol = np.where(up == ids[:, 0], q_ab, -q_ab)
    rho = spearmanr(np.abs(q_sim), np.abs(q_sol)).correlation
    same = float(np.sum(np.abs(q_sol) * (np.sign(q_sim) == np.sign(q_sol))) / np.sum(np.abs(q_sol)))
    return dict(t_used=t, rho_absflux=rho, same_direction=same, flux_ratio=float(np.abs(q_sim).sum() / np.abs(q_sol).sum()),
                gravity=gravity_on(d)), sol, xyz


def within_band(values, vol, mask):
    """volume-weighted percentile of each cell's value within the cells of mask"""
    idx = np.where(mask)[0]; o = idx[np.argsort(values[idx], kind="stable")]
    pv = vol[o] / vol[o].sum(); pct = np.empty(len(vol)); pct[:] = np.nan
    pct[o] = 100 * (np.cumsum(pv) - pv / 2)
    return pct


def check2(d, sol, xyz):
    vol = np.asarray(sol["vol"], float)
    age = np.where(np.asarray(sol["stagnant"], bool) | ~np.isfinite(sol["age"]), np.inf, np.asarray(sol["age"], float))
    t, p, k = cb.snapshots(d)[-1]; mv = cb.molar_volumes(d); mol = np.zeros(len(vol))
    with h5py.File(p, "r") as f:
        g = f[k]
        for m in cb.CARB:
            n = next((n for n in g if " VF" in n and n.split(" VF")[0] == m), None)
            if n is not None:
                mol += np.clip(np.asarray(g[n][:], float).ravel() - cb.SEED, 0, None) * vol / mv[m]
    w = mol / mol.sum()
    x = xyz[:, 0]; xpct = within_band(x, vol, np.ones(len(vol), bool))
    band = np.minimum((xpct / 10).astype(int), 9)
    apct = np.full(len(vol), np.nan)
    for b in range(10):
        m = band == b
        if m.sum() > 1:
            apct[m] = within_band(age, vol, m)[m]
    agepct_all = within_band(age, vol, np.ones(len(vol), bool))
    return dict(age_pct_carb=float(np.nansum(w * agepct_all)), dist_pct_carb=float(np.sum(w * xpct)),
                age_pct_within_band=float(np.nansum(w * apct)))


def partial_spearman(x, y, controls):
    rx, ry = rankdata(x), rankdata(y)
    if np.ptp(rx) == 0 or np.ptp(ry) == 0:
        return float("nan"), float("nan")             # a constant variable has no correlation
    C = np.column_stack([np.ones(len(x))] + [rankdata(c) for c in controls])
    ex = rx - C @ np.linalg.lstsq(C, rx, rcond=None)[0]; ey = ry - C @ np.linalg.lstsq(C, ry, rcond=None)[0]
    if ex.std() < 1e-9 * rx.std() or ey.std() < 1e-9 * ry.std():
        return float("nan"), float("nan")             # fully explained by the controls
    r = float(np.corrcoef(ex, ey)[0, 1]); dof = len(x) - 2 - len(controls)
    tval = r * np.sqrt(dof / max(1e-12, 1 - r * r))
    return r, float(2 * tdist.sf(abs(tval), dof))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs_gravityoff")
    ap.add_argument("--blocks", nargs="+", default=["B", "A"])
    ap.add_argument("--old-root", default=None, help="circulating runs, for check (1)")
    ap.add_argument("--old-prefix", default="A_")
    ap.add_argument("--old-n", type=int, default=10)
    a = ap.parse_args()
    R = lambda r: r if os.path.isabs(r) else os.path.join(REV, r)
    out = []
    if a.old_root:
        dirs = sorted(d for d in glob.glob(os.path.join(R(a.old_root), a.old_prefix + "p32_*"))
                      if os.path.isfile(os.path.join(d, "boundary_right_e.ex")) and cb.snapshots(d))[:a.old_n]
        print(f"(1) circulating runs, {len(dirs)} networks: simulated flux vs solve() flux")
        for d in dirs:
            try:
                c1, _, _ = check1(d)
            except Exception as e:
                print(f"  skipped {os.path.basename(d)}: {e}"); continue
            print(f"  {os.path.basename(d):32s} gravity {'on ' if c1['gravity'] else 'off'}  rho(|q|) {c1['rho_absflux']:+.2f}  "
                  f"same direction {100 * c1['same_direction']:5.1f}%  total |q| ratio {c1['flux_ratio']:9.3g}  (t = {c1['t_used']:g} y)")
            out.append(dict(set="circulating", network=os.path.basename(d), **c1))
    for blk in a.blocks:
        eff = {}
        p = os.path.join(REV, "results", f"{blk}_coupled_networks.csv")
        if os.path.isfile(p):
            eff = {r["network"]: float(r["efficiency_pct"]) for r in csv.DictReader(open(p))}
        print(f"\nBlock {blk} (corrected runs)")
        rows = []
        for d in sorted(glob.glob(os.path.join(R(a.root), f"{blk}_feedback__p32_*"))):
            name = os.path.basename(d)
            try:
                c1, sol, xyz = check1(d); c2 = check2(d, sol, xyz)
            except Exception as e:
                print(f"  skipped {name}: {e}"); continue
            r = dict(set=f"corrected_{blk}", network=name, shape99=float(shape99(sol["age"])),
                     pore_volume_m3=sol["pore_volume_m3"], n_cells=sol["n_cells"], efficiency_pct=eff.get(name, np.nan), **c1, **c2)
            rows.append(r)
            print(f"  {name:30s} (1) rho(|q|) {c1['rho_absflux']:+.2f}  same dir {100 * c1['same_direction']:5.1f}%  |q| ratio {c1['flux_ratio']:6.3f}"
                  f"   (2) age pct {c2['age_pct_carb']:5.1f}, within distance band {c2['age_pct_within_band']:5.1f}, distance pct {c2['dist_pct_carb']:5.1f}")
        if rows:
            med = lambda k: np.nanmedian([r[k] for r in rows])
            print(f"  medians: (1) rho(|q|) {med('rho_absflux'):+.2f}, same direction {100 * med('same_direction'):.1f}%, |q| ratio {med('flux_ratio'):.3f}")
            print(f"           (2) carbonate-weighted age percentile {med('age_pct_carb'):.1f}; within distance bands {med('age_pct_within_band'):.1f}"
                  f" (50 = no age preference at fixed distance); distance percentile {med('dist_pct_carb'):.1f}")
            ok = [r for r in rows if np.isfinite(r["efficiency_pct"])]
            if len(ok) >= 6:
                x = [r["shape99"] for r in ok]; y = [r["efficiency_pct"] for r in ok]
                pv = [r["pore_volume_m3"] for r in ok]; nc = [r["n_cells"] for r in ok]
                rs = spearmanr(x, y); rp, pp = partial_spearman(x, y, [pv, nc])
                print(f"           (3) shape99 vs efficiency: rho {rs.correlation:+.3f} (p = {rs.pvalue:.2g}); "
                      f"controlling pore volume and cell count: {rp:+.3f} (p = {pp:.2g}); n = {len(ok)}")
                print(f"               shape99 vs pore volume rho {spearmanr(x, pv).correlation:+.2f}, vs cell count {spearmanr(x, nc).correlation:+.2f}")
            out += rows
    if out:
        keys = sorted({k for r in out for k in r}, key=lambda k: (k not in ("set", "network"), k))
        p = os.path.join(REV, "results", "flow_field_validation.csv")
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(out)
        print(f"\nwritten: results/{os.path.basename(p)}")


if __name__ == "__main__":
    main()
