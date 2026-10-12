#!/usr/bin/env python3
"""
Local cation balance of the cells that form carbonate: are their cations released by local
dissolution, or advected from upstream?

    python3 src/local_cation_balance.py [--root runs_gravityoff] [--prefix A_feedback__]

For each network and interval (as in colocation_mapping.py), the precipitating cells are those
whose carbonate volume fraction increases by more than --min-increase during the interval.
In those cells, in moles:
  supply   Ca and Mg released by the net loss of anorthite (1 Ca), diopside (1 Ca, 1 Mg)
           and forsterite (2 Mg) during the interval; regrowth of a primary mineral
           reduces the supply
  demand   Ca in new calcite and Mg in new magnesite
Reported (medians of networks, per interval):
  coverage        sum over cells of min(supply, demand) / sum of demand: the share of the
                  cations in new carbonate that local dissolution in the same cell provides
                  (100 % = fully local; the rest must arrive by transport)
  supply/demand   summed over all precipitating cells
  Ca and Mg separately
  carbonate formed  share of the total carbonate at 50 y formed in the interval
With --stages, the intervals run consecutively from 0.01 to 50 y (early, transport-controlled stage
included). Writes results/local_cation_balance.csv (or _stages.csv).

Molar volumes (m3/mol) must match the thermodynamic database of the decks; edit MOLAR_VOL if not.
"""
import argparse, csv, glob, os, sys
import warnings
import numpy as np, h5py
warnings.filterwarnings("ignore", message="All-NaN")

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb

INTERVALS = [(1.0, 2.0), (5.0, 7.0), (10.0, 15.0), (20.0, 30.0)]
# --stages: consecutive intervals over the whole history, from the snapshot times of the decks
STAGES = [(0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 5.0),
          (5.0, 7.0), (7.0, 10.0), (10.0, 15.0), (15.0, 20.0), (20.0, 30.0), (30.0, 50.0)]
MOLAR_VOL = {"Anorthite": 100.79e-6, "Diopside": 66.09e-6, "Forsterite": 43.79e-6,
             "Calcite": 36.934e-6, "Magnesite": 28.018e-6}
RELEASE = {"Anorthite": (1, 0), "Diopside": (1, 1), "Forsterite": (0, 2)}      # (Ca, Mg) per mole


def vf(g, mineral, n):
    k = next((x for x in g if x.startswith(f"{mineral} VF")), None)
    return np.asarray(g[k][:], float).ravel()[:n] if k else np.zeros(n)


def at(snaps, t):
    s = min(snaps, key=lambda x: abs(x[0] - t))
    if abs(s[0] - t) > 0.01:
        raise ValueError(f"no snapshot at {t:g} y (nearest {s[0]:g})")
    return s


def network(d, min_inc):
    snaps = cb.snapshots(d); vol = cb.cell_volumes(d); n = len(vol); rows = []
    t0, p0, k0 = snaps[0]; tl, pl, kl = at(snaps, 50.0)
    with h5py.File(p0, "r") as f0, h5py.File(pl, "r") as fl:
        total = sum(((vf(fl[kl], m, n) - vf(f0[k0], m, n)) * vol / MOLAR_VOL[m]).sum() for m in ("Calcite", "Magnesite"))
    for ta, tb in INTERVALS:
        try:
            (_, pa, ka), (_, pb, kb) = at(snaps, ta), at(snaps, tb)
        except ValueError:
            continue
        with h5py.File(pa, "r") as fa, h5py.File(pb, "r") as fb:
            ga, gb = fa[ka], fb[kb]
            cc = (vf(gb, "Calcite", n) - vf(ga, "Calcite", n)) * vol / MOLAR_VOL["Calcite"]
            mg = (vf(gb, "Magnesite", n) - vf(ga, "Magnesite", n)) * vol / MOLAR_VOL["Magnesite"]
            dca = np.zeros(n); dmg = np.zeros(n)
            for m, (sca, smg) in RELEASE.items():
                mol = (vf(ga, m, n) - vf(gb, m, n)) * vol / MOLAR_VOL[m]        # net loss, mol
                dca += sca * mol; dmg += smg * mol
        carb_inc = (cc * MOLAR_VOL["Calcite"] + mg * MOLAR_VOL["Magnesite"]) / vol
        prec = carb_inc > min_inc
        if not prec.any():
            continue
        dem_ca, dem_mg = np.clip(cc, 0, None), np.clip(mg, 0, None)
        sup_ca, sup_mg = np.clip(dca, 0, None), np.clip(dmg, 0, None)
        p = prec
        dem, sup = dem_ca[p] + dem_mg[p], sup_ca[p] + sup_mg[p]
        cover = lambda s, d_: 100 * np.minimum(s, d_).sum() / d_.sum() if d_.sum() > 0 else float("nan")
        rows.append(dict(network=os.path.basename(d), t_start=ta, t_end=tb, cells=int(p.sum()),
                         demand_mol=float(dem.sum()), supply_mol=float(sup.sum()),
                         share_of_total_pct=float(100 * (np.clip(cc, 0, None).sum() + np.clip(mg, 0, None).sum()) / total) if total > 0 else float('nan'),
                         coverage_pct=cover(sup, dem), supply_over_demand=float(sup.sum() / dem.sum()) if dem.sum() > 0 else float("nan"),
                         coverage_ca_pct=cover(sup_ca[p], dem_ca[p]), coverage_mg_pct=cover(sup_mg[p], dem_mg[p]),
                         supply_over_demand_ca=float(sup_ca[p].sum() / dem_ca[p].sum()) if dem_ca[p].sum() > 0 else float("nan"),
                         supply_over_demand_mg=float(sup_mg[p].sum() / dem_mg[p].sum()) if dem_mg[p].sum() > 0 else float("nan")))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs_gravityoff"); ap.add_argument("--prefix", default="A_feedback__")
    ap.add_argument("--stages", action="store_true", help="consecutive intervals over the whole history (0.01 to 50 y)")
    ap.add_argument("--min-increase", type=float, default=1e-9, help="carbonate volume fraction increase that defines a precipitating cell")
    a = ap.parse_args()
    global INTERVALS
    if a.stages:
        INTERVALS = STAGES
    root = a.root if os.path.isabs(a.root) else os.path.join(REV, a.root)
    out = []
    for d in sorted(glob.glob(os.path.join(root, a.prefix + "p32_*"))):
        try:
            out += network(d, a.min_increase)
        except Exception as e:
            print(f"  skipped {os.path.basename(d)}: {e}")
    if not out:
        sys.exit("no networks")
    print(f"local cation balance of the precipitating cells ({len({r['network'] for r in out})} networks; medians, with the range across networks)")
    print(f"{'interval':>12} {'carbonate formed':>17} {'coverage (local share)':>24} {'supply/demand':>15} {'Ca coverage':>12} {'Mg coverage':>12}")
    for ta, tb in INTERVALS:
        rs = [r for r in out if r["t_start"] == ta]
        if not rs:
            continue
        f = lambda k: (np.nanmedian([r[k] for r in rs]), np.nanmin([r[k] for r in rs]), np.nanmax([r[k] for r in rs]))
        c, s, cca, cmg, sh = f("coverage_pct"), f("supply_over_demand"), f("coverage_ca_pct"), f("coverage_mg_pct"), f("share_of_total_pct")
        print(f"{ta:>5g}-{tb:<5g}y  {sh[0]:7.2f}% of total  {c[0]:6.1f}% ({c[1]:5.1f}-{c[2]:5.1f})   {s[0]:7.3g} ({s[1]:.3g}-{s[2]:.3g})  {cca[0]:6.1f}%      {cmg[0]:6.1f}%")
    p = os.path.join(REV, "results", "local_cation_balance_stages.csv" if a.stages else "local_cation_balance.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    print(f"written: {os.path.relpath(p, REV)}")


if __name__ == "__main__":
    main()
