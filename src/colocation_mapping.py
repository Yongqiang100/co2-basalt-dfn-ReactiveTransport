#!/usr/bin/env python3
"""
Co-location test on the corrected runs: does carbonate form where dissolved
calcium or magnesium, dissolved carbon and a high pH coincide?

    python3 src/colocation_mapping.py --list runs_gravityoff/A_feedback__p32_100_s1181   # dataset names
    python3 src/colocation_mapping.py --root runs_gravityoff --prefix A_feedback__

For each network and each interval (1-2, 5-7, 10-15, 20-30 y, nearest snapshots):
  new carbonate   carbonate formed in each cell during the interval (mol, positive part)
  conditions      at the START of the interval (before that carbonate forms):
                  cations  --cations complexes (default): Ca/Mg complexes, CaCO3(aq) + CaHCO3+ +
                           MgCO3(aq) + MgHCO3+, as in generate_figures.py. --cations freeion: the
                           free-ion proxies [CaHCO3+]/[HCO3-] and [MgHCO3+]/[HCO3-] (mass action,
                           activity coefficients taken as equal); high if Ca or Mg is above its mean
                  carbon   --carbon total (default): carbon complexes, HCO3- + CO3--, as in
                           generate_figures.py. --carbon carbonate: the carbonate ion CO3-- alone
                  pH
  rank            where the new carbonate sits in each condition, as a volume-weighted
                  percentile within the network (50 = no preference)
  all three high  share of the new carbonate in cells where all three conditions are
                  above the network's volume-weighted mean, and the volume share of
                  those cells
  one alone       share of the new carbonate in cells where that condition is high and
                  at least one other is low
  overlap         share of the supersaturated volume outside the cells where all three
                  conditions are high, and the share of new carbonate in both sets
  absolute        the thresholds above are relative to each network. Two absolute measures,
                  at the start of the interval: the volume where a carbonate precipitates
                  (positive Calcite, Magnesite, Siderite or Dawsonite Rate, meaning the water is
                  supersaturated with it), and the volume with pH above 6, with the share of
                  the new carbonate in each
Writes results/colocation_mapping_<carbon>.csv.
"""
import argparse, csv, glob, os, re, sys
import numpy as np, h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb

INTERVALS = ((1.0, 2.0), (5.0, 7.0), (10.0, 15.0), (20.0, 30.0))
# --stages: consecutive intervals over the whole history, from the snapshot times of the decks
STAGES = ((0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 5.0),
          (5.0, 7.0), (7.0, 10.0), (10.0, 15.0), (15.0, 20.0), (20.0, 30.0), (30.0, 50.0))


def find(g, patterns):
    names = list(g)
    for p in patterns:
        n = next((n for n in names if re.match(p, n)), None)
        if n:
            return n
    return None


def arr(g, n):
    return np.asarray(g[n][:], float).ravel()


CARBON = "total"
CATIONS = "complexes"


def conditions(g):
    """{'Ca', 'Mg', 'carbon', 'pH'} arrays and a description of what was used"""
    ph = find(g, [r"^pH"]); hco3 = find(g, [r"^HCO3- "]); co3 = find(g, [r"^CO3-- "])
    ca = find(g, [r"^Total Ca\+\+", r"^Ca\+\+ "]); mg = find(g, [r"^Total Mg\+\+", r"^Mg\+\+ "])
    if CATIONS == "complexes":
        names = [n for n in (find(g, [r"^CaCO3\(aq\)"]), find(g, [r"^CaHCO3\+"]), find(g, [r"^MgCO3\(aq\)"]), find(g, [r"^MgHCO3\+"])) if n]
        if not names:
            raise RuntimeError("no Ca/Mg complexes written (run with --list)")
        cx = sum(arr(g, n) for n in names)
        C = dict(Ca=cx, Mg=cx)
        how = "cations = Ca/Mg complexes (" + " + ".join(n.split(" [")[0] for n in names) + ")"
        ca = mg = None
    if not (ph and hco3):
        raise RuntimeError("pH or HCO3- not written (run with --list)")
    if CARBON == "carbonate":
        if not co3:
            raise RuntimeError("CO3-- not written (use --carbon total)")
        carbon = arr(g, co3)
    else:
        carbon = arr(g, hco3) + (arr(g, co3) if co3 else 0.0)
    if CATIONS == "complexes":
        pass
    elif ca and mg:
        C = dict(Ca=arr(g, ca), Mg=arr(g, mg)); how = f"cations {ca} and {mg}"
    else:
        cah, mgh = find(g, [r"^CaHCO3\+"]), find(g, [r"^MgHCO3\+"])
        if not (cah and mgh):
            raise RuntimeError("no Ca/Mg species written (run with --list)")
        h = np.where(arr(g, hco3) > 0, arr(g, hco3), np.nan)
        C = dict(Ca=np.nan_to_num(arr(g, cah) / h), Mg=np.nan_to_num(arr(g, mgh) / h))
        how = "cations from [CaHCO3+]/[HCO3-] and [MgHCO3+]/[HCO3-] (free-ion proxies)"
    C.update(carbon=carbon, pH=arr(g, ph))
    return C, how + ("; carbon = CO3-- (carbonate ion)" if CARBON == "carbonate" else "; carbon = carbon complexes, HCO3- + CO3--")


def wpct(x, vol):
    """volume-weighted percentile of each cell; tied values share the middle of their band"""
    o = np.argsort(x, kind="stable"); xs = x[o]; pv = vol[o] / vol.sum(); cum = np.cumsum(pv)
    p = np.empty(len(x)); i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        lo = cum[i] - pv[i]; hi = cum[j]
        p[o[i:j + 1]] = 100 * (lo + hi) / 2
        i = j + 1
    return p


def carbonate_mol(g, vol, mv):
    return sum(np.clip(np.asarray(g[n][:], float).ravel() - cb.SEED, 0, None) * vol / mv[m]
               for m in cb.CARB for n in [next((n for n in g if " VF" in n and n.split(" VF")[0] == m), None)] if n)


def network(d):
    vol, mv = cb.cell_volumes(d), cb.molar_volumes(d)
    snaps = cb.snapshots(d); out = []
    for ta, tb in INTERVALS:
        a = min(snaps, key=lambda s: abs(s[0] - ta)); b = min(snaps, key=lambda s: abs(s[0] - tb))
        if abs(a[0] - ta) > 0.5 or abs(b[0] - tb) > 1.0 or b[0] <= a[0]:
            continue
        with h5py.File(a[1], "r") as fa, h5py.File(b[1], "r") as fb:
            ga, gb = fa[a[2]], fb[b[2]]
            C, how = conditions(ga)
            rates = [np.asarray(ga[n][:], float).ravel() for m in cb.CARB for n in ga if n.startswith(f"{m} Rate")]
            supersat = np.any(np.vstack(rates) > 0, axis=0) if rates else np.zeros(len(vol), bool)
            w = np.clip(carbonate_mol(gb, vol, mv) - carbonate_mol(ga, vol, mv), 0, None)
        if w.sum() <= 0:
            continue
        w = w / w.sum()
        mean = lambda v: (v * vol).sum() / vol.sum()
        pCa, pMg = wpct(C["Ca"], vol), wpct(C["Mg"], vol)
        P = {"cations": np.maximum(pCa, pMg), "carbon": wpct(C["carbon"], vol), "pH": wpct(C["pH"], vol)}
        high = {"cations": (C["Ca"] > mean(C["Ca"])) | (C["Mg"] > mean(C["Mg"])),
                "carbon": C["carbon"] > mean(C["carbon"]), "pH": C["pH"] > mean(C["pH"])}
        allh = high["cations"] & high["carbon"] & high["pH"]
        r = dict(t_start=a[0], t_end=b[0], carbon_dataset=how,
                 **{f"rank_{k}": float((w * P[k]).sum()) for k in P},
                 share_all_three=float(100 * w[allh].sum()), volume_all_three=float(100 * vol[allh].sum() / vol.sum()),
                 volume_supersaturated=float(100 * vol[supersat].sum() / vol.sum()), share_supersaturated=float(100 * w[supersat].sum()),
                 volume_pH_above_6=float(100 * vol[C["pH"] > 6].sum() / vol.sum()), share_pH_above_6=float(100 * w[C["pH"] > 6].sum()),
                 supersat_outside_allhigh_pct=float(100 * vol[supersat & ~allh].sum() / vol[supersat].sum()) if vol[supersat].sum() > 0 else float("nan"),
                 share_supersat_and_allhigh=float(100 * w[supersat & allh].sum()))
        for k in P:
            others = [o for o in P if o != k]
            alone = high[k] & ~(high[others[0]] & high[others[1]])
            r[f"share_{k}_without_others"] = float(100 * w[alone].sum())
        out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs_gravityoff"); ap.add_argument("--prefix", default="A_feedback__")
    ap.add_argument("--list", default=None, help="print the dataset names of one run's first snapshot and stop")
    ap.add_argument("--stages", action="store_true", help="consecutive intervals over the whole history (0.01 to 50 y)")
    ap.add_argument("--carbon", choices=["total", "carbonate"], default="total")
    ap.add_argument("--cations", choices=["complexes", "freeion"], default="complexes")
    a = ap.parse_args()
    global INTERVALS
    if a.stages:
        INTERVALS = STAGES
    global CARBON, CATIONS; CARBON, CATIONS = a.carbon, a.cations
    if a.list:
        t, p, k = cb.snapshots(a.list)[0]
        with h5py.File(p, "r") as f:
            try:
                print("\n".join(sorted(f[k])))
            except BrokenPipeError:
                pass
        return
    root = a.root if os.path.isabs(a.root) else os.path.join(REV, a.root)
    rows = []
    for d in sorted(glob.glob(os.path.join(root, a.prefix + "p32_*"))):
        try:
            for r in network(d):
                rows.append(dict(network=os.path.basename(d), **r))
        except Exception as e:
            print(f"  skipped {os.path.basename(d)}: {e}")
    if not rows:
        sys.exit("no results")
    print(f"measured as: {rows[0]['carbon_dataset']}")
    if CATIONS == "freeion":
        print("rank for cations: the higher of the Ca and Mg ranks")
    print(f"{'interval':>10} {'n':>3}  {'rank: cations':>13} {'carbon':>7} {'pH':>6}   {'new carbonate in cells with all three high':>42}"
          f"   {'with one high, another low: cations / carbon / pH':>50}")
    for ta, tb in INTERVALS:
        rs = [r for r in rows if abs(r["t_start"] - ta) <= 1e-4 + 0.005 * ta]
        if not rs:
            continue
        m = lambda k: np.median([r[k] for r in rs])
        print(f"{ta:>4g}-{tb:<4g}y {len(rs):>3}  {m('rank_cations'):13.1f} {m('rank_carbon'):7.1f} {m('rank_pH'):6.1f}   "
              f"{m('share_all_three'):14.1f}% of it, in {m('volume_all_three'):5.1f}% of the volume      "
              f"{m('share_cations_without_others'):6.1f}% / {m('share_carbon_without_others'):5.1f}% / {m('share_pH_without_others'):5.1f}%")
    print(f"\noverlap of the supersaturated volume with the cells where all three conditions are high (medians of networks)")
    print(f"{'interval':>10}   {'supersaturated volume outside those cells':>42}   {'new carbonate in both':>22}")
    for ta, tb in INTERVALS:
        rs = [r for r in rows if abs(r["t_start"] - ta) <= 1e-4 + 0.005 * ta]
        if rs:
            m = lambda k: np.nanmedian([r[k] for r in rs]); mx = lambda k: np.nanmax([r[k] for r in rs])
            print(f"{ta:>4g}-{tb:<4g}y   {m('supersat_outside_allhigh_pct'):8.2f}% (max {mx('supersat_outside_allhigh_pct'):6.2f}%)                  {m('share_supersat_and_allhigh'):8.1f}%")
    print(f"\nabsolute measures at the start of each interval (medians of networks)")
    print(f"{'interval':>10}   {'supersaturated with a carbonate':>32}   {'pH above 6':>28}")
    for ta, tb in INTERVALS:
        rs = [r for r in rows if abs(r["t_start"] - ta) <= 1e-4 + 0.005 * ta]
        if not rs:
            continue
        m = lambda k: np.median([r[k] for r in rs])
        print(f"{ta:>4g}-{tb:<4g}y   {m('volume_supersaturated'):6.2f}% of the volume, {m('share_supersaturated'):5.1f}% of new carbonate"
              f"   {m('volume_pH_above_6'):6.2f}% of the volume, {m('share_pH_above_6'):5.1f}% of new carbonate")
    keys = list(rows[0])
    p = os.path.join(REV, "results", f"colocation_mapping_{CATIONS}_{CARBON}.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"written: results/{os.path.basename(p)}")


if __name__ == "__main__":
    main()
