#!/usr/bin/env python3
"""
Basic results for a block: one row per network, one per intensity group, as CSV.

    python3 src/summarize_block.py --prefix A_ --name A_fixed
    python3 src/summarize_block.py --prefix A_feedback__ --name A_coupled

Per network (last snapshot, normally 50 y), all volume-weighted with the mesh's
cell volumes:
  co2_kg, efficiency_pct      CO2 mineralized as carbonate, and / CO2 injected
                              (as carbon_budget.py: moles = net VF x volume / molar volume)
  overfilled_pct              share of that CO2 in cells whose carbonate VF exceeds 1
  co2_valid_kg                CO2 in cells whose carbonate VF is at most 1 (lenient)
  co2_strict_kg               CO2 in cells whose minerals together (primary, carbonate, clay)
                              are at most 1, i.e. non-negative porosity (the physical limit)
  calcite/magnesite/siderite/dawsonite_pct   shares of the carbonate, by moles
  precip_cells, precip_vol_pct               cells with net carbonate > 0, and their share of the volume
  <primary>_dissolved_pct     1 - final / initial mineral volume over the domain, all six primaries
  <primary>_dissolved_mol     moles dissolved (dissolved volume / molar volume from hanford.dat),
                              comparable with the carbonate formed
  pH_mean, pH_median          at the last snapshot, volume-weighted mean and median
Per intensity: n, mean / median / CV of co2_kg and efficiency, median shares and
dissolution; overall: Spearman(P32 factor, co2_kg).
Writes results/<name>_networks.csv and results/<name>_by_intensity.csv.
"""
import argparse, csv, glob, os, re, sys
import numpy as np, h5py
from scipy.stats import spearmanr

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb

CARB = cb.CARB
PRIM = ("Anorthite", "Albite", "Diopside", "Forsterite", "Fayalite", "Enstatite")


def tval(k):
    return float(k.split("Time")[1].split()[0])


def vf(g, mineral):
    for n in g:
        if " VF" in n and n.split(" VF")[0] == mineral:
            return np.asarray(g[n][:], float).ravel()
    return None


def primary_molar_volumes(d):
    """m3/mol for the primary minerals, from the run's hanford.dat (NaN where missing)."""
    db = os.path.join(d, "hanford.dat")
    if not os.path.isfile(db):
        db = os.path.join(os.path.dirname(REV), "database", "hanford.dat")
    mv = {}
    for line in open(db, errors="replace"):
        t = line.split()
        if len(t) > 1 and t[0].startswith("'") and t[0].strip("'") in PRIM and t[0].strip("'") not in mv:
            try:
                mv[t[0].strip("'")] = float(t[1].replace("d", "e")) * 1e-6
            except ValueError:
                pass
    return mv


def wmedian(x, w):
    o = np.argsort(x); c = np.cumsum(w[o])
    return float(x[o][np.searchsorted(c, 0.5 * c[-1])])


def network(d):
    deck = next(iter(sorted(glob.glob(os.path.join(d, "*.in")))), None)
    snaps = cb.snapshots(d) if deck else []          # all parts of a continued run, in time order
    if not snaps:
        return None
    vol, mv = cb.cell_volumes(d), cb.molar_volumes(d)
    pmv = primary_molar_volumes(d)
    (_, p0, k0), (t, p1, k1) = snaps[0], snaps[-1]
    with h5py.File(p0, "r") as f0, h5py.File(p1, "r") as f:
        g0, g1 = f0[k0], f[k1]
        net = {m: np.clip(vf(g1, m) - cb.SEED, 0, None) for m in CARB if vf(g1, m) is not None}
        mol = {m: net[m] * vol / mv[m] for m in net}
        prim, pmol = {}, {}
        for m in PRIM:
            a0, a1 = vf(g0, m), vf(g1, m)
            if a0 is not None and a1 is not None and (a0 * vol).sum() > 0:
                prim[m] = 100 * (1 - (a1 * vol).sum() / (a0 * vol).sum())
                if m in pmv:
                    pmol[m] = float(((a0 - a1) * vol).sum()) / pmv[m]
        ph = next((np.asarray(g1[n][:], float).ravel() for n in g1 if n.strip().lower().startswith("ph")), None)
        allvf = sum(np.asarray(g1[n][:], float).ravel() for n in g1 if " VF" in n)
    tot_mol_cell = sum(mol.values()); tot = float(tot_mol_cell.sum())
    carb_vf = sum(net.values())
    inj = cb.injected_mol(deck, t)
    row = dict(network=os.path.basename(d), years=t, co2_kg=tot * cb.M_CO2,
               efficiency_pct=100 * tot / inj if inj else np.nan,
               overfilled_pct=100 * float(tot_mol_cell[carb_vf > 1].sum()) / tot if tot else 0.0,
               co2_valid_kg=float(tot_mol_cell[carb_vf <= 1].sum()) * cb.M_CO2,
               co2_strict_kg=float(tot_mol_cell[allvf <= 1].sum()) * cb.M_CO2,
               precip_cells=int((carb_vf > 0).sum()), precip_vol_pct=100 * float(vol[carb_vf > 0].sum() / vol.sum()))
    for m in CARB:
        row[f"{m.lower()}_pct"] = 100 * float(mol[m].sum()) / tot if tot and m in mol else 0.0
    for m in PRIM:
        row[f"{m.lower()}_dissolved_pct"] = prim.get(m, np.nan)
        row[f"{m.lower()}_dissolved_mol"] = pmol.get(m, np.nan)
    row["pH_mean"] = float((ph * vol).sum() / vol.sum()) if ph is not None else np.nan
    row["pH_median"] = wmedian(ph, vol) if ph is not None else np.nan
    m = re.search(r"p32_(\d{3})_s(\d+)", row["network"])
    row["p32_factor"] = int(m.group(1)) / 100 if m else np.nan
    row["seed"] = int(m.group(2)) if m else -1
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--root", default="runs_gravityoff")
    a = ap.parse_args()
    root = a.root if os.path.isabs(a.root) else os.path.join(REV, a.root)
    dirs = sorted(d for d in glob.glob(os.path.join(root, a.prefix + "p32_*")) if os.path.isdir(d))
    rows, skipped = [], []
    for d in dirs:
        try:
            r = network(d)
        except Exception as e:
            skipped.append(f"{os.path.basename(d)} ({e})"); continue
        if r is None or r["years"] < 49.99:
            skipped.append(os.path.basename(d) + ("" if r is None else f" (at {r['years']:.3g} y)")); continue
        rows.append(r)
    if not rows:
        sys.exit("no finished runs")
    os.makedirs(os.path.join(REV, "results"), exist_ok=True)
    out1 = os.path.join(REV, "results", f"{a.name}_networks.csv")
    keys = ["network", "p32_factor", "seed", "years", "co2_kg", "efficiency_pct", "overfilled_pct",
            "co2_valid_kg", "co2_strict_kg",
            "calcite_pct", "magnesite_pct", "siderite_pct", "dawsonite_pct", "precip_cells", "precip_vol_pct"] + \
           [f"{m.lower()}_dissolved_pct" for m in PRIM] + [f"{m.lower()}_dissolved_mol" for m in PRIM] + ["pH_mean", "pH_median"]
    with open(out1, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in sorted(rows, key=lambda r: (r["p32_factor"], r["seed"])):
            w.writerow({k: (f"{r[k]:.6g}" if isinstance(r[k], float) else r[k]) for k in keys})
    groups = sorted({r["p32_factor"] for r in rows})
    gk = ["p32_factor", "n", "co2_kg_mean", "co2_kg_median", "co2_kg_cv_pct", "co2_valid_kg_median",
          "co2_strict_kg_median", "co2_strict_kg_cv_pct", "efficiency_pct_mean",
          "efficiency_pct_median", "efficiency_cv_pct", "overfilled_pct_median", "calcite_pct_median",
          "magnesite_pct_median", "siderite_pct_median", "dawsonite_pct_median"] + \
         [f"{m.lower()}_dissolved_pct_median" for m in PRIM] + ["pH_mean_median"]
    grows = []
    for p in groups:
        rs = [r for r in rows if r["p32_factor"] == p]
        c = np.array([r["co2_kg"] for r in rs]); e = np.array([r["efficiency_pct"] for r in rs])
        med = lambda k: float(np.nanmedian([r[k] for r in rs]))
        cv = lambda x: 100 * x.std(ddof=1) / x.mean() if len(x) > 1 and x.mean() else np.nan
        cs = np.array([r["co2_strict_kg"] for r in rs])
        grows.append(dict(p32_factor=p, n=len(rs), co2_kg_mean=c.mean(), co2_kg_median=np.median(c), co2_kg_cv_pct=cv(c),
                          co2_valid_kg_median=med("co2_valid_kg"), co2_strict_kg_median=np.median(cs),
                          co2_strict_kg_cv_pct=cv(cs),
                          efficiency_pct_mean=e.mean(), efficiency_pct_median=np.median(e), efficiency_cv_pct=cv(e),
                          overfilled_pct_median=med("overfilled_pct"), calcite_pct_median=med("calcite_pct"),
                          magnesite_pct_median=med("magnesite_pct"), siderite_pct_median=med("siderite_pct"),
                          dawsonite_pct_median=med("dawsonite_pct"), pH_mean_median=med("pH_mean"),
                          **{f"{m.lower()}_dissolved_pct_median": med(f"{m.lower()}_dissolved_pct") for m in PRIM}))
    out2 = os.path.join(REV, "results", f"{a.name}_by_intensity.csv")
    with open(out2, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=gk); w.writeheader()
        for g in grows:
            w.writerow({k: (f"{g[k]:.6g}" if isinstance(g[k], float) else g[k]) for k in gk})
    print(f"{a.name}: {len(rows)} networks at 50 y" + (f"; left out: {', '.join(skipped)}" if skipped else ""))
    print(f"{'P32':>5}{'n':>4}{'CO2 kg median':>15}{'CV %':>7}{'strict kg':>11}{'CV %':>7}{'eff % median':>14}{'overfilled %':>14}"
          f"{'cal/mag/sid/daw %':>22}{'anor diss %':>13}{'pH':>6}")
    for g in grows:
        mix = f"{g['calcite_pct_median']:.0f}/{g['magnesite_pct_median']:.0f}/{g['siderite_pct_median']:.0f}/{g['dawsonite_pct_median']:.0f}"
        print(f"{g['p32_factor']:>5.2f}{g['n']:>4}{g['co2_kg_median']:>15.4g}{g['co2_kg_cv_pct']:>7.0f}"
              f"{g['co2_strict_kg_median']:>11.4g}{g['co2_strict_kg_cv_pct']:>7.0f}"
              f"{g['efficiency_pct_median']:>14.3g}{g['overfilled_pct_median']:>14.1f}{mix:>22}"
              f"{g['anorthite_dissolved_pct_median']:>13.1f}{g['pH_mean_median']:>6.2f}")
    print(f"\nprimary minerals dissolved by 50 y (median % of the initial volume)")
    print(f"{'P32':>5}" + "".join(f"{m:>12}" for m in PRIM))
    for g in grows:
        print(f"{g['p32_factor']:>5.2f}" + "".join(f"{g[m.lower() + '_dissolved_pct_median']:>12.1f}" for m in PRIM))
    c = np.array([r["co2_kg"] for r in rows]); p = np.array([r["p32_factor"] for r in rows])
    rho = spearmanr(p, c)
    allc = np.array([r["co2_kg"] for r in rows])
    print(f"all: CO2 median {np.median(allc):.4g} kg, CV {100 * allc.std(ddof=1) / allc.mean():.0f}%, "
          f"Spearman(P32, CO2) = {rho.correlation:.3f} (p = {rho.pvalue:.2g})")
    st = np.array([r["co2_strict_kg"] for r in rows]); rs2 = spearmanr(p, st)
    print(f"strict (non-negative porosity only): CO2 median {np.median(st):.4g} kg, "
          f"CV {100 * st.std(ddof=1) / st.mean():.0f}%, Spearman(P32, CO2) = {rs2.correlation:.3f} (p = {rs2.pvalue:.2g})")
    print(f"written: {os.path.relpath(out1, REV)}, {os.path.relpath(out2, REV)}")


if __name__ == "__main__":
    main()
