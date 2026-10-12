#!/usr/bin/env python3
"""
Dissolution location, capture of released cations in carbonate, secondary silicates, and the
minimum carbonate per set.

    python3 src/dissolution_capture_analysis.py --root runs_gravityoff

(1) Where does dissolution happen, compared with carbonate? For each coupled Block A
    network and each time (2, 10, 20, 50 y), the dissolved volume of forsterite,
    diopside and anorthite in each cell, and its position in the water-age range
    (volume-weighted percentile, 50 = no preference, low = young water).
    The carbonate position is given at the same times.
(2) What share of the released calcium, magnesium and iron ends up in carbonate?
    Released: anorthite + diopside (Ca), 2 x forsterite + diopside + enstatite (Mg),
    2 x fayalite (Fe), in moles, from results/<set>_networks.csv. Carbonate: moles of CO2.
(3) How much kaolinite and chalcedony form? Net moles at 50 y from PFLOTRAN's mass balance,
    and their volume as a share of the initial pore volume.
(4) Does every network form carbonate? The lowest carbonate per set.
Writes results/dissolution_capture_analysis.csv.
"""
import argparse, csv, glob, os, re, statistics as st, sys
import numpy as np, h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb
import flowfield

FAST = ("Forsterite", "Diopside", "Anorthite")
SILIC = ("Kaolinite", "Chalcedony")
TIMES = (2.0, 10.0, 20.0, 50.0)


def vf(g, m):
    n = next((n for n in g if " VF" in n and n.split(" VF")[0] == m), None)
    return np.asarray(g[n][:], float).ravel() if n else None


def pct_weighted(age, vol, w):
    o = np.argsort(age, kind="stable"); pv = vol[o] / vol.sum()
    pct = np.empty(len(vol)); pct[o] = 100 * (np.cumsum(pv) - pv / 2)
    return float((w * pct).sum() / w.sum()) if w.sum() > 0 else float("nan")


def zones(d):
    sol = flowfield.solve(d, verbose=False)
    vol = np.asarray(sol["vol"], float)
    age = np.where(np.asarray(sol["stagnant"], bool) | ~np.isfinite(sol["age"]), np.inf, np.asarray(sol["age"], float))
    snaps = cb.snapshots(d); mv = cb.molar_volumes(d); out = {}
    t0, p0, k0 = snaps[0]
    with h5py.File(p0, "r") as f:
        a0 = {m: vf(f[k0], m) for m in FAST}
    for T in TIMES:
        t, p, k = min(snaps, key=lambda s: abs(s[0] - T))
        if abs(t - T) > 0.5:
            continue
        with h5py.File(p, "r") as f:
            g = f[k]
            diss = sum(np.clip(a0[m] - vf(g, m), 0, None) * vol for m in FAST if a0[m] is not None and vf(g, m) is not None)
            carb = sum(np.clip(vf(g, m) - cb.SEED, 0, None) * vol / mv[m] for m in cb.CARB if vf(g, m) is not None)
        out[T] = (pct_weighted(age, vol, diss), pct_weighted(age, vol, carb))
    return out


def silicates(d):
    ps = cb._parts(d, "pflotran_co2-mas", ".dat")
    if not ps:
        return None
    (c0, r0), (c1, r1) = cb._mas_rows(ps[0]), cb._mas_rows(ps[-1])
    res = {}
    for m in SILIC:
        i0 = next((i for i, c in enumerate(c0) if f" {m} Total Mass" in c), None)
        i1 = next((i for i, c in enumerate(c1) if f" {m} Total Mass" in c), None)
        res[m] = float(r1[-1][i1]) - float(r0[0][i0]) if i0 is not None and i1 is not None else float("nan")
    iw = next((i for i, c in enumerate(c0) if "Global Water Mass" in c), None)
    res["pore_m3"] = float(r0[0][iw]) / 1000.0 if iw is not None else float("nan")
    return res


def silicate_molar_volumes(d):
    db = os.path.join(d, "hanford.dat"); mv = {}
    for line in open(db, errors="replace"):
        t = line.split()
        if len(t) > 1 and t[0].strip("'") in SILIC and t[0].strip("'") not in mv:
            try:
                mv[t[0].strip("'")] = float(t[1].replace("d", "e")) * 1e-6
            except ValueError:
                pass
    return mv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs_gravityoff")
    ap.add_argument("--zones-prefix", default="A_feedback__")
    a = ap.parse_args()
    root = a.root if os.path.isabs(a.root) else os.path.join(REV, a.root)
    rows = []

    print("(1) Where dissolution and carbonate sit in the water-age range (50 = no preference, low = young water)")
    Z = {T: ([], []) for T in TIMES}
    for d in sorted(glob.glob(os.path.join(root, a.zones_prefix + "p32_*"))):
        try:
            z = zones(d)
        except Exception as e:
            print(f"  skipped {os.path.basename(d)}: {e}"); continue
        for T, (dp, cp) in z.items():
            Z[T][0].append(dp); Z[T][1].append(cp)
            rows.append(dict(check="zones", network=os.path.basename(d), t=T, dissolution_age_pct=dp, carbonate_age_pct=cp))
    for T in TIMES:
        if Z[T][0]:
            print(f"  {T:4g} y ({len(Z[T][0])} networks): dissolution at age percentile {np.nanmedian(Z[T][0]):5.1f} "
                  f"({np.nanmin(Z[T][0]):.1f} to {np.nanmax(Z[T][0]):.1f}),  carbonate at {np.nanmedian(Z[T][1]):5.1f}")

    print("\n(2) Share of the released Ca, Mg and Fe that ends up in carbonate")
    for name in ("A_coupled", "B_coupled", "E_coupled"):
        p = os.path.join(REV, "results", f"{name}_networks.csv")
        if not os.path.isfile(p):
            print(f"  {name}: no results file"); continue
        q = []
        for r in csv.DictReader(open(p)):
            g = lambda k: float(r[k]) if r.get(k) not in (None, "", "nan") else 0.0
            rel = g("anorthite_dissolved_mol") + 2 * g("diopside_dissolved_mol") + 2 * g("forsterite_dissolved_mol") \
                + g("enstatite_dissolved_mol") + 2 * g("fayalite_dissolved_mol")
            if rel > 0:
                q.append(100 * g("co2_kg") / cb.M_CO2 / rel)
                rows.append(dict(check="capture", network=r["network"], capture_pct=q[-1]))
        if q:
            print(f"  {name}: median {st.median(q):.3g}% (range {min(q):.3g} to {max(q):.3g}%), {len(q)} networks")

    print("\n(3) Kaolinite and chalcedony at 50 y (net, from the mass balance)")
    for pre in ("A_feedback__", "E_feedback__"):
        k, c, fr = [], [], []
        for d in sorted(glob.glob(os.path.join(root, pre + "p32_*"))):
            s = silicates(d)
            if not s:
                continue
            mv = silicate_molar_volumes(d)
            vol = sum(s[m] * mv.get(m, np.nan) for m in SILIC)
            k.append(s["Kaolinite"]); c.append(s["Chalcedony"]); fr.append(100 * vol / s["pore_m3"])
            rows.append(dict(check="silicates", network=os.path.basename(d), kaolinite_mol=s["Kaolinite"],
                             chalcedony_mol=s["Chalcedony"], silicate_pct_of_pore=fr[-1]))
        if k:
            print(f"  {pre[:-2]}: kaolinite median {st.median(k):.4g} mol, chalcedony median {st.median(c):.4g} mol; "
                  f"together {st.median(fr):.3g}% of the initial pore volume (range {min(fr):.3g} to {max(fr):.3g}%), {len(k)} networks")

    print("\n(4) Lowest carbonate per set (every network forms carbonate?)")
    for name in ("A_coupled", "B_coupled", "E_coupled", "D30_coupled", "D40_coupled"):
        p = os.path.join(REV, "results", f"{name}_networks.csv")
        if os.path.isfile(p):
            v = sorted((float(r["co2_kg"]), r["network"]) for r in csv.DictReader(open(p)))
            print(f"  {name}: lowest {v[0][0]:.3g} kg CO2 ({v[0][1]}), median {st.median(x for x, _ in v):.3g} kg, "
                  f"networks below 1% of the median: {sum(x < 0.01 * st.median(y for y, _ in v) for x, _ in v)}")

    if rows:
        keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("check", "network"), k))
        p = os.path.join(REV, "results", "dissolution_capture_analysis.csv")
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
        print(f"\nwritten: results/{os.path.basename(p)}")


if __name__ == "__main__":
    main()
