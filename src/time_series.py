#!/usr/bin/env python3
"""
Trace pH, carbonate and dissolution through time, for every network of a block.

    python3 src/time_series.py --prefix A_ --name A_fixed
    python3 src/time_series.py --prefix A_feedback__ --name A_coupled

At every snapshot, volume-weighted with the mesh's cell volumes:
  pH_mean, pH_p10, pH_p90       mean pH and its 10th / 90th percentiles over the volume
  vol_pH_gt6_pct                share of the volume above pH 6 (where carbonate can be stable)
  co2_kg, co2_strict_kg         CO2 in carbonate, and in cells with non-negative porosity only
  precip_cells, overfilled_pct  cells with net carbonate > 0; share of the CO2 in cells with carbonate VF > 1
  <mineral>_dissolved_pct       anorthite, diopside, forsterite
  <carbonate>_pct, <carbonate>_kg   each carbonate's share of the carbonate (moles), and its own
                                CO2 in kg: a decrease over time indicates that the mineral re-dissolves
Writes results/<name>_timeseries.csv (one row per network and time),
results/<name>_timeseries_by_intensity.csv (median, 25th and 75th percentiles
per intensity and time) and results/<name>_timeseries.png (if matplotlib is present).
"""
import argparse, csv, glob, os, re, sys
import numpy as np, h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb

CARB, FAST = cb.CARB, ("Anorthite", "Diopside", "Forsterite")
Q = ("pH_mean", "pH_p10", "pH_p90", "vol_pH_gt6_pct", "co2_kg", "co2_strict_kg", "precip_cells",
     "overfilled_pct") + tuple(f"{m.lower()}_dissolved_pct" for m in FAST) \
    + tuple(f"{m.lower()}_pct" for m in CARB) + tuple(f"{m.lower()}_kg" for m in CARB)


def tval(k):
    return float(k.split("Time")[1].split()[0])


def get(g, mineral):
    for n in g:
        if " VF" in n and n.split(" VF")[0] == mineral:
            return np.asarray(g[n][:], float).ravel()
    return None


def wq(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]); return float(x[o][min(len(x) - 1, np.searchsorted(c, q * c[-1]))])


def trace(d):
    vol, mv = cb.cell_volumes(d), cb.molar_volumes(d)
    out = []
    snaps = cb.snapshots(d)                          # all parts of a continued run, in time order
    files = {p: h5py.File(p, "r") for p in {p for _, p, _ in snaps}}
    try:
        a0 = {m: get(files[snaps[0][1]][snaps[0][2]], m) for m in FAST}
        for t, p, k in snaps:
            g = files[p][k]
            net = {m: np.clip(get(g, m) - cb.SEED, 0, None) for m in CARB if get(g, m) is not None}
            mol = sum(net[m] * vol / mv[m] for m in net); carb = sum(net.values())
            allvf = sum(np.asarray(g[n][:], float).ravel() for n in g if " VF" in n)
            ph = next((np.asarray(g[n][:], float).ravel() for n in g if n.strip().lower().startswith("ph")), None)
            tot = float(mol.sum())
            row = dict(t=t, co2_kg=tot * cb.M_CO2, co2_strict_kg=float(mol[allvf <= 1].sum()) * cb.M_CO2,
                       precip_cells=int((carb > 0).sum()),
                       overfilled_pct=100 * float(mol[carb > 1].sum()) / tot if tot else 0.0)
            if ph is not None:
                row.update(pH_mean=float((ph * vol).sum() / vol.sum()), pH_p10=wq(ph, vol, 0.10), pH_p90=wq(ph, vol, 0.90),
                           vol_pH_gt6_pct=100 * float(vol[ph > 6].sum() / vol.sum()))
            for m in CARB:
                mm = float((net[m] * vol / mv[m]).sum()) if m in net else 0.0
                row[f"{m.lower()}_kg"] = mm * cb.M_CO2
                row[f"{m.lower()}_pct"] = 100 * mm / tot if tot else np.nan
            for m in FAST:
                a1 = get(g, m)
                row[f"{m.lower()}_dissolved_pct"] = (100 * (1 - (a1 * vol).sum() / (a0[m] * vol).sum())
                                                    if a0[m] is not None and a1 is not None and (a0[m] * vol).sum() else np.nan)
            out.append(row)
    finally:
        for fh in files.values():
            fh.close()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True); ap.add_argument("--name", required=True)
    ap.add_argument("--root", default="runs_gravityoff")
    a = ap.parse_args()
    root = a.root if os.path.isabs(a.root) else os.path.join(REV, a.root)
    dirs = sorted(d for d in glob.glob(os.path.join(root, a.prefix + "p32_*"))
                  if os.path.isfile(os.path.join(d, "pflotran_co2.h5")))
    rows = []
    for d in dirs:
        m = re.search(r"p32_(\d{3})_s(\d+)", os.path.basename(d))
        try:
            tr = trace(d)
        except Exception as e:
            print(f"  skipped {os.path.basename(d)}: {e}"); continue
        if not tr or tr[-1]["t"] < 49.99:
            print(f"  skipped {os.path.basename(d)}: not at 50 y"); continue
        for r in tr:
            r.update(network=os.path.basename(d), p32_factor=int(m.group(1)) / 100 if m else np.nan,
                     seed=int(m.group(2)) if m else -1)
            rows.append(r)
    if not rows:
        sys.exit("no finished runs")
    os.makedirs(os.path.join(REV, "results"), exist_ok=True)
    keys = ["network", "p32_factor", "seed", "t"] + list(Q)
    p1 = os.path.join(REV, "results", f"{a.name}_timeseries.csv")
    with open(p1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r[k]:.6g}" if isinstance(r.get(k), float) else r.get(k, "")) for k in keys})
    groups, times = sorted({r["p32_factor"] for r in rows}), sorted({round(r["t"], 6) for r in rows})
    summ = []
    for p in groups:
        for t in times:
            rs = [r for r in rows if r["p32_factor"] == p and round(r["t"], 6) == t]
            if not rs:
                continue
            s = dict(p32_factor=p, t=t, n=len(rs))
            for q in Q:
                v = np.array([r.get(q, np.nan) for r in rs], float)
                s[f"{q}_median"], s[f"{q}_p25"], s[f"{q}_p75"] = (np.nanmedian(v), np.nanpercentile(v, 25), np.nanpercentile(v, 75)) \
                    if np.isfinite(v).any() else (np.nan,) * 3
            summ.append(s)
    p2 = os.path.join(REV, "results", f"{a.name}_timeseries_by_intensity.csv")
    sk = ["p32_factor", "t", "n"] + [f"{q}_{x}" for q in Q for x in ("median", "p25", "p75")]
    with open(p2, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sk); w.writeheader()
        for s in summ:
            w.writerow({k: (f"{s[k]:.6g}" if isinstance(s[k], float) else s[k]) for k in sk})
    nnet = len({r["network"] for r in rows})
    print(f"{a.name}: {nnet} networks, {len(times)} snapshots")
    show = [t for t in times if t in (0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0)]
    for label, q, fmt in (("mean pH", "pH_mean", "{:6.2f}"), ("volume above pH 6 (%)", "vol_pH_gt6_pct", "{:6.1f}"),
                          ("CO2, strict (kg)", "co2_strict_kg", "{:6.3g}"), ("overfilled share (%)", "overfilled_pct", "{:6.1f}"),
                          ("calcite share of the carbonate (%)", "calcite_pct", "{:6.1f}"),
                          ("magnesite share of the carbonate (%)", "magnesite_pct", "{:6.1f}"),
                          ("CO2 in magnesite (kg); a decrease over time indicates re-dissolution", "magnesite_kg", "{:6.3g}"),
                          ("CO2 in calcite (kg)", "calcite_kg", "{:6.3g}")):
        print(f"\n{label}, median over networks\n{'P32':>5}" + "".join(f"{t:>7g}y" for t in show))
        for p in groups:
            line = f"{p:>5.2f}"
            for t in show:
                s = next((s for s in summ if s["p32_factor"] == p and s["t"] == t), None)
                line += " " + (fmt.format(s[f"{q}_median"]) if s and np.isfinite(s[f"{q}_median"]) else "     -") + " "
            print(line)
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
        panels = (("pH_mean", "mean pH"), ("vol_pH_gt6_pct", "volume above pH 6 (%)"),
                  ("co2_strict_kg", "CO$_2$ in carbonate, non-negative porosity (kg)"), ("forsterite_dissolved_pct", "dissolved (%)"))
        cmap = plt.get_cmap("viridis")
        for i, (q, lab) in enumerate(panels):
            A = ax.flat[i]
            for j, p in enumerate(groups):
                ss = [s for s in summ if s["p32_factor"] == p and s["t"] > 0]
                tt = [s["t"] for s in ss]; col = cmap(j / max(1, len(groups) - 1))
                if q == "forsterite_dissolved_pct":
                    for mn, ls in (("forsterite", "-"), ("diopside", "--"), ("anorthite", ":")):
                        A.plot(tt, [s[f"{mn}_dissolved_pct_median"] for s in ss], ls, color=col, lw=1.2,
                               label=f"{mn}" if j == len(groups) - 1 else None)
                else:
                    A.plot(tt, [s[f"{q}_median"] for s in ss], color=col, lw=1.5, label=f"x{p:.2f}")
                    A.fill_between(tt, [s[f"{q}_p25"] for s in ss], [s[f"{q}_p75"] for s in ss], color=col, alpha=0.15)
            A.set_xscale("log"); A.set_xlabel("time (y)"); A.set_ylabel(lab)
            A.legend(fontsize=7, frameon=False)
        fig.suptitle(f"{a.name}: median over networks (bands: 25th to 75th percentile)", fontsize=10)
        fig.tight_layout()
        p3 = os.path.join(REV, "results", f"{a.name}_timeseries.png"); fig.savefig(p3, dpi=150)
        print(f"\nwritten: results/{os.path.basename(p1)}, results/{os.path.basename(p2)}, results/{os.path.basename(p3)}")
    except ImportError:
        print(f"\nwritten: results/{os.path.basename(p1)}, results/{os.path.basename(p2)} (no matplotlib: no figure)")


if __name__ == "__main__":
    main()
