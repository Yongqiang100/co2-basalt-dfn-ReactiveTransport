#!/usr/bin/env python3
"""
Is dissolution uniform? Map of the Ca and Mg released by dissolution between two times, for the
networks of Figure 9. The map shows the dissolution alone; the precipitation is not drawn.

    python3 src/dissolution_map.py --root runs_gravityoff \
        --runs A_feedback__p32_100_s1397 A_feedback__p32_100_s1063 --t0 5 --t1 7

Per network it prints
  CV of release      coefficient of variation of the release per unit volume over the cells
  top-10% share      share of the total release in the 10% of the network volume that releases
                     most (10% = uniform dissolution)
  prec/rest          release per unit volume in the precipitating cells divided by that in
                     the other cells (1 = the precipitating cells are not dissolution hot spots)
and writes figures/dissolution_map.png (x-z view, as Figure 9, dissolution only).
Release: net loss of anorthite (1 Ca), diopside (1 Ca, 1 Mg), forsterite (2 Mg), from the
volume fractions, the cell volumes and the molar volumes of local_cation_balance.py.
"""
import argparse, os, sys
import numpy as np, h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
sys.path.insert(0, os.path.join(REV, "src"))
import carbon_budget as cb
from local_cation_balance import MOLAR_VOL, RELEASE, vf, at


def coords(run_dir):
    u = os.path.join(run_dir, "full_mesh.uge"); n = int(open(u).readline().split()[1])
    a = np.loadtxt(u, skiprows=1, max_rows=n)
    return a[:, 1:4], a[:, 4]


def network(run_dir, t0, t1, min_inc):
    snaps = cb.snapshots(run_dir); xyz, vol = coords(run_dir); n = len(vol)
    (_, pa, ka), (_, pb, kb) = at(snaps, t0), at(snaps, t1)
    with h5py.File(pa, "r") as fa, h5py.File(pb, "r") as fb:
        ga, gb = fa[ka], fb[kb]
        rel = np.zeros(n)
        for m, (sca, smg) in RELEASE.items():
            rel += (sca + smg) * (vf(ga, m, n) - vf(gb, m, n)) * vol / MOLAR_VOL[m]
        carb = sum((vf(gb, m, n) - vf(ga, m, n)) for m in ("Calcite", "Magnesite"))
    rel = np.clip(rel, 0, None); dens = rel / vol; prec = carb > min_inc
    order = np.argsort(dens)[::-1]; cumv = np.cumsum(vol[order]) / vol.sum()
    top = rel[order][cumv <= 0.10].sum() / rel.sum() if rel.sum() > 0 else float("nan")
    w = vol / vol.sum(); mu = (w * dens).sum(); cv = np.sqrt((w * (dens - mu) ** 2).sum()) / mu if mu > 0 else float("nan")
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = (rel[prec].sum() / vol[prec].sum()) / (rel[~prec].sum() / vol[~prec].sum()) if prec.any() and (~prec).any() else float("nan")
    return dict(xyz=xyz, dens=dens, prec=prec, cv=100 * cv, top=100 * top, ratio=ratio, nprec=int(prec.sum()), n=n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs_gravityoff")
    ap.add_argument("--runs", nargs="+", default=["A_feedback__p32_100_s1397", "A_feedback__p32_100_s1063"])
    ap.add_argument("--t0", type=float, default=5.0); ap.add_argument("--t1", type=float, default=7.0)
    ap.add_argument("--min-increase", type=float, default=1e-9)
    ap.add_argument("--no-figure", action="store_true")
    a = ap.parse_args()
    root = a.root if os.path.isabs(a.root) else os.path.join(REV, a.root)
    res = {}
    print(f"dissolution between {a.t0:g} and {a.t1:g} years (release of Ca + Mg)")
    for r in a.runs:
        d = network(os.path.join(root, r), a.t0, a.t1, a.min_increase); res[r] = d
        print(f"  {r}: CV of release {d['cv']:.0f}%, top-10% share {d['top']:.1f}%, "
              f"prec/rest {d['ratio']:.2f} ({d['nprec']} of {d['n']} cells precipitate)")
    if a.no_figure:
        return
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    fig, axes = plt.subplots(1, len(res), figsize=(3.2 * len(res), 3.0), squeeze=False)
    allpos = np.concatenate([d["dens"][d["dens"] > 0] for d in res.values()])
    norm = LogNorm(vmin=np.percentile(allpos, 5), vmax=np.percentile(allpos, 99.5)) if allpos.size else None
    for ax, (r, d) in zip(axes[0], res.items()):
        x, z = d["xyz"][:, 0], d["xyz"][:, 2]; v = np.where(d["dens"] > 0, d["dens"], np.nan)
        sc = ax.scatter(x, z, s=0.6, c=v, cmap="viridis", norm=norm, linewidths=0, rasterized=True)
        ax.set_title(f"{r.split('__')[-1]}\nCV {d['cv']:.0f}%, top 10% {d['top']:.0f}%", fontsize=8)
        ax.set_xlabel("x (m)", fontsize=8); ax.set_ylabel("z (m)", fontsize=8); ax.tick_params(labelsize=7); ax.set_aspect("equal")
    cb_ = fig.colorbar(sc, ax=axes[0].tolist(), shrink=0.8); cb_.set_label(f"Ca + Mg released, {a.t0:g}-{a.t1:g} y (mol m$^{{-3}}$)", fontsize=7); cb_.ax.tick_params(labelsize=6)
    out = os.path.join(REV, "figures", "dissolution_map.png"); os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight"); print("written:", os.path.relpath(out, REV))


if __name__ == "__main__":
    main()
