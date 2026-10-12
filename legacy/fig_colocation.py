#!/usr/bin/env python3
"""
Flow, dissolution and precipitation on the same networks (AE-2, R2-6).

The associate editor asked for a figure showing the flow paths, the dissolution
zones and the precipitation patterns together, so that the flow-reactive
co-location argument can be assessed visually rather than only through a
correlation coefficient.

WHAT IS PLOTTED
---------------
Two realisations at the SAME fracture intensity, one that traps and one that
does not, in two rows. Four columns:

  (a) flow        cell flux magnitude, computed from the simulated pressure
                  field and the mesh connectivity by Darcy's law. Log colour
                  scale, because fluxes span several orders of magnitude.
  (b) dissolution anorthite volume fraction consumed between t=0 and t=50 yr,
                  i.e. where cations are being released.
  (c) carbonate   net carbonate volume fraction at 50 yr, plotted as markers
                  rather than a field: precipitation occupies 7 to 740 cells of
                  18,000 to 253,000, so a continuous colour map renders it
                  invisible.
  (d) overlay     the principal flow paths in grey with the precipitating cells
                  in red on top. This is the panel that carries the argument.

HOW THE FLUX IS OBTAINED
------------------------
Rather than re-solving the flow problem, the figure uses the pressure field
PFLOTRAN actually simulated. For each connection in the .uge file,

    q_ij = (k/mu) * A_ij * (P_i - P_j) / d_ij

with A_ij the interface area from the .uge (aperture-corrected) and d_ij the
centroid separation. The per-cell flux magnitude is half the sum of |q| over a
cell's connections, which for a cell with one inlet and one outlet recovers the
throughflow. Only relative magnitudes matter for the figure, so k/mu is set to
1 and the colour bar is labelled in arbitrary units.

This has the advantage of showing the field the reactive transport actually
experienced, including any effect of the corrected volumes and areas.

Usage
-----
    python3 src/fig_colocation.py --runs runs \
        --high C_baseline__p32_150_s117 --low C_baseline__p32_150_s259 \
        --out figures/fig_colocation.pdf
    python3 src/fig_colocation.py --runs runs --list      # what is available
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEED_VF = 1e-6


# ---------------------------------------------------------------- mesh reading
def read_uge(path):
    """Cell centroids and volumes, plus the connection list.

    .uge cell line       : id  x  y  z  volume
    .uge connection line : id1 id2 x  y  z  area
    """
    with open(path) as f:
        n = int(f.readline().split()[1])
        xyz = np.empty((n, 3)); vol = np.empty(n)
        for i in range(n):
            t = f.readline().split()
            xyz[i] = (float(t[1]), float(t[2]), float(t[3]))
            vol[i] = float(t[4])
        line = f.readline().split()
        m = int(line[1])
        ids = np.empty((m, 2), dtype=np.int64); area = np.empty(m)
        for j in range(m):
            t = f.readline().split()
            ids[j] = (int(t[0]) - 1, int(t[1]) - 1)      # to 0-based
            area[j] = float(t[5])
    return xyz, vol, ids, area


def cell_flux(xyz, ids, area, pressure):
    """Per-cell flux magnitude from the simulated pressure field.

    q_ij = A_ij * (P_i - P_j) / d_ij, with k/mu = 1 (relative units only).
    The per-cell value is half the sum of |q| over its connections, which for
    a cell with one inlet and one outlet equals the throughflow.
    """
    i, j = ids[:, 0], ids[:, 1]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    d[d <= 0] = np.nan
    q = area * (pressure[i] - pressure[j]) / d
    q = np.nan_to_num(q, nan=0.0, posinf=0.0, neginf=0.0)
    acc = np.zeros(len(xyz))
    np.add.at(acc, i, np.abs(q))
    np.add.at(acc, j, np.abs(q))
    return 0.5 * acc


# ------------------------------------------------------------- field extraction
def load_case(run_dir, final_year=50.0):
    import h5py
    from validate_run import time_groups

    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5")))
          if "dfn_properties" not in f]
    uge = os.path.join(run_dir, "full_mesh.uge")
    if not h5 or not os.path.isfile(uge):
        return None, "no output"
    xyz, vol, ids, area = read_uge(uge)

    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg:
            return None, "no time groups"
        if tg[-1][0] < final_year * 0.999:
            return None, f"incomplete ({tg[-1][0]:g} yr)"
        g0, g1 = f[tg[0][1]], f[tg[-1][1]]

        def fld(grp, *prefixes):
            for p in prefixes:
                k = next((x for x in grp.keys() if x.startswith(p)), None)
                if k is not None:
                    return np.asarray(grp[k][:], dtype=float)
            return None

        press = fld(g1, "Liquid Pressure", "Liquid_Pressure", "Pressure")
        anor0 = fld(g0, "Anorthite VF")
        anor1 = fld(g1, "Anorthite VF")
        carb = np.zeros(len(xyz))
        for m in CARB:
            a = fld(g1, f"{m} VF")
            if a is not None:
                carb += np.clip(a - SEED_VF, 0.0, None)
        available = sorted(g1.keys())

    if press is None:
        return None, ("no pressure field in the output; add "
                      "LIQUID_PRESSURE to the OUTPUT block, or the flux "
                      "panel cannot be drawn")
    if anor0 is None or anor1 is None:
        return None, "no anorthite field"
    for nm, arr in (("pressure", press), ("anorthite", anor1), ("carbonate", carb)):
        if arr.size != len(xyz):
            return None, f"{nm} length {arr.size} != {len(xyz)} cells"

    return dict(xyz=xyz, vol=vol,
                flux=cell_flux(xyz, ids, area, press),
                diss=anor0 - anor1,
                carb=carb,
                n_carb=int((carb > 0).sum()),
                carb_per_cell=float(carb.sum() / len(xyz)),
                available=available), None


# --------------------------------------------------------------------- plotting
def panel(ax, xyz, c, *, log=False, cmap="viridis", size=0.6,
          vmin=None, vmax=None, label="", markers=False, mask=None):
    """Project onto the flow plane (x horizontal, z vertical)."""
    x, z = xyz[:, 0], xyz[:, 2]
    if markers:
        sel = mask if mask is not None else (c > 0)
        ax.scatter(x, z, s=0.15, c="0.86", linewidths=0, rasterized=True)
        if sel.any():
            sc = ax.scatter(x[sel], z[sel], s=14, c=c[sel],
                            cmap=cmap, norm=_norm(c[sel], log, vmin, vmax),
                            linewidths=0.3, edgecolors="k", rasterized=True)
        else:
            sc = None
    else:
        sc = ax.scatter(x, z, s=size, c=c, cmap=cmap,
                        norm=_norm(c, log, vmin, vmax),
                        linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.4)
    if label:
        ax.set_title(label, fontsize=8, pad=3)
    return sc


def _norm(c, log, vmin, vmax):
    from matplotlib.colors import LogNorm, Normalize
    finite = c[np.isfinite(c)]
    if log:
        pos = finite[finite > 0]
        if pos.size == 0:
            return Normalize()
        lo = vmin if vmin is not None else np.percentile(pos, 5)
        hi = vmax if vmax is not None else pos.max()
        return LogNorm(vmin=max(lo, hi * 1e-6), vmax=hi)
    return Normalize(vmin=vmin if vmin is not None else finite.min(),
                     vmax=vmax if vmax is not None else finite.max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--high", default="C_baseline__p32_150_s117")
    ap.add_argument("--low", default="C_baseline__p32_150_s259")
    ap.add_argument("--out", default="figures/fig_colocation.pdf")
    ap.add_argument("--list", action="store_true",
                    help="report trapping per case so a contrasting pair can be chosen")
    a = ap.parse_args()

    if a.list:
        rows = []
        for d in sorted(glob.glob(os.path.join(a.runs, "C_baseline__*"))):
            r, why = load_case(d)
            nm = os.path.basename(d).replace("C_baseline__", "")
            if r is None:
                print(f"  {nm:<18} {why}")
            else:
                rows.append((nm, r["carb_per_cell"], r["n_carb"]))
        print(f"\n{'case':<18}{'carb/cell':>12}{'cells':>8}")
        for nm, v, n in sorted(rows, key=lambda t: -t[1]):
            print(f"{nm:<18}{v:>12.3e}{n:>8}")
        print("\nPick two at the SAME P32: one near the top, one near the bottom.")
        return

    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 8, "axes.linewidth": 0.4})
    import matplotlib.pyplot as plt

    cases = []
    for tag, name in (("traps", a.high), ("does not trap", a.low)):
        d = os.path.join(a.runs, name)
        r, why = load_case(d)
        if r is None:
            sys.exit(f"{name}: {why}")
        r["name"] = name.replace("C_baseline__", "")
        r["tag"] = tag
        cases.append(r)
        print(f"  {r['name']:<18} {r['n_carb']:>5} carbonate cells, "
              f"{r['carb_per_cell']:.3e} per cell")

    # common colour limits so the two rows are comparable
    fl = np.concatenate([c["flux"][c["flux"] > 0] for c in cases])
    fl_lo, fl_hi = np.percentile(fl, 25), np.percentile(fl, 99.9)
    ds = np.concatenate([c["diss"] for c in cases])
    ds_lo, ds_hi = np.percentile(ds, 2), np.percentile(ds, 98)
    cb = np.concatenate([c["carb"][c["carb"] > 0] for c in cases
                         if (c["carb"] > 0).any()])
    cb_lo, cb_hi = (cb.min(), cb.max()) if cb.size else (1e-9, 1e-6)

    fig, axes = plt.subplots(2, 4, figsize=(11.0, 5.4))
    titles = ("(a) flow: cell flux magnitude",
              "(b) dissolution: anorthite consumed",
              "(c) carbonate at 50 yr",
              "(d) flow paths with precipitation")
    for row, c in enumerate(cases):
        s0 = panel(axes[row, 0], c["xyz"], c["flux"], log=True, cmap="Blues",
                   vmin=fl_lo, vmax=fl_hi,
                   label=titles[0] if row == 0 else "")
        s1 = panel(axes[row, 1], c["xyz"], c["diss"], cmap="Oranges",
                   vmin=ds_lo, vmax=ds_hi,
                   label=titles[1] if row == 0 else "")
        s2 = panel(axes[row, 2], c["xyz"], c["carb"], log=True, cmap="Reds",
                   vmin=cb_lo, vmax=cb_hi, markers=True,
                   label=titles[2] if row == 0 else "")
        # (d) the principal flow paths, with precipitation on top
        ax = axes[row, 3]
        x, z = c["xyz"][:, 0], c["xyz"][:, 2]
        ax.scatter(x, z, s=0.15, c="0.9", linewidths=0, rasterized=True)
        thr = np.percentile(c["flux"][c["flux"] > 0], 90) if (c["flux"] > 0).any() else 0
        fast = c["flux"] >= thr
        ax.scatter(x[fast], z[fast], s=0.5, c="0.35", linewidths=0, rasterized=True)
        sel = c["carb"] > 0
        if sel.any():
            ax.scatter(x[sel], z[sel], s=16, c="#c1272d", linewidths=0.3,
                       edgecolors="k", rasterized=True, zorder=3)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.4)
        if row == 0:
            ax.set_title(titles[3], fontsize=8, pad=3)
        axes[row, 0].set_ylabel(f"{c['name']}\n({c['tag']})", fontsize=8)

    for col, (sc, lab) in enumerate(((s0, "flux (arb. units)"),
                                     (s1, "$\\Delta$VF anorthite"),
                                     (s2, "carbonate VF"))):
        if sc is None:
            continue
        cax = fig.add_axes([0.055 + col * 0.2455, 0.07, 0.14, 0.016])
        cbar = fig.colorbar(sc, cax=cax, orientation="horizontal")
        cbar.set_label(lab, fontsize=7); cbar.ax.tick_params(labelsize=6)
    fig.text(0.80, 0.075, "grey: fastest 10% of cells\nred: precipitating cells",
             fontsize=7, ha="center")

    fig.subplots_adjust(left=0.045, right=0.99, top=0.94, bottom=0.13,
                        wspace=0.06, hspace=0.06)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300)
    fig.savefig(a.out.rsplit(".", 1)[0] + ".png", dpi=200)
    print(f"\nwrote {a.out}")
    print(f"      {a.out.rsplit('.',1)[0]}.png")


if __name__ == "__main__":
    main()
