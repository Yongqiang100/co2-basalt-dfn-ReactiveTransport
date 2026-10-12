#!/usr/bin/env python3
"""
Two presentations of flow-reactive co-location (AE-2, R2-6).

The associate editor asked for the flow paths, dissolution zones and
precipitation patterns to be shown together, and for an account of how
negligible-carbonate cases were handled. Two approaches are provided because
they answer slightly different questions and it is not obvious in advance which
reads better on real data.

    --mode profile   Every quantity binned by distance from the injection face
                     and plotted as a median with an interquartile band. This
                     asks what the fluid encounters as it moves downstream, in
                     aggregate over all cells rather than along one streamline.
    --mode overlay   One panel per realisation with the network, the fastest
                     flow paths, the dissolution field and the precipitating
                     cells superimposed, so spatial coincidence is visible
                     directly.

WHY NOT A TRACED STREAMLINE
---------------------------
Following the largest flux from cell to cell was tried first and does not work
on these meshes: the walk becomes trapped in local recirculation, and the paths
it produces span only a few metres of a 20 m domain, never reaching either the
injection region or the outflow. On one realisation the traced path lay entirely
in unreacted formation water (pH 7.54 throughout) and on another pH decreased
downstream, which is the wrong sense for a path leaving an acid injector. Both
modes below avoid that failure: the profile mode uses the injection-to-outflow
coordinate directly, and the overlay mode makes no path assumption at all.

WHAT THE PROFILE MODE SHOWS
---------------------------
Five rows, top to bottom, as the causal chain the manuscript claims:

  1  anorthite dissolution rate  -- where cations enter the fluid
  2  Ca and Mg carbonate complexes -- whether they accumulate downstream
  3  pH -- whether the proton load has been consumed
  4  carbonate precipitation rate -- where precipitation is active at 50 yr
  5  accumulated carbonate volume fraction -- what has formed

If the trapping realisation shows cations rising and pH recovering upstream of
its carbonate, while the non-trapping realisation shows neither, the mechanism
is demonstrated rather than asserted. If not, that is equally worth knowing.

Usage
-----
    python3 src/fig_coloc2.py --runs runs --mode profile \\
        --high C_baseline__p32_100_s383 --low C_baseline__p32_100_s42 \\
        --out figures/fig_profile.pdf

    python3 src/fig_coloc2.py --runs runs --mode overlay \\
        --high C_baseline__p32_100_s383 --low C_baseline__p32_100_s42 \\
        --out figures/fig_overlay.pdf
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEED_VF = 1e-6
PHYS_VF_MAX = 1.0          # a cell cannot hold more carbonate than its volume


def read_uge(path):
    with open(path) as f:
        n = int(f.readline().split()[1])
        xyz = np.empty((n, 3)); vol = np.empty(n)
        for i in range(n):
            t = f.readline().split()
            xyz[i] = (float(t[1]), float(t[2]), float(t[3])); vol[i] = float(t[4])
        m = int(f.readline().split()[1])
        ids = np.empty((m, 2), dtype=np.int64); area = np.empty(m)
        for j in range(m):
            t = f.readline().split()
            ids[j] = (int(t[0]) - 1, int(t[1]) - 1); area[j] = float(t[5])
    return xyz, vol, ids, area


def gather(run_dir, final_year=50.0):
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

        def fld(grp, *pre):
            for p in pre:
                k = next((x for x in grp.keys() if x.startswith(p)), None)
                if k is not None:
                    return np.asarray(grp[k][:], dtype=float)
            return None

        def tot(*pre):
            out = None
            for p in pre:
                v = fld(g1, p)
                if v is not None:
                    out = v if out is None else out + v
            return out

        press = fld(g1, "Liquid Pressure", "Pressure")
        ph = fld(g1, "pH")
        anor0, anor1 = fld(g0, "Anorthite VF"), fld(g1, "Anorthite VF")
        diss_rate = fld(g1, "Anorthite Rate")
        prec_rate = tot(*[f"{m} Rate" for m in CARB])
        cations = tot("CaCO3(aq)", "CaHCO3+", "MgCO3(aq)", "MgHCO3+")
        carb = np.zeros(len(xyz))
        for m in CARB:
            a = fld(g1, f"{m} VF")
            if a is not None:
                carb += np.clip(a - SEED_VF, 0.0, None)

    if press is None:
        return None, "no pressure field in the output"

    # cell flux magnitude from the simulated pressure field
    i, j = ids[:, 0], ids[:, 1]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1); d[d <= 0] = np.nan
    q = np.nan_to_num(area * (press[i] - press[j]) / d, nan=0.0,
                      posinf=0.0, neginf=0.0)
    flux = np.zeros(len(xyz))
    np.add.at(flux, i, np.abs(q)); np.add.at(flux, j, np.abs(q)); flux *= 0.5

    # cells whose carbonate exceeds their own volume are unphysical -- the
    # fixed-porosity model does not prevent it -- and are excluded from the
    # plotted fields rather than allowed to set the colour or axis limits
    unphys = carb > PHYS_VF_MAX
    return dict(name=os.path.basename(run_dir).replace("C_baseline__", ""),
                xyz=xyz, vol=vol, flux=flux, ph=ph,
                diss=(anor0 - anor1) if anor0 is not None else None,
                diss_rate=diss_rate, prec_rate=prec_rate,
                cations=cations, carb=carb, unphys=unphys,
                n_carb=int((carb > 0).sum()), n_unphys=int(unphys.sum()),
                carb_per_cell=float(carb.sum() / len(xyz))), None


# ------------------------------------------------------------------ mode: profile
def binned(x, y, edges, keep=None):
    """Median and interquartile range of y in each x bin."""
    if y is None:
        return None
    ok = np.isfinite(y)
    if keep is not None:
        ok &= keep
    idx = np.digitize(x, edges) - 1
    med = np.full(len(edges) - 1, np.nan)
    lo = np.full(len(edges) - 1, np.nan)
    hi = np.full(len(edges) - 1, np.nan)
    for b in range(len(edges) - 1):
        sel = ok & (idx == b)
        if sel.sum() >= 3:
            v = y[sel]
            med[b], lo[b], hi[b] = np.median(v), *np.percentile(v, [25, 75])
    return med, lo, hi


def mode_profile(cases, out, nbins=40):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 8, "axes.linewidth": 0.5,
                                "xtick.labelsize": 7, "ytick.labelsize": 7})
    import matplotlib.pyplot as plt

    COL = {"traps": "#1b5e9c", "does not trap": "#b8532a"}
    rows = [("anorthite dissolution rate", "|rate|\n(mol m$^{-3}$ s$^{-1}$)",
             "diss_rate", True, True),
            ("Ca and Mg carbonate complexes", "concentration\n(M)",
             "cations", True, False),
            ("pH", "pH", "ph", False, False),
            ("carbonate precipitation rate", "rate\n(mol m$^{-3}$ s$^{-1}$)",
             "prec_rate", True, True),
            ("accumulated carbonate", "carbonate VF", "carb", True, False)]
    rows = [r for r in rows if any(c.get(r[2]) is not None for c in cases)]

    x0 = min(c["xyz"][:, 0].min() for c in cases)
    x1 = max(c["xyz"][:, 0].max() for c in cases)
    edges = np.linspace(x0, x1, nbins + 1)
    ctr = 0.5 * (edges[:-1] + edges[1:]) - x0        # distance from inlet face

    fig, axes = plt.subplots(len(rows), 1, figsize=(6.4, 1.35 * len(rows)),
                             sharex=True)
    if len(rows) == 1:
        axes = [axes]
    for ax, (title, ylab, key, logy, absval) in zip(axes, rows):
        for c in cases:
            y = c.get(key)
            if y is None:
                continue
            if absval:
                y = np.abs(y)
            keep = ~c["unphys"] if key in ("carb", "prec_rate") else None
            b = binned(c["xyz"][:, 0], y, edges, keep)
            if b is None:
                continue
            med, lo, hi = b
            col = COL[c["tag"]]
            ax.fill_between(ctr, lo, hi, color=col, alpha=0.18, lw=0)
            ax.plot(ctr, med, lw=1.1, color=col,
                    label=f"{c['name']} ({c['tag']})")
        if logy:
            # Test the plotted collections as well as the lines: in a bin where
            # most cells do not precipitate the median is zero while the upper
            # quartile is not, so checking lines alone mislabels the panel.
            vals = [l.get_ydata() for l in ax.get_lines()]
            for coll in ax.collections:
                try:
                    vals.append(np.asarray(coll.get_paths()[0].vertices)[:, 1])
                except Exception:
                    pass
            pos = [np.nanmax(v) for v in vals if len(v) and np.isfinite(v).any()]
            if pos and max(pos) > 0:
                ax.set_yscale("log")
                # a median of zero cannot appear on a log axis; say so rather
                # than leaving the reader to infer it from a missing line
                for c_ in cases:
                    y_ = c_.get(key)
                    if y_ is None:
                        continue
                    b_ = binned(c_["xyz"][:, 0], np.abs(y_) if absval else y_,
                                edges, ~c_["unphys"] if key in ("carb", "prec_rate") else None)
                    if b_ is not None and np.nansum(b_[0]) == 0:
                        ax.text(0.5, 0.06,
                                f"{c_['name']}: bin medians are zero "
                                f"(precipitation is confined to few cells)",
                                transform=ax.transAxes, ha="center",
                                fontsize=6, color="0.45")
                        break
            else:
                ax.text(0.5, 0.5, "zero throughout", transform=ax.transAxes,
                        ha="center", va="center", fontsize=7, color="0.5")
        ax.set_ylabel(ylab, fontsize=7)
        ax.text(0.006, 0.90, title, transform=ax.transAxes, fontsize=8,
                va="top", fontweight="bold")
        ax.grid(alpha=0.25, lw=0.4)
        ax.tick_params(length=2.5)
        ax.axvspan(0, 0.2 * (x1 - x0), color="#4a7fb5", alpha=0.06, lw=0)

    axes[0].legend(fontsize=6.5, loc="lower right", framealpha=0.9)
    axes[-1].set_xlabel("distance from the injection face (m)", fontsize=8)
    fig.text(0.5, 0.012,
             "Line: bin median. Band: interquartile range over all cells in "
             "the bin. Shaded strip: injection region.",
             fontsize=6.5, ha="center", color="0.3")
    fig.subplots_adjust(left=0.135, right=0.985, top=0.99, bottom=0.105,
                        hspace=0.12)
    save(fig, out)


# ------------------------------------------------------------------ mode: overlay
def mode_overlay(cases, out):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 8, "axes.linewidth": 0.5})
    import matplotlib.pyplot as plt

    ds = np.concatenate([c["diss"] for c in cases if c["diss"] is not None])
    ds_lo, ds_hi = np.percentile(ds, [2, 98])

    fig, axes = plt.subplots(1, len(cases), figsize=(4.6 * len(cases), 4.9))
    if len(cases) == 1:
        axes = [axes]
    for ax, c in zip(axes, cases):
        x, z = c["xyz"][:, 0], c["xyz"][:, 2]
        # dissolution field as the backdrop
        sc = ax.scatter(x, z, s=0.7, c=c["diss"], cmap="Oranges",
                        vmin=ds_lo, vmax=ds_hi, linewidths=0, rasterized=True)
        # the fastest tenth of cells, i.e. the principal flow paths
        thr = np.percentile(c["flux"][c["flux"] > 0], 90)
        fast = c["flux"] >= thr
        ax.scatter(x[fast], z[fast], s=1.6, c="#12406e", linewidths=0,
                   alpha=0.55, rasterized=True, label="fastest 10% of cells")
        # precipitation, excluding cells holding more than their own volume
        sel = (c["carb"] > 0) & ~c["unphys"]
        if sel.any():
            ax.scatter(x[sel], z[sel], s=26, facecolors="none",
                       edgecolors="#c1272d", linewidths=0.9, zorder=4,
                       label=f"precipitating cells (n={int(sel.sum())})")
        if c["unphys"].any():
            ax.scatter(x[c["unphys"]], z[c["unphys"]], s=26, marker="x",
                       c="0.35", linewidths=0.8, zorder=4,
                       label=f"excluded, VF>1 (n={c['n_unphys']})")
        ax.axvspan(x.min(), x.min() + 0.2 * (x.max() - x.min()),
                   color="#4a7fb5", alpha=0.07, lw=0)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{c['name']} — {c['tag']}\n"
                     f"{c['n_carb']} precipitating cells, "
                     f"{c['carb_per_cell']:.2e} per cell", fontsize=8)
        ax.legend(fontsize=6.5, loc="upper left", framealpha=0.92,
                  markerscale=1.2, borderpad=0.35, handletextpad=0.5)
        # scale bar
        sb = 5.0
        bx, bz = x.max() - 0.08 * (x.max() - x.min()) - sb, z.min() + 0.05 * (z.max() - z.min())
        ax.plot([bx, bx + sb], [bz, bz], lw=1.5, c="k", solid_capstyle="butt")
        ax.text(bx + sb / 2, bz + 0.025 * (z.max() - z.min()), f"{sb:g} m",
                fontsize=6.5, ha="center")

    cax = fig.add_axes([0.32, 0.085, 0.36, 0.016])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("anorthite volume fraction consumed over 50 yr", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    fig.text(0.5, 0.018,
             "Projection: $x$ (injection to outflow) horizontal, $z$ vertical. "
             "Shaded strip: injection region.",
             fontsize=6.5, ha="center", color="0.3")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.155, wspace=0.04)
    save(fig, out)


def save(fig, out):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, dpi=300)
    fig.savefig(out.rsplit(".", 1)[0] + ".png", dpi=200)
    print(f"\nwrote {out}")
    print(f"      {out.rsplit('.',1)[0]}.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--mode", choices=("profile", "overlay"), default="profile")
    ap.add_argument("--high", default="C_baseline__p32_100_s383")
    ap.add_argument("--low", default="C_baseline__p32_100_s42")
    ap.add_argument("--out", default=None)
    ap.add_argument("--bins", type=int, default=40)
    a = ap.parse_args()
    out = a.out or f"figures/fig_{a.mode}.pdf"

    cases = []
    for tag, name in (("traps", a.high), ("does not trap", a.low)):
        c, why = gather(os.path.join(a.runs, name))
        if c is None:
            sys.exit(f"{name}: {why}")
        c["tag"] = tag
        cases.append(c)
        print(f"  {c['name']:<18} {c['n_carb']:>5} precipitating cells, "
              f"{c['carb_per_cell']:.3e} per cell"
              + (f", {c['n_unphys']} excluded with VF>1" if c["n_unphys"] else ""))

    (mode_profile(cases, out, a.bins) if a.mode == "profile"
     else mode_overlay(cases, out))


if __name__ == "__main__":
    main()
