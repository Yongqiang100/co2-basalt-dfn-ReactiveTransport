#!/usr/bin/env python3
"""
make_all_manuscript_figures.py
==============================
ONE self-contained driver that regenerates ALL TEN figures of the Cooper Basin
Li-ISR manuscript in a single run, on the seed-42 base case + the ten-member
network ensemble.

It depends only on things that are reliably present in the analysis tree:
  * the run outputs              (pflotran_field_results/...)
  * the seed-42 mesh             (full_mesh.inp, for fig:domain only)
  * the core modules            postprocess_field.py, channelization.py
  * generate_revised_figures.py  (ensemble engine -> 4 data figures)
  * ph_sweep_summary.py, co2_sweep_summary.py   (chemistry figures)

The three figures that previously failed when a helper file / input was missing
-- the mechanism schematic, the kinetics-area grid, and the domain projection --
are now built INSIDE this file. The schematic and the verified kinetics grid
need no run data at all; the domain projection reads only the mesh and uses a
pure-numpy convex hull, so it needs neither scipy nor a separate script. Nothing
this driver must have can therefore be "missing on Setonix".

    figure         file in the .tex              built by
    ------         ----------------              --------
    fig:scaling    scaling_loglog.pdf           \\  generate_revised_figures.py
    fig:channel    channelization_composite.pdf  >  (called with empty --mesh /
    fig:bimodal    bimodal_composite.pdf        /    --sweep-csv so it produces
    fig:headline   master_summary.pdf          /     ONLY these four figures)
    fig:ph         ph_sweep_summary.pdf          ph_sweep_summary.py
    fig:co2        co2_sweep_summary.pdf         co2_sweep_summary.py
    fig:ksweep     ksweep_grid.pdf              <inline; verified grid, --sweep-csv override>
    fig:domain     domain_dfn_proj.pdf          <inline; reads mesh, scipy-free>
    fig:calcite    calcite_regime.pdf           <inline; seed-42 mass balance>
    fig:mechanism  mechanism_schematic.pdf      <inline; schematic, no data>

Environment (Setonix):
  module load python/3.11.6 py-numpy/2.1.2 py-h5py/3.12.1 py-matplotlib

Usage:
  python make_all_manuscript_figures.py \\
      --results   pflotran_field_results \\
      --mesh      dfn_library_field/habanero_field_s42/full_mesh.inp \\
      --sweep-csv sweep_kinetics_area_grid.csv \\
      --base 42 --outdir figs
(--sweep-csv is optional; without it the embedded verified grid is used.
 --mesh is only needed for fig:domain; the other nine are produced regardless.)
"""

import os
import time
import sys
import csv as _csv
import glob
import json
import argparse
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["text.usetex"] = False
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
RATES = [1, 5, 10, 25, 50]
INJ_C, PROD_C = "#1A5FB4", "#C01C28"

sys.path.insert(0, HERE)
try:
    import nature_style
    nature_style.apply()
    _HAVE_NATURE = True
except Exception:
    _HAVE_NATURE = False


# ---------------------------------------------------------------- helpers
_RUN_T0 = time.time()


def _script(name):
    """Prefer a copy sitting next to this driver; else rely on PATH/cwd."""
    cand = os.path.join(HERE, name)
    return cand if os.path.exists(cand) else name


def _run(cmd, label):
    print(f"\n>>> [{label}]\n    {' '.join(map(str, cmd))}")
    try:
        subprocess.run(cmd, check=True)
        print(f"    [{label}] OK")
        return True
    except Exception as exc:                                       # noqa: BLE001
        print(f"    [{label}] FAILED: {exc}")
        return False


def _save(fig, outdir, name, dpi=200):
    if _HAVE_NATURE:
        nature_style.finalize(fig)
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig:calcite
def fig_calcite(results, base, outdir):
    """Whole-domain calcite change vs rate, 20 y + 50 y, seed-42 base case."""
    sys.path.insert(0, HERE)
    sys.path.insert(0, os.getcwd())
    try:
        from postprocess_field import analyze_run
    except Exception as exc:                                       # noqa: BLE001
        print(f"    [calcite] cannot import postprocess_field: {exc}")
        return False

    def cal(run_name):
        a = analyze_run(os.path.join(results, run_name))
        if not a:
            print(f"    [calcite] missing run: {run_name}")
            return np.nan
        return a["minerals"]["Calcite"]["delta_pct"]

    c20 = np.array([cal(f"habanero_field_s{base}_rate{r}") for r in RATES])
    c50 = np.array([cal(f"habanero_field_s{base}_rate{r}_50y") for r in RATES])
    finite = np.concatenate([c20[np.isfinite(c20)], c50[np.isfinite(c50)]])
    if finite.size == 0:
        print("    [calcite] no calcite data found -- skipped")
        return False
    ymax = max(0.6, float(np.nanmax(finite)) * 1.25)
    ymin = min(-0.6, float(np.nanmin(finite)) * 1.15)

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.axhspan(0.0, ymax, color="#2e7d32", alpha=0.07)
    ax.axhspan(ymin, 0.0, color="#c62828", alpha=0.07)
    ax.axhline(0.0, color="0.35", lw=1.0)
    ax.text(1.05, ymax * 0.10, "net precipitation", color="#2e7d32", fontsize=10, va="bottom")
    ax.text(1.05, ymin * 0.10, "net dissolution", color="#c62828", fontsize=10, va="top")
    ax.plot(RATES, c20, "o-", color="#1f4e79", lw=1.9, ms=6, label="20 y")
    if np.isfinite(c50).any():
        ax.plot(RATES, c50, "s--", color="#b8731b", lw=1.9, ms=6, label="50 y")
    else:
        print("    [calcite] WARNING: no _50y runs found -- plotting 20 y only")
    # mark the 20-y precip->dissol crossover (Da* ~ 1)
    for i in range(len(RATES) - 1):
        a0, a1 = c20[i], c20[i + 1]
        if np.isfinite(a0) and np.isfinite(a1) and a0 > 0.0 >= a1:
            xc = RATES[i] + (RATES[i + 1] - RATES[i]) * a0 / (a0 - a1)
            ax.axvline(xc, color="0.5", ls=":", lw=1.2)
            # label tied directly to the crossover line (no arrow), so the marker
            # cannot appear to point at both the 20 y and 50 y curves
            ax.text(xc / 1.07, ymin * 0.60, r"$\mathrm{Da}^{\ast}\!\sim\!1$",
                    rotation=90, ha="center", va="center", fontsize=10, color="0.3")
            break
    ax.set_xscale("log")
    ax.set_xticks(RATES)
    ax.set_xticklabels([str(r) for r in RATES])
    ax.set_xlim(0.8, 62)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"Injection rate $Q$  (kg s$^{-1}$)")
    ax.set_ylabel("Whole-domain calcite change  (%)")
    ax.set_title(f"Calcite buffering regime (base case, seed {base})", loc="left", fontsize=10)
    ax.legend(frameon=False, loc="upper right")
    ax.grid(alpha=0.25, which="both")
    _save(fig, outdir, "calcite_regime")
    print(f"    [calcite] wrote calcite_regime.pdf/.png   "
          f"20y={np.array2string(c20, precision=2)}  50y={np.array2string(c50, precision=2)}")
    return True


# ---------------------------------------------------------------- fig:ksweep
# verified k0 x reactive-area grid (= tab:ksweep); base cell (1e-12, 100) = 8.49
_KSWEEP_K0 = [1e-12, 1e-11, 1e-10]
_KSWEEP_A = [50.0, 100.0, 200.0, 500.0]
_KSWEEP_GRID = np.array([
    [4.83,  8.49, 14.96, 27.86],    # k0 = 1e-12
    [27.86, 39.22, 45.03, 46.35],   # k0 = 1e-11
    [46.35, 46.78, 47.10, 47.44],   # k0 = 1e-10
])


def _load_ksweep(csv=None):
    """Return (k0s, As, grid, source) for the rate-constant by area sweep."""
    k0s, As, grid, src = _KSWEEP_K0, _KSWEEP_A, _KSWEEP_GRID, "embedded verified grid"
    if csv and os.path.exists(csv):
        cells = {}
        try:
            with open(csv) as f:
                for row in _csv.DictReader(f):
                    cells[(float(row["k0_mol_m2_s"]), float(row["area_cm2_cm3"]))] = \
                        float(row["avg_li_mg_L"])
        except Exception as exc:                                   # noqa: BLE001
            print(f"    [ksweep] CSV unreadable ({exc}); using embedded grid")
            cells = {}
        if cells:
            k0s = sorted({k for k, _ in cells})
            As = sorted({a for _, a in cells})
            grid = np.array([[cells.get((k, a), np.nan) for a in As] for k in k0s])
            src = csv
    return k0s, As, grid, src


def _draw_ksweep(ax, k0s, As, grid, fs=10, cb=True, fig=None):
    """Draw the heatmap on a supplied axis. Returns the image for a colorbar."""
    im = ax.imshow(grid, cmap="viridis", aspect="auto", origin="upper",
                   vmin=0, vmax=max(50.0, float(np.nanmax(grid))))
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            v = grid[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                        color="white" if v < 30 else "black", fontsize=fs, fontweight="bold")
    ax.set_xticks(range(len(As)), [f"{int(a)}" for a in As])
    ax.set_yticks(range(len(k0s)), [f"$10^{{{int(round(np.log10(k)))}}}$" for k in k0s])
    ax.set_xlabel(r"Reactive surface area $A$  (cm$^2$ cm$^{-3}$)")
    ax.set_ylabel(r"Rate constant $k_0$  (mol m$^{-2}$ s$^{-1}$)")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    if cb and fig is not None:
        bar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        bar.set_label("Produced Li  (mg L$^{-1}$)")
        bar.outline.set_visible(False)
    return im


def fig_ksweep(outdir, csv=None):
    """fig:ksweep -- produced-Li heatmap over k0 x reactive area, standalone."""
    k0s, As, grid, src = _load_ksweep(csv)
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    _draw_ksweep(ax, k0s, As, grid, fs=10, cb=True, fig=fig)
    _save(fig, outdir, "ksweep_grid")
    print(f"    [ksweep] wrote ksweep_grid.pdf/.png   source={src}  "
          f"range {np.nanmin(grid):.1f}-{np.nanmax(grid):.1f}  base(1e-12,100)={grid[k0s.index(1e-12), As.index(100.0)]:.2f}"
          if (1e-12 in k0s and 100.0 in As) else
          f"    [ksweep] wrote ksweep_grid.pdf/.png   source={src}")
    return True


# ---------------------------------------------------------------- fig:domain
def _hull2d(P2):
    """2-D convex hull (Andrew's monotone chain) -- pure numpy, no scipy."""
    pts = sorted(set(map(tuple, np.asarray(P2).tolist())))
    if len(pts) < 3:
        return np.asarray(pts, float)

    def cr(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lo = []
    for p in pts:
        while len(lo) >= 2 and cr(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    up = []
    for p in reversed(pts):
        while len(up) >= 2 and cr(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    return np.asarray(lo[:-1] + up[:-1], float)


def _hull(P2):
    try:
        h = _hull2d(P2)
        if len(h) >= 3:
            return np.vstack([h, h[:1]])
        raise ValueError("degenerate")
    except Exception:                                              # noqa: BLE001
        d = P2 - P2.mean(0)
        _, _, Vt = np.linalg.svd(d, full_matrices=False)
        t = d @ Vt[0]
        return np.vstack([P2[t.argmin()], P2[t.argmax()]])


def fig_domain(mesh, outdir):
    """Three orthographic projections of the real DFN mesh (full_mesh.inp)."""
    if not mesh or not os.path.exists(mesh):
        print(f"    [domain] skipped (mesh not found: {mesh!r})")
        return False
    import matplotlib.cm as cm
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    INJ = np.array([-125.0, 0.0])
    PROD = np.array([125.0, 0.0])
    Z_TOP, Z_BOT, RBOX = -160.0, -210.0, 40.0

    with open(mesh) as fh:
        nnode, nelem = map(int, fh.readline().split()[:2])
        coords = np.empty((nnode, 3))
        for i in range(nnode):
            p = fh.readline().split()
            coords[i] = (float(p[1]), float(p[2]), float(p[3]))
        tris = np.empty((nelem, 3), dtype=np.int64)
        mat = np.empty(nelem, dtype=np.int64)
        for e in range(nelem):
            p = fh.readline().split()
            mat[e] = int(p[1])
            tris[e] = (int(p[3]) - 1, int(p[4]) - 1, int(p[5]) - 1)
    frac_ids = np.unique(mat)
    nfrac = frac_ids.size
    xmin, ymin, zmin = coords.min(0)
    xmax, ymax, zmax = coords.max(0)
    print(f"    [domain] nfrac={nfrac} bounds x[{xmin:.0f},{xmax:.0f}] "
          f"y[{ymin:.0f},{ymax:.0f}] z[{zmin:.0f},{zmax:.0f}]")
    rng = np.random.default_rng(7)
    COLS = cm.turbo(np.linspace(0.05, 0.95, nfrac))[rng.permutation(nfrac)]
    FN = [coords[np.unique(tris[mat == f].ravel())] for f in frac_ids]

    def draw_frac(ax, cols):
        for k, P in enumerate(FN):
            tr = _hull(P[:, cols])
            ax.fill(tr[:, 0], tr[:, 1], facecolor=COLS[k], edgecolor="none", alpha=0.07, zorder=1)
            ax.plot(tr[:, 0], tr[:, 1], color=COLS[k], lw=0.5, alpha=0.6, zorder=2)

    PAD = 10
    # Right column twice the width of the left so the spanning plan-view square
    # aligns top-to-bottom with the two stacked cross-section squares;
    # constrained_layout then packs the equal-aspect axes without dead space.
    fig = plt.figure(figsize=(11.0, 7.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 2.0], height_ratios=[1, 1])
    axXZ = fig.add_subplot(gs[0, 0])
    axYZ = fig.add_subplot(gs[1, 0])
    axXY = fig.add_subplot(gs[:, 1])

    # (c) X-Y plane (plan)
    draw_frac(axXY, [0, 1])
    axXY.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,
                             fill=False, edgecolor="0.4", lw=1.2))
    for (wx, wy), c, lab in [(INJ, INJ_C, "Injector"), (PROD, PROD_C, "Producer")]:
        axXY.add_patch(Rectangle((wx - RBOX, wy - RBOX), 2 * RBOX, 2 * RBOX,
                                 fill=False, edgecolor=c, lw=1.8, zorder=6))
        axXY.scatter([wx], [wy], color=c, s=55, zorder=7, edgecolor="white", linewidth=0.9)
        axXY.annotate(lab, (wx, wy + RBOX + 12), color=c, ha="center",
                      fontsize=9, fontweight="bold", zorder=7)
    axXY.annotate("", xy=(PROD[0], -150), xytext=(INJ[0], -150),
                  arrowprops=dict(arrowstyle="<->", color="0.15", lw=1.4), zorder=7)
    axXY.text(0, -162, "250 m", ha="center", va="top", fontsize=9, color="0.15", zorder=7)
    axXY.annotate("N", xy=(xmax - 38, ymax - 30), xytext=(xmax - 38, ymax - 88),
                  arrowprops=dict(arrowstyle="->", color="0.15", lw=1.4),
                  ha="center", fontsize=10, fontweight="bold", color="0.15", zorder=7)
    axXY.set_xlim(xmin - PAD, xmax + PAD); axXY.set_ylim(ymin - PAD, ymax + PAD)
    axXY.set_aspect("equal")
    axXY.set_xlabel("$x$  (m)"); axXY.set_ylabel("$y$  (m)")
    axXY.set_title(r"(c)  $x$–$y$ plane, plan view along $z$", fontsize=11)

    # (a) X-Z plane
    draw_frac(axXZ, [0, 2])
    for wx, c in [(INJ[0], INJ_C), (PROD[0], PROD_C)]:
        axXZ.plot([wx, wx], [Z_BOT, Z_TOP], color=c, lw=5, zorder=6, solid_capstyle="round")
    axXZ.set_xlim(xmin - PAD, xmax + PAD); axXZ.set_ylim(zmin - PAD, zmax + PAD)
    axXZ.set_aspect("equal")
    axXZ.set_xlabel("$x$  (m)"); axXZ.set_ylabel("$z$  (m)")
    axXZ.set_title(r"(a)  $x$–$z$ plane, view along $y$", fontsize=11)

    # (b) Y-Z plane
    draw_frac(axYZ, [1, 2])
    axYZ.plot([0, 0], [Z_BOT, Z_TOP], color="0.2", lw=5, zorder=6, solid_capstyle="round")
    axYZ.set_xlim(ymin - PAD, ymax + PAD); axYZ.set_ylim(zmin - PAD, zmax + PAD)
    axYZ.set_aspect("equal")
    axYZ.set_xlabel("$y$  (m)"); axYZ.set_ylabel("$z$  (m)")
    axYZ.set_title(r"(b)  $y$–$z$ plane, view along $x$", fontsize=11)

    hh = [Line2D([0], [0], color=INJ_C, lw=6, label="Injector"),
          Line2D([0], [0], color=PROD_C, lw=6, label="Producer"),
          Line2D([0], [0], color="0.6", lw=1.5, label="Domain outline (panel c)"),
          Line2D([0], [0], marker='s', ls='none', mfc='none', mec='0.4', ms=10,
                 label="40 m drainage box")]
    fig.legend(handles=hh, loc="outside lower center", ncol=4,
               fontsize=9, frameon=False, handletextpad=0.7, columnspacing=1.8)
    if _HAVE_NATURE:
        nature_style.finalize(fig)
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"domain_dfn_proj.{ext}"), dpi=300)
    plt.close(fig)
    print(f"    [domain] wrote domain_dfn_proj.pdf/.png  ({nfrac} fractures)")
    return True


# ---------------------------------------------------------------- fig:mechanism
def fig_study_design(outdir):
    """fig:design -- study-design schematic. Needs no run data.

    Everything is contained in this function and drawn inside an rc_context,
    because two of its settings would corrupt the other figures if applied
    globally:
      * the sans-serif (Arial) family and stixsans mathtext, and
      * savefig.bbox = "standard", i.e. cropping DISABLED.
    The second matters most: the canvas is exactly \\textwidth (390 pt) wide with
    the axes filling it, so \\includegraphics[width=\\textwidth] reproduces it 1:1
    and every point size prints at its nominal value. The tight bbox used by
    _save() elsewhere in this file would crop it and destroy that calibration.
    This is why the figure used to be a separate subprocess.
    """
    RC = {"font.family": "sans-serif",
          "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                              "TeX Gyre Heros", "Nimbus Sans", "DejaVu Sans"],
          "mathtext.fontset": "stixsans",
          "savefig.dpi": 300, "savefig.bbox": "standard"}
    with matplotlib.rc_context(RC):
        FIG_W = 5.40                       # inches; == \textwidth -> printed 1:1

        ORANGE = ("#c9781f", "#fbf0e2", "#c9781f")
        BLUE   = ("#2c6fbb", "#eaf1fb", "#2c6fbb")
        GREEN  = ("#2e8b57", "#e8f5ee", "#2e8b57")
        PURPLE = ("#7b4fa3", "#f3ecf9", "#7b4fa3")
        GREY_F, GREY_E, INK = "#eef1f4", "#8794a2", "#333333"
        COLS = [ORANGE, BLUE, GREEN, PURPLE]

        LBL_X = 0.055
        X_LEFT, X_RIGHT, GAP = 0.120, 0.955, 0.016
        BW = (X_RIGHT - X_LEFT - 3 * GAP) / 4
        CX = [X_LEFT + BW / 2 + i * (BW + GAP) for i in range(4)]

        M_TOP, M_BOT = 4.0, 4.0
        H_BAN, G_BC, H_CARD, G_CR, H_RES, G_RB = 14.0, 6.0, 45.0, 6.0, 51.0, 11.0
        FS_LBL, FS_BAN, FS_CARD, FS_RUN, FS_RES = 6.5, 7.0, 6.5, 6.5, 6.3
        FS_BAND_TITLE, FS_BULLET = 7.5, 6.3
        BAND_PAD, BAND_HANG = 0.014, 0.012

        TITLE = ["Injection rate", "Fracture intensity",
                 "Injectate chemistry", "Reaction kinetics"]
        CARD = [
            "five rates, 1$-$50 kg s$^{-1}$\nten network realizations\n(+ 50-year extension)",
            "$P_{32}$ = 0.05$-$0.15 m$^{-1}$\nfive intensities,\ntwo networks",
            "injection pH 3$-$8\nCO$_2$-saturated,\nfive rates",
            "$k_0A$ over four decades\n12-cell grid plus\nlow-reactivity cases",
        ]
        NRUN = ["55 simulations", "9 simulations", "9 simulations", "16 simulations"]
        FIND = [
            "Lithium amount rises,\nLi-concentration falls:\n"
            "$M_{\\mathrm{Li}}\\propto Q^{0.41}$\n$C\\propto Q^{-0.59}$",
            "Comparable\nfirst-order control\non Li-concentration",
            "Secondary:\n$\\leq$10 % increase\n(reagent confined\nto channels)",
            "Li-concentration\nincreases with the\nrate$-$area product:\n"
            "$C\\propto (k_0A)^{0.4}$\n(1.2$-$47 mg L$^{-1}$)",
        ]
        BAND_TITLE = ("Recovery is transport-limited: flow channelization rather than\n"
                      "injectate chemistry governs the Li-concentration")
        BULLETS = [
            "Flow is confined to a limited zone of high-transmissivity fractures "
            "($\\mathrm{Da}<1$); fluid is produced before it can approach spodumene "
            "saturation",
            "The producer recovers preferentially from the fast flow channels; 96 % of "
            "the dissolved lithium resides in stagnant zones already close to saturation",
            "Improved lithium recovery depends on enlarging the swept volume rather than "
            "on further chemical enhancement",
        ]

        # -- measure the bullet wrapping on a throwaway canvas of the same width.
        #    The wrap depends only on the figure WIDTH, which is fixed, so the band
        #    height can be derived from the resulting line count and the canvas made
        #    exactly as tall as its content.
        def _wrap(text, max_pt, fontsize, fig):
            r = fig.canvas.get_renderer()
            probe = fig.text(0, 0, "", fontsize=fontsize)
            scale = 72.0 / fig.dpi
            lines, cur = [], ""
            for w in text.split():
                trial = (cur + " " + w).strip()
                probe.set_text(trial)
                if cur and probe.get_window_extent(renderer=r).width * scale > max_pt:
                    lines.append(cur); cur = w
                else:
                    cur = trial
            if cur:
                lines.append(cur)
            probe.remove()
            return lines

        band_w_pt = (X_RIGHT - X_LEFT) * FIG_W * 72.0
        avail = band_w_pt - (2 * BAND_PAD + BAND_HANG) * FIG_W * 72.0
        tmp = plt.figure(figsize=(FIG_W, 4.0)); tmp.canvas.draw()
        WRAPPED = [_wrap(b, avail, FS_BULLET, tmp) for b in BULLETS]
        # widest line actually produced: greedy wrapping stops short of `avail` by
        # however much the next word would have overflowed, so the real extent is
        # ragged. Measuring it lets the block be centred -> equal gaps both sides.
        _r = tmp.canvas.get_renderer(); _p = tmp.text(0, 0, "", fontsize=FS_BULLET)
        _sc = 72.0 / tmp.dpi
        TEXT_W = 0.0
        for blk in WRAPPED:
            for line in blk:
                _p.set_text(line)
                TEXT_W = max(TEXT_W, _p.get_window_extent(renderer=_r).width * _sc)
        _p.remove(); plt.close(tmp)

        n_title = len(BAND_TITLE.split("\n"))
        title_h = n_title * FS_BAND_TITLE * 1.30
        lead = FS_BULLET * 1.45
        H_BAND = (7.0 + title_h + 7.0 + sum(len(w) for w in WRAPPED) * lead
                  + (len(WRAPPED) - 1) * lead * 0.45 + 8.0)

        FIG_H_PT = (M_TOP + H_BAN + G_BC + H_CARD + G_CR + H_RES + G_RB
                    + H_BAND + M_BOT)
        FIG_H = FIG_H_PT / 72.0

        def y(pt):
            return 1.0 - pt / FIG_H_PT

        B1, B0 = y(M_TOP), y(M_TOP + H_BAN)
        C1 = y(M_TOP + H_BAN + G_BC)
        C0 = y(M_TOP + H_BAN + G_BC + H_CARD)
        F1 = y(M_TOP + H_BAN + G_BC + H_CARD + G_CR)
        F0 = y(M_TOP + H_BAN + G_BC + H_CARD + G_CR + H_RES)
        V1 = y(M_TOP + H_BAN + G_BC + H_CARD + G_CR + H_RES + G_RB)
        V0 = y(FIG_H_PT - M_BOT)
        TIERS = [("Control\nvaried", (B0 + B1) / 2),
                 ("Parameter\nswept", (C0 + C1) / 2),
                 ("Principal\nresult", (F0 + F1) / 2),
                 ("Conclusions", (V0 + V1) / 2)]

        fig = plt.figure(figsize=(FIG_W, FIG_H))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
        checks = []

        def rbox(cx, y0, y1, fc, ec, w=BW, lw=0.9):
            ax.add_patch(FancyBboxPatch((cx - w / 2, y0), w, y1 - y0,
                         boxstyle="round,pad=0.002,rounding_size=0.010",
                         linewidth=lw, facecolor=fc, edgecolor=ec, mutation_aspect=1.0))

        def vconn(cx, ytop, ybot, color, lw=0.9):
            ax.add_patch(FancyArrowPatch((cx, ytop), (cx, ybot), arrowstyle="-|>",
                         mutation_scale=6, lw=lw, color=color, shrinkA=0, shrinkB=0))

        for lab, yy in TIERS:
            ax.text(LBL_X, yy, lab, ha="center", va="center", fontsize=FS_LBL,
                    color="#555555", weight="bold", linespacing=1.2)

        for cx, title, card, nrun, find, (ban, tint, ec) in zip(
                CX, TITLE, CARD, NRUN, FIND, COLS):
            rbox(cx, B0, B1, ban, ban)
            t = ax.text(cx, (B0 + B1) / 2, title, ha="center", va="center",
                        color="white", fontsize=FS_BAN, weight="bold")
            checks.append((t, (cx - BW / 2, B0, cx + BW / 2, B1)))

            vconn(cx, B0 - 0.002, C1 + 0.002, ec)
            rbox(cx, C0, C1, "white", ec)
            t = ax.text(cx, C1 - 16.0 / FIG_H_PT, card, ha="center", va="center",
                        color=INK, fontsize=FS_CARD, linespacing=1.4)
            checks.append((t, (cx - BW / 2, C0, cx + BW / 2, C1)))
            t = ax.text(cx, C0 + 7.0 / FIG_H_PT, nrun, ha="center", va="center",
                        color=ec, fontsize=FS_RUN, weight="bold", style="italic")
            checks.append((t, (cx - BW / 2, C0, cx + BW / 2, C1)))

            vconn(cx, C0 - 0.002, F1 + 0.002, ec)
            rbox(cx, F0, F1, tint, ec)
            t = ax.text(cx, (F0 + F1) / 2, find, ha="center", va="center",
                        color=ec, fontsize=FS_RES, weight="bold", linespacing=1.3)
            checks.append((t, (cx - BW / 2, F0, cx + BW / 2, F1)))

            ax.add_patch(FancyArrowPatch((cx, F0 - 0.002), (cx, V1 + 0.002),
                         arrowstyle="-|>", mutation_scale=6, lw=0.7, color=GREY_E,
                         linestyle=(0, (2.5, 1.8))))

        bx0, bx1 = CX[0] - BW / 2, CX[-1] + BW / 2
        ax.add_patch(FancyBboxPatch((bx0, V0), bx1 - bx0, V1 - V0,
                     boxstyle="round,pad=0.002,rounding_size=0.010",
                     linewidth=0.9, facecolor=GREY_F, edgecolor=GREY_E,
                     mutation_aspect=1.0))
        t = ax.text((bx0 + bx1) / 2, V1 - (7.0 + title_h / 2) / FIG_H_PT, BAND_TITLE,
                    ha="center", va="center", fontsize=FS_BAND_TITLE, color=INK,
                    weight="bold", linespacing=1.3)
        checks.append((t, (bx0, V0, bx1, V1)))

        # Centre the block (hanging bullet + widest text line) so the leftover
        # space splits evenly instead of all collecting on the right.
        block_w = BAND_HANG + TEXT_W / (FIG_W * 72.0)
        x0 = bx0 + ((bx1 - bx0) - block_w) / 2.0
        ypt = 7.0 + title_h + 7.0 + lead * 0.75
        for blk in WRAPPED:
            for k, line in enumerate(blk):
                if k == 0:
                    t = ax.text(x0, V1 - ypt / FIG_H_PT, "\u2022", ha="left",
                                va="center", fontsize=FS_BULLET, color=INK)
                    checks.append((t, (bx0, V0, bx1, V1)))
                t = ax.text(x0 + BAND_HANG, V1 - ypt / FIG_H_PT, line, ha="left",
                            va="center", fontsize=FS_BULLET, color=INK)
                checks.append((t, (bx0, V0, bx1, V1)))
                ypt += lead
            ypt += lead * 0.45

        # every text artist must stay inside the box it belongs to; the card and
        # result boxes are FIXED height, so longer wording there overflows and must
        # be reported rather than silently clipped in the compiled PDF
        fig.canvas.draw()
        rr = fig.canvas.get_renderer(); inv = ax.transAxes.inverted()
        bad = []
        for t, (qx0, qy0, qx1, qy1) in checks:
            bb = t.get_window_extent(renderer=rr)
            (a0, b0), (a1, b1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
            if a0 < qx0 - 2e-3 or a1 > qx1 + 2e-3 or b0 < qy0 - 2e-3 or b1 > qy1 + 2e-3:
                bad.append(t.get_text().split("\n")[0][:52])

        os.makedirs(outdir, exist_ok=True)
        fig.savefig(os.path.join(outdir, "fig_study_design.pdf"))
        fig.savefig(os.path.join(outdir, "fig_study_design.png"), dpi=200)
        plt.close(fig)

    if bad:
        print(f"    [study_design] WARNING: {len(bad)} text block(s) overflow:")
        for b in bad:
            print(f"      - {b!r}")
    else:
        print(f"    [study_design] wrote fig_study_design.pdf/.png  "
              f"(all {len(checks)} text blocks fit, band {H_BAND:.0f} pt, "
              f"figure {FIG_H_PT:.0f} pt)")
    return True


def fig_mechanism(outdir):
    """Cell Li concentration against Damkohler number (schematic; no data)."""
    CEILING = 314.0
    cmap = plt.cm.YlOrRd
    C_CHAN, C_POCK = cmap(0.16), cmap(0.62)
    fig, ax = plt.subplots(figsize=(6.6, 4.7))

    # panel (b) concentration vs Damkohler number
    Da = np.logspace(-2, 2, 400)
    C = CEILING * (1.0 - np.exp(-Da))
    ax.semilogx(Da, C, color="0.2", lw=1.8, zorder=4)
    ax.axhline(CEILING, color="0.45", ls=(0, (5, 4)), lw=1.1, zorder=3)
    ax.text(95, CEILING - 9, "Closed-system saturation\nlimit ($\\approx$314 mg L$^{-1}$)",
             ha="right", va="top", fontsize=9.0, color="0.35")
    ax.axvline(1.0, color="0.75", lw=0.9, ls=":", zorder=1)
    # Faint tints separate the two transport-geochemical regimes either side of
    # Da = 1 (drawn behind everything else).
    # Faint tints separate the two regimes either side of Da = 1 (behind all).
    # Faint tints separate the two regimes either side of Da = 1 (behind all).
    ax.axvspan(1e-2, 1.0, color="#4f7fbf", alpha=0.11, zorder=0)
    ax.axvspan(1.0, 1e2, color="#4fae4f", alpha=0.11, zorder=0)
    # Channel-mode concentration window, labelled inside the band itself (the
    # band is curve-free for Da > ~0.1, so the text clears the curve).
    ax.axhspan(1, 27, color=C_CHAN, alpha=0.6, zorder=2)
    ax.text(0.34, 13.5, "channel mode\n1\u201327 mg L$^{-1}$",
             ha="center", va="center", fontsize=8.0, color="0.2", zorder=5)
    # Two representative parcels, named in plain terms; the small legend sits
    # far into the lower-right corner, clear of the curve.
    ax.scatter([0.03], [9], color="0.15", edgecolor=PROD_C, linewidth=1.6,
                s=58, zorder=6, label="Producer output")
    ax.scatter([0.21], [61], color=C_POCK, edgecolor="0.3", s=52, zorder=6,
                label="Bypassed region")
    leg = ax.legend(loc="lower right", fontsize=10, framealpha=0.95,
                     handletextpad=0.4, borderpad=0.6, labelspacing=0.5)
    leg.set_zorder(8)
    # Regime labels, one per half, with the Damkohler condition and its effect.
    ax.text(0.0125, 250, "kinetics / transport-limited\n$\\mathrm{Da}\\ll1$ \u2192 dilute",
             fontsize=9, color="0.42", va="center", zorder=5)
    ax.text(13, 162, "equilibrium-limited\n$\\mathrm{Da}\\gg1$ \u2192 near saturation",
             fontsize=9, color="0.42", ha="center", va="center", zorder=5)
    ax.set_xlim(1e-2, 1e2); ax.set_ylim(0, 342)
    ax.set_xlabel("Damk\u00f6hler number  $\\mathrm{Da}=\\tau_{\\mathrm{res}}/\\tau_{\\mathrm{rxn}}$", fontsize=14)
    ax.set_ylabel("Cell Li concentration  (mg L$^{-1}$)", fontsize=14)
    ax.tick_params(labelsize=12)
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

    fig.tight_layout()
    _save(fig, outdir, "mechanism_schematic", dpi=300)
    print("    [mechanism] wrote mechanism_schematic.pdf/.png")
    return True


# ---------------------------------------------------------------- fig:p32
_BASE_P32_TOTAL = 0.083   # sum of the three family P32 targets at scale 1.0 (fallback only)


def _rate_band(results, base):
    """Base-case produced-concentration range over the rate sweep."""
    from postprocess_field import analyze_run
    g = []
    for r in RATES:
        a = analyze_run(os.path.join(results, f"habanero_field_s{base}_rate{r}"))
        if a is not None:
            g.append(a["avg_li_mg_L"])
    return (min(g), max(g)) if g else (2.0, 37.0)


def _p32_total(results, scale, seed):
    """Total P32 (sum of family targets) from the run or library dfn_summary.json."""
    for cand in (os.path.join(results, f"habanero_field_p{scale}_s{seed}_rate10", "dfn_summary.json"),
                 os.path.join("dfn_library_field", f"habanero_field_p{scale}_s{seed}", "dfn_summary.json")):
        if os.path.exists(cand):
            d = json.load(open(cand))
            v = d.get("p32_values")
            if isinstance(v, (list, tuple)) and v:
                return float(sum(v))
    return None


def fig_p32(results, base, outdir, scales=("0.6", "0.8", "1.0", "1.4", "1.8"), seeds=(42, 101)):
    """fig:p32 -- produced concentration and Gini against fracture intensity at 10 kg/s.
    Reads the P32-sweep runs live; skips cleanly if they are not present."""
    sys.path.insert(0, HERE); sys.path.insert(0, os.getcwd())
    try:
        from postprocess_field import analyze_run
        from channelization import load_run, gini
    except Exception as exc:                                       # noqa: BLE001
        print(f"    [p32] cannot import core modules: {exc}"); return False
    P32, grade, gsd, gtr, grx = [], [], [], [], []
    for sc in scales:
        gg, tt, rr, pp32 = [], [], [], []
        for s in seeds:
            d = os.path.join(results, f"habanero_field_p{sc}_s{s}_rate10")
            a = analyze_run(d)
            if a is not None:
                gg.append(a["avg_li_mg_L"])
            try:
                run = load_run(d)
                tt.append(gini(run["li_mass_per_cell_g"]))
                rr.append(gini(run["spod_lost_mol_per_cell"]))
            except Exception:
                pass
            pv = _p32_total(results, sc, s)
            if pv is not None:
                pp32.append(pv)
        if gg:
            P32.append(np.mean(pp32) if pp32 else float(sc) * _BASE_P32_TOTAL)
            grade.append(np.mean(gg)); gsd.append(np.std(gg, ddof=1) if len(gg) > 1 else 0.0)
            gtr.append(np.mean(tt) if tt else np.nan); grx.append(np.mean(rr) if rr else np.nan)
    if not grade:
        print(f"    [p32] skipped (no P32 runs under {results})"); return False
    P32 = np.array(P32); grade = np.array(grade); gsd = np.array(gsd)
    gtr = np.array(gtr); grx = np.array(grx)
    lo, hi = _rate_band(results, base)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    axA.axhspan(lo, hi, color="#777777", alpha=0.16, lw=0)
    axA.text(P32.max(), hi * 0.96, f"rate-driven range\n({RATES[0]}\u2013{RATES[-1]} kg s$^{{-1}}$)",
             ha="right", va="top", fontsize=7, color="#555555")
    axA.errorbar(P32, grade, yerr=gsd, fmt="o-", color="#1A5FB4", capsize=2.5, ecolor="#1A5FB4",
                 elinewidth=0.9, mfc="white", mec="#1A5FB4", mew=1.1, zorder=3)
    axA.set_xlabel(r"Fracture intensity $P_{32}$  (m$^{-1}$)")
    axA.set_ylabel(r"Produced Li concentration  (mg L$^{-1}$)")
    axA.set_ylim(0, max(hi, grade.max()) * 1.06)
    axA.set_title("(a)  Concentration against intensity at 10 kg s$^{-1}$", loc="left", fontsize=10)
    axB.plot(P32, gtr, "s-", color="#C01C28", mfc="white", mec="#C01C28", mew=1.1, label="Transport (fluid Li)")
    axB.plot(P32, grx, "o-", color="#1A5FB4", mfc="white", mec="#1A5FB4", mew=1.1, label="Reaction (spodumene)")
    axB.set_xlabel(r"Fracture intensity $P_{32}$  (m$^{-1}$)")
    axB.set_ylabel("Gini coefficient"); axB.set_ylim(0.45, 0.80)
    axB.legend(fontsize=8, frameon=False, loc="upper right")
    axB.set_title("(b)  Channelization against intensity", loc="left", fontsize=10)
    _save(fig, outdir, "fig_p32")
    print(f"    [p32] wrote fig_p32.pdf/.png   P32 {P32.min():.3f}-{P32.max():.3f} m^-1  "
          f"concentration {grade.min():.2f}-{grade.max():.2f} mg/L")
    return True


# ---------------------------------------------------------------- fig:ksweep-extended
def _load_kinetics_extended(kin_dirs):
    """Return (kA, concentration) arrays for every sweep case found, or (None, None)."""
    sys.path.insert(0, HERE); sys.path.insert(0, os.getcwd())
    try:
        from postprocess_field import analyze_run
    except Exception as exc:                                       # noqa: BLE001
        print(f"    [kinext] cannot import postprocess_field: {exc}")
        return None, None
    rows = {}
    for dd in kin_dirs:
        for case in sorted(glob.glob(os.path.join(dd, "spd_*"))):
            cj = os.path.join(case, "sweep_case.json")
            if not os.path.exists(cj):
                continue
            c = json.load(open(cj)); a = analyze_run(case)
            if a is None:
                continue
            rows[c["spd_k0_mol_m2_s"] * c["spd_area_cm2_cm3"]] = a["avg_li_mg_L"]
    if not rows:
        print(f"    [kinext] skipped (no sweep_case.json under {', '.join(kin_dirs)})")
        return None, None
    kA = np.array(sorted(rows))
    return kA, np.array([rows[k] for k in kA])


def _draw_kinetics_extended(ax, kA, gr, fs=12, title=True):
    """Draw the produced concentration against k0*A on a supplied axis.

    Returns the saturation limit and the base-case value.
    """
    from matplotlib.ticker import LogLocator
    ceil = float(gr.max()); low_grade, low_kA = gr[0], kA[0]
    base_grade = None
    for k, g in zip(kA, gr):
        if abs(k - 5.0e-11) < 1e-13:
            base_grade = float(g)
    ax.set_xscale("log")
    ax.tick_params(axis="both", labelsize=fs)
    ax.axhline(ceil, color="#2e7d32", ls="--", lw=1.2)
    ax.text(kA.max(), ceil * 1.03, f"Saturation limit $\\approx${ceil:.0f} mg L$^{{-1}}$",
            ha="right", va="bottom", fontsize=12, color="#2e7d32")
    ax.plot(kA, gr, "o-", color="#1A5FB4", mfc="white", mec="#1A5FB4", mew=1.3, ms=6, zorder=3)
    # markers on the curve, with colour-matched labels (no connector lines)
    ax.plot([low_kA], [low_grade], "v", color="#8172B3", ms=9, zorder=4)
    ax.text(2.7e-12, 11, f"Extended low ({low_grade:.2f} mg L$^{{-1}}$)",
            color="#8172B3", fontsize=12, ha="left", va="bottom")
    if base_grade is not None:
        ax.plot([5.0e-11], [base_grade], "s", color="#C01C28", ms=9, zorder=4)
        ax.text(6.0e-11, 2.6, f"Base case ({base_grade:.2f} mg L$^{{-1}}$)",
                color="#C01C28", fontsize=12, ha="left", va="center")
    ax.set_xlabel(r"$k_0\,A$   (mol m$^{-2}$ s$^{-1}$ $\times$ cm$^2$ cm$^{-3}$)", fontsize=fs + 2)
    ax.set_ylabel(r"Produced Li concentration  (mg L$^{-1}$)", fontsize=fs + 2)
    ax.set_ylim(0, ceil * 1.18)
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=9))
    if title:
        ax.set_title("Produced concentration against spodumene rate $\\times$ area",
                     loc="left", fontsize=fs + 2)
    return ceil, base_grade


def fig_kinetics_extended(kin_dirs, outdir):
    """fig:ksweep-extended -- produced concentration against k0*A, standalone."""
    kA, gr = _load_kinetics_extended(kin_dirs)
    if kA is None:
        return False
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ceil, base_grade = _draw_kinetics_extended(ax, kA, gr, fs=12, title=True)
    _save(fig, outdir, "fig_kinetics_extended")
    print(f"    [kinext] wrote fig_kinetics_extended.pdf/.png   "
          f"concentration {gr.min():.2f}-{ceil:.2f} mg/L  base(5e-11)={'n/a' if base_grade is None else f'{base_grade:.2f}'}")
    return True


# ------------------------------------------------- fig:kinetics (combined a, b)
def fig_kinetics_combined(kin_dirs, outdir, csv=None):
    """One figure with both kinetics panels.

    (a) the cases together with the sub-floor extension, collapsed onto the
        product k0*A
    (b) the rate-constant by reactive-area grid, as a heatmap

    The two panels show different things: (a) shows that the concentration depends
    on the product rather than on the two factors separately, and adds the four
    cases below the grid; (b) gives the value at each point of the grid. They are
    drawn together because a reader compares them, and because a single figure
    spares a manuscript one float.

    Falls back to writing nothing if the extension cases cannot be read, in which
    case the two standalone figures remain the ones to use.
    """
    kA, gr = _load_kinetics_extended(kin_dirs)
    if kA is None:
        print("    [kincomb] skipped: the extension cases could not be read")
        return False
    k0s, As, grid, src = _load_ksweep(csv)

    # width ratio follows the standalone sizes, 7.4 inches for the collapse and
    # 5.4 for the grid
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 3.9),
                             gridspec_kw={"width_ratios": [7.4, 5.4], "wspace": 0.28})
    ceil, base_grade = _draw_kinetics_extended(axes[0], kA, gr, fs=10, title=False)
    _draw_ksweep(axes[1], k0s, As, grid, fs=8.5, cb=True, fig=fig)

    for ax, letter in zip(axes, "ab"):
        ax.text(-0.14, 1.06, f"({letter})", transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top", ha="left")
    _save(fig, outdir, "fig_kinetics_combined")
    print(f"    [kincomb] wrote fig_kinetics_combined.pdf/.png   "
          f"(a) concentration {gr.min():.2f}-{ceil:.2f} mg/L across {len(kA)} cases; "
          f"(b) grid {np.nanmin(grid):.1f}-{np.nanmax(grid):.1f} from {src}")
    return True



# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="pflotran_field_results")
    ap.add_argument("--mesh", default="dfn_library_field/habanero_field_s42/full_mesh.inp")
    ap.add_argument("--sweep-csv", default="sweep_kinetics_area_grid.csv")
    ap.add_argument("--kinetics-dirs", nargs="+",
                    default=["pflotran_sweep_kinetics_area", "pflotran_sweep_kinetics_extend"],
                    help="sweep dirs for fig_kinetics_extended (original + sub-floor)")
    ap.add_argument("--base", type=int, default=42)
    ap.add_argument("--only", nargs="+", metavar="TAG",
                    help="build only these figures, e.g. --only design p32 . Tags: "
                         "revised ph co2 ksweep domain calcite mechanism design "
                         "p32 kinext kincomb")
    ap.add_argument("--outdir", default="figs")
    args = ap.parse_args()
    only = set(args.only or [])

    def _want(tag):
        return not only or tag in only
    os.makedirs(args.outdir, exist_ok=True)
    py = sys.executable
    status = {}

    # 1) Four ensemble/base-case data figures via the verified engine.
    #    Empty --mesh/--sweep-csv make generate_revised_figures.py skip its OWN
    #    ksweep/domain (we build those below), so it emits exactly
    #    scaling/channel/bimodal/headline + the revised_numbers.txt dump.
    if _want("revised"):
        status["generate_revised_figures.py  (scaling, channel, bimodal, headline)"] = _run(
        [py, _script("generate_revised_figures.py"),
         "--results", args.results, "--base", str(args.base),
         "--outdir", args.outdir, "--mesh", "", "--sweep-csv", ""],
        "revised-4")

    # 2) Chemistry sweeps (self-contained scripts).
    if _want("ph"):
        status["ph_sweep_summary.py  (fig:ph)"] = _run(
        [py, _script("ph_sweep_summary.py"),
         "--rate-sweep-dir", args.results, "--ph-sweep-dir", args.results,
         "--output-dir", args.outdir], "ph")
    if _want("co2"):
        status["co2_sweep_summary.py  (fig:co2)"] = _run(
        [py, _script("co2_sweep_summary.py"),
         "--plain-dir", args.results, "--co2-dir", args.results,
         "--output-dir", args.outdir], "co2")


    # 3) The four figures built INSIDE this file (no external script needed).
    print("\n>>> [inline figures]")
    if _want("ksweep"):
        status["ksweep_grid          (inline, fig:ksweep)"] = fig_ksweep(args.outdir, csv=args.sweep_csv)
    if _want("domain"):
        status["domain_dfn_proj      (inline, fig:domain)"] = fig_domain(args.mesh, args.outdir)
    if _want("calcite"):
        status["calcite_regime       (inline, fig:calcite)"] = fig_calcite(args.results, args.base, args.outdir)
    if _want("mechanism"):
        status["mechanism_schematic  (inline, fig:mechanism)"] = fig_mechanism(args.outdir)
    if _want("design"):
        status["fig_study_design      (inline, fig:design)"] = fig_study_design(args.outdir)
    if _want("p32"):
        status["fig_p32              (inline, fig:p32)"] = fig_p32(args.results, args.base, args.outdir)
    if _want("kinext"):
        status["fig_kinetics_extended(inline, fig:ksweep-extended)"] = fig_kinetics_extended(args.kinetics_dirs, args.outdir)
    if _want("kincomb"):
        status["fig_kinetics_combined (inline, fig:kinetics)"] = fig_kinetics_combined(
            args.kinetics_dirs, args.outdir, csv=args.sweep_csv)

    expect = ["scaling_loglog", "channelization_composite", "bimodal_composite",
              "master_summary", "ph_sweep_summary", "co2_sweep_summary",
              "ksweep_grid", "domain_dfn_proj", "calcite_regime", "mechanism_schematic",
              "fig_p32", "fig_kinetics_extended", "fig_study_design"]
    print("\n" + "=" * 70 + "\nFIGURE-GENERATION STATUS\n" + "=" * 70)
    for k, ok in status.items():
        print(f"  {'OK ' if ok else 'XX '} {k}")
    print("-" * 70)
    # Credit a figure only if its PDF was written during THIS run. Counting
    # files that merely exist lets stale output from an earlier run mask a
    # total failure -- e.g. invoking the driver without --results, where every
    # data-dependent figure fails but the old PDFs are still on disk.
    def _fresh(n):
        f = os.path.join(args.outdir, n + ".pdf")
        try:
            return os.path.getmtime(f) >= _RUN_T0 - 1.0
        except OSError:
            return False

    if only:                       # judge only what was asked for
        _map = {"revised": ["scaling_loglog", "channelization_composite",
                            "bimodal_composite", "master_summary"],
                "ph": ["ph_sweep_summary"], "co2": ["co2_sweep_summary"],
                "ksweep": ["ksweep_grid"], "domain": ["domain_dfn_proj"],
                "calcite": ["calcite_regime"], "mechanism": ["mechanism_schematic"],
                "design": ["fig_study_design"], "p32": ["fig_p32"],
                "kinext": ["fig_kinetics_extended"],
                "kincomb": ["fig_kinetics_combined"]}
        expect = [n for t in only for n in _map.get(t, [])]
    have = [n for n in expect if _fresh(n)]
    miss = [n for n in expect if n not in have]
    stale = [n for n in miss
             if os.path.exists(os.path.join(args.outdir, n + ".pdf"))]
    failed = [k for k, v in status.items() if not v]

    print(f"  {len(have)}/{len(expect)} manuscript PDFs written by this run "
          f"in '{args.outdir}/'")
    if stale:
        print("  STALE (present but NOT regenerated -- do not use): "
              + ", ".join(m + ".pdf" for m in stale))
    if miss:
        print("  MISSING: " + ", ".join(m + ".pdf" for m in miss))
        print("  (domain_dfn_proj needs --mesh to point at full_mesh.inp;")
        print("   every other figure runs without a mesh.)")
    if failed:
        print("  GENERATORS THAT FAILED: " + ", ".join(sorted(failed)))
    if not miss and not failed:
        print(f"  All {len(expect)} manuscript figures generated. Copy "
              f"'{args.outdir}/*.pdf' into the manuscript figs/ directory.")
        # The combined kinetics figure is an alternative to the pair, not a
        # fourteenth manuscript figure, so it is reported separately rather than
        # counted above. Using it means one \includegraphics in place of two, and
        # one caption covering both panels.
        if _fresh("fig_kinetics_combined"):
            print("  Also written: fig_kinetics_combined.pdf, panels (a) the grid and")
            print("  (b) the collapse onto k0*A. It replaces ksweep_grid.pdf and")
            print("  fig_kinetics_extended.pdf if the manuscript uses a single float.")
    else:
        print("  RUN INCOMPLETE -- fix the errors above and re-run before using figs/.")
        print("  Common cause: paths are relative to the CURRENT directory. From")
        print("  pipeline/, pass --results ../pflotran_field_results and")
        print("  --mesh ../dfn_library_field/habanero_field_s42/full_mesh.inp .")
    print("=" * 70)
    return 0 if (not miss and not failed) else 1


if __name__ == "__main__":
    sys.exit(main())
