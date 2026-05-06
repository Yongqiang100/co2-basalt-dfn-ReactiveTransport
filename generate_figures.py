#!/usr/bin/env python3
"""
generate_figures.py — Single entry point for all manuscript figures.

This script reproduces every figure in:

    Chen, Y., Xie, Q., & Regenauer-Lieb, K. (2026). Fracture network
    connectivity controls on CO2 mineral trapping efficiency in basalt.
    Water Resources Research.

It replaces the previous multi-script workflow (generate_all_figures.py,
generate_new_figures.py, stress_test_analysis.py) with a single command.

Generated figures (all in paper_figures/):

  Fig 1 (a)  fig_dfn_geometry             3D DFN realizations at 4 P32 levels
  Fig 2      fig_timeseries               pH and mineral time series
  Fig 3      fig_carbonate_budget         Carbonate composition by P32
  Fig 4      fig_connectivity_trapping    Trapping vs P32, pH vs P32, CV vs P32
  Fig 5      fig_dissolution              Primary mineral dissolution box plots
  Fig 6      fig_p32_spatial_comparison   3D pH and carbonate at t = 50 yr
  Fig 7      fig_topology_trapping        Intersection density vs trapping
  Fig 8      fig_trapping_efficiency      CO2 trapping efficiency over time
  Fig 9      fig_stagnation_zones         Velocity-decile precipitation

Required inputs (relative to current directory):
  pflotran_results/<dfn_name>/pflotran_co2.h5     PFLOTRAN HDF5 outputs
  pflotran_results/<dfn_name>/simulation_params.json
  dfn_library/<dfn_name>/full_mesh.uge            DFN meshes
  (Fig 1 generates DFNs from seeds; does not need pre-existing meshes.)

Usage:
  python generate_figures.py                       # generate all figures
  python generate_figures.py --only 2 4 8          # only specific figures
  python generate_figures.py --skip-3d             # skip slow 3D figures (1, 6)
  python generate_figures.py --output-dir myfigs   # custom output directory

Dependencies:
  numpy, scipy, h5py, matplotlib
"""

import os
import sys
import glob
import json
import argparse
from typing import Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LogNorm, Normalize, to_rgba
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import stats

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False
    print("WARNING: h5py not installed. Only Fig 1 (DFN geometry) will work.")
    print("Install with: pip install h5py")


# ============================================================
# CONFIGURATION
# ============================================================
RESULTS_ROOT = "pflotran_results"
DFN_ROOT = "dfn_library"
FIG_DIR = "paper_figures"  # overwritten by --output-dir

SEC_PER_YEAR = 3.156e7
CO2_MOLAL = 0.82  # mol/kg of injected CO2

P32_COLORS = {0.75: "#4477AA", 1.0: "#228833", 1.25: "#CCBB44",
              1.5: "#EE6677", 2.0: "#AA3377"}
P32_LABELS = {0.75: "0.75×", 1.0: "1.00×", 1.25: "1.25×",
              1.5: "1.50×", 2.0: "2.00×"}
CARB_COLORS = {"calcite": "#2166AC", "magnesite": "#1A9850",
               "siderite": "#D6394C", "dawsonite": "#F5A623"}
FAM_COLORS = ["#2166AC", "#D6394C", "#1A9850"]
FAM_LABELS = ["Columnar joints", "Cooling fractures", "Conjugate joints"]
MOLAR_VOL = {"Calcite": 3.693e-5, "Magnesite": 2.803e-5,
             "Siderite": 2.938e-5, "Dawsonite": 5.840e-5}

# Publication style — applied at module load
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "axes.linewidth": 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.4,
    "ytick.major.width": 0.4,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "legend.frameon": False,
    "lines.linewidth": 1.0,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "mathtext.default": "regular",
})


# ============================================================
# OUTPUT HELPERS
# ============================================================
def save_fig(fig, name):
    for ext in ("png", "pdf"):
        path = os.path.join(FIG_DIR, f"{name}.{ext}")
        fig.savefig(path, dpi=600 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    print(f"    Saved: {name}.png and {name}.pdf")


def style_3d_ax(ax, x, y, z, margin=1.5):
    ax.set_xlim(x.min() - margin, x.max() + margin)
    ax.set_ylim(y.min() - margin, y.max() + margin)
    ax.set_zlim(z.min() - margin, z.max() + margin)
    ax.view_init(elev=25, azim=-50)
    ax.set_xlabel("x (m)", fontsize=6.5, labelpad=-2)
    ax.set_ylabel("y (m)", fontsize=6.5, labelpad=-2)
    ax.set_zlabel("z (m)", fontsize=6.5, labelpad=-2)
    ax.tick_params(pad=-3, labelsize=5.5)
    for pn in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pn.fill = False
        pn.set_edgecolor((0.85, 0.85, 0.85, 0.3))
    ax.grid(False)


# ============================================================
# HDF5 + MESH I/O HELPERS
# ============================================================
def read_uge_centroids(path):
    xs, ys, zs = [], [], []
    with open(path) as f:
        n = int(f.readline().strip().split()[1])
        for _ in range(n):
            p = f.readline().strip().split()
            xs.append(float(p[1])); ys.append(float(p[2])); zs.append(float(p[3]))
    return np.array(xs), np.array(ys), np.array(zs)


def parse_time_groups(h5f):
    """Return sorted list of (group_key, time_yr) tuples."""
    out = []
    for k in h5f.keys():
        if "Time" not in k or "failure" in k.lower() or "cut" in k.lower():
            continue
        p = k.strip().split()
        try:
            i = p.index("Time"); t = float(p[i + 1])
            t_yr = t if (len(p) > i + 2 and p[i + 2] == "y") else t / SEC_PER_YEAR
            out.append((k, t_yr))
        except Exception:
            continue
    out.sort(key=lambda x: x[1])
    return out


def find_h5_var(grp, prefix):
    for k in grp.keys():
        if k.startswith(prefix):
            return k
    return None


def find_closest_time(tg, target):
    return min(tg, key=lambda x: abs(x[1] - target))


def get_total_carbonate(grp, n):
    """Return per-cell total carbonate VF (raw, not seed-subtracted)."""
    t = np.zeros(n)
    for pf in ("Calcite VF", "Magnesite VF", "Siderite VF", "Dawsonite VF"):
        vk = find_h5_var(grp, pf)
        if vk:
            d = grp[vk][:].flatten()
            t[:min(len(d), n)] += d[:min(len(d), n)]
    return t


def get_velocity_magnitude(h5f, grp_name, n_cells):
    """Extract or estimate velocity magnitude per cell.

    Tries (in order):
    1. Direct velocity variables (Liquid/Darcy Velocity)
    2. Component fluxes (Liquid X/Y/Z-Flux)
    3. Pressure-gradient proxy (last resort)
    """
    grp = h5f[grp_name]
    for prefix in ("Liquid Velocity", "Darcy Velocity", "Liquid X-Velocity"):
        vk = find_h5_var(grp, prefix)
        if vk:
            return np.abs(grp[vk][:].flatten()[:n_cells])
    vx = find_h5_var(grp, "Liquid X-Flux")
    vy = find_h5_var(grp, "Liquid Y-Flux")
    vz = find_h5_var(grp, "Liquid Z-Flux")
    if vx and vy and vz:
        fx = grp[vx][:].flatten()[:n_cells]
        fy = grp[vy][:].flatten()[:n_cells]
        fz = grp[vz][:].flatten()[:n_cells]
        return np.sqrt(fx**2 + fy**2 + fz**2)
    pk = find_h5_var(grp, "Liquid Pressure") or find_h5_var(grp, "Pressure")
    if pk:
        p = grp[pk][:].flatten()[:n_cells]
        dp = np.abs(p - np.mean(p))
        return dp / (np.max(dp) + 1e-30)
    return None


def extract_all_data():
    """Build per-DFN dict of time series for Figs 2–8."""
    D = {}
    for d in sorted(glob.glob(f"{RESULTS_ROOT}/*")):
        nm = os.path.basename(d)
        h5 = os.path.join(d, "pflotran_co2.h5")
        if not os.path.exists(h5):
            continue
        p32 = int(nm.split("_")[1]) / 100.0
        with h5py.File(h5, "r") as f:
            tg = parse_time_groups(f)
            if len(tg) < 2:
                continue
            nc = f[tg[0][0]]["pH"].shape[0]
            ti, pH = [], []
            ca, mg, si, dw = [], [], [], []
            an, fo, di = [], [], []
            for k, t in tg:
                g = f[k]
                ti.append(t)
                pH.append(np.nanmean(g["pH"][:]) if "pH" in g else np.nan)
                for pf, store in [("Calcite VF", ca), ("Magnesite VF", mg),
                                   ("Siderite VF", si), ("Dawsonite VF", dw),
                                   ("Anorthite VF", an), ("Forsterite VF", fo),
                                   ("Diopside VF", di)]:
                    v = find_h5_var(g, pf)
                    store.append(np.nansum(g[v][:]) if v else 0.0)
        tc = np.array(ca) + np.array(mg) + np.array(si) + np.array(dw)
        ft = tc[-1]
        D[nm] = {
            "p32": p32, "ncells": nc, "times": np.array(ti),
            "pH_mean": np.array(pH),
            "calcite": np.array(ca), "magnesite": np.array(mg),
            "siderite": np.array(si), "dawsonite": np.array(dw),
            "anorthite": np.array(an), "forsterite": np.array(fo),
            "diopside": np.array(di),
            "total_carb": tc, "total_carb_per_cell": tc / nc,
            "pct_calc": ca[-1] / ft * 100 if ft > 0 else 0,
            "pct_mag": mg[-1] / ft * 100 if ft > 0 else 0,
            "pct_sid": si[-1] / ft * 100 if ft > 0 else 0,
            "pct_daw": dw[-1] / ft * 100 if ft > 0 else 0,
            "anor_dissolved_pct": (an[0] - an[-1]) / an[0] * 100 if an[0] > 0 else 0,
            "forst_dissolved_pct": (fo[0] - fo[-1]) / fo[0] * 100 if fo[0] > 0 else 0,
            "diop_dissolved_pct": (di[0] - di[-1]) / di[0] * 100 if di[0] > 0 else 0,
        }
    return D


def discover_dfns_with_mesh():
    """List DFNs that have both an HDF5 result and a UGE mesh (for Fig 6)."""
    L = []
    for d in sorted(glob.glob(f"{RESULTS_ROOT}/*")):
        nm = os.path.basename(d)
        h5 = os.path.join(d, "pflotran_co2.h5")
        if not os.path.exists(h5):
            continue
        p32 = int(nm.split("_")[1]) / 100.0
        uge = os.path.join(DFN_ROOT, nm, "full_mesh.uge")
        if not os.path.exists(uge):
            uge = os.path.join(d, "full_mesh.uge")
        if not os.path.exists(uge):
            continue
        with h5py.File(h5, "r") as f:
            tg = parse_time_groups(f)
            if len(tg) < 2:
                continue
            g = f[tg[-1][0]]
            nc = g["pH"].shape[0] if "pH" in g else 0
            tc = get_total_carbonate(g, nc)
            n_precip = int(np.sum(tc > 1e-8))
        L.append({"name": nm, "p32": p32, "h5": h5, "uge": uge,
                  "ncells": nc, "n_precip": n_precip})
    L.sort(key=lambda d: (d["p32"], d["name"]))
    return L


# ============================================================
# DFN GENERATION (used only by Fig 1, no PFLOTRAN run needed)
# ============================================================
def fisher_orientation(t0, p0, k, n, rng):
    xi = rng.uniform(0, 1, n)
    cos_psi = 1 + np.log(xi + (1 - xi) * np.exp(-2 * k)) / k
    cos_psi = np.clip(cos_psi, -1, 1)
    psi = np.arccos(cos_psi)
    chi = rng.uniform(0, 2 * np.pi, n)
    out = np.zeros((n, 3))
    st, ct = np.sin(t0), np.cos(t0); sp, cp = np.sin(p0), np.cos(p0)
    R = np.array([[ct*cp, -sp, st*cp], [ct*sp, cp, st*sp], [-st, 0, ct]])
    for i in range(n):
        sP, cP = np.sin(psi[i]), np.cos(psi[i])
        sc, cc = np.sin(chi[i]), np.cos(chi[i])
        out[i] = R @ np.array([sP * cc, sP * sc, cP])
    return out


def disc_verts(c, n, r, np_=24):
    n = n / np.linalg.norm(n)
    ref = np.array([0, 0, 1]) if abs(n[2]) < 0.9 else np.array([1, 0, 0])
    u = np.cross(n, ref); u /= np.linalg.norm(u); v = np.cross(n, u)
    a = np.linspace(0, 2 * np.pi, np_, endpoint=False)
    return np.array([c + r * (np.cos(t) * u + np.sin(t) * v) for t in a])


def clip_poly(verts, L):
    poly = verts.tolist()
    for ax in range(3):
        for sign, lim in ((1, L), (-1, 0)):
            if len(poly) < 3:
                return np.array(poly) if poly else np.empty((0, 3))
            cl = []
            for i in range(len(poly)):
                a, b = np.array(poly[i]), np.array(poly[(i + 1) % len(poly)])
                ai, bi = sign * a[ax] <= sign * lim, sign * b[ax] <= sign * lim
                if ai:
                    cl.append(a)
                if ai != bi:
                    e = b - a
                    t = (lim - a[ax]) / e[ax] if abs(e[ax]) > 1e-12 else 0
                    cl.append(a + t * e)
            poly = cl
    return np.array(poly) if poly else np.empty((0, 3))


def gen_dfn(pm, L=20, seed=42):
    rng = np.random.default_rng(seed)
    base = [0.8, 0.6, 0.6]
    fams = [(np.pi/2, 0, 20, 1.6, 0.3, 3, 12),
            (0, 0, 15, 1.8, 0.4, 3, 15),
            (np.pi/2, np.pi/2, 20, 1.6, 0.3, 3, 12)]
    out = []
    for fi, (th, ph, kp, lm, ls, rmin, rmax) in enumerate(fams):
        p32 = base[fi] * pm
        mr = np.exp(lm + 0.5 * ls**2)
        n = max(1, int(p32 * L**3 / (np.pi * mr**2)))
        norms = fisher_orientation(th, ph, kp, n, rng)
        rs = np.clip(np.exp(rng.normal(lm, ls, n)), rmin, rmax)
        cs = rng.uniform(0, L, (n, 3))
        for j in range(n):
            v = clip_poly(disc_verts(cs[j], norms[j], rs[j]), L)
            if len(v) >= 3:
                out.append((v, fi))
    return out


# ============================================================
# FIGURE 1 — DFN geometry at four P32 levels
# ============================================================
def fig1_geometry():
    print("\n  Fig 1: DFN geometry...")
    cfgs = [(0.75, 117), (1.00, 42), (1.50, 259), (2.00, 383)]
    L = 20
    AL = [0.55, 0.65, 0.50]
    dfns = [gen_dfn(p, seed=s) for p, s in cfgs]
    subs = [f"P$_{{32}}$ × {p:.2f}  ({len(d)} fractures)"
            for (p, _), d in zip(cfgs, dfns)]

    fig = plt.figure(figsize=(7.2, 7.6), facecolor="white")
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.05],
                          hspace=0.10, wspace=0.02,
                          left=0.01, right=0.99, top=0.97, bottom=0.03)
    for i, (dfn, sub) in enumerate(zip(dfns, subs)):
        r, c = divmod(i, 2)
        ax = fig.add_subplot(gs[r, c], projection="3d", computed_zorder=False)
        cs = np.array([[0,0,0],[L,0,0],[L,L,0],[0,L,0],
                        [0,0,L],[L,0,L],[L,L,L],[0,L,L]], float)
        for a, b in [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
                     (0,4),(1,5),(2,6),(3,7)]:
            ax.plot(*zip(cs[a], cs[b]), color="#444", lw=0.6, alpha=0.5)
        az_r, el_r = np.radians(-52), np.radians(28)
        def dep(v):
            cx, cy, cz = v.mean(0)
            return ((cx*np.cos(az_r) + cy*np.sin(az_r)) * np.cos(el_r)
                    + cz * np.sin(el_r))
        draw = sorted(dfn, key=lambda f: dep(f[0]))
        mx = [60, 80, 100, 120][i]
        if len(draw) > mx:
            idx = np.random.default_rng(0).choice(len(draw), mx, replace=False)
            draw = [draw[j] for j in sorted(idx)]
        ps, fc, ec = [], [], []
        for v, fi in draw:
            ps.append(v)
            fc.append(to_rgba(FAM_COLORS[fi], AL[fi]))
            ec.append(to_rgba("#222", 0.6))
        ax.add_collection3d(Poly3DCollection(
            ps, facecolors=fc, edgecolors=ec,
            linewidths=0.4 if len(draw) < 80 else 0.15, zsort='average'))
        ax.text2D(0.02, 0.96, "abcd"[i], transform=ax.transAxes,
                  fontsize=11, fontweight="bold", va="top",
                  path_effects=[pe.withStroke(linewidth=3, foreground="white")])
        ax.text2D(0.02, 0.87, sub, transform=ax.transAxes,
                  fontsize=7.5, va="top", color="#222",
                  path_effects=[pe.withStroke(linewidth=2, foreground="white")])
        ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
        ax.view_init(elev=28, azim=-52)
        ax.set_xlabel("x (m)", fontsize=7, labelpad=-2)
        ax.set_ylabel("y (m)", fontsize=7, labelpad=-2)
        ax.set_zlabel("z (m)", fontsize=7, labelpad=-2)
        ax.tick_params(pad=-3, labelsize=6)
        for pn in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
            pn.fill = False
            pn.set_edgecolor((0.85, 0.85, 0.85, 0.25))
        ax.grid(False)

    al = fig.add_subplot(gs[2, :])
    al.set_axis_off()
    hs = [Line2D([0], [0], color=FAM_COLORS[i], lw=0, marker='s',
                  markersize=9, markerfacecolor=to_rgba(FAM_COLORS[i], 0.7),
                  markeredgecolor="#333", markeredgewidth=0.4,
                  label=f"Family {i+1}: {FAM_LABELS[i]}") for i in range(3)]
    al.legend(handles=hs, loc="center", ncol=3, frameon=False, fontsize=10)
    save_fig(fig, "fig_dfn_geometry")
    plt.close(fig)


# ============================================================
# FIGURE 2 — Time series
# ============================================================
def fig2_timeseries(data):
    print("\n  Fig 2: Time series...")
    p32v = sorted(set(d["p32"] for d in data.values()))
    grp = {p: [d for d in data.values() if d["p32"] == p] for p in p32v}
    tc = list(data.values())[0]["times"]

    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.8))
    for p in p32v:
        g = grp[p]; c = P32_COLORS[p]
        pH = np.mean([np.interp(tc, d["times"], d["pH_mean"]) for d in g], 0)
        ca = np.mean([np.interp(tc, d["times"], d["calcite"]/d["ncells"]) for d in g], 0)
        mg = np.mean([np.interp(tc, d["times"], d["magnesite"]/d["ncells"]) for d in g], 0)
        an = np.mean([np.interp(tc, d["times"], d["anorthite"]/d["ncells"]) for d in g], 0)
        fo = np.mean([np.interp(tc, d["times"], d["forsterite"]/d["ncells"]) for d in g], 0)
        tb = np.mean([np.interp(tc, d["times"], d["total_carb_per_cell"]) for d in g], 0)
        kw = dict(color=c, lw=1.5)
        axes[0, 0].plot(tc, pH, **kw)
        axes[0, 1].plot(tc, ca, **kw)
        axes[0, 2].plot(tc, mg, **kw)
        axes[1, 0].plot(tc, an, **kw)
        axes[1, 1].plot(tc, fo, **kw)
        axes[1, 2].plot(tc, tb, label=P32_LABELS[p], **kw)
    titles = ["Mean pH", "Calcite VF/cell", "Magnesite VF/cell",
              "Anorthite VF/cell", "Forsterite VF/cell", "Total carbonate VF/cell"]
    yls = ["pH"] + ["Volume fraction"] * 5
    for i, (ax, yl, ti, lb) in enumerate(zip(axes.flat, yls, titles, "abcdef")):
        ax.set_xlabel("")
        if i // 3 == 0:
            ax.set_xticklabels([])
        ax.set_ylabel(yl, fontsize=7)
        ax.set_title(ti, fontsize=8, fontweight="bold", pad=12)
        ax.text(-0.12, 1.15, lb, transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")
    for ax in (axes[0, 1], axes[0, 2], axes[1, 2]):
        ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))
        ax.yaxis.get_offset_text().set_fontsize(6)
    fig.tight_layout(w_pad=1.8, h_pad=2.5, rect=[0, 0.10, 1, 1])
    fig.text(0.5, 0.06, "Time (years)", ha="center", fontsize=8)
    h, l = axes[1, 2].get_legend_handles_labels()
    if axes[1, 2].get_legend():
        axes[1, 2].get_legend().remove()
    fig.legend(h, l, loc="lower center", ncol=5, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, 0.0))
    save_fig(fig, "fig_timeseries")
    plt.close(fig)


# ============================================================
# FIGURE 3 — Carbonate budget
# ============================================================
def fig3_carbonate_budget(data):
    print("\n  Fig 3: Carbonate budget...")
    p32v = sorted(set(d["p32"] for d in data.values()))
    grp = {p: [d for d in data.values() if d["p32"] == p] for p in p32v}

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5))
    x = np.arange(len(p32v)); w = 0.55
    mins = ["calcite", "magnesite", "siderite", "dawsonite"]
    cols = [CARB_COLORS[m] for m in mins]
    lbls = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]
    for ai, (ax, mode) in enumerate(zip(axes, ["abs", "pct"])):
        bot = np.zeros(len(p32v))
        for mn, cl, lb in zip(mins, cols, lbls):
            if mode == "abs":
                v = np.array([np.mean([d[mn][-1] / d["ncells"] for d in grp[p]])
                               for p in p32v])
            else:
                v = np.array([np.mean([d[mn][-1] / d["total_carb"][-1] * 100
                                         if d["total_carb"][-1] > 0 else 0
                                         for d in grp[p]]) for p in p32v])
            ax.bar(x, v, w, bottom=bot, color=cl, label=lb,
                   edgecolor="white", linewidth=0.3)
            bot += v
        ax.set_xticks(x)
        ax.set_xticklabels([P32_LABELS[p] for p in p32v])
        ax.set_ylabel("Carbonate VF per cell" if ai == 0 else "% of total carbonate")
        ax.text(-0.12, 1.15, "ab"[ai], transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")
    fig.text(0.5, 0.02, "P$_{32}$ multiplier", ha="center", fontsize=8)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, -0.05))
    fig.tight_layout(rect=[0, 0.06, 1, 1], w_pad=2.0)
    save_fig(fig, "fig_carbonate_budget")
    plt.close(fig)


# ============================================================
# FIGURE 4 — Connectivity–trapping (3 panels)
# ============================================================
def fig4_connectivity(data):
    print("\n  Fig 4: Connectivity–trapping...")
    p32v = sorted(set(d["p32"] for d in data.values()))
    grp = {p: [d for d in data.values() if d["p32"] == p] for p in p32v}

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.2))

    # Panel a — trapping vs P32
    for d in data.values():
        axes[0].scatter(d["p32"], d["total_carb_per_cell"][-1],
                        color=P32_COLORS[d["p32"]], s=30, alpha=0.8,
                        edgecolor="black", linewidth=0.3, zorder=3)
    ms = [np.mean([d["total_carb_per_cell"][-1] for d in grp[p]]) for p in p32v]
    ss = [np.std([d["total_carb_per_cell"][-1] for d in grp[p]]) for p in p32v]
    axes[0].errorbar(p32v, ms, yerr=ss, fmt="k-o", lw=1, capsize=3,
                     markersize=5, zorder=4)
    xf, yf = np.array(p32v[1:]), np.array(ms[1:])
    r2 = 0
    if len(xf) > 2:
        sl, ic, r, _, _ = stats.linregress(xf, yf); r2 = r**2
        xx = np.linspace(min(xf), max(xf), 50)
        axes[0].plot(xx, sl * xx + ic, "k--", lw=0.7)
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Total carbonate VF/cell")
    axes[0].text(-0.12, 1.15, "a", transform=axes[0].transAxes,
                 fontsize=10, fontweight="bold", va="top")

    # Panel b — pH vs P32
    for d in data.values():
        axes[1].scatter(d["p32"], d["pH_mean"][-1], color=P32_COLORS[d["p32"]],
                        s=30, alpha=0.8, edgecolor="black", linewidth=0.3, zorder=3)
    pm = [np.mean([d["pH_mean"][-1] for d in grp[p]]) for p in p32v]
    ps = [np.std([d["pH_mean"][-1] for d in grp[p]]) for p in p32v]
    axes[1].errorbar(p32v, pm, yerr=ps, fmt="k-o", lw=1, capsize=3,
                     markersize=5, zorder=4)
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Mean pH at t = 50 yr")
    axes[1].text(-0.12, 1.15, "b", transform=axes[1].transAxes,
                 fontsize=10, fontweight="bold", va="top")

    # Panel c — CV vs P32
    cvs = []
    for p in p32v:
        means = np.mean([d["total_carb_per_cell"][-1] for d in grp[p]])
        stds = np.std([d["total_carb_per_cell"][-1] for d in grp[p]])
        cvs.append(stds / means * 100 if means > 0 else 0)
    bars = axes[2].bar(range(len(p32v)), cvs,
                        color=[P32_COLORS[p] for p in p32v],
                        edgecolor="black", linewidth=0.3, width=0.6)
    for b, cv in zip(bars, cvs):
        axes[2].text(b.get_x() + b.get_width()/2, b.get_height() + 2,
                     f"{cv:.0f}%", ha="center", va="bottom", fontsize=6)
    axes[2].set_xticks(range(len(p32v)))
    axes[2].set_xticklabels([P32_LABELS[p] for p in p32v])
    axes[2].set_xlabel("")
    axes[2].set_ylabel("CV (%)")
    axes[2].text(-0.12, 1.15, "c", transform=axes[2].transAxes,
                 fontsize=10, fontweight="bold", va="top")

    fig.tight_layout(w_pad=1.5, rect=[0, 0.12, 1, 1])
    fig.text(0.5, 0.08, "P$_{32}$ multiplier", ha="center", fontsize=8)
    hs = [Line2D([0], [0], color=P32_COLORS[p], lw=0, marker="o", markersize=5,
                  markerfacecolor=P32_COLORS[p], markeredgecolor="black",
                  markeredgewidth=0.3, label=f"P$_{{32}}$ {P32_LABELS[p]}")
          for p in p32v]
    hs.append(Line2D([0], [0], color="black", lw=0.9, marker="o",
                     markersize=4, markerfacecolor="black", label="Mean ± s.d."))
    if r2 > 0:
        hs.append(Line2D([0], [0], color="black", ls="--", lw=0.7,
                         label=f"R² = {r2:.3f}"))
    fig.legend(handles=hs, loc="lower center", ncol=len(hs), frameon=False,
               fontsize=6.5, bbox_to_anchor=(0.5, 0.0))
    save_fig(fig, "fig_connectivity_trapping")
    plt.close(fig)


# ============================================================
# FIGURE 5 — Dissolution box plots
# ============================================================
def fig5_dissolution(data):
    print("\n  Fig 5: Dissolution box plots...")
    p32v = sorted(set(d["p32"] for d in data.values()))
    grp = {p: [d for d in data.values() if d["p32"] == p] for p in p32v}

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.2))
    for ax, (key, title), lb in zip(
            axes,
            [("anor_dissolved_pct", "Anorthite"),
             ("forst_dissolved_pct", "Forsterite"),
             ("diop_dissolved_pct", "Diopside")],
            "abc"):
        vals = [[d[key] for d in grp[p]] for p in p32v]
        bp = ax.boxplot(vals, tick_labels=[P32_LABELS[p] for p in p32v],
                         patch_artist=True, widths=0.5,
                         medianprops=dict(color="black", lw=0.8),
                         whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                         flierprops=dict(markersize=3))
        for patch, p in zip(bp["boxes"], p32v):
            patch.set_facecolor(P32_COLORS[p])
            patch.set_alpha(0.4)
            patch.set_edgecolor("black")
            patch.set_linewidth(0.5)
        for j, (p, v) in enumerate(zip(p32v, vals)):
            jit = np.random.default_rng(42).uniform(-0.12, 0.12, len(v))
            ax.scatter(np.array([j+1]*len(v)) + jit, v,
                       color=P32_COLORS[p], s=22, edgecolor="black",
                       linewidth=0.3, zorder=3, alpha=0.85)
        ax.set_ylabel("Dissolved (%)")
        ax.set_title(title, fontsize=8, fontweight="bold", pad=12)
        ax.set_xlabel("")
        ax.text(-0.12, 1.15, lb, transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top")
    fig.tight_layout(w_pad=1.5, rect=[0, 0.12, 1, 1])
    fig.text(0.5, 0.08, "P$_{32}$ multiplier", ha="center", fontsize=8)
    hs = [Line2D([0], [0], color=P32_COLORS[p], lw=0, marker="o", markersize=5,
                  markerfacecolor=P32_COLORS[p], markeredgecolor="black",
                  markeredgewidth=0.3, label=f"P$_{{32}}$ {P32_LABELS[p]}")
          for p in p32v]
    fig.legend(handles=hs, loc="lower center", ncol=5, frameon=False,
               fontsize=7, bbox_to_anchor=(0.5, 0.0))
    save_fig(fig, "fig_dissolution")
    plt.close(fig)


# ============================================================
# FIGURE 6 — 3D spatial comparison
# ============================================================
def fig6_spatial(dfn_list):
    print("\n  Fig 6: 3D spatial comparison...")
    targets = [0.75, 1.0, 1.5, 2.0]
    reps = {}
    for d in dfn_list:
        if d["p32"] in targets:
            if d["p32"] not in reps or d["n_precip"] > reps[d["p32"]]["n_precip"]:
                reps[d["p32"]] = d
    sel = [reps[p] for p in targets if p in reps]
    nc = len(sel)
    if nc < 2:
        print("    Not enough DFNs with meshes — skipping.")
        return

    fig = plt.figure(figsize=(11, 7.5), facecolor="white")
    gs = fig.add_gridspec(3, nc, height_ratios=[1, 1, 0.06], hspace=-0.25,
                          wspace=0.08, left=0.01, right=0.99, top=0.97,
                          bottom=0.04)
    sc_ph = sc_carb = None
    for ci, dfn in enumerate(sel):
        x, y, z = read_uge_centroids(dfn["uge"])
        n = dfn["ncells"]
        sz = max(0.3, 12000 / n)
        with h5py.File(dfn["h5"], "r") as f:
            tg = parse_time_groups(f)
            tk, _ = find_closest_time(tg, 50.0)
            g = f[tk]
            ph = g["pH"][:].flatten()[:n]
            tc = get_total_carbonate(g, n)

        ax = fig.add_subplot(gs[0, ci], projection="3d", computed_zorder=False)
        sc_ph = ax.scatter(x[:n], y[:n], z[:n], c=ph, cmap="RdYlBu",
                            norm=Normalize(3, 7.5), s=sz, alpha=0.7,
                            edgecolors="none", rasterized=True)
        ax.text2D(0.02, 0.98, "abcd"[ci], transform=ax.transAxes,
                  fontsize=11, fontweight="bold", va="top",
                  path_effects=[pe.withStroke(linewidth=3, foreground="white")])
        ax.text2D(0.02, 0.88, f"P$_{{32}}$ × {dfn['p32']:.2f}",
                  transform=ax.transAxes, fontsize=8, va="top", color="#222",
                  path_effects=[pe.withStroke(linewidth=2, foreground="white")])
        style_3d_ax(ax, x, y, z)

        ax2 = fig.add_subplot(gs[1, ci], projection="3d", computed_zorder=False)
        tc_plot = np.copy(tc); tc_plot[tc_plot < 1e-12] = 1e-12
        vm = np.percentile(tc_plot[tc > 1e-8], 99) if np.any(tc > 1e-8) else 1e-5
        carb_norm = LogNorm(1e-10, max(vm, 1e-4))
        pt_sz = max(1.0, 14000 / n)
        sc_carb = ax2.scatter(x, y, z, c=tc_plot, cmap="RdYlBu_r",
                               norm=carb_norm, s=pt_sz, alpha=0.65,
                               edgecolors="none", linewidths=0,
                               rasterized=True, zorder=2)
        mask_hi = tc > 1e-6
        if np.any(mask_hi):
            ax2.scatter(x[mask_hi], y[mask_hi], z[mask_hi], c=tc[mask_hi],
                        cmap="RdYlBu_r", norm=carb_norm,
                        s=max(2.0, 22000/n), alpha=0.95,
                        edgecolors="black", linewidths=0.15,
                        rasterized=True, zorder=5)
        ax2.text2D(0.02, 0.98, "efgh"[ci], transform=ax2.transAxes,
                   fontsize=11, fontweight="bold", va="top",
                   path_effects=[pe.withStroke(linewidth=3, foreground="white")])
        ax2.text2D(0.02, 0.88, f"P$_{{32}}$ × {dfn['p32']:.2f}",
                   transform=ax2.transAxes, fontsize=8, va="top", color="#222",
                   path_effects=[pe.withStroke(linewidth=2, foreground="white")])
        style_3d_ax(ax2, x, y, z)

    gs_bot = gs[2, :].subgridspec(1, 5, width_ratios=[1.5, 0.15, 1.5, 0.15, 1],
                                    wspace=0.1)
    if sc_ph:
        cax = fig.add_subplot(gs_bot[0, 0])
        cb = fig.colorbar(sc_ph, cax=cax, orientation="horizontal")
        cb.set_label("pH (top row, a–d)", fontsize=7.5)
        cb.ax.tick_params(labelsize=6.5)
    if sc_carb:
        cax = fig.add_subplot(gs_bot[0, 2])
        cb = fig.colorbar(sc_carb, cax=cax, orientation="horizontal")
        cb.set_label("Total carbonate VF (bottom row, e–h)", fontsize=7.5)
        cb.ax.tick_params(labelsize=6.5)
    at = fig.add_subplot(gs_bot[0, 4])
    at.set_axis_off()
    at.text(0.5, 0.5, "All panels at t = 50 yr", transform=at.transAxes,
            fontsize=9, va="center", ha="center", fontstyle="italic", color="#444")
    save_fig(fig, "fig_p32_spatial_comparison")
    plt.close(fig)


# ============================================================
# FIGURE 7 — Topology–trapping (intersection density)
# ============================================================
def _count_intersections_from_dfnworks(dfn_dir):
    candidates = [
        os.path.join(dfn_dir, "intersection_list.dat"),
        os.path.join(dfn_dir, "intersections.dat"),
        os.path.join(dfn_dir, "dfnGen_output", "intersection_list.dat"),
    ]
    candidates += glob.glob(os.path.join(dfn_dir, "intersect*.dat"))
    candidates += glob.glob(os.path.join(dfn_dir, "intersect*.inp"))
    for path in candidates:
        if not os.path.exists(path):
            continue
        try:
            with open(path) as f:
                lines = f.readlines()
            count = sum(1 for line in lines
                        if line.strip() and not line.startswith("#")
                        and not line.startswith("//"))
            first = lines[0].strip().split()
            if len(first) == 1 and first[0].isdigit():
                count = int(first[0])
            if count > 0:
                return count
        except Exception:
            continue
    return None


def _count_intersections_from_inp(dfn_dir):
    inp_path = os.path.join(dfn_dir, "full_mesh.inp")
    if not os.path.exists(inp_path):
        return None
    try:
        with open(inp_path) as f:
            header = f.readline().strip().split()
            n_nodes, n_cells = int(header[0]), int(header[1])
            for _ in range(n_nodes):
                f.readline()
            node_materials = {}
            for _ in range(n_cells):
                parts = f.readline().strip().split()
                mat_id = int(parts[1])
                for nid in [int(x) for x in parts[3:]]:
                    if nid not in node_materials:
                        node_materials[nid] = set()
                    node_materials[nid].add(mat_id)
        return sum(1 for mats in node_materials.values() if len(mats) >= 2)
    except Exception:
        return None


def _count_intersections_from_params(dfn_dir, results_dir):
    for path in [os.path.join(results_dir, "simulation_params.json"),
                 os.path.join(dfn_dir, "params.json")]:
        if not os.path.exists(path):
            continue
        try:
            with open(path) as f:
                params = json.load(f)
            for key in ("n_intersections", "num_intersections", "intersection_count"):
                if key in params:
                    return int(params[key])
        except Exception:
            continue
    return None


def _count_intersections(dfn_name):
    dfn_dir = os.path.join(DFN_ROOT, dfn_name)
    results_dir = os.path.join(RESULTS_ROOT, dfn_name)
    for func, label in [
        (lambda: _count_intersections_from_dfnworks(dfn_dir), "dfnWorks file"),
        (lambda: _count_intersections_from_inp(dfn_dir), "LaGriT .inp"),
        (lambda: _count_intersections_from_params(dfn_dir, results_dir), "params.json")]:
        n = func()
        if n is not None:
            return n, label
    return None, None


def fig7_topology():
    print("\n  Fig 7: Topology–trapping (intersection density)...")
    records, methods_used = [], set()
    for d in sorted(glob.glob(f"{RESULTS_ROOT}/*")):
        nm = os.path.basename(d)
        h5 = os.path.join(d, "pflotran_co2.h5")
        pj = os.path.join(d, "simulation_params.json")
        if not os.path.exists(h5):
            continue
        p32 = int(nm.split("_")[1]) / 100.0
        domain_vol = 8000.0
        if os.path.exists(pj):
            with open(pj) as f:
                params = json.load(f)
            domain_vol = params.get("domain_volume_m3", domain_vol)
        n_int, method = _count_intersections(nm)
        if n_int is None:
            continue
        methods_used.add(method)
        with h5py.File(h5, "r") as f:
            tg = parse_time_groups(f)
            if len(tg) < 2:
                continue
            g = f[tg[-1][0]]
            nc = g["pH"].shape[0] if "pH" in g else 0
            if nc == 0:
                continue
            tc = get_total_carbonate(g, nc)
            carb_per_cell = np.sum(tc) / nc
        records.append({"name": nm, "p32": p32, "n_intersections": n_int,
                        "int_density": n_int / domain_vol,
                        "carb_per_cell": carb_per_cell})
    if len(records) < 3:
        print(f"    Only {len(records)} DFNs with intersection data — skipping.")
        return
    print(f"    {len(records)} DFNs, methods: {', '.join(methods_used)}")

    int_dens = np.array([r["int_density"] for r in records])
    carb_pc = np.array([r["carb_per_cell"] for r in records])
    rho_s, p_s = stats.spearmanr(int_dens, carb_pc)
    r_p, _ = stats.pearsonr(int_dens, carb_pc)
    r2 = r_p ** 2
    slope, intercept, _, _, _ = stats.linregress(int_dens, carb_pc)
    print(f"    Spearman ρ = {rho_s:.3f} (p = {p_s:.4f}), R² = {r2:.3f}")

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 4.0))
    for r in records:
        ax.scatter(r["int_density"], r["carb_per_cell"],
                   color=P32_COLORS[r["p32"]], s=35, alpha=0.85,
                   edgecolor="black", linewidth=0.3, zorder=3)
    x_fit = np.linspace(int_dens.min(), int_dens.max(), 100)
    ax.plot(x_fit, slope * x_fit + intercept, "k--", lw=0.8, zorder=2)
    ax.text(0.05, 0.95, f"R² = {r2:.3f}\nSpearman ρ = {rho_s:.3f}",
            transform=ax.transAxes, fontsize=7, va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#ccc", alpha=0.9))
    ax.set_xlabel("Intersection density (m$^{-3}$)")
    ax.set_ylabel("Total carbonate VF per cell")
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_fontsize(6)

    p32v = sorted(set(r["p32"] for r in records))
    handles = [Line2D([0], [0], color=P32_COLORS[p], lw=0, marker="o",
                       markersize=5.5, markerfacecolor=P32_COLORS[p],
                       markeredgecolor="black", markeredgewidth=0.3,
                       label=f"P$_{{32}}$ {P32_LABELS[p]}") for p in p32v]
    handles.append(Line2D([0], [0], color="black", ls="--", lw=0.8,
                          label=f"Linear fit (R² = {r2:.3f})"))
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=6.5, handletextpad=0.3, columnspacing=1.0,
               bbox_to_anchor=(0.55, 0.02))
    save_fig(fig, "fig_topology_trapping")
    plt.close(fig)


# ============================================================
# FIGURE 8 — Trapping efficiency
# ============================================================
def _cum_co2_inj(times, rate):
    """Mid-point integration of injection rate accounting for ramp-up."""
    rt = np.array([0, 1e-4, 1e-2])
    rr = np.array([0, 0.1 * rate, rate])
    c = np.zeros(len(times))
    for i in range(1, len(times)):
        tm = 0.5 * (times[i - 1] + times[i])
        r = np.interp(tm, rt, rr) if tm < 1e-2 else rate
        c[i] = c[i - 1] + r * (times[i] - times[i - 1]) * SEC_PER_YEAR * CO2_MOLAL
    return c


def _co2_in_carb(grp, nc):
    """Net CO2 in carbonates (mol), seed-subtracted."""
    t = 0.0
    for mn, mv in MOLAR_VOL.items():
        vk = find_h5_var(grp, f"{mn} VF")
        if vk:
            t += max(np.nansum(grp[vk][:]) - 1e-6 * nc, 0) / mv
    return t


def fig8_efficiency():
    print("\n  Fig 8: Trapping efficiency...")
    recs = []
    for d in sorted(glob.glob(f"{RESULTS_ROOT}/*")):
        nm = os.path.basename(d)
        h5 = os.path.join(d, "pflotran_co2.h5")
        pj = os.path.join(d, "simulation_params.json")
        if not os.path.exists(h5):
            continue
        p32 = int(nm.split("_")[1]) / 100.0
        rate = 0.01
        if os.path.exists(pj):
            with open(pj) as f:
                pr = json.load(f)
            rate = pr.get("water_rate_kg_s_scaled", pr.get("water_rate_kg_s", 0.01))
        with h5py.File(h5, "r") as f:
            tg = parse_time_groups(f)
            if len(tg) < 3:
                continue
            nc = f[tg[0][0]]["pH"].shape[0]
            ti = np.array([t for _, t in tg])
            inj = _cum_co2_inj(ti, rate)
            cb = np.array([_co2_in_carb(f[k], nc) for k, _ in tg])
        eff = np.where(inj > 0, cb / inj * 100, 0.0)
        recs.append({"name": nm, "p32": p32, "times": ti,
                     "eff": eff, "final": eff[-1]})
    if not recs:
        return
    p32v = sorted(set(r["p32"] for r in recs))
    grp = {p: [r for r in recs if r["p32"] == p] for p in p32v}
    tc = recs[0]["times"]

    fig = plt.figure(figsize=(7.2, 3.5), facecolor="white")
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.04],
                          width_ratios=[1.4, 1], hspace=0.30, wspace=0.35,
                          left=0.08, right=0.97, top=0.96, bottom=0.02)
    ax = fig.add_subplot(gs[0, 0])
    for p in p32v:
        ae = np.mean([np.interp(tc, r["times"], r["eff"]) for r in grp[p]], 0)
        ax.plot(tc, ae, color=P32_COLORS[p], lw=1.8)
    ax.axvline(x=2, color="#666", ls="--", lw=0.6)
    ax.text(2, 1.03, "t = 2 yr", fontsize=6, color="#666",
            ha="center", va="bottom", transform=ax.get_xaxis_transform())
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("CO$_2$ trapping efficiency (%)")
    ax.set_xlim(0, 52)
    ax.text(-0.12, 1.06, "a", transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="top")

    ax2 = fig.add_subplot(gs[0, 1])
    rng = np.random.default_rng(42)
    for p in p32v:
        vs = [r["final"] for r in grp[p]]
        jit = rng.uniform(-0.04, 0.04, len(vs))
        ax2.scatter(np.array([p]*len(vs)) + jit, vs,
                    color=P32_COLORS[p], s=28, alpha=0.75,
                    edgecolor="black", linewidth=0.3, zorder=3)
    ms = [np.mean([r["final"] for r in grp[p]]) for p in p32v]
    ss = [np.std([r["final"] for r in grp[p]]) for p in p32v]
    ax2.errorbar(p32v, ms, yerr=ss, fmt="k-o", lw=0.9, capsize=3,
                 markersize=4.5, zorder=4, capthick=0.6)
    ax2.set_xlabel("P$_{32}$ multiplier")
    ax2.set_ylabel("Efficiency at t = 50 yr (%)")
    ax2.set_xticks(p32v)
    ax2.set_xticklabels([P32_LABELS[p] for p in p32v])
    ax2.text(-0.12, 1.06, "b", transform=ax2.transAxes,
             fontsize=10, fontweight="bold", va="top")

    al = fig.add_subplot(gs[1, :])
    al.set_axis_off()
    hs = []
    for p in p32v:
        hs.append(Line2D([0], [0], color=P32_COLORS[p], lw=2,
                         marker="none", label=P32_LABELS[p]))
    for p in p32v:
        hs.append(Line2D([0], [0], color=P32_COLORS[p], lw=0,
                         marker="o", markersize=4.5,
                         markerfacecolor=P32_COLORS[p],
                         markeredgecolor="black", markeredgewidth=0.3,
                         label=P32_LABELS[p]))
    hs.append(Line2D([0], [0], color="black", lw=0.9, marker="o",
                     markersize=4, markerfacecolor="black", label="Mean ± s.d."))
    al.legend(handles=hs, loc="center", ncol=11, frameon=False,
              fontsize=7, handlelength=1.4, handletextpad=0.3,
              columnspacing=0.8)
    save_fig(fig, "fig_trapping_efficiency")
    plt.close(fig)


# ============================================================
# FIGURE 9 — Stagnation zones (velocity-decile precipitation)
# ============================================================
def fig9_stagnation():
    print("\n  Fig 9: Stagnation zones (velocity-decile precipitation)...")
    decile_records = []
    for d in sorted(glob.glob(f"{RESULTS_ROOT}/*")):
        nm = os.path.basename(d)
        h5 = os.path.join(d, "pflotran_co2.h5")
        if not os.path.exists(h5):
            continue
        with h5py.File(h5, "r") as f:
            tg = parse_time_groups(f)
            if len(tg) < 2:
                continue
            grp_name = tg[-1][0]
            grp = f[grp_name]
            nc = grp["pH"].shape[0] if "pH" in grp else 0
            if nc == 0:
                continue
            tc = get_total_carbonate(grp, nc)
            vel = get_velocity_magnitude(f, grp_name, nc)
        if vel is None:
            print(f"    No velocity data for {nm}, skipping.")
            continue

        carb_total = np.sum(tc)
        if carb_total < 1e-15:
            continue  # no precipitation — exclude from average

        decile_fracs = []
        for lo_pct, hi_pct in [(0,10), (10,20), (20,30), (30,40), (40,50),
                                (50,60), (60,70), (70,80), (80,90), (90,100)]:
            lo = np.percentile(vel, lo_pct)
            hi = np.percentile(vel, hi_pct)
            if lo_pct == 0:
                mask = vel <= hi
            elif hi_pct == 100:
                mask = vel > lo
            else:
                mask = (vel > lo) & (vel <= hi)
            decile_fracs.append(np.sum(tc[mask]) / carb_total)
        decile_records.append(decile_fracs)

    if not decile_records:
        print("    No DFNs with velocity data and precipitation — skipping.")
        return

    mean_deciles = np.mean(decile_records, axis=0)
    n_active = len(decile_records)
    print(f"    Averaged over {n_active} DFNs with measurable precipitation")

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.5))
    decile_labels = [f"{i*10}–{(i+1)*10}" for i in range(10)]
    ax.bar(range(10), mean_deciles * 100, color="#4477AA",
           edgecolor="black", linewidth=0.3, width=0.7)
    ax.set_xticks(range(10))
    ax.set_xticklabels(decile_labels, rotation=45, ha="right", fontsize=5.5)
    ax.set_xlabel(r"Flow velocity percentile (low $\rightarrow$ high)")
    ax.set_ylabel("Carbonate fraction (%)")
    ax.axhline(y=10, color="#999", ls=":", lw=0.5, label="Uniform (10%)")
    ax.legend(fontsize=6, loc="upper left")
    ax.text(0.02, 0.88,
            f"n = {n_active} DFNs\nwith precipitation",
            transform=ax.transAxes, fontsize=6, ha="left", va="top",
            bbox=dict(boxstyle="round", facecolor="white",
                      edgecolor="#ccc", alpha=0.9))
    fig.tight_layout()
    save_fig(fig, "fig_stagnation_zones")
    plt.close(fig)


# ============================================================
# MAIN
# ============================================================
ALL_FIGURES = [
    (1, "fig_dfn_geometry",            "DFN geometry",                     fig1_geometry,    False),
    (2, "fig_timeseries",              "pH and mineral time series",       fig2_timeseries,  True),
    (3, "fig_carbonate_budget",        "Carbonate budget",                  fig3_carbonate_budget, True),
    (4, "fig_connectivity_trapping",   "Connectivity–trapping (3 panels)",  fig4_connectivity, True),
    (5, "fig_dissolution",             "Dissolution box plots",             fig5_dissolution,  True),
    (6, "fig_p32_spatial_comparison",  "3D spatial pH and carbonate",       fig6_spatial,      "mesh"),
    (7, "fig_topology_trapping",       "Intersection density vs trapping",  fig7_topology,     False),
    (8, "fig_trapping_efficiency",     "CO2 trapping efficiency",           fig8_efficiency,   False),
    (9, "fig_stagnation_zones",        "Velocity-decile precipitation",     fig9_stagnation,   False),
]


def main():
    global FIG_DIR, RESULTS_ROOT, DFN_ROOT

    parser = argparse.ArgumentParser(
        description="Generate all manuscript figures.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Figures generated in paper_figures/:\n" +
            "\n".join(f"  Fig {n}: {fname:<35s} {desc}"
                       for n, fname, desc, _, _ in ALL_FIGURES) +
            "\n"
        ))
    parser.add_argument("--only", type=int, nargs="+", default=None,
                        help="Generate only specified figures (e.g., --only 2 4 8)")
    parser.add_argument("--skip-3d", action="store_true",
                        help="Skip slow 3D figures (Fig 1 and Fig 6)")
    parser.add_argument("--output-dir", type=str, default=FIG_DIR,
                        help=f"Output directory (default: {FIG_DIR})")
    parser.add_argument("--results-dir", type=str, default=RESULTS_ROOT,
                        help=f"PFLOTRAN results directory (default: {RESULTS_ROOT})")
    parser.add_argument("--dfn-dir", type=str, default=DFN_ROOT,
                        help=f"DFN library directory (default: {DFN_ROOT})")
    args = parser.parse_args()

    # Apply CLI overrides
    FIG_DIR = args.output_dir
    RESULTS_ROOT = args.results_dir
    DFN_ROOT = args.dfn_dir
    os.makedirs(FIG_DIR, exist_ok=True)

    # Decide which figures to run
    requested = set(args.only) if args.only else set(n for n, *_ in ALL_FIGURES)
    skip_3d = {1, 6} if args.skip_3d else set()
    to_run = [(n, fn, desc, func, needs)
              for n, fn, desc, func, needs in ALL_FIGURES
              if n in requested and n not in skip_3d]

    if not to_run:
        print("No figures selected. Exiting.")
        return

    print("="*60)
    print(f"Generating {len(to_run)} figure(s) → {FIG_DIR}/")
    print("="*60)

    needs_data = any(needs is True for _, _, _, _, needs in to_run)
    needs_mesh = any(needs == "mesh" for _, _, _, _, needs in to_run)

    data = None
    dfn_list = None

    if needs_data or needs_mesh:
        if not HAS_H5PY:
            print("\nERROR: h5py is required for figures 2–9.")
            print("Install with: pip install h5py")
            sys.exit(1)
        if not os.path.isdir(RESULTS_ROOT):
            print(f"\nERROR: results directory not found: {RESULTS_ROOT}")
            sys.exit(1)

    if needs_data:
        print(f"\n  Loading data from {RESULTS_ROOT}/...")
        data = extract_all_data()
        print(f"  Loaded {len(data)} DFN realizations")
        if len(data) == 0:
            print("\nERROR: no DFN data found. Check --results-dir.")
            sys.exit(1)

    if needs_mesh:
        print(f"\n  Discovering DFNs with meshes...")
        dfn_list = discover_dfns_with_mesh()
        print(f"  Found {len(dfn_list)} DFN(s) with both HDF5 and mesh")

    # Run each requested figure
    for n, fname, desc, func, needs in to_run:
        if needs is True:
            func(data)
        elif needs == "mesh":
            if dfn_list and len(dfn_list) >= 2:
                func(dfn_list)
            else:
                print(f"\n  Fig {n} skipped: needs mesh files in {DFN_ROOT}/")
        else:
            func()

    print("\n" + "="*60)
    print(f"Done. {len(to_run)} figure(s) saved to {FIG_DIR}/")
    print("="*60)


if __name__ == "__main__":
    main()
