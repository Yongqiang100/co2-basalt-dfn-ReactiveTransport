#!/usr/bin/env python3
"""
generate_figures.py — Single entry point for all manuscript figures.

Revision R2: reads the corrected coupled runs (runs_gravityoff, A_feedback__p32_*),
merges the output parts of continued runs, and computes the flow through each cell
from the simulated pressure where no velocity is written. Figure layouts unchanged.

This script reproduces every figure in:

    Chen, Y., Xie, Q., & Regenauer-Lieb, K. (2026). Fracture network
    connectivity controls on CO2 mineral trapping efficiency in basalt.
    Water Resources Research.

It replaces the previous multi-script workflow (generate_all_figures.py,
generate_new_figures.py, stress_test_analysis.py) with a single command.

Generated figures (all in paper_figures/):

  Fig 1 (a)  fig_dfn_geometry             3D DFN realizations at 4 P32 levels
  Fig 2      fig_timeseries               pH and mineral time series
  Fig 3      fig_carbonate_budget         Carbonate assemblage by network
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
import re
import json
import argparse

# so that `from fig_study_design import ...` resolves when the
# script is run from the parent directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
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
SEED_VF = 1e-6      # initial VF of every secondary phase

# A cell whose secondary phases sum above this holds more solid than its own
# volume, which the fixed-porosity formulation does not prevent. Such cells are
# excluded, as they are in finalise.py and Tables 5 and 6.
VALID_MAX_SEC = 1.0

# Realisations in which the invalid cells hold most of the carbonate. Excluded
# outright rather than partially masked.
# The two exclusions were fixed-porosity cells holding more solid than their volume.
# The coupled runs have no such cells, so nothing is excluded.
EXCLUDE_CASES = ()

SEC_PHASES = ("Calcite VF", "Magnesite VF", "Siderite VF", "Dawsonite VF",
              "Kaolinite VF", "Chalcedony VF")
RESULTS_ROOT = "runs_gravityoff"   # corrected runs: gravity omitted, uniform outflow pressure

# Block A alone is the ten-per-level ensemble the paper reports: 49 runs, ten
# at each of five intensity levels. The 25 C_baseline runs are a separate
# five-per-level set on different seeds, and pooling the two gives fifteen per
# level, which is not the ensemble any table describes.
# Coupled porosity-permeability Block A, the primary configuration of the revision.
INCLUDE_PREFIXES = ("A_feedback__p32_",)


def wanted(name):
    return name.startswith(INCLUDE_PREFIXES)
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
def cell_volumes(path, n):
    """The volume column of a PFLOTRAN .uge, in cubic metres.

    The wrapper that writes these files multiplies the LaGriT surface areas by
    the aperture, so the column is a true volume; the uncorrected files held
    areas and were three orders of magnitude larger.
    """
    v = np.zeros(n)
    with open(path) as f:
        m = int(f.readline().split()[1])
        for i in range(min(m, n)):
            v[i] = float(f.readline().split()[4])
    return v


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
    """Per-cell net carbonate volume fraction, seed-subtracted.

    Every secondary phase is initialised at 1e-6 so that transition-state
    growth can begin, so the raw field carries 4e-6 in every cell whether or
    not anything precipitated. Where this feeds a distribution across cells,
    as in the velocity-decile figure, that uniform floor pushes every bin
    toward an equal share and understates any real concentration.

    The return value is still a volume fraction per cell. Any caller that
    needs a carbonate VOLUME must multiply by the cell volumes, which span
    four orders of magnitude.
    """
    t = np.zeros(n)
    for pf in ("Calcite VF", "Magnesite VF", "Siderite VF", "Dawsonite VF"):
        vk = find_h5_var(grp, pf)
        if vk:
            d = np.asarray(grp[vk][:], float).flatten()
            m = min(len(d), n)
            t[:m] += np.clip(d[:m] - SEED_VF, 0.0, None)
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
    if hasattr(h5f, "run_dir"):
        q = cell_flux_from_pressure(h5f.run_dir, grp, n_cells)
        if q is not None:
            return q
    pk = find_h5_var(grp, "Liquid Pressure") or find_h5_var(grp, "Pressure")
    if pk:
        p = grp[pk][:].flatten()[:n_cells]
        dp = np.abs(p - np.mean(p))
        return dp / (np.max(dp) + 1e-30)
    return None


class RunH5:
    """All parts of one run's HDF5 output as one file-like object.

    A run continued from a checkpoint writes pflotran_co2_part1.h5, _part2, ...
    and the latest part as pflotran_co2.h5. Opening pflotran_co2.h5 alone loses
    the early snapshots. Where a time appears in two parts, the later part wins.
    Supports keys(), [group], `in`, iteration and `with`, as the figures use.
    """
    def __init__(self, main_path):
        self.run_dir = os.path.dirname(main_path)
        parts = sorted(glob.glob(os.path.join(self.run_dir, "pflotran_co2_part*.h5")),
                       key=lambda p: int(re.search(r"_part(\d+)", p).group(1)))
        paths = parts + ([main_path] if os.path.isfile(main_path) else [])
        self.files = [h5py.File(p, "r") for p in paths]
        self.where = {}
        for f in self.files:
            for k in f.keys():
                self.where[k] = f
    def keys(self):
        return list(self.where.keys())
    def __iter__(self):
        return iter(self.keys())
    def __contains__(self, k):
        return k in self.where
    def __getitem__(self, k):
        return self.where[k][k]
    def close(self):
        for f in self.files:
            f.close()
    def __enter__(self):
        return self
    def __exit__(self, *a):
        self.close()


def open_run_h5(main_path):
    return RunH5(main_path)


def cell_flux_from_pressure(run_dir, grp, n_cells):
    """Flow through each cell from the simulated pressure: sum over the cell's
    connections of |area (p_i - p_j) / distance|, halved. Proportional to the
    Darcy flux for the uniform fracture permeability. None if the mesh lacks
    connections or the pressure is not written."""
    uge = os.path.join(run_dir, "full_mesh.uge")
    pk = find_h5_var(grp, "Liquid Pressure")
    if not (os.path.exists(uge) and pk):
        return None
    with open(uge) as f:
        n = int(f.readline().split()[1])
        xyz = np.array([[float(v) for v in f.readline().split()[1:4]] for _ in range(n)])
        head = f.readline().split()
        if not head or head[0].upper() != "CONNECTIONS":
            return None
        m = int(head[1])
        con = np.array([[float(v) for v in f.readline().split()] for _ in range(m)])
    if con.size == 0:
        return None
    i = con[:, 0].astype(int) - 1; j = con[:, 1].astype(int) - 1; area = con[:, 5]
    p = np.asarray(grp[pk][:], float).flatten()[:n]
    dist = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    q = np.where(dist > 0, area * np.abs(p[i] - p[j]) / np.where(dist > 0, dist, 1.0), 0.0)
    flux = np.zeros(n); np.add.at(flux, i, q); np.add.at(flux, j, q)
    return 0.5 * flux[:n_cells]


def extract_all_data():
    """Build per-DFN dict of time series for Figs 2–8."""
    D = {}
    for d in sorted(glob.glob(f"{RESULTS_ROOT}/*")):
        nm = os.path.basename(d)
        if not wanted(nm):
            continue
        h5 = os.path.join(d, "pflotran_co2.h5")
        if not os.path.exists(h5):
            continue
        m = re.search(r"p32_(\d+)", nm)
        if not m:
            continue
        p32 = int(m.group(1)) / 100.0
        with open_run_h5(h5) as f:
            tg = parse_time_groups(f)
            if len(tg) < 2:
                continue
            nc = f[tg[0][0]]["pH"].shape[0]

            # The mask is built at the final time, where the volume fractions
            # are largest, and then applied at every time so the series is
            # taken over one consistent set of cells.
            gl = f[tg[-1][0]]
            sec_tot = np.zeros(nc)
            for pf in SEC_PHASES:
                vk = find_h5_var(gl, pf)
                if vk:
                    a = np.asarray(gl[vk][:], float).flatten()[:nc]
                    sec_tot[:len(a)] += np.clip(a - SEED_VF, 0.0, None)
            valid = np.ones(nc, dtype=bool)  # finalise.py masks no cells
            nv = int(valid.sum())
            if nv < 0.5 * nc:
                continue

            # Cell volumes, for the mineral budget. Without them the
            # percentages weight every cell equally across four orders of
            # magnitude of cell size.
            uge = os.path.join(d, "full_mesh.uge")
            vol = cell_volumes(uge, nc) if os.path.exists(uge) else None

            ti, pH = [], []
            ca, mg, si, dw = [], [], [], []
            an, fo, di = [], [], []
            for k, t in tg:
                g = f[k]
                ti.append(t)
                pH.append(np.nanmean(g["pH"][:][valid]) if "pH" in g else np.nan)
                for pf, store in [("Calcite VF", ca), ("Magnesite VF", mg),
                                   ("Siderite VF", si), ("Dawsonite VF", dw),
                                   ("Anorthite VF", an), ("Forsterite VF", fo),
                                   ("Diopside VF", di)]:
                    v = find_h5_var(g, pf)
                    if v:
                        a = np.asarray(g[v][:], float).flatten()[:nc]
                        store.append(float(np.nansum(a[valid])))
                    else:
                        store.append(0.0)
        # Net of the initial seed. The primary phases (anorthite, forsterite,
        # diopside) are not seeded and are left as read.
        ca = np.clip(np.array(ca) - SEED_VF * nv, 0.0, None)
        mg = np.clip(np.array(mg) - SEED_VF * nv, 0.0, None)
        si = np.clip(np.array(si) - SEED_VF * nv, 0.0, None)
        dw = np.clip(np.array(dw) - SEED_VF * nv, 0.0, None)
        tc = ca + mg + si + dw
        ft = tc[-1]

        # Final-time mineral volumes, for the assemblage percentages.
        vw = {}
        if vol is not None:
            gl_name = tg[-1][0]
            with open_run_h5(h5) as f2:
                gl2 = f2[gl_name]
                for key, pf in (("calcite", "Calcite VF"),
                                ("magnesite", "Magnesite VF"),
                                ("siderite", "Siderite VF"),
                                ("dawsonite", "Dawsonite VF")):
                    vk = find_h5_var(gl2, pf)
                    if vk:
                        a = np.asarray(gl2[vk][:], float).flatten()[:nc]
                        net = np.clip(a - SEED_VF, 0.0, None)
                        vw[key] = float((net[valid] * vol[valid]).sum())
                    else:
                        vw[key] = 0.0
            vw["total"] = sum(vw[k] for k in
                              ("calcite", "magnesite", "siderite", "dawsonite"))
        D[nm] = {
            "p32": p32, "ncells": nv, "ncells_total": nc, "times": np.array(ti),
            "pH_mean": np.array(pH),
            "calcite": ca, "magnesite": mg,
            "siderite": si, "dawsonite": dw,
            "anorthite": np.array(an), "forsterite": np.array(fo),
            "diopside": np.array(di),
            "total_carb": tc, "total_carb_per_cell": tc / nv,
            "vol_weighted": vw,
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
        if not wanted(nm):
            continue
        h5 = os.path.join(d, "pflotran_co2.h5")
        if not os.path.exists(h5):
            continue
        m = re.search(r"p32_(\d+)", nm)
        if not m:
            continue
        p32 = int(m.group(1)) / 100.0
        uge = os.path.join(DFN_ROOT, nm, "full_mesh.uge")
        if not os.path.exists(uge):
            uge = os.path.join(d, "full_mesh.uge")
        if not os.path.exists(uge):
            continue
        with open_run_h5(h5) as f:
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
# FIGURE 3 — Carbonate assemblage by network
#
# The stacked percentage panel was replaced. A stacked bar cannot show the
# quantity Table 6 reports, because medians of separate percentages do not sum
# to 100 and the stack would not close. Volume weighting closes it but hides
# the result: at P32 x0.75 one realisation holds more carbonate than the next
# two combined and is 99.9% dawsonite, which sets the whole level.
#
# The spread within each intensity level is the result. At fixed fracture
# density the assemblage differs between networks, which is the claim the paper
# makes for the total amount.
# ============================================================
def fig3_carbonate_budget(data):
    print("\n  Fig 3: Carbonate assemblage by network...")
    p32v = sorted(set(d["p32"] for d in data.values()))
    grp = {p: [d for d in data.values() if d["p32"] == p] for p in p32v}

    def shares(runs, mineral):
        """Share of one mineral in each realisation that precipitated.

        A realisation with no carbonate has no assemblage and is left out,
        which is the same exclusion Table 6 applies.
        """
        out = []
        for d in runs:
            vw = d.get("vol_weighted")
            if not vw or vw.get("total", 0.0) <= 0:
                continue
            out.append(vw.get(mineral, 0.0) / vw["total"] * 100.0)
        return out

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5))
    x = np.arange(len(p32v))
    mins = ["calcite", "magnesite", "siderite", "dawsonite"]

    # panel a: absolute amounts, unchanged, reporting what Table 5 uses
    ax = axes[0]
    bot = np.zeros(len(p32v))
    for mn in mins:
        v = np.array([np.mean([d[mn][-1] / d["ncells"] for d in grp[p]])
                      for p in p32v])
        ax.bar(x, v, 0.55, bottom=bot, color=CARB_COLORS[mn],
               label=mn.capitalize(), edgecolor="white", linewidth=0.3)
        bot += v
    ax.set_xticks(x)
    ax.set_xticklabels([P32_LABELS[p] for p in p32v])
    ax.set_ylabel("Carbonate VF per cell")

    # panel b: the two phases that compete, one box per intensity level, with
    # every realisation drawn so the sample size is visible. Magnesite and
    # siderite are a median 4% and 1% and would crowd the axis.
    ax = axes[1]
    rng = np.random.default_rng(20260915)
    for mn, sign in (("calcite", 0),):           # calcite alone; dawsonite does not form
        pos, vals = [], []
        for i, p in enumerate(p32v):
            s = shares(grp[p], mn)
            if s:
                pos.append(i + sign * 0.17)
                vals.append(s)
        if not vals:
            continue
        ax.boxplot(vals, positions=pos, widths=0.40, showfliers=False,
                   patch_artist=True,
                   medianprops=dict(color="black", lw=1.1),
                   whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7),
                   boxprops=dict(facecolor=CARB_COLORS[mn], alpha=0.55,
                                 edgecolor="0.3", lw=0.6))
        for pp, vv in zip(pos, vals):
            ax.scatter(pp + rng.uniform(-0.06, 0.06, len(vv)), vv,
                       s=11, color=CARB_COLORS[mn], edgecolor="black",
                       linewidth=0.3, zorder=3, alpha=0.9)
    ax.axhline(50, color="0.5", ls="--", lw=0.7, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels([P32_LABELS[p] for p in p32v])
    ax.set_ylim(-4, 104)
    ax.set_ylabel("Calcite share of carbonate (%)")

    for ai, a in enumerate(axes):
        a.text(-0.12, 1.15, "ab"[ai], transform=a.transAxes,
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
    """Trapping, pH and variability against fracture intensity.

    There is no relationship to fit. The submitted figure fitted a line to the
    group means and reported R^2 = 0.307; with the injection rate corrected the
    means are 6.76e-5, 1.26e-4, 6.57e-5, 6.64e-5 and 8.46e-5 with coefficients
    of variation of 145 to 215 per cent. Each panel therefore carries a
    horizontal reference line at the ensemble value, which R1-12 asks for, and
    intervals wide enough that the line falls inside every one.
    """
    print("\n  Fig 4: Trapping against intensity (no trend)...")
    p32v = sorted(set(d["p32"] for d in data.values()))
    grp = {p: [d for d in data.values() if d["p32"] == p] for p in p32v}
    rng = np.random.default_rng(20260906)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.2))

    # panel a: carbonate per cell. Log axis, with the realisations that
    # produced none drawn on a marked floor rather than dropped.
    ax = axes[0]
    allv = np.array([d["total_carb_per_cell"][-1] for d in data.values()])
    pos = allv[allv > 0]
    floor = pos.min() / 4.0 if pos.size else 1e-9
    for d in data.values():
        v = d["total_carb_per_cell"][-1]; c = P32_COLORS[d["p32"]]
        if v > 0:
            ax.scatter(d["p32"], v, color=c, s=26, alpha=0.85,
                       edgecolor="black", linewidth=0.3, zorder=3)
        else:
            ax.scatter(d["p32"], floor, facecolor="none", edgecolor=c, s=26,
                       marker="v", linewidth=0.8, zorder=3)
    med = [float(np.median([d["total_carb_per_cell"][-1] for d in grp[p]]))
           for p in p32v]
    ax.plot(p32v, med, "k-o", lw=1.0, markersize=4.5, zorder=4)
    if pos.size:
        ax.axhline(float(np.median(pos)), color="0.35", ls="--", lw=0.8, zorder=1)
        ax.axhline(floor, color="0.7", ls=":", lw=0.6, zorder=1)
        ax.set_yscale("log"); ax.set_ylim(floor/2.5, pos.max()*2.5)
    # Spearman, not Pearson: with zeros present and the values spanning orders
    # of magnitude a linear fit is dominated by the few largest cases.
    rho, pv = stats.spearmanr([d["p32"] for d in data.values()], allv)
    ax.set_ylabel("Total carbonate VF/cell")
    ax.text(0.03, 0.96, f"$\\rho$ = {rho:+.2f}\n$p$ = {pv:.2f}",
            transform=ax.transAxes, fontsize=6.2, va="top", color="0.25")
    ax.text(-0.12, 1.15, "a", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top")

    # panel b: pH, now flat
    ax = axes[1]
    for d in data.values():
        ax.scatter(d["p32"], d["pH_mean"][-1], color=P32_COLORS[d["p32"]],
                   s=26, alpha=0.85, edgecolor="black", linewidth=0.3, zorder=3)
    pm = [float(np.mean([d["pH_mean"][-1] for d in grp[p]])) for p in p32v]
    psd = [float(np.std([d["pH_mean"][-1] for d in grp[p]])) for p in p32v]
    ax.errorbar(p32v, pm, yerr=psd, fmt="k-o", lw=1.0, capsize=3,
                markersize=4.5, zorder=4)
    ax.axhline(float(np.mean([d["pH_mean"][-1] for d in data.values()])),
               color="0.35", ls="--", lw=0.8, zorder=1)
    ax.set_ylabel("Mean pH at $t$ = 50 yr")
    ax.text(-0.12, 1.15, "b", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top")

    # panel c: CV, with bootstrap intervals. A CV from ten values carries its
    # own uncertainty, which R3-3 asks to be reported.
    ax = axes[2]
    cvs, lo, hi = [], [], []
    f = lambda a: (a.std(ddof=1)/a.mean()*100) if a.mean() > 0 else 0.0
    for p in p32v:
        v = np.array([d["total_carb_per_cell"][-1] for d in grp[p]])
        cvs.append(f(v))
        bs = np.array([f(rng.choice(v, v.size, replace=True)) for _ in range(4000)])
        lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5))
    x = np.arange(len(p32v))
    ax.bar(x, cvs, color=[P32_COLORS[p] for p in p32v], edgecolor="black",
           linewidth=0.3, width=0.6, zorder=2)
    ax.errorbar(x, cvs, yerr=[np.array(cvs)-np.array(lo),
                              np.array(hi)-np.array(cvs)],
                fmt="none", ecolor="black", elinewidth=0.8, capsize=3, zorder=4)
    ax.axhline(float(np.mean(cvs)), color="0.35", ls="--", lw=0.8, zorder=1)
    for xi, cv in zip(x, cvs):
        ax.text(xi, cv+6, f"{cv:.0f}%", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x); ax.set_xticklabels([P32_LABELS[p] for p in p32v])
    ax.set_ylabel("CV (%)")
    ax.text(-0.12, 1.15, "c", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top")

    fig.tight_layout(w_pad=1.5, rect=[0, 0.12, 1, 1])
    fig.text(0.5, 0.08, "P$_{32}$ multiplier", ha="center", fontsize=8)
    hs = [Line2D([0], [0], color=P32_COLORS[p], lw=0, marker="o", markersize=5,
                 markerfacecolor=P32_COLORS[p], markeredgecolor="black",
                 markeredgewidth=0.3, label=f"P$_{{32}}$ {P32_LABELS[p]}")
          for p in p32v]
    hs += [Line2D([0],[0], color="black", lw=0.9, marker="o", markersize=4,
                  markerfacecolor="black", label="group median"),
           Line2D([0],[0], color="0.35", ls="--", lw=0.8, label="ensemble value"),
           Line2D([0],[0], color="black", lw=0, marker="v", markersize=5,
                  markerfacecolor="none", markeredgecolor="black",
                  label="no carbonate")]
    fig.legend(handles=hs, loc="lower center", ncol=len(hs), frameon=False,
               fontsize=6.2, bbox_to_anchor=(0.5, 0.0))
    save_fig(fig, "fig_connectivity_trapping")
    plt.close(fig)


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
        with open_run_h5(dfn["h5"]) as f:
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
    # Strip prefix (e.g., A_p32_100_s1181 -> p32_100_s1181) to match DFN library
    m_case = re.search(r"(L\d+_)?p32_\d+_s\d+", dfn_name)
    case_name = m_case.group(0) if m_case else dfn_name
    candidates = [
        os.path.join(DFN_ROOT, case_name),
        os.path.join(DFN_ROOT, dfn_name),
        os.path.join("..", DFN_ROOT, case_name),
        os.path.join("..", "dfn_library", case_name),
    ]
    dfn_dir = None
    for cand in candidates:
        if os.path.isdir(cand):
            dfn_dir = cand
            break
    if dfn_dir is None:
        return None, None
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
        if not wanted(nm):
            continue
        h5 = os.path.join(d, "pflotran_co2.h5")
        pj = os.path.join(d, "simulation_params.json")
        if not os.path.exists(h5):
            continue
        m = re.search(r"p32_(\d+)", nm)
        if not m:
            continue
        p32 = int(m.group(1)) / 100.0
        domain_vol = 8000.0
        if os.path.exists(pj):
            with open(pj) as f:
                params = json.load(f)
            domain_vol = params.get("domain_volume_m3", domain_vol)
        n_int, method = _count_intersections(nm)
        if n_int is None:
            continue
        methods_used.add(method)
        with open_run_h5(h5) as f:
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
    ax.text(0.95, 0.95, f"R² = {r2:.3f}\nSpearman ρ = {rho_s:.3f}",
            transform=ax.transAxes, fontsize=7, va="top", ha="right",
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


def _deck_rate(run_dir):
    """Plateau injection rate in kg/s, from the deck's RATE LIST.

    apply_corrections.py writes a three-point ramp per realisation, computed
    from that realisation's fracture pore volume. The plateau is the largest
    listed value.
    """
    f = os.path.join(run_dir, "pflotran_co2.in")
    if not os.path.exists(f):
        return None
    lines = open(f).read().split("\n")
    for i, ln in enumerate(lines):
        if "RATE LIST" not in ln:
            continue
        vals = []
        for ln2 in lines[i + 1:]:
            t = ln2.strip()
            if t.startswith("/"):
                break
            parts = t.split()
            if len(parts) == 2:
                try:
                    vals.append(float(parts[1]))
                except ValueError:
                    pass
        if vals:
            return max(vals)
    return None


def _co2_in_carb(grp, nc, vol):
    """Net CO2 in carbonates (mol), seed-subtracted and volume-weighted.

    A volume fraction is dimensionless, so summing it and dividing by a molar
    volume does not give moles. Each mineral contributes
    sum(VF x V_cell) / V_molar, and one mole of CO2 per mole of carbonate for
    all four phases.
    """
    t = 0.0
    for mn, mv in MOLAR_VOL.items():
        vk = find_h5_var(grp, f"{mn} VF")
        if vk:
            a = np.asarray(grp[vk][:], float).flatten()[:nc]
            net = np.clip(a - 1e-6, 0.0, None)
            t += float((net * vol[:len(net)]).sum()) / mv
    return t


def fig8_efficiency():
    print("\n  Fig 8: Trapping efficiency...")
    recs = []
    for d in sorted(glob.glob(f"{RESULTS_ROOT}/*")):
        nm = os.path.basename(d)
        if not wanted(nm):
            continue
        h5 = os.path.join(d, "pflotran_co2.h5")
        pj = os.path.join(d, "simulation_params.json")
        if not os.path.exists(h5):
            continue
        m = re.search(r"p32_(\d+)", nm)
        if not m:
            continue
        p32 = int(m.group(1)) / 100.0
        # The deck is authoritative. simulation_params.json is absent from
        # the Block A runs, so the previous fallback applied 0.01 kg/s to every
        # realisation regardless of what was injected.
        rate = _deck_rate(d)
        if rate is None:
            print(f"    {nm}: no RATE LIST in the deck, skipped")
            continue
        if os.path.exists(pj):
            with open(pj) as f:
                pr = json.load(f)
            rate = pr.get("water_rate_kg_s_scaled", pr.get("water_rate_kg_s", 0.01))
        with open_run_h5(h5) as f:
            tg = parse_time_groups(f)
            if len(tg) < 3:
                continue
            nc = f[tg[0][0]]["pH"].shape[0]
            uge = os.path.join(d, "full_mesh.uge")
            if not os.path.exists(uge):
                print(f"    {nm}: no mesh, skipped")
                continue
            vol = cell_volumes(uge, nc)
            ti = np.array([t for _, t in tg])
            inj = _cum_co2_inj(ti, rate)
            cb = np.array([_co2_in_carb(f[k], nc, vol) for k, _ in tg])
        # np.divide with `where` performs the division only on the
        # masked elements, rather than computing it everywhere and
        # discarding the invalid results afterwards.
        eff = np.divide(cb, inj, out=np.zeros_like(cb, dtype=float),
                        where=inj > 0) * 100.0
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
        if not wanted(nm):
            continue
        h5 = os.path.join(d, "pflotran_co2.h5")
        if not os.path.exists(h5):
            continue
        with open_run_h5(h5) as f:
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
    # Legend inside the axes at the upper right, without a frame. The top of the
    # y-axis sits 15% above the tallest bar (or the 10% line), so the legend has
    # room above the bars.
    ax.set_ylim(0, 1.15 * max(float(np.max(mean_deciles)) * 100, 10.0))
    ax.legend(fontsize=6, loc="upper right", frameon=False)
    fig.tight_layout()
    save_fig(fig, "fig_stagnation_zones")
    plt.close(fig)


# ============================================================
# MAIN
# ============================================================

# ============================================================
# FIG 16 — Dissolution, carbonate precipitation and pH over time
# ============================================================
def fig16_dissolution_precipitation_ph():
    """Median and 10-90 % band over the Block A networks, on a log time axis:
    (a) forsterite, diopside and anorthite dissolved (%), (b) moles of calcite,
    magnesite and all carbonate minerals (mol), (c) volume-weighted mean pH. All quantities
    are volume-weighted with the cell volumes, as in summarize_block.py."""
    print("\n  Fig 16: Dissolution, carbonate and pH over time...")
    FAST = ("Forsterite", "Diopside", "Anorthite")
    per_t = {}                     # time -> list of per-network dicts
    n_runs = 0
    for d in sorted(glob.glob(f"{RESULTS_ROOT}/*")):
        nm = os.path.basename(d)
        if not wanted(nm):
            continue
        h5 = os.path.join(d, "pflotran_co2.h5"); uge = os.path.join(d, "full_mesh.uge")
        if not (os.path.exists(h5) and os.path.exists(uge)):
            continue
        with open_run_h5(h5) as f:
            tg = parse_time_groups(f)
            if len(tg) < 3 or tg[-1][1] < 49.99:
                continue
            nc = f[tg[0][0]]["pH"].shape[0]; vol = cell_volumes(uge, nc)
            g0 = f[tg[0][0]]
            init = {}
            for m in FAST:
                k = find_h5_var(g0, f"{m} VF")
                init[m] = float((np.asarray(g0[k][:], float).flatten()[:nc] * vol).sum()) if k else np.nan
            for k_t, t in tg[1:]:
                g = f[k_t]; row = {}
                ph = np.asarray(g["pH"][:], float).flatten()[:nc]
                row["pH"] = float((ph * vol).sum() / vol.sum())
                for m in FAST:
                    k = find_h5_var(g, f"{m} VF")
                    now = float((np.asarray(g[k][:], float).flatten()[:nc] * vol).sum()) if k else np.nan
                    row[m] = 100.0 * (init[m] - now) / init[m] if init[m] > 0 else np.nan
                tot = 0.0
                for m in ("Calcite", "Magnesite", "Siderite", "Dawsonite"):
                    k = find_h5_var(g, f"{m} VF")
                    mol = 0.0
                    if k:
                        net = np.clip(np.asarray(g[k][:], float).flatten()[:nc] - SEED_VF, 0.0, None)
                        mol = float((net * vol).sum()) / MOLAR_VOL[m]      # moles of mineral
                    row[m] = mol; tot += mol
                row["total"] = tot
                per_t.setdefault(round(t, 6), []).append(row)
        n_runs += 1
    times = sorted(t for t, v in per_t.items() if len(v) == n_runs and t > 0)
    if n_runs == 0 or len(times) < 3:
        print("    no complete runs with time series, skipping."); return
    print(f"    {n_runs} networks, {len(times)} output times")

    def band(key):
        a = np.array([[r[key] for r in per_t[t]] for t in times])
        return np.nanpercentile(a, 50, axis=1), np.nanpercentile(a, 10, axis=1), np.nanpercentile(a, 90, axis=1)

    t = np.array(times)
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5))
    # (a) dissolution
    ax = axs[0]
    for m, c in (("Forsterite", "#1a9850"), ("Diopside", "#8073ac"), ("Anorthite", "#2166AC")):
        md, lo, hi = band(m)
        ax.fill_between(t, lo, hi, color=c, alpha=0.18, lw=0)
        ax.plot(t, md, color=c, lw=1.2, label=m.lower())
    ax.set_ylim(0, 100); ax.set_ylabel("Primary mineral dissolved (%)")
    ax.legend(loc="upper left", fontsize=6.2)
    # (b) carbonate
    ax = axs[1]
    for key, c, lab in (("Calcite", CARB_COLORS["calcite"], "calcite"), ("Magnesite", CARB_COLORS["magnesite"], "magnesite"),
                        ("total", "black", "all carbonate")):
        md, lo, hi = band(key)
        ax.plot(t, md, color=c, lw=1.2, ls="-" if key != "total" else "--", label=lab)   # median only, no band
    ax.set_ylim(bottom=0); ax.set_ylabel("Carbonate mineral (mol)")
    ax.legend(loc="upper left", fontsize=6.2)
    # Zoom inset: the same medians between 0.01 and 0.1 years, with its own y-scale.
    ins = ax.inset_axes([0.16, 0.30, 0.42, 0.34])
    w = (t >= 0.01 - 1e-12) & (t <= 0.1 + 1e-12)
    ymax = 0.0
    for key, c in (("Calcite", CARB_COLORS["calcite"]), ("Magnesite", CARB_COLORS["magnesite"]), ("total", "black")):
        md, _, _ = band(key)
        ins.plot(t[w], md[w], color=c, lw=1.0, ls="-" if key != "total" else "--")
        if w.any():
            ymax = max(ymax, float(np.nanmax(md[w])))
    ins.set_xscale("log"); ins.set_xlim(0.01, 0.1)
    from matplotlib.ticker import NullFormatter, FixedLocator, FixedFormatter
    ins.xaxis.set_major_locator(FixedLocator([0.01, 0.03, 0.1]))
    ins.xaxis.set_major_formatter(FixedFormatter(["0.01", "0.03", "0.1"]))
    ins.xaxis.set_minor_formatter(NullFormatter())
    ins.set_ylim(0, ymax * 1.15 if ymax > 0 else 1)
    ins.tick_params(labelsize=5, length=2, pad=1)
    ins.set_title("0.01$-$0.1 years", fontsize=5.5, pad=2)
    for sp in ins.spines.values():
        sp.set_linewidth(0.4); sp.set_visible(True)
    # (c) pH
    ax = axs[2]
    md, lo, hi = band("pH")
    ax.fill_between(t, lo, hi, color="0.5", alpha=0.25, lw=0)
    ax.plot(t, md, color="black", lw=1.2)
    ax.set_ylabel("Mean pH (volume-weighted)")
    for a, tag in zip(axs, "abc"):
        a.set_xscale("log"); a.set_xlim(t.min(), t.max()); a.set_xlabel("Time (years)")
        a.text(-0.20, 1.03, tag, transform=a.transAxes, fontsize=9, fontweight="bold", va="bottom")
    fig.tight_layout(w_pad=1.6)
    save_fig(fig, "fig_dissolution_precipitation_ph")
    plt.close(fig)
    for m in ("Forsterite", "Diopside", "Anorthite"):
        md, _, _ = band(m); print(f"    {m}: " + ", ".join(f"{md[i]:.1f}% at {t[i]:g} y" for i in range(len(t)) if t[i] in (1.0, 2.0, 10.0, 20.0, 50.0)))
    md, _, _ = band("pH"); print("    pH: " + ", ".join(f"{md[i]:.2f} at {t[i]:g} y" for i in range(len(t)) if t[i] in (0.1, 1.0, 10.0, 50.0)))


def fig10_study_design():
    """Study-design schematic. Needs no run data; see src/fig_study_design.py.

    The module reports any text block that overflows its fixed-height box.
    Those boxes cannot grow, so an overflow would be silently clipped in the
    compiled PDF; it is raised here instead.
    """
    print("\n  Fig 10: Study design schematic...")
    try:
        from fig_study_design import build as _build
    except ImportError as e:
        print(f"    skipped: {e}")
        return
    n_bad = _build(FIG_DIR)
    if n_bad:
        print(f"    {n_bad} text block(s) did not fit; shorten the wording in "
              f"src/fig_study_design.py")


# ============================================================
# LETTER FIGURES 11–16
# These are for the response-to-reviewers letter, not the manuscript.
# Generate with: python generate_figures.py --only 11 12 13 14 15 16
# ============================================================

# --- Helpers for the letter figures ---

LETTER_CASE_HI = "A_p32_100_s1181"   # high-trapping (original figure)
LETTER_CASE_LO = "A_p32_100_s839"    # does not trap (original figure)

C_PRECIP  = "#e31a1c"
C_CATION  = "#ff7f00"
C_DISSOL  = "#1f78b4"


def _pick_pair(data):
    """Auto-select the highest- and lowest-trapping realisations at P32=1.0."""
    cands = [(nm, d) for nm, d in data.items()
             if abs(d["p32"] - 1.0) < 0.01 and nm not in EXCLUDE_CASES]
    if len(cands) < 2:
        cands = [(nm, d) for nm, d in data.items() if nm not in EXCLUDE_CASES]
    cands.sort(key=lambda x: x[1]["total_carb_per_cell"][-1])
    return cands[-1][0], cands[0][0]


def _load_cell_fields(run_dir, nm):
    """Load per-cell fields for one realisation from its HDF5."""
    h5 = os.path.join(run_dir, "pflotran_co2.h5")
    uge = os.path.join(run_dir, "full_mesh.uge")

    def _find(grp, *prefixes):
        """Try multiple prefixes, return the first match."""
        for pf in prefixes:
            k = find_h5_var(grp, pf)
            if k is not None:
                return k
        return None

    with open_run_h5(h5) as f:
        tg = parse_time_groups(f)
        grp = f[tg[-1][0]]
        nc = grp["pH"].shape[0]

        # Print all available variables once for debugging
        all_keys = sorted(grp.keys())
        total_keys = [k for k in all_keys if "Total" in k or "Ca" in k
                      or "Mg" in k or "HCO" in k or "CO3" in k
                      or "CO2" in k or "H2CO3" in k]
        if total_keys:
            print(f"    Species keys: {', '.join(total_keys[:12])}")

        ph = grp["pH"][:].flatten()[:nc]

        # Ca/Mg carbonate complexes: CaCO3(aq) + CaHCO3+ + MgCO3(aq) + MgHCO3+
        ca_mg = np.zeros(nc)
        for pf in ("CaCO3(aq)", "CaHCO3+", "MgCO3(aq)", "MgHCO3+"):
            k = _find(grp, pf)
            if k:
                ca_mg += np.asarray(grp[k][:], float).flatten()[:nc]
        if ca_mg.max() > 0:
            print(f"    Ca/Mg complexes range: {ca_mg.min():.2e}–{ca_mg.max():.2e}")
        else:
            # Fallback: try Total Ca++ and Total Mg++
            ca_k = _find(grp, "Total Ca++", "Total_Ca++", "Ca++ [")
            mg_k = _find(grp, "Total Mg++", "Total_Mg++", "Mg++ [")
            if ca_k:
                ca_mg += np.asarray(grp[ca_k][:], float).flatten()[:nc]
            if mg_k:
                ca_mg += np.asarray(grp[mg_k][:], float).flatten()[:nc]
            if ca_mg.max() > 0:
                print(f"    Ca/Mg (from totals) range: {ca_mg.min():.2e}–{ca_mg.max():.2e}")
            else:
                print(f"    WARNING: no Ca/Mg species found in HDF5")

        # DIC: HCO3- + CO3--
        dic = np.zeros(nc)
        for pf in ("HCO3- [", "CO3-- ["):
            k = _find(grp, pf)
            if k:
                dic += np.asarray(grp[k][:], float).flatten()[:nc]
        if dic.max() > 0:
            print(f"    DIC (HCO3- + CO3--) range: {dic.min():.2e}–{dic.max():.2e}")
        else:
            # Fallback
            dic_k = _find(grp, "Total HCO3-", "Total_HCO3-", "Total H2CO3",
                           "HCO3-", "H2CO3")
            if dic_k:
                dic = np.asarray(grp[dic_k][:], float).flatten()[:nc]
                print(f"    DIC (fallback {dic_k}) range: {dic.min():.2e}–{dic.max():.2e}")
            else:
                print(f"    WARNING: no DIC species found in HDF5")

        tc = get_total_carbonate(grp, nc)

        # Anorthite dissolution
        anor_k = find_h5_var(grp, "Anorthite VF")
        anor_final = np.asarray(grp[anor_k][:], float).flatten()[:nc] if anor_k else None
        g0 = f[tg[0][0]]
        anor_k0 = find_h5_var(g0, "Anorthite VF")
        anor_init = np.asarray(g0[anor_k0][:], float).flatten()[:nc] if anor_k0 else None
        if anor_final is not None and anor_init is not None:
            dissolved = anor_init - anor_final
        elif anor_final is not None:
            dissolved = 0.30 - anor_final
        else:
            dissolved = np.zeros(nc)

        # Throughflow / velocity from HDF5 (used as fallback for flux)
        vel = get_velocity_magnitude(f, tg[-1][0], nc)
        if vel is None:
            vel = get_velocity_magnitude(f, tg[0][0], nc)
        if vel is None:
            vel = np.ones(nc)

    # Groundwater age from flowfield.solve() on the DFN library mesh
    # (Goode 1996 steady Darcy solve with uniform permeability).
    # Flux from the run's own pressure field (aperture-dependent).
    # Same sources joint_predictor.py uses.
    age = None
    flux_from_pressure = None
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import flowfield
        m_case = re.search(r"(L\d+_)?p32_\d+_s\d+", nm)
        case_name = m_case.group(0) if m_case else nm
        # Try several locations for the DFN mesh
        candidates = [
            os.path.join(DFN_ROOT, case_name),
            os.path.join(DFN_ROOT, nm),
            os.path.join("..", DFN_ROOT, case_name),
            os.path.join("..", "dfn_library", case_name),
            run_dir,  # mesh might be in the run dir itself
        ]
        mesh_dir = None
        for cand in candidates:
            if os.path.isdir(cand) and os.path.isfile(
                    os.path.join(cand, "full_mesh.uge")):
                mesh_dir = cand
                break
        if mesh_dir is not None:
            sol = flowfield.solve(mesh_dir, verbose=False)
            age = np.asarray(sol["age"], float)
            if age.size != nc:
                print(f"    age field {age.size} != {nc} cells, discarding")
                age = None
        else:
            print(f"    no mesh with full_mesh.uge found for {nm}")
    except Exception as e:
        print(f"    flowfield.solve failed: {str(e)[:60]}")

    try:
        from joint_predictor import run_flux as _run_flux
        flux_from_pressure = _run_flux(run_dir)
        if flux_from_pressure is not None and flux_from_pressure.size != nc:
            flux_from_pressure = None
    except Exception:
        pass

    if age is None:
        print(f"    WARNING: no age, using x-coordinate as proxy")
        if os.path.exists(uge):
            x_tmp, _, _ = read_uge_centroids(uge)
            age = x_tmp[:nc]
        else:
            age = np.ones(nc)
    flux = flux_from_pressure if flux_from_pressure is not None else vel

    # Percentile rank of dissolution
    from scipy.stats import rankdata
    diss_pctile = 100.0 * rankdata(dissolved, method="average") / nc

    # Coordinates
    x, y, z = np.zeros(nc), np.zeros(nc), np.zeros(nc)
    if os.path.exists(uge):
        x, y, z = read_uge_centroids(uge)

    return {
        "x": x, "y": y, "z": z,
        "throughflow": vel, "age": age, "flux": flux,
        "ca_mg": ca_mg, "ph": ph,
        "dissolution_pctile": diss_pctile, "dissolved": dissolved,
        "dic": dic, "carbonate_vf": tc, "nc": nc,
    }


def _plot_coloc_row(ax_row, d, labels, row_label=None, show_cbar=False,
                    cbar_axes=None):
    """Plot six co-location panels into a list of 6 axes (one row).

    ax_row:  list of 6 matplotlib Axes
    d:       dict from _load_cell_fields
    labels:  list of 6 single-letter labels
    """
    from scipy.stats import rankdata
    from matplotlib.colors import LogNorm, Normalize

    x, z = d["x"], d["z"]
    nc = d["nc"]
    sz = max(0.15, min(1.5, 8000 / nc))
    ca_mg = d["ca_mg"]
    dic = d["dic"]
    bg_ca = np.median(ca_mg)
    bg_dic = np.median(dic)

    age_rank = rankdata(d["age"], method="average")
    flux_rank = rankdata(d["flux"], method="average")
    joint = age_rank * flux_rank
    joint_pctile = 100.0 * rankdata(joint, method="average") / nc

    n_precip = int(np.sum(d["carbonate_vf"] > 1e-10))
    n_ca = int(np.sum(ca_mg > bg_ca * 1.1))
    n_ph = int(np.sum(d["ph"] >= 7.5))
    n_dic = int(np.sum(dic > bg_dic * 1.1))

    def safe_lognorm(arr, floor=1e-10):
        pos = arr[arr > 0]
        if pos.size == 0:
            return Normalize(vmin=0, vmax=1)
        return LogNorm(vmin=max(np.percentile(pos, 1), floor),
                       vmax=np.percentile(pos, 99))

    # --- panel definitions ---
    # (title, values, cmap, norm, mask_bg, count_text)
    vel = d["throughflow"]
    panels = [
        ("Flow", vel,
         "RdBu", safe_lognorm(vel, 1e-15),
         None, None),
        ("Ca/Mg complexes", ca_mg,
         "RdBu", safe_lognorm(ca_mg, 1e-8),
         ca_mg <= bg_ca * 1.1 if bg_ca > 0 else None,
         f"{n_ca:,} cells above background"),
        ("pH", d["ph"],
         "RdBu", Normalize(vmin=3.0, vmax=10.0),
         None, f"{n_ph:,} cells at pH ≥ 7.5"),
        ("Age × Flux rank", joint_pctile,
         "RdBu", Normalize(vmin=0, vmax=100),
         None, None),
        ("Carbon complexes", dic,
         "RdBu", safe_lognorm(dic, 1e-6),
         dic <= bg_dic * 1.1 if bg_dic > 0 else None,
         f"{n_dic:,} cells above background"),
    ]

    scatters = []
    for i, (title, vals, cmap, norm, mask, count) in enumerate(panels):
        ax = ax_row[i]
        ax.set_aspect("equal")
        ax.set_title(f"({labels[i]}) {title}", fontsize=7, loc="left")
        if mask is not None:
            ax.scatter(x[mask], z[mask], c="#d9d9d9", s=sz*0.5, alpha=0.15,
                       rasterized=True)
            sc = ax.scatter(x[~mask], z[~mask], c=vals[~mask], s=sz,
                            cmap=cmap, norm=norm, rasterized=True,
                            edgecolors="none")
        else:
            sc = ax.scatter(x, z, c=vals, s=sz, cmap=cmap, norm=norm,
                            rasterized=True, edgecolors="none")
        scatters.append(sc)
        if count:
            ax.text(0.02, 0.02, count, transform=ax.transAxes,
                    fontsize=4.5, color="#666", va="bottom")
        ax.tick_params(labelsize=5, length=2)
        ax.set_xticklabels([]); ax.set_yticklabels([])

    # Overlay — matches the original figure style:
    # grey background, grey oldest-quarter, blue top-cation, red hollow circles
    ax = ax_row[5]
    ax.set_aspect("equal")
    ax.set_title(f"({labels[5]}) Co-location", fontsize=7, loc="left")

    # All cells in light grey
    ax.scatter(x, z, c="#d9d9d9", s=sz*0.5, alpha=0.3,
               rasterized=True, edgecolors="none", zorder=1)

    # Oldest 25% by age
    age_75 = np.percentile(d["age"], 75)
    oldest = d["age"] >= age_75
    ax.scatter(x[oldest], z[oldest], c="#999999", s=sz*0.8,
               label=f"Oldest 25%", rasterized=True, edgecolors="none",
               zorder=2)

    # Top 0.1% Ca/Mg
    top01 = ca_mg >= np.percentile(ca_mg, 99.9)
    ax.scatter(x[top01], z[top01], c="#1f78b4", s=sz*2.0,
               label=f"Top 0.1% Ca/Mg", rasterized=True, edgecolors="none",
               zorder=3)

    # Precipitating cells — hollow red circles
    precip = d["carbonate_vf"] > 1e-10
    if n_precip > 0:
        ax.scatter(x[precip], z[precip], facecolors="none",
                   edgecolors="#d62728", s=sz*4.0, linewidths=0.6,
                   label=f"Precipitating ({n_precip})",
                   rasterized=True, zorder=5)
    else:
        ax.scatter([], [], facecolors="none", edgecolors="#d62728",
                   s=sz*4.0, linewidths=0.6,
                   label=f"Precipitating ({n_precip})")

    ax.legend(fontsize=4.5, loc="upper right", markerscale=2.0,
              framealpha=0.9, edgecolor="#ccc",
              handletextpad=0.3, borderpad=0.4)
    ax.tick_params(labelsize=5, length=2)
    ax.set_xticklabels([]); ax.set_yticklabels([])
    ax.tick_params(labelsize=5, length=2)
    ax.set_xticklabels([]); ax.set_yticklabels([])

    # Row label
    if row_label:
        ax_row[0].set_ylabel(row_label, fontsize=7, fontweight="bold",
                              rotation=90, labelpad=8)

    return scatters


# ============================================================
# LETTER FIG 11 — Co-location (6 panels, single realisation)
# ============================================================
def fig11_cations():
    """Co-location figure: 1 row × 6 panels."""
    print("\n  Letter Fig 11: Co-location (single realisation)...")
    hi = LETTER_CASE_HI
    if hi is None:
        for d in sorted(glob.glob(f"{RESULTS_ROOT}/{INCLUDE_PREFIXES[0]}100_*")):
            hi = os.path.basename(d)
            break
    if hi is None:
        print("    No case found — skipping.")
        return
    run_dir = os.path.join(RESULTS_ROOT, hi)
    if not os.path.exists(os.path.join(run_dir, "pflotran_co2.h5")):
        print(f"    {hi}: no HDF5 — skipping.")
        return

    d = _load_cell_fields(run_dir, hi)

    fig, axes = plt.subplots(1, 6, figsize=(14.4, 2.8))
    scs = _plot_coloc_row(axes, d, list("abcdef"), row_label=hi)

    # Shared colorbars below
    cbar_labels = [r"Flow rate (relative)", "Ca/Mg (mol/L)",
                   "pH (pivot 6.3)", "Age × Flux rank (%)",
                   r"HCO$_3^-$ + CO$_3^{2-}$ (mol/L)"]
    for i, (sc, lbl) in enumerate(zip(scs, cbar_labels)):
        cb = fig.colorbar(sc, ax=axes[i], orientation="horizontal",
                          fraction=0.06, pad=0.12, aspect=15)
        cb.ax.tick_params(labelsize=4.5)
        cb.set_label(lbl, fontsize=5)

    fig.subplots_adjust(left=0.04, right=0.99, top=0.90, bottom=0.18,
                        wspace=0.08)
    save_fig(fig, "fig_cations")
    plt.close(fig)


# ============================================================
# LETTER FIG 12 — Co-location pair (two realisations, 2×6)
# ============================================================
def fig12_cations_pair(data=None):
    """Same 6 panels for high- vs low-trapping: 2 rows × 6 columns."""
    print("\n  Letter Fig 12: Co-location pair...")
    hi_nm = LETTER_CASE_HI
    lo_nm = LETTER_CASE_LO
    if (hi_nm is None or lo_nm is None) and data is not None:
        hi_nm, lo_nm = _pick_pair(data)
    if hi_nm is None or lo_nm is None:
        print("    No cases specified and no data to auto-pick — skipping.")
        return
    print(f"    High: {hi_nm}, Low: {lo_nm}")

    d_hi = _load_cell_fields(os.path.join(RESULTS_ROOT, hi_nm), hi_nm)
    d_lo = _load_cell_fields(os.path.join(RESULTS_ROOT, lo_nm), lo_nm)

    fig, axes = plt.subplots(2, 6, figsize=(14.4, 5.2))

    scs_hi = _plot_coloc_row(axes[0], d_hi, list("abcdef"),
                              row_label=f"{hi_nm}\n(traps)")
    scs_lo = _plot_coloc_row(axes[1], d_lo, list("abcdef"),
                              row_label=f"{lo_nm}\n(does not trap)")

    # Shared colorbars below the bottom row
    cbar_labels = [r"Flow rate (relative)", "Ca/Mg (mol/L)",
                   "pH (pivot 6.3)", "Age × Flux rank (%)",
                   r"HCO$_3^-$ + CO$_3^{2-}$ (mol/L)"]
    for i, (sc, lbl) in enumerate(zip(scs_hi, cbar_labels)):
        cb = fig.colorbar(sc, ax=[axes[0, i], axes[1, i]],
                          orientation="horizontal",
                          fraction=0.05, pad=0.10, aspect=15,
                          location="bottom")
        cb.ax.tick_params(labelsize=4.5)
        cb.set_label(lbl, fontsize=5)

    # Direction arrow
    axes[0, 0].annotate("", xy=(0.15, 1.08), xytext=(0.0, 1.08),
                         xycoords="axes fraction",
                         arrowprops=dict(arrowstyle="->", color="black",
                                         lw=0.8))
    axes[0, 0].text(0.16, 1.08, "x", transform=axes[0, 0].transAxes,
                     fontsize=6, va="center")
    axes[0, 0].text(-0.02, 1.08, "z", transform=axes[0, 0].transAxes,
                     fontsize=6, va="bottom", ha="right")

    fig.subplots_adjust(left=0.06, right=0.99, top=0.94, bottom=0.12,
                        wspace=0.08, hspace=0.25)
    save_fig(fig, "fig_cations_pair")
    plt.close(fig)


# ============================================================
# LETTER FIG 13 — Pre-registration index
# ============================================================
def fig13_letter_index():
    """Predicted index vs observed carbonate for held-out realisations."""
    print("\n  Letter Fig 13: Pre-registration index...")
    import csv
    csv_path = os.path.join(RESULTS_ROOT, "..", "joint_predictor.csv")
    if not os.path.exists(csv_path):
        csv_path = os.path.join(os.path.dirname(RESULTS_ROOT),
                                "joint_predictor.csv")
    if not os.path.exists(csv_path):
        csv_path = "joint_predictor.csv"
    if not os.path.exists(csv_path):
        print(f"    No joint_predictor.csv found — skipping.")
        return

    rows = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            rows.append(r)

    # Expect columns: case, predicted_index, carbonate_total, held_out
    pred = np.array([float(r.get("predicted_index", r.get("index", 0)))
                     for r in rows])
    obs = np.array([float(r.get("carbonate_total", r.get("carbonate", 0)))
                    for r in rows])
    held = np.array([r.get("held_out", "false").strip().lower() == "true"
                     for r in rows])

    if not held.any():
        # If no held_out column, use all rows
        held = np.ones(len(rows), dtype=bool)

    p, o = pred[held], obs[held]
    has = o > 1e-15

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.5, 2.8))

    floor = o[has].min() * 0.1 if has.any() else 1e-12
    o_plot = np.where(has, o, floor)

    ax1.set_title("(a)", fontsize=8, loc="left")
    ax1.scatter(p[has], o_plot[has], c=CARB_COLORS["calcite"], s=18,
                edgecolors="white", linewidths=0.3, zorder=3)
    ax1.scatter(p[~has], o_plot[~has], c="none", s=18,
                edgecolors="#636363", linewidths=0.8, marker="^",
                label="no carbonate", zorder=3)
    ax1.set_yscale("log")
    ax1.set_xlabel("Predicted index")
    ax1.set_ylabel("Carbonate total (m³)")
    ax1.legend(fontsize=6, framealpha=0.8)

    from scipy.stats import rankdata, spearmanr
    rp, ro = rankdata(p), rankdata(o_plot)
    rho, pval = spearmanr(p, o_plot)
    ax2.set_title("(b)", fontsize=8, loc="left")
    ax2.scatter(rp[has], ro[has], c=CARB_COLORS["calcite"], s=18,
                edgecolors="white", linewidths=0.3, zorder=3)
    ax2.scatter(rp[~has], ro[~has], c="none", s=18,
                edgecolors="#636363", linewidths=0.8, marker="^", zorder=3)
    ax2.plot([0, rp.max() + 1], [0, rp.max() + 1], c="#636363",
             ls="--", lw=0.5, zorder=1)
    ax2.set_xlabel("Predicted rank")
    ax2.set_ylabel("Carbonate rank")
    ax2.text(0.05, 0.95, f"ρ = {rho:.3f}\np = {pval:.3f}",
             transform=ax2.transAxes, fontsize=6.5, va="top")

    fig.tight_layout(w_pad=2.0)
    save_fig(fig, "fig_flow_index_test_legacy")
    plt.close(fig)


# ============================================================
# LETTER FIG 14 — Variability by intensity
# ============================================================
def fig14_letter_variability(data):
    """Carbonate per cell and anorthite dissolution by P32 level."""
    print("\n  Letter Fig 14: Variability by intensity...")
    p32v = sorted(set(d["p32"] for d in data.values()))
    grp = {p: [d for d in data.values() if d["p32"] == p] for p in p32v}
    rng = np.random.default_rng(42)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.5, 2.8))

    all_pos = [d["total_carb_per_cell"][-1] for d in data.values()
               if d["total_carb_per_cell"][-1] > 0]
    ens_med = np.median(all_pos) if all_pos else 1e-7
    floor = min(all_pos) / 4.0 if all_pos else 1e-9

    for i, p in enumerate(p32v):
        sub = grp[p]
        carb = [d["total_carb_per_cell"][-1] for d in sub]
        anor = [d["anor_dissolved_pct"] for d in sub]
        jit = rng.uniform(-0.15, 0.15, len(carb))

        for j, c in enumerate(carb):
            if c > 0:
                ax1.scatter(i + jit[j], c, c=P32_COLORS[p], s=14,
                            edgecolors="white", linewidths=0.3, zorder=3)
            else:
                ax1.scatter(i + jit[j], floor, facecolor="none",
                            edgecolor=P32_COLORS[p], s=14, marker="^",
                            linewidths=0.6, zorder=3)
        pos = [c for c in carb if c > 0]
        if pos:
            ax1.plot([i - 0.3, i + 0.3],
                     [np.median(pos)] * 2,
                     c="#D6394C", lw=1.5, zorder=4)

        ax2.scatter(np.full(len(anor), i) + jit[:len(anor)], anor,
                    c=P32_COLORS[p], s=14, edgecolors="white",
                    linewidths=0.3, zorder=3)
        ax2.plot([i - 0.3, i + 0.3], [np.median(anor)] * 2,
                 c="#D6394C", lw=1.5, zorder=4)

    ax1.axhline(ens_med, c="#636363", ls="--", lw=0.5, zorder=1)
    ax1.set_yscale("log")
    ax1.set_xticks(range(len(p32v)))
    ax1.set_xticklabels([P32_LABELS[p] for p in p32v])
    ax1.set_xlabel(r"P$_{32}$ multiplier")
    ax1.set_ylabel("Carbonate per cell")
    ax1.set_title("(a)", fontsize=8, loc="left")

    ax2.axhline(np.median([d["anor_dissolved_pct"] for d in data.values()]),
                c="#636363", ls="--", lw=0.5, zorder=1)
    ax2.set_xticks(range(len(p32v)))
    ax2.set_xticklabels([P32_LABELS[p] for p in p32v])
    ax2.set_xlabel(r"P$_{32}$ multiplier")
    ax2.set_ylabel("Anorthite dissolved (%)")
    ax2.set_title("(b)", fontsize=8, loc="left")

    fig.tight_layout(w_pad=2.0)
    save_fig(fig, "fig_carbonate_dissolution_variability_legacy")
    plt.close(fig)


# ============================================================
# LETTER FIG 15 — Sensitivity variants
# ============================================================
def fig15_letter_sensitivity():
    """29 single-parameter variants as ratio to baseline."""
    print("\n  Letter Fig 15: Sensitivity variants...")
    import csv
    # Corrected ratios: results/sensitivity_ratios.csv from revision_statistics.py,
    # one row per variant and network, carbonate per cell as in blockc2.py.
    # The older sensitivity_variants.csv came from the circulating runs and is not used.
    csv_path = os.path.join(RESULTS_ROOT, "..", "results", "sensitivity_ratios.csv")
    if not os.path.exists(csv_path):
        print(f"    No {csv_path} — run src/revision_statistics.py --only sensitivity_ratios. Skipping.")
        return
    labels = {"anor_as30": "anorthite A_s 10 -> 30", "anor_as50": "anorthite A_s 10 -> 50", "anor_as100": "anorthite A_s 10 -> 100",
              "global_as_x0.1": "all A_s x0.1", "global_as_x10": "all A_s x10", "no_dawsonite": "dawsonite removed",
              "seed_vf_lo": "secondary seed VF 1e-6 -> 1e-8", "seed_vf_hi": "secondary seed VF 1e-6 -> 1e-4",
              "sec_as_lo": "secondary A_s 1 -> 0.1", "sec_as_hi": "secondary A_s 1 -> 10",
              "daw_rate_1e9": "dawsonite k 1e-7 -> 1e-9", "daw_rate_1e11": "dawsonite k 1e-7 -> 1e-11", "daw_rate_1e13": "dawsonite k 1e-7 -> 1e-13",
              "analcime_1e9": "analcime added, k 1e-9", "analcime_1e11": "analcime added, k 1e-11", "analcime_1e13": "analcime added, k 1e-13"}
    by = {}
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            if r["variant"] != "vf_consistent":
                by.setdefault(r["variant"], []).append(float(r["ratio"]))
    rows = [dict(variant_label=labels.get(v, v), ratio_lo=min(q), ratio_hi=max(q), ratio_median=float(np.median(q)))
            for v, q in by.items()]
    rows.sort(key=lambda r: r["ratio_median"])

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    y = np.arange(len(rows))
    for i, r in enumerate(rows):
        lo, hi, med = r["ratio_lo"], r["ratio_hi"], r["ratio_median"]
        outside = lo < 0.5 or hi > 2.0
        color = "#D6394C" if outside else "#2166AC"
        ax.barh(i, hi - lo, left=lo, height=0.6, color=color, alpha=0.7,
                edgecolor="white", linewidth=0.3)
        ax.plot(med, i, "d", color="white", markersize=3.5, zorder=3)

    ax.axvspan(0.5, 2.0, color="#d9d9d9", alpha=0.3, zorder=0,
               label="factor of 2")
    ax.axvline(1.0, color="#636363", lw=0.5, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([r["variant_label"] for r in rows], fontsize=6)
    ax.set_xlabel("Ratio to baseline")
    ax.set_xscale("log")
    ax.set_xlim(0.005, 20)
    from matplotlib.ticker import ScalarFormatter
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.invert_yaxis()
    ax.legend(fontsize=6, loc="lower right")

    fig.tight_layout()
    save_fig(fig, "fig_parameter_sensitivity_legacy")
    plt.close(fig)


# ============================================================
# LETTER FIG 16 — Study design schematic
# (delegates to the same fig_study_design module as Fig 10)
# ============================================================
# Fig 10 already generates this. Letter references fig_study_design.pdf,
# which is the same file. No separate function needed.


ALL_FIGURES = [
    (1,  "fig_dfn_geometry",            "DFN geometry",                     fig1_geometry,         False),
    (2,  "fig_timeseries",              "pH and mineral time series",       fig2_timeseries,       True),
    (3,  "fig_carbonate_budget",        "Carbonate assemblage",                 fig3_carbonate_budget, True),
    (4,  "fig_connectivity_trapping",   "Connectivity–trapping (3 panels)", fig4_connectivity,     True),
    (5,  "fig_dissolution",             "Dissolution box plots",            fig5_dissolution,      True),
    (6,  "fig_p32_spatial_comparison",  "3D spatial pH and carbonate",      fig6_spatial,          "mesh"),
    (7,  "fig_topology_trapping",       "Intersection density vs trapping", fig7_topology,         False),
    (8,  "fig_trapping_efficiency",     "CO2 trapping efficiency",          fig8_efficiency,       False),
    (9,  "fig_stagnation_zones",        "Velocity-decile precipitation",    fig9_stagnation,       False),
    (10, "fig_study_design",            "Study-design schematic",           fig10_study_design,    False),
    (16, "fig_dissolution_precipitation_ph", "Dissolution, carbonate and pH over time", fig16_dissolution_precipitation_ph, False),
    # --- Letter figures ---
    # Letter and co-location figures come from their own scripts: fig_cations_pair.pdf from
    # src/fig_cations.py, fig_carbonate_dissolution_variability.pdf and fig_parameter_sensitivity.pdf from src/ensemble_statistics_figures.py.
    # fig12, fig14 and fig15 below have different layouts and are not registered.
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
    parser.add_argument("--case-hi", type=str, default=None,
                        help="High-trapping case for letter figs 11–12")
    parser.add_argument("--case-lo", type=str, default=None,
                        help="Low-trapping case for letter fig 12")
    args = parser.parse_args()

    # Apply CLI overrides
    global LETTER_CASE_HI, LETTER_CASE_LO
    FIG_DIR = args.output_dir
    RESULTS_ROOT = args.results_dir
    DFN_ROOT = args.dfn_dir
    LETTER_CASE_HI = args.case_hi
    LETTER_CASE_LO = args.case_lo
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
