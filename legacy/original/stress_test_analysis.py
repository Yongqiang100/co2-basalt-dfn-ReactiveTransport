#!/usr/bin/env python3
"""
stress_test_analysis.py — Post-processing analyses addressing reviewer stress test points

  Analysis 1: Backbone network extraction and correlation with trapping
              (backbone fraction, backbone tortuosity, backbone alignment)

  Analysis 2: Stagnation zone mapping — velocity vs. precipitation co-location
              (low-velocity cells as precipitation sites)

  Analysis 3: Dead-end fraction — topological metric from DFN graph
              (dead-end fractures / total fractures)

  Analysis 4: Finite-size scaling — CV vs. domain volume
              (does variance decrease with domain size?)

  Figures:
    fig_backbone_trapping.pdf     — backbone metrics vs carbonate VF/cell
    fig_stagnation_zones.pdf      — velocity percentile vs precipitation fraction
    fig_deadend_trapping.pdf      — dead-end fraction vs carbonate VF/cell
    fig_finitesize_scaling.pdf    — CV vs domain volume

Usage:
    python stress_test_analysis.py
    python stress_test_analysis.py --backbone
    python stress_test_analysis.py --stagnation
    python stress_test_analysis.py --deadend
    python stress_test_analysis.py --finitesize
    python stress_test_analysis.py --all

Requirements:
    - HDF5 simulation outputs in pflotran_results/
    - UGE mesh files in dfn_library/ (or pflotran_results/)
    - For backbone/dead-end: dfnWorks connectivity files or LaGriT .inp meshes
"""

import os, sys, glob, argparse, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from collections import defaultdict

try:
    import h5py
except ImportError:
    print("ERROR: h5py required. Install with: pip install h5py")
    sys.exit(1)

# ============================================================
# CONFIGURATION
# ============================================================
RESULTS_ROOT = "pflotran_results"
DFN_ROOT = "dfn_library"
FIG_DIR = "paper_figures"
OUT_DIR = "stress_test_results"
SEC_PER_YEAR = 3.156e7

P32_COLORS = {0.75: "#4477AA", 1.0: "#228833", 1.25: "#CCBB44",
              1.5: "#EE6677", 2.0: "#AA3377"}
P32_LABELS = {0.75: "0.75×", 1.0: "1.00×", 1.25: "1.25×",
              1.5: "1.50×", 2.0: "2.00×"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 8, "axes.titlesize": 8,
    "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.4, "ytick.major.width": 0.4,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 6.5, "legend.frameon": False,
    "lines.linewidth": 1.0,
    "figure.dpi": 150, "savefig.dpi": 600,
    "pdf.fonttype": 42, "mathtext.default": "regular",
})


def save_fig(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"),
                    dpi=600 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    print(f"    Saved: {name}.png/pdf")


def save_csv(data, headers, name):
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, f"{name}.csv")
    with open(path, "w") as f:
        f.write(",".join(headers) + "\n")
        for row in data:
            f.write(",".join(str(x) for x in row) + "\n")
    print(f"    Saved: {path}")


def safe_spearman(x, y):
    """Spearman correlation with constant-input protection."""
    if np.std(x) < 1e-30 or np.std(y) < 1e-30:
        return np.nan, 1.0
    return stats.spearmanr(x, y)


def safe_pearson(x, y):
    """Pearson correlation with constant-input protection."""
    if np.std(x) < 1e-30 or np.std(y) < 1e-30:
        return np.nan, 1.0
    return stats.pearsonr(x, y)


def safe_linregress(x, y):
    """Linear regression with constant-input protection. Returns result or None."""
    if np.std(x) < 1e-30:
        return None
    return stats.linregress(x, y)


# ============================================================
# HDF5 HELPERS
# ============================================================
def parse_time_groups(h5f):
    r = []
    for k in h5f.keys():
        if "Time" not in k or "failure" in k.lower() or "cut" in k.lower():
            continue
        p = k.strip().split()
        try:
            i = p.index("Time")
            t = float(p[i + 1])
            r.append((k, t if (len(p) > i + 2 and p[i + 2] == "y") else t / SEC_PER_YEAR))
        except:
            continue
    r.sort(key=lambda x: x[1])
    return r


def find_h5_var(grp, prefix):
    for k in grp.keys():
        if k.startswith(prefix):
            return k
    return None


def get_total_carbonate(grp, n):
    t = np.zeros(n)
    for pf in ["Calcite VF", "Magnesite VF", "Siderite VF", "Dawsonite VF"]:
        vk = find_h5_var(grp, pf)
        if vk:
            d = grp[vk][:].flatten()
            t[:min(len(d), n)] += d[:min(len(d), n)]
    return t


def read_uge_centroids(path):
    xs, ys, zs = [], [], []
    with open(path) as f:
        n = int(f.readline().strip().split()[1])
        for _ in range(n):
            p = f.readline().strip().split()
            xs.append(float(p[1])); ys.append(float(p[2])); zs.append(float(p[3]))
    return np.array(xs), np.array(ys), np.array(zs)


def read_uge_connections(path):
    """Read cell-cell connections from .uge file.
    Returns list of (cell_i, cell_j) pairs and connection areas."""
    conns = []
    areas = []
    with open(path) as f:
        # First line: CELLS n
        line = f.readline().strip().split()
        n_cells = int(line[1])
        # Skip cell centroid lines
        for _ in range(n_cells):
            f.readline()
        # Connection header: CONNECTIONS n
        line = f.readline().strip().split()
        if len(line) < 2 or line[0].upper() != "CONNECTIONS":
            return conns, areas
        n_conn = int(line[1])
        for _ in range(n_conn):
            parts = f.readline().strip().split()
            if len(parts) >= 2:
                ci, cj = int(parts[0]), int(parts[1])
                area = float(parts[2]) if len(parts) >= 3 else 1.0
                conns.append((ci, cj))
                areas.append(area)
    return conns, areas


def get_velocity_magnitude(h5f, grp_name, n_cells):
    """Extract or estimate velocity magnitude per cell.
    
    PFLOTRAN may store Liquid Velocity [m/yr] or Darcy Velocity.
    If neither is available, estimate from Liquid Flux.
    """
    grp = h5f[grp_name]
    
    # Try direct velocity
    for prefix in ["Liquid Velocity", "Darcy Velocity", "Liquid X-Velocity"]:
        vk = find_h5_var(grp, prefix)
        if vk:
            v = grp[vk][:].flatten()[:n_cells]
            return np.abs(v)
    
    # Try flux-based estimate: look for liquid flux variables
    vx = find_h5_var(grp, "Liquid X-Flux")
    vy = find_h5_var(grp, "Liquid Y-Flux")
    vz = find_h5_var(grp, "Liquid Z-Flux")
    if vx and vy and vz:
        fx = grp[vx][:].flatten()[:n_cells]
        fy = grp[vy][:].flatten()[:n_cells]
        fz = grp[vz][:].flatten()[:n_cells]
        return np.sqrt(fx**2 + fy**2 + fz**2)
    
    # Try material ID based estimate using pressure gradient
    pk = find_h5_var(grp, "Liquid Pressure")
    if pk is None:
        pk = find_h5_var(grp, "Pressure")
    if pk:
        p = grp[pk][:].flatten()[:n_cells]
        # Return pressure as proxy (higher pressure gradient = higher velocity)
        # This is a rough proxy; actual velocity requires mesh connectivity
        dp = np.abs(p - np.mean(p))
        return dp / (np.max(dp) + 1e-30)  # normalized 0-1
    
    return None


def get_material_ids(dfn_dir):
    """Extract material IDs from dfnWorks output or LaGriT mesh.
    Material ID typically corresponds to fracture index.
    Returns dict: cell_id -> material_id."""
    
    # Try materialid.dat (dfnWorks standard output)
    for fname in ["materialid.dat", "material_id.dat", "cellid.dat"]:
        path = os.path.join(dfn_dir, fname)
        if os.path.exists(path):
            try:
                ids = np.loadtxt(path, dtype=int)
                return {i+1: mid for i, mid in enumerate(ids)}
            except:
                pass
    
    # Try full_mesh.inp (LaGriT AVS-UCD format)
    inp_path = os.path.join(dfn_dir, "full_mesh.inp")
    if os.path.exists(inp_path):
        try:
            with open(inp_path) as f:
                header = f.readline().strip().split()
                n_nodes, n_cells = int(header[0]), int(header[1])
                for _ in range(n_nodes):
                    f.readline()
                cell_mats = {}
                for _ in range(n_cells):
                    parts = f.readline().strip().split()
                    cell_id = int(parts[0])
                    mat_id = int(parts[1])
                    cell_mats[cell_id] = mat_id
            return cell_mats
        except:
            pass
    
    return None


# ============================================================
# ANALYSIS 1: BACKBONE NETWORK
# ============================================================
def build_cell_graph(uge_path):
    """Build adjacency graph from UGE mesh connections."""
    conns, areas = read_uge_connections(uge_path)
    graph = defaultdict(set)
    for ci, cj in conns:
        graph[ci].add(cj)
        graph[cj].add(ci)
    return graph, conns, areas


def find_backbone(graph, centroids_x, n_cells, inj_frac=0.2):
    """Identify backbone cells: cells on shortest paths from injection to outflow.
    
    Backbone = union of all cells on shortest paths from injection zone
    (left 20% of domain) to outflow zone (right boundary).
    Uses BFS flood from injection cells.
    """
    x = centroids_x
    x_min, x_max = np.min(x), np.max(x)
    x_range = x_max - x_min
    
    # Injection zone: left 20%
    inj_threshold = x_min + inj_frac * x_range
    # Outflow zone: right 10%
    out_threshold = x_max - 0.1 * x_range
    
    inj_cells = set(i for i in range(1, n_cells + 1) if x[i-1] <= inj_threshold)
    out_cells = set(i for i in range(1, n_cells + 1) if x[i-1] >= out_threshold)
    
    if not inj_cells or not out_cells:
        return set(), 0.0, np.inf
    
    # BFS from injection cells to find distances
    from collections import deque
    dist_from_inj = {}
    queue = deque()
    for c in inj_cells:
        dist_from_inj[c] = 0
        queue.append(c)
    
    while queue:
        curr = queue.popleft()
        for neigh in graph.get(curr, []):
            if neigh not in dist_from_inj:
                dist_from_inj[neigh] = dist_from_inj[curr] + 1
                queue.append(neigh)
    
    # BFS from outflow cells
    dist_from_out = {}
    queue = deque()
    for c in out_cells:
        dist_from_out[c] = 0
        queue.append(c)
    
    while queue:
        curr = queue.popleft()
        for neigh in graph.get(curr, []):
            if neigh not in dist_from_out:
                dist_from_out[neigh] = dist_from_out[curr] + 1
                queue.append(neigh)
    
    # Find shortest path length from inj to out
    min_path = np.inf
    for c in out_cells:
        if c in dist_from_inj:
            min_path = min(min_path, dist_from_inj[c])
    
    if min_path == np.inf:
        # No connected path
        return set(), 0.0, np.inf
    
    # Backbone = cells on any shortest (or near-shortest) path
    # A cell is on a shortest path if dist_from_inj[c] + dist_from_out[c] <= min_path + tolerance
    tolerance = max(1, int(0.05 * min_path))  # 5% tolerance
    backbone = set()
    for c in range(1, n_cells + 1):
        if c in dist_from_inj and c in dist_from_out:
            if dist_from_inj[c] + dist_from_out[c] <= min_path + tolerance:
                backbone.add(c)
    
    # Backbone fraction
    bb_fraction = len(backbone) / n_cells if n_cells > 0 else 0
    
    # Backbone tortuosity: ratio of graph-distance path length to Euclidean distance
    # Euclidean distance from injection centroid to outflow centroid
    if backbone:
        bb_x = centroids_x[np.array(list(backbone)) - 1]
        eucl_dist = x_max - x_min
        graph_dist = min_path  # in cell hops
        # Approximate physical path length: graph_dist * mean_cell_spacing
        mean_spacing = x_range / max(1, len(set(np.round(x, 2))))
        tortuosity = (graph_dist * mean_spacing) / eucl_dist if eucl_dist > 0 else 1.0
    else:
        tortuosity = np.inf
    
    return backbone, bb_fraction, tortuosity


def compute_backbone_alignment(backbone, centroids_x, centroids_y, centroids_z):
    """Compute alignment of backbone with injection-outflow axis (x-axis).
    
    Returns: alignment score (0 = perpendicular, 1 = perfectly aligned with x)
    """
    if len(backbone) < 2:
        return 0.0
    
    bb_idx = np.array(list(backbone)) - 1
    bx = centroids_x[bb_idx]
    by = centroids_y[bb_idx]
    bz = centroids_z[bb_idx]
    
    # PCA on backbone cell positions
    coords = np.column_stack([bx - bx.mean(), by - by.mean(), bz - bz.mean()])
    cov = np.cov(coords.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    
    # Principal direction (largest eigenvalue)
    principal = eigenvectors[:, np.argmax(eigenvalues)]
    
    # Alignment = |cos(angle)| between principal direction and x-axis
    alignment = abs(principal[0])  # dot product with [1, 0, 0]
    
    return alignment


def analysis_backbone():
    """Extract backbone metrics for all 25 DFNs and correlate with trapping."""
    print("\n=== Analysis 1: Backbone Network ===")
    
    records = []
    for d in sorted(glob.glob(f"{RESULTS_ROOT}/*")):
        nm = os.path.basename(d)
        h5 = os.path.join(d, "pflotran_co2.h5")
        if not os.path.exists(h5):
            continue
        
        p32 = int(nm.split("_")[1]) / 100.0
        
        # Find UGE mesh
        uge = os.path.join(DFN_ROOT, nm, "full_mesh.uge")
        if not os.path.exists(uge):
            uge = os.path.join(d, "full_mesh.uge")
        if not os.path.exists(uge):
            print(f"    WARNING: No UGE mesh for {nm}, skipping.")
            continue
        
        # Read mesh
        x, y, z = read_uge_centroids(uge)
        n_cells = len(x)
        
        # Build graph
        graph, conns, areas = build_cell_graph(uge)
        
        if not graph:
            print(f"    WARNING: No connections in UGE for {nm}, skipping.")
            continue
        
        # Find backbone
        backbone, bb_frac, tortuosity = find_backbone(graph, x, n_cells)
        alignment = compute_backbone_alignment(backbone, x, y, z)
        
        # Get carbonate data
        with h5py.File(h5, "r") as f:
            tg = parse_time_groups(f)
            if len(tg) < 2:
                continue
            grp = f[tg[-1][0]]
            nc = grp["pH"].shape[0] if "pH" in grp else 0
            if nc == 0:
                continue
            tc = get_total_carbonate(grp, nc)
            carb_per_cell = np.sum(tc) / nc
            
            # Carbonate in backbone vs non-backbone
            if backbone:
                bb_idx = np.array([c - 1 for c in backbone if c - 1 < nc])
                non_bb_idx = np.array([i for i in range(nc) if (i + 1) not in backbone])
                carb_in_bb = np.sum(tc[bb_idx]) if len(bb_idx) > 0 else 0
                carb_in_nonbb = np.sum(tc[non_bb_idx]) if len(non_bb_idx) > 0 else 0
                bb_carb_fraction = carb_in_bb / (carb_in_bb + carb_in_nonbb) if (carb_in_bb + carb_in_nonbb) > 0 else 0
            else:
                bb_carb_fraction = 0
        
        records.append({
            "name": nm, "p32": p32, "ncells": n_cells,
            "bb_fraction": bb_frac, "bb_tortuosity": tortuosity,
            "bb_alignment": alignment, "bb_carb_fraction": bb_carb_fraction,
            "carb_per_cell": carb_per_cell, "bb_size": len(backbone),
        })
        
        print(f"    {nm}: backbone={len(backbone)}/{n_cells} ({bb_frac:.3f}), "
              f"tortuosity={tortuosity:.2f}, alignment={alignment:.3f}, "
              f"carb/cell={carb_per_cell:.2e}")
    
    if len(records) < 3:
        print("    ERROR: Not enough DFNs with backbone data.")
        return records
    
    # Save CSV
    save_csv(
        [[r["name"], r["p32"], r["ncells"], r["bb_fraction"], r["bb_tortuosity"],
          r["bb_alignment"], r["bb_carb_fraction"], r["carb_per_cell"]]
         for r in records],
        ["DFN", "P32", "Cells", "Backbone_Fraction", "Backbone_Tortuosity",
         "Backbone_Alignment", "Carb_in_Backbone_Frac", "Carb_VF_per_Cell"],
        "backbone_metrics"
    )
    
    # Correlations
    bb_frac = np.array([r["bb_fraction"] for r in records])
    bb_align = np.array([r["bb_alignment"] for r in records])
    carb = np.array([r["carb_per_cell"] for r in records])
    
    for metric_name, metric_vals in [("Backbone fraction", bb_frac),
                                      ("Backbone alignment", bb_align)]:
        rho, p_val = safe_spearman(metric_vals, carb)
        r_p, p_p = safe_pearson(metric_vals, carb)
        r2 = r_p**2 if not np.isnan(r_p) else np.nan
        print(f"\n    {metric_name} vs carb/cell:")
        print(f"      Spearman ρ = {rho:.3f} (p = {p_val:.4f})" if not np.isnan(rho) else f"      Spearman ρ = N/A (constant input)")
        print(f"      Pearson R² = {r2:.3f} (p = {p_p:.4f})" if not np.isnan(r2) else f"      Pearson R² = N/A (constant input)")
    
    # Plot: 3-panel figure
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0))
    
    metrics = [
        ("bb_fraction", "Backbone fraction", "Backbone fraction (cells on shortest path / total)"),
        ("bb_alignment", "Backbone alignment", "Backbone alignment with flow axis"),
        ("bb_carb_fraction", "Carb. in backbone (%)", "Fraction of carbonate in backbone cells"),
    ]
    
    for ax, (key, ylabel, xlabel), lb in zip(axes, metrics, "abc"):
        vals = np.array([r[key] for r in records])
        for r in records:
            ax.scatter(r[key], r["carb_per_cell"], color=P32_COLORS[r["p32"]],
                       s=30, alpha=0.85, edgecolor="black", linewidth=0.3, zorder=3)
        
        # Correlation (safe)
        rho, p_val = safe_spearman(vals, carb)
        r_p, _ = safe_pearson(vals, carb)
        r2 = r_p**2 if not np.isnan(r_p) else np.nan
        
        # Trend line (safe)
        reg = safe_linregress(vals, carb)
        if reg is not None:
            x_fit = np.linspace(vals.min(), vals.max(), 50)
            ax.plot(x_fit, reg.slope * x_fit + reg.intercept, "k--", lw=0.7, zorder=2)
        
        r2_s = f"{r2:.3f}" if not np.isnan(r2) else "N/A"
        rho_s = f"{rho:.3f}" if not np.isnan(rho) else "N/A"
        ax.text(0.05, 0.95, f"R² = {r2_s}\nρ = {rho_s}",
                transform=ax.transAxes, fontsize=6.5, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="#ccc", alpha=0.9))
        ax.set_xlabel(xlabel, fontsize=7)
        ax.set_ylabel("Total carbonate VF/cell" if lb == "a" else "")
        ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))
        ax.yaxis.get_offset_text().set_fontsize(6)
        ax.text(-0.12, 1.12, lb, transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="top")
    
    # Legend
    p32v = sorted(set(r["p32"] for r in records))
    handles = [Line2D([0], [0], color=P32_COLORS[p], lw=0, marker="o", markersize=5,
                      markerfacecolor=P32_COLORS[p], markeredgecolor="black",
                      markeredgewidth=0.3, label=f"P$_{{32}}$ {P32_LABELS[p]}")
               for p in p32v]
    fig.tight_layout(rect=[0, 0.10, 1, 1], w_pad=1.5)
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=6.5, bbox_to_anchor=(0.5, 0.0))
    
    save_fig(fig, "fig_backbone_trapping")
    plt.close(fig)
    
    return records


# ============================================================
# ANALYSIS 2: STAGNATION ZONES
# ============================================================
def analysis_stagnation():
    """Map velocity field to precipitation locations."""
    print("\n=== Analysis 2: Stagnation Zone Mapping ===")
    
    records = []
    all_vel_bins = []
    all_carb_fracs = []
    
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
            
            grp_name = tg[-1][0]
            grp = f[grp_name]
            nc = grp["pH"].shape[0] if "pH" in grp else 0
            if nc == 0:
                continue
            
            tc = get_total_carbonate(grp, nc)
            vel = get_velocity_magnitude(f, grp_name, nc)
        
        if vel is None:
            print(f"    WARNING: No velocity data for {nm}, skipping.")
            continue
        
        # Bin cells by velocity percentile
        vel_pcts = np.percentile(vel[vel > 0], [10, 25, 50, 75, 90]) if np.any(vel > 0) else [0]*5
        
        # For each velocity quartile, compute fraction of total carbonate
        carb_total = np.sum(tc)
        if carb_total < 1e-15:
            carb_frac_low = 0
            carb_frac_high = 0
        else:
            p25 = np.percentile(vel, 25)
            p75 = np.percentile(vel, 75)
            low_vel_mask = vel <= p25
            high_vel_mask = vel >= p75
            carb_frac_low = np.sum(tc[low_vel_mask]) / carb_total
            carb_frac_high = np.sum(tc[high_vel_mask]) / carb_total
        
        # Detailed: fraction of carbonate in each velocity decile
        decile_fracs = []
        for lo_pct, hi_pct in [(0, 10), (10, 20), (20, 30), (30, 40), (40, 50),
                                (50, 60), (60, 70), (70, 80), (80, 90), (90, 100)]:
            lo = np.percentile(vel, lo_pct)
            hi = np.percentile(vel, hi_pct)
            if lo_pct == 0:
                mask = vel <= hi
            elif hi_pct == 100:
                mask = vel > lo
            else:
                mask = (vel > lo) & (vel <= hi)
            
            frac = np.sum(tc[mask]) / carb_total if carb_total > 1e-15 else 0
            decile_fracs.append(frac)
        
        records.append({
            "name": nm, "p32": p32, "ncells": nc,
            "carb_frac_low25": carb_frac_low,
            "carb_frac_high25": carb_frac_high,
            "decile_fracs": decile_fracs,
            "carb_per_cell": np.sum(tc) / nc,
            "has_carb": carb_total > 1e-10,
        })
        
        print(f"    {nm}: carb in low-vel (Q1)={carb_frac_low:.1%}, "
              f"high-vel (Q4)={carb_frac_high:.1%}")
    
    if not records:
        print("    ERROR: No velocity data available.")
        return records
    
    # Aggregate: mean decile fractions across DFNs with precipitation
    active = [r for r in records if r["has_carb"]]
    if not active:
        print("    No DFNs with measurable carbonate for stagnation analysis.")
        return records
    
    mean_deciles = np.mean([r["decile_fracs"] for r in active], axis=0)
    
    # Save CSV
    save_csv(
        [[r["name"], r["p32"], r["carb_frac_low25"], r["carb_frac_high25"],
          r["carb_per_cell"]] + r["decile_fracs"] for r in records],
        ["DFN", "P32", "Carb_in_LowVel_Q1", "Carb_in_HighVel_Q4", "Carb_VF_per_Cell"] +
        [f"Decile_{i+1}" for i in range(10)],
        "stagnation_zone_analysis"
    )
    
    # Plot: single-panel figure
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.5))
    
    decile_labels = [f"{i*10}–{(i+1)*10}" for i in range(10)]
    colors = "#4477AA"
    ax.bar(range(10), mean_deciles * 100, color=colors, edgecolor="black",
           linewidth=0.3, width=0.7)
    ax.set_xticks(range(10))
    ax.set_xticklabels(decile_labels, rotation=45, ha="right", fontsize=5.5)
    ax.set_xlabel("Flow velocity percentile (low $\\rightarrow$ high)")
    ax.set_ylabel("Carbonate fraction (%)")
    ax.axhline(y=10, color="#999", ls=":", lw=0.5, label="Uniform (10%)")
    ax.legend(fontsize=6, loc="upper left")
    ax.text(0.02, 0.88, f"n = {len(active)} DFNs\nwith precipitation",
            transform=ax.transAxes, fontsize=6, ha="left", va="top",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="#ccc", alpha=0.9))
    
    fig.tight_layout()
    
    save_fig(fig, "fig_stagnation_zones")
    plt.close(fig)
    
    return records


# ============================================================
# ANALYSIS 3: DEAD-END FRACTION
# ============================================================
def compute_dead_end_fraction(graph, n_cells):
    """Compute fraction of cells that are topological dead-ends.
    
    A dead-end cell has only 1 connection (degree = 1).
    Also compute cells with degree <= 2 (near-dead-end).
    """
    degrees = {c: len(neighbors) for c, neighbors in graph.items()}
    
    dead_end_cells = sum(1 for c, d in degrees.items() if d == 1)
    near_dead_end = sum(1 for c, d in degrees.items() if d <= 2)
    
    n_in_graph = len(degrees)
    if n_in_graph == 0:
        return 0.0, 0.0, 0.0
    
    de_frac = dead_end_cells / n_in_graph
    nde_frac = near_dead_end / n_in_graph
    mean_degree = np.mean(list(degrees.values()))
    
    return de_frac, nde_frac, mean_degree


def analysis_deadend():
    """Compute dead-end metrics and correlate with trapping."""
    print("\n=== Analysis 3: Dead-End Fraction ===")
    
    records = []
    for d in sorted(glob.glob(f"{RESULTS_ROOT}/*")):
        nm = os.path.basename(d)
        h5 = os.path.join(d, "pflotran_co2.h5")
        if not os.path.exists(h5):
            continue
        
        p32 = int(nm.split("_")[1]) / 100.0
        
        # Find UGE mesh
        uge = os.path.join(DFN_ROOT, nm, "full_mesh.uge")
        if not os.path.exists(uge):
            uge = os.path.join(d, "full_mesh.uge")
        if not os.path.exists(uge):
            continue
        
        x, y, z = read_uge_centroids(uge)
        n_cells = len(x)
        graph, _, _ = build_cell_graph(uge)
        
        if not graph:
            continue
        
        de_frac, nde_frac, mean_deg = compute_dead_end_fraction(graph, n_cells)
        
        with h5py.File(h5, "r") as f:
            tg = parse_time_groups(f)
            if len(tg) < 2:
                continue
            grp = f[tg[-1][0]]
            nc = grp["pH"].shape[0] if "pH" in grp else 0
            if nc == 0:
                continue
            tc = get_total_carbonate(grp, nc)
            carb_per_cell = np.sum(tc) / nc
        
        records.append({
            "name": nm, "p32": p32, "ncells": n_cells,
            "dead_end_frac": de_frac, "near_dead_end_frac": nde_frac,
            "mean_degree": mean_deg, "carb_per_cell": carb_per_cell,
        })
        
        print(f"    {nm}: dead-end={de_frac:.3f}, near-dead-end={nde_frac:.3f}, "
              f"mean_degree={mean_deg:.1f}, carb/cell={carb_per_cell:.2e}")
    
    if len(records) < 3:
        print("    ERROR: Not enough DFNs.")
        return records
    
    # Save CSV
    save_csv(
        [[r["name"], r["p32"], r["ncells"], r["dead_end_frac"],
          r["near_dead_end_frac"], r["mean_degree"], r["carb_per_cell"]]
         for r in records],
        ["DFN", "P32", "Cells", "Dead_End_Fraction", "Near_Dead_End_Fraction",
         "Mean_Degree", "Carb_VF_per_Cell"],
        "dead_end_metrics"
    )
    
    # Correlations
    de_frac = np.array([r["dead_end_frac"] for r in records])
    mean_deg = np.array([r["mean_degree"] for r in records])
    carb = np.array([r["carb_per_cell"] for r in records])
    
    for mname, mvals in [("Dead-end fraction", de_frac), ("Mean degree", mean_deg)]:
        rho, p_val = safe_spearman(mvals, carb)
        r_p, p_p = safe_pearson(mvals, carb)
        print(f"\n    {mname} vs carb/cell:")
        print(f"      Spearman ρ = {rho:.3f} (p = {p_val:.4f})")
        print(f"      Pearson R² = {r_p**2:.3f} (p = {p_p:.4f})")
    
    # Plot: 2-panel
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.5))
    
    for r in records:
        ax1.scatter(r["dead_end_frac"], r["carb_per_cell"],
                    color=P32_COLORS[r["p32"]], s=30, alpha=0.85,
                    edgecolor="black", linewidth=0.3, zorder=3)
        ax2.scatter(r["mean_degree"], r["carb_per_cell"],
                    color=P32_COLORS[r["p32"]], s=30, alpha=0.85,
                    edgecolor="black", linewidth=0.3, zorder=3)
    
    for ax, vals, xlabel, lb in [(ax1, de_frac, "Dead-end cell fraction", "a"),
                                  (ax2, mean_deg, "Mean cell connectivity (degree)", "b")]:
        rho, p_val = safe_spearman(vals, carb)
        r_p, _ = safe_pearson(vals, carb)
        r2 = r_p**2 if not np.isnan(r_p) else np.nan
        reg = safe_linregress(vals, carb)
        if reg is not None:
            x_fit = np.linspace(vals.min(), vals.max(), 50)
            ax.plot(x_fit, reg.slope * x_fit + reg.intercept, "k--", lw=0.7, zorder=2)
        r2_s = f"{r2:.3f}" if not np.isnan(r2) else "N/A"
        rho_s = f"{rho:.3f}" if not np.isnan(rho) else "N/A"
        ax.text(0.05, 0.95, f"R² = {r2_s}\nρ = {rho_s}",
                transform=ax.transAxes, fontsize=6.5, va="top",
                bbox=dict(boxstyle="round", facecolor="white", edgecolor="#ccc", alpha=0.9))
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Total carbonate VF/cell" if lb == "a" else "")
        ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))
        ax.yaxis.get_offset_text().set_fontsize(6)
        ax.text(-0.12, 1.12, lb, transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="top")
    
    p32v = sorted(set(r["p32"] for r in records))
    handles = [Line2D([0], [0], color=P32_COLORS[p], lw=0, marker="o", markersize=5,
                      markerfacecolor=P32_COLORS[p], markeredgecolor="black",
                      markeredgewidth=0.3, label=f"P$_{{32}}$ {P32_LABELS[p]}")
               for p in p32v]
    fig.tight_layout(rect=[0, 0.10, 1, 1], w_pad=2.0)
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=6.5, bbox_to_anchor=(0.5, 0.0))
    
    save_fig(fig, "fig_deadend_trapping")
    plt.close(fig)
    
    return records


# ============================================================
# ANALYSIS 4: FINITE-SIZE SCALING
# ============================================================
def analysis_finitesize():
    """Examine whether CV depends on domain volume."""
    print("\n=== Analysis 4: Finite-Size Scaling ===")
    
    records = []
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
        
        with h5py.File(h5, "r") as f:
            tg = parse_time_groups(f)
            if len(tg) < 2:
                continue
            grp = f[tg[-1][0]]
            nc = grp["pH"].shape[0] if "pH" in grp else 0
            if nc == 0:
                continue
            tc = get_total_carbonate(grp, nc)
            carb_per_cell = np.sum(tc) / nc
        
        records.append({
            "name": nm, "p32": p32, "ncells": nc,
            "domain_vol": domain_vol, "carb_per_cell": carb_per_cell,
        })
    
    if len(records) < 5:
        print("    ERROR: Not enough DFNs.")
        return
    
    # Group by P32
    p32v = sorted(set(r["p32"] for r in records))
    grp = {p: [r for r in records if r["p32"] == p] for p in p32v}
    
    # Compute CV and mean domain volume per group
    group_stats = []
    for p in p32v:
        vals = [r["carb_per_cell"] for r in grp[p]]
        vols = [r["domain_vol"] for r in grp[p]]
        mu = np.mean(vals)
        sigma = np.std(vals)
        cv = (sigma / mu * 100) if mu > 0 else 0
        mean_vol = np.mean(vols)
        std_vol = np.std(vols)
        group_stats.append({
            "p32": p, "cv": cv, "mean_vol": mean_vol, "std_vol": std_vol,
            "n": len(vals), "mean_carb": mu, "std_carb": sigma,
        })
        print(f"    P32 {P32_LABELS[p]}: CV={cv:.1f}%, "
              f"domain_vol={mean_vol:.0f}±{std_vol:.0f} m³, n={len(vals)}")
    
    # Save CSV
    save_csv(
        [[g["p32"], g["cv"], g["mean_vol"], g["std_vol"], g["n"],
          g["mean_carb"], g["std_carb"]] for g in group_stats],
        ["P32", "CV_pct", "Mean_Domain_Vol_m3", "Std_Domain_Vol_m3",
         "N_realizations", "Mean_Carb_per_Cell", "Std_Carb_per_Cell"],
        "finite_size_scaling"
    )
    
    # Correlation: CV vs mean domain volume
    vols = np.array([g["mean_vol"] for g in group_stats])
    cvs = np.array([g["cv"] for g in group_stats])
    rho, p_val = safe_spearman(vols, cvs)
    print(f"\n    CV vs domain volume: Spearman ρ = {rho:.3f} (p = {p_val:.4f})")
    
    # Plot: 2-panel
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2))
    
    # Panel a: CV vs mean domain volume
    for g in group_stats:
        ax1.scatter(g["mean_vol"], g["cv"], color=P32_COLORS[g["p32"]],
                    s=60, edgecolor="black", linewidth=0.5, zorder=3)
        ax1.errorbar(g["mean_vol"], g["cv"], xerr=g["std_vol"],
                     fmt="none", ecolor="#666", capsize=3, capthick=0.5, zorder=2)
    
    ax1.axhline(y=100, color="#999", ls=":", lw=0.5)
    ax1.text(0.05, 0.95, f"Spearman ρ = {rho:.3f}\np = {p_val:.3f}",
             transform=ax1.transAxes, fontsize=6.5, va="top",
             bbox=dict(boxstyle="round", facecolor="white", edgecolor="#ccc", alpha=0.9))
    ax1.set_xlabel("Mean domain volume (m$^3$)")
    ax1.set_ylabel("Coefficient of variation (%)")
    ax1.text(-0.12, 1.12, "a", transform=ax1.transAxes, fontsize=10,
             fontweight="bold", va="top")
    
    # Panel b: Individual carbonate VF/cell vs domain volume (all 25 points)
    for r in records:
        ax2.scatter(r["domain_vol"], r["carb_per_cell"],
                    color=P32_COLORS[r["p32"]], s=25, alpha=0.75,
                    edgecolor="black", linewidth=0.3, zorder=3)
    
    ax2.set_xlabel("Domain volume (m$^3$)")
    ax2.set_ylabel("Total carbonate VF/cell")
    ax2.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))
    ax2.yaxis.get_offset_text().set_fontsize(6)
    ax2.text(-0.12, 1.12, "b", transform=ax2.transAxes, fontsize=10,
             fontweight="bold", va="top")
    
    p32v_plot = sorted(set(r["p32"] for r in records))
    handles = [Line2D([0], [0], color=P32_COLORS[p], lw=0, marker="o", markersize=5,
                      markerfacecolor=P32_COLORS[p], markeredgecolor="black",
                      markeredgewidth=0.3, label=f"P$_{{32}}$ {P32_LABELS[p]}")
               for p in p32v_plot]
    fig.tight_layout(rect=[0, 0.10, 1, 1], w_pad=2.0)
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=6.5, bbox_to_anchor=(0.5, 0.0))
    
    save_fig(fig, "fig_finitesize_scaling")
    plt.close(fig)


# ============================================================
# SUMMARY TABLE
# ============================================================
def print_metric_comparison(bb_records, de_records):
    """Print comparison of all topological metrics vs P32 and intersection density."""
    print("\n" + "=" * 70)
    print("METRIC COMPARISON SUMMARY")
    print("=" * 70)
    print(f"{'Metric':<35s} {'Spearman ρ':>12s} {'p-value':>10s} {'R²':>8s}")
    print("-" * 70)
    
    # P32 (from all records)
    if bb_records:
        p32s = np.array([r["p32"] for r in bb_records])
        carb = np.array([r["carb_per_cell"] for r in bb_records])
        rho, p = safe_spearman(p32s, carb)
        r_p, _ = safe_pearson(p32s, carb)
        r2 = r_p**2 if not np.isnan(r_p) else np.nan
        rho_s = f"{rho:.3f}" if not np.isnan(rho) else "N/A"
        r2_s = f"{r2:.3f}" if not np.isnan(r2) else "N/A"
        print(f"{'P32 (fracture intensity)':<35s} {rho_s:>12s} {p:>10.4f} {r2_s:>8s}")
    
    # Backbone metrics
    if bb_records:
        for key, label in [("bb_fraction", "Backbone fraction"),
                            ("bb_alignment", "Backbone alignment"),
                            ("bb_carb_fraction", "Carb. in backbone (frac)")]:
            vals = np.array([r[key] for r in bb_records])
            rho, p = safe_spearman(vals, carb)
            r_p, _ = safe_pearson(vals, carb)
            r2 = r_p**2 if not np.isnan(r_p) else np.nan
            rho_s = f"{rho:.3f}" if not np.isnan(rho) else "N/A"
            r2_s = f"{r2:.3f}" if not np.isnan(r2) else "N/A"
            print(f"{label:<35s} {rho_s:>12s} {p:>10.4f} {r2_s:>8s}")
    
    # Dead-end metrics
    if de_records:
        carb_de = np.array([r["carb_per_cell"] for r in de_records])
        for key, label in [("dead_end_frac", "Dead-end cell fraction"),
                            ("mean_degree", "Mean cell degree")]:
            vals = np.array([r[key] for r in de_records])
            rho, p = safe_spearman(vals, carb_de)
            r_p, _ = safe_pearson(vals, carb_de)
            r2 = r_p**2 if not np.isnan(r_p) else np.nan
            rho_s = f"{rho:.3f}" if not np.isnan(rho) else "N/A"
            r2_s = f"{r2:.3f}" if not np.isnan(r2) else "N/A"
            print(f"{label:<35s} {rho_s:>12s} {p:>10.4f} {r2_s:>8s}")
    
    print("=" * 70)
    print("Note: P32 R² = 0.307 and intersection density R² = 0.037 from")
    print("      the main manuscript analysis (generate_all_figures.py)")
    print("=" * 70)


# ============================================================
# MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="Stress test analyses for WRR manuscript")
    parser.add_argument("--backbone", action="store_true",
                        help="Run backbone network analysis")
    parser.add_argument("--stagnation", action="store_true",
                        help="Run stagnation zone mapping")
    parser.add_argument("--deadend", action="store_true",
                        help="Run dead-end fraction analysis")
    parser.add_argument("--finitesize", action="store_true",
                        help="Run finite-size scaling analysis")
    parser.add_argument("--all", action="store_true",
                        help="Run all analyses")
    args = parser.parse_args()
    
    # Default: run all
    if not any([args.backbone, args.stagnation, args.deadend, args.finitesize, args.all]):
        args.all = True
    
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)
    
    print("=" * 60)
    print("Stress Test Post-Processing Analyses")
    print(f"Results: {RESULTS_ROOT}/")
    print(f"DFN library: {DFN_ROOT}/")
    print(f"Figures: {FIG_DIR}/")
    print(f"CSV output: {OUT_DIR}/")
    print("=" * 60)
    
    bb_records = None
    de_records = None
    
    if args.backbone or args.all:
        bb_records = analysis_backbone()
    
    if args.stagnation or args.all:
        analysis_stagnation()
    
    if args.deadend or args.all:
        de_records = analysis_deadend()
    
    if args.finitesize or args.all:
        analysis_finitesize()
    
    # Summary comparison
    if bb_records or de_records:
        print_metric_comparison(bb_records or [], de_records or [])
    
    print(f"\n{'=' * 60}")
    print("Done. Check figures in paper_figures/ and data in stress_test_results/")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()