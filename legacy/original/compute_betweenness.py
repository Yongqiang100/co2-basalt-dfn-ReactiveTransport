#!/usr/bin/env python3
"""
compute_betweenness.py — Cell-graph betweenness centrality vs. CO2 trapping

Tests whether betweenness-based topological metrics predict carbonate
mineral trapping per cell, where simpler scalar metrics (P32, intersection
density, backbone fraction, dead-end fraction) failed.

Computes BOTH edge and node betweenness on the cell-connectivity graph
extracted from PFLOTRAN .uge files, then reduces each to several scalar
candidates for inclusion in Table 8 of the manuscript:

  R1: flow-weighted reactive co-location
       sum over cells of (node_betweenness * dissolution_rate_proxy)
  R2: top-k betweenness path trapping
       fraction of total carbonate VF on cells in top 10% betweenness
  R3: mean betweenness of actively precipitating cells
       mean node_betweenness over cells where carb VF > eps
  R4: edge-betweenness Gini concentration
       Gini coefficient of the edge-betweenness distribution
  R5: source-sink edge-betweenness sum
       sum of edge betweenness for source-target paths from injection
       boundary to outflow boundary (edge load proxy)

For each candidate, computes Spearman rho, Pearson r, p-values, and R^2
versus per-cell carbonate volume fraction at t=50 yr. Output is a table
the user picks from for Table 8.

Performance:
  - Cell graphs reach ~250k nodes for the densest DFNs.
  - Full all-pairs betweenness is O(VE), prohibitive at this scale.
  - Default uses k-source approximation: betweenness restricted to paths
    starting at a sample of injection-boundary cells (--k 50 by default).
  - This is more physically meaningful anyway: we care about flow paths
    *from the injection zone*, not arbitrary source-target pairs.

Requirements:
  pip install h5py networkx numpy scipy

Usage:
  python compute_betweenness.py                    # all DFNs, default k=50
  python compute_betweenness.py --k 100            # higher accuracy, slower
  python compute_betweenness.py --dfn p32_100_s42  # single DFN for testing
  python compute_betweenness.py --no-edge          # skip edge betweenness
  python compute_betweenness.py --csv out.csv      # save per-DFN values
"""

import os
import sys
import glob
import time
import argparse
import json
from collections import defaultdict

import numpy as np
import h5py
import networkx as nx
from scipy import stats


RESULTS_ROOT = "pflotran_results"
DFN_ROOT = "dfn_library"


# ============================================================
# I/O helpers — copied from stress_test_analysis.py for self-containment
# ============================================================
def read_uge_centroids(path):
    xs, ys, zs = [], [], []
    with open(path) as f:
        n = int(f.readline().strip().split()[1])
        for _ in range(n):
            p = f.readline().strip().split()
            xs.append(float(p[1]))
            ys.append(float(p[2]))
            zs.append(float(p[3]))
    return np.array(xs), np.array(ys), np.array(zs)


def read_uge_connections(path):
    """Return (connections, areas) where connections is list of (i, j)
    pairs with 1-based cell IDs as in PFLOTRAN .uge files."""
    conns = []
    areas = []
    with open(path) as f:
        line = f.readline().strip().split()
        n_cells = int(line[1])
        for _ in range(n_cells):
            f.readline()
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


def parse_time_groups(h5f):
    out = []
    for k in h5f.keys():
        if "Time" not in k or "failure" in k.lower() or "cut" in k.lower():
            continue
        parts = k.strip().split()
        try:
            ti = parts.index("Time")
            t = float(parts[ti + 1])
            t_yr = t if (len(parts) > ti + 2 and parts[ti + 2] == "y") \
                else t / 3.156e7
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


def get_total_carbonate(grp, n, seed_vf=1e-6):
    """Return per-cell NET carbonate VF (current minus initial seed value).

    Each of the 4 carbonate phases is initialized at seed_vf = 1e-6 per
    cell (set in run_pflotran.py to allow nucleation). The raw VF stored
    in HDF5 is therefore (seed + actual precipitation). We subtract the
    seed and clip at zero, matching the convention used by
    plot_trapping_efficiency.py and the manuscript figures.
    """
    total = np.zeros(n)
    for prefix in ["Calcite VF", "Magnesite VF", "Siderite VF", "Dawsonite VF"]:
        vk = find_h5_var(grp, prefix)
        if vk:
            d = grp[vk][:].flatten()
            net = np.maximum(d[:min(len(d), n)] - seed_vf, 0.0)
            total[:min(len(d), n)] += net
    return total


def get_carbonate_history(h5_path, n_cells, seed_vf=1e-6):
    """Return (final_carb_per_cell, peak_carb_per_cell) over the simulation.

    final = net carbonate VF/cell at t=50 yr (after possible re-dissolution)
    peak  = maximum net carbonate VF/cell ever achieved during simulation

    Distinguishes "never precipitated" from "precipitated then re-dissolved".
    """
    with h5py.File(h5_path, "r") as f:
        tg = parse_time_groups(f)
        if len(tg) < 2:
            return 0.0, 0.0
        # Final
        last_grp = f[tg[-1][0]]
        carb_final = get_total_carbonate(last_grp, n_cells, seed_vf=seed_vf)
        # Peak: scan all time groups for max sum
        peak = 0.0
        for tg_key, _ in tg[1:]:  # skip t=0
            grp = f[tg_key]
            c = get_total_carbonate(grp, n_cells, seed_vf=seed_vf)
            s = float(c.sum() / n_cells)
            if s > peak:
                peak = s
    return float(carb_final.sum() / n_cells), peak


# ============================================================
# Build the cell graph as a networkx Graph
# ============================================================
def build_nx_graph(uge_path):
    conns, _ = read_uge_connections(uge_path)
    G = nx.Graph()
    G.add_edges_from(conns)
    return G


# ============================================================
# Identify injection / outflow cells from centroid x-coordinates
# ============================================================
def injection_cells(xs, frac=0.2):
    """1-based cell IDs in the leftmost frac of the domain."""
    xmin, xmax = xs.min(), xs.max()
    cutoff = xmin + frac * (xmax - xmin)
    return np.where(xs <= cutoff)[0] + 1


def outflow_cells(xs, frac=0.2):
    """1-based cell IDs in the rightmost frac of the domain."""
    xmin, xmax = xs.min(), xs.max()
    cutoff = xmax - frac * (xmax - xmin)
    return np.where(xs >= cutoff)[0] + 1


# ============================================================
# Betweenness — k-source approximation
# ============================================================
def compute_node_betweenness(G, sources, k_max=50, seed=42):
    """Approximate node betweenness using sources as Brandes 'k' sample.
    networkx.betweenness_centrality with k samples from given seed list."""
    rng = np.random.default_rng(seed)
    sources_in_g = [s for s in sources if G.has_node(s)]
    if len(sources_in_g) == 0:
        return None
    if len(sources_in_g) > k_max:
        sources_in_g = list(rng.choice(sources_in_g, size=k_max, replace=False))
    # networkx accepts 'k' (number of pivots) but picks them randomly
    # internally; we want OUR pivots, so we compute by hand using
    # single-source shortest paths and Brandes accumulation.
    return _brandes_node(G, sources_in_g)


def _brandes_node(G, sources):
    """Brandes betweenness restricted to given source set, normalized by
    1 / (n_sources * (|V| - 1)) to be comparable across DFNs."""
    n = G.number_of_nodes()
    bet = {v: 0.0 for v in G.nodes()}

    for s in sources:
        # Single-source shortest paths (BFS, unweighted)
        S = []
        P = {v: [] for v in G.nodes()}
        sigma = dict.fromkeys(G.nodes(), 0.0)
        sigma[s] = 1.0
        d = dict.fromkeys(G.nodes(), -1)
        d[s] = 0
        Q = [s]
        while Q:
            v = Q.pop(0)
            S.append(v)
            for w in G.neighbors(v):
                if d[w] < 0:
                    Q.append(w)
                    d[w] = d[v] + 1
                if d[w] == d[v] + 1:
                    sigma[w] += sigma[v]
                    P[w].append(v)
        delta = dict.fromkeys(G.nodes(), 0.0)
        while S:
            w = S.pop()
            for v in P[w]:
                delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
            if w != s:
                bet[w] += delta[w]

    # Normalize: divide by number of source nodes used (so values are
    # average per-source contribution, comparable across DFNs)
    norm = 1.0 / max(len(sources), 1)
    for v in bet:
        bet[v] *= norm
    return bet


def compute_edge_betweenness(G, sources):
    """Edge betweenness restricted to given source set."""
    bet = {e: 0.0 for e in G.edges()}
    # Make undirected canonical edge keys
    def ek(u, v):
        return (u, v) if u < v else (v, u)
    bet = defaultdict(float)
    for s in sources:
        if not G.has_node(s):
            continue
        S = []
        P = {v: [] for v in G.nodes()}
        sigma = dict.fromkeys(G.nodes(), 0.0)
        sigma[s] = 1.0
        d = dict.fromkeys(G.nodes(), -1)
        d[s] = 0
        Q = [s]
        while Q:
            v = Q.pop(0)
            S.append(v)
            for w in G.neighbors(v):
                if d[w] < 0:
                    Q.append(w)
                    d[w] = d[v] + 1
                if d[w] == d[v] + 1:
                    sigma[w] += sigma[v]
                    P[w].append(v)
        delta = dict.fromkeys(G.nodes(), 0.0)
        while S:
            w = S.pop()
            coeff = (1.0 + delta[w]) / sigma[w]
            for v in P[w]:
                c = sigma[v] * coeff
                bet[ek(v, w)] += c
                delta[v] += c
    norm = 1.0 / max(len(sources), 1)
    return {e: b * norm for e, b in bet.items()}


# ============================================================
# Reduction candidates — turn betweenness into per-DFN scalars
# ============================================================
def gini(arr):
    """Gini coefficient of a non-negative array."""
    a = np.asarray(arr, dtype=float)
    a = a[a >= 0]
    if a.size == 0 or np.all(a == 0):
        return 0.0
    a = np.sort(a)
    n = a.size
    idx = np.arange(1, n + 1)
    return (2.0 * np.sum(idx * a) / (n * np.sum(a))) - (n + 1) / n


def reduce_metrics(node_bet, edge_bet, carb_per_cell, dissolution_proxy,
                   xs, inj_ids, out_ids):
    """Return dict of candidate scalar reductions for one DFN.

    Inputs use 1-based cell IDs (PFLOTRAN convention).
    """
    n_cells = len(carb_per_cell)
    bet_arr = np.zeros(n_cells)
    if node_bet is not None:
        for cid, b in node_bet.items():
            if 1 <= cid <= n_cells:
                bet_arr[cid - 1] = b

    out = {}

    # R1: flow-weighted reactive co-location
    # sum_i (node_bet_i * dissolution_proxy_i)
    out["R1_flow_react_colocation"] = float(np.sum(bet_arr * dissolution_proxy))

    # R2: top-10% betweenness fraction of carbonate
    if bet_arr.sum() > 0:
        thresh = np.percentile(bet_arr, 90)
        mask = bet_arr >= thresh
        total_carb = carb_per_cell.sum()
        if total_carb > 0:
            out["R2_top10pct_carb_fraction"] = \
                float(carb_per_cell[mask].sum() / total_carb)
        else:
            out["R2_top10pct_carb_fraction"] = 0.0
    else:
        out["R2_top10pct_carb_fraction"] = 0.0

    # R3: mean betweenness of actively precipitating cells
    eps = 1e-8
    active = carb_per_cell > eps
    if active.sum() > 0 and bet_arr.sum() > 0:
        out["R3_mean_bet_active_cells"] = float(bet_arr[active].mean())
    else:
        out["R3_mean_bet_active_cells"] = 0.0

    # R4: Gini of edge betweenness distribution
    if edge_bet is not None and len(edge_bet) > 0:
        eb = np.array(list(edge_bet.values()))
        out["R4_edge_bet_gini"] = float(gini(eb))
    else:
        out["R4_edge_bet_gini"] = np.nan

    # R5: edge-betweenness load on the injection->outflow corridor
    # sum of edge betweenness for edges where at least one endpoint is
    # within 30% of the injection-to-outflow line
    if edge_bet is not None and len(edge_bet) > 0:
        xmin, xmax = xs.min(), xs.max()
        # Top 30% of cells by node betweenness anchor the corridor
        if node_bet is not None and bet_arr.sum() > 0:
            thresh = np.percentile(bet_arr, 70)
            corridor = set(np.where(bet_arr >= thresh)[0] + 1)
            load = 0.0
            for (u, v), b in edge_bet.items():
                if u in corridor or v in corridor:
                    load += b
            out["R5_corridor_edge_load"] = float(load)
        else:
            out["R5_corridor_edge_load"] = np.nan
    else:
        out["R5_corridor_edge_load"] = np.nan

    return out


# ============================================================
# Dissolution proxy — anorthite VF loss per cell as fraction
# ============================================================
def get_dissolution_proxy(h5_path, n_cells):
    """Return per-cell anorthite dissolution magnitude at t=50 yr.

    Anorthite is the dominant Ca-supplier and its dissolution rate at
    the local pH is what governs whether enough Ca is released for
    calcite supersaturation. We use the absolute VF decrease, capped
    at zero (no precipitation expected for primary minerals).
    """
    with h5py.File(h5_path, "r") as f:
        tg = parse_time_groups(f)
        if len(tg) < 2:
            return np.zeros(n_cells)
        first_grp = f[tg[0][0]]
        last_grp = f[tg[-1][0]]
        ak = find_h5_var(first_grp, "Anorthite VF")
        if ak is None:
            return np.zeros(n_cells)
        a0 = first_grp[ak][:].flatten()[:n_cells]
        a1 = last_grp[find_h5_var(last_grp, "Anorthite VF")][:].flatten()[:n_cells]
        loss = np.maximum(a0 - a1, 0.0)
    return loss


# ============================================================
# Main per-DFN routine
# ============================================================
def process_dfn(name, k_pivots=50, do_edge=True, verbose=True):
    """Compute all betweenness reductions for one DFN.

    Returns dict with the reduction values plus ncells, p32, etc.
    """
    res_dir = os.path.join(RESULTS_ROOT, name)
    h5_path = os.path.join(res_dir, "pflotran_co2.h5")
    if not os.path.exists(h5_path):
        if verbose:
            print(f"  [{name}] no h5, skip")
        return None

    # Find UGE
    uge_path = os.path.join(DFN_ROOT, name, "full_mesh.uge")
    if not os.path.exists(uge_path):
        uge_path = os.path.join(res_dir, "full_mesh.uge")
    if not os.path.exists(uge_path):
        if verbose:
            print(f"  [{name}] no .uge, skip")
        return None

    parts = name.split("_")
    p32 = int(parts[1]) / 100.0

    t0 = time.time()
    if verbose:
        print(f"  [{name}] reading mesh...", end="", flush=True)
    xs, ys, zs = read_uge_centroids(uge_path)
    n_cells = len(xs)

    G = build_nx_graph(uge_path)
    if verbose:
        print(f" {n_cells} cells, {G.number_of_edges()} edges "
              f"({time.time()-t0:.1f}s)")

    # Sources = injection cells
    inj = injection_cells(xs, frac=0.2)
    out_ids = outflow_cells(xs, frac=0.2)

    # Carbonate at final time AND peak over history
    with h5py.File(h5_path, "r") as f:
        tg = parse_time_groups(f)
        last = f[tg[-1][0]]
        carb = get_total_carbonate(last, n_cells)
    diss = get_dissolution_proxy(h5_path, n_cells)

    # Also track peak carbonate (distinguishes never-precipitated vs
    # precipitated-then-redissolved)
    carb_final_total, carb_peak_total = get_carbonate_history(h5_path, n_cells)

    t0 = time.time()
    if verbose:
        print(f"  [{name}] node betweenness (k={k_pivots})...",
              end="", flush=True)
    nb = compute_node_betweenness(G, inj, k_max=k_pivots)
    if verbose:
        print(f" {time.time()-t0:.1f}s")

    eb = None
    if do_edge:
        t0 = time.time()
        if verbose:
            print(f"  [{name}] edge betweenness (k={k_pivots})...",
                  end="", flush=True)
        # Subsample sources for edge betweenness too
        rng = np.random.default_rng(42)
        srcs = [s for s in inj if G.has_node(s)]
        if len(srcs) > k_pivots:
            srcs = list(rng.choice(srcs, k_pivots, replace=False))
        eb = compute_edge_betweenness(G, srcs)
        if verbose:
            print(f" {time.time()-t0:.1f}s")

    metrics = reduce_metrics(nb, eb, carb, diss, xs, inj, out_ids)

    return {
        "name": name,
        "p32": p32,
        "ncells": n_cells,
        "n_edges": G.number_of_edges(),
        "carb_per_cell_final": carb_final_total,
        "carb_per_cell_peak": carb_peak_total,
        "n_inj_cells": int(len(inj)),
        **metrics,
    }


# ============================================================
# Cross-DFN correlation analysis
# ============================================================
def correlate_metrics(records, target_key="carb_per_cell_final"):
    """For each candidate reduction, compute Spearman/Pearson against
    target. Return dataframe-like list of dicts."""
    metric_keys = ["R1_flow_react_colocation",
                   "R2_top10pct_carb_fraction",
                   "R3_mean_bet_active_cells",
                   "R4_edge_bet_gini",
                   "R5_corridor_edge_load"]
    target = np.array([r[target_key] for r in records], dtype=float)

    out = []
    for k in metric_keys:
        x = np.array([r.get(k, np.nan) for r in records], dtype=float)
        valid = np.isfinite(x) & np.isfinite(target)
        if valid.sum() < 3:
            out.append({"metric": k, "n": int(valid.sum()),
                        "spearman_rho": np.nan, "spearman_p": np.nan,
                        "pearson_r": np.nan, "pearson_p": np.nan,
                        "r_squared": np.nan})
            continue
        rho, prho = stats.spearmanr(x[valid], target[valid])
        r, pr = stats.pearsonr(x[valid], target[valid])
        out.append({"metric": k,
                    "n": int(valid.sum()),
                    "spearman_rho": float(rho),
                    "spearman_p": float(prho),
                    "pearson_r": float(r),
                    "pearson_p": float(pr),
                    "r_squared": float(r * r)})
    return out


# ============================================================
# Main
# ============================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dfn", default=None,
                    help="Single DFN name (e.g. p32_100_s42); default is all")
    ap.add_argument("--k", type=int, default=50,
                    help="Number of source pivots for betweenness "
                         "approximation (default 50)")
    ap.add_argument("--no-edge", action="store_true",
                    help="Skip edge betweenness (faster, R4/R5 will be NaN)")
    ap.add_argument("--csv", default="betweenness_results.csv",
                    help="Output CSV with per-DFN values")
    ap.add_argument("--json", default="betweenness_results.json",
                    help="Output JSON with per-DFN values + correlations")
    args = ap.parse_args()

    if args.dfn:
        names = [args.dfn]
    else:
        names = sorted([
            os.path.basename(d)
            for d in glob.glob(f"{RESULTS_ROOT}/*")
            if os.path.exists(os.path.join(d, "pflotran_co2.h5"))
        ])

    print(f"Processing {len(names)} DFN(s) with k={args.k} "
          f"{'(node only)' if args.no_edge else '(node + edge)'}\n")

    records = []
    for i, name in enumerate(names):
        print(f"[{i+1}/{len(names)}] {name}")
        try:
            rec = process_dfn(name, k_pivots=args.k,
                              do_edge=not args.no_edge, verbose=True)
            if rec is not None:
                records.append(rec)
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
        print()

    if not records:
        print("No records produced.")
        return

    # ----- Save raw values -----
    keys = list(records[0].keys())
    with open(args.csv, "w") as f:
        f.write(",".join(keys) + "\n")
        for r in records:
            f.write(",".join(str(r.get(k, "")) for k in keys) + "\n")
    print(f"Saved per-DFN values: {args.csv}")

    # ----- Sample diagnostics -----
    n_zero_final = sum(1 for r in records if r.get("carb_per_cell_final", 0) <= 0)
    n_zero_peak = sum(1 for r in records if r.get("carb_per_cell_peak", 0) <= 0)
    n_redissolved = sum(1 for r in records
                       if r.get("carb_per_cell_peak", 0) > 0
                       and r.get("carb_per_cell_final", 0) <= 0)
    print("\n" + "=" * 78)
    print("SAMPLE DIAGNOSTICS")
    print("=" * 78)
    print(f"  Total DFNs:                       {len(records)}")
    print(f"  Zero final carbonate:             {n_zero_final}")
    print(f"  Zero peak carbonate (never any):  {n_zero_peak}")
    print(f"  Re-dissolved (peak>0, final=0):   {n_redissolved}")
    print()
    print("  Per-DFN final / peak carbonate per cell:")
    for r in sorted(records, key=lambda x: (x["p32"], x["name"])):
        final = r.get("carb_per_cell_final", 0)
        peak = r.get("carb_per_cell_peak", 0)
        flag = ""
        if final <= 0 and peak > 0:
            flag = "  <- re-dissolved"
        elif final <= 0 and peak <= 0:
            flag = "  <- never any"
        print(f"    {r['name']:<22s}  final={final:.3e}  peak={peak:.3e}{flag}")

    # ----- Correlations against TWO targets -----
    for target in ["carb_per_cell_final", "carb_per_cell_peak"]:
        corr = correlate_metrics(records, target_key=target)

        print("\n" + "=" * 78)
        print(f"CORRELATION OF REDUCTION CANDIDATES WITH: {target}")
        print("=" * 78)
        print(f"{'metric':<32s} {'n':>3s} {'rho':>7s} {'p_rho':>8s} "
              f"{'r':>7s} {'p_r':>8s} {'R^2':>7s}")
        print("-" * 78)
        for c in corr:
            rho_s = f"{c['spearman_rho']:>7.3f}" if not np.isnan(c['spearman_rho']) else "    NaN"
            prho_s = f"{c['spearman_p']:>8.3f}" if not np.isnan(c['spearman_p']) else "     NaN"
            r_s = f"{c['pearson_r']:>7.3f}" if not np.isnan(c['pearson_r']) else "    NaN"
            pr_s = f"{c['pearson_p']:>8.3f}" if not np.isnan(c['pearson_p']) else "     NaN"
            r2_s = f"{c['r_squared']:>7.3f}" if not np.isnan(c['r_squared']) else "    NaN"
            print(f"{c['metric']:<32s} {c['n']:>3d} {rho_s} {prho_s} "
                  f"{r_s} {pr_s} {r2_s}")

        valid_corr = [c for c in corr if not np.isnan(c['r_squared'])]
        if valid_corr:
            best = max(valid_corr, key=lambda c: c['r_squared'])
            print("-" * 78)
            print(f"BEST for {target}: {best['metric']}  "
                  f"R^2 = {best['r_squared']:.3f}  "
                  f"rho = {best['spearman_rho']:.3f}")

    # Save JSON with both correlation sets
    with open(args.json, "w") as f:
        json.dump({"records": records,
                   "correlations_final": correlate_metrics(records, "carb_per_cell_final"),
                   "correlations_peak": correlate_metrics(records, "carb_per_cell_peak"),
                   "args": vars(args)}, f, indent=2, default=str)
    print(f"\nSaved correlation summary: {args.json}")


if __name__ == "__main__":
    main()
