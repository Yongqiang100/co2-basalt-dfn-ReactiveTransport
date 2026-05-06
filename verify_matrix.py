"""
verify_matrix.py — Full quality verification of the DFN matrix

Checks all DFNs in dfn_library/ and produces:
  1. Summary table with pass/fail for each case
  2. Statistical consistency across realizations
  3. Mesh quality analysis
  4. Boundary file verification
  5. Domain coverage check
  6. Publication-ready summary figure (matrix_verification.png)

Usage:
    python verify_matrix.py
    python verify_matrix.py --verbose       # detailed per-DFN output
    python verify_matrix.py --fix           # attempt to regenerate .ex files
"""

import os
import sys
import json
import glob
import argparse
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

DFN_ROOT = os.path.join(os.getcwd(), "dfn_library")

# Quality thresholds
THRESH = {
    "min_fractures": 10,
    "min_connections": 5,
    "min_conn_per_frac": 1.5,
    "max_ar5_pct": 5.0,         # % elements with aspect ratio > 5
    "max_small_angle_pct": 5.0,  # % elements with min angle < 10 deg
    "max_obtuse_pct": 10.0,      # % elements with max angle > 120 deg
    "domain_coverage_frac": 0.7,  # must cover 70% of domain extent
    "min_ex_files": 6,
    "min_right_ex_lines": 2,     # need outflow boundary
}


# =================================================================
# Load all DFN summaries
# =================================================================
def load_summaries():
    summaries = []
    for entry in sorted(os.listdir(DFN_ROOT)):
        path = os.path.join(DFN_ROOT, entry)
        if not os.path.isdir(path):
            continue
        sp = os.path.join(path, "dfn_summary.json")
        if os.path.exists(sp):
            with open(sp) as f:
                s = json.load(f)
            s["_dir"] = path
            summaries.append(s)
        else:
            summaries.append({
                "name": entry,
                "_dir": path,
                "_error": "No dfn_summary.json"
            })
    return summaries


# =================================================================
# Mesh quality check for one DFN
# =================================================================
def check_mesh_quality(dfn_dir):
    inp_path = os.path.join(dfn_dir, "full_mesh.inp")
    if not os.path.exists(inp_path):
        return {"error": "full_mesh.inp not found"}

    nodes = []
    elements = []
    mat_ids = []

    with open(inp_path) as f:
        header = f.readline().strip().split()
        n_nodes = int(header[0])
        n_elements = int(header[1])

        for _ in range(n_nodes):
            parts = f.readline().strip().split()
            nodes.append([float(parts[1]), float(parts[2]), float(parts[3])])

        for _ in range(n_elements):
            parts = f.readline().strip().split()
            mat_ids.append(int(parts[1]))
            node_ids = [int(p) - 1 for p in parts[3:]]
            elements.append(node_ids)

    nodes = np.array(nodes)

    result = {
        "n_nodes": n_nodes,
        "n_elements": n_elements,
        "x_range": [float(nodes[:, 0].min()), float(nodes[:, 0].max())],
        "y_range": [float(nodes[:, 1].min()), float(nodes[:, 1].max())],
        "z_range": [float(nodes[:, 2].min()), float(nodes[:, 2].max())],
        "n_materials": len(np.unique(mat_ids)),
    }

    # Triangle quality
    tris = [e for e in elements if len(e) == 3]
    if not tris:
        result["error"] = "No triangles"
        return result

    aspect_ratios = []
    min_angles = []
    max_angles = []

    for nids in tris:
        p0, p1, p2 = nodes[nids[0]], nodes[nids[1]], nodes[nids[2]]
        e01 = np.linalg.norm(p1 - p0)
        e12 = np.linalg.norm(p2 - p1)
        e20 = np.linalg.norm(p0 - p2)
        eds = [e01, e12, e20]
        mn, mx = min(eds), max(eds)
        aspect_ratios.append(mx / mn if mn > 0 else 999)

        def ang(a, b, c):
            v = np.clip((a**2 + b**2 - c**2) / (2 * a * b + 1e-30), -1, 1)
            return np.degrees(np.arccos(v))

        if all(x > 0 for x in eds):
            a1 = ang(e01, e12, e20)
            a2 = ang(e12, e20, e01)
            a3 = ang(e20, e01, e12)
            min_angles.append(min(a1, a2, a3))
            max_angles.append(max(a1, a2, a3))

    ar = np.array(aspect_ratios)
    mina = np.array(min_angles)
    maxa = np.array(max_angles)

    result["ar_mean"] = float(np.mean(ar))
    result["ar_max"] = float(np.max(ar))
    result["ar5_pct"] = float(100 * np.mean(ar > 5))
    result["min_angle_min"] = float(np.min(mina)) if len(mina) > 0 else 0
    result["min_angle_mean"] = float(np.mean(mina)) if len(mina) > 0 else 0
    result["small_angle_pct"] = float(100 * np.mean(mina < 10)) if len(mina) > 0 else 0
    result["obtuse_pct"] = float(100 * np.mean(maxa > 120)) if len(maxa) > 0 else 0

    return result


# =================================================================
# Boundary file check
# =================================================================
def check_boundaries(dfn_dir):
    result = {}
    ex_files = glob.glob(os.path.join(dfn_dir, "*.ex"))
    result["n_ex"] = len(ex_files)

    for ef in ex_files:
        name = os.path.basename(ef)
        with open(ef) as f:
            n_lines = len(f.readlines())
        result[name] = n_lines

    right_ex = os.path.join(dfn_dir, "boundary_right_e.ex")
    result["has_outflow"] = (
        os.path.exists(right_ex) and
        result.get("boundary_right_e.ex", 0) > 1
    )
    return result


# =================================================================
# Domain coverage check
# =================================================================
def check_domain_coverage(mesh_result, domain_size=20.0):
    half = domain_size / 2.0
    coverage = {}
    for axis, key in [("x", "x_range"), ("y", "y_range"), ("z", "z_range")]:
        if key in mesh_result:
            rng = mesh_result[key]
            extent = rng[1] - rng[0]
            coverage[axis] = extent / domain_size
        else:
            coverage[axis] = 0.0
    coverage["min"] = min(coverage.values()) if coverage else 0.0
    return coverage


# =================================================================
# Run all checks
# =================================================================
def verify_all(verbose=False):
    summaries = load_summaries()
    if not summaries:
        print("No DFNs found in dfn_library/")
        return []

    results = []

    print(f"\n{'='*100}")
    print(f"DFN MATRIX VERIFICATION — {len(summaries)} cases")
    print(f"{'='*100}\n")

    for s in summaries:
        name = s.get("name", "?")
        dfn_dir = s.get("_dir", "")
        r = {"name": name, "issues": [], "warnings": []}

        # Check for generation error
        if "_error" in s:
            r["issues"].append(s["_error"])
            r["status"] = "FAIL"
            results.append(r)
            continue

        r["num_fractures"] = s.get("num_fractures", 0)
        r["n_elements"] = s.get("n_elements", 0)
        r["n_connections"] = s.get("n_connections", 0)
        r["conn_per_frac"] = s.get("conn_per_frac", 0)
        r["p32_mult"] = s.get("p32_mult", 0)
        r["seed"] = s.get("seed", 0)

        # Check 1: Minimum fractures
        if r["num_fractures"] < THRESH["min_fractures"]:
            r["issues"].append(
                f"Only {r['num_fractures']} fractures "
                f"(need {THRESH['min_fractures']})")

        # Check 2: Connectivity
        if r["conn_per_frac"] < THRESH["min_conn_per_frac"]:
            r["warnings"].append(
                f"Low connectivity: {r['conn_per_frac']:.1f} conn/frac")

        # Check 3: Mesh quality
        mesh = check_mesh_quality(dfn_dir)
        r["mesh"] = mesh
        if "error" in mesh:
            r["issues"].append(f"Mesh: {mesh['error']}")
        else:
            if mesh.get("ar5_pct", 0) > THRESH["max_ar5_pct"]:
                r["issues"].append(
                    f"AR>5: {mesh['ar5_pct']:.1f}% "
                    f"(limit {THRESH['max_ar5_pct']}%)")
            if mesh.get("small_angle_pct", 0) > THRESH["max_small_angle_pct"]:
                r["issues"].append(
                    f"Angle<10°: {mesh['small_angle_pct']:.1f}% "
                    f"(limit {THRESH['max_small_angle_pct']}%)")
            if mesh.get("obtuse_pct", 0) > THRESH["max_obtuse_pct"]:
                r["warnings"].append(
                    f"Obtuse: {mesh['obtuse_pct']:.1f}%")

        # Check 4: Domain coverage
        if "error" not in mesh:
            cov = check_domain_coverage(mesh)
            r["coverage"] = cov
            if cov["min"] < THRESH["domain_coverage_frac"]:
                r["warnings"].append(
                    f"Domain coverage: "
                    f"X={cov['x']:.0%} Y={cov['y']:.0%} Z={cov['z']:.0%}")

        # Check 5: Boundary files
        bc = check_boundaries(dfn_dir)
        r["boundaries"] = bc
        if bc["n_ex"] < THRESH["min_ex_files"]:
            r["issues"].append(
                f"Only {bc['n_ex']} .ex files "
                f"(need {THRESH['min_ex_files']})")
        if not bc.get("has_outflow", False):
            r["issues"].append("No outflow boundary (right_e.ex)")

        # Check 6: UGE file
        uge_path = os.path.join(dfn_dir, "full_mesh.uge")
        if os.path.exists(uge_path):
            with open(uge_path) as f:
                header = f.readline().strip().split()
            if len(header) < 2 or header[0].upper() != "CELLS":
                r["issues"].append("Invalid UGE header")
        else:
            r["issues"].append("full_mesh.uge missing")

        # Status
        if r["issues"]:
            r["status"] = "FAIL"
        elif r["warnings"]:
            r["status"] = "WARN"
        else:
            r["status"] = "PASS"

        results.append(r)

    # ---- Print summary table ----
    print(f"{'Name':<20s} {'P32x':>5s} {'Seed':>5s} {'Frac':>5s} "
          f"{'Elems':>8s} {'C/F':>5s} {'AR5%':>5s} {'<10°':>5s} "
          f"{'BCs':>3s} {'Out':>3s} {'Status':>6s}")
    print("-" * 100)

    n_pass, n_warn, n_fail = 0, 0, 0
    for r in results:
        status = r["status"]
        if status == "PASS":
            n_pass += 1
            sym = "OK"
        elif status == "WARN":
            n_warn += 1
            sym = "WARN"
        else:
            n_fail += 1
            sym = "FAIL"

        mesh = r.get("mesh", {})
        bc = r.get("boundaries", {})

        print(f"{r['name']:<20s} "
              f"{r.get('p32_mult', 0):>5.2f} "
              f"{r.get('seed', '?'):>5} "
              f"{r.get('num_fractures', 0):>5} "
              f"{r.get('n_elements', 0):>8,} "
              f"{r.get('conn_per_frac', 0):>5.1f} "
              f"{mesh.get('ar5_pct', 0):>5.1f} "
              f"{mesh.get('small_angle_pct', 0):>5.1f} "
              f"{bc.get('n_ex', 0):>3} "
              f"{'Y' if bc.get('has_outflow', False) else 'N':>3} "
              f"{sym:>6s}")

        if verbose and (r["issues"] or r["warnings"]):
            for issue in r["issues"]:
                print(f"  ISSUE: {issue}")
            for warn in r["warnings"]:
                print(f"  WARN:  {warn}")

    print("-" * 100)
    print(f"Total: {len(results)} | "
          f"PASS: {n_pass} | WARN: {n_warn} | FAIL: {n_fail}")

    # ---- Statistical consistency ----
    print(f"\n{'='*80}")
    print("STATISTICAL CONSISTENCY ACROSS REALIZATIONS")
    print(f"{'='*80}\n")

    by_p32 = {}
    for r in results:
        if r["status"] == "FAIL" and "num_fractures" not in r:
            continue
        mult = f"{r.get('p32_mult', 0):.2f}"
        if mult not in by_p32:
            by_p32[mult] = []
        by_p32[mult].append(r)

    print(f"{'P32x':<8s} {'N':>3s} {'Fractures':>16s} {'Elements':>18s} "
          f"{'Conn/frac':>14s} {'Coverage':>10s}")
    print("-" * 80)

    for mult in sorted(by_p32.keys()):
        group = by_p32[mult]
        n = len(group)

        def stats(key):
            vals = [r[key] for r in group if key in r and r[key] > 0]
            if not vals:
                return "—"
            m, s = np.mean(vals), np.std(vals)
            return f"{m:.0f} +/- {s:.0f}"

        covs = []
        for r in group:
            if "coverage" in r:
                covs.append(r["coverage"]["min"])
        cov_str = f"{np.mean(covs):.0%}" if covs else "—"

        print(f"{mult:<8s} {n:>3d} {stats('num_fractures'):>16s} "
              f"{stats('n_elements'):>18s} "
              f"{stats('conn_per_frac'):>14s} "
              f"{cov_str:>10s}")

        # Flag issues
        fracs = [r["num_fractures"] for r in group
                 if "num_fractures" in r]
        if fracs and min(fracs) < THRESH["min_fractures"]:
            print(f"  *** Some realizations have <{THRESH['min_fractures']} "
                  f"fractures — may need higher P32")
        cpfs = [r["conn_per_frac"] for r in group
                if "conn_per_frac" in r and r["conn_per_frac"] > 0]
        if cpfs and np.mean(cpfs) < THRESH["min_conn_per_frac"]:
            print(f"  *** Low average connectivity — may be below "
                  f"percolation threshold")

    # ---- PFLOTRAN readiness ----
    print(f"\n{'='*80}")
    print("PFLOTRAN READINESS")
    print(f"{'='*80}\n")

    ready = [r for r in results
             if r["status"] in ("PASS", "WARN")
             and r.get("boundaries", {}).get("has_outflow", False)]
    not_ready = [r for r in results if r not in ready]

    print(f"  Ready for PFLOTRAN: {len(ready)}/{len(results)}")
    if not_ready:
        print(f"  Not ready:")
        for r in not_ready:
            reasons = r.get("issues", []) + r.get("warnings", [])
            print(f"    {r['name']}: {'; '.join(reasons[:3])}")

    print(f"\n  To run all ready cases:")
    print(f"    python run_pflotran.py --dfn all --nprocs 4")

    return results


# =================================================================
# Generate verification figure
# =================================================================
def make_figure(results):
    valid = [r for r in results
             if "num_fractures" in r and r.get("n_elements", 0) > 0]
    if not valid:
        print("No valid results to plot")
        return

    fig = plt.figure(figsize=(22, 16))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)

    # Color by P32 multiplier
    p32_vals = sorted(set(r.get("p32_mult", 0) for r in valid))
    cmap = plt.cm.viridis(np.linspace(0.15, 0.85, max(len(p32_vals), 1)))
    p32_color = {v: cmap[i] for i, v in enumerate(p32_vals)}

    def get_color(r):
        return p32_color.get(r.get("p32_mult", 0), "gray")

    # 1. Fractures vs P32
    ax = fig.add_subplot(gs[0, 0])
    for r in valid:
        ax.scatter(r.get("p32_mult", 0), r["num_fractures"],
                   color=get_color(r), s=60, edgecolor="k", lw=0.5)
    ax.set_xlabel("P32 multiplier")
    ax.set_ylabel("Connected fractures")
    ax.set_title("Network size vs intensity")
    ax.axhline(THRESH["min_fractures"], color="red", ls="--", lw=1,
               alpha=0.5, label=f"Min = {THRESH['min_fractures']}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # 2. Elements vs P32
    ax = fig.add_subplot(gs[0, 1])
    for r in valid:
        ax.scatter(r.get("p32_mult", 0), r.get("n_elements", 0),
                   color=get_color(r), s=60, edgecolor="k", lw=0.5)
    ax.set_xlabel("P32 multiplier")
    ax.set_ylabel("Mesh elements")
    ax.set_title("Mesh resolution vs intensity")
    ax.grid(alpha=0.3)

    # 3. Connectivity vs P32
    ax = fig.add_subplot(gs[0, 2])
    for r in valid:
        ax.scatter(r.get("p32_mult", 0), r.get("conn_per_frac", 0),
                   color=get_color(r), s=60, edgecolor="k", lw=0.5)
    ax.set_xlabel("P32 multiplier")
    ax.set_ylabel("Connections per fracture")
    ax.set_title("Connectivity vs intensity")
    ax.axhline(THRESH["min_conn_per_frac"], color="red", ls="--", lw=1,
               alpha=0.5, label=f"Min = {THRESH['min_conn_per_frac']}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # 4. Aspect ratio quality
    ax = fig.add_subplot(gs[1, 0])
    ar5 = [r.get("mesh", {}).get("ar5_pct", 0) for r in valid]
    names = [r["name"] for r in valid]
    colors = [get_color(r) for r in valid]
    x = np.arange(len(valid))
    ax.bar(x, ar5, color=colors, edgecolor="k", lw=0.3, alpha=0.8)
    ax.axhline(THRESH["max_ar5_pct"], color="red", ls="--", lw=1, alpha=0.5)
    ax.set_ylabel("% elements with AR > 5")
    ax.set_title("Mesh quality: aspect ratio")
    ax.set_xticks([])
    ax.grid(alpha=0.3, axis="y")

    # 5. Small angle quality
    ax = fig.add_subplot(gs[1, 1])
    sa = [r.get("mesh", {}).get("small_angle_pct", 0) for r in valid]
    ax.bar(x, sa, color=colors, edgecolor="k", lw=0.3, alpha=0.8)
    ax.axhline(THRESH["max_small_angle_pct"], color="red", ls="--",
               lw=1, alpha=0.5)
    ax.set_ylabel("% elements with angle < 10°")
    ax.set_title("Mesh quality: minimum angle")
    ax.set_xticks([])
    ax.grid(alpha=0.3, axis="y")

    # 6. Domain coverage
    ax = fig.add_subplot(gs[1, 2])
    cov = [r.get("coverage", {}).get("min", 0) * 100 for r in valid]
    ax.bar(x, cov, color=colors, edgecolor="k", lw=0.3, alpha=0.8)
    ax.axhline(THRESH["domain_coverage_frac"] * 100, color="red",
               ls="--", lw=1, alpha=0.5)
    ax.set_ylabel("Min axis coverage [%]")
    ax.set_title("Domain coverage")
    ax.set_xticks([])
    ax.grid(alpha=0.3, axis="y")

    # 7. Boundary cells (right_e for outflow)
    ax = fig.add_subplot(gs[2, 0])
    right_cells = [
        r.get("boundaries", {}).get("boundary_right_e.ex", 0)
        for r in valid
    ]
    ax.bar(x, right_cells, color=colors, edgecolor="k", lw=0.3, alpha=0.8)
    ax.axhline(THRESH["min_right_ex_lines"], color="red", ls="--",
               lw=1, alpha=0.5)
    ax.set_ylabel("Outflow boundary cells")
    ax.set_title("Right boundary (outflow)")
    ax.set_xticks([])
    ax.grid(alpha=0.3, axis="y")

    # 8. Pass/Warn/Fail summary
    ax = fig.add_subplot(gs[2, 1])
    status_counts = {"PASS": 0, "WARN": 0, "FAIL": 0}
    for r in results:
        status_counts[r.get("status", "FAIL")] += 1
    bars = ax.bar(
        status_counts.keys(),
        status_counts.values(),
        color=["#2ca02c", "#ff7f0e", "#d62728"],
        edgecolor="k", lw=0.5, alpha=0.8
    )
    for bar, val in zip(bars, status_counts.values()):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                str(val), ha="center", fontsize=14, fontweight="bold")
    ax.set_ylabel("Count")
    ax.set_title("Overall status")
    ax.grid(alpha=0.3, axis="y")

    # 9. Fractures vs connectivity (scatter)
    ax = fig.add_subplot(gs[2, 2])
    for r in valid:
        ax.scatter(r["num_fractures"], r.get("conn_per_frac", 0),
                   color=get_color(r), s=60, edgecolor="k", lw=0.5)
    ax.set_xlabel("Connected fractures")
    ax.set_ylabel("Connections per fracture")
    ax.set_title("Network structure")
    ax.grid(alpha=0.3)

    # Legend for P32 colors
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=p32_color[v], edgecolor="k", lw=0.5,
              label=f"P32 = {v:.2f}x")
        for v in p32_vals
    ]
    fig.legend(handles=legend_elements, loc="upper center",
               ncol=len(p32_vals), fontsize=9,
               bbox_to_anchor=(0.5, 0.98))

    fig.suptitle("DFN Matrix Verification", fontsize=16, y=1.01)
    fig.tight_layout()

    path = os.path.join(DFN_ROOT, "matrix_verification.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {path}")


# =================================================================
# Attempt to fix missing .ex files
# =================================================================
def fix_missing_ex(results):
    from pydfnworks import DFNWORKS

    missing = [r for r in results
               if r.get("boundaries", {}).get("n_ex", 0) < 6
               and r.get("n_elements", 0) > 0]

    if not missing:
        print("All DFNs have boundary files.")
        return

    print(f"\nAttempting to generate .ex files for {len(missing)} DFNs...")

    cwd = os.getcwd()
    for r in missing:
        name = r["name"]
        dfn_dir = os.path.join(DFN_ROOT, name)
        print(f"\n  {name}...")

        if not os.path.exists(os.path.join(dfn_dir, "full_mesh.uge")):
            print(f"    Skip: no UGE file")
            continue

        os.chdir(dfn_dir)
        try:
            DFN = DFNWORKS(jobname=dfn_dir, ncpu=1)
            DFN.inp_file = "full_mesh.inp"
            DFN.uge_file = "full_mesh.uge"
            DFN.flow_solver = "PFLOTRAN"

            with open("params.txt") as f:
                parts = f.readline().strip().split()
                DFN.h = float(parts[1]) if len(parts) >= 2 else 0.2

            DFN.zone2ex(zone_file="all", boundary_cell_area=1.e-1)

            n_ex = len(glob.glob("*.ex"))
            print(f"    Generated {n_ex} .ex files")
        except Exception as e:
            print(f"    Failed: {e}")
        os.chdir(cwd)


# =================================================================
# Main
# =================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Verify all DFNs in the matrix")
    parser.add_argument("--verbose", action="store_true",
        help="Show detailed issues per DFN")
    parser.add_argument("--fix", action="store_true",
        help="Attempt to regenerate missing .ex files")
    args = parser.parse_args()

    if not os.path.exists(DFN_ROOT):
        print(f"No dfn_library/ found. Run prepare_dfn.py first.")
        sys.exit(1)

    results = verify_all(verbose=args.verbose)

    if args.fix:
        fix_missing_ex(results)
        print("\nRe-verifying after fix...")
        results = verify_all(verbose=args.verbose)

    make_figure(results)

    # Final verdict
    n_pass = sum(1 for r in results if r["status"] == "PASS")
    n_warn = sum(1 for r in results if r["status"] == "WARN")
    n_fail = sum(1 for r in results if r["status"] == "FAIL")
    n_ready = sum(1 for r in results
                  if r["status"] in ("PASS", "WARN")
                  and r.get("boundaries", {}).get("has_outflow", False))

    print(f"\n{'='*60}")
    print(f"FINAL VERDICT")
    print(f"{'='*60}")
    print(f"  Total DFNs:       {len(results)}")
    print(f"  PASS:             {n_pass}")
    print(f"  WARN:             {n_warn}")
    print(f"  FAIL:             {n_fail}")
    print(f"  PFLOTRAN-ready:   {n_ready}")
    print(f"  Verification fig: {DFN_ROOT}/matrix_verification.png")

    if n_fail > 0:
        print(f"\n  Action needed: {n_fail} DFNs failed verification.")
        print(f"  Run with --verbose for details.")
        print(f"  Run with --fix to regenerate missing .ex files.")
    elif n_ready == len(results):
        print(f"\n  All {len(results)} DFNs are ready for PFLOTRAN.")
        print(f"  Next: python run_pflotran.py --dfn all --nprocs 4")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()