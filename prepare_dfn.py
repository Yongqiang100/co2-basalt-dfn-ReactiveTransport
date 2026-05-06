"""
prepare_dfn.py — Generate 5×5 DFN matrix for WRR study

Matrix: 5 P32 levels × 5 stochastic realizations = 25 DFNs
Each DFN includes mesh + .ex boundary files for PFLOTRAN.

Usage:
    conda activate dfnworks
    python prepare_dfn.py single           # 1 baseline DFN
    python prepare_dfn.py matrix           # full 5×5 = 25 DFNs
    python prepare_dfn.py matrix --dry     # show matrix, don't run
    python prepare_dfn.py list             # list all generated DFNs
"""

import os
import sys
import json
import shutil
import argparse
import time as timer
import numpy as np
from pydfnworks import DFNWORKS

OUTPUT_ROOT = os.path.join(os.getcwd(), "dfn_library")

# =================================================================
# Base configuration (20×20×20 m, h=0.2)
# =================================================================
BASE_CONFIG = {
    "domain": [20.0, 20.0, 20.0],
    "h": 0.2,
    "family1": {
        "shape": "ell",
        "distribution": "log_normal",
        "kappa": 20.0,
        "theta": 1.5708,
        "phi": 0.0,
        "log_mean": 1.6,
        "log_std": 0.3,
        "min_radius": 3.0,
        "max_radius": 12.0,
        "p32": 0.8,
        "aspect": 1.5,
        "aperture_mu": -6.9,
        "aperture_sigma": 0.5,
    },
    "family2": {
        "shape": "ell",
        "distribution": "log_normal",
        "kappa": 15.0,
        "theta": 0.0,
        "phi": 0.0,
        "log_mean": 1.8,
        "log_std": 0.4,
        "min_radius": 3.0,
        "max_radius": 15.0,
        "p32": 0.6,
        "aspect": 1.0,
        "aperture_mu": -7.5,
        "aperture_sigma": 0.4,
    },
    "family3": {
        "shape": "ell",
        "distribution": "log_normal",
        "kappa": 20.0,
        "theta": 1.5708,
        "phi": 1.5708,
        "log_mean": 1.6,
        "log_std": 0.3,
        "min_radius": 3.0,
        "max_radius": 12.0,
        "p32": 0.6,
        "aspect": 1.5,
        "aperture_mu": -6.9,
        "aperture_sigma": 0.5,
    },
}

# =================================================================
# Study matrix definition
# =================================================================
P32_MULTIPLIERS = [0.75, 1.00, 1.25, 1.50, 2.00]
SEEDS = [42, 117, 259, 383, 501]

P32_LABELS = {
    0.75: "p32_075",
    1.00: "p32_100",
    1.25: "p32_125",
    1.50: "p32_150",
    2.00: "p32_200",
}


# =================================================================
# Core DFN generation
# =================================================================
def make_config(p32_mult):
    """Create config with scaled P32 values"""
    config = {
        "domain": list(BASE_CONFIG["domain"]),
        "h": BASE_CONFIG["h"],
    }
    for fk in ["family1", "family2", "family3"]:
        config[fk] = dict(BASE_CONFIG[fk])
        config[fk]["p32"] = BASE_CONFIG[fk]["p32"] * p32_mult
    return config


def generate_single_dfn(config, output_name, seed):
    """Generate one DFN with mesh and .ex boundary files"""
    jobname = os.path.join(OUTPUT_ROOT, output_name)
    if os.path.exists(jobname):
        shutil.rmtree(jobname)

    t0 = timer.perf_counter()
    print(f"\n{'='*60}")
    print(f"  DFN: {output_name}")
    print(f"  P32: {config['family1']['p32']:.2f} / "
          f"{config['family2']['p32']:.2f} / "
          f"{config['family3']['p32']:.2f}")
    print(f"  Seed: {seed}")
    print(f"{'='*60}")

    DFN = DFNWORKS(jobname=jobname, ncpu=1)
    DFN.params["domainSize"]["value"] = config["domain"]
    DFN.params["h"]["value"] = config["h"]
    DFN.params["orientationOption"]["value"] = None
    DFN.params["keepOnlyLargestCluster"]["value"] = True
    DFN.params["seed"]["value"] = seed

    for fam_key in ["family1", "family2", "family3"]:
        fam = config[fam_key]
        DFN.add_fracture_family(
            shape=fam["shape"],
            distribution=fam["distribution"],
            kappa=fam["kappa"],
            orientation_distribution="fisher",
            theta=fam["theta"],
            phi=fam["phi"],
            log_mean=fam["log_mean"],
            log_std=fam["log_std"],
            min_radius=fam["min_radius"],
            max_radius=fam["max_radius"],
            number_of_points=8,
            p32=fam["p32"],
            aspect=fam["aspect"],
            beta_distribution=1,
            beta=0.0,
            hy_variable="aperture",
            hy_function="log-normal",
            hy_params={"mu": fam["aperture_mu"],
                       "sigma": fam["aperture_sigma"]},
        )

    try:
        DFN.dfn_gen(output=False)
    except SystemExit:
        print(f"  FAILED: dfn_gen exited for {output_name}")
        return None

    elapsed = timer.perf_counter() - t0

    # Collect summary
    summary = {
        "name": output_name,
        "seed": seed,
        "p32_mult": config["family1"]["p32"] / BASE_CONFIG["family1"]["p32"],
        "p32_values": [config[f]["p32"] for f in ["family1","family2","family3"]],
        "domain": config["domain"],
        "h": config["h"],
        "num_fractures": DFN.num_frac,
        "elapsed_seconds": round(elapsed, 1),
    }

    inp_path = os.path.join(jobname, "full_mesh.inp")
    if os.path.exists(inp_path):
        with open(inp_path) as f:
            header = f.readline().strip().split()
            summary["n_nodes"] = int(header[0])
            summary["n_elements"] = int(header[1])

    uge_path = os.path.join(jobname, "full_mesh.uge")
    if os.path.exists(uge_path):
        summary["uge_size_mb"] = round(os.path.getsize(uge_path) / 1e6, 1)

    conn_path = os.path.join(jobname, "dfnGen_output", "connectivity.dat")
    if os.path.exists(conn_path):
        with open(conn_path) as f:
            lines = f.readlines()
        n_conn = sum(int(line.strip().split()[0]) for line in lines
                     if line.strip()) // 2
        summary["n_connections"] = n_conn
        summary["conn_per_frac"] = round(
            n_conn / max(summary["num_fractures"], 1), 2)

    # Generate .ex boundary files
    print("  Generating .ex boundary files...")
    cwd_save = os.getcwd()
    os.chdir(jobname)
    try:
        DFN2 = DFNWORKS(jobname=jobname, ncpu=1)
        DFN2.inp_file = 'full_mesh.inp'
        DFN2.uge_file = 'full_mesh.uge'
        DFN2.flow_solver = 'PFLOTRAN'
        DFN2.h = config["h"]
        DFN2.zone2ex(zone_file='all', boundary_cell_area=1.e-1)

        import glob
        ex_files = glob.glob("*.ex")
        summary["ex_files"] = len(ex_files)
        ex_sizes = {}
        for ef in ex_files:
            with open(ef) as f:
                n = len(f.readlines())
            ex_sizes[ef] = n
        summary["boundary_cells"] = ex_sizes
    except Exception as e:
        print(f"  .ex generation failed: {e}")
        summary["ex_files"] = 0
    os.chdir(cwd_save)

    # Save summary
    summary_path = os.path.join(jobname, "dfn_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Fractures: {summary.get('num_fractures', '?')}")
    print(f"  Nodes: {summary.get('n_nodes', '?'):,}")
    print(f"  Elements: {summary.get('n_elements', '?'):,}")
    print(f"  Connections: {summary.get('n_connections', '?')} "
          f"({summary.get('conn_per_frac', '?')}/frac)")
    print(f"  Boundary .ex: {summary.get('ex_files', 0)}")
    print(f"  Time: {elapsed:.0f} s")
    return summary


# =================================================================
# Study modes
# =================================================================
def run_single():
    """One baseline DFN at P32×1.0, seed=42"""
    config = make_config(1.0)
    return [generate_single_dfn(config, "baseline", seed=42)]


def run_matrix(dry=False):
    """Full 5×5 matrix: 5 P32 levels × 5 seeds"""
    total = len(P32_MULTIPLIERS) * len(SEEDS)

    print("\n" + "=" * 60)
    print("SIMULATION MATRIX")
    print("=" * 60)
    print(f"  P32 levels:    {P32_MULTIPLIERS}")
    print(f"  Seeds:         {SEEDS}")
    print(f"  Total DFNs:    {total}")
    print(f"  Domain:        {BASE_CONFIG['domain']}")
    print(f"  Mesh h:        {BASE_CONFIG['h']}")
    print()

    # Show matrix
    print(f"  {'P32 mult':<10s}", end="")
    for s in SEEDS:
        print(f"  seed_{s:<5d}", end="")
    print()
    print(f"  {'-'*70}")

    cases = []
    for mult in P32_MULTIPLIERS:
        label = P32_LABELS[mult]
        print(f"  {mult:<10.2f}", end="")
        for seed in SEEDS:
            name = f"{label}_s{seed}"
            cases.append((mult, seed, name))
            print(f"  {name:<11s}", end="")
        print()

    print(f"\n  Total: {len(cases)} DFN realizations")
    print(f"  Est. time: ~{len(cases) * 70 / 60:.0f} min "
          f"(~70s each at 94k elements)")

    if dry:
        print("\n  --dry flag set. No DFNs generated.")
        return []

    print(f"\n  Starting generation...")

    results = []
    failed = []
    for i, (mult, seed, name) in enumerate(cases):
        print(f"\n  [{i+1}/{len(cases)}] {name}")
        config = make_config(mult)
        r = generate_single_dfn(config, name, seed=seed)
        if r is not None:
            results.append(r)
        else:
            failed.append(name)

    # Summary
    print(f"\n\n{'='*60}")
    print(f"MATRIX GENERATION COMPLETE")
    print(f"{'='*60}")
    print(f"  Succeeded: {len(results)}/{len(cases)}")
    if failed:
        print(f"  Failed: {failed}")

    return results


# =================================================================
# List and summary
# =================================================================
def list_dfns():
    if not os.path.exists(OUTPUT_ROOT):
        print("No DFN library. Run: python prepare_dfn.py single")
        return

    # Collect all summaries
    summaries = []
    for entry in sorted(os.listdir(OUTPUT_ROOT)):
        sp = os.path.join(OUTPUT_ROOT, entry, "dfn_summary.json")
        if os.path.exists(sp):
            with open(sp) as f:
                summaries.append(json.load(f))

    if not summaries:
        print("No DFNs found.")
        return

    print(f"\n{'='*90}")
    print(f"DFN Library: {OUTPUT_ROOT}")
    print(f"{'='*90}")
    print(f"{'Name':<20s} {'P32x':>5s} {'Seed':>5s} {'Frac':>5s} "
          f"{'Nodes':>8s} {'Elems':>8s} {'Conn':>5s} {'C/F':>5s} "
          f"{'BCs':>3s} {'Time':>5s}")
    print(f"{'-'*90}")

    # Group by P32 multiplier
    by_p32 = {}
    for s in summaries:
        mult = s.get("p32_mult", 1.0)
        key = f"{mult:.2f}"
        if key not in by_p32:
            by_p32[key] = []
        by_p32[key].append(s)

    for key in sorted(by_p32.keys()):
        for s in sorted(by_p32[key], key=lambda x: x.get("seed", 0)):
            print(f"{s.get('name','?'):<20s} "
                  f"{s.get('p32_mult',0):>5.2f} "
                  f"{s.get('seed','?'):>5} "
                  f"{s.get('num_fractures','?'):>5} "
                  f"{s.get('n_nodes','?'):>8,} "
                  f"{s.get('n_elements','?'):>8,} "
                  f"{s.get('n_connections','?'):>5} "
                  f"{s.get('conn_per_frac','?'):>5} "
                  f"{s.get('ex_files',0):>3} "
                  f"{s.get('elapsed_seconds',0):>4.0f}s")
        print()

    print(f"Total: {len(summaries)} DFNs")

    # Statistics per P32 level
    print(f"\n{'='*70}")
    print(f"Statistics by P32 level")
    print(f"{'='*70}")
    print(f"{'P32x':>6s} {'N':>3s} {'Fracs':>12s} {'Elems':>16s} "
          f"{'Conn/frac':>12s}")
    print(f"{'-'*70}")

    for key in sorted(by_p32.keys()):
        group = by_p32[key]
        n = len(group)
        fracs = [s["num_fractures"] for s in group if "num_fractures" in s]
        elems = [s["n_elements"] for s in group if "n_elements" in s]
        cpf = [s["conn_per_frac"] for s in group if "conn_per_frac" in s]

        def fmt_range(arr):
            if not arr:
                return "?"
            return f"{np.mean(arr):.0f} [{np.min(arr):.0f}-{np.max(arr):.0f}]"

        print(f"{key:>6s} {n:>3d} {fmt_range(fracs):>12s} "
              f"{fmt_range(elems):>16s} {fmt_range(cpf):>12s}")

    print(f"{'='*70}")
    print(f"\nPFLOTRAN meshes at: {OUTPUT_ROOT}/<name>/full_mesh.uge")


def plot_matrix_summary(results):
    """Generate summary comparison plots"""
    if not results:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Organize by P32 multiplier
    by_p32 = {}
    for r in results:
        mult = r.get("p32_mult", 1.0)
        if mult not in by_p32:
            by_p32[mult] = []
        by_p32[mult].append(r)

    mults = sorted(by_p32.keys())
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(mults)))

    # 1. Fractures vs P32
    ax = axes[0, 0]
    for i, mult in enumerate(mults):
        fracs = [r["num_fractures"] for r in by_p32[mult]]
        x = [mult] * len(fracs)
        ax.scatter(x, fracs, color=colors[i], s=60, zorder=5)
    means = [np.mean([r["num_fractures"] for r in by_p32[m]]) for m in mults]
    ax.plot(mults, means, "k--", lw=1.5, alpha=0.5)
    ax.set_xlabel("P32 multiplier")
    ax.set_ylabel("Connected fractures")
    ax.set_title("Network size vs intensity")
    ax.grid(True, alpha=0.3)

    # 2. Elements vs P32
    ax = axes[0, 1]
    for i, mult in enumerate(mults):
        elems = [r.get("n_elements", 0) for r in by_p32[mult]]
        x = [mult] * len(elems)
        ax.scatter(x, elems, color=colors[i], s=60, zorder=5)
    means = [np.mean([r.get("n_elements", 0) for r in by_p32[m]]) for m in mults]
    ax.plot(mults, means, "k--", lw=1.5, alpha=0.5)
    ax.set_xlabel("P32 multiplier")
    ax.set_ylabel("Mesh elements")
    ax.set_title("Mesh resolution vs intensity")
    ax.grid(True, alpha=0.3)

    # 3. Connections per fracture vs P32
    ax = axes[1, 0]
    for i, mult in enumerate(mults):
        cpf = [r.get("conn_per_frac", 0) for r in by_p32[mult]]
        x = [mult] * len(cpf)
        ax.scatter(x, cpf, color=colors[i], s=60, zorder=5)
    means = [np.mean([r.get("conn_per_frac", 0) for r in by_p32[m]]) for m in mults]
    ax.plot(mults, means, "k--", lw=1.5, alpha=0.5)
    ax.set_xlabel("P32 multiplier")
    ax.set_ylabel("Connections per fracture")
    ax.set_title("Connectivity vs intensity")
    ax.grid(True, alpha=0.3)

    # 4. Generation time vs elements
    ax = axes[1, 1]
    for i, mult in enumerate(mults):
        for r in by_p32[mult]:
            ax.scatter(r.get("n_elements", 0), r.get("elapsed_seconds", 0),
                       color=colors[i], s=60, label=f"{mult:.2f}×" if r == by_p32[mult][0] else "",
                       zorder=5)
    ax.set_xlabel("Mesh elements")
    ax.set_ylabel("Generation time [s]")
    ax.set_title("Computational cost")
    ax.legend(fontsize=8, title="P32×")
    ax.grid(True, alpha=0.3)

    fig.suptitle("DFN Matrix Summary: 5 P32 levels × 5 realizations",
                 fontsize=14)
    fig.tight_layout()
    path = os.path.join(OUTPUT_ROOT, "matrix_summary.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {path}")


# =================================================================
# Main
# =================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Generate DFN matrix for WRR study")
    parser.add_argument("mode",
        choices=["single", "matrix", "list"],
        help="single=1 baseline, matrix=5×5, list=show library")
    parser.add_argument("--dry", action="store_true",
        help="Show matrix plan without generating")
    args = parser.parse_args()

    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    if args.mode == "list":
        list_dfns()
        return

    results = []

    if args.mode == "single":
        results = run_single()

    elif args.mode == "matrix":
        results = run_matrix(dry=args.dry)

    if results:
        plot_matrix_summary(results)
        list_dfns()

    print("\nDone.")


if __name__ == "__main__":
    main()