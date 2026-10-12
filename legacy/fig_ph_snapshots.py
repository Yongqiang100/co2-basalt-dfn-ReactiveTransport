#!/usr/bin/env python3
"""
pH snapshots on the fracture network at 4 timesteps,
matching the x-z projection style of fig_cations.

Run on hpc01:
  python3 src/fig_ph_snapshots.py \
    --case backup_production_vf085_20260909/production/A_p32_100_s613 \
    --mesh runs/A_p32_100_s613
"""
import argparse, glob, os, sys
import numpy as np

PH_PIVOT = 6.3
PH_EMPHASIS = 7.5
CMAP = "RdBu_r"


def read_uge(path):
    with open(path) as f:
        n = int(f.readline().split()[1])
        xyz = np.empty((n, 3))
        for i in range(n):
            t = f.readline().split()
            xyz[i] = (float(t[1]), float(t[2]), float(t[3]))
    return xyz


def load_ph_at_times(run_dir, target_times):
    import h5py
    h5_files = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5")))
                if "dfn_properties" not in f]
    if not h5_files:
        sys.exit(f"No HDF5 output in {run_dir}")

    results = {}
    with h5py.File(h5_files[0], "r") as f:
        all_t, all_g = [], []
        for g in sorted(f.keys()):
            try:
                parts = g.strip().split()
                idx = parts.index("Time")
                all_t.append(float(parts[idx + 1]))
                all_g.append(g)
            except:
                continue
        all_t = np.array(all_t)

        for tt in target_times:
            i = np.argmin(np.abs(all_t - tt))
            grp = f[all_g[i]]
            ph = None
            for k in grp.keys():
                if "pH" in k and "Free" not in k:
                    ph = grp[k][:].flatten()
                    break
            if ph is None:
                for k in grp.keys():
                    if "H+" in k:
                        ph = -np.log10(np.clip(grp[k][:].flatten(), 1e-14, None))
                        break
            if ph is not None:
                results[tt] = (all_t[i], ph)
                print(f"  t={tt} -> {all_t[i]:.4f} yr, pH [{ph.min():.2f}, {ph.max():.2f}]")
    return results


def style(ax, x, z, title=None):
    ax.axvspan(x.min(), x.min() + 0.2 * (x.max() - x.min()),
               color="#4a7fb5", alpha=0.06, lw=0, zorder=0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.4)
    if title:
        ax.set_title(title, fontsize=12, pad=4)


def make_figure(xyz, ph_data, target_times, outfile):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 11, "axes.linewidth": 0.5})
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    x, z = xyz[:, 0], xyz[:, 2]

    ph_min = min(ph.min() for _, ph in ph_data.values())
    ph_max = max(ph.max() for _, ph in ph_data.values())
    norm = TwoSlopeNorm(vmin=max(ph_min, 3.0), vcenter=PH_PIVOT,
                        vmax=min(ph_max, 8.5))

    ncols = len(target_times)
    fig, axes = plt.subplots(1, ncols, figsize=(4.5 * ncols, 4.5))
    if ncols == 1:
        axes = [axes]

    labels = [f"({chr(97 + i)})" for i in range(ncols)]

    sc = None
    for ax, tt, lab in zip(axes, target_times, labels):
        if tt not in ph_data:
            ax.set_title(f"t = {tt} yr — no data")
            ax.axis("off")
            continue

        actual_t, ph = ph_data[tt]

        # Base layer: all cells small
        s_base = np.full(len(ph), 0.6)
        # Emphasise cells near initial pH (un-acidified)
        s_base[ph >= PH_EMPHASIS] = 3.0

        sc = ax.scatter(x, z, c=ph, s=s_base, cmap=CMAP, norm=norm,
                       linewidths=0, rasterized=True)

        if actual_t < 1:
            title = f"t = {actual_t:.2f} yr"
        elif actual_t < 10:
            title = f"t = {actual_t:.1f} yr"
        else:
            title = f"t = {actual_t:.0f} yr"

        style(ax, x, z, title)

        # Panel label and mean pH
        ax.text(0.02, 0.97, lab, transform=ax.transAxes,
               fontsize=13, fontweight="bold", va="top")
        ax.text(0.98, 0.03, f"mean pH {ph.mean():.2f}",
               transform=ax.transAxes, fontsize=8, ha="right",
               bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.85))

        # Axes arrows on first panel
        if ax == axes[0]:
            ax.annotate("", xy=(0.25, -0.04), xytext=(0.04, -0.04),
                        xycoords="axes fraction", textcoords="axes fraction",
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle="-|>", lw=0.9, color="0.15"))
            ax.text(0.285, -0.06, "$x$", transform=ax.transAxes,
                    fontsize=11, va="center", ha="left", color="0.15",
                    clip_on=False)
            ax.annotate("", xy=(-0.045, 0.30), xytext=(-0.045, 0.04),
                        xycoords="axes fraction", textcoords="axes fraction",
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle="-|>", lw=0.9, color="0.15"))
            ax.text(-0.045, 0.325, "$z$", transform=ax.transAxes,
                    fontsize=11, ha="center", va="bottom", color="0.15",
                    clip_on=False)

    fig.subplots_adjust(left=0.04, right=0.97, top=0.93, bottom=0.08,
                        wspace=0.06)

    # Shared colour bar
    if sc is not None:
        cax = fig.add_axes([0.25, 0.02, 0.50, 0.015])
        cbar = fig.colorbar(sc, cax=cax, orientation="horizontal")
        cbar.set_label(f"pH  (pivot {PH_PIVOT})", fontsize=10)
        cbar.ax.tick_params(labelsize=9, length=2, pad=1.5)

    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    fig.savefig(outfile, dpi=300)
    fig.savefig(outfile.rsplit(".", 1)[0] + ".png", dpi=200)
    print(f"\nSaved: {outfile}")
    plt.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="runs/A_p32_100_s613")
    ap.add_argument("--mesh", default=None)
    ap.add_argument("--out", default="figures/letter_ph_snapshots.pdf")
    ap.add_argument("--times", nargs="+", type=float,
                    default=[0.1, 1.0, 10.0, 50.0])
    args = ap.parse_args()

    mesh_dir = args.mesh or args.case
    print(f"Loading mesh from {mesh_dir} ...")
    xyz = read_uge(os.path.join(mesh_dir, "full_mesh.uge"))
    print(f"  {xyz.shape[0]} cells")

    print(f"Loading pH from {args.case} ...")
    ph_data = load_ph_at_times(args.case, args.times)

    make_figure(xyz, ph_data, args.times, args.out)
