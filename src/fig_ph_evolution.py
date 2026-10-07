#!/usr/bin/env python3
"""
pH spatiotemporal evolution figure.
Panel (a): Space-time heatmap — distance along flow vs time, colour = pH
Panel (b): Spatial pH snapshots at t = 0.1, 1, 10, 50 years

Run on hpc01:
  python3 fig_ph_evolution.py --case runs/A_p32_100_s613
"""
import argparse, os, sys, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm

def load_case(run_dir, mesh_dir=None):
    """Load pH data from PFLOTRAN HDF5 output."""
    import h5py
    
    if mesh_dir is None:
        mesh_dir = run_dir
    
    h5_files = sorted(glob.glob(os.path.join(run_dir, '*.h5')))
    h5_files = [f for f in h5_files if 'dfn_properties' not in f]
    if not h5_files:
        sys.exit(f"No HDF5 output in {run_dir}")
    
    # Read mesh geometry (cell centres)
    uge = os.path.join(mesh_dir, 'full_mesh.uge')
    with open(uge) as f:
        n = int(f.readline().split()[1])
        coords = []
        for _ in range(n):
            parts = f.readline().split()
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    coords = np.array(coords)
    
    # Read pH at all timesteps
    with h5py.File(h5_files[0], 'r') as f:
        groups = sorted(f.keys())
        
        times = []
        ph_all = []
        for g in groups:
            grp = f[g]
            # Parse time from group name: "  N Time  X.XXE+XX y"
            try:
                parts = g.strip().split()
                t_idx = parts.index('Time')
                t = float(parts[t_idx + 1])
            except:
                continue
            
            # Find pH field
            ph_key = None
            for k in grp.keys():
                kl = k.lower()
                if 'ph' in kl and 'free' not in kl:
                    ph_key = k
                    break
            
            if ph_key is None:
                # Compute from H+
                for k in grp.keys():
                    if 'Total H+' in k or 'Free H+' in k:
                        h_conc = grp[k][:]
                        ph = -np.log10(np.clip(h_conc, 1e-14, None))
                        times.append(t)
                        ph_all.append(ph)
                        break
                continue
            
            times.append(t)
            ph_all.append(grp[ph_key][:].flatten())
    
    times = np.array(times)
    ph_all = np.array(ph_all)
    
    # Sort by time
    order = np.argsort(times)
    times = times[order]
    ph_all = ph_all[order]
    
    return coords, times, ph_all


def make_figure(coords, times, ph_all, outfile):
    """Create the space-time pH figure."""
    
    # Project onto flow direction (x-axis)
    x = coords[:, 0]
    x_min, x_max = x.min(), x.max()
    x_range = x_max - x_min
    x_norm = (x - x_min) / x_range  # 0 = inlet, 1 = outlet
    
    # Bin cells by normalised distance
    nbins = 60
    bin_edges = np.linspace(0, 1, nbins + 1)
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    cell_bin = np.digitize(x_norm, bin_edges) - 1
    cell_bin = np.clip(cell_bin, 0, nbins - 1)
    
    # Build space-time pH matrix (ntimes x nbins)
    ph_matrix = np.full((len(times), nbins), np.nan)
    for ti in range(len(times)):
        for bi in range(nbins):
            mask = cell_bin == bi
            if mask.sum() > 0:
                ph_matrix[ti, bi] = np.mean(ph_all[ti, mask])
    
    fig = plt.figure(figsize=(14, 5.5))
    gs = GridSpec(1, 2, width_ratios=[1.2, 1], wspace=0.35)
    
    # === Panel (a): Space-time heatmap ===
    ax1 = fig.add_subplot(gs[0])
    
    # Use diverging colourmap centred at ~5.5 (approximate carbonate stability threshold)
    vmin, vmax = 3.0, 8.0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=5.5, vmax=vmax)
    
    # Pcolormesh: x = distance, y = time
    T, X = np.meshgrid(times, bin_centres, indexing='ij')
    im = ax1.pcolormesh(X, T, ph_matrix, cmap='RdYlBu', norm=norm,
                        shading='auto', rasterized=True)
    
    cb = fig.colorbar(im, ax=ax1, label='pH', shrink=0.9)
    ax1.set_xlabel('Normalised distance along flow direction', fontsize=11)
    ax1.set_ylabel('Time (years)', fontsize=11)
    ax1.set_title('(a) pH space-time evolution', fontsize=12)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, times.max())
    
    # Mark snapshot times
    snapshot_times = [0.1, 1.0, 10.0, 50.0]
    for st in snapshot_times:
        if st <= times.max():
            ax1.axhline(st, color='white', linewidth=0.5, linestyle='--', alpha=0.7)
    
    # === Panel (b): Spatial snapshots ===
    ax2 = fig.add_subplot(gs[1])
    
    colours = ['#e41a1c', '#ff7f00', '#4daf4a', '#377eb8']
    linestyles = ['-', '-', '-', '-']
    
    for st, col, ls in zip(snapshot_times, colours, linestyles):
        idx = np.argmin(np.abs(times - st))
        actual_t = times[idx]
        
        # Get binned pH profile
        profile = ph_matrix[idx, :]
        ax2.plot(bin_centres, profile, ls, color=col, linewidth=1.8,
                label=f't = {actual_t:.1f} yr')
    
    ax2.axhline(y=7.5, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)
    ax2.axhline(y=3.4, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)
    ax2.text(0.95, 7.6, 'initial pH = 7.5', fontsize=8, ha='right', color='gray')
    ax2.text(0.95, 3.5, 'injectate pH = 3.4', fontsize=8, ha='right', color='gray')
    
    ax2.set_xlabel('Normalised distance along flow direction', fontsize=11)
    ax2.set_ylabel('pH', fontsize=11)
    ax2.set_title('(b) Spatial pH profiles', fontsize=12)
    ax2.legend(fontsize=9, loc='center right')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(2.5, 8.5)
    
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"Saved: {outfile}")
    plt.close()


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='runs/A_p32_100_s613')
    ap.add_argument('--mesh', default=None,
                    help='Directory containing full_mesh.uge (defaults to --case)')
    ap.add_argument('--out', default='figures/letter_ph_evolution.pdf')
    args = ap.parse_args()
    
    mesh_dir = args.mesh or args.case
    print(f"Loading {args.case} (mesh from {mesh_dir}) ...")
    coords, times, ph_all = load_case(args.case, mesh_dir=mesh_dir)
    print(f"  {len(times)} timesteps, {coords.shape[0]} cells")
    print(f"  Time range: {times.min():.4f} to {times.max():.1f} yr")
    
    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    make_figure(coords, times, ph_all, args.out)
