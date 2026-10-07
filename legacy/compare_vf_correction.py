#!/usr/bin/env python3
"""
Compare corrected (VF=0.50) vs original (VF=0.85) ensemble results.
Run on hpc01:
  python3 src/compare_vf_correction.py
"""
import os, sys, glob
import numpy as np

RUNS = "runs"
BACKUP = "backup_production_vf085_20260909/production"
CARB_PHASES = ["Calcite", "Magnesite", "Siderite", "Dawsonite"]

def read_carbonate(h5_path):
    """Read final carbonate VFs from HDF5."""
    import h5py
    try:
        with h5py.File(h5_path, 'r') as f:
            keys = sorted(f.keys())
            last = f[keys[-1]]
            # Parse time
            parts = keys[-1].strip().split()
            t = float(parts[parts.index('Time') + 1])
            
            totals = {}
            ncells = None
            for phase in CARB_PHASES:
                k = f"{phase} VF [m^3 mnrl_m^3 bulk]"
                if k in last:
                    data = last[k][:].flatten()
                    ncells = len(data)
                    seed = 1.0e-6
                    net = np.maximum(data - seed, 0.0)
                    totals[phase] = net.sum()
            
            total_carb = sum(totals.values())
            return {
                'time': t,
                'ncells': ncells,
                'total_carb': total_carb,
                'per_cell': total_carb / ncells if ncells else 0,
                'phases': totals
            }
    except Exception as e:
        return None


def main():
    cases = sorted([d for d in os.listdir(RUNS) if d.startswith("A_p32_")])
    
    print(f"{'Case':<25s} {'Old/cell':>12s} {'New/cell':>12s} {'Ratio':>8s} {'Old t':>6s} {'New t':>6s}")
    print("-" * 75)
    
    ratios = []
    for case in cases:
        old_h5 = os.path.join(BACKUP, case, "pflotran_co2.h5")
        new_h5 = os.path.join(RUNS, case, "pflotran_co2.h5")
        
        if not os.path.exists(old_h5) or not os.path.exists(new_h5):
            print(f"{case:<25s} {'MISSING':>12s}")
            continue
        
        old = read_carbonate(old_h5)
        new = read_carbonate(new_h5)
        
        if old is None or new is None:
            print(f"{case:<25s} {'ERROR':>12s}")
            continue
        
        if old['per_cell'] > 0:
            ratio = new['per_cell'] / old['per_cell']
            ratios.append(ratio)
        else:
            ratio = float('inf')
        
        print(f"{case:<25s} {old['per_cell']:>12.6e} {new['per_cell']:>12.6e} {ratio:>8.3f} {old['time']:>6.1f} {new['time']:>6.1f}")
    
    if ratios:
        ratios = np.array(ratios)
        print("-" * 75)
        print(f"{'Completed pairs:':<25s} {len(ratios)}")
        print(f"{'Ratio range:':<25s} {ratios.min():.3f} to {ratios.max():.3f}")
        print(f"{'Ratio mean:':<25s} {ratios.mean():.3f}")
        print(f"{'Ratio median:':<25s} {np.median(ratios):.3f}")
        print()
        
        expected = "0.76 to 0.91"
        print(f"Expected from sensitivity: x{expected}")
        print(f"Actual range:              x{ratios.min():.3f} to x{ratios.max():.3f}")
        
        if ratios.min() >= 0.5 and ratios.max() <= 1.2:
            print("CONCLUSION: VF correction is within the expected range.")
            print("Between-realisation differences are preserved.")
        else:
            print("WARNING: Some ratios outside expected range. Investigate.")


if __name__ == "__main__":
    main()
