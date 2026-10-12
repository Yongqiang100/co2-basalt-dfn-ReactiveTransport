#!/usr/bin/env python3
"""
Domain-average pH at 50 years, across the production ensemble.

The letter tells R2 that the domain-average pH is 6.07 and that the distinction
between that and the near-injection value of 3.4 is stated in the revised text.
Neither 6.07 nor "domain-average pH" appears in the manuscript, so the figure
has to be reproduced before the claim can stand.

Reported two ways, because the answer depends on the weighting: an unweighted
mean over cells, and a mean weighted by cell volume. The volume-weighted value
is the physically meaningful one.

    python3 src/domain_mean_ph.py
"""
import glob
import os
import sys

import numpy as np
import h5py


def cell_volumes(run_dir, n):
    """Cell volumes from the UGE file, or None if it is unavailable."""
    uge = os.path.join(run_dir, "full_mesh.uge")
    if not os.path.isfile(uge):
        return None
    with open(uge) as f:
        head = f.readline().split()
        if len(head) < 2:
            return None
        ncell = int(head[1])
        vol = np.empty(ncell)
        for i in range(ncell):
            vol[i] = float(f.readline().split()[4])
    return vol[:n] if len(vol) >= n else None


def mean_ph(run_dir):
    with h5py.File(os.path.join(run_dir, "pflotran_co2.h5"), "r") as f:
        g = f[sorted(f.keys())[-1]]
        k = next((x for x in g if x.startswith("pH")), None)
        if k is None:
            return None
        ph = np.asarray(g[k][:], float).flatten()
    vol = cell_volumes(run_dir, len(ph))
    if vol is None or vol.sum() <= 0:
        return float(ph.mean()), None
    return float(ph.mean()), float((ph * vol).sum() / vol.sum())


def main():
    runs = sorted(d for d in glob.glob("runs/A_p32_*") if os.path.isdir(d))
    unw, vw = [], []
    for d in runs:
        try:
            r = mean_ph(d)
        except (OSError, KeyError):
            continue
        if r is None:
            continue
        unw.append(r[0])
        if r[1] is not None:
            vw.append(r[1])

    print(f"{len(unw)} realization(s) read\n")
    if unw:
        print("unweighted mean over cells")
        print(f"  ensemble mean   {np.mean(unw):.2f}")
        print(f"  ensemble median {np.median(unw):.2f}")
        print(f"  range           {min(unw):.2f} to {max(unw):.2f}")
    if vw:
        print("\nvolume-weighted")
        print(f"  ensemble mean   {np.mean(vw):.2f}")
        print(f"  ensemble median {np.median(vw):.2f}")
        print(f"  range           {min(vw):.2f} to {max(vw):.2f}")
    else:
        print("\nno UGE files found, so no volume weighting")

    print("\nThe letter states 6.07. Compare the values above before that number")
    print("goes into the manuscript.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
