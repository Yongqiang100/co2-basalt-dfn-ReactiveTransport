#!/usr/bin/env python3
"""
Trapping efficiency as a function of time, across the production ensemble.

Efficiency is the carbon fixed as carbonate divided by the carbon injected, and
both terms grow during the simulation. The denominator grows at a constant rate
once the injection ramp ends; the numerator grows at a rate that falls as the
rock upstream is stripped. Efficiency is therefore a curve rather than a single
figure, and reporting only the value at fifty years reports the endpoint of a
decline rather than what the rock achieved.

On the highest-carbonate network the efficiency rises through the first two
years, holds from two to seven, and then falls by 42 per cent to the value at
fifty years. This script computes the same curve for every production
realization and reports, for each, the peak efficiency, the time at which it
occurs, and the value at fifty years.

The carbonate comes from the volume fractions in the HDF5 output, converted
through the molar volumes read from hanford.dat, with the seeded volume
fraction removed. That route reproduces the efficiency column in finalised.csv
to seven figures. The injected carbon is the water injection rate multiplied by
the elapsed time and the injectate concentration, which agrees with the
injected mass in the PFLOTRAN mass-balance file to one part in a thousand.

    python3 src/efficiency_time.py
    python3 src/efficiency_time.py --block A --csv efficiency_time.csv
    python3 src/efficiency_time.py --detail A_p32_100_s1397
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import sys

import numpy as np
import h5py

MOLAR_VOL_CM3 = {
    "Calcite": 36.9340,
    "Magnesite": 28.0180,
    "Siderite": 29.3780,
    "Dawsonite": 58.5366,
}
SEED = 1.0e-6
YEAR = 3.156e7
CONC = 0.82                     # mol of carbon per kg of water

# The injection rate is scaled to each network's pore volume and varies by a
# factor of three across the ensemble, from 5.75e-03 to 1.80e-02 kg/s. An
# earlier version of this script used a single constant and therefore
# misstated the efficiency for every network but one. The injected mass is now
# read from each run's own mass-balance file, and the rate is recovered from
# the deck only if that file is missing.


def read_volumes(uge):
    with open(uge) as f:
        n = int(f.readline().split()[1])
        v = np.empty(n)
        for i in range(n):
            v[i] = float(f.readline().split()[4])
    return v


def carbonate_series(h5, vol):
    """Carbonate in moles at every snapshot."""
    out = []
    with h5py.File(h5, "r") as f:
        tg = sorted((float(m.group(1)), k) for k in f.keys()
                    if (m := re.search(r"([-\d.]+[eE][-+]?\d+|[-\d.]+)\s*y", k)))
        for t, k in tg:
            g = f[k]
            tot = 0.0
            seen = False
            for p, mv in MOLAR_VOL_CM3.items():
                kk = next((x for x in g if x.startswith(p + " VF")), None)
                if kk is None:
                    continue
                seen = True
                v = np.maximum(np.asarray(g[kk][:], float).flatten() - SEED, 0.0)
                tot += float((v * vol).sum()) / (mv * 1e-6)
            if seen:
                out.append((t, tot))
    return out


def injection_rate(run_dir):
    """The plateau injection rate in kg/s, from this network's own deck."""
    deck = os.path.join(run_dir, "pflotran_co2.in")
    if not os.path.isfile(deck):
        return None
    with open(deck) as f:
        txt = f.read()
    m = re.search(r"RATE LIST(.*?)/", txt, re.S)
    if not m:
        return None
    rates = []
    for line in m.group(1).splitlines():
        p_ = line.split()
        if len(p_) == 2:
            try:
                rates.append((float(p_[0]), float(p_[1])))
            except ValueError:
                pass
    return max((r for _, r in rates), default=None)


def injected_series(run_dir, times):
    """Carbon injected by each snapshot time, in moles.

    Uses the final injected mass from the mass-balance file to fix the rate,
    which accounts for the ramp over the first hundredth of a year, and scales
    linearly to earlier times.
    """
    mas = os.path.join(run_dir, "pflotran_co2-mas.dat")
    kg_final = t_final = None
    if os.path.isfile(mas):
        with open(mas) as f:
            header = [h.strip().strip('"') for h in f.readline().split(",")]
            last = None
            for line in f:
                if line.strip():
                    last = line
        if last:
            vals = [float(x) for x in last.split()]
            t_final = vals[0]
            for n_, v_ in zip(header, vals):
                if "injector water [kg]" in n_.lower():
                    kg_final = v_
                    break
    if kg_final and t_final:
        rate = kg_final / (t_final * YEAR)
    else:
        rate = injection_rate(run_dir)
        if rate is None:
            return None
    return {t: rate * t * YEAR * CONC for t in times}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--block", default="A")
    ap.add_argument("--detail", help="print the full curve for this network")
    ap.add_argument("--csv")
    a = ap.parse_args()

    dirs = sorted(glob.glob(os.path.join(a.runs, f"{a.block}_p32_*")))
    if a.detail:
        dirs = [d for d in dirs if os.path.basename(d) == a.detail] or \
               [os.path.join(a.runs, a.detail)]
    if not dirs:
        print(f"no {a.block}_p32_* directories under {a.runs}")
        return 1

    rows = []
    for d in dirs:
        h5 = os.path.join(d, "pflotran_co2.h5")
        uge = os.path.join(d, "full_mesh.uge")
        if not (os.path.isfile(h5) and os.path.isfile(uge)):
            continue
        try:
            vol = read_volumes(uge)
            ser = carbonate_series(h5, vol)
        except (OSError, KeyError, ValueError):
            continue
        inj = injected_series(d, [t for t, _ in ser])
        if inj is None:
            continue
        pts = [(t, c, 100.0 * c / inj[t])
               for t, c in ser if t > 0 and c > 0 and inj.get(t, 0) > 0]
        if not pts:
            rows.append((os.path.basename(d), np.nan, np.nan, np.nan, np.nan))
            continue
        if a.detail:
            print(f"{os.path.basename(d)}")
            print("%8s %12s %14s %12s" %
                  ("t (y)", "carb (mol)", "C injected", "eff (%)"))
            for t, c, e in pts:
                print("%8.2f %12.3f %14.4e %12.6f"
                      % (t, c, inj[t], e))
            print()
        eff = np.array([e for _, _, e in pts])
        tt = np.array([t for t, _, _ in pts])
        k = int(np.argmax(eff))
        final = eff[-1]
        rows.append((os.path.basename(d), tt[k], eff[k], final,
                     eff[k] / final if final > 0 else np.nan))

    ok = [r for r in rows if np.isfinite(r[2])]
    if not ok:
        print("no realization produced carbonate")
        return 1
    ok.sort(key=lambda r: -r[2])

    print(f"{len(ok)} realizations with carbonate, of {len(rows)}")
    print()
    print("%-22s %10s %12s %12s %8s" %
          ("network", "peak at y", "peak eff %", "50 y eff %", "ratio"))
    print("-" * 70)
    for nm, tp, ep, ef, r in ok[:12]:
        print("%-22s %10.2f %12.6f %12.6f %8.2f" % (nm, tp, ep, ef, r))
    if len(ok) > 12:
        print(f"... {len(ok)-12} more")

    tp = np.array([r[1] for r in ok])
    ep = np.array([r[2] for r in ok])
    ef = np.array([r[3] for r in ok])
    rr = np.array([r[4] for r in ok])
    print()
    print(f"peak efficiency   : highest {ep.max():.6f} %, "
          f"median {np.median(ep):.6f} %")
    print(f"efficiency at 50 y: highest {ef.max():.6f} %, "
          f"median {np.median(ef):.6f} %")
    print(f"time of the peak  : median {np.median(tp):.1f} y, "
          f"range {tp.min():.2f} to {tp.max():.1f} y")
    print(f"peak over final   : median {np.median(rr):.2f}, "
          f"max {rr.max():.2f}")
    print()
    late = int((tp >= 50).sum())
    print(f"realizations whose efficiency is still rising at 50 y: "
          f"{late} of {len(ok)}")
    if late < len(ok) / 2:
        print("  In most realizations the efficiency peaks before the end of")
        print("  the simulation and declines thereafter, because the injected")
        print("  carbon continues to accumulate while the rate of carbonate")
        print("  formation falls.")

    if a.csv:
        import csv as _csv
        with open(a.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["run_id", "peak_time_y", "peak_efficiency_pct",
                        "efficiency_50y_pct", "peak_over_final"])
            w.writerows(rows)
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
