#!/usr/bin/env python3
"""
Carbon budget per run: CO2 mineralized as carbonate, CO2 injected, efficiency.

    python3 src/carbon_budget.py --prefix A_
    python3 src/carbon_budget.py --prefix A_feedback__ --root runs

Mineralized: over all cells, net carbonate volume fraction at the last snapshot
(final VF minus the 1e-6 seed, never below 0) x cell volume / molar volume,
summed over calcite, magnesite, siderite and dawsonite (one C each). Molar
volumes are read from the run's own hanford.dat. Also reported: the share of
that CO2 in cells whose carbonate volume fraction exceeds 1 (physically
impossible: more mineral than the cell holds).

Injected: the water PFLOTRAN actually injected, from its mass-balance file
(pflotran_co2-mas.dat, "carbonated_water_injector Water [kg]", cumulative), times
the CO2(aq) of the injected water (constraint co2_rich_water, mol/kg). Only if
that file is missing is the deck's RATE LIST integrated instead, as a step
function, which underestimates short injections that end inside the start-up
ramp. Treating 1 kg of injected water as 1 kg of solvent is a small approximation.
"""
import argparse, glob, os, re, sys
import numpy as np, h5py

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
CARB = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
SEED, M_CO2, YEAR = 1.0e-6, 44.0095e-3, 365.25 * 86400.0


def num(t):
    return float(t.replace("d", "e").replace("D", "e"))


def molar_volumes(d):
    db = os.path.join(d, "hanford.dat")
    if not os.path.isfile(db):
        db = os.path.join(os.path.dirname(REV), "database", "hanford.dat")
    mv = {}
    for line in open(db, errors="replace"):
        t = line.split()
        if len(t) > 1 and t[0].strip("'") in CARB and t[0].startswith("'") and t[0].strip("'") not in mv:
            try:
                mv[t[0].strip("'")] = num(t[1]) * 1e-6          # cm3/mol -> m3/mol
            except ValueError:
                pass
    missing = [m for m in CARB if m not in mv]
    if missing:
        raise RuntimeError(f"molar volume not found in {db}: {missing}")
    return mv


def cell_volumes(d):
    with open(os.path.join(d, "full_mesh.uge")) as f:
        n = int(f.readline().split()[1])
        return np.array([float(f.readline().split()[4]) for _ in range(n)])


def _parts(d, stem, ext):
    """A continued run's files in time order: <stem>_part1<ext>, _part2, ..., then <stem><ext>."""
    parts = sorted(glob.glob(os.path.join(d, f"{stem}_part*{ext}")),
                   key=lambda p: int(re.search(r"_part(\d+)", p).group(1)))
    main = os.path.join(d, stem + ext)
    return parts + ([main] if os.path.isfile(main) else [])


def snapshots(d):
    """[(t, h5 path, group name)] over all parts of a run, in time order; where a time
    appears in two parts, the later part's copy is kept."""
    out = {}
    for p in _parts(d, "pflotran_co2", ".h5"):
        with h5py.File(p, "r") as f:
            for k in f:
                if "Time" in k:
                    out[round(float(k.split("Time")[1].split()[0]), 9)] = (p, k)
    return [(t, p, k) for t, (p, k) in sorted(out.items())]


def _mas_rows(p):
    with open(p, errors="replace") as f:
        cols = [c.strip().strip('"') for c in f.readline().split(",")]
        rows = [l.split() for l in f if l.strip()]
    return cols, rows


def carbonate_mas(d):
    """Net carbonate (mol) from PFLOTRAN's mass balance, across all parts of a run: the
    last row of the last part minus the first row of the first part. None if unavailable."""
    ps = _parts(d, "pflotran_co2-mas", ".dat")
    if not ps:
        return None
    (c0, r0), (c1, r1) = _mas_rows(ps[0]), _mas_rows(ps[-1])
    i0 = [i for i, c in enumerate(c0) if any(f" {m} Total Mass" in c for m in CARB)]
    i1 = [i for i, c in enumerate(c1) if any(f" {m} Total Mass" in c for m in CARB)]
    if not (i0 and i1 and r0 and r1):
        return None
    return sum(float(r1[-1][i]) for i in i1) - sum(float(r0[0][i]) for i in i0)


def injected_water_mas(d):
    """Cumulative injected water (kg) from PFLOTRAN's mass-balance file, or None."""
    total, prev_last = 0.0, 0.0
    ps = _parts(d, "pflotran_co2-mas", ".dat")
    if not ps:
        return None
    for p in ps:
        cols, rows = _mas_rows(p)
        k = next((i for i, c in enumerate(cols) if "injector" in c.lower() and "water [kg]" in c.lower()), None)
        if k is None or not rows:
            return None
        try:
            first, last = abs(float(rows[0][k])), abs(float(rows[-1][k]))
        except (IndexError, ValueError):
            return None
        # a continuation either carries the cumulative total on, or starts again from zero
        total = last if (prev_last and first >= 0.999 * prev_last) or not prev_last else total + last
        prev_last = last
    return total


def injected_mol(deck, t_end):
    s = open(deck).read()
    c = re.search(r"CONSTRAINT\s+co2_rich_water.*?CO2\(aq\)\s+([0-9.eEdD+-]+)", s, re.S | re.I)
    w = injected_water_mas(os.path.dirname(deck))
    if w is not None and c:
        return w * num(c.group(1))
    m = re.search(r"RATE\s+LIST(.*?)\n\s*/", s, re.S | re.I)
    rows = [(num(a), num(b)) for a, b in re.findall(r"^\s*([0-9.eEdD+-]+)\s+([0-9.eEdD+-]+)\s*$", m.group(1), re.M)] if m else []
    c = re.search(r"CONSTRAINT\s+co2_rich_water.*?CO2\(aq\)\s+([0-9.eEdD+-]+)", s, re.S | re.I)
    if not rows or not c:
        return None
    kg = 0.0
    for i, (t, r) in enumerate(rows):
        t_next = rows[i + 1][0] if i + 1 < len(rows) else t_end
        kg += r * max(0.0, min(t_next, t_end) - min(t, t_end)) * YEAR
    return kg * num(c.group(1))


def budget(d):
    deck = next(iter(sorted(glob.glob(os.path.join(d, "*.in")))), None)
    snaps = snapshots(d) if deck else []
    if not snaps:
        return None
    t_end, p, k = snaps[-1]
    with h5py.File(p, "r") as f:
        g = f[k]
        vf = {m: np.clip(np.asarray(g[n][:], float).ravel() - SEED, 0, None)
              for n in g for m in CARB if n.split(" VF")[0] == m and " VF" in n}
    vol, mv = cell_volumes(d), molar_volumes(d)
    mol_cell = sum(vf[m] * vol / mv[m] for m in vf)
    total_vf = sum(vf.values())
    mol = float(mol_cell.sum())
    over = float(mol_cell[total_vf > 1].sum())
    inj = injected_mol(deck, t_end)
    return dict(t=t_end, kg=mol * M_CO2, over=100 * over / mol if mol else 0.0,
                inj=inj * M_CO2 if inj else None, eff=100 * mol / inj if inj else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", default="A_")
    ap.add_argument("--root", default="runs_gravityoff")
    a = ap.parse_args()
    root = a.root if os.path.isabs(a.root) else os.path.join(REV, a.root)
    dirs = sorted(d for d in glob.glob(os.path.join(root, a.prefix + "p32_*")) if os.path.isdir(d))
    if not dirs:
        sys.exit(f"no {a.prefix}p32_* directories in {root}")
    W = max(22, max(len(os.path.basename(d)) for d in dirs) + 2)
    print(f"{'run':<{W}}{'years':>6}{'CO2 mineralized kg':>20}{'in VF>1 %':>11}{'CO2 injected kg':>17}{'efficiency %':>14}")
    res = []
    for d in dirs:
        try:
            b = budget(d)
        except Exception as e:
            print(f"{os.path.basename(d):<{W}}  error: {e}"); continue
        if b is None:
            print(f"{os.path.basename(d):<{W}}  no output"); continue
        res.append(b)
        inj = f"{b['inj']:.4g}" if b["inj"] else "-"
        eff = f"{b['eff']:.4g}" if b["eff"] is not None else "-"
        print(f"{os.path.basename(d):<{W}}{b['t']:>6.3g}{b['kg']:>20.4g}{b['over']:>11.1f}{inj:>17}{eff:>14}")
    fin = [b for b in res if b["t"] >= 49.99]
    if fin:
        kg = np.array([b["kg"] for b in fin]); eff = np.array([b["eff"] for b in fin if b["eff"] is not None])
        ov = np.array([b["over"] for b in fin])
        print(f"\n{len(fin)} runs at 50 y: CO2 mineralized per run median {np.median(kg):.4g} kg (range {kg.min():.3g} to {kg.max():.3g}); "
              f"efficiency median {np.median(eff):.3g}% (range {eff.min():.3g} to {eff.max():.3g}%)")
        print(f"share of the mineralized CO2 in cells with carbonate VF above 1: median {np.median(ov):.1f}%")


if __name__ == "__main__":
    main()
