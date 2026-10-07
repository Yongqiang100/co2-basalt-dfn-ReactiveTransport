#!/usr/bin/env python3
"""Close the outstanding numerical gaps in the revised manuscript.

Reports, per realisation and by P32 level: pH at 50 yr, primary mineral
dissolution fractions, trapping efficiency, carbonate amount and assemblage,
and pore volumes flushed.

Nothing is assumed: molar volumes are parsed from pflotran_co2.out, injected
mass is read from the mass balance file rather than reconstructed from the
deck rate, and runs that did not reach the final time are skipped.
"""
from __future__ import annotations
import sys, os, re, csv, glob, argparse
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

CARB    = ("Calcite", "Magnesite", "Siderite", "Dawsonite")
PRIMARY = ("Anorthite", "Albite", "Diopside", "Forsterite", "Fayalite", "Enstatite")
SEED_VF, PHI, CO2_MOLAL = 1e-6, 0.50, 0.82
CO2_PER_FORMULA = {"Calcite": 1, "Magnesite": 1, "Siderite": 1, "Dawsonite": 1}

def molar_volumes(out_path):
    if not os.path.isfile(out_path): return {}
    mv, pending = {}, None
    with open(out_path, errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s in CARB or s in PRIMARY:
                pending = s
            elif pending and s.startswith("Molar Volume:"):
                m = re.search(r"([0-9.]+E[+-][0-9]+)", s)
                if m: mv[pending] = float(m.group(1))
                pending = None
    return mv

def injected_kg(mas_path):
    if not mas_path or not os.path.isfile(mas_path): return None
    txt = open(mas_path, errors="ignore").read().splitlines()
    if len(txt) < 2: return None
    hdr = [h.strip().strip('"') for h in txt[0].split(",")]
    try:
        col = next(i for i, h in enumerate(hdr) if "injector Water [kg]" in h)
    except StopIteration:
        return None
    rows = [l.split() for l in txt[1:] if l.strip()]
    if not rows: return None
    try: return float(rows[-1][col])
    except (IndexError, ValueError): return None

def cell_volumes(uge_path):
    with open(uge_path) as f:
        n = int(f.readline().split()[1])
        return np.array([float(f.readline().split()[4]) for _ in range(n)])

def analyse(run_dir, final_year=50.0):
    import h5py
    from validate_run import time_groups
    h5 = [f for f in sorted(glob.glob(os.path.join(run_dir, "*.h5")))
          if "dfn_properties" not in f]
    uge = os.path.join(run_dir, "full_mesh.uge")
    if not h5 or not os.path.isfile(uge): return None, "no output"
    vol = cell_volumes(uge)
    mv  = molar_volumes(os.path.join(run_dir, "pflotran_co2.out"))
    mas = glob.glob(os.path.join(run_dir, "*-mas.dat"))
    inj = injected_kg(mas[0] if mas else None)

    with h5py.File(h5[0], "r") as f:
        tg = time_groups(f)
        if not tg: return None, "no time groups"
        if tg[-1][0] < final_year * 0.999:
            return None, f"incomplete ({tg[-1][0]:g} yr)"
        g0, g1 = f[tg[0][1]], f[tg[-1][1]]
        def field(grp, prefix):
            k = next((x for x in grp.keys() if x.startswith(prefix)), None)
            return np.asarray(grp[k][:], dtype=float) if k else None
        ph = field(g1, "pH")
        if ph is None:
            h = field(g1, "Total H+")
            ph = -np.log10(np.clip(h, 1e-30, None)) if h is not None else None
        carb_vf, carb_vol, carb_mol, n = {}, {}, 0.0, len(vol)
        for m in CARB:
            a = field(g1, f"{m} VF")
            if a is None:
                carb_vf[m] = carb_vol[m] = 0.0; continue
            if a.size != vol.size:
                return None, f"length mismatch {a.size} vs {vol.size}"
            net = np.clip(a - SEED_VF, 0.0, None)
            carb_vf[m]  = float(net.sum())
            carb_vol[m] = float((net * vol).sum())
            if m in mv: carb_mol += carb_vol[m] / mv[m] * CO2_PER_FORMULA[m]
        diss = {}
        for m in PRIMARY:
            a0, a1 = field(g0, f"{m} VF"), field(g1, f"{m} VF")
            if a0 is None or a1 is None: diss[m] = None; continue
            v0, v1 = float((a0*vol).sum()), float((a1*vol).sum())
            diss[m] = 100.0*(v0-v1)/v0 if v0 > 0 else None

    pore = float((vol * PHI).sum())
    tot_vf, tot_vol = sum(carb_vf.values()), sum(carb_vol.values())
    inj_mol = inj * CO2_MOLAL if inj else None
    r = dict(run_id=os.path.basename(run_dir.rstrip("/")), n_cells=n,
             fracture_volume_m3=float(vol.sum()), pore_volume_m3=pore,
             injected_kg=inj, injected_CO2_mol=inj_mol,
             pore_volumes_flushed=(inj/988.0/pore if inj else None),
             carb_per_cell=tot_vf/n, carb_intensity=tot_vol/pore,
             carb_CO2_mol=carb_mol if mv else None,
             efficiency_pct=(100.0*carb_mol/inj_mol if (mv and inj_mol) else None),
             pH_mean=(float(np.nanmean(ph)) if ph is not None else None),
             pH_volwt=(float(np.nansum(ph*vol)/vol.sum()) if ph is not None else None))
    for m in CARB:
        r[f"pct_{m.lower()}"] = 100*carb_vol[m]/tot_vol if tot_vol > 0 else 0.0
    for m in PRIMARY:
        r[f"diss_{m.lower()}_pct"] = diss[m]
    return r, None

def level_of(run_id):
    m = re.search(r"p32_(\d{3})", run_id)
    return int(m.group(1))/100 if m else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_gravityoff")
    ap.add_argument("--prefix", default="")
    ap.add_argument("--out", default="finalised.csv")
    a = ap.parse_args()
    dirs = [d for d in sorted(glob.glob(os.path.join(a.runs, a.prefix+"*")))
            if os.path.isdir(d)]
    rows, skipped = [], []
    for d in dirs:
        r, why = analyse(d)
        rows.append(r) if r else skipped.append((os.path.basename(d), why))
    if skipped:
        print(f"skipped {len(skipped)}:")
        for nm, why in skipped[:8]: print(f"  {nm}: {why}")
    if not rows: sys.exit("nothing analysed")
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {a.out} ({len(rows)} runs)\n")

    print("="*74); print("SANITY CHECKS"); print("="*74)
    for label, key, msg in (("molar volumes not parsed","carb_CO2_mol","efficiency unavailable"),
                            ("no mass balance file","injected_kg",""),
                            ("pH field absent","pH_mean","add pH to the OUTPUT block")):
        bad = [r["run_id"] for r in rows if r[key] is None]
        if bad: print(f"  {label} for {len(bad)} run(s) {msg}, e.g. {bad[0]}")
    pv = [r["pore_volumes_flushed"] for r in rows if r["pore_volumes_flushed"]]
    if pv:
        print(f"  pore volumes flushed: {min(pv):.3g} to {max(pv):.3g}")
        print(f"    (>10^3 = strongly flushed; >10^6 = rate still wrong)")
    eff = [r["efficiency_pct"] for r in rows if r["efficiency_pct"] is not None]
    if eff:
        over = [r["run_id"] for r in rows
                if r["efficiency_pct"] is not None and r["efficiency_pct"] > 100]
        print(f"  efficiency: {min(eff):.4g}% to {max(eff):.4g}%"
              + (f"  ** {len(over)} EXCEED 100% -- check **" if over else ""))

    grp = defaultdict(list)
    for r in rows:
        lv = level_of(r["run_id"])
        if lv is not None: grp[lv].append(r)
    def agg(lv, key):
        v = [r[key] for r in grp[lv] if r[key] is not None]
        return np.array(v) if v else None

    print("\n"+"="*74); print("TABLE 5  --  outcomes at t = 50 yr, by P32 level"); print("="*74)
    print(f"{'P32':>6}{'n':>4}{'pH':>15}{'carb/cell':>13}{'CV %':>7}{'eff %':>10}{'zeros':>7}")
    for lv in sorted(grp):
        ph, pc, ef = agg(lv,"pH_volwt"), agg(lv,"carb_per_cell"), agg(lv,"efficiency_pct")
        if ph is None: phs = "n/a"
        elif ph.size > 1: phs = f"{ph.mean():.2f} +/- {ph.std(ddof=1):.2f}"
        else: phs = f"{ph.mean():.2f}"
        cv = (100*pc.std(ddof=1)/pc.mean()
              if pc is not None and pc.size > 1 and pc.mean() > 0 else float("nan"))
        efs = f"{ef.mean():.4g}" if ef is not None else "n/a"
        print(f"{lv:>6.2f}{len(grp[lv]):>4}{phs:>15}{pc.mean():>13.3e}"
              f"{cv:>7.0f}{efs:>10}{int((pc==0).sum()):>7}")

    print("\n"+"="*74); print("PRIMARY MINERAL DISSOLUTION  --  % of initial volume consumed"); print("="*74)
    print(f"{'P32':>6}" + "".join(f"{m[:9]:>12}" for m in PRIMARY))
    for lv in sorted(grp):
        line = f"{lv:>6.2f}"
        for m in PRIMARY:
            v = agg(lv, f"diss_{m.lower()}_pct")
            if v is None: line += f"{'n/a':>12}"
            elif v.size > 1: line += f"{v.mean():>7.1f}+-{v.std(ddof=1):<4.1f}"
            else: line += f"{v.mean():>12.1f}"
        print(line)

    print("\n"+"="*74); print("TABLE 6  --  carbonate assemblage, volume-weighted by level"); print("="*74)
    print(f"{'P32':>6}" + "".join(f"{m[:9]:>11}" for m in CARB))
    tot = defaultdict(float)
    for lv in sorted(grp):
        # weight by carbonate VOLUME, not by intensity: intensity is volume
        # divided by pore volume, so weighting by it biases the assemblage
        # toward realisations with small pore volume.
        acc = {m: sum(r[f"pct_{m.lower()}"] * r["carb_intensity"] * r["pore_volume_m3"]
                      for r in grp[lv]) for m in CARB}
        t = sum(acc.values())
        for m in CARB: tot[m] += acc[m]
        print(f"{lv:>6.2f}" + "".join(f"{(100*acc[m]/t if t else 0):>10.0f}%" for m in CARB))
    gt = sum(tot.values())
    print(f"{'ens':>6}" + "".join(f"{(100*tot[m]/gt if gt else 0):>10.0f}%" for m in CARB))

if __name__ == "__main__":
    main()
