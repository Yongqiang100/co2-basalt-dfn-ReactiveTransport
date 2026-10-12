#!/usr/bin/env python3
"""
Generate the ADDITIONAL DFN realisations for the WRR revision.
Wraps prepare_dfn.py without modifying its 5x5 matrix logic, so the
original 25 realisations stay reproducible from the archived code.

Requires the one-line make_config patch in 01_matrix_ext.patch
(only needed for --block D; A and B work unpatched).

Usage
-----
  python 01b_extend_matrix.py --block A --dry     # +10 seeds x 5 levels  = 50
  python 01b_extend_matrix.py --block B --dry     # 10 seeds x x0.90,x1.75 = 20
  python 01b_extend_matrix.py --block D --dry     # 8 seeds x L=30,40 m    = 16
  python 01b_extend_matrix.py --block A --index 7 # generate ONE case (Slurm array)

Drop --dry to generate. Use --index for one case per array task.
"""
import argparse, os, sys, json

# prepare_dfn.py lives in the deposit root and computes OUTPUT_ROOT from cwd at
# IMPORT time, so both sys.path and cwd must point there before importing.
_HERE = os.path.dirname(os.path.abspath(__file__))
_DEPOSIT = os.environ.get("DEPOSIT_ROOT", os.path.abspath(os.path.join(_HERE, "..", "..")))
if not os.path.exists(os.path.join(_DEPOSIT, "prepare_dfn.py")):
    sys.exit(f"FATAL: prepare_dfn.py not found in {_DEPOSIT}\n"
             f"       set DEPOSIT_ROOT to the co2-basalt clone root")
os.chdir(_DEPOSIT)
sys.path.insert(0, _DEPOSIT)
import prepare_dfn as P

# 10 new seeds, disjoint from the published set [42,117,259,383,501]
NEW_SEEDS = [613, 727, 839, 941, 1063, 1181, 1289, 1397, 1481, 1597]

BLOCKS = {
    # (label, p32_mult, seed, domain_L)  -- domain None => BASE_CONFIG default
    "A": [(f"p32_{int(m*100):03d}_s{s}", m, s, None)
          for m in [0.75, 1.00, 1.25, 1.50, 2.00] for s in NEW_SEEDS],
    "B": [(f"p32_{int(m*100):03d}_s{s}", m, s, None)
          for m in [0.90, 1.75] for s in NEW_SEEDS],
    "D": [(f"L{L:02d}_p32_100_s{s}", 1.00, s, L)
          for L in [30, 40] for s in NEW_SEEDS[:8]],
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block", required=True, choices=sorted(BLOCKS))
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--index", type=int, default=None,
                    help="0-based: generate only this case (for Slurm arrays)")
    a = ap.parse_args()

    cases = BLOCKS[a.block]
    if a.index is not None:
        if not 0 <= a.index < len(cases):
            sys.exit(f"index {a.index} out of range 0..{len(cases)-1}")
        cases = [cases[a.index]]

    print(f"Block {a.block}: {len(cases)} case(s)")
    for name, m, s, L in cases:
        dom = [float(L)]*3 if L else P.BASE_CONFIG["domain"]
        est = P.BASE_CONFIG["h"]
        print(f"  {name:<24s} P32x{m:<5.2f} seed={s:<5d} domain={dom} h={est}")

    if a.dry:
        # Warn loudly about the cost of the large-domain block
        if a.block == "D":
            base = 81_000
            for L in [30, 40]:
                print(f"  NOTE L={L} m at h=0.2 -> ~{base*(L/20)**3/1e3:.0f}k cells/realisation")
        print("\n--dry set; nothing generated.")
        return

    os.makedirs(P.OUTPUT_ROOT, exist_ok=True)
    manifest = []
    for name, m, s, L in cases:
        print(f"\n=== {name} ===", flush=True)
        # make_config() returns a plain dict, so override the domain on the
        # RESULT rather than patching prepare_dfn.py. Keeps the deposit
        # byte-identical to the archive, which the reproducibility claim needs.
        cfg = P.make_config(m)
        if L:
            cfg["domain"] = [float(L)] * 3
        r = P.generate_single_dfn(cfg, name, seed=s)
        manifest.append({"name": name, "p32_mult": m, "seed": s,
                         "domain_L": L, "ok": r is not None})
    mf = f"manifest_block{a.block}" + (f"_{a.index}" if a.index is not None else "") + ".json"
    with open(mf, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nWrote {mf}")

if __name__ == "__main__":
    main()
