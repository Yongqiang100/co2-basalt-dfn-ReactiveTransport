#!/usr/bin/env python3
"""Resolve a run set block + array index to a concrete run specification.

Output (tab-separated): run_id <TAB> dfn_case <TAB> run_dir <TAB> variant <TAB> variant_case_index
variant is "-" for plain runs (Blocks A, B, D).
"""
import os, sys, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
DEPOSIT = os.environ.get("DEPOSIT_ROOT", os.path.abspath(os.path.join(HERE, "..", "..")))
RUNS = os.path.join(HERE, "..", "runs")
NEW_SEEDS = [613, 727, 839, 941, 1063, 1181, 1289, 1397, 1481, 1597]

def block_runs(block):
    out = []
    if block == "A":
        for m in (75, 100, 125, 150, 200):
            for s in NEW_SEEDS:
                c = f"p32_{m:03d}_s{s}"
                out.append((f"A_{c}", c, os.path.join(RUNS, f"A_{c}"), "-", ""))
    elif block == "B":
        for m in (90, 175):
            for s in NEW_SEEDS:
                c = f"p32_{m:03d}_s{s}"
                out.append((f"B_{c}", c, os.path.join(RUNS, f"B_{c}"), "-", ""))
    elif block in ("D30", "D40"):
        L = 30 if block == "D30" else 40
        for s in NEW_SEEDS[:8]:
            c = f"L{L}_p32_100_s{s}"
            out.append((f"{block}_{c}", c, os.path.join(RUNS, f"{block}_{c}"), "-", ""))
    elif block == "C":
        sys.path.insert(0, os.path.join(HERE, "..", "config"))
        import variants as V
        for vname, spec in V.VARIANTS.items():
            for i, case in enumerate(spec["cases"]):
                rid = f"C_{vname}__{case}"
                out.append((rid, case, os.path.join(RUNS, rid), vname, str(i)))
    else:
        sys.exit(f"unknown block {block!r} (A|B|C|D30|D40)")
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block", required=True)
    ap.add_argument("--index", type=int)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--count", action="store_true")
    a = ap.parse_args()
    runs = block_runs(a.block)
    if a.count:
        print(len(runs)); return
    if a.list or a.index is None:
        miss = 0
        for i, r in enumerate(runs):
            mesh = os.path.join(DEPOSIT, "dfn_library", r[1], "full_mesh.uge")
            ok = os.path.isfile(mesh); miss += (not ok)
            print(f"  {i:>3} {r[0]:<34} {r[1]:<20} {'' if ok else 'MESH MISSING'}")
        print(f"\n  {len(runs)} run(s) in block {a.block}"
              + (f", {miss} with NO MESH -- generate DFNs first" if miss else ""))
        print(f"  submit: sbatch --array=0-{len(runs)-1}%6 --export=ALL,BLOCK={a.block} slurm/run_rt.sh")
        return
    if not 0 <= a.index < len(runs):
        sys.exit(f"index {a.index} out of range 0..{len(runs)-1}")
    print("\t".join(runs[a.index]))

if __name__ == "__main__":
    main()
