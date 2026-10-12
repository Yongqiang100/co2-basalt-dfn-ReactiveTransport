#!/usr/bin/env python3
"""
Verification runs for the dissolved carbon at shut-in (three shut-in networks, 0 to 10 years).

    python3 src/prepare_verify_shutin.py                     # default networks below
    python3 src/prepare_verify_shutin.py --networks p32_075_s1597 p32_200_s1289 p32_100_s1063

For each network, copies runs_gravityoff/E_feedback__<network> to runs_verify_shutin/ and edits
the deck:
  - CHEMISTRY OUTPUT: adds PRIMARY_SPECIES, so that TOTAL writes the total concentration of each
    primary species (Total CO2(aq) covers all dissolved carbon)
  - TIME: FINAL_TIME 10 y; SNAPSHOT_FILE times up to and including 10 y
  - removes a RESTART line, so the run starts at t = 0
Simulation outputs of the original run (h5, mass balance, logs, checkpoints) are not copied.
The original runs are not modified.
"""
import argparse, os, re, shutil, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
SKIP = re.compile(r"(\.h5$|-mas.*\.dat$|\.out$|\.log$|checkpoint|restart|\.chk|slurm-)", re.I)


def edit(deck):
    s = deck
    # chemistry output: PRIMARY_SPECIES inside the CHEMISTRY ... OUTPUT block
    a = s.index("\nCHEMISTRY"); b = s.index("\nEND", a)
    chem = s[a:b]
    if "PRIMARY_SPECIES\n" not in chem.split("OUTPUT", 1)[-1]:
        chem = re.sub(r"(\n\s*OUTPUT\s*\n)", r"\1    PRIMARY_SPECIES\n", chem, count=1)
    s = s[:a] + chem + s[b:]
    # final time and snapshot times
    s, n = re.subn(r"FINAL_TIME\s+\S+\s+y", "FINAL_TIME 10.0d0 y", s)
    if n != 1:
        raise ValueError("FINAL_TIME line not found once")
    def times(m):
        vals = [float(v) for v in m.group(2).split() if float(v) <= 10.0 + 1e-9]
        if 10.0 not in vals:
            vals.append(10.0)
        return m.group(1) + " ".join(f"{v:g}" for v in vals)
    s, n = re.subn(r"(TIMES\s+y\s+)([0-9.eEdD+\s-]+?)(?=\n)", times, s, count=1)
    if n != 1:
        raise ValueError("SNAPSHOT_FILE TIMES line not found")
    s = re.sub(r"\n\s*RESTART[^\n]*", "", s)
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--networks", nargs="+", default=["p32_075_s1597", "p32_200_s1289", "p32_100_s1063"])
    ap.add_argument("--root", default="runs_gravityoff"); ap.add_argument("--out", default="runs_verify_shutin")
    a = ap.parse_args()
    for net in a.networks:
        src = os.path.join(REV, a.root, f"E_feedback__{net}"); dst = os.path.join(REV, a.out, f"E_feedback__{net}")
        if not os.path.isdir(src):
            sys.exit(f"missing {src}")
        if os.path.exists(dst):
            sys.exit(f"{dst} exists; remove it first")
        os.makedirs(dst)
        for f in os.listdir(src):
            if SKIP.search(f) or os.path.isdir(os.path.join(src, f)):
                continue
            shutil.copy2(os.path.join(src, f), dst)
        p = os.path.join(dst, "pflotran_co2.in")
        with open(p) as fh:
            deck = fh.read()
        with open(p, "w") as fh:
            fh.write(edit(deck))
        print(f"  prepared {os.path.relpath(dst, REV)}")
    print("Run each deck with the usual PFLOTRAN launcher, then: python3 src/verify_shutin_dic.py")


if __name__ == "__main__":
    main()
