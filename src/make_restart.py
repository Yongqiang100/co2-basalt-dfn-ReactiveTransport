#!/usr/bin/env python3
"""
Prepare a run that stopped early (wall-clock stop, time limit, node failure) to
continue from its restart or checkpoint file, then resubmit it with run_dirs.sh.

    python3 src/make_restart.py runs/D40_feedback__L40_p32_100_s1181 [...]

For each run directory:
  1. picks the restart file: '*-restart*' if present (written at a graceful
     stop), otherwise the newest periodic checkpoint;
  2. saves this part's outputs under new names, since the continued run
     rewrites them:  pflotran_co2.h5 -> pflotran_co2_part1.h5,
     pflotran.log -> pflotran_part1.log, pflotran_co2-mas.dat -> ..._part1.dat
     (part2, part3 on later continuations);
  3. puts  RESTART / FILENAME <file> /  in the SIMULATION block (replacing an
     earlier RESTART block).
run_dirs.sh sees the RESTART block and keeps the files instead of clearing them.
"""
import glob, os, re, shutil, sys

PREFIX = "pflotran_co2"


def restart_file(d):
    rs = [f for f in glob.glob(os.path.join(d, "*restart*")) if not f.endswith(".in")]
    if rs:
        return max(rs, key=os.path.getmtime)
    ck = glob.glob(os.path.join(d, "*.chk")) + [f for f in glob.glob(os.path.join(d, PREFIX + "-*.h5"))
                                                  if re.search(r"-\d[\d.]*[a-z]+\.h5$", f)]
    return max(ck, key=os.path.getmtime) if ck else None


def prepare(d):
    deck = os.path.join(d, PREFIX + ".in")
    if not os.path.isfile(deck):
        return f"NOT DONE: no {PREFIX}.in"
    rf = restart_file(d)
    if rf is None:
        return "NOT DONE: no restart or checkpoint file (was CHECKPOINT in the deck?)"
    n = 1
    while os.path.exists(os.path.join(d, f"{PREFIX}_part{n}.h5")):
        n += 1
    saved = []
    for src, dst in ((f"{PREFIX}.h5", f"{PREFIX}_part{n}.h5"), ("pflotran.log", f"pflotran_part{n}.log"),
                     (f"{PREFIX}-mas.dat", f"{PREFIX}-mas_part{n}.dat")):
        if os.path.isfile(os.path.join(d, src)):
            shutil.move(os.path.join(d, src), os.path.join(d, dst)); saved.append(dst)
    lines = open(deck).read().split("\n")
    # drop an earlier RESTART block
    out, skip = [], False
    for l in lines:
        if re.fullmatch(r"\s*RESTART\s*", l, re.I):
            skip = True; continue
        if skip:
            if l.strip() == "/":
                skip = False
            continue
        out.append(l)
    sim = next(k for k, l in enumerate(out) if re.fullmatch(r"\s*SIMULATION\s*", l, re.I))
    end = next(k for k in range(sim + 1, len(out)) if re.fullmatch(r"\s*END\s*", out[k], re.I))
    out[end:end] = ["  RESTART", f"    FILENAME {os.path.basename(rf)}", "  /"]
    open(deck, "w").write("\n".join(out))
    return f"continues from {os.path.basename(rf)}; part {n} saved ({', '.join(saved) or 'no outputs found'})"


def main():
    dirs = [x for x in sys.argv[1:] if os.path.isdir(x)]
    if not dirs:
        sys.exit(__doc__)
    bad = 0
    for d in dirs:
        r = prepare(d); bad += r.startswith("NOT DONE")
        print(f"  {d}: {r}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
