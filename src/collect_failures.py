#!/usr/bin/env python3
"""
Collect every failed task across all sets and report how each one stopped.

    python3 src/collect_failures.py              # report, and write lists/failures.tsv
    python3 src/collect_failures.py --archive    # also copy each failed run's outputs to runs_failed/

Sets and job IDs come from lists/jobs.txt (the table status.sh uses). Each
failed task's run directory is read from its own job log (rt-<job>_<i>.out:
" dir     :", rtD-<job>_<i>.out: " run     :"), so no job-to-list mapping is
needed. Per failure: set, run, Slurm state, simulated years reached, PFLOTRAN's
stop reason, coupled / floored, and whether a later rerun has since reached
50 years (then the failure is superseded, and nothing needs doing).
Archiving copies only outputs (h5, logs); decks and meshes are left in place.
"""
import glob, os, re, shutil, subprocess, sys

_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE
FAILED = ("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL")      # plus COMPLETED runs short of 50 y


def sacct(ids):
    # standard = full ISO dates; some sites default to a relative format that drops the date for today
    r = subprocess.run(["sacct", "-X", "-n", "-P", "-j", ids, "-o", "JobID,State,Elapsed,End"],
                       capture_output=True, text=True, env=dict(os.environ, SLURM_TIME_FORMAT="standard"))
    for line in r.stdout.splitlines():
        t = line.split("|")
        if len(t) == 4 and "_" in t[0] and "[" not in t[0]:
            yield t[0], t[1].split()[0], t[2], t[3]


def floor_added():
    """{run id: time the floor was added}, from lists/minpor_applied.txt."""
    p, out = os.path.join(REV, "lists", "minpor_applied.txt"), {}
    if os.path.isfile(p):
        for line in open(p):
            t = line.rstrip("\n").split("\t")
            if len(t) >= 2:
                out[t[0]] = t[1].replace(" ", "T")
    return out


def run_dir_from_log(task):
    job, idx = task.split("_")
    for pat, key in ((f"rtD-{job}_{idx}.out", "run     :"), (f"rt-{job}_{idx}.out", "dir     :")):
        p = os.path.join(REV, pat)
        if os.path.isfile(p):
            for line in open(p, errors="replace"):
                if key in line:
                    v = line.split(":", 1)[1].strip()
                    if os.path.isabs(v):
                        return os.path.normpath(v), p
                    return (os.path.join(REV, v) if "/" in v else os.path.join(REV, "runs", v)), p
            return None, p
    return None, None


def last_time(logfile):
    t = None
    if os.path.isfile(logfile):
        for line in open(logfile, errors="replace"):
            m = re.search(r"\bStep\s+\d+\s+Time=\s*([0-9.Ee+-]+)", line)
            if m:
                t = float(m.group(1))
    return t


def reason(logfile):
    r = ""
    if os.path.isfile(logfile):
        for line in open(logfile, errors="replace"):
            if re.search(r"Stopping:|ERROR", line):
                r = line.strip()
    return r[:70]


H5_ERROR = []


def reached_50(d):
    """True if the run's output reaches 50 y: from the HDF5 file, or, if that cannot be
    read, from the last time step in its PFLOTRAN log (the error is reported once)."""
    h5 = os.path.join(d, "pflotran_co2.h5")
    if os.path.isfile(h5):
        try:
            import h5py
            with h5py.File(h5, "r") as f:
                return max((float(k.split("Time")[1].split()[0]) for k in f if "Time" in k), default=0) >= 49.99
        except Exception as e:
            if not H5_ERROR:
                H5_ERROR.append(f"{type(e).__name__}: {e}")
                print(f"  note: cannot read {h5} ({H5_ERROR[0][:120]}); using the PFLOTRAN log instead", file=sys.stderr)
    t = last_time(os.path.join(d, "pflotran.log"))
    return t is not None and t >= 49.99


def main():
    archive = "--archive" in sys.argv
    table = os.path.join(REV, "lists", "jobs.txt")
    if not os.path.isfile(table):
        sys.exit(f"no job table: {table}")
    rows, added = [], floor_added()
    for line in open(table):
        t = line.split()
        if len(t) < 2 or t[0].startswith("#"):
            continue
        name, ids = t[0], t[1]
        for task, state, elapsed, end in sacct(ids):
            if state not in FAILED and state != "COMPLETED":
                continue
            d, log = run_dir_from_log(task)
            if state == "COMPLETED":             # exited cleanly: fine, unless it stopped short of 50 y
                if not d or reached_50(d):
                    continue
                state = "STOPPED_EARLY"
            rid = os.path.basename(d) if d else "?"
            deck = next(iter(sorted(glob.glob(os.path.join(d, "*.in")))), None) if d else None
            txt = open(deck).read() if deck else ""
            coupled = bool(re.search(r"^\s*UPDATE_POROSITY", txt, re.M | re.I))
            floored = bool(re.search(r"^\s*MINIMUM_POROSITY", txt, re.M | re.I))
            if floored and rid in added and end[:2] == "20" and added[rid] > end:
                floored = False          # the floor went in after this task ended: it ran without it
            plog = os.path.join(d, "pflotran.log") if d else ""
            done = bool(d) and reached_50(d)
            why = reason(plog) if d else ""
            if not why:
                why = {"TIMEOUT": "(Slurm time limit reached)", "OUT_OF_MEMORY": "(Slurm: out of memory)",
                       "NODE_FAIL": "(Slurm: node failure)",
                       "STOPPED_EARLY": "(clean exit before 50 y: wall-clock stop? continue with make_restart.py)"
                       }.get(state, "(no stop message in pflotran.log)")
            rows.append(dict(set=name, task=task, run=rid, state=state, elapsed=elapsed,
                             years=last_time(plog) if d else None, why=why,
                             config=("coupled+floor" if floored else "coupled") if coupled else "fixed",
                             now="rerun reached 50 y" if done else "", dir=d))
    if not rows:
        print("no failed tasks"); return
    W = max(24, max(len(r["run"]) for r in rows) + 2)
    print(f"{'set':<14}{'run':<{W}}{'state':<15}{'config':<15}{'years':>7}  stop reason / status")
    for r in rows:
        y = f"{r['years']:.3g}" if r["years"] is not None else "-"
        print(f"{r['set']:<14}{r['run']:<{W}}{r['state']:<15}{r['config']:<15}{y:>7}  {r['now'] or r['why']}")
    open_ = [r for r in rows if not r["now"]]
    print(f"\n{len(rows)} failed tasks; {len(rows) - len(open_)} since completed by a rerun; {len(open_)} outstanding")
    by = {}
    for r in open_:
        by.setdefault((r["set"], r["config"]), []).append(r["run"])
    for (s, c), runs in sorted(by.items()):
        print(f"  {s} ({c}): {len(runs)}  {', '.join(runs)}")
    out = os.path.join(REV, "lists", "failures.tsv")
    with open(out, "w") as f:
        f.write("set\ttask\trun\tstate\telapsed\tconfig\tyears_reached\tstop_reason\tstatus\n")
        for r in rows:
            f.write("\t".join(str(r[k] if r[k] is not None else "") for k in
                              ("set", "task", "run", "state", "elapsed", "config", "years", "why", "now")) + "\n")
    print(f"written: {os.path.relpath(out, REV)}")
    if archive:
        n = 0
        for r in open_:
            if not r["dir"] or not os.path.isdir(r["dir"]):
                continue
            dst = os.path.join(REV, "runs_failed", r["run"])
            if os.path.exists(dst):
                continue
            os.makedirs(dst)
            for f in glob.glob(os.path.join(r["dir"], "pflotran_co2*.h5")) + glob.glob(os.path.join(r["dir"], "*.log")) \
                     + glob.glob(os.path.join(r["dir"], "pflotran_co2.out")) + glob.glob(os.path.join(r["dir"], "*.in")):
                shutil.copy2(f, dst)
            n += 1
        print(f"archived outputs of {n} outstanding failures to runs_failed/ (existing copies kept)")


if __name__ == "__main__":
    main()
