#!/usr/bin/env python3
"""
Add checkpointing and a graceful wall-clock stop to a deck, so a run that would
reach the queue's time limit stops cleanly, writes a restart file, and can be
continued with make_restart.py (PFLOTRAN v6 CHECKPOINT and WALLCLOCK_STOP cards).

  SIMULATION block:   CHECKPOINT
                        PERIODIC TIME <every> y     (safety copies, e.g. after a node failure)
                        FORMAT HDF5
                      /
  SUBSURFACE block:   WALLCLOCK_STOP <stop> h       (stop before the job's time limit)

    python3 src/add_checkpoint.py runs/D40_feedback__*/pflotran_co2.in --stop 23 --every 5
Defaults: --stop 23 (hours; for a 24 h job), --every 5 (years). Safe to run twice.
A deck whose SIMULATION or SUBSURFACE block cannot be found once is left unwritten.
"""
import re, sys


def edit(path, stop, every):
    lines = open(path).read().split("\n")
    sim = [k for k, l in enumerate(lines) if re.fullmatch(r"\s*SIMULATION\s*", l, re.I)]
    sub = [k for k, l in enumerate(lines) if re.fullmatch(r"\s*SUBSURFACE\s*", l, re.I)]
    if len(sim) != 1 or len(sub) != 1:
        return f"NOT WRITTEN: need one SIMULATION and one SUBSURFACE line (found {len(sim)}, {len(sub)})"
    end = next((k for k in range(sim[0] + 1, len(lines)) if re.fullmatch(r"\s*END\s*", lines[k], re.I)), None)
    if end is None or end > sub[0]:
        return "NOT WRITTEN: no END closing the SIMULATION block before SUBSURFACE"
    done = []
    has_ckpt = any(re.match(r"\s*CHECKPOINT\b", lines[k], re.I) for k in range(sim[0], end))
    has_stop = any(re.match(r"\s*WALLCLOCK_STOP\b", l, re.I) for l in lines)
    # later insertion first, so the earlier index stays valid
    if not has_stop:
        lines[sub[0] + 1:sub[0] + 1] = [f"  WALLCLOCK_STOP {stop:g} h"]; done.append(f"WALLCLOCK_STOP {stop:g} h")
    if not has_ckpt:
        lines[end:end] = ["  CHECKPOINT", f"    PERIODIC TIME {every:g}.d0 y", "    FORMAT HDF5", "  /"]
        done.append(f"CHECKPOINT every {every:g} y (HDF5)")
    if not done:
        return "already has CHECKPOINT and WALLCLOCK_STOP; no change"
    open(path, "w").write("\n".join(lines))
    return "added " + " and ".join(done)


def main():
    a = sys.argv[1:]
    stop = float(a[a.index("--stop") + 1]) if "--stop" in a else 23.0
    every = float(a[a.index("--every") + 1]) if "--every" in a else 5.0
    decks = [x for x in a if x.endswith(".in")]
    if not decks:
        sys.exit(__doc__)
    bad = 0
    for p in decks:
        r = edit(p, stop, every); bad += r.startswith("NOT WRITTEN")
        print(f"  {p}: {r}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
