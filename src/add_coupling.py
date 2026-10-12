#!/usr/bin/env python3
"""
Enable porosity-permeability coupling in a PFLOTRAN deck: the five lines the
coupled ensemble (runs_aperture) carries, plus a porosity floor:

  CHEMISTRY block          UPDATE_POROSITY
                           MINIMUM_POROSITY 1.d-2   (floor = the critical porosity;
                                                     without it a clogged cell reaches
                                                     zero porosity and the solver breaks down)
                           UPDATE_PERMEABILITY
  MATERIAL_PROPERTY block  PERMEABILITY_POWER 3.d0
                           PERMEABILITY_CRITICAL_POROSITY 0.01d0
                           PERMEABILITY_MIN_SCALE_FACTOR 1.d-6

The chemistry lines go after the DATABASE line, the permeability lines after
the '/' that closes the PERMEABILITY sub-block holding PERM_ISO: the same
blocks as in the runs_aperture decks (position within a block does not matter
to PFLOTRAN).

    python3 src/add_coupling.py runs_gravityoff/C_feedback__*/pflotran_co2.in
A deck that is already coupled but has no floor gets only the floor.
Safe to run twice. A deck in which either place cannot be found exactly once
is left unwritten and reported.
"""
import re, sys

FLOOR = "MINIMUM_POROSITY 1.d-2"
CHEM = ["UPDATE_POROSITY", FLOOR, "UPDATE_PERMEABILITY"]
PERM = ["PERMEABILITY_POWER 3.d0", "PERMEABILITY_CRITICAL_POROSITY 0.01d0",
        "PERMEABILITY_MIN_SCALE_FACTOR 1.d-6"]


def indent(line):
    return line[:len(line) - len(line.lstrip())]


def edit(path):
    lines = open(path).read().split("\n")
    up = [k for k, l in enumerate(lines) if re.match(r"\s*UPDATE_POROSITY\b", l, re.I)]
    if up:
        if any(re.match(r"\s*MINIMUM_POROSITY\b", l, re.I) for l in lines):
            return "already coupled, with floor; no change"
        lines[up[0] + 1:up[0] + 1] = [indent(lines[up[0]]) + FLOOR]
        open(path, "w").write("\n".join(lines))
        return "already coupled; porosity floor added"
    chem = [k for k, l in enumerate(lines) if re.match(r"\s*CHEMISTRY\s*$", l, re.I)]
    db = [k for k, l in enumerate(lines) if re.match(r"\s*DATABASE\b", l, re.I)]
    iso = [k for k, l in enumerate(lines) if re.match(r"\s*PERM_ISO\b", l, re.I)]
    if len(chem) != 1 or len(db) != 1 or db[0] < chem[0]:
        return f"NOT WRITTEN: need one CHEMISTRY block with one DATABASE line (found {len(chem)}, {len(db)})"
    if len(iso) != 1:
        return f"NOT WRITTEN: need one PERM_ISO line (found {len(iso)})"
    close = next((k for k in range(iso[0] + 1, min(iso[0] + 4, len(lines))) if lines[k].strip() == "/"), None)
    if close is None:
        return "NOT WRITTEN: no '/' closing the PERMEABILITY sub-block after PERM_ISO"
    ci = indent(lines[db[0]])
    pi = indent(lines[iso[0] - 1]) if re.match(r"\s*PERMEABILITY\s*$", lines[iso[0] - 1], re.I) else "  "
    # insert the later block first, so the earlier index stays valid
    for at, block, ind in sorted([(db[0] + 1, CHEM, ci), (close + 1, PERM, pi)], reverse=True):
        lines[at:at] = [ind + b for b in block]
    open(path, "w").write("\n".join(lines))
    return "coupling added (6 lines, with the porosity floor)"


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: add_coupling.py <deck.in> [...]")
    bad = 0
    for p in sys.argv[1:]:
        r = edit(p)
        bad += r.startswith("NOT WRITTEN")
        print(f"  {p}: {r}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
