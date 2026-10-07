#!/usr/bin/env python3
"""
Prepare a PFLOTRAN deck for option B, immediately before PFLOTRAN starts.

1. Gravity off. The flow is solved in Richards mode, where fluid density does
   not depend on the dissolved species, so gravity adds only a hydrostatic
   offset and drives no flow of its own. With gravity on, the uniform-pressure
   outflow boundary is out of equilibrium and drives circulation through that
   face; with gravity off it is in equilibrium and the injection alone drives
   the flow.
2. Databases. Every DATABASE the deck names must be reachable from the run
   directory. A missing file, or a link that does not resolve, is copied from
   the project's own database first: <co2-basalt>/database next to revision/,
   then ~/co2-basalt/database, and only then beside $PFLOTRAN_DB (on Setonix
   that is the stock pflotran-v6 database, which may differ from the project's).

    python3 src/fix_gravity_off.py pflotran_co2.in          # both steps
    python3 src/fix_gravity_off.py pflotran_co2.in --no-db  # gravity only

Safe to run twice. Exits non-zero, before writing anything, unless the deck
runs in Richards mode, has exactly one GRAVITY line, and every database it
names can be found.
"""
import os, re, shutil, sys

GRAV = re.compile(r"^([ \t]*)GRAVITY[ \t]+(\S+)[ \t]+(\S+)[ \t]+(\S+)[ \t]*$", re.M | re.I)
DB = re.compile(r"^[ \t]*DATABASE[ \t]+(\S+)", re.M | re.I)
_HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(_HERE) if os.path.basename(_HERE) == "src" else _HERE


def num(t):
    return float(t.replace("d", "e").replace("D", "E"))


def db_candidates(name):
    """Reference copies of a database, the project's own first."""
    cands = [os.path.join(os.path.dirname(REV), "database", name),
             os.path.join(os.path.expanduser("~"), "co2-basalt", "database", name)]
    env = os.environ.get("PFLOTRAN_DB", "")
    if env:
        cands += [env] if os.path.basename(env) == name else []
        cands.append(os.path.join(os.path.dirname(env), name))
    return cands


def db_sources(name):
    return next((c for c in db_candidates(name) if os.path.isfile(c)), None)


def plan_databases(deck_dir, text):
    """[(target, source or None)] for every database that is not reachable."""
    todo = []
    for raw in DB.findall(text):
        name = raw.strip("'\"")
        if os.path.isabs(name):
            if not os.path.isfile(name):
                sys.exit(f"FATAL: database {name} not found. Deck NOT written")
            continue
        target = os.path.join(deck_dir, name)
        if not os.path.isfile(target):
            src = db_sources(name)
            if src is None:
                sys.exit(f"FATAL: database {name} is named in the deck but found neither in the run "
                         f"directory nor in any reference directory. Deck NOT written")
            todo.append((target, src))
    return todo


def main():
    args = [a for a in sys.argv[1:] if a != "--no-db"]
    if not args:
        sys.exit("usage: fix_gravity_off.py <deck.in> [--no-db]")
    path = args[0]
    s = open(path).read()
    if not re.search(r"^[ \t]*MODE[ \t]+RICHARDS\b", s, re.M | re.I):
        sys.exit("FATAL: flow mode is not RICHARDS; gravity off is only equivalent there. Deck NOT written")
    g = list(GRAV.finditer(s))
    if len(g) != 1:
        sys.exit(f"FATAL: expected one GRAVITY line, found {len(g)} "
                 f"(none means PFLOTRAN's default, which is on). Deck NOT written")
    dbs = [] if "--no-db" in sys.argv else plan_databases(os.path.dirname(os.path.abspath(path)), s)

    m = g[0]
    if all(num(v) == 0 for v in m.groups()[1:]):
        print("  gravity  : already off; no change")
    else:
        open(path, "w").write(s[:m.start()] + f"{m.group(1)}GRAVITY 0.d0 0.d0 0.d0" + s[m.end():])
        print(f"  gravity  : switched off (was {' '.join(m.groups()[1:])})")
    for target, src in dbs:
        if os.path.islink(target) or os.path.exists(target):
            os.remove(target)                       # a link that does not resolve
        shutil.copy(src, target)
        print(f"  database : copied {os.path.basename(target)} from {os.path.dirname(src)}")


if __name__ == "__main__":
    main()
