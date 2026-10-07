#!/usr/bin/env python3
"""
Let mineral surface areas evolve (PFLOTRAN v6 SURFACE_AREA_FUNCTION), for the
surface-area sensitivity test.

  primary minerals  SURFACE_AREA_FUNCTION VOLUME_FRACTION_RATIO,  SURFACE_AREA_VOL_FRAC_POWER p
                    (shrinking grains: a = a0 (phi_m / phi_m0)^p)
  secondary         SURFACE_AREA_FUNCTION POROSITY_RATIO,  SURFACE_AREA_POROSITY_POWER s
                    (a = a0 (phi / phi0)^s: less surface as the pores fill)

No chemistry-block switch is needed: PFLOTRAN 6.0 rejects the old
UPDATE_MINERAL_SURFACE_AREA ("updated to SURFACE_AREA_FUNCTION within each
mineral block"), and this script removes it from any deck that still has it.

The volume-fraction form is for primary minerals only: secondary minerals start
at the 1e-6 seed, so their ratio would grow to ~1e5 and so would their surface.

    python3 src/add_surface_area.py runs/S1_feedback__*/pflotran_co2.in [--prim 0.666666667] [--sec 0.666666667]
Safe to run twice. A deck whose MINERAL_KINETICS block or mineral entries cannot
be found is left unwritten and reported.
"""
import re, sys

PRIM = ("Anorthite", "Albite", "Diopside", "Forsterite", "Fayalite", "Enstatite")
SEC = ("Calcite", "Magnesite", "Siderite", "Dawsonite", "Kaolinite", "Chalcedony")


def indent(l):
    return l[:len(l) - len(l.lstrip())]


def edit(path, p, s):
    lines = open(path).read().split("\n")
    n0 = len(lines)
    lines = [l for l in lines if not re.fullmatch(r"\s*UPDATE_MINERAL_SURFACE_AREA\s*", l, re.I)]
    removed = n0 - len(lines)
    if any(re.match(r"\s*SURFACE_AREA_FUNCTION\b", l, re.I) for l in lines):
        if removed:
            open(path, "w").write("\n".join(lines))
            return "already has surface-area functions; removed the obsolete UPDATE_MINERAL_SURFACE_AREA"
        return "already has surface-area functions; no change"
    mk = [k for k, l in enumerate(lines) if re.fullmatch(r"\s*MINERAL_KINETICS\s*", l, re.I)]
    if len(mk) != 1:
        return f"NOT WRITTEN: need one MINERAL_KINETICS block (found {len(mk)})"
    base = indent(lines[mk[0]])
    end = next((k for k in range(mk[0] + 1, len(lines))
                if lines[k].strip() == "/" and indent(lines[k]) == base), None)
    if end is None:
        return "NOT WRITTEN: no '/' closing MINERAL_KINETICS"
    heads = {}
    for k in range(mk[0] + 1, end):
        name = lines[k].strip()
        if name in PRIM + SEC and name not in heads:
            heads[name] = k
    missing = [m for m in PRIM + SEC if m not in heads]
    if missing:
        return f"NOT WRITTEN: minerals not found in MINERAL_KINETICS: {', '.join(missing)}"
    ins = []
    for m, k in heads.items():
        ind = indent(lines[k]) + "  "
        if m in PRIM:
            ins.append((k + 1, [f"{ind}SURFACE_AREA_FUNCTION VOLUME_FRACTION_RATIO", f"{ind}SURFACE_AREA_VOL_FRAC_POWER {p:.9g}"]))
        else:
            ins.append((k + 1, [f"{ind}SURFACE_AREA_FUNCTION POROSITY_RATIO", f"{ind}SURFACE_AREA_POROSITY_POWER {s:.9g}"]))
    for at, block in sorted(ins, reverse=True):
        lines[at:at] = block
    open(path, "w").write("\n".join(lines))
    return f"surface areas evolve: {len(PRIM)} primary (VF power {p:g}), {len(SEC)} secondary (porosity power {s:g})"


def main():
    a = sys.argv[1:]
    p = float(a[a.index("--prim") + 1]) if "--prim" in a else 2 / 3
    s = float(a[a.index("--sec") + 1]) if "--sec" in a else 2 / 3
    decks = [x for x in a if x.endswith(".in")]
    if not decks:
        sys.exit(__doc__)
    bad = 0
    for d in decks:
        r = edit(d, p, s); bad += r.startswith("NOT WRITTEN")
        print(f"  {d}: {r}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
