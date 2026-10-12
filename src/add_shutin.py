#!/usr/bin/env python3
"""
Stop the injection at time T (years): the steady rate at T, zero just after
(T x (1 + offset)). The steady rate is the rate of the last listed row.

  T after the ramp   the two rows are appended (E decks: 10 y -> 1.000001e+01, offset 1e-6)
  T inside the ramp  rows at or after T are dropped first, then the two rows
                     added: the old F_1d deck (2.74e-3 y -> 2.740274e-3, offset 1e-4)

    python3 src/add_shutin.py [--offset REL] <deck.in> <years>          # edit in place
    python3 src/add_shutin.py --check [--offset REL] <deck.in> <years>  # exit 0 if it stops at T

Offset defaults to 1e-6 (the E decks); the old F decks use 1e-4. Refuses,
writing nothing, if the RATE LIST is not in years or already ends at zero.
"""
import re, sys

NUM = r"[-+]?(\d+\.?\d*|\.\d+)([dDeE][-+]?\d+)?"


def val(t):
    return float(t.replace("d", "e").replace("D", "e"))


def rows_of(lines):
    heads = [k for k, l in enumerate(lines) if re.match(r"\s*RATE\s+LIST\s*$", l, re.I)]
    if len(heads) != 1:
        return None, None, f"expected one RATE LIST, found {len(heads)}"
    units, rows, k = None, [], heads[0] + 1
    while k < len(lines):
        t = lines[k].split()
        if not t:
            k += 1; continue
        if t[0].upper() == "TIME_UNITS":
            units = t[1].lower() if len(t) > 1 else None; k += 1; continue
        if t[0].upper() == "DATA_UNITS":
            k += 1; continue
        if len(t) >= 2 and re.fullmatch(NUM, t[0]) and re.fullmatch(NUM, t[1]):
            rows.append(k); k += 1; continue
        break
    if units not in ("y", "yr", "year", "years"):
        return None, None, f"RATE LIST time units are '{units}', not years"
    if not rows:
        return None, None, "RATE LIST has no data rows"
    return rows, units, ""


def edit(path, T, off):
    lines = open(path).read().split("\n")
    rows, _, err = rows_of(lines)
    if err:
        return "NOT WRITTEN: " + err
    t_last, r_last = lines[rows[-1]].split()[:2]
    if val(r_last) == 0:
        return "NOT WRITTEN: the schedule already ends at zero"
    keep = [k for k in rows if val(lines[k].split()[0]) < T]
    if not keep:
        return f"NOT WRITTEN: no listed row before {T:g} y"
    dropped = [k for k in rows if k not in keep]
    ind = lines[keep[-1]][:len(lines[keep[-1]]) - len(lines[keep[-1]].lstrip())]
    new = [f"{ind}{T:.9e}   {r_last}", f"{ind}{T * (1 + off):.9e}   0.000000e+00"]
    for k in sorted(dropped, reverse=True):
        del lines[k]
    lines[keep[-1] + 1:keep[-1] + 1] = new
    open(path, "w").write("\n".join(lines))
    how = f"; {len(dropped)} ramp row(s) at or after T dropped, as in the old F_1d deck" if dropped else ""
    return f"injection stops at {T:g} y (rate {r_last} until then){how}"


def check(path, T, off):
    lines = open(path).read().split("\n")
    rows, _, err = rows_of(lines)
    if err or len(rows) < 2:
        return False
    (t1, r1), (t2, r2) = (lines[rows[-2]].split()[:2], lines[rows[-1]].split()[:2])
    close = lambda a, b: abs(a - b) <= 1e-9 * max(1.0, abs(b))
    return val(r2) == 0 and val(r1) > 0 and close(val(t1), T) and close(val(t2), T * (1 + off))


def main():
    a = sys.argv[1:]
    chk = a[:1] == ["--check"]
    if chk:
        a = a[1:]
    off = 1e-6
    if a[:1] == ["--offset"]:
        off = float(a[1]); a = a[2:]
    if len(a) != 2:
        sys.exit("usage: add_shutin.py [--check] [--offset REL] <deck.in> <years>")
    if chk:
        sys.exit(0 if check(a[0], float(a[1]), off) else 1)
    r = edit(a[0], float(a[1]), off)
    print(f"  {a[0]}: {r}")
    sys.exit(1 if r.startswith("NOT WRITTEN") else 0)


if __name__ == "__main__":
    main()
