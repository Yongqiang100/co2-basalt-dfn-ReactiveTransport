"""
deckmod — declarative edits to a PFLOTRAN deck.

WHY THIS EXISTS
---------------
run_pflotran.py builds the deck as one ~500-line f-string (line 260) and exposes
only six parameters: temperature, pressure, water_rate_kg_s, co2_molal,
inject_years, total_years. Everything the revision needs to vary — mineral
volume fractions, reactive surface areas, the secondary-phase list, porosity,
UPDATE_POROSITY, dual-continuum — is hard-coded.

Rather than fork or reimplement that validated template, this module post-
processes the deck it emits. Consequences:

  * The baseline deck stays byte-identical to the published one.
  * Every sensitivity variant is a small, reviewable diff.
  * Provenance is exact: the applied ops ARE the description of the variant.

Every op records what it matched. `apply()` raises if an op matches nothing, so
a silently-ineffective sensitivity run is impossible.

Usage
-----
    d = Deck.from_file("carbfix.in")
    d.set_mineral_field("Anorthite", "SPECIFIC_SURFACE_AREA", "30.d0 cm^2/cm^3")
    d.remove_mineral("Dawsonite")
    d.write("carbfix_variant.in")
    print(d.diff())          # eyeball once per variant
    print(d.provenance())    # -> manifest
"""
from __future__ import annotations
import re, difflib, json


class DeckError(RuntimeError):
    pass


class Deck:
    def __init__(self, text: str, name: str = "deck"):
        self.original = text
        self.text = text
        self.name = name
        self.ops: list[dict] = []

    # ---------------------------------------------------------------- io
    @classmethod
    def from_file(cls, path):
        with open(path) as f:
            return cls(f.read(), name=str(path))

    def write(self, path):
        with open(path, "w") as f:
            f.write(self.text)
        return path

    # ------------------------------------------------------- block finding
    @staticmethod
    def _indent(line):
        return len(line) - len(line.lstrip())

    def _block_span(self, header_re, text=None, base=0):
        """Return (start, end) char offsets of a BLOCK ... END/'/' region.

        Nesting is resolved by INDENTATION, not by a keyword whitelist: the
        block ends at the first following line that is END or a bare '/' at an
        indent no deeper than the header. PFLOTRAN decks emitted by
        run_pflotran.py are consistently indented, and this handles arbitrary
        sub-blocks (mineral names, PREFACTOR, PERMEABILITY) without needing to
        enumerate them.
        """
        src = self.text if text is None else text
        m = re.search(header_re, src, re.M)
        if not m:
            raise DeckError(f"[{self.name}] no block matching {header_re!r}")
        # `^\s*` in caller-supplied patterns can swallow preceding blank lines,
        # because \s matches newlines. Anchor to the line holding the KEYWORD,
        # not to m.start(), or parts[0] is a blank line and the indent
        # calculation in insert_in_block is off by one.
        g = m.group(0)
        kw = m.start() + (len(g) - len(g.lstrip()))
        head_start = src.rfind("\n", 0, kw) + 1
        lines = src[head_start:].splitlines(keepends=True)
        h_ind = self._indent(lines[0])
        closer = re.compile(r"^\s*(END|/)\s*$", re.I)
        off = len(lines[0])
        for ln in lines[1:]:
            if closer.match(ln) and self._indent(ln) <= h_ind:
                return base + head_start, base + head_start + off + len(ln)
            off += len(ln)
        raise DeckError(
            f"[{self.name}] unterminated block {header_re!r} "
            f"(header indent {h_ind}; check deck indentation is consistent)")

    def _sub_in_block(self, header_re, pattern, repl, op):
        s, e = self._block_span(header_re)
        block = self.text[s:e]
        new, n = re.subn(pattern, repl, block, flags=re.M)
        if n == 0:
            raise DeckError(
                f"[{self.name}] op {op['op']}: pattern {pattern!r} not found "
                f"in block {header_re!r} — the deck template has changed; "
                f"fix the op rather than ignoring this")
        self.text = self.text[:s] + new + self.text[e:]
        op["matches"] = n
        self.ops.append(op)
        return n

    # -------------------------------------------------------------- ops
    def set_block_value(self, block_header, key, value):
        """MATERIAL_PROPERTY fracture / POROSITY 0.50d0 -> new value."""
        return self._sub_in_block(
            block_header, rf"^(\s*{re.escape(key)}\s+).*$", "\\g<1>" + value.replace("\\", "\\\\"),
            {"op": "set_block_value", "block": block_header,
             "key": key, "value": value})

    def set_mineral_field(self, mineral, field, value):
        """Set a field inside MINERAL_KINETICS / <mineral>.

        Scoped to the named mineral's own sub-block, so setting Anorthite's
        SPECIFIC_SURFACE_AREA cannot silently hit Albite's.
        """
        ks, ke = self._block_span(r"^\s*MINERAL_KINETICS\s*$")
        kin = self.text[ks:ke]
        ms, me = self._block_span(rf"^\s*{re.escape(mineral)}\s*$", text=kin, base=0)
        body, n = re.subn(rf"^(\s*{re.escape(field)}\s+).*$",
                          "\\g<1>" + value.replace("\\", "\\\\"),
                          kin[ms:me], flags=re.M)
        if n == 0:
            raise DeckError(
                f"[{self.name}] field {field!r} not found under {mineral!r}")
        self.text = self.text[:ks] + kin[:ms] + body + kin[me:] + self.text[ke:]
        self.ops.append({"op": "set_mineral_field", "mineral": mineral,
                         "field": field, "value": value, "matches": n})
        return n

    def set_volume_fraction(self, mineral, value):
        """MINERAL_VOLUME_FRACTIONS entry in the initial CONSTRAINT."""
        return self._sub_in_block(
            r"^\s*MINERAL_VOLUME_FRACTIONS\s*$",
            rf"^(\s*{re.escape(mineral)}\s+)\S+(.*)$",
            "\\g<1>" + value.replace("\\", "\\\\") + "\\g<2>",
            {"op": "set_volume_fraction", "mineral": mineral, "value": value})

    def _constraint_minerals_span(self):
        """Locate the CONSTRAINT's MINERALS block, not CHEMISTRY's name list.

        The deck has THREE blocks headed MINERALS (run_pflotran.py:356, 508, 598).
        Only the CONSTRAINT one carries data:
            Anorthite     0.30d0    10.d0  cm^2/cm^3
        i.e. name, volume fraction, specific surface area, units. Discriminate on
        that: the block whose entries have three or more fields.
        """
        for m in re.finditer(r"^[ \t]*MINERALS[ \t]*$", self.text, re.M):
            try:
                s0, e0 = self._block_span(r"^[ \t]*MINERALS[ \t]*$",
                                          text=self.text[m.start():], base=m.start())
            except DeckError:
                continue
            body = self.text[s0:e0].splitlines()[1:-1]
            if any(len(l.split()) >= 3 for l in body if l.strip()):
                return s0, e0
        raise DeckError(f"[{self.name}] no CONSTRAINT MINERALS block "
                        f"(entries with >=3 fields) found")

    def _set_constraint_field(self, mineral, field_no, value, label):
        # A value containing whitespace injects extra tokens into a positional
        # line: 'Anorthite 0.30d0 30.d0 cm^2/cm^3  cm^2/cm^3'. Fields 2 and 3
        # are bare numbers; units are preserved from the existing line.
        if re.search(r"\s", str(value)):
            raise DeckError(
                f"[{self.name}] {label}({mineral!r}) value {value!r} contains "
                f"whitespace. Fields 2 and 3 are bare numbers -- do not append units.")
        s, e = self._constraint_minerals_span()
        blk = self.text[s:e]
        # groups: 1 = indent+name+spaces, 2 = volume fraction, 3 = spaces,
        #         4 = specific surface area, 5 = trailing (units)
        pat = rf"^(\s*{re.escape(mineral)}\s+)(\S+)(\s+)(\S+)(.*)$"
        def rep(m):
            pre, vf, gap, area, rest = m.groups()
            if field_no == 1:
                vf = value
            else:
                area = value
            return pre + vf + gap + area + rest
        new, n = re.subn(pat, rep, blk, flags=re.M)
        if n == 0:
            raise DeckError(f"[{self.name}] {mineral!r} not found in the "
                            f"CONSTRAINT MINERALS block")
        self.text = self.text[:s] + new + self.text[e:]
        self.ops.append({"op": label, "mineral": mineral, "value": value, "matches": n})
        return n

    def set_constraint_vf(self, mineral, value):
        """Volume fraction: field 2 of the CONSTRAINT MINERALS line."""
        return self._set_constraint_field(mineral, 1, value, "set_constraint_vf")

    def set_constraint_area(self, mineral, value):
        """Specific surface area: field 3. NOT in MINERAL_KINETICS -- that block
        holds only PREFACTOR / RATE_CONSTANT / ACTIVATION_ENERGY."""
        return self._set_constraint_field(mineral, 2, value, "set_constraint_area")

    def remove_mineral(self, mineral):
        """Drop a phase from MINERALS, MINERAL_KINETICS and every CONSTRAINT.

        Used for the dawsonite-suppression sensitivity. Removing the phase is
        cleaner than setting its rate to zero: a suppressed-but-present phase
        still participates in the equilibrium system.
        """
        hits, off_scan = 0, 0
        # 1. MINERALS list
        s, e = self._block_span(r"^\s*MINERALS\s*$")
        blk, n = re.subn(rf"^\s*{re.escape(mineral)}\s*$\n", "", self.text[s:e], flags=re.M | re.I)
        self.text, hits = self.text[:s] + blk + self.text[e:], hits + n
        # 2. MINERAL_KINETICS sub-block
        ks, ke = self._block_span(r"^\s*MINERAL_KINETICS\s*$")
        kin = self.text[ks:ke]
        try:
            ms, me = self._block_span(rf"^\s*{re.escape(mineral)}\s*$", text=kin, base=0)
            kin, n = kin[:ms] + kin[me:], 1
        except DeckError:
            n = 0
        self.text, hits = self.text[:ks] + kin + self.text[ke:], hits + n
        # 3. every MINERAL_VOLUME_FRACTIONS block -- SCOPED, not deck-wide.
        # A deck-wide substitution would delete any line beginning with the
        # mineral name, including in OUTPUT or comment blocks.
        n = 0
        while True:
            try:
                vs, ve = self._constraint_minerals_span()
                if vs < off_scan:
                    break
            except DeckError:
                break
            seg, k = re.subn(rf"^\s*{re.escape(mineral)}\s+\S+.*$\n", "",
                             self.text[vs:ve], flags=re.M | re.I)
            self.text = self.text[:vs] + seg + self.text[ve:]
            n += k
            off_scan = vs + len(seg)
        hits += n
        if hits < 2:
            raise DeckError(
                f"[{self.name}] remove_mineral({mineral!r}) only matched "
                f"{hits} place(s); expected >=2 (MINERALS + MINERAL_KINETICS)")
        self.ops.append({"op": "remove_mineral", "mineral": mineral,
                         "matches": hits})
        return hits

    def set_shutin(self, inject_years, rate=None):
        """Stop injection at inject_years by extending the RATE LIST.

        PFLOTRAN holds the last tabulated value indefinitely, which is why the
        published deck's three entries (ending at 1.d-2 yr) inject for the full
        50 yr. Appending (t, rate) then (t*1.0001, 0) makes the source stop.

        rate defaults to the deck's own final tabulated rate, which
        run_pflotran.py scales by domain volume -- so it differs per case and
        must not be hardcoded.
        """
        s0, e0 = self._block_span(r"^\s*RATE LIST\s*$")
        blk = self.text[s0:e0]
        lines = blk.splitlines(keepends=True)
        if rate is None:
            for l in reversed(lines[:-1]):
                p_ = l.split()
                if len(p_) == 2:
                    try:
                        rate = float(p_[1].replace("d", "e").replace("D", "e"))
                        break
                    except ValueError:
                        continue
            if rate is None:
                raise DeckError(f"[{self.name}] could not read the final rate "
                                f"from the RATE LIST")
        ind = self._indent(lines[1]) if len(lines) > 1 else 4
        add = (f"{' '*ind}{inject_years:.6e}   {rate:.6e}\n"
               f"{' '*ind}{inject_years*1.0001:.6e}   0.d0\n")
        self.text = self.text[:s0] + "".join(lines[:-1]) + add + lines[-1] + self.text[e0:]
        self.ops.append({"op": "set_shutin", "inject_years": inject_years,
                         "rate": rate, "matches": 1})
        return 1

    def set_rate_type(self, itype="TOTAL_MASS_RATE"):
        """Switch the injection FLOW_CONDITION's rate type.

        PFLOTRAN condition.F90:1260-1276 distinguishes:
            MASS_RATE        applied PER CELL over the source region
            TOTAL_MASS_RATE  applied as a total across the region
        The published deck uses MASS_RATE on a multi-cell REGION, so the
        effective rate is (deck value) x (cells in region) -- 2.36 to 7820 kg/s
        across the ensemble against a nominal 0.01.
        """
        s0, e0 = self._block_span(r"^\s*FLOW_CONDITION carbonated_water_injection\s*$")
        blk = self.text[s0:e0]
        # the line is "RATE MASS_RATE": sub-condition name then type
        new, n = re.subn(r"^(\s*RATE\s+)(?:TOTAL_|SCALED_)?MASS_RATE\s*$",
                         rf"\g<1>{itype}", blk, flags=re.M)
        if n == 0:
            raise DeckError(f"[{self.name}] no MASS_RATE line in the injection "
                            f"FLOW_CONDITION")
        self.text = self.text[:s0] + new + self.text[e0:]
        self.ops.append({"op": "set_rate_type", "itype": itype, "matches": n})
        return n

    def set_rate_schedule(self, entries):
        """Replace the RATE LIST table. entries: [(time_yr, rate_kg_s), ...]."""
        t = [e[0] for e in entries]
        if t != sorted(t):
            raise DeckError(f"[{self.name}] rate schedule times not ascending: {t}")
        s0, e0 = self._block_span(r"^\s*RATE LIST\s*$")
        lines = self.text[s0:e0].splitlines(keepends=True)
        head = [l for l in lines[1:-1] if "UNITS" in l.upper()]
        ind = self._indent(lines[1]) if len(lines) > 1 else 4
        body = "".join(f"{' '*ind}{a:.6e}   {b:.6e}\n" for a, b in entries)
        self.text = (self.text[:s0] + lines[0] + "".join(head) + body
                     + lines[-1] + self.text[e0:])
        self.ops.append({"op": "set_rate_schedule", "entries": entries, "matches": 1})
        return 1

    def set_permeability_dataset(self, h5="dfn_properties.h5",
                                 dsname="Permeability", label="permeability"):
        """Replace PERM_ISO with the aperture-derived per-cell field.

        lagrit2pflotran() writes dfn_properties.h5 with per-cell permeability
        from the fracture apertures (cubic law). The published deck discarded it
        via a single PERM_ISO 8.33d-8, so Table 1's aperture distributions
        affected nothing. Apertures span 2.3e-4 to 2.4e-3 m -- a ~100x
        permeability range -- so this makes aperture heterogeneity active and
        answers R2-4/R3-4.
        """
        s0, e0 = self._block_span(r"^\s*PERMEABILITY\s*$")
        blk = self.text[s0:e0]
        new, n = re.subn(r"^(\s*)PERM_ISO\s+\S+\s*$", rf"\g<1>DATASET {label}",
                         blk, flags=re.M)
        if n == 0:
            raise DeckError(f"[{self.name}] no PERM_ISO inside PERMEABILITY")
        self.text = self.text[:s0] + new + self.text[e0:]
        self.append_block(f"DATASET {label}\n  FILENAME {h5}\n"
                          f"  HDF5_DATASET_NAME {dsname}\nEND")
        self.ops.append({"op": "set_permeability_dataset", "h5": h5,
                         "dataset": dsname, "matches": n})
        return n

    def set_solver_value(self, method, block, key, value):
        """Set a value in NUMERICAL_METHODS <FLOW|TRANSPORT> / <block>.

        The published deck gives flow 25 Newton iterations and transport 250.
        With corrected cell volumes (1e-6 to 2.8e-3 m3 rather than 3.9e-3 to
        1.66) the initial pressure solve exceeds 25 and reports
        SNES_DIVERGED_MAX_IT with FLOW TS SNES steps = 0. The ATOL of 1.d-6 is
        also absolute on a residual whose scale fell with the volumes.
        """
        ms, me = self._block_span(rf"^\s*NUMERICAL_METHODS {method}\s*$")
        seg = self.text[ms:me]
        bs, be = self._block_span(rf"^\s*{block}\s*$", text=seg, base=0)
        body, n = re.subn(rf"^(\s*{re.escape(key)}\s+).*$", rf"\g<1>{value}",
                          seg[bs:be], flags=re.M)
        if n == 0:
            raise DeckError(f"[{self.name}] {key!r} not found in "
                            f"NUMERICAL_METHODS {method}/{block}")
        self.text = self.text[:ms] + seg[:bs] + body + seg[be:] + self.text[me:]
        self.ops.append({"op": "set_solver_value", "method": method,
                         "block": block, "key": key, "value": value,
                         "matches": n})
        return n

    def set_mineral_rate(self, mineral, value, which=0):
        """Set a RATE_CONSTANT inside a mineral's MINERAL_KINETICS entry.

        Table 2 attributes every rate to Palandri and Kharaka (2004), but that
        compilation does not include dawsonite: its rate of 1.0e-7 mol/m2/s is
        a round number, and at 50 C that is 7.1e-7, comparable to calcite and
        about seven orders above the only other secondary aluminosilicate in
        the set. Dawsonite is 40% of the ensemble carbonate budget, so how much
        of that share depends on the assumed rate needs bounding.

        `which` selects the prefactor when a mineral has more than one (an acid
        and a neutral mechanism); dawsonite has only one.
        """
        ms, me = self._block_span(r"^\s*MINERAL_KINETICS\s*$")
        seg = self.text[ms:me]
        bs, be = self._block_span(rf"^\s*{re.escape(mineral)}\s*$",
                                  text=seg, base=0)
        blk = seg[bs:be]
        hits = list(re.finditer(r"^(\s*RATE_CONSTANT\s+)(\S+)(\s+\S+)\s*$",
                                blk, re.M))
        if not hits:
            raise DeckError(f"[{self.name}] no RATE_CONSTANT for {mineral}")
        if which >= len(hits):
            raise DeckError(f"[{self.name}] {mineral} has {len(hits)} "
                            f"prefactor(s), asked for index {which}")
        h = hits[which]
        blk = blk[:h.start()] + f"{h.group(1)}{value}{h.group(3)}" + blk[h.end():]
        self.text = self.text[:ms] + seg[:bs] + blk + seg[be:] + self.text[me:]
        self.ops.append({"op": "set_mineral_rate", "mineral": mineral,
                         "value": value, "which": which, "matches": 1})
        return 1

    def insert_in_block(self, block_header, lines):
        """Add lines immediately before a block's OWN terminator.

        The terminator is the last line of the span returned by _block_span --
        NOT the first END/'/' found by searching, which would be an inner
        sub-block's closer. Getting this wrong put UPDATE_POROSITY inside the
        PERMEABILITY sub-block, where PFLOTRAN ignores it and the sensitivity
        run silently reproduces the baseline.
        """
        s, e = self._block_span(block_header)
        blk = self.text[s:e]
        parts = blk.splitlines(keepends=True)
        # A genuine block has at least one child line between header and
        # terminator. Exactly two lines means either an empty block or -- far
        # more likely -- a leaf statement misidentified as a block header
        # (e.g. `DATABASE hanford.dat`, whose span runs to the enclosing END).
        # Both are mis-specifications; refuse rather than insert somewhere odd.
        if len(parts) < 3:
            raise DeckError(
                f"[{self.name}] {block_header!r} has no child lines "
                f"({len(parts)} line(s) in span) -- not a block, or empty. "
                f"Refusing to insert.")
        # Guard against treating a NON-block keyword as a block. _block_span
        # will happily return a span running to the enclosing block's
        # terminator, e.g. `DATABASE hanford.dat` -> ... -> `END`, which would
        # insert lines in a meaningless place. A real block's first child is
        # strictly more indented than its header.
        if self._indent(parts[1]) <= self._indent(parts[0]):
            raise DeckError(
                f"[{self.name}] {block_header!r} does not look like a block: "
                f"first child is not indented deeper than the header. "
                f"Refusing to insert.")
        head, term = "".join(parts[:-1]), parts[-1]
        if not re.match(r"^\s*(END|/)\s*$", term, re.I):
            raise DeckError(
                f"[{self.name}] last line of {block_header!r} is not a terminator: {term!r}")
        ind = self._indent(parts[0]) + 2          # child indent, not hardcoded
        add = "".join(" " * ind + l.strip() + "\n" for l in lines)
        self.text = self.text[:s] + head + add + term + self.text[e:]
        self.ops.append({"op": "insert_in_block", "block": block_header,
                         "lines": lines, "indent": ind, "matches": 1})
        return 1

    def append_block(self, text):
        """Append a top-level block, e.g. MULTIPLE_CONTINUUM. Inserted before
        the final END of SUBSURFACE if present, else at end of file."""
        anchor = list(re.finditer(r"^END_SUBSURFACE\s*$", self.text, re.M))
        if anchor:
            i = anchor[-1].start()
            self.text = self.text[:i] + text.rstrip() + "\n\n" + self.text[i:]
        else:
            self.text = self.text.rstrip() + "\n\n" + text.rstrip() + "\n"
        self.ops.append({"op": "append_block",
                         "head": text.strip().splitlines()[0], "matches": 1})
        return 1

    # ------------------------------------------------------- audit trail
    def diff(self, n=3):
        return "".join(difflib.unified_diff(
            self.original.splitlines(keepends=True),
            self.text.splitlines(keepends=True),
            fromfile=f"{self.name} (baseline)", tofile=f"{self.name} (variant)",
            n=n))

    def provenance(self):
        return {"deck": self.name, "n_ops": len(self.ops), "ops": self.ops}

    def __repr__(self):
        return f"<Deck {self.name} ops={len(self.ops)}>"


# ------------------------------------------------------------------ variants
def apply_variant(deck: Deck, ops: list) -> Deck:
    """Apply a variant spec (list of (method, args, kwargs)) to a deck."""
    for spec in ops:
        method, args = spec[0], spec[1] if len(spec) > 1 else ()
        kwargs = spec[2] if len(spec) > 2 else {}
        getattr(deck, method)(*args, **kwargs)
    return deck


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        sys.exit("usage: deckmod.py <deck.in>   # self-test on a real deck")
    d = Deck.from_file(sys.argv[1])
    print("block spans found:")
    for h in [r"^\s*MINERALS\s*$", r"^\s*MINERAL_KINETICS\s*$",
              r"^\s*MATERIAL_PROPERTY fracture\b",
              r"^\s*MINERAL_VOLUME_FRACTIONS\s*$"]:
        try:
            s, e = d._block_span(h)
            print(f"  OK   {h:<42s} {e-s:6d} chars")
        except DeckError as ex:
            print(f"  FAIL {h:<42s} {ex}")
