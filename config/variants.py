"""
Block C sensitivity variants, as declarative deck edits.

Each entry is a list of (method, args) tuples applied by deckmod to the BASELINE
deck that run_pflotran.py emits. The list IS the definition of the variant --
there is no second copy of the deck to diverge.

Add a variant here, not by editing run_pflotran.py.
"""

# Cases the sensitivity suite runs on. Geochemical tests do not need the
# largest mesh: the x1.00 pair is ~84k cells vs ~166k for x1.50, which roughly
# halves Block C's cost with no loss of inferential value.
PAIR_CHEAP = ["p32_100_s383", "p32_100_s42"]    # high- and low-trapping
PAIR_FULL  = ["p32_150_s117", "p32_150_s42"]    # the pair quoted in the paper

VARIANTS = {
    # -- R3-6, R1-2c: anorthite reactive surface area -----------------------
    **{f"anor_as{v}": {
        "comment": f"Anorthite SPECIFIC_SURFACE_AREA = {v} cm^2/cm^3 "
                   f"(published value 10; R3 asked for 30-50)",
        "cases": PAIR_FULL,
        "ops": [("set_constraint_area", ("Anorthite", f"{v}.d0"))],
       } for v in (30, 50, 100)},

    # -- R1-2c: the reduction argument is not anorthite-specific ------------
    **{f"global_as_x{tag}": {
        "comment": f"All primary surface areas scaled x{tag} "
                   f"-- reframes the tuning as a Damkohler sensitivity",
        "cases": PAIR_FULL,
        "ops": [("set_constraint_area", (m, f"{a}.d0"))
                for m, a in zip(
                    ["Anorthite","Albite","Diopside","Forsterite","Fayalite","Enstatite"],
                    ([1,10,10,10,10,10] if tag == "0.1" else [100,1000,1000,1000,1000,1000]))],
       } for tag in ("0.1", "10")},

    # -- R2-1: mineral VFs constrained by 1 - phi ---------------------------
    "vf_consistent": {
        "comment": "Primary VFs rescaled by 0.5/0.85 so solid + void = 1. "
                   "Partially-filled-fracture conceptual model.",
        "cases": PAIR_FULL,
        "ops": [("set_constraint_vf", (m, f"{v*0.5/0.85:.4f}d0"))
                for m, v in [("Anorthite",0.30), ("Diopside",0.25), ("Albite",0.20),
                             ("Enstatite",0.05), ("Forsterite",0.03), ("Fayalite",0.02)]],
    },

    # -- R2-4, R3-4: aperture feedback --------------------------------------
    "feedback": {
        "comment": "Porosity-permeability feedback enabled. Tests whether "
                   "dissolution-driven opening amplifies channelling.",
        "cases": PAIR_FULL + ["p32_200_s117", "p32_200_s383"],
        # reaction.F90:987 -- "UPDATE_POROSITY must be listed under CHEMISTRY".
        # The full set at reaction.F90:1130 is UPDATE_POROSITY,
        # UPDATE_TORTUOSITY, UPDATE_PERMEABILITY.
        "ops": [("insert_in_block",
                 (r"^\s*CHEMISTRY\s*$",
                  ["UPDATE_POROSITY", "UPDATE_PERMEABILITY"]))],
    },

    # -- R1-2a, R1-8: dawsonite suppression --------------------------------
    "no_dawsonite": {
        "comment": "Dawsonite removed entirely (not merely rate-suppressed: a "
                   "present-but-suppressed phase still enters the equilibrium "
                   "system). One case per P32 level.",
        "cases": ["p32_075_s117","p32_100_s383","p32_125_s383",
                  "p32_150_s117","p32_200_s117"],
        "ops": [("remove_mineral", ("Dawsonite",))],
    },

    # -- R3-5, R1-2b: nucleation surrogate sensitivity ----------------------
    **{f"seed_vf_{tag}": {
        "comment": f"Secondary seed volume fraction {val} (published 1e-6). "
                   f"Brackets nucleation-kinetics uncertainty within TST.",
        "cases": PAIR_CHEAP,
        "ops": [("set_constraint_vf", (m, val)) for m in
                ["Calcite","Magnesite","Siderite","Dawsonite","Kaolinite","Chalcedony"]],
       } for tag, val in [("lo","1.d-8"), ("hi","1.d-4")]},

    **{f"sec_as_{tag}": {
        "comment": f"Secondary reactive surface area {val} cm^2/cm^3 (published 1).",
        "cases": PAIR_CHEAP,
        "ops": [("set_constraint_area", (m, f"{val}"))
                for m in ["Calcite","Magnesite","Siderite","Dawsonite"]],
       } for tag, val in [("lo","0.1d0"), ("hi","10.d0")]},

    # -- R1-10: shut-in rather than WAG. Inject for 2 yr, then stop and let the
    # system react for 48 yr. Tests the prediction in the paper's own section
    # 4.3/4.5 that a finite pulse raises conversion and suppresses the
    # re-dissolution seen under continuous injection.
    "shutin_2yr": {
        "comment": "R1-10  inject 2 yr then shut in for 48 yr",
        "cases": ["p32_075_s117", "p32_100_s383", "p32_125_s383",
                  "p32_150_s117", "p32_150_s42", "p32_200_s117"],
        "ops": [("set_shutin", (2.0,))],
    },

    # -- R2-1 + R2-4/R3-4 combined. Feedback alone diverges because the
    # published volume fractions (0.85) against porosity 0.50 drive phi past 1.0
    # once ~44% of the solid dissolves. With VFs rescaled to sum to 1-phi,
    # dissolving 44% gives phi 0.50 -> 0.72 and the coupled solve should hold.
    "feedback_vfc": {
        "comment": "R2-4, R3-4, R2-1  feedback WITH corrected volume fractions",
        "cases": ["p32_150_s117", "p32_150_s42"],
        "ops": [("set_constraint_vf", (m, f"{v*0.5/0.85:.4f}d0")) for m, v in
                [("Anorthite",0.30),("Diopside",0.25),("Albite",0.20),
                 ("Enstatite",0.05),("Forsterite",0.03),("Fayalite",0.02)]]
               + [("insert_in_block", (r"^\s*CHEMISTRY\s*$",
                                       ["UPDATE_POROSITY", "UPDATE_PERMEABILITY"]))],
    },

    # -- R1-10, R3-8. CarbFix1 Phase I injected 175 t CO2 over 45 days
    # (24 Jan - 9 Mar 2012); >95% mineralised within ~2 yr of monitoring.
    # The 2-yr pulse is ~16x longer than the field injection, so this matches
    # the schedule the 95% figure refers to.
    "shutin_carbfix": {
        "comment": "R1-10, R3-8  45-day pulse (CarbFix1 Phase I) then shut in",
        "cases": ["p32_075_s117", "p32_100_s383", "p32_150_s117", "p32_150_s42"],
        "ops": [("set_shutin", (0.123,))],
    },

    # -- field-scaled: 2% of domain pore volume over a 45-day pulse, then
    # shut in. TOTAL_MASS_RATE so the deck value is the actual total rate.
    # CarbFix1 filled ~1.6% of its monitored pore volume, so 2% is comparable.
    "fieldscale_002pct": {
        "comment": "R1-10, R3-8  2% pore-volume fill, 45-day pulse, TOTAL_MASS_RATE",
        "cases": ["p32_100_s383"],
        "ops": [("set_rate_type", ("SCALED_MASS_RATE VOLUME",)),
                ("set_rate_schedule", ([(0.0, 0.0), (1e-4, 1.919078e-03),
                                        (1e-2, 1.919078e-02), (0.123, 1.919078e-02),
                                        (0.1230123, 0.0)],))],
    },
    # -- field-scaled: 10% of domain pore volume over a 45-day pulse, then
    # shut in. TOTAL_MASS_RATE so the deck value is the actual total rate.
    # CarbFix1 filled ~1.6% of its monitored pore volume, so 2% is comparable.
    "fieldscale_010pct": {
        "comment": "R1-10, R3-8  10% pore-volume fill, 45-day pulse, TOTAL_MASS_RATE",
        "cases": ["p32_100_s383"],
        "ops": [("set_rate_type", ("SCALED_MASS_RATE VOLUME",)),
                ("set_rate_schedule", ([(0.0, 0.0), (1e-4, 9.595391e-03),
                                        (1e-2, 9.595391e-02), (0.123, 9.595391e-02),
                                        (0.1230123, 0.0)],))],
    },
    # -- field-scaled: 25% of domain pore volume over a 45-day pulse, then
    # shut in. TOTAL_MASS_RATE so the deck value is the actual total rate.
    # CarbFix1 filled ~1.6% of its monitored pore volume, so 2% is comparable.
    "fieldscale_025pct": {
        "comment": "R1-10, R3-8  25% pore-volume fill, 45-day pulse, TOTAL_MASS_RATE",
        "cases": ["p32_100_s383"],
        "ops": [("set_rate_type", ("SCALED_MASS_RATE VOLUME",)),
                ("set_rate_schedule", ([(0.0, 0.0), (1e-4, 2.398848e-02),
                                        (1e-2, 2.398848e-01), (0.123, 2.398848e-01),
                                        (0.1230123, 0.0)],))],
    },
    # -- field-scaled: 100% of domain pore volume over a 45-day pulse, then
    # shut in. TOTAL_MASS_RATE so the deck value is the actual total rate.
    # CarbFix1 filled ~1.6% of its monitored pore volume, so 2% is comparable.
    "fieldscale_100pct": {
        "comment": "R1-10, R3-8  100% pore-volume fill, 45-day pulse, TOTAL_MASS_RATE",
        "cases": ["p32_100_s383"],
        "ops": [("set_rate_type", ("SCALED_MASS_RATE VOLUME",)),
                ("set_rate_schedule", ([(0.0, 0.0), (1e-4, 9.595391e-02),
                                        (1e-2, 9.595391e-01), (0.123, 9.595391e-01),
                                        (0.1230123, 0.0)],))],
    },

    # -- CORRECTED CONFIGURATION -------------------------------------------
    # Requires meshes regenerated with lagrit2pflotran(): volumes are physical
    # (6.997 m3 vs 7550.9 m2 for p32_100_s383) and dfn_properties.h5 carries
    # aperture-derived per-cell permeability (apertures 2.3e-4 to 2.4e-3 m).
    "corrected_carbfix": {
        "comment": "corrected volumes + aperture permeability + 45-day pulse",
        "cases": ["p32_100_s383", "p32_150_s117", "p32_150_s42"],
        "ops": [("set_permeability_dataset", ()),
                ("set_rate_type", ("SCALED_MASS_RATE VOLUME",)),
                ("set_rate_schedule", ([(0.0, 0.0), (1e-4, 3.6e-7),
                                        (1e-2, 3.6e-6), (0.123, 3.6e-6),
                                        (0.1230123, 0.0)],))],
    },
    "corrected_continuous": {
        "comment": "corrected volumes + aperture permeability + continuous",
        "cases": ["p32_100_s383", "p32_150_s117", "p32_150_s42"],
        "ops": [("set_permeability_dataset", ()),
                ("set_rate_type", ("SCALED_MASS_RATE VOLUME",))],
    },

    # Corrected volumes, but PERM_ISO retained. Isolates the volume correction
    # from the permeability heterogeneity: corrected_continuous changes both
    # and fails at the first flow solve (SNES steps = 0).
    "corrected_uniformperm": {
        "comment": "corrected volumes, uniform PERM_ISO, continuous injection",
        "cases": ["p32_100_s383"],
        "ops": [("set_rate_type", ("SCALED_MASS_RATE VOLUME",))],
    },

    # Corrected volumes + aperture permeability + relaxed FLOW solver limits.
    # The published deck allows flow 25 Newton iterations against transport's
    # 250; on the corrected geometry the initial pressure solve exceeds 25 and
    # returns SNES_DIVERGED_MAX_IT with zero completed steps.
    "corrected_solver": {
        "comment": "corrected mesh + relaxed flow solver, 45-day pulse",
        "cases": ["p32_100_s383", "p32_150_s117"],
        "ops": [("set_permeability_dataset", ()),
                ("set_rate_type", ("SCALED_MASS_RATE VOLUME",)),
                ("set_solver_value", ("FLOW", "NEWTON_SOLVER",
                                      "MAXIMUM_NUMBER_OF_ITERATIONS", "250")),
                ("set_solver_value", ("FLOW", "NEWTON_SOLVER", "ATOL", "1.d-10")),
                ("set_rate_schedule", ([(0.0, 0.0), (1e-4, 3.6e-7),
                                        (1e-2, 3.6e-6), (0.123, 3.6e-6),
                                        (0.1230123, 0.0)],))],
    },

    # THE CONFIGURATION FOR THE RE-RUNS.
    # Corrected volumes (lagrit2pflotran) + SCALED_MASS_RATE + uniform PERM_ISO.
    # The aperture-derived permeability field (112x contrast across cells whose
    # volumes span 2700x) prevents the initial flow solve from converging even
    # at 250 Newton iterations; uniform permeability converges with the
    # published solver settings. Aperture heterogeneity therefore remains a
    # stated limitation, as in the original manuscript.
    "corrected_baseline": {
        "comment": "corrected volumes, uniform perm, continuous injection",
        "cases": ["p32_075_s117", "p32_100_s383", "p32_125_s383",
                  "p32_150_s117", "p32_150_s42", "p32_200_s117"],
        "ops": [("set_rate_type", ("SCALED_MASS_RATE VOLUME",))],
    },

    # THE CONFIGURATION FOR THE RE-RUNS.
    # Corrected volumes + uniform PERM_ISO + field-scaled 45-day pulse.
    # Rationale: at the published 0.01 kg/s the corrected 3.5 m3 pore volume is
    # flushed ~4500 times in 50 yr, which is both physically implausible and
    # numerically very slow (119 timesteps for 0.01 yr). 3.6e-6 kg/s over 45
    # days gives ~0.02 pore volumes, comparable to CarbFix1's ~1.6% fill.
    # Uniform permeability because the aperture-derived field (112x contrast)
    # prevents the initial flow solve converging even at 250 Newton iterations.
    "corrected_field": {
        "comment": "corrected volumes, uniform perm, field-scaled 45-day pulse",
        "cases": ["p32_075_s117", "p32_100_s383", "p32_125_s383",
                  "p32_150_s117", "p32_150_s42", "p32_200_s117"],
        # The published deck allows flow 25 Newton iterations against
        # transport's 250. On corrected geometry the initial pressure solve
        # exceeds 25 above ~126k cells (55k and 84k converge; 126k, 162k, 172k
        # and 231k return SNES_DIVERGED_MAX_IT with zero completed steps). ATOL
        # is absolute on a residual whose scale fell with the cell volumes.
        "ops": [("set_solver_value", ("FLOW", "NEWTON_SOLVER",
                                      "MAXIMUM_NUMBER_OF_ITERATIONS", "250")),
                # ATOL relaxed from 1.d-6 to 1.d-4. The flow residual stagnates
                # at 2r ~ 1.39e-6 against a solution norm 2x ~ 2.01e9, i.e. a
                # RELATIVE residual of 7e-16 -- machine precision. The residual
                # infinity norm (ir 6.04e-08) is already well inside 1.d-6. So
                # the solve converges and fails only an absolute test whose
                # floor scales with the solution norm, which is why 55k and 84k
                # cells pass while 126k+ do not.
                ("set_solver_value", ("FLOW", "NEWTON_SOLVER", "ATOL", "1.d-4")),
                ("set_rate_type", ("SCALED_MASS_RATE VOLUME",)),
                ("set_rate_schedule", ([(0.0, 0.0), (1e-4, 3.6e-7),
                                        (1e-2, 3.6e-6), (0.123, 3.6e-6),
                                        (0.1230123, 0.0)],))],
    },

    # Seeding sensitivity on the CORRECTED, field-scaled configuration. The
    # earlier test (1e-8 to 1e-4) changed trapping by <=1%, but that was under
    # ~57 million pore volumes of flushing where precipitation was flux-limited.
    # At the field-scaled rate the system is dissolution-limited and the
    # mineralogy has flipped to magnesite, so the seeding dependence may differ.
    "corrected_seed_lo": {
        "comment": "corrected field-scaled + secondary seed VF 1.d-10",
        "cases": ["p32_100_s383", "p32_150_s117"],
        "ops": [("set_rate_type", ("SCALED_MASS_RATE VOLUME",)),
                ("set_rate_schedule", ([(0.0, 0.0), (1e-4, 3.6e-7),
                                        (1e-2, 3.6e-6), (0.123, 3.6e-6),
                                        (0.1230123, 0.0)],))]
               + [("set_constraint_vf", (m, "1.d-10")) for m in
                  ("Calcite","Magnesite","Siderite","Dawsonite",
                   "Kaolinite","Chalcedony")],
    },

    # Seeding sensitivity on the CORRECTED, field-scaled configuration. The
    # earlier test (1e-8 to 1e-4) changed trapping by <=1%, but that was under
    # ~57 million pore volumes of flushing where precipitation was flux-limited.
    # At the field-scaled rate the system is dissolution-limited and the
    # mineralogy has flipped to magnesite, so the seeding dependence may differ.
    "corrected_seed_off": {
        "comment": "corrected field-scaled + secondary seed VF 0.d0",
        "cases": ["p32_100_s383", "p32_150_s117"],
        "ops": [("set_rate_type", ("SCALED_MASS_RATE VOLUME",)),
                ("set_rate_schedule", ([(0.0, 0.0), (1e-4, 3.6e-7),
                                        (1e-2, 3.6e-6), (0.123, 3.6e-6),
                                        (0.1230123, 0.0)],))]
               + [("set_constraint_vf", (m, "0.d0")) for m in
                  ("Calcite","Magnesite","Siderite","Dawsonite",
                   "Kaolinite","Chalcedony")],
    },

    # The 25 published cases with NO variant applied -- only the mandatory
    # corrections from apply_corrections.py. Needed as the reference for every
    # Block C ratio: comparing a corrected variant against the uncorrected
    # published value in betweenness_results.csv is not a valid comparison.
    "baseline": {
        "comment": "the published 25, corrected volumes, no variant",
        "cases": ['p32_075_s42', 'p32_075_s117', 'p32_075_s259', 'p32_075_s383', 'p32_075_s501', 'p32_100_s42', 'p32_100_s117', 'p32_100_s259', 'p32_100_s383', 'p32_100_s501', 'p32_125_s42', 'p32_125_s117', 'p32_125_s259', 'p32_125_s383', 'p32_125_s501', 'p32_150_s42', 'p32_150_s117', 'p32_150_s259', 'p32_150_s383', 'p32_150_s501', 'p32_200_s42', 'p32_200_s117', 'p32_200_s259', 'p32_200_s383', 'p32_200_s501'],
        "ops": [],
    },

    # Dawsonite rate sensitivity (two orders slower). Its published rate of 1.0e-7
    # mol/m2/s is not from Palandri and Kharaka and is seven orders above the
    # only other secondary aluminosilicate in the set, yet dawsonite is 40% of
    # the ensemble carbonate budget. This bounds how much of that share is
    # rate-controlled. It also answers R2-2: zeolites are present in
    # hanford.dat, but published analcime dissolution rates are 5 to 6 orders
    # below the rate assumed for dawsonite, so analcime cannot compete for Na
    # and Al on a 50-year timescale within a kinetic framework.
    "daw_rate_1e9": {
        "comment": "dawsonite RATE_CONSTANT -> 1.0e-9 (two orders slower)",
        "cases": ["p32_100_s383", "p32_125_s383", "p32_150_s501"],
        "ops": [("set_mineral_rate", ("Dawsonite", "1.0e-9"))],
    },

    # Dawsonite rate sensitivity (four orders slower, near albite neutral). Its published rate of 1.0e-7
    # mol/m2/s is not from Palandri and Kharaka and is seven orders above the
    # only other secondary aluminosilicate in the set, yet dawsonite is 40% of
    # the ensemble carbonate budget. This bounds how much of that share is
    # rate-controlled. It also answers R2-2: zeolites are present in
    # hanford.dat, but published analcime dissolution rates are 5 to 6 orders
    # below the rate assumed for dawsonite, so analcime cannot compete for Na
    # and Al on a 50-year timescale within a kinetic framework.
    "daw_rate_1e11": {
        "comment": "dawsonite RATE_CONSTANT -> 1.0e-11 (four orders slower, near albite neutral)",
        "cases": ["p32_100_s383", "p32_125_s383", "p32_150_s501"],
        "ops": [("set_mineral_rate", ("Dawsonite", "1.0e-11"))],
    },

    # Dawsonite rate sensitivity (six orders slower, near kaolinite and published analcime rates). Its published rate of 1.0e-7
    # mol/m2/s is not from Palandri and Kharaka and is seven orders above the
    # only other secondary aluminosilicate in the set, yet dawsonite is 40% of
    # the ensemble carbonate budget. This bounds how much of that share is
    # rate-controlled. It also answers R2-2: zeolites are present in
    # hanford.dat, but published analcime dissolution rates are 5 to 6 orders
    # below the rate assumed for dawsonite, so analcime cannot compete for Na
    # and Al on a 50-year timescale within a kinetic framework.
    "daw_rate_1e13": {
        "comment": "dawsonite RATE_CONSTANT -> 1.0e-13 (six orders slower, near kaolinite and published analcime rates)",
        "cases": ["p32_100_s383", "p32_125_s383", "p32_150_s501"],
        "ops": [("set_mineral_rate", ("Dawsonite", "1.0e-13"))],
    },
}

# Variants that need more than a deck edit -- flagged rather than silently absent
NEEDS_MORE_WORK = {
    "matrix_diffusion": "Requires a MULTIPLE_CONTINUUM block AND a matrix "
                        "mineral assemblage. Draft with append_block, then "
                        "verify against the PFLOTRAN manual before trusting it.",
    "expanded_secondary": "Adding zeolite/smectite/Fe-sink phases requires "
                          "thermodynamic data. CHECK hanford.dat CONTAINS THEM "
                          "before designing the run -- it may not.",
    "redox": "Needs the O2(aq) constraint changed in the CONSTRAINT block and a "
             "decision on the buffering assemblage.",
    "pulsed": "Needs injection scheduling, not a deck edit. run_pflotran.py "
              "exposes inject_years; may be usable directly.",
}
