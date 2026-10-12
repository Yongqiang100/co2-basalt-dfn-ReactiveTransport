# Results summary for the revision of 2026WR044716

Corrected simulations, compiled 28 September 2026.

**Where the data are.** Everything is on hpc02 in `/scratch/chen/co2-basalt/revision`. Paths below start from there.

| Folder | Contents |
|---|---|
| `runs_gravityoff/` | the 362 corrected runs with results |
| `runs_gravityoff_setonix_overlap/` | Setonix copies of 77 runs that also ran on hpc02, including all of fixed Block C |
| `runs_gravityoff_prefloor/` | 8 runs from before the porosity limit was added |
| `runs_surface_area_hpc01/` | the 20 surface-area runs (`S0`, `S1`) |
| `runs_dirichlet_20260924_setonix/` | the old runs with the boundary error (source of the old numbers) |
| `results/`, `archive_setonix/results/`, `archive_hpc01/results/` | result files from hpc02, Setonix and hpc01 |

**Terms used below.**
- *Coupled model*: porosity and permeability change as minerals dissolve and form. This is the model we report.
- *Fixed model*: porosity stays constant. We report it as a limiting case.
- *Efficiency*: the share of the injected carbon that ends up in carbonate minerals.
- *Water age*: how long water has been in the network since it was injected, calculated from the flow field.
- *Flow index* (shape99): the age of the oldest 1% of the water divided by the median age. A high value means the network has zones where water moves slowly and stays long.
- *Old runs*: the runs with the boundary error, before gravity was switched off.

**Source column.** Each value names the script and the results file that produce it. "Inline" marks the values still computed by a session command (see section 9).

---

## 1. What changes in the paper

| Finding | Evidence | What it means for the paper |
|---|---|---|
| Different minerals dissolve at different times, and each forms its own carbonate | Forsterite dissolves in the first 2 years. Its magnesium forms magnesite (49% of the carbonate at 2 years). Diopside and anorthite dissolve over decades. Their calcium forms calcite (85% at 50 years). Injections of 45 days or less form 74 to 86% magnesite. Injections of 2 years or more form 86 to 91% calcite | The carbonate type shows which mineral supplied the element. The carbonate mix depends on time and on injection duration. This sequence belongs in the Results and the Abstract |
| Denser networks form more carbonate because they contain more pore space. Their efficiency is the same | Carbonate increases with fracture intensity (ρ = 0.66). Pore space increases in the same way (ρ = 0.97). Efficiency does not change (ρ = −0.06, p = 0.67) | The claim that fracture intensity controls efficiency is not supported |
| The co-location of three conditions controls mineralization: dissolved calcium or magnesium, carbonate ion, and a high pH | Calcium and magnesium are abundant: 98.8% of the released amount is not captured. The injected water supplies the carbon. A high pH is rare. Cells where all three conditions exceed the network mean contain 94.7 to 99.3% of the new carbonate and represent 4.0 to 7.7% of the volume (section 2.1) | The co-location mapping identifies the mineralization condition. It explains the low efficiency: during injection, the three conditions coincide in a small part of the network |
| The flow-field test requested by the reviewers gives a statistical association | The test specified in advance gives ρ = 0.84 between the flow index and the efficiency (20 networks, p = 1.5 × 10⁻⁶). The index does not depend on fracture intensity | The flow field carries no geochemical information and does not identify the mineralization condition. The letter reports the test result, as the reviewers requested, without using it as evidence for the mechanism |
| Calcium and magnesium never limit the carbonate | During continuous injection, 1.25% of the released calcium, magnesium and iron forms carbonate (median, 50 networks). After shut-in, 13.7% forms carbonate, and the carbonate equals 108% of the dissolved carbon | During injection, the low pH limits the carbonate. After shut-in, the dissolved carbon limits it. Reverses the letter's claim that the supply of calcium and magnesium limits the carbonate |
| Kaolinite fills more pore space than the carbonate | Kaolinite and chalcedony fill 24.1% of the initial pore volume at 50 years (median 11,630 mol of kaolinite against 642 mol of carbonate) | Kaolinite belongs in the answers on secondary minerals (R2-2) and on porosity and clogging |
| Stopping injection lets most dissolved carbon turn into carbonate | Carbonate at 50 years equals 108% of the carbon dissolved in the water at shut-in (median). Stopping gives 6.2 times the carbonate of continuous injection | Reverses the old shut-in result (1.1 times) |
| The rock's reactivity sets the amount of carbonate formed | Ten times less reactive surface: 8% of the carbonate. Ten times more: 202% | Flow produces no carbonate in rock without reactive minerals. All model cells start with the same minerals. The flow results apply to a rock of uniform composition |

**Title.** The current title says fracture network connectivity controls efficiency. The results support a different statement. Carbonate forms where dissolved calcium or magnesium, carbonate ion and a high pH coincide. The co-location mapping identifies this condition and explains the low efficiency. During injection, the three conditions coincide in a small part of the network. After injection stops, the pH increases across the network. The co-located volume then expands.

---

## 2. The mineralization condition: co-location of three conditions (AE-1, R2-5, R3-2, R1-11, R1-15)

### 2.1 Carbonate forms where three conditions coincide

Carbonate forms where three conditions coincide in the same cell:

1. Dissolved calcium or magnesium, released by the dissolving minerals.
2. Carbonate ion, from the injected carbon.
3. A pH high enough for carbonate to be stable.

The co-location mapping of the three conditions identifies the mineralization condition. The mapping also explains the low efficiency:

- The minerals release calcium and magnesium across the network (section 3.4). During continuous injection, 98.8% of the released calcium, magnesium and iron is not captured in carbonate (section 6.3).
- The injected water supplies the carbon.
- A high pH is rare during injection. From five weeks on, 0.7 to 3.3% of the pore space has a pH above 6 (section 3.5).
- The three conditions coincide in this small part of the network. The carbonate forms there.
- After injection stops, the minerals neutralize the pore water. The pH increases to 9.1. The co-located volume extends across the network. The carbonate increases until the dissolved carbon is used up (section 6.1).

**Co-location mapping on the corrected runs** (coupled Block A, medians of 50 networks). The three conditions are measured at the start of each interval, before the carbonate of that interval forms.

| Interval | Rank of new carbonate in Ca or Mg | in carbonate ion | in pH | New carbonate in cells with all three high | Volume of these cells | New carbonate where one condition is high and another low: Ca or Mg / carbonate ion / pH |
|---|---|---|---|---|---|---|
| 1 to 2 y | 96.8 | 98.5 | 97.4 | 96.6% | 6.6% | 0.8 / 2.2 / 3.4% |
| 5 to 7 y | 98.7 | 98.5 | 97.2 | 94.7% | 4.0% | 4.7 / 1.1 / 3.3% |
| 10 to 15 y | 98.5 | 98.1 | 96.5 | 96.8% | 6.0% | 2.6 / 1.5 / 1.7% |
| 20 to 30 y | 98.6 | 97.9 | 97.6 | 99.3% | 7.7% | 0.6 / 0.2 / 0.5% |

Rank: 50 means no preference, 100 means the highest values. "High" means above the network's volume-weighted mean. Source: `src/colocation_mapping.py` → `results/colocation_mapping_carbonate.csv`, `results/colocation_mapping.log`.

The outputs do not contain free calcium and magnesium ions. The script uses [CaHCO₃⁺]/[HCO₃⁻] and [MgHCO₃⁺]/[HCO₃⁻]. By mass action, these ratios are proportional to the free ions, with activity coefficients taken as equal across cells. The old letter mapped the three conditions in the old runs (median rank 99.1% for carbonate ion). The corrected runs replace this mapping.

**Letter text:**

> Carbonate forms where three conditions coincide: dissolved calcium or magnesium, carbonate ion, and a pH high enough for carbonate to be stable. We mapped the three conditions at the start of each time interval and the carbonate formed during the interval. In each condition, the new carbonate ranks at 96.5 to 98.7, where 50 means no preference (medians of 50 networks). Cells where all three conditions exceed the network mean contain 94.7 to 99.3% of the new carbonate. These cells represent 4.0 to 7.7% of the network volume. Cells where one condition is high and another is low contain 0.2 to 4.7% of the new carbonate. The co-location of the three conditions explains the low efficiency. The dissolving minerals release calcium and magnesium across the network. During continuous injection, 98.8% of the released calcium, magnesium and iron is not captured in carbonate. The injected water supplies the carbon. A high pH is rare during injection. From five weeks on, 0.7 to 3.3% of the pore space has a pH above 6. The three conditions coincide in this small part of the network. After injection stops, the minerals neutralize the pore water. The pH increases to 9.1. The co-located volume extends across the network. The carbonate increases until the dissolved carbon is used up.


### 2.2 Flow-field tests requested by the reviewers

The reviewers asked for a test with a metric from the flow field alone (AE-1, R2-5, R3-2). The flow field carries no geochemical information. The tests below give statistical associations between the flow field and the carbonate. The tests do not identify the mineralization condition.

**Index test specified in advance.** The flow index is the age of the oldest 1% of the water divided by the median age. We defined the index and the test before running the 20 networks of Block B.

| Test | Block B (20 networks, test set up in advance) | Block A (50 networks) | Source |
|---|---|---|---|
| Index vs efficiency | ρ = 0.844, p = 1.5 × 10⁻⁶ (one-sided) | ρ = 0.767, p = 8.2 × 10⁻¹¹ | `src/flow_index_test.py` → `results/flow_index_test.csv` |
| Index vs carbonate per cell (the measure named in advance) | ρ = 0.744, p = 8.4 × 10⁻⁵ (one-sided) | ρ = 0.812, p = 8.3 × 10⁻¹³ | same |
| Index vs fracture intensity | ρ = −0.28, p = 0.24 | ρ = 0.03, p = 0.86 | same |
| Old runs | ρ = 0.398, p = 0.041 | | current letter, R2-5 |

**Carbonate and water age.** The 10% of the pore space with the oldest water contains a median 93.8% of the carbonate (Block A, range 51.8 to 99.9%, 50 networks). In Block B, the share is 92.7% (range 50.5 to 99.2%). Source: `src/carbonate_water_age.py` → `results/carbonate_water_age.csv`.

**Checks of these associations.**

| Check | Old runs (Block A, 20) | Corrected, Block B (20) | Corrected, Block A (50) |
|---|---|---|---|
| Share of the flow in the direction of the injection-driven flow (random directions give 50%) | 53.1% (32.2 to 78.5%) | 99.6% (98.9 to 99.8%) | 99.6% (97.8 to 99.8%) |
| Total flow ÷ injection-driven flow | 2,605 (171 to 4,960) | 2.08 (0.58 to 3.12) | 2.04 (0.32 to 3.55) |
| Where carbonate sits in the age range, comparing cells at the same distance from the inlet | | 95.4 | 94.5 |
| Where carbonate sits along the distance from the inlet (50 = no preference) | | 58.9 | 61.7 |
| Index vs efficiency, after correcting for pore volume and cell number | | ρ = 0.840 (p = 1.3 × 10⁻⁵) | ρ = 0.767 (p = 2 × 10⁻¹⁰) |

Source: `src/flow_field_validation.py`. Corrected runs: `results/flow_field_validation.log`. Old runs: hpc02 `runs/A_p32_*`, `results/flow_field_validation_circulating.csv`. The old runs had the wrong flow: 53.1% of their flow followed the direction of the injection-driven flow, and random directions give 50%.

**Letter text (R2-5, R3-2, AE-1):**

> We tested a metric from the flow field alone, as requested. We defined the index and the test before running the 20 networks of Block B. In these networks, the index correlates with the efficiency (ρ = 0.84, one-sided p = 1.5 × 10⁻⁶). The flow field carries no geochemical information. The correlation is a statistical association and does not identify the mineralization condition. The co-location mapping of calcium or magnesium, carbonate ion and pH identifies this condition.

### 2.3 The model rock has a uniform composition

All model cells start with the same minerals and the same reactive surface area. The rock's reactive surface area sets the amount of carbonate formed: ten times less reactive surface leaves 8% of the carbonate (section 5.1). In real basalt, the minerals vary from place to place. The places where the three conditions coincide then depend on the local minerals.

**Letter text (limitation):**

> All cells in the model start with the same minerals. The reactive surface area of the rock sets the amount of carbonate formed. A tenfold reduction of this area reduces the carbonate by 92%. In real basalt, the minerals vary from place to place. The places where the three conditions coincide then depend on the local minerals.

### 2.4 Letter statements to replace

- Summary (ii): "Mineralization does not correlate with the flow field (AE)." Replace with the result of the test specified in advance, reported as an association (section 2.2).
- R2-5: "treats the 16% as the upper limit of prediction from the flow field".
- AE-1: the co-location ranks of the three conditions from the old runs (repeat on the corrected runs, section 2.1).
- R1-6 (letter line 642): "the least-flushed portion of cells contains no carbonate at all, whereas the most-flushed portion contains 99.6%". The corrected runs give 93.8% of the carbonate in the 10% of pore space with the oldest water.
- The figure caption for `fig_cations_pair.pdf`: "Delivered carbon and cation supply are not sufficient without co-location in the same cell." Keep the statement and redraw the figure from the corrected runs.

## 3. Different minerals dissolve at different times, and different carbonates form (AE-2, R2-6, R1-3, R1-9)

The rock does not dissolve all at once. Forsterite dissolves first, within two years. Diopside and anorthite dissolve over decades. Each mineral releases a different element. Magnesium from forsterite forms magnesite. Calcium from diopside and anorthite forms calcite. The carbonate that forms changes over time with the element supplied.

### 3.1 Four stages

| Period | Mineral dissolving | Element released | Carbonate forming | Mean pH |
|---|---|---|---|---|
| First 5 weeks | forsterite starts (6% dissolved) | magnesium, in small amounts | small amounts (0.015 kg CO₂) | 8.5 → 4.9 |
| 5 weeks to 2 years | forsterite, 6 → 95% dissolved | magnesium | magnesite, from 8 to 49% of the carbonate | 4.9 → 4.5 |
| 2 to 20 years | diopside 8 → 85%, anorthite 3 → 41% | calcium and magnesium | calcite share increases from 51 to 74% | 4.5 → 4.3 |
| 20 to 50 years | anorthite 41 → 96%, dissolving faster as the pH decreases | calcium | calcite share increases from 74 to 85% | 4.3 → 3.8 |

Fayalite dissolves by 1.0% and albite by 0.2% or less. Iron and sodium stay scarce. Siderite and dawsonite do not form.

Source: coupled model (`results/A_coupled_timeseries_by_intensity.csv`). Carbonate and pH at P₃₂ × 1.00.

The order of dissolution does not depend on fracture intensity. At each time, the five intensities agree within 1.9 percentage points:

| Dissolved (%), coupled, range over the five intensities | 5 weeks | 1 y | 2 y | 5 y | 10 y | 20 y | 50 y |
|---|---|---|---|---|---|---|---|
| Forsterite | 6.0 to 6.1 | 66.0 to 67.1 | 95.1 to 96.2 | 98.7 to 99.6 | 99.4 to 99.9 | 99.5 to 99.9 | 99.6 to 99.9 |
| Diopside | 0.3 | 3.1 to 3.2 | 7.6 to 7.8 | 22.9 to 23.5 | 48.3 to 49.1 | 84.6 to 86.5 | 96.6 to 98.2 |
| Anorthite | 0.1 | 1.6 | 3.4 | 9.0 | 18.6 to 18.7 | 40.9 to 41.1 | 95.0 to 95.8 |

### 3.2 The amount of each carbonate over time

Coupled Block A, P₃₂ × 1.00, CO₂ in each carbonate (kg).

| | 5 weeks | 6 months | 1 y | 2 y | 5 y | 10 y | 20 y | 50 y |
|---|---|---|---|---|---|---|---|---|
| Magnesite | 0.0013 | 0.097 | 0.46 | 1.06 | 2.43 | 3.40 | 3.98 | 4.99 |
| Calcite | 0.012 | 0.23 | 0.52 | 1.08 | 2.43 | 4.49 | 9.32 | 22.7 |

Magnesite forms fastest in the first 5 years, while forsterite dissolves. Afterwards it grows slowly, from 3.4 to 5.0 kg over 40 years, with magnesium from diopside. Calcite grows fastest in the last 30 years, from 9.3 to 22.7 kg, as anorthite dissolution speeds up. Source: `results/A_coupled_timeseries_by_intensity.csv`.

### 3.3 Injection duration changes the carbonate mix in the same way

Network s1181, P₃₂ × 1.00, coupled model, values at 50 years.

| Injection | Forsterite dissolved | Diopside dissolved | Anorthite dissolved | Calcium share of the released Ca + Mg | Calcite share of the carbonate |
|---|---|---|---|---|---|
| 1 day | 3.2% | −0.8% | 1.2% | 28% | 26% |
| 10 days | 50.8% | −4.7% | 8.9% | 14% | 14% |
| 30 days | 61.6% | −5.2% | 10.1% | 14% | 14% |
| 45 days | 65.4% | −5.1% | 10.4% | 14% | 14% |
| 2 years | 99.2% | 9.1% | 22.9% | 38% | 86% |
| 5 years | 99.4% | 25.5% | 29.0% | 44% | 88% |
| 10 years | 99.4% | 51.6% | 39.7% | 49% | 87% |
| 50 years, no stop | 99.5% | 98.0% | 95.9% | 56% | 91% |

Source: `src/revision_statistics.py` → `results/injection_duration.csv`.network()` on the `F_*`, `E_*` and `A_*` runs of `p32_100_s1181`, dissolved moles from the results rows. Released calcium: anorthite + diopside. Released magnesium: 2 × forsterite + diopside + enstatite.

Two patterns follow from the table:

- **Short injections (45 days or less):** the carbonate mix equals the mix of released elements (calcite 26% against calcium 28%, and 14% against 14%). The released calcium and magnesium turn into carbonate in proportion. Forsterite supplies most of the magnesium.
- **Long injections (2 years or more):** calcite represents 86 to 91% of the carbonate, and calcium represents 38 to 56% of the released elements. Calcite forms in preference to magnesite. Most of the released magnesium does not form carbonate and leaves the network in solution.

Negative values mean that the mineral grew. Section 5.4 describes this effect.

**Letter text:**

> The injection duration changes the carbonate mix. After injections of 45 days or less, the carbonate mix matches the mix of released elements. Calcite represents 14 to 26% of the carbonate. Calcium represents 14 to 28% of the released calcium and magnesium. After injections of 2 years or more, calcite represents 86 to 91% of the carbonate. Calcium represents 38 to 56% of the released elements. Calcite forms in preference to magnesite during long injections. Most of the released magnesium leaves the network in solution.

### 3.4 Dissolution happens throughout the network, carbonate in the oldest water

AE-2 asks for dissolution zones, flow paths and precipitation shown together. The rock dissolves throughout each network, with a slight preference for younger water. Carbonate forms in the oldest water.

| Time | Position of dissolution in the water-age range | Position of new carbonate | Networks |
|---|---|---|---|
| 2 years | 43.3 (41.8 to 44.2) | 95.3 | 50 |
| 10 years | 41.6 (40.3 to 42.4) | 95.9 | 50 |
| 20 years | 44.1 (42.5 to 45.2) | 95.6 | 50 |
| 50 years | 48.7 (47.6 to 49.3) | 95.1 | 50 |

Scale: 50 means no preference, low values mean young water, high values mean old water. Dissolution is the volume of forsterite, diopside and anorthite dissolved in each cell since the start. Source: `src/dissolution_capture_analysis.py` (part 1) → `results/dissolution_capture_analysis.csv`, coupled Block A.

The letter's claim (iii) describes a dissolution front with carbonate forming downstream. The corrected runs show dissolution across the whole network and a weak downstream preference of the carbonate (59 to 62 on the distance scale, section 2.2).

**Letter text (AE-2):**

> The rock dissolves throughout each network. On a scale of water age where 50 means no preference, the dissolved volume sits at 42 to 44 between 2 and 20 years. Dissolution shows a slight preference for younger water. The carbonate sits at 95 to 96 on the same scale. The minerals release calcium and magnesium across the network. Carbonate forms in the slowly renewed water, where the pH is high enough.

### 3.5 pH and the space where carbonate forms

Coupled Block A, P₃₂ × 1.00.

| | 4 days | 5 weeks | 1 y | 2 y | 5 y | 10 y | 20 y | 50 y |
|---|---|---|---|---|---|---|---|---|
| Mean pH | 8.54 | 4.88 | 4.73 | 4.51 | 4.45 | 4.40 | 4.26 | 3.75 |
| Pore space with pH above 6 (%) | 81.8 | 3.3 | 1.6 | 1.7 | 1.2 | 0.8 | 0.7 | 1.1 |
| CO₂ in carbonate (kg) | 0.0018 | 0.0152 | 0.983 | 2.14 | 4.86 | 7.9 | 13.2 | 27.9 |
| Magnesite share (%) | 0.0 | 7.9 | 46.1 | 49.4 | 47.7 | 38.8 | 25.5 | 14.5 |
| Calcite share (%) | 63.9 | 82.4 | 53.9 | 50.5 | 52.2 | 61.2 | 74.4 | 85.3 |

Source: `src/time_series.py` → `results/A_coupled_timeseries_by_intensity.csv`. At every intensity, magnesite peaks at 38 to 52% at 2 years and calcite reaches 84 to 87% at 50 years. No cell contains more mineral than its volume at any time.

### 3.6 Carbonate that dissolves again (R1-9)

In the fixed model, one of 50 networks lost more than 5% of its peak carbonate (9%). In the old runs, five networks lost 52 to 100%. In the coupled model, the median magnesite amount increases at every intensity through the 50 years. The network-by-network check needs running on the coupled model (section 9).

AE-2 also asks how networks with negligible or fully dissolved carbonate enter the analysis. In the corrected coupled runs, each network forms carbonate.

| Set | Lowest carbonate (kg CO₂) | Median (kg CO₂) | Networks below 1% of the median |
|---|---|---|---|
| Block A | 5.09 (s1481 × 0.75) | 28.3 | 0 |
| Block B | 4.81 (s1063 × 0.90) | 24.3 | 0 |
| Shut-in (E) | 33.6 (s941 × 0.75) | 213 | 0 |
| 30 m | 57.4 | 106 | 0 |
| 40 m | 168 | 233 | 0 |

Source: `src/dissolution_capture_analysis.py` (part 4).

**Letter text (AE-2):**

> Each network forms carbonate. In Block A, the lowest amount is 5.09 kg of CO₂, 18% of the median of 28.3 kg. Each network enters each analysis.

**Letter text (dissolution and carbonate sequence):**

> The minerals that dissolve change over time. The carbonates that form change with them. Forsterite dissolves first. By 2 years, 95% of the forsterite has dissolved. The released magnesium forms magnesite. Magnesite represents 49% of the carbonate at 2 years. Diopside and anorthite dissolve over the following decades. By 20 years, 85% of the diopside and 41% of the anorthite have dissolved. Anorthite dissolution speeds up as the pH decreases and reaches 96% at 50 years. The released calcium forms calcite. Calcite represents 85% of the carbonate at 50 years. Fayalite and albite dissolve by 1% or less. Siderite and dawsonite do not form.

**Letter text (pH):**

> The pH changes in three stages. The pH increases to 8.5 in the first four days. The injected acidic water then reaches the whole network. The pH decreases to 4.9 by five weeks. The pH then decreases slowly, to 4.4 at 10 years and 3.8 at 50 years (P₃₂ × 1.00). From five weeks on, 0.7 to 3.3% of the pore space has a pH above 6. Carbonate forms in this space throughout the 50 years.

---

## 4. Number of networks, fracture intensity and domain size (AE-3, R2-3, R3-3, R3-7, R1-12)

### 4.1 Coupled Block A, 10 networks per intensity

| P₃₂ | CO₂ in carbonate, median (kg) | Spread between networks, CV (%) | Efficiency, median (%) | Calcite / magnesite (%) | pH at 50 y |
|---|---|---|---|---|---|
| ×0.75 | 9.85 | 69 | 0.0057 | 86 / 14 | 3.79 |
| ×1.00 | 27.9 | 53 | 0.0055 | 85 / 14 | 3.75 |
| ×1.25 | 25.6 | 38 | 0.0039 | 86 / 14 | 3.71 |
| ×1.50 | 46.7 | 51 | 0.0054 | 87 / 12 | 3.72 |
| ×2.00 | 58.4 | 35 | 0.0056 | 84 / 16 | 3.71 |

Source: `src/summarize_block.py` → `results/A_coupled_networks.csv`. Block B: ×0.90 gives 16.9 kg (efficiency 0.0036%), ×1.75 gives 36.4 kg (0.0038%).

### 4.2 Denser networks form more carbonate because they contain more pore space

| Correlation with fracture intensity | Block A (50) | Block B (20) | Source |
|---|---|---|---|
| Carbonate | 0.66 (p = 1.8 × 10⁻⁷) | 0.62 (p = 0.0033) | `summarize_block.py` |
| Pore space | 0.97 (p = 7 × 10⁻³²) | 0.87 (p = 8 × 10⁻⁷) | `src/revision_statistics.py` → `results/volume_effect.csv` |
| Efficiency | −0.06 (p = 0.67) | 0.03 (p = 0.88) | `src/revision_statistics.py` → `results/volume_effect.csv` |

**Letter text:**

> Denser fracture networks form more carbonate because they contain more pore space. The carbonate increases with fracture intensity (ρ = 0.66, p = 2 × 10⁻⁷, 50 networks). The pore space increases with fracture intensity in the same way (ρ = 0.97). The injection rate scales with the pore space. The efficiency does not change with fracture intensity (ρ = −0.06, p = 0.67).

### 4.3 Domain size

Networks at 40 m vary less from one to another than networks at 20 or 30 m. The networks at each size are different networks. They share seed numbers but not fractures. The comparison shows the variation between similar networks at each size. It does not show how one network changes with size.

| Domain | Networks | Efficiency, median (%) | Spread, CV (%) (95% range) |
|---|---|---|---|
| 20 m | 10 | 0.0055 | 51 (27 to 72) |
| 30 m | 8 | 0.0056 | 71 (26 to 87) |
| 40 m | 8 | 0.0054 | 19 (10 to 24) |

| Comparison | Ratio of spreads (95% range) | Probability of a ratio this low by chance (permutation test) |
|---|---|---|
| 40 m vs 20 m | 0.37 (0.18 to 0.68) | 0.012 |
| 40 m vs 30 m | 0.26 (0.14 to 0.71) | 0.012 |
| 30 m vs 20 m | 1.38 (0.48 to 2.56) | 0.69 |

Source: `src/revision_statistics.py` → `results/domain_size_variability.csv`. The spread does not decrease steadily (51, 71, 19%), and each size has 8 to 10 networks. An earlier paired analysis treated equal seeds as the same network and is withdrawn.

**Letter text:**

> The efficiency is similar at 20, 30 and 40 m (median 0.0055, 0.0056 and 0.0054%). Networks at 40 m vary less from one to another than networks at 20 or 30 m (coefficient of variation 19%, against 51% and 71%, permutation p = 0.012). The networks at each size are different networks with the same fracture statistics. The comparison shows the variation between similar networks at each size. Showing how one network changes with size requires smaller domains cut from the same network. We did not simulate nested domains.

### 4.4 Letter statements to replace

- AE-3: "The pH at 50 years remains consistent across all five intensity levels (5.97 to 6.10)." The corrected pH is 3.71 to 3.79.
- Summary (v): the spread trend 147 → 112 → 99% for 20, 30 and 40 m.

---

## 5. Model assumptions (AE-4, R1-1, R1-2a to c, R1-7, R1-8, R2-1, R2-2, R2-4, R3-4, R3-5, R3-6)

### 5.1 Which parameters matter

The reactive surface area of the rock changes the carbonate more than any other parameter tested. Dawsonite, analcime and the seed amount do not change it.

| Change | Networks | Carbonate compared with the baseline | Anorthite dissolved (baseline) |
|---|---|---|---|
| All reactive surface areas ×0.1 | 2 | 0.08 (0.07 to 0.09) | 17% (96.5%) |
| All reactive surface areas ×10 | 1 | 2.02 | 99% (95.8%) |
| Anorthite surface area 30 / 50 / 100 | 2 each | 1.11 / 1.14 / 1.19 | 99.7% (96.5%) |
| Surface area of new minerals, higher / lower | 2 each | 1.10 / 0.94 | 96.5% |
| Dawsonite rate 10⁻⁹ to 10⁻¹³ (baseline 10⁻⁷ mol m⁻² s⁻¹) | 3 each | 1.00 | 96.0% |
| No dawsonite | 5 | 1.00 | 95.8% |
| Analcime allowed, rate 10⁻⁹ to 10⁻¹³ | 2 each | 1.00 | 96.5% |
| Seed amount higher / lower | 2 each | 1.01 / 1.00 | 96.5% |
| `vf_consistent` (same as the baseline since the correction) | 2 | 1.00 | drop from Table R5 |

Source: `src/revision_statistics.py` → `results/sensitivity_ratios.csv`. We checked each variant deck against the baseline deck. Each differs in the intended lines. Analcime at the fastest rate reaches a mean volume fraction of 8 × 10⁻¹⁰. The anorthite surface-area values (baseline 10, variants 30 to 100 cm² cm⁻³) need confirming in the decks.

**Letter text:**

> The reactive surface area of the rock changes the carbonate more than any other parameter tested. A tenfold smaller surface area reduces the carbonate by 92%. A tenfold increase doubles it. Changing the anorthite surface area alone changes the carbonate by 11 to 19%. Changing the surface area of the new minerals changes it by 6 to 10%.

> Dawsonite and analcime need sodium. Albite supplies the sodium and dissolves by 0.2% or less. Dawsonite rate constants from 10⁻⁷ to 10⁻¹³ mol m⁻² s⁻¹ do not change the carbonate. Removing dawsonite does not change the carbonate. At the fastest rate tested, analcime reaches a mean volume fraction below 10⁻⁹. Siderite does not form. Its iron source, fayalite, dissolves by 1.0%.

### 5.2 Shrinking mineral surfaces reduce the carbonate by 27 to 50%

We let each mineral's surface shrink as the mineral dissolves and as pores fill (hpc01, 10 networks).

| P₃₂ | Constant surface (kg) | Shrinking surface (kg) | Ratio | Anorthite dissolved | pH at 50 y |
|---|---|---|---|---|---|
| ×0.75 | 11.2 | 5.6 | 0.50 | 95.5 → 78.6% | 3.75 → 3.90 |
| ×1.00 | 17.0 | 12.3 | 0.72 | 95.6 → 78.9% | 3.73 → 3.92 |
| ×1.25 | 27.8 | 19.4 | 0.70 | 95.4 → 78.9% | 3.74 → 3.92 |
| ×1.50 | 55.1 | 39.3 | 0.71 | 94.4 → 78.2% | 3.78 → 3.94 |
| ×2.00 | 64.6 | 47.0 | 0.73 | 95.1 → 78.8% | 3.75 → 3.92 |

Source: `summarize_block.py` on hpc01 → `archive_hpc01/results/S0_*.csv`, `S1_*.csv`. Settings: `src/add_surface_area.py`. The ratios compare medians of two networks. The network-by-network ratios need calculating (section 9).

**Letter text:**

> Letting mineral surfaces shrink as minerals dissolve and pores fill reduces the carbonate by 27 to 50% (ten networks at five intensities). Anorthite dissolution decreases from 95% to 79%. The dependence on fracture intensity and the mineral mix stay the same. The association between the flow index and the efficiency stays positive with shrinking surfaces (ρ = 0.64, p = 0.048, ten networks).

### 5.3 Dissolution opens the fractures (R2-4, R3-4)

| Measure | Value | Source |
|---|---|---|
| Porosity at 50 y, median per network | 0.70 to 0.72 (start 0.50), maximum 0.84 | inline on Setonix, part of Block A |
| Permeability at 50 y ÷ start | median 2.8 to 3.0, maximum 4.9 | same |
| Injection pressure, early ÷ at 50 y | median 3.2 (2.1 to 4.1), 32 networks | same |
| Fixed model, cells with more carbonate than their volume | 45 of 49 networks, up to 145 times | `check_block.py` |

These numbers come from Setonix before Block A finished. They need recalculating on the full set (section 9).

**Letter text:**

> Dissolution opens the fractures. The porosity increases from 0.50 to 0.70 by 50 years (median). The permeability increases 2.8 times. The pressure needed to keep injecting decreases 3.2 times. Cells filled with carbonate clog. The injected water flows around them.

### 5.4 Diopside and albite grow in the model after injection stops (limitation)

After short injections, diopside grows by up to 5.2% of its initial volume (section 3.3). Albite grows by 0.1 to 0.5% in several sets. The pH increases to 9 to 10 after the injection stops, and the pore water becomes supersaturated with these minerals. The rate law in PFLOTRAN acts in both directions, and supersaturated minerals grow. Diopside and albite form at igneous temperatures and do not form at the model temperature. The growth takes up calcium and magnesium. Without the growth, these elements form carbonate or remain in solution.

PFLOTRAN v6 lists a separate `PRECIPITATION_RATE_CONSTANT` for a mineral, paired with a dissolution constant. Whether it combines with the `PREFACTOR` mechanisms in the decks is not known. A test on the 1-day and 45-day injections, with a negligible precipitation constant for diopside and albite, measures the effect. These runs are short.

**Letter text (limitation):**

> After injection stops, the pore water becomes supersaturated with diopside and albite. The rate law in the model allows these minerals to grow. Diopside grows by up to 5.2% of its initial volume after injections of 45 days or less. Diopside and albite do not form at the model temperature. The carbonate formed after short injections is a lower estimate.

### 5.5 Kaolinite fills more pore space than the carbonate

Kaolinite and chalcedony take up the aluminium and silica from the dissolving silicates. Kaolinite forms far more volume than the carbonate.

| Set | Kaolinite (mol, median) | Chalcedony (mol, median) | Kaolinite and chalcedony, share of the initial pore volume | Carbonate (mol CO₂, median) |
|---|---|---|---|---|
| Block A, continuous | 11,630 | 28.65 | 24.1% (23.4 to 24.7%) | 642 |
| Shut-in (E) | 6,563 | 735.5 | 14.1% (13.5 to 14.3%) | 4,840 |

Source: `src/dissolution_capture_analysis.py` (part 3), net change from PFLOTRAN's mass balance. Carbonate from `results/A_coupled_networks.csv` and `results/E_coupled_networks.csv` (median CO₂ ÷ 0.0440 kg mol⁻¹).

Carbonate closes the pores where clogging occurs, and clogging is rare. At 50 years, 222 of 6.4 million cells (0.0035%) have a porosity of 0.02 or less, across the 50 networks of Block A. In these cells, the new mineral volume contains 66.2% carbonate, 18.0% kaolinite and 15.8% chalcedony. Kaolinite spreads across the network. Carbonate concentrates in the cells with old water and closes the pores there. Source: `src/revision_statistics.py` → `results/clogging_cells.csv`.

**Letter text (R2-2):**

> Kaolinite and chalcedony take up the aluminium and silica released by the dissolving silicates. At 50 years, they fill 24.1% of the initial pore volume (median of 50 networks, range 23.4 to 24.7%). Kaolinite accounts for 11,630 mol, compared with 642 mol of carbonate (medians). Kaolinite fills more pore space than the carbonate across the network. Clogging is rare. At 50 years, 222 of 6.4 million cells have a porosity of 0.02 or less. In these cells, the new mineral volume contains 66.2% carbonate, 18.0% kaolinite and 15.8% chalcedony.

### 5.6 Letter statements to replace

- Summary (iv): the dawsonite rate correction. It is no longer needed.
- Summary (vi): "Aperture evolution reduces mineralization by 31%".
- AE-4: sensitivity numbers from the old runs.
- R1-1 and R2-8: anorthite dissolution of 69.2 to 74.2%. The corrected value is 95.0 to 96.3%.

---

## 6. Injection time and field comparison (AE-5, R1-5, R1-10, R3-8)

### 6.1 Stopping injection increases the carbonate 6.2 times

After injection stops at 10 years, the rock neutralizes the acidic water. The pH increases to 9.1. The carbon dissolved in the water then turns into carbonate.

| Measure (coupled, 50 networks) | Value | Source |
|---|---|---|
| Carbonate, injection stopped ÷ continuous, same network | median 6.24 (3.08 to 30.58) | `src/revision_statistics.py` → `results/shutin_porewater_carbon.csv` |
| Carbonate ÷ carbon dissolved in the water at shut-in | median 1.08 (1.02 to 1.17) | `src/revision_statistics.py` → `results/shutin_porewater_carbon.csv` |
| pH at 50 y | 9.05 to 9.09 | `results/E_coupled_*.csv` |
| Anorthite / diopside dissolved | 39 to 40% / 51% | same |
| Old runs | 1.1 times (45 pairs) | current letter |

**Letter text:**

> Stopping the injection increases the carbonate 6.2 times (median of 50 networks, range 3.1 to 30.6). After injection stops at 10 years, the rock neutralizes the acidic water. The pH increases to 9.1 by 50 years. The carbon dissolved in the water turns into carbonate. By 50 years, the carbonate equals 108% of the carbon dissolved at shut-in (median, range 102 to 117%).

### 6.2 Short injections convert a higher share of their carbon

All values come from PFLOTRAN's mass balance (network s1181, P₃₂ × 1.00).

| Injection | Carbon injected (mol) | Carbonate formed (mol) | Efficiency |
|---|---|---|---|
| 1 day | 72.9 | 71.2 | 97.6% |
| 10 days | 5,076 | 2,947 | 58.1% |
| 30 days | 20,202 | 3,378 | 16.7% |
| 45 days | 31,463 | 3,444 | 11.0% |
| 2 years | 5.50 × 10⁵ | 3,715 | 0.68% |
| 5 years | 1.38 × 10⁶ | 3,921 | 0.28% |
| 10 years | 2.76 × 10⁶ | 4,270 | 0.15% |
| 50 years, no stop | 1.38 × 10⁷ | 354 | 0.0026% |

Source: `src/revision_statistics.py` → `results/injection_duration.csv`. The network's pore space contains 4,258 kg of water. Filled with injected water, this water contains 3,492 mol of carbon. After a stop, the carbonate formed approaches this amount.

**Letter text:**

> Short injections convert a higher share of their carbon into carbonate. A one-day injection converts 97.6%. A 45-day injection converts 11.0%. Continuous injection over 50 years converts 0.0026%. After injection stops, the carbonate formed approaches the carbon in one pore volume of injected water. Carbon injected beyond one pore volume leaves the network before injection stops. The efficiency after a stop is close to the pore volume divided by the injected volume.

> At CarbFix, the pore volume reached by the injected water exceeded the injected volume. The short injections in our simulations reach efficiencies close to the CarbFix value of 95%.

### 6.3 Calcium and magnesium never limit the carbonate

| Set | Share of the released Ca, Mg and Fe in carbonate, median (range) | Networks |
|---|---|---|
| Block A, continuous injection | 1.25% (0.23 to 2.7%) | 50 |
| Block B, continuous injection | 0.89% (0.24 to 4.2%) | 20 |
| Shut-in at 10 years (E) | 13.7% (12.8 to 15.0%) | 50 |

Released elements in moles: calcium from anorthite and diopside, magnesium from forsterite (×2), diopside and enstatite, iron from fayalite (×2). Source: `src/dissolution_capture_analysis.py` (part 2), from `results/*_coupled_networks.csv`.

During continuous injection, 98.8% of the released elements leave the network or stay dissolved. The low pH limits the carbonate. After shut-in, the carbonate equals 108% of the dissolved carbon (section 6.1), and 86.3% of the released elements stay in excess. The dissolved carbon limits the carbonate.

**Letter text (R1-10, R3-8, AE-5):**

> Calcium and magnesium do not limit the carbonate. During continuous injection, 1.25% of the calcium, magnesium and iron released by dissolution forms carbonate (median of 50 networks). The rest leaves the network or stays dissolved. The low pH limits the carbonate during injection. After injection stops, 13.7% of the released elements forms carbonate. The carbonate then equals 108% of the carbon dissolved at shut-in. The dissolved carbon limits the carbonate after injection stops.

### 6.4 Letter statements to replace

- Summary: the efficiency peak of 0.0020% at 3 to 5 years, and the gap of five orders of magnitude with field figures.
- AE-5: "The simulated mineralization efficiency peaks at 0.0020%".
- The statements that calcium and magnesium supply limits the carbonate, and that dilution by the injected water is excluded.
- Summary (iii): "Precipitation is spatially decoupled from dissolution", occurring downstream of the dissolution front. Replace after the test in section 3.4.
- The cube-root scaling of carbonate with injected carbon (R1-4, Summary).
- Table R9, and the shut-in ratio of 1.1.

---

## 7. Other comments

| Comment | Result |
|---|---|
| R1-3, R2-7 (pH) | Section 3.1. The injected water has pH 3.4 |
| R2-8 (anorthite and fracture intensity) | Anorthite dissolves 95.0 to 96.3% at every intensity |
| R1-6, R3-9 (dead ends, diffusion) | No pore space without flow. Carbonate forms in slow-flowing water (section 2.2). Diffusion into the rock matrix is not simulated |
| R1-14, R3-1 (glass, brine) | No new results. Update the numbers cited |
| R1-4 (contributions) | Section 1 |

---

## 8. Methods statements

| Topic | Text | Source |
|---|---|---|
| Gravity | We switched gravity off. With gravity on, the uniform-pressure outflow boundary made water circulate 400 to 6,600 times faster than the injection | `src/fix_gravity_off.py` |
| Porosity limit | We set a lower porosity limit of 0.01. Without the limit, 13 networks stopped when carbonate filled a cell. With the limit, each of the 13 completed | `src/add_coupling.py`, `lists/failures.tsv` |
| Effect of the limit | We reran the 8 completed networks that reached the limit. The limit changed their carbonate by less than 0.2% (median −0.01%) | `src/revision_statistics.py` → `results/porosity_limit_pairs.csv` |
| Computers | We ran the simulations on three computers. Identical runs on different computers agreed within 0.06% (29 pairs) | `src/revision_statistics.py` → `results/platform_comparison.csv` |
| Carbon accounting | We took the injected carbon and the carbonate from PFLOTRAN's mass balance | `src/carbon_budget.py` |
| Continued runs | Six runs continued from saved states. We combined their output files | `src/make_restart.py`, `src/carbon_budget.py` |
| Unfinished runs | Coupled: one network (s42 at ×1.50 with ten times the surface area) reached 0.17 years. Fixed: B s941 at ×1.75 (1.1 years) and the same s42 case (0.2 years) | `lists/failures.tsv` |
| Replaced network | Network s941 at ×2.00 has no connected path across the domain. Network s117 replaces it | Block A definition |

---

## 9. Remaining work

**Scripts for the values.** `src/revision_statistics.py` computes the values from session commands. Each section writes `results/<section>.csv`: `volume_effect`, `domain_size_variability`, `sensitivity_ratios`, `injection_duration`, `shutin_porewater_carbon`, `porosity_limit_pairs`, `platform_comparison`, `flow_index_surface_area`, `clogging_cells`, `shrinking_surface_ratios`, `anorthite_surface_area`, `silicate_regrowth`. The other scripts: `src/flow_index_test.py`, `src/carbonate_water_age.py`, `src/flow_field_validation.py`, `src/colocation_mapping.py`, `src/dissolution_capture_analysis.py`, `src/plot_mineralization_figures.py`.

**Run on the complete coupled data (hpc02):**

1. Done: coupled dissolution over time (3.1) and in the injection-duration runs (3.3).
2a. Measure diopside and albite growth after shut-in in the coupled E ensemble, and test a negligible precipitation constant on the 1-day and 45-day injections (5.4).
2. Carbonate that dissolves again, network by network (3.6).
3. Porosity, permeability and injection pressure on all 50 Block A networks (5.3).
4. Network-by-network ratios for shrinking surfaces (5.2).
5. Anorthite surface-area values in the Block C decks (5.1).
6. Whether the test plan written in advance names the coupled or the fixed model.
7. Rerun `src/flow_field_validation.py` on the corrected runs and save the output as `results/flow_field_validation_corrected.csv`.
8. Run `src/colocation_mapping.py` on coupled Block A: the direct test of the three-condition mechanism (2.1).
9. Done: `src/dissolution_capture_analysis.py` (sections 3.4, 3.6, 5.5 and 6.3).
10. Done: the clogging cells contain 66.2% carbonate by new mineral volume (5.5).
11. Done: co-location mapping on the corrected runs (2.1).

**Redraw the figures from the corrected runs.** The letter uses three figures made from the old runs:

| Old figure | Used in the letter | New figure (`src/plot_mineralization_figures.py`) | Content |
|---|---|---|---|
| `figures/fig_cations_pair.pdf` | AE-1, R1-11, R2-6, R3-2 | `figures/fig_mineralization_colocation.pdf` | two networks at P₃₂ × 1.00 (highest and lowest efficiency): flow, dissolution, Ca or Mg, carbonate ion, pH and new carbonate, with the co-located cells outlined |
| `figures/letter_variability.pdf` | AE-3, R1-12 | `figures/fig_efficiency_intensity_domain.pdf` | carbonate and efficiency of each network against fracture intensity, and efficiency at 20, 30 and 40 m |
| `figures/letter_sensitivity.pdf` | AE-4 | `figures/fig_parameter_sensitivity.pdf` | carbonate of each Block C variant relative to the baseline, by network |

Update the `\includegraphics` lines in the letter to the new file names. `figures/fig_study_design.pdf` shows the study design and does not depend on the results.

**Letter.** Run the prose checklist on each response once it is in the letter.
