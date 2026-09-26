# Regression suite for hot-companion gravity

Every test the law has passed or failed since round 3, in one place, scored the same way, so
that any change to the law, its constants or the code can be checked against everything at once.

```bash
cd research_work/results/hot-companion/regression
python run_suite.py                               # the adopted law (round 12), quick tier (about 1 minute)
python run_suite.py --tier full                   # + the colliding clusters (about 20 minutes)
python run_suite.py --law no_hold                 # a candidate law from candidates/
python run_suite.py --law round11                 # round 11's constants and distance scale
python run_suite.py --law round9                  # the same law with the round-3 constants (before round 11)
python run_suite.py --law round3                  # the law before round 9 (released at once)
python run_suite.py --only dwarfs,precision       # some groups only
python run_suite.py --tier full --save-baseline   # make this run the new baseline
```

The runner prints a scoreboard and writes `runs/<law>-<tier>/results.json` and `report.md`
(both ignored by git). The exit code is 1 when any check got worse in status than in
`baseline.json`, or crashed. Each run records the commit it came from, marked `+changes` when
tracked files differ from it.

**The baseline** (`baseline.json`) is the adopted law, round 12, on the full tier: 90 checks, 59 pass,
11 close, 7 fail and 13 tracked (round 12, step 3: the distance law's scale fitted jointly, α x0.95,
SPARC's Hubble-flow galaxies in the static law too, and the constants refitted there; results README
§22.3). It changes only when a change is adopted or a test is added (RULES.md §12). Earlier baselines:
round 9, 62 / 8 / 6; round 10 (the far collisions and SLACS in the static distances), 58 / 11 / 7;
round 11, 57 / 9 / 10; round 12 step 1 (MACS J0025 at the age its shock fronts give, 0.3 Gyr; §22.1),
58 / 9 / 9; step 2 (the KiDS and Mistele lenses at their measured heat, and the Sérsic gap graded;
§22.2), 59 / 10 / 8.

**The law's distances.** A law carries its distance scale (`alpha_per_Mpc`) and how SPARC is placed
(`sparc_distances`); `common.apply_distances` sets both before the tests run, so a candidate can move the
distance scale like any other constant (`run_suite.py --law round11` runs the previous scale).

## The law under test

`law_config.py` turns a candidate file into one dictionary that every test reads. Nothing in the
suite reads constants from anywhere else. The adopted law (round 12) is the round-3 law with
gradual release over 30,000 AU (0.15 pc), adopted for Cassini in round 9, with every data set in the
project's own (static) distances: the distance scale α = 2.3645 × 10⁻⁴ Mpc⁻¹ (x0.95 of rounds 10–11,
fitted jointly to Pantheon+ supernovae and SPARC's Hubble-flow galaxies), SPARC's Hubble-flow galaxies
placed there too, X-COP's stellar profiles deprojected, and the three constants refitted:
a = 6.298 × 10⁻¹¹ m/s², g_d = 2.027 × 10⁻¹⁰ m/s², u = 169.4 km/s (`run-distance-scale-v12`; results
README §22.3). `--law round11` loads round 11's (a = 6.547 × 10⁻¹¹, g_d = 2.107 × 10⁻¹⁰, u = 162.6 km/s,
α x1, SPARC at its published distances); `--law round9` the round-3 constants (fitted in the release's
ΛCDM units on projected stars). A candidate is a small JSON file in
`candidates/`, a change to a base law:

| Key | Meaning | Default |
|---|---|---|
| `gd_scale` | the release scale g_d multiplied by this | 1 |
| `release_length_au` | the companion is released gradually over this length around each emitter, R(r) = 1 − e^(−r/L) | 30,000 AU (adopted in round 9; `round3` has 0, released at once) |
| `external_hold` | how strongly a subsystem's companion follows an outside galaxy's pull (dwarfs, the Sun, wide binaries) | 1 (fully) |
| `heat_exponent` | the heat weight k = 3 (σ/u)^p (round 14). Every heat weight in the code goes through `law.heat_weight` / `law.k_from_sig2`, set from this key by `common.apply_distances`; at p = 2 the arithmetic is the scripts' own, so the suite reproduces its baseline to the last digit | 2 |
| `hot_geometry`, `stream_kappa_per_Mpc` | how the hot matter's extra glow is heard (round 19; `law.HOT_GEOMETRY`, `law.STREAM_KAPPA`): 'two_way' (the law: every shell, inside and outside the receiver), 'one_way' (inner shells only), 'one_way_vector' (the inner shells' net flux), or 'stream' (a stream absorbing inward-travelling waves at κ per Mpc, which also filters the clusters' cold glow in `run_v3.cluster_M3`). Spherical sums only: the collision maps' 3D sums stay two-way. The stream's weights come from a table of the part it removes (`law._absorbed_table`, cached in `cache/`), exact to 0.3% | 'two_way', 0 |
| `refit` | constants refitted on their home data after the change: `"a"` on the 149 SPARC galaxies (g_d held), `"u"` on the 12 X-COP clusters (for bases 'round11' and 'round12', since round 14, the sample the clusters test grades: stars deprojected, static distances, as those constants were fitted; before, the round-3 sample, which gave u = 192 km/s for the adopted law instead of its 169.4) | none |
| `a_SI`, `lam`, `u_kms`, `base` | explicit constants; `base` 'round3' (the round-3 constants, as the rounds 7–9 candidates use), 'round11', 'round12' (the adopted ones) or another `results.json` | 'round3' |
| `alpha_per_Mpc`, `sparc_distances` | the static distance law's scale, and 'published' or 'static' SPARC distances | the base's ('round12': 2.3645 × 10⁻⁴, 'static'; others 2.4890 × 10⁻⁴, 'published') |
| `distance_variant`, `eta_path` | round 21, registered comparisons: the distance geometry, 'fixed' (the adopted law, D_A = D) or 'metric' (D_A = D/(1 + z), so D_L = (1 + z)² D_A), and a path factor √(1 + η z/(1 + z)) on D (`collisions_v10.VARIANT`, `ETA_PATH`). They reach every conversion that goes through `collisions_v10.static` (KiDS, Mistele, the Bullet, the far collisions); SLACS keeps its own conversion, and SPARC's Hubble-flow distances take α but not the path factor (under 1% at their redshifts) | 'fixed', 0 |
| `collision_field` | round 24, registered comparisons: the field on the collision maps. 'law' (the adopted law's `bullet_v4.kappa_map_v4`) or a mode of `code/eft_field_v24.py`, installed by `common.apply_distances` for every collision code (they all call `bullet_v4.kappa_map_v4` at run time): 'eft' (the Casimir-EFT note's local equation, heat as extra source mass, today's matter), 'eft_memory' (the same with the law's memory), 'eft_excess' / 'eft_excess_local' (the note's source plus the hot glow's non-flowing brightness S − \|g_hot\|, with / without memory), 'eft_scalar' / 'eft_scalar_local' (brightness only), and '_aqual' variants solved exactly (curl field kept; standalone runs only). Mode 'law' reproduces the law's map to 2 × 10⁻¹⁶ | 'law' |

Candidates now in the folder:

| File | What it changes | Why it was proposed |
|---|---|---|
| `gd_x1p5.json` | g_d × 1.5 | the Sun's rotation speed (round 7) |
| `gd_x1p5_refit.json` | g_d × 1.5, then a and u refitted | the same, with the constants re-balanced |
| `weak_hold.json` | external hold 10% | the six faint dwarf galaxies (round 7) |
| `no_hold.json` | external hold 0, release over 200,000 AU (about 1 pc) | the dwarfs, once the release length protects the binaries (round 9) |
| `no_hold_r12.json` | the same on the round-12 law (round 20, a registered comparison, not adopted) | the review's request to rerun it under today's inputs |
| `dist_metric_eta_r12.json` | the supplied distance law: path factor with η = 1/2, metric geometry, α from the supernovae at η = 1/2 (H₀-like 73.2); a and u refitted (round 21) | the owner's feedback: round 12's beam-area term read as a path correction, with distance duality |
| `dist_fixed_eta_r12.json` | the path factor with η = 1/2 in the adopted fixed geometry, same α; a and u refitted | separates the path factor from the geometry |
| `dist_metric_r12.json` | the metric geometry alone at the adopted α; a and u refitted | separates the geometry from the path factor |
| `gd_x1p25.json` | g_d × 1.25 | the smallest useful step toward the Sun's speed |
| `combined.json` | all three, a and u refitted | together |
| `heat_p175.json` | the heat weight k = 3 (σ/u)^1.75, a and u refitted (6.181 × 10⁻¹¹, 132.2 km/s) | the KiDS early/late gap allows p = 1.75–2 and p = 1.75 removes KiDS's level (round 13) |
| `heat_p2_refit.json` | the adopted law through the same refit | the control: reproduces the baseline exactly |
| `eft_r12.json` | the Casimir-EFT note's local equation (round 24): heat as extra source mass (Gauss in round systems, `hot_geometry` 'one_way_vector'; collision maps 'eft'), hold kept | test A of round24-casimir-eft.md |
| `eft_r12_refit.json` | the same, a and u refitted (5.98 × 10⁻¹¹, 98.6 km/s) | whether the constants can absorb the change |
| `eft_nohold_r12.json` | the same without the strong-field hold (lam = 10⁹): the note's equation exactly as written | why the hold is needed |
| `eft_memory_r12.json` | the note's equation with the law's memory on the collision maps | separates memory from how the heat adds up |
| `eft_excess_r12.json` | the proposed repair: the note's source plus the hot glow's non-flowing brightness S_ex = S − \|g_hot\| (round systems: the law exactly), with memory | the repair (round 24, results README §35.3) |
| `eft_excess_local_r12.json` | the repair with no memory (today's matter only) | a purely local field equation |
| `eft_excess_aqual_r12.json` | the repair solved exactly on the collision maps (curl field kept; damped iteration, relaxation 0.3, 60 steps) | whether the field form's shortcut matters |
| `eft_scalar_r12.json`, `eft_scalar_local_r12.json` | brightness only (the heat enters through S, the direction from the cold flow), with / without memory | an alternative repair |

`gradual_release.json` (release over 30,000 AU) was adopted in round 9 and is now the default.

Which tests each amendment reaches: g_d reaches everything; the release length reaches only the
precision tests (below 0.15 pc no galaxy-scale quantity changes); the external hold reaches the
dwarfs, Cassini and the wide binaries.

## How a check is scored

Each check compares one number from the law with one measurement (`checks.py`):

* **pass**: within 2 standard errors of the measurement, or inside the required range;
* **close**: within 3;
* **fail**: further out;
* **info**: tracked but not graded (disputed data, or a number with no measurement);
* **error**: the test crashed.

Every check also carries a **score** (lower is better): its distance from the measurement in
standard errors, or its value over a limit. Against the baseline each check is marked:
* **regressed** / **improved**: its status changed;
* **worse** / **better**: its score moved by more than 0.1 or 5%;
* **known**: a fail already in the baseline;
* **changed**: its number moved by more than 0.1% while its grade and score did not (a tracked
  number, or a value moving inside a passing range);
* unchanged otherwise.

## What is tested (89 checks in the full tier)

| Group | Checks | Data |
|---|---|---|
| machinery | the standing rule (not MOND, Newton or dark matter: the formula guard); the switch-off in strong pulls; the constants and energy budget (tracked) | `research_work/tools/formula_guard.py` |
| galaxies | SPARC typical speed miss (all, test, validation) against MOND's; median residual overall and where the switch-off acts; bulge-dominated galaxies | Lelli et al. 2016 |
| clusters | X-COP mass miss and its trend with radius; the u X-COP prefers; held-out miss over 924 half splits (full tier) | Ettori et al. 2019; Ghizzardi et al. 2021 |
| lensing | KiDS-1000 lensing pull for six samples and the early/late gap; lensing circular speeds of spirals and ellipticals; SLACS light = matter (the project's static distances since round 10); Einstein Cross; microlensing | Brouwer et al. 2021; Mistele et al. 2024; Auger et al. 2009 |
| milky_way | speed at the Sun; 15–27 kpc curve (four Gaia analyses); vertical pull at 1.1 kpc; mass inside 20/50/100/200 kpc; escape speed; inner Galaxy | `data/mw_literature_v7.json` |
| dwarfs | speed spread of ten dwarf spheroidals, with the Galaxy's pull and heat | `data/mw_dwarfs.json` |
| precision | planets, S2, the Double Pulsar, light bending, Cassini's Q2; wide binaries (pass between the two published analyses, 1.0–1.5) | Hees et al. 2014; Chae 2023–24; Banik et al. 2024 |
| collisions (full) | the Bullet Cluster (round-5 case): outer stars, lensing strengths, gas residuals, peak positions, masses inside 250 kpc; the 72-collision stack (β); MACS J0025.4−1222 (lensing inside 300 kpc, peak positions, galaxy speeds); Abell 520 (six clumps inside 150 kpc, 710 kpc, galaxy speeds per clump); El Gordo (aperture lensing masses inside 0.5 and 1 Mpc, galaxy speeds; the SE peak's offset from the cool core tracked); the three in the project's static distances since round 10 | Clowe et al. 2006; Harvey et al. 2015; Bradač et al. 2008; Jee et al. 2014; Clowe et al. 2012; Mahdavi et al. 2007; Girardi et al. 2008; Menanteau et al. 2012; Kim et al. 2021; `code/collisions_v8.py`, `code/collisions_v10.py` |

Every test calls the same code that produced the published numbers (`../code/`); with
`--law round9` the suite reproduces the round-9 numbers for the tests whose conventions did not
change (SPARC 15.85 km/s, collision stack β 0.027). The round-11 baseline gives SPARC 15.94 km/s,
X-COP 0.221 (static distances, deprojected stars), KiDS +0.063 / +0.077 / +0.040 (static
distances), Bullet κ 0.715 and 0.259, collision stack β 0.016.

## Adding a test

Write a function in the right `t_*.py` module that returns `checks.make(...)` objects: an id, a
plain-language title, the law's value, and a criterion from `checks.py` (`z_check`,
`range_check`, `floor_check`, `at_most`, `rms_z`). Then run the full tier and save the baseline.
The Milky Way's Newtonian fields are cached in `cache/` (they do not depend on the law).

## Adding a candidate

Copy a file in `candidates/`, change the keys, and run
`python run_suite.py --law <name> --tier full`. The scoreboard lists what the change fixes
(improved), what it breaks (regressed), and what moved without changing status.

## Round 14: the heat exponent (full tier, against the round-12 baseline)

| Law | Constants | Pass / close / fail | What moves |
|---|---|---|---|
| round 12 (adopted) | a 6.298 × 10⁻¹¹, u 169.4 | 59 / 11 / 7 | |
| `heat_p2_refit` (control) | the same, refitted: a 6.298 × 10⁻¹¹, u 169.44 | 59 / 11 / 7 | nothing: 83 checks the same |
| `heat_p175` | a 6.181 × 10⁻¹¹, u 132.2 | **60 / 10 / 7** | better: KiDS all 0.065 → 0.012 dex (pass), red 0.078 → 0.018 (fail → pass), disks 0.066 → 0.038 (pass), Abell 520 P6 (pass); worse: Mistele's ellipticals rms z 2.6 → 5.2 (fail), SLACS light = matter −0.028 → −0.054 (close), SPARC bulges 29.4 → 30.5 km/s (close; MOND 30.35), Abell 520 P2 (close) |

p = 1.75 trades the KiDS level for the lensing of massive ellipticals: Brouwer et al.'s KiDS relation wants a little
more heat at the red lenses' 140 km/s, Mistele et al.'s circular speeds from the same survey, and SLACS, a little less.
Not adopted; the law keeps p = 2 (results README §24.1).

## Round 19: how the hot glow is heard (quick tier, against the round-12 baseline)

The results README §29.3; the runs are kept in `../run-hot-shell-v19/suite/`.

| Law | Constants | Pass / close / fail | What moves |
|---|---|---|---|
| round 12 (adopted) | a 6.298 × 10⁻¹¹, u 169.4 | 36 / 7 / 6 | |
| `stream_k3` (no refit) | the same | 35 / 8 / 6 | only the X-COP radial trend, 0.239 → 0.285 (close) |
| `stream_k3_refit` | a 6.206 × 10⁻¹¹, u 139.7 | 35 / 8 / 6 | better: KiDS all 0.065 → 0.024, red 0.078 → 0.029 (pass); worse: Mistele's ellipticals 2.55 → 5.15 (fail), SLACS −0.028 → −0.060 (close), SPARC bulges 29.4 → 30.4 km/s (close), X-COP trend 0.251 (close): round 14's trade, since a lower u makes every hot star louder |
| `stream_k10_refit` | a 6.116 × 10⁻¹¹, u 119.9 | 32 / 10 / 7 | the same, further |
| `stream_k30_refit` | a 6.020 × 10⁻¹¹, u 104.4 | 30 / 9 / 10 | the same, further |
| `one_way_refit` | a 6.083 × 10⁻¹¹, u 113.8 | 32 / 6 / 11 | X-COP trend 0.399 (fail), bulges, KiDS bulges and GAMA, Mistele, SLACS fail |

The data want the hot glow heard from all around: an absorbing stream only with an absorption length of about 300 kpc
or more. None adopted.

## Round 20: the no-hold rule on today's law (full tier, against the round-12 baseline)

The results README §31.7; the run is kept in `../run-no-hold-v20/`. A registered comparison, held to the same standard as
the release length (neither is derived); not adopted.

| Law | Pass / close / fail | What moves |
|---|---|---|
| round 12 (adopted) | 59 / 11 / 7 | |
| `no_hold_r12` | **61 / 12 / 4** | better: Carina 3.47 → 4.39 km/s (pass), Antlia 2 1.11 → 4.14 (pass), Sextans 2.10 → 4.45 (close), Crater II 1.00 → 3.35 (close), the ten dwarfs' χ² 137.7 → 60.4; unchanged failures: Draco 2.79 → 4.10 and Ursa Minor 3.44 → 4.18 (against 9.1 and 9.5), the Bullet's smaller half, KiDS's red lenses; changed predictions: Cassini's Q2 → 0 (the Galaxy's pull no longer enters the Sun's law), wide binaries 1.076 → 1.182 at 20,000 AU |

## Round 21: the supplied distance law, in parts (full tier, against the round-12 baseline)

The results README §32.3; the reports are kept in `../run-distance-eta-v21/suite/`. Registered comparisons, a and u
refitted on their home data; none adopted.

| Law | Pass / close / fail | What moves |
|---|---|---|
| round 12 (adopted) | 59 / 11 / 7 | |
| `dist_fixed_eta_r12` (η = 1/2, H₀-like 73.2, D_A = D*) | **62 / 8 / 7** | improved: KiDS all, red and disk lenses (+0.035, +0.048, +0.036 dex), El Gordo NW speeds; masses up in MACS J0025, El Gordo and the Bullet's smaller half (1.50 against 2.45–2.81); regressed: Mistele's ellipticals 2.55 → 3.71 (fail), Abell 520 P2 (close) |
| `dist_metric_eta_r12` (η = 1/2, H₀-like 73.2, D_A = D*/(1 + z): the proposal) | **62 / 8 / 7** | improved: KiDS all, red and disk lenses (+0.027, +0.042, +0.029 dex); regressed: Mistele's ellipticals 4.03 (fail); the Bullet's smaller half 0.91 against 1.89–2.17 |
| `dist_metric_r12` (the metric geometry alone, η = 0, H₀-like 70.9) | 57 / 13 / 7 | regressed: MACS J0025's galaxy speeds (631 km/s, fail), El Gordo's masses (close); KiDS red fail → close |

The path factor (with its own scale) is what helps; the metric geometry alone costs the far clusters.

## Round 24: the Casimir-EFT candidate (full tier, against the round-12 baseline)

The results README §35; the reports are kept in `../run-eft-v24/suite/` (and `../run-eft-v24/suite_table_v24.json`,
from `code/eft_suite_table_v24.py`). Registered comparisons; none adopted. (Cassini's Q2 shows as "worse" for every
law, the adopted one included: its target became Park et al. 2026 in round 18, after this baseline was saved.)

| Law | Pass / close / fail | What moves |
|---|---|---|
| round 12 (adopted) | 59 / 11 / 7 | |
| `eft_r12` (the note's local equation, hold kept) | 47 / 13 / 17 | regressed: X-COP's miss and radial trend (0.486, 0.715), the Bullet (outer stars, smaller half's κ 0.067 and peak 845 kpc, main mass 2.39), the stack (β 0.180), MACS J0025 (NW peak 249 kpc, speeds 571 km/s), Abell 520 (P4, 710 kpc, speeds), El Gordo (masses, NW speeds) |
| `eft_r12_refit` (a 5.98 × 10⁻¹¹, u 98.6) | 50 / 12 / 15 | the collisions mostly back (21 / 3 / 3); regressed: SPARC (16.7, close; bulges fail), X-COP trend (0.419), KiDS bulge and GAMA, both Mistele samples (ellipticals z 11.3), SLACS (−0.108) |
| `eft_nohold_r12` (no hold) | 41 / 12 / 24 | regressed: the planets (3 × 10⁻³ of the Sun's pull), S2, SPARC (18.3), microlensing, SLACS, the Milky Way's vertical pull, and the clusters as `eft_r12`; improved: the Sun's speed 234 km/s |
| `eft_memory_r12` | 50 / 12 / 15 | as `eft_r12` outside the collisions; the Bullet's smaller half recovers its peak, the rest does not |
| `eft_excess_r12` (the repair) | **59 / 11 / 7** | nothing changes status; the collision numbers move by 1% or less (peaks by up to 4 kpc) |
| `eft_excess_local_r12` (the repair, no memory) | **58 / 11 / 8** | regressed: MACS J0025's NW peak (259 kpc); the Bullet's smaller half κ 0.114 (still passes), stack β 0.061 |
| `eft_excess_aqual_r12` (the repair solved exactly, curl field kept; collisions only, the rest is the law) | **58 / 11 / 8** | regressed: MACS J0025's NW peak (213 kpc); the Bullet's peaks 27 and 50 kpc, gas lensing 0.064 / 0.084, κ 0.695 / 0.222 (all pass); stack β 0.021 |
| `eft_scalar_r12` (brightness only) | 59 / 10 / 8 | regressed: MACS J0025's NW peak (223 kpc); improved: Abell 520 P6 |
| `eft_scalar_local_r12` | 51 / 11 / 15 | the lensing follows the gas: Bullet peaks 265 and 617 kpc, stack β 0.98, MACS J0025 both peaks |

## What the candidates do now (quick tier, against the round-9 baseline)

| Candidate | Fixes (fail or close → pass) | Breaks (pass → close or fail) | Moved without changing grade |
|---|---|---|---|
| no hold, release over about 1 pc | Carina and Antlia 2 → pass; Sextans and Crater II → close | nothing | Draco 2.8 → 4.1, Ursa Minor 3.5 → 4.2 km/s (still fail); dwarfs χ² 135 → 60; Cassini Q2 → 0; wide binaries 1.09 → 1.19; Fornax 12.36 → 12.58 (11.7 ± 0.9) |
| external hold 10% | Carina | wide binaries: 1.72 × Newton at 20,000 AU (close) | every dwarf rises (χ² 82) |
| g_d × 1.25 | nothing (the Sun: 211 → 215 km/s, still close) | pull above the disk: 74.1 → 76.4 | strong lenses −0.017 → −0.037 dex |
| g_d × 1.5 | the Sun: 217.6 km/s | strong lenses (−0.057 dex); pull above the disk (78.1) | SPARC median residual 0.031 → 0.015 dex |

Quick-tier totals: no_hold 41 / 5 / 2 against the baseline's quick subset of 39 / 4 / 5. With the
collisions, which read neither the hold nor the release, that is 64 / 9 / 3 in all. The refitted
candidates were not rerun: the release changes only Cassini and the binaries, so their round-8
results below stand, with Cassini now passing for `gd_x1p5_refit`. (`combined` already had the
release.)

## What the candidates did in round 8 (quick tier, against the round-3 baseline)

| Candidate | Fixes (fail or close → pass) | Breaks (pass → close or fail) | Moved without changing grade |
|---|---|---|---|
| gradual release, L = 30,000 AU | Cassini Q2: 3.1 × 10⁻²⁶ → 4.6 × 10⁻²⁷ s⁻² | nothing | wide binaries at 20,000 AU: 19% → 9% extra pull (4% at 7,000 AU) |
| external hold 10% | Cassini Q2 (3.8 × 10⁻²⁷); the dwarf Carina | wide binaries: 2.5× Newton's pull at 20,000 AU, beyond both published analyses | every dwarf rises (χ² 135 → 82); the five faint ones still fail |
| g_d × 1.25 | nothing (the Sun: 211 → 215 km/s, still close) | pull above the disk: 74.1 → 76.4 (just past 2σ) | strong lenses −0.017 → −0.037 dex; SPARC median residual halves |
| g_d × 1.5 | the Sun: 217.6 km/s | strong lenses (−0.057 dex); pull above the disk (78.1) | SPARC median residual 0.031 → 0.015 dex |
| g_d × 1.5, a and u refitted | the Sun (216.8) | the same two, and the lensing speeds of spirals (Mistele et al.) | a 4.6% lower, u 196.4 km/s |
| all three, refitted | the Sun, Cassini, Carina | the same three, and wide binaries (1.76× Newton, close) | Cassini Q2 1.2 × 10⁻²⁷ |

Read across (round 8): gradual release costs nothing anywhere else (adopted in round 9); a 10% external hold helps the dwarfs
and Cassini but makes wide binaries far stronger than either published analysis allows, so the
hold has to weaken for the dwarfs without weakening for binaries (they differ: a dwarf moves past
the Galaxy's companion at 100–300 km/s, a binary's stars move with it); moving the switch-off
later trades the Sun's speed against the pull above the disk and the strong lenses. The Milky
Way's own visible matter is the other lever: with Bovy & Rix's (2013) shorter disk instead of
McMillan's, our law gives 217 km/s at the Sun (round 7), so the disk's shape comes next.

Round 9 adds that the adopted release length, if it is 0.5 pc or more, protects the binaries
whatever the hold. So the dwarfs' hold no longer has to depend on speed to spare the binaries,
although it still needs a physical reason (results README §19.2).

**Known borrowed assumptions** (results README §19.4, §20.3):
* The SLACS check moved to the project's static distances in round 10: light = matter −0.012 dex;
  the stars needed are 1.44–1.95 × Salpeter.
* The three collisions moved to the static distances in round 10 (`code/collisions_v10.py`), each
  lensing mass converted with its paper's own cosmology and source redshifts, El Gordo against
  aperture masses. The suite's tally went from 62 / 8 / 6 to 58 / 11 / 7; every grade that moved
  is a far cluster whose stars came out 20–29% lighter. Their star masses still carry the
  Big-Bang age cap or fixed light-to-mass ratios; with 1.4 times the stars the tally would be
  63 / 7 / 6 (results README §20.6).
* Round 11 moved the rest (results README §21):
  * X-COP (`t_clusters.static_clusters`): profiles converted to the static distances
    (`code/xcop_static_v11.py`), with the stellar profiles deprojected. The release's profiles are
    projected (cylinder) masses, which rounds 1–10 used as spherical ones;
  * KiDS and Mistele (`t_lensing`): g_bar × stars/size², g_obs × Σ_crit(static)/Σ_crit(ΛCDM), at
    z_l = 0.25 with sources at 0.75 (`code/kids_static_v11.py`);
  * the Bullet Cluster (`t_collisions.bullet`): inputs and targets converted at z = 0.296
    (`code/bullet_static_v11.py`);
  * MACS J0025's star masses put on X-COP's (Chabrier) basis, × 10^−0.25
    (`collisions_v10.chabrier_basis`).

  With the constants refitted in these conventions (the adopted law), the tally is 57 / 9 / 10
  against round 10's 58 / 11 / 7. Five grades improved: Abell 520's P4, its mass inside 710 kpc
  and its galaxy speeds; El Gordo inside 0.5 and 1 Mpc. Six regressed: KiDS all (close), blue
  and disks (fail) and the early/late gap (close); Mistele's spirals (fail); MACS J0025's NW
  lensing peak (fail; `code/macs_peak_scan_v11.py`: at 0.5 Gyr with Chabrier-basis stars the NW
  side is a flat ridge from the galaxies to the gas, highest at the gas; at ≤ 0.35 Gyr, or with
  1.33 times the stars, the peak returns to the galaxies).
  The KiDS level scales with the static distance law's α (Σ_crit ∝ α; g_bar does not depend on
  α), and the gap with u.
