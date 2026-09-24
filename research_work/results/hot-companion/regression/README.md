# Regression suite for hot-companion gravity

Every test the law has passed or failed since round 3, in one place, scored the same way, so
that any change to the law, its constants or the code can be checked against everything at once.

```bash
cd research_work/results/hot-companion/regression
python run_suite.py                               # the adopted law (round 9), quick tier (about 1 minute)
python run_suite.py --tier full                   # + the colliding clusters (about 20 minutes)
python run_suite.py --law no_hold                 # a candidate law from candidates/
python run_suite.py --law round3                  # the law before round 9 (released at once)
python run_suite.py --only dwarfs,precision       # some groups only
python run_suite.py --tier full --save-baseline   # make this run the new baseline
```

The runner prints a scoreboard and writes `runs/<law>-<tier>/results.json` and `report.md`
(both ignored by git). The exit code is 1 when any check got worse in status than in
`baseline.json`, or crashed. Each run records the commit it came from, marked `+changes` when
tracked files differ from it.

**The baseline** (`baseline.json`) is the adopted law, round 9, on the full tier: 89 checks,
62 pass, 8 close, 6 fail and 13 tracked. It changes only when a change is adopted or a test is
added (RULES.md §12).

## The law under test

`law_config.py` turns a candidate file into one dictionary that every test reads. Nothing in the
suite reads constants from anywhere else. The adopted law (round 9) is the round-3 constants with
gradual release over 30,000 AU (0.15 pc), adopted for Cassini. A candidate is a small JSON file in
`candidates/`, a change to the adopted law:

| Key | Meaning | Default (the adopted law, round 9) |
|---|---|---|
| `gd_scale` | the release scale g_d multiplied by this | 1 |
| `release_length_au` | the companion is released gradually over this length around each emitter, R(r) = 1 − e^(−r/L) | 30,000 AU (adopted in round 9; `round3` has 0, released at once) |
| `external_hold` | how strongly a subsystem's companion follows an outside galaxy's pull (dwarfs, the Sun, wide binaries) | 1 (fully) |
| `refit` | constants refitted on their home data after the change: `"a"` on the 149 SPARC galaxies (g_d held), `"u"` on the 12 X-COP clusters | none |
| `a_SI`, `lam`, `u_kms`, `base` | explicit constants, or another `results.json` to start from | run-v3 |

Candidates now in the folder:

| File | What it changes | Why it was proposed |
|---|---|---|
| `gd_x1p5.json` | g_d × 1.5 | the Sun's rotation speed (round 7) |
| `gd_x1p5_refit.json` | g_d × 1.5, then a and u refitted | the same, with the constants re-balanced |
| `weak_hold.json` | external hold 10% | the six faint dwarf galaxies (round 7) |
| `no_hold.json` | external hold 0, release over 200,000 AU (about 1 pc) | the dwarfs, once the release length protects the binaries (round 9) |
| `gd_x1p25.json` | g_d × 1.25 | the smallest useful step toward the Sun's speed |
| `combined.json` | all three, a and u refitted | together |

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

Every test calls the same code that produced the published numbers (`../code/`); the baseline
reproduces them to the last digit (SPARC 15.85 km/s, X-COP 0.227, KiDS +0.024 / −0.005 / +0.021,
Bullet κ 0.675 and 0.14, collision stack β 0.027).

## Adding a test

Write a function in the right `t_*.py` module that returns `checks.make(...)` objects: an id, a
plain-language title, the law's value, and a criterion from `checks.py` (`z_check`,
`range_check`, `floor_check`, `at_most`, `rms_z`). Then run the full tier and save the baseline.
The Milky Way's Newtonian fields are cached in `cache/` (they do not depend on the law).

## Adding a candidate

Copy a file in `candidates/`, change the keys, and run
`python run_suite.py --law <name> --tier full`. The scoreboard lists what the change fixes
(improved), what it breaks (regressed), and what moved without changing status.

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
* The Bullet Cluster, KiDS, Mistele and X-COP comparisons still use the papers' ΛCDM
  conversions.
