# Regression suite for hot-companion gravity

Every test the law has passed or failed since round 3, in one place, scored the same way, so
that any change to the law, its constants or the code can be checked against everything at once.

```bash
cd research_work/results/hot-companion/regression
python run_suite.py                               # the round-3 law, quick tier (about 2 minutes)
python run_suite.py --tier full                   # + the colliding clusters (about 12 minutes)
python run_suite.py --law gradual_release         # a candidate law from candidates/
python run_suite.py --only dwarfs,precision       # some groups only
python run_suite.py --tier full --save-baseline   # make this run the new baseline
```

The runner prints a scoreboard and writes `runs/<law>-<tier>/results.json` and `report.md`
(both ignored by git). The exit code is 1 when any check got worse in status than in
`baseline.json`, or crashed.

## The law under test

`law_config.py` turns a candidate file into one dictionary that every test reads. Nothing in the
suite reads constants from anywhere else. A candidate is a small JSON file in `candidates/`:

| Key | Meaning | Default (the round-3 law) |
|---|---|---|
| `gd_scale` | the release scale g_d multiplied by this | 1 |
| `release_length_au` | the companion is released gradually over this length around each emitter, R(r) = 1 − e^(−r/L) | 0 (at once) |
| `external_hold` | how strongly a subsystem's companion follows an outside galaxy's pull (dwarfs, the Sun, wide binaries) | 1 (fully) |
| `refit` | constants refitted on their home data after the change: `"a"` on the 149 SPARC galaxies (g_d held), `"u"` on the 12 X-COP clusters | none |
| `a_SI`, `lam`, `u_kms`, `base` | explicit constants, or another `results.json` to start from | run-v3 |

Candidates now in the folder:

| File | What it changes | Why it was proposed |
|---|---|---|
| `gd_x1p5.json` | g_d × 1.5 | the Sun's rotation speed (round 7) |
| `gd_x1p5_refit.json` | g_d × 1.5, then a and u refitted | the same, with the constants re-balanced |
| `gradual_release.json` | release over 30,000 AU (0.15 pc) | Cassini's Q2 limit (round 7) |
| `weak_hold.json` | external hold 10% | the six faint dwarf galaxies (round 7) |
| `combined.json` | all three, a and u refitted | together |

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
standard errors, or its value over a limit. Against the baseline each check is marked
**regressed** / **improved** (status changed), **worse** / **better** (score moved by more than
0.1 or 5%), **known** (a fail already in the baseline) or unchanged.

## What is tested (71 checks in the full tier)

| Group | Checks | Data |
|---|---|---|
| machinery | the standing rule (not MOND, Newton or dark matter: the formula guard); the switch-off in strong pulls; the constants and energy budget (tracked) | `research_work/tools/formula_guard.py` |
| galaxies | SPARC typical speed miss (all, test, validation) against MOND's; median residual overall and where the switch-off acts; bulge-dominated galaxies | Lelli et al. 2016 |
| clusters | X-COP mass miss and its trend with radius; the u X-COP prefers; held-out miss over 924 half splits (full tier) | Ettori et al. 2019; Ghizzardi et al. 2021 |
| lensing | KiDS-1000 lensing pull for six samples and the early/late gap; lensing circular speeds of spirals and ellipticals; SLACS light = matter; Einstein Cross; microlensing | Brouwer et al. 2021; Mistele et al. 2024; Auger et al. 2009 |
| milky_way | speed at the Sun; 15–27 kpc curve (four Gaia analyses); vertical pull at 1.1 kpc; mass inside 20/50/100/200 kpc; escape speed; inner Galaxy | `data/mw_literature_v7.json` |
| dwarfs | speed spread of ten dwarf spheroidals, with the Galaxy's pull and heat | `data/mw_dwarfs.json` |
| precision | planets, S2, the Double Pulsar, light bending, Cassini's Q2; wide binaries (tracked) | Hees et al. 2014 and others |
| collisions (full) | the Bullet Cluster (round-5 case): outer stars, lensing strengths, gas residuals, peak positions, masses inside 250 kpc; the 72-collision stack (β); MACS J0025.4−1222, Abell 520 and El Gordo (round 8) | Clowe et al. 2006; Harvey et al. 2015; see `code/collisions_v8.py` |

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
