# Standing rule: no MOND, no plain Newton, no dark matter, and check every formula

**Recorded 22 September 2026, at the owner's instruction:**

> "make a note that we are not to use mond or any derivation of it, make sure we always
> check every formula against it. same with normal newtonian physics and dark matter.
> then go back to our creative solutions that don't work anything like these and make
> them work."

This rule applies to every candidate from now on. It overrides earlier notes where they
conflict.

---

## 1. What is not allowed as our answer

| Excluded | What that means in practice |
|---|---|
| **MOND and anything derived from it** | We do not adopt, as our law, any of the following: Milgrom's MOND; any interpolation function ν or μ; AQUAL or QUMOND field equations; TeVeS; AeST; EMOND; or any other "MOND plus a correction" construction. This includes the law this project used until today, `g = g_N/2 + √(g_N²/4 + g_N·a₀)`, which is MOND's "simple" function. |
| **Plain Newtonian gravity** | A law that is only Newton, including Newton with rescaled masses or a changed constant, is not a solution. |
| **Dark matter** | We do not add an invisible mass component with its own freedom: halos fitted per galaxy, new particles, or any extra mass that is not fixed by the ordinary matter we see. |

## 2. "Check every formula against it": the procedure

Every candidate formula goes through two checks before it is reported.

1. **The guard.** `research_work/tools/formula_guard.py` runs these tests:
   * **MOND equivalence, isolated mass.** For one isolated, cold mass, does the candidate
     reduce to a single function of the Newtonian acceleration? If so, which published
     interpolation function does it match, and how closely?
   * **MOND equivalence, real galaxies.** Across 3,150 SPARC measurements, how much of the
     candidate's prediction cannot be written as a function of the local Newtonian
     acceleration alone?
   * **Newton equivalence.** Is the candidate's boost a constant in every galaxy? That is
     just Newton with rescaled masses.
   * **Dark-matter equivalence.** Does the candidate carry free per-object parameters for
     an extra component, or an extra substance with its own conserved amount?
2. **The benchmark.** The candidate is scored on the same data as Newton (ordinary matter
   only), MOND (one fitted constant) and dark matter (an NFW halo fitted to each object).
   The number of adjustable numbers is reported beside every score.

A formula is reported as ours only if the guard finds it differs from all three somewhere.
That difference must also be tested against data.

## 3. What the data themselves force: stated honestly

Disk-galaxy rotation obeys an observed law, the *radial acceleration relation*, very
tightly. On 149 SPARC galaxies we tested three laws built on different ideas from MOND
(22 September 2026):

| Law | Mean velocity error (lower is better) |
|---|---:|
| Companion pool: pull set by radius and the galaxy's total mass | 20.35 km/s |
| Stream density: pull set by the total of all contributions, no cancellation | 19.40 km/s |
| Stream density, weighted by net flow | 18.68 km/s |
| MOND-family context laws | 16.13 – 17.00 km/s |

All three came out worse. There is also a mathematical reason. For one isolated mass,
suppose a law uses nothing but gravity's strength, one acceleration constant and the mass.
Dimensional analysis then forces it to be a function of the Newtonian acceleration alone,
and so to have MOND's form.

**What follows.** In cold, isolated disk galaxies, any law that fits the data agrees with
MOND, just as any law must agree with Newton in the Solar System. That agreement is not
"using MOND", provided:

* our equations come from our own mechanism and do not take MOND's formulas as inputs;
* the guard reports any place where our result coincides with a published MOND function;
* our law differs from MOND where MOND is silent or wrong, and that difference is tested
  on data.

## 4. What this retires

* The dynamics equation in [THEORY.md](THEORY.md), which is Milgrom's QUMOND, and the AeST
  "relativistic home" built on it. **Both are MOND. They stay on the record, but they are
  no longer our solution.**
* The gravitational-slip results (η) were measured *relative to that MOND law*. They must
  be re-derived against the new law before they are quoted as ours.

## 5. Where this led

The first candidate that passes the rule is recorded in
[research_work/results/hot-companion/](research_work/results/hot-companion/README.md).

## 6. MOND as an output, not an input (23 September 2026)

Round 2 derived MOND's deep law, its constant (`a₀ = 2ℓ/u`) and an interpolating function
from the hot-companion mechanism, as the limit of cold matter
([hot-companion §8.1](research_work/results/hot-companion/README.md)). This is consistent
with the rule. It is the case §3 anticipated:
* no MOND formula went in;
* the guard still flags the cold limit as MOND-like, and it should;
* every result we claim as ours rests on where the law is **not** MOND: hot gas in clusters,
  hot haloes around ellipticals, and the release of the companion in strong fields.

Our derived function, `ν(y) = 1 + e^(−y/λ)/√y`, is not any published one. The closest is
"simple", 0.031 dex away. It must not be swapped for a published function to simplify
calculations.

## 7. Round 3 passes the check (23 September 2026)

The round-3 law (stars and galaxies carry the heat; gas collides and counts as cold; the
pull follows the companion's net flow) was run through the guard:
* **Cold, isolated mass.** A MOND form: `ν = 1 + e^(−y/λ)/√y`, λ = 3.45. The nearest
  published function is RAR-exponential, 0.030 dex away. As §6 says, this is an output, not
  an input.
* **Same Newtonian pull, random speeds of 0–1000 km/s.** Predictions differ by 0.93 dex, so
  **not a function of g_N alone**.
* **SPARC.** 0.022 dex of the prediction is beyond any local function of g_N.
* **No per-object parameters, and no extra conserved substance.** The companion's amount is
  fixed by the visible matter.

Details: [research_work/results/hot-companion/README.md](research_work/results/hot-companion/README.md) §10.

## 8. Round 4 passes the check (23 September 2026)

Round 4 adds no formula. It makes explicit that the companion travels at u relative to the
matter that emitted it, which Galilean invariance requires, and follows the consequence after
collisions:
* **Not dark matter.** The pre-collision companion around a collided cluster is fixed by the
  visible matter's own history. It has no adjustable amount, it fades as the fresh sphere
  grows at u, and it is absent in settled systems.
* **No new parameter.** The pre-collision mass ratio (1:6–1:10) comes from merger
  reconstructions, and the star mass-to-light ratio (1–1.5) stays inside the published
  0.5–3. Neither is fitted to the lensing.
* **Settled systems are unchanged,** so round 3's guard verdict stands: the cold limit is a
  derived MOND form (RAR-exponential, 0.030 dex away), and predictions differ by 0.93 dex at
  fixed g_N.

Details: [research_work/results/hot-companion/README.md](research_work/results/hot-companion/README.md) §12.

## 9. Round 5 passes the check (23 September 2026)

Round 5 adds no formula. It rebuilds the Bullet main cluster's outskirts and reads the lensing
data correctly:
* **Not dark matter.** The two numbers inferred are both amounts of visible matter, checkable by
  counting stars and galaxies:
  * the main cluster's outer stars, from the galaxy speed (star/gas 0.048 at R500, inside the
    X-COP range);
  * the subcluster's pre-collision size, about 1:3, from its lensing mass.
  Neither is an invisible or free component.
* **Settled systems unchanged; the guard verdict of §7 stands.**
* **The blog's new §4** explains each piece of the law with hypotheses for its origin. The cold
  limit is still reported as MOND-like, as it must be.

Details: [research_work/results/hot-companion/README.md](research_work/results/hot-companion/README.md) §14.

## 10. Round 6 passes the check (23 September 2026)

No new formula. Round 6 tests round 5's inferred amounts of visible matter against data:
* the main cluster's outer stars: confirmed by a Legacy Survey star count;
* the subcluster's lost galaxies: not seen, so they are withdrawn as an explanation.

This is the rule working as intended: inferred matter must be found, or the inference is
dropped.

The wide-binary prediction uses the law unchanged, in the cold limit. It is reported against
MOND (43%) and Newton (0%) and sits at 19%.

Details: [research_work/results/hot-companion/README.md](research_work/results/hot-companion/README.md) §16.

## 11. Round 7: the full check against the standard tests (23 September 2026)

No new formula. Round 7 runs the round-3 law, unchanged, against every standard Galactic and
lensing measurement we could trace to its paper. It adds one standing practice:

* **The coverage table** (BLOG §7; results README §17.8) is re-run whenever the law changes.
  Every row is a published measurement, with MOND, Newton and dark matter alongside.
* **A failed row triggers the owner's standing instruction**: we are looking for solutions,
  so first say what would have to change to solve the problem, then test it.
  Round 7's three candidate amendments are recorded as candidates, not adopted:
  * gradual release, for Cassini;
  * a larger g_d, for the Sun's speed;
  * a weaker external hold, for dwarfs.

  Each must pass the formula guard and the refit before it enters the law.

Details: [research_work/results/hot-companion/README.md](research_work/results/hot-companion/README.md) §17.

## 12. Round 8: the regression suite (23 September 2026)

No new formula. Round 8 turns the coverage table into code and adds one standing practice:

* **Every change is run through the regression suite before it is reported**:
  a change to the law, to its constants, or to the code that computes any test.
  The suite is `research_work/results/hot-companion/regression/`
  (`python run_suite.py --tier full --law <candidate>`).
  * Each check compares one number with one published measurement and is graded the same way
    everywhere: pass (within 2 standard errors), close (within 3), fail.
  * The report states what the change fixes (improved), what it breaks (regressed), and what
    moved without changing grade.
* **A regression triggers the owner's standing instruction**: say what would have to change to
  recover it, then test that.
* **The baseline (`baseline.json`) records the adopted law.** It changes only when a change is
  adopted or a new test is added, in the same commit.
* **Candidate amendments are JSON files in `candidates/`**, never edits to the law's code.
  The suite's README lists which tests each amendment reaches.

Details: [research_work/results/hot-companion/README.md](research_work/results/hot-companion/README.md) §18.
