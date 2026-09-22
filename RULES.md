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
