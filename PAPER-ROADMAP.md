# What it would take to publish this

**22 September 2026.** An honest assessment of the gap between where the work is and
a paper that survives peer review, with the tasks ordered by what blocks what.

---

## The strategic read, before the list

**The galaxy-scale result is not publishable as new.** Our `a₀ = 1.171 × 10⁻¹⁰ m/s²`
and our interpolation function are Milgrom's, and the acceleration relation we
reproduce was published by McGaugh, Lelli and Schombert in 2016 from the same SPARC
data. A paper whose headline is "we fit rotation curves with one constant" will be
desk-rejected as a rediscovery, correctly.

**The novelty is τ, and only τ.** The claim worth a paper is: *the residual factor of
two that has blocked MOND in clusters for thirty years is the Tolman active-mass
factor of the gravitational field itself, whose two values general relativity fixes
exactly.* Everything else in this project is scaffolding for that one sentence.

**Nature is the wrong first target.** Nature wants a result that is both novel and
settled. Ours is novel and *suggestive*: eleven clusters, a factor that needs a
hydrostatic-bias correction sitting at the edge of the measured range, and a theory
that is not yet written down as field equations. The right sequence is a strong
MNRAS or Physical Review D paper first, then Nature only if someone independent
confirms it. Aiming at Nature now costs a year and produces a rejection.

**One strategic amputation.** The "no expanding universe" position is not required by
any result we have — every test here is local, at z < 0.3, and none of them touches
cosmology. Carrying it into the paper means owing the referee a replacement account
of the CMB, big-bang nucleosynthesis and structure formation. That is a decade of
work attached to a result that does not need it. **Drop it from the paper.** Keep it
in the notebook if you like, but a paper that claims it will not be read past the
abstract.

---

## Tier 0 — Blockers. There is no paper until these are done.

### T0.1 Resolve the contradiction at the centre of the argument

This is the single most important task on the list, and a referee will find it in
five minutes.

Section 3.2 of the notebook argues the extra gravity is **not matter** — that there
is no reservoir, that the force is the source's own field with a longer reach, and
that this is why the energy catastrophe dissolves. Section 2 then uses **τ**, which
is the Tolman active mass of a **stress-energy tensor**.

Those cannot both be true as written. Either:

* **(a) the extra gravity carries a real T_μν** — in which case τ applies cleanly,
  but the 50,000× energy problem comes straight back and must be answered; or
* **(b) it is a modified propagator** — in which case the energy problem stays
  dissolved, but τ has to be *re-derived* as a property of the modified field
  equations rather than borrowed from general relativity's Tolman formula.

(b) is the more promising branch and the one consistent with the rest of the work.
It is real theory, not bookkeeping: what plays the role of "streaming versus
standing" when there is no substance to stream?

**Done when:** the paper can state in one paragraph what the extra gravity *is*,
and τ follows from that statement rather than from analogy.

### T0.2 Write down field equations

A covariant action, or at minimum a set of field equations, that reduces to:

* the interpolation function in the static weak-field limit;
* τ = 1 and τ = 2 as its two limiting stress states;
* a transition between them.

Without this the work is a fitting formula with a physical story attached. PRD will
not take it; MNRAS might, but the referee will ask anyway, and the answer "we don't
have one" caps the paper's reach permanently.

**Note:** T0.2 done properly may hand us T0.3 for free — if the field equations
contain the switch, we stop hunting for it empirically.

### T0.3 Derive what sets τ, or demote it honestly

Four candidates tested and ruled out: local starlight flux, galaxy mass/size/
luminosity/gas fraction, surface density, and support mode (blocked by a stellar-mass
convention). Only *scale* survives as a discriminator, and we cannot yet say why.

Two acceptable outcomes. Derive the switch from T0.2. Or present τ honestly as a
**measured two-valued parameter** — weaker, but still publishable if the measurement
itself is clean, which is what Tier 1 is for.

Not acceptable: fitting a crossover function to the gap between two populations and
presenting it as a prediction. We did that once already (notebook §3.5) and had to
retract it when raw data replaced a smoothed model.

---

## Tier 1 — Observational work the claim rests on

### T1.1 The Milky Way, properly — *partly done, and it works*

**Already run today, with nothing refitted:**

| | RMS vs 542-Cepheid rotation proxy |
|---|---:|
| Visible matter alone | 75.24 km/s |
| **This law, `a₀` from SPARC only** | **15.71 km/s** |
| Repo's earlier model, *fitted to this data* | 8.76 km/s |

Mean fractional error 6.4%. But **every one of the twelve bins is under-predicted**,
by 7 to 22 km/s. That is a systematic, not scatter, and it needs explaining before a
referee does it for us — most likely the baryonic mass model or the Jeans proxy,
possibly real.

**Still to do:** a modern Gaia DR3 rotation curve rather than a Cepheid proxy; the
**vertical force K_z(z)** at the solar radius, which is a genuinely different test of
the same law and one the literature takes seriously; and a baryonic model with
propagated uncertainties.

### T1.2 Cluster masses from weak lensing, not X-rays — *not optional*

The headline cluster result currently needs a hydrostatic bias of b ≈ 0.17–0.195
against a measured 0.10–0.15. That is the same size, but it sits at the edge, and the
whole claim leans on it. Weak-lensing masses for the same eleven clusters remove the
correction instead of applying it. **This is the task that most improves the paper per
unit effort.**

### T1.3 Groups — the decisive measurement

The only systems between the two regimes. With four candidate variables eliminated we
now know what to measure against.

### T1.4 Dwarf spheroidals — the test we might fail

Deepest into the low-acceleration regime, and where MOND-like laws have known
tensions. If the law breaks here we need to know before the referee does. Run it even
though — especially though — it may hurt.

### T1.5 The reserved SLACS lenses

The KCWI release carries 14; only six are usable here. The rest are a genuine holdout
for the redshift trend, and the way to keep them one is to declare the prediction in
writing before computing a single angle.

### T1.6 SPARC with standard quality cuts

Our 16.25 km/s uses no inclination or quality flags, deliberately, so the comparison
against our own old model stayed like for like. For publication it has to be run the
way the field runs it, so the numbers are comparable to what is already in print.

---

## Tier 2 — The four things that kill no-dark-matter papers

Referees ask all of these. Silence on any one is fatal.

### T2.1 The Bullet Cluster
Mandatory. The canonical argument against modified gravity. The paper needs a
position, even if the position is "this is outside our scope and here is why."

### T2.2 Gravitational wave speed — GW170817
The simultaneous arrival of gravitational waves and light killed a large class of
modified-gravity theories overnight. Whatever comes out of T0.2 must be checked
against it. Do this *early* — it can invalidate a theory before it is written.

### T2.3 Solar system and binary pulsars
The law must return to Newton and to standard general relativity where those are
tested to high precision. The interpolation function does this by construction; the
τ mechanism needs checking.

### T2.4 Cosmology — scope it out
See the strategic amputation above. The honest move is to state plainly that the
paper addresses galaxy and cluster dynamics and takes no position on cosmology.

---

## Tier 3 — Rigor the reviewers will check

* **T3.1** Full posteriors with propagated uncertainties, not point fits. Every
  headline number needs a credible interval.
* **T3.2** A blind analysis on reserved data, declared before unblinding.
* **T3.3** A formal prior-art section comparing against MOND, TeVeS, EMOND,
  emergent gravity and superfluid dark matter — stating exactly what is ours. The
  register at `research_plan/prior-art/` is the right home and needs updating.
* **T3.4** Independent reproduction: a clean public repository someone else can run.
* **T3.5** Fix the stellar-mass convention (notebook §3.9). It currently moves
  results by a factor of 2.4 and is the largest known systematic in the programme.

---

## Answers to the three questions asked

**Do we need first principles?** Yes, and it is the top of the list. T0.1 and T0.2
are the difference between "a formula that works" and "a theory." Without them the
ceiling is a solid MNRAS paper; with them, PRD and a real claim.

**Do we need more data?** Yes, but less than you might think, and specifically:
weak-lensing cluster masses (T1.2) and groups (T1.3). Everything else on the list is
confirmation rather than foundation.

**Do we need the Milky Way?** Yes — and it partly works already, with zero parameters
refitted. Finishing it properly with Gaia and the vertical force is one of the
cheaper high-value tasks here.

---

## Suggested sequencing

1. **T0.1** — resolve the contradiction. Everything downstream depends on which
   branch is taken, and it costs thinking rather than compute.
2. **T2.2** — check the gravitational wave constraint early, before investing in a
   theory it might rule out.
3. **T1.2 and T1.3** in parallel — the two measurements that most strengthen the
   cluster claim.
4. **T0.2** — field equations, informed by 1 and 2.
5. **Paper 1** to MNRAS or ApJ: the cluster τ result, scoped to galaxies and
   clusters, cosmology explicitly out of scope.
6. **Paper 2** to PRD: the theory, if T0.2 succeeds.
7. **Nature** only after someone independent reproduces the cluster result.
