# What it would take to publish this

**22 September 2026.** An honest assessment of the gap between where the work is and
a paper that survives peer review, with the tasks ordered by what blocks what.

---

## Status — 23 September 2026 (hot-companion rows updated after round 6)

| Task | Status | Result |
|---|---|---|
| **RULE** No MOND, no Newton, no dark matter; check every formula | **In force** | Recorded in [RULES.md](RULES.md). The automatic check lives in `research_work/tools/formula_guard.py`. It correctly flags MOND itself. The MOND law used until now, and everything built on it (T0.2's dynamics equation, the AeST home, the η values), is retired as our answer. |
| **NEW** Hot-companion gravity | **Works on galaxies, clusters and the Solar System** | Heat feeds the companion without cancelling, and strong fields hold it back. **Galaxies:** 15.93 km/s, against MOND's 16.13, better on train, validation and test. **Clusters:** error 0.223, against MOND's 1.062; held-out cross-check 0.234. MOND refit on clusters only reaches 0.338, and needs a constant 9.7 times its galaxy value. **Solar System:** zero anomalous pull. Three universal constants. **Round 2 (rows below):** MOND is derived as its cold limit, it has a field equation and a momentum-conserving action, and the KiDS, lens and slip problems are resolved or made testable. **Open:** merging clusters. Details in [research_work/results/hot-companion/](research_work/results/hot-companion/README.md). |
| **NEW-2** MOND from our law | **Done** | Cold matter emits in step, so companion flows add like Newton's field (identity to 1.8 × 10⁻¹⁵). Heat Doppler-scrambles them into a plain sum, shown numerically. Energy balance gives the deep law with **a₀ = 2ℓ/u** (ℓ = 2.9 × 10⁻⁵ W/kg). The resulting interpolating function, `1 + e^(−y/λ)/√y`, is new: 0.031 dex from "simple". MOND is an output, not an input. Test: the Sun's extra mass loss, 1.4 × 10⁻¹⁴ per year. `code/derive_mond.py`. |
| **T0.2′** Field equation and action for the new law | **Done (non-relativistic)** | Conservative potential form; the cold limit is QUMOND with our ν. The heat term as first written broke momentum (25%); an action adds a reaction force on hot matter and restores it (1.3% at the largest box, shrinking with box size). `code/field_equation.py`. Relativistic form open. |
| **KiDS** Early/late lensing offset | **Reproduced** | Hot gas haloes of about one stellar mass (0.6–1.0 keV) give +0.20 to +0.22 dex beyond 300 kpc, against the observed ≥ 0.2. MOND with the same haloes gives 0.14. Prediction: the offset grows with radius. `code/kids_haloes.py`. |
| **T3.5** Stellar-mass convention | **Done** | Standard distances plus published Chabrier masses remove 0.13–0.16 dex. Our law then needs stars 1.2–1.6× Salpeter in the six SLACS lenses; dark-matter models need about 1.0. Spectroscopic IMFs are the test. `code/lenses_t35.py`. |
| **Slip** Re-derived against the new law | **Not needed** | Lensing and resolved kinematics agree to +0.040 ± 0.026 dex. The earlier η ≠ 1 was an artefact of measuring against the MOND law. Cluster X-ray masses sit 22% below our prediction (b = 0.22), in the direction of the known hydrostatic bias. |
| **ROUND 3** Collision rule + flow direction | **Done; law locked in** | Gas collides, so its companion stays in step (Dicke narrowing): only stars and galaxies carry heat. The companion pulls along its net energy flow. Refit: u = 197 km/s; galaxies 15.85 km/s (MOND 16.13); clusters 0.227 (held-out 0.244; law-only stellar speeds 0.329). KiDS early/late from early-type stars: 0.17–0.27 dex (obs 0.17/0.27). SLACS: stars 1.05–1.35 × Salpeter; no slip. Momentum is carried by the companion; a matter-only action would push cluster stars outward by 30–70% of gravity, so it is ruled out. Scripts: `code/dicke_toy.py`, `code/run_v3.py`, `code/kids_v3.py`. |
| **T2.1** Bullet Cluster (round 6) | **Main cluster settled; subcluster lensing mass open** | Legacy Survey DR10 star count (1,652 cluster galaxies by photometric redshift): light rises to 2.5–3 Mpc with 23–28% beyond 1.5 Mpc; total 85–103% of the round-5 model inferred from the galaxy speed. Subcluster: the 1:3 solution's lost galaxies are not in Barrena et al.'s 71 non-core velocities (about 3σ for a compact population, 1.4σ for a wide fast one) or in the starlight within 500 kpc (0.46 of the main's central light; core alone 0.32; 1:3 would give 1.55). In our law its lensing needs its hot stars, so the 2× lensing gap is open. Next: the heat gained in the crossing, carried by the memory. Scripts: `code/bullet_light_v6.py`, `code/bullet_members_v6.py`. |
| **WB** Wide binaries (round 6) | **Predicted** | Exact orientation-averaged solution of our field equation in the Galaxy's field at the Sun (1.58 × 10⁻¹⁰ m/s², half the companion released): 19% more pull than Newton beyond 7,000 AU; 1.14–1.26 for ±20% in the Galaxy's field. MOND (simple) 43%. Data disputed (Chae ~1.4; Banik et al. Newton). `code/wide_binaries_v6.py`. |
| **T2.1** Bullet Cluster (round 5, subcluster size superseded by round 6) | **Pattern, speeds and lensing masses reproduced; subcluster size to test** | Main cluster: half of round 4's speed gap was the model cutting its galaxies at 1.5 Mpc; outer stars at star/gas 0.048 at R500 (X-COP 0.035–0.073) give 1,249 km/s and a lensing mass of 2.42–2.54 × 10¹⁴ inside 250 kpc (obs 2.5 ± 0.1, 2.8 ± 0.2). Clowe et al.'s κ are lower bounds, so round 4's "within 1σ" and its lighter stars are withdrawn. Subcluster: its strong-lensing mass (2.0–2.3) needs a pre-collision size of about 1:3 in visible matter (1.93 at 1:3; 1.04 at 1:8); its original galaxies then moved at about 800 km/s (Barrena et al.'s X-ray estimate: about 700). Predictions: about 7 × 10¹² M☉ of its stars now around it; star/gas about 0.05 in the main cluster's outskirts. `code/bullet_main_v5.py`. |
| **T2.1** Bullet Cluster (round 4, strength readings superseded by round 5) | **Pattern and strengths reproduced; main-cluster galaxy speeds open** | The companion is slow (u = 197 km/s) and keeps its emitter's velocity, so each cluster is still wrapped in the companion of the settled cluster it was before the collision, riding with the galaxies. Only a sphere of radius u·t (30 kpc) has been rebuilt around the stopped gas. No new rule; nothing fitted to the lensing. Subcluster: 0.175 (0.13–0.25 for pre-collision ratios 1:6–1:10), against 0.20 ± 0.05. At star M/L 1–1.5 all four strengths agree within about 1σ. Barrena et al. 2002 independently argued the subcluster is the stripped core of a cluster with a pre-merger ratio of about 1:6. Collision stack with memory: β = 0.03 (−0.01 to 0.12), against −0.04 ± 0.07. Prediction: the subcluster's stars move at 460–610 km/s (7 galaxies now give 212 ± 60). **Open:** the main cluster's galaxies come out about 20% slower than the 1,249 ± 100 km/s measured (2.5σ). Scripts: `code/bullet_v4.py`, `code/collisions_v4.py`, `code/bullet_speeds_v4.py`, `code/stream_tidal_v4.py`. |
| **T2.1** Bullet Cluster (round 3, strengths superseded by round 4) | **Pattern reproduced; strength open** | On Clowe et al. 2006's published masses, both lensing peaks sit on the galaxies (8 and 34 kpc from the BCGs) and the gas residuals match (0.04 vs 0.05 ± 0.06 and 0.02 ± 0.06). The old direction rule, hot gas, or no heat each put the peaks on the gas. Strengths: main 0.51 (obs 0.36 ± 0.06; 0.31–0.41 at M/L 1–1.5), sub 0.07 (obs 0.20 ± 0.05). **The subcluster strength is the open item.** 20 modelled collision pieces all keep lensing with the galaxies (Harvey et al. 2015: 5.8 ± 8.2 kpc). Scripts: `code/bullet_v3.py`, `code/collisions_v3.py`. |
| **T2.1** Bullet Cluster (round 2, superseded) | **Partly; top open problem** | Companions that keep their emitter's velocity move the lensing peak 42–75% of the way from gas toward galaxies in a toy merger. But Harvey et al. 2015 (72 collisions) keep lensing on the stars, which needs a much slower companion rebuild on shocked gas than u = 874 km/s gives. `code/bullet_toy.py`. |
| **T0.1** Resolve the contradiction | **Done** | Branch (a) fails on energy by 1.5 × 10⁴; the field's own energy is 1.7 × 10⁻⁴ of what it would need. Branch (b), a modified propagator, survives. τ splits into two objects: in lensing it is the **gravitational slip η** (legitimate, no energy cost); in cluster dynamics it has **no mechanism**. The cluster "closed within scatter" claim is **retracted**. |
| **T2.2** GW170817 | **Done** | Survives, with a construction constraint: the modification and the slip must live in the scalar sector. TeVeS-type completions are excluded; a surviving class exists. |
| **T0.2** Field equations | **Done**; relativistic home identified | Three equations. The dynamics equation reproduces the fitted law to 4.4 × 10⁻¹⁶ and is Milgrom's QUMOND, attributed. The slip equation is ours; JR-10's γ_χ is identically η. **Relativistic:** slip is exactly a radial stress of the extra gravity (the radial Einstein equation). Of four places the extra gravity could live, one passes GW170817 and our data: the aether-scalar-tensor theory (which gives η = 1, verified from the paper) plus a conformal coupling κ of matter to its scalar, giving η = (1−κ)/(1+κ). Clusters need κ ≈ −0.12 to −0.17. What sets κ is open. |
| **T1.2** Cluster lensing masses | **Done, via the published bias** | The X-COP hydrostatic bias implies η = 1.27–1.42, overlapping the under-bent SLACS lenses (1.36–1.54). Degenerate with non-thermal pressure in one comparison; a relaxed-versus-disturbed split at matched f separates them. A published slip measurement from galaxy orbits, which gas pressure does not affect (Pizzuti et al., MACS J1206), gives 1.01 +0.31/−0.28 against our predicted 1.21–1.37. That is consistent at about 1σ. |
| **T1.3** Groups | **Done, inconclusive** | Three SL2S groups with real weak-lensing masses point to η > 1 at face value, but a 30% dispersion bias the source paper itself reports erases it. **Correction:** an earlier "lensing masses ~50% above dynamical" attribution was not in the paper. |
| **T2.3** Solar system (slip part) | **Done** | Slip moves the Cassini γ by at most 3.7 × 10⁻¹², against a bound of 2.3 × 10⁻⁵. Still to check, separately: the simple interpolation function's high-acceleration tail against Cassini (Hees et al. 2016 rule out several popular choices). |

Full working in [THEORY.md](THEORY.md); relativistic slots in
[research_work/results/relativistic-slip/](research_work/results/relativistic-slip/README.md).

**Next, after round 6 (hot-companion gravity):**

1. **The Bullet subcluster's lensing mass** (twice our law's): model the heat its galaxies gained
   while crossing the main cluster, carried with them by the memory; look for diffuse
   starlight moving with it; a few hundred more velocities around it.
2. **A field theory for the companion,** with a momentum flux, its travel time (a retarded
   source) and a relativistic form, so lensing is derived.
3. **Microphysics:** derive `u` and `g_d`; pin the companion wavelength and locking time;
   check that cluster gas scatters faster than about 0.1–1 Myr.
4. **Wide binaries:** a proper external-field prediction.
5. **Tests others can run now:**
   * spectroscopic IMFs of the six SLACS lenses (1.05–1.35 × Salpeter);
   * the KiDS offset flat beyond 100 kpc;
   * galaxy and gas dynamics equal in relaxed clusters;
   * a larger spectroscopic sample of the Bullet's subcluster (460–610 km/s predicted);
   * extra lensing around the smaller clump in older collisions.

**Superseded list (22 September, built on the retired MOND law and slip), kept for the record.
Next, in order of value per effort:**

1. **Stack galaxy-orbit against lensing slip in CLASH-VLT clusters.** This is pressure-free,
   builds on a published single-cluster result, and tests η directly.
2. **Relaxed-versus-disturbed hydrostatic bias at matched f.** This separates slip from
   non-thermal pressure using public catalogues.
3. **The weak-lensing acceleration relation by galaxy type (KiDS-1000).** Brouwer et al. report a
   ≥6σ early/late difference that a universal law cannot make. A type-dependent κ can,
   and this is the most direct handle on what sets κ.
4. **T3.5, the stellar-mass convention**, before re-reading the three over-bent lenses.

**What changed about the paper.** The headline is no longer "τ solves the cluster
problem." It is: **the gravitational slip of the extra gravity is not 1**, measured
object by object with every constant frozen. That is smaller, contradiction-free, not a
rediscovery, and falsifiable. The cluster dynamical deficit stays open and is reported
as open.

---

## The strategic read, before the list

**The galaxy-scale result is not publishable as new.** Our `a₀ = 1.171 × 10⁻¹⁰ m/s²`
and our interpolation function are Milgrom's, and the acceleration relation we
reproduce was published by McGaugh, Lelli and Schombert in 2016 from the same SPARC
data. A paper whose headline is "we fit rotation curves with one constant" will be
desk-rejected as a rediscovery, correctly.

**The novelty is the slip η, and only η.** ~~The claim worth a paper is that the
residual factor of two blocking MOND in clusters is the Tolman active-mass factor of
the gravitational field.~~ *Superseded by T0.1 — see the status table above.* The
claim worth a paper is now: *the extra gravity bends light and moves matter through
potentials that differ, by a slip η ≠ 1 measurable object by object, where every
relativistic MOND theory and every dark-matter halo predicts η = 1.*

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
