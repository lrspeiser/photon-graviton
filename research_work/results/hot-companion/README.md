# Hot-companion gravity: a law that is not MOND, not Newton and not dark matter

**22 September 2026.** Built under the standing rule in [RULES.md](../../../RULES.md):
* no MOND or anything derived from it;
* no plain Newton;
* no dark matter;
* every formula checked against all three.

Reproduce with:

```
cd code
python run.py --output-dir ../run-v2
```

The run takes about 15 seconds and needs no network. All numbers below come from
`run-v1/results.json`.

---

## 1. How we got here: the creative ideas, revisited

We went back to the project's own ideas: matter feeds a companion field; the companion
spreads, attaches and detaches; energy converts to gravity. We then asked what would have
to change for each of them to work.

| Idea | What happened on the data | What it taught us |
|---|---|---|
| **Companion pool.** Each system has a companion whose strength scales with √(its mass), released as a 1/distance pull from every piece of matter | 149 galaxies: **20.35 km/s**, against MOND's 16.13 | Faint galaxies need extra pull right from the centre. A pull that grows with radius can't give it. |
| **Stream density.** Every piece of matter adds companion that thins as 1/distance², with no cancellation | **18.7 – 19.4 km/s** | In cold disks the extra pull tracks the *cancelling* Newtonian pull, not the non-cancelling total. |
| **Heat feeds the companion.** Hot, fast-moving matter feeds it more | Cluster gap against MOND grows with gas temperature (Spearman 0.52, p ≈ 0.09) | Clusters are the one place where matter is hot. |
| **The clusters' own requirement** | The extra source clusters need is 27× their baryons at 0.1 R500 but only 4× at R500 | The need is concentrated, while the gas temperature is flat. |
| **Hot matter adds without cancelling** | Gas outside a point adds 4–12× at 0.1 R500 and about 2× at R500, the same steep shape | **This is the mechanism.** |

## 2. The mechanism, from first principles

*Proposed here; originality unverified; prior art in §6.*

1. **All ordinary matter feeds a companion field.**
2. **Cold, orderly matter feeds it coherently.** Contributions keep a fixed relative
   phase and add as arrows, so they cancel exactly as Newton's pull does. The coherent
   intensity at a point is |g_N|. This is Gauss's law, established physics.
3. **Hot, randomly moving matter feeds it incoherently.** Its contribution is in
   proportion to its random kinetic energy measured against the companion's own speed `u`,
   with weight `k = 3σ²/u²`. Incoherent contributions add as intensities and never cancel:
   `S_hot(x) = G ∫ k ρ(x′) / |x − x′|² dV′`.
4. **The companion pulls with its amplitude.** The pull is the square root of its total
   intensity times a constant `a`, directed along the net pull.
5. **Strong fields hold the companion attached.** It is released only where the ordinary
   pull is weak, like a thermally activated escape. The released fraction is
   `f = exp(−|g_N| / g_d)`, with `g_d = λa`.

```
g  =  g_N  +  exp(−|g_N|/g_d) · √( a · ( |g_N| + S_hot ) )
```

**Three universal constants**, all fitted to data:

| Constant | Value | Fitted on |
|---|---:|---|
| `a` | 6.58 × 10⁻¹¹ m/s² | 149 SPARC galaxies |
| `g_d = λa` | 2.82 × 10⁻¹⁰ m/s² (λ = 4.28) | 149 SPARC galaxies |
| `u`, the companion speed | **874 km/s** | 12 X-COP clusters |

Nothing is fitted to any individual object.

**In plain words:** ordinary gravity is the orderly part, where pulls from opposite sides
cancel. Our extra gravity is a companion that heat feeds without cancelling. Cold disks
barely feed it, which is why disk galaxies look MOND-like. Hot clusters feed it strongly,
which is why MOND fails there and this law does not.

## 3. Results on real data

### Galaxies: 149 SPARC rotation curves, 3,150 points

Mean velocity error, lower is better:

| Law | All | Train | Validation | Test | Adjustable numbers |
|---|---:|---:|---:|---:|---:|
| **Hot companion (ours)** | **15.93** | **15.73** | **19.57** | **13.09** | 3 universal |
| MOND, simple function (context) | 16.13 | 15.85 | 19.80 | 13.52 | 1 universal |
| Newton, visible matter only | 45.58 | 45.94 | 48.64 | 41.68 | 0 |
| NFW dark matter (context) | 7.52 | — | — | — | 298 (2 per galaxy) |

### Clusters: 12 X-COP hydrostatic mass profiles, measured gas density and temperature

`ln(M_measured / M_predicted)`, where 0 is perfect:

| Law | 0.1 R500 | 0.2 | 0.35 | 0.5 | 0.7 | 1.0 R500 | rms |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Hot companion (ours)** | +0.14 | +0.13 | +0.07 | +0.01 | −0.10 | −0.25 | **0.223** |
| Ours without the heat term | +1.36 | +1.26 | +1.10 | +0.96 | +0.77 | +0.54 | 1.048 |
| MOND, simple function | +1.38 | +1.28 | +1.12 | +0.97 | +0.78 | +0.54 | 1.062 |
| MOND with its constant refit **on clusters** (×9.7) | +0.39 | +0.28 | +0.11 | −0.04 | −0.24 | −0.49 | 0.338 |
| Newton | +2.45 | +2.38 | +2.26 | +2.16 | +2.03 | +1.89 | 2.211 |
| NFW (the release's own fit) | | | | | | | 0.101 (24 numbers) |

* **The fit is predictive, not over-fitted.** Across all 924 ways of fitting `u` on 6
  clusters and predicting the other 6, the held-out error is **0.234** (5–95%:
  0.186–0.274). The fitted speed stays at 874 km/s (805–949).
* **Robust.** Treating member galaxies as cold gives 0.230. Cutting the gas at 2 R500
  instead of 3 gives 0.220.
* **The one clear outlier is expected.** A2255, a known merger, is over-predicted at every
  radius. That is what an X-ray mass that reads low should look like.
* **The −0.25 at R500** is the size of the known non-thermal pressure bias there. With a
  typical bias profile added (context only, not fitted), the edge comes to −0.06.

### Solar System

For every planet from Mercury to Neptune, the released companion adds **0 m/s²** (below
10⁻³⁰⁰). The strong-field release is what achieves this. Without it, the law predicted an
extra 6 × 10⁻⁸ m/s² at Saturn, which planetary tracking rules out.

### Lensing

* **Coma weak lensing, 1.7–14 Mpc.** At these distances most of Coma's hot gas lies inside
  the radius, so the law's pull falls as 1/R. That is the shape that already fitted the
  reconstructed shear in the earlier run (χ² = 4.9 for 5 points). We have not recomputed it
  with a Coma gas model, and the amplitude needs source redshifts we don't have.
* **Six SLACS strong lenses: blocked, and a sharp prediction.**
  * Inside their Einstein radii these massive ellipticals are in strong fields, so the law
    makes them nearly Newtonian.
  * The observed rings then need stars 3–5σ heavier than a Salpeter population.
  * Put differently, ≈ 1.5–2.2× Salpeter mass per unit light: a "bottom-heavy" population.
  * Some published studies of galaxies this massive (σ ≈ 255–280 km/s) report exactly
    that; others do not. The known stellar-mass and distance convention problem (roadmap
    T3.5) also enters here.
  * The heat term moves all six lenses in the right direction, by about 0.1σ each.
* **KiDS-1000 ellipticals versus spirals.** Brouwer et al. 2021 measure ≥ 0.2 dex more
  pull around ellipticals than around spirals of equal stellar mass. Their paper says no
  current modified-gravity theory can produce this.
  * Our law predicts the **right direction**, because ellipticals are hot. The size is
    0.02–0.07 dex for σ = 150–300 km/s, **about a tenth to a third** of what is observed.
  * Hot gas haloes around ellipticals, which in this law also feed the companion, may
    supply more. Not yet computed.

## 4. The formula guard

Every check comes from `research_work/tools/formula_guard.py`.

| Law | Verdict |
|---|---|
| MOND simple (control) | **MOND-CLASS**, correctly flagged |
| **Hot companion (ours)** | Only for perfectly cold, isolated matter does it reduce to a MOND-shaped function. Even then it is **not identical to any published one**: the closest is the "simple" function, 0.031 dex away. For the same Newtonian pull, predictions differ by **0.34 dex** between random speeds of 0 and 1,000 km/s, so it is **not a function of g_N alone**. On cold SPARC disks it differs from a pure g_N law by only 0.0016 dex, as the data force. |
| Companion pool (failed) | Point-mass limit identical to the Bekenstein toy function; 0.13 dex beyond local g_N on SPARC |
| Newton test | Not Newtonian: the boost varies inside every galaxy |
| Dark-matter test | No per-object parameters; no independent substance. The companion is fixed by the ordinary matter and its temperature. **Not dark matter.** |

## 5. Distinct predictions, where it parts company with MOND and dark matter

1. **Cluster masses track gas temperature and the extended hot envelope.** The companion
   is concentrated where much hot gas lies outside the point. Colder clusters and groups
   get proportionally less.
2. **Hot systems pull harder than cold ones at the same Newtonian pull.** Ellipticals beat
   spirals, as KiDS sees, and groups sit in between.
3. **Massive ellipticals are nearly Newtonian inside their Einstein radii.** Their lensing
   then requires heavy stellar populations. Spectroscopic IMF measurements can falsify
   this.
4. **Reach.** A companion that spreads at 874 km/s for ~13.8 Gyr reaches ~12 Mpc, so extra
   pull around isolated galaxies should persist to ~10 Mpc and then fall. KiDS sees it
   flat to at least 1 Mpc.
5. **Time.** Recently heated gas (merger shocks) has had less time to spread its
   companion. This is a lever for the Bullet Cluster, not yet computed.

## 6. What is borrowed, and what may be ours

* **Established:**
  * Gauss's law;
  * shell-averaged 1/d² geometry;
  * the SPARC, X-COP, SLACS and KiDS data;
  * the NFW profile, which is dark-matter context.
* **Closest relatives we know of:**
  * **EMOND** (Zhao & Famaey 2012; Hodson & Zhao 2017) raises MOND's constant in deep
    potentials, to rescue clusters.
  * **Superfluid dark matter** (Berezhiani & Khoury 2015) switches behaviour with the dark
    matter's temperature.

  Ours uses neither MOND's constant-changing nor dark matter. It uses the *ordinary matter's*
  random motions, fed in without cancellation. We have not found this combination
  published; a full check is roadmap T3.3.
* **The strong-field release factor** is a new choice here. It makes the cold limit
  MOND-shaped (§4), as the galaxy data force on every theory.

## 7. Open, and next

1. **A field equation.** The law is written for the static weak field. A field equation is
   needed to settle how stars and planets move as whole bodies, to guarantee energy
   conservation, and to compute lensing without assuming light and matter feel the same
   pull.
2. **Why u ≈ 874 km/s and g_d ≈ 2.8 × 10⁻¹⁰ m/s².** Both are measured, not derived.
3. **The KiDS offset.** We get a tenth to a third of it. Compute hot haloes around
   ellipticals in this law.
4. **SLACS lenses.** Settle the stellar-mass convention (T3.5), then test the heavy-IMF
   prediction.
5. **The Bullet Cluster**, through companion attachment and spreading time.
6. **The earlier slip (η) results** were measured against the retired MOND law. Re-derive
   them against this law.

---

## 8. Round 2, 23 September 2026: MOND derived, and the next steps run

Scripts: `code/derive_mond.py`, `code/field_equation.py`, `code/kids_haloes.py`,
`code/lenses_t35.py` and `code/bullet_toy.py`. Their outputs are in `run-derive/`,
`run-field/`, `run-kids/`, `run-lenses/` and `run-bullet/`.

### 8.1 MOND is the cold limit of this law

**The direction matters: MOND comes out, it does not go in.** We start from the companion
mechanism and let the matter cool.

1. **Emission.** Every kilogram of ordinary matter emits companion energy at a steady rate
   ℓ. The energy streams away at speed `u`, conserved.
2. **Cold matter emits in step,** so the energy flows add as arrows. Their sum has exactly
   the shape of Newton's field: `F = (ℓ/4πG) g_N`. This is a geometric identity, checked
   numerically to 1.8 × 10⁻¹⁵.
3. **Why hot matter is different.** Random motion Doppler-shifts every emitter, so the
   cross terms average away.
   * Cold emitters, all in step, give the squared vector sum: 160.50 against 160.49.
   * Hot emitters (σ/u = 0.3) give the scalar sum: 210.55 against 209.20, within 0.6%.
   * So *heat means incoherence* is demonstrated here, not assumed.
4. **The companion's energy density is A²/8πG, and its pull equals its amplitude A.**
   Energy balance through a sphere, `4πr² u A²/8πG = ℓM`, gives exactly
   **A = √(G M a)/r with a = 2ℓ/u** (checked symbolically).
5. **The result satisfies all three of MOND's defining properties:**
   * the deep law g = √(a g_N), which gives flat rotation curves and **v⁴ = GMa**;
   * the deep-regime scale symmetry (the residual is exactly 0);
   * the Newtonian limit, where ν − 1 = 3 × 10⁻¹² at y = 100.

So **MOND is what this law becomes for cold matter**, and it explains MOND's two unexplained
pieces:

| MOND's ingredient | In MOND | In hot-companion gravity |
|---|---|---|
| The constant a₀ | A new constant of nature | **a₀ = 2ℓ/u**: twice the companion power per kilogram, divided by the companion speed. With a = 6.58 × 10⁻¹¹ m/s² and u = 874 km/s, ℓ = 2.9 × 10⁻⁵ W/kg |
| The interpolating function | Chosen by hand | Fixed by the release physics: `ν(y) = 1 + exp(−y/λ)/√y`. It is new: the closest published function is "simple", 0.031 dex away; RAR-exponential 0.034; standard 0.090 |
| Why galaxies obey it | Postulated | Disks are cold |
| Why clusters don't | Unexplained (needs extra matter) | Cluster gas is hot |

**A concrete, testable consequence.** At ℓ = 2.9 × 10⁻⁵ W/kg, all matter converts
1.4 × 10⁻⁴ of its rest energy over 13.8 Gyr. The Sun, hot inside, would emit 7.7 × 10²⁵ W
(0.20 L☉) into the companion. That is an extra mass loss of **1.4 × 10⁻¹⁴ per year**, about
15% on top of sunlight and solar wind, which is at the edge of what planetary ranging
measures. The number holds if companion energy emitted inside the attachment radius is not
re-absorbed.

### 8.2 The field equation, and why it needs an action

```
∇²Φ_N = 4πGρ
S_hot = G ∫ (3σ²/u²) ρ / d²
∇²Φ  = −∇·h,   h = g_N + exp(−|g_N|/g_d) √(a(|g_N| + S_hot)) ĝ_N,   g = −∇Φ
```

* **The pull is conservative.** It is the gradient of Φ. In spherical symmetry it equals
  the algebraic law used on the data.
* **The cold limit is Milgrom's QUMOND, with our derived ν.**
* **3D check.** A disk galaxy with a hot bulge on a 160³ grid: between 2 and 30 kpc, the
  in-plane pull from the field equation is −9% to +2% from the algebraic law.
* **Momentum: a real flaw, found and fixed.** Written as above, the heat term breaks
  Newton's third law: an isolated hot-plus-cold pair pushes itself along.

| Box size | Cold | Hot, as written | Hot, with the action's reaction force |
|---:|---:|---:|---:|
| 160 kpc | 0.035 | 0.224 | 0.074 |
| 240 kpc | 0.017 | 0.251 | 0.032 |
| 360 kpc | 0.009 | 0.258 | **0.013** |

The table gives the net force on the pair divided by the force on one body. A residual that
shrinks with the box is truncation; a violation does not shrink.

**The fix is an action,** QUMOND-structured:
`L = −(1/8πG)[2∇Φ·∇ψ − W(|∇ψ|, S)] + ρ(v²/2 − Φ)`, with
`W = q² + 2∫₀^q e^(−q′/g_d) √(a(q′+S)) dq′`.
* It reproduces the field equation above.
* It hands hot matter an extra reaction force, `(k/8π) ∇ ∫ (∂W/∂S)/d²`.
* With that force included, momentum is conserved (last column).
* In real clusters the reaction is small: about 0 at the centre and **+5% (median) extra
  inward pull on the gas at R500**.

### 8.3 KiDS ellipticals versus spirals: hot haloes close the gap

The model gives each elliptical a hot gas halo (β = 0.5, 300 kpc) in addition to its own
stellar heat. The offset is log g(elliptical) − log g(spiral) at equal stars and cold gas:

| Halo | at 30 kpc | 100 kpc | 300 kpc | 1 Mpc |
|---|---:|---:|---:|---:|
| none (stars' heat only) | +0.020 | +0.025 | +0.027 | +0.027 |
| 1 stellar mass at 0.6 keV | +0.031 | +0.084 | **+0.197** | **+0.186** |
| 1 stellar mass at 1.0 keV | +0.034 | +0.099 | **+0.220** | **+0.205** |
| MOND with the same halo | +0.006 | +0.035 | +0.139 | +0.137 |

* **Hot haloes of about one stellar mass** at the X-ray temperatures typical of these
  galaxies reproduce the observed ≥ 0.2 dex beyond 300 kpc.
* **MOND counts only the halo's mass, not its heat.** At 100 kpc it needs 8.3 stellar
  masses for 0.2 dex, where ours needs 3.0–5.0.
* **Prediction:** the offset grows with radius.

### 8.4 SLACS lenses: the bookkeeping fix, and the heavy-star test

We recomputed the lenses with standard ΛCDM distances and the published SLACS Chabrier
masses (Auger et al. 2009), instead of the project's static-universe convention (T3.5):

| | Project convention | Standard convention |
|---|---:|---:|
| Extra stellar mass needed by our law, lensing | +0.46 to +0.61 dex | **+0.34 to +0.45 dex** |
| Relative to Salpeter | — | **1.23 – 1.57×** |
| Lensing-minus-kinematics gap (slip test) | +0.051 ± 0.026 | **+0.040 ± 0.026 dex** |

* **The convention alone accounted for about 0.13–0.16 dex (a factor of about 1.4)** of the
  earlier lens tension.
* **The published dark-matter-based IMF** for these σ (Posacki, Cappellari & Treu 2015,
  `log α = 0.38 log(σ/200) − 0.06`) is α ≈ 0.95–0.99. It assumes a halo supplies part of the
  mass.
* **Spectroscopic measurements that don't assume dark matter** find bottom-heavy,
  super-Salpeter IMFs in massive ellipticals (van Dokkum & Conroy 2010; Conroy & van Dokkum
  2017).
* **Clean test:** spectroscopic IMFs of these six lenses. We predict 1.2–1.6× Salpeter;
  dark-matter models predict about 1.0.

### 8.5 The slip, re-derived against the new law

* **Lenses.** The stellar mass that lensing needs and the mass the resolved kinematics
  need differ by +0.040 ± 0.026 dex, consistent with zero. Light and matter feel the same
  pull (η = 1).
* **Clusters.** At R500 the X-ray mass sits 22% below our law's prediction, a bias of
  b = 0.22. That is the direction of the known X-ray bias, and it is where weak-lensing
  masses usually sit, though b = 0.22 is at the high end of published estimates. Part of it
  may be a genuine over-prediction.

**The earlier "η ≠ 1" was an artefact of measuring against the MOND law.** The new law needs
no slip.

### 8.6 The Bullet Cluster: partly solved, and the top open problem

This is a toy geometry, not a fit to the maps: 150 Myr after a 4,700 km/s collision, with
gas outweighing stars 7:1.

**Idea:** companions keep the velocity of the matter that emitted them. When the gas is
stopped, its earlier companion flies on with the galaxies, and only its fresh companion
(spread to u·t) sits on the gas.

| | Lensing peak (galaxies −200, gas −80 kpc) | Share of the way to the galaxies |
|---|---:|---:|
| No memory (companion follows the matter now) | −70 | 0% (on the gas) |
| Fresh companion spread 31–63 kpc | −170 | **75%** |
| 94 kpc | −150 | 58% |
| 134 kpc (u × 150 Myr) | −130 | 42% |
| ≥ 268 kpc | −90 to −70 | back on the gas |

* **The mechanism moves lensing off the gas toward the galaxies**, which the famous
  argument says no force law can do, because this one has memory.
* **But the stacked measurement pushes back.** Across 72 collisions, lensing sits within
  5.8 ± 8.2 kpc of the stars (Harvey et al. 2015, Science 347, 1462). Our model predicts
  lensing drifting back to the gas within about 0.3 Gyr, which those data don't show.
* **What would have to change:** the gas must rebuild its companion far more slowly than
  u = 874 km/s allows. For example, shocked gas might emit less for a while. This is now
  the top open problem.

## 9. Open, and next (updated)

1. **Mergers.** The companion rebuild time after a collision (§8.6).
2. **Derive `u` and `g_d`.** They are still measured.
3. **Solar mass loss** of 1.4 × 10⁻¹⁴ per year, against planetary ranging.
4. **Spectroscopic IMFs** for the six SLACS lenses.
5. **The KiDS offset profile** against radius.
6. **The full action**, written relativistically, for lensing beyond the weak field.

## 10. Round 3, 23 September 2026: colliding clusters, and the law locked in

Scripts: `code/dicke_toy.py`, `code/run_v3.py`, `code/bullet_v3.py`, `code/collisions_v3.py`,
`code/kids_v3.py`, and `code/lenses_t35.py --constants`. Their outputs are in `run-dicke/`,
`run-v3/`, `run-bullet-v3/`, `run-collisions-v3/`, `run-kids-v3/` and `run-lenses-v3/`.

Colliding clusters are the classic argument for dark matter, so passing them was required
before this law could stand as an alternative. Two changes, each from first principles, did
it. Neither changes the law for round, relaxed systems.

### 10.1 Change 1: colliding matter is cold to the companion (Dicke narrowing)

The heat term rests on Doppler scrambling (§8.1): a moving emitter's companion phase drifts
by `k_c x(t)`, the companion wavenumber times the emitter's displacement. Laboratory
spectroscopy has a known exception. An emitter that changes direction many times before
it travels one wavelength drifts only diffusively, and its Doppler effect switches off
(R. H. Dicke, Phys. Rev. 89, 472, 1953). The confined-emitter version is the recoil-free
Mössbauer line.

`dicke_toy.py` checks this for the companion. Units are k_c = 1 and σ = 1, with phase locking
at rate Γ.

* **Coherence time.** Free streaming gives 1.40, against 1.41 predicted. Colliding with
  k_c ℓ = 0.1 gives 10.0 (theory 10.1). With k_c ℓ = 0.01 it gives 102 (theory 100).
* **Steady-state phase variance.** It matches `k_c²σ²τ_c / (Γ(1 + Γτ_c))` at every τ_c,
  for example 0.0200 against 0.0200 at τ_c = 0.001. Free streaming gives `k_c²σ²/Γ²`,
  proportional to σ². That is the heat weight's form, `k = 3σ²/u²` with `u = √3 Γ/k_c`.
  Collisions multiply it by `Γτ_c/(1 + Γτ_c)`.
* **The summed field**, the round-2 step-2 set-up, at one temperature:
  * free-streaming emitters give 0.998 of the scalar sum;
  * colliding emitters (k_c ℓ = 0.001) give 1.20 × the squared vector sum, like cold matter
    (11.2 × for free streamers).

**Rule:** stars and galaxies, which never collide, carry the heat. Gas and plasma, whose
particles collide or gyrate far faster, add coherently like cold matter. They still count
in full as ordinary mass. This needs the hot cluster gas to scatter its ions faster than
about the locking time (0.1–1 Myr is the natural window). That is likely in magnetised
cluster plasma but not established; it is flagged as an assumption.

### 10.2 The refit: stars now carry the cluster heat

The stars' random speeds come from the isotropic Jeans equation in the measured X-ray
acceleration. The stellar mass profiles are the X-COP release's (Ghizzardi et al. 2021,
seven clusters). The other five use the median star-to-gas ratio. Stars dominate cluster
centres: M*/M_gas is 2–11 at 0.02 R500, about 1 at 0.05, 0.3–0.5 at 0.1, and 0.02–0.07 at
R500.

| | Round 3 | Round 1 | MOND | Newton | NFW dark matter |
|---|---:|---:|---:|---:|---:|
| Constants | a = 6.561 × 10⁻¹¹ m/s², g_d = 2.262 × 10⁻¹⁰ m/s² (λ = 3.448), **u = 197.4 km/s** | u = 874 | 1 | 0 | 2 per object |
| 149 galaxies, velocity error | **15.85 km/s** | 15.93 | 16.13 | 45.58 | 7.52 |
| 25 bulge-dominated galaxies | **29.18** | 29.48 | 30.35 | | |
| 12 clusters, rms ln(M_HSE/M_pred) | **0.227** | 0.223 | 1.062 | 2.211 | 0.101 |

The 25 bulge-dominated galaxies were the worry: their bulges are now 19× hotter per
kilogram. They improve.

Cluster residuals by radius (0.1 → 1 R500) are +0.06, +0.11, +0.08, +0.03, −0.07, −0.19.
That is flatter than round 1's +0.14 → −0.25.

Variants and robustness:

| Variant | u (km/s) | rms |
|---|---:|---:|
| Stars moving at the gas's speed | 224 | 0.202 |
| Seven clusters with measured stellar profiles only | — | 0.200 |
| **Stars' speeds from the law's own gravity (no X-ray input)** | 216 | **0.329** |
| Held-out, all 924 half-splits | 197 (178–218) | 0.244 (0.197–0.285) |

A free gas share of the heat scores 0.217–0.228 at every weight. The clusters alone cannot
tell gas from stars; colliding clusters can (§10.4).

### 10.3 Change 2: the companion pulls along its net flow

Energy fluxes add as vectors even when the waves are out of step. So the companion's net
flow is along `g_N + g_hot`, where `g_hot = G ∫ k ρ_free (x′ − x)/|x′ − x|³` is the
heat-weighted Newtonian field of the free-streaming matter. The extra pull takes that
direction, with the two parts mixed in proportion to their sizes. It fades smoothly to
zero where the parts cancel:

```
h = g_N + e^(−|g_N|/g_d) √(a(|g_N| + S)) · (g_N + g_hot) / (|g_N| + |g_hot|)
∇²Φ = −∇·h ,   g = −∇Φ   (the pull on matter and light)
S = G ∫ k ρ_free / d² ,   k = 3σ²/u² for free-streaming matter, 0 for colliding matter
```

For spherical systems, and for disks with bulges, `g_hot ∥ g_N`. The pull is then exactly
the fitted law, so none of §10.2 changes.

A first version used the pure unit vector of the flow. It flips abruptly where the two
parts oppose, which would plant sheets of negative lensing. The mixed form is continuous;
the minimum κ in every Bullet map is positive.

### 10.4 The Bullet Cluster, on its published numbers

**Data.** Clowe et al. 2006 (Table 2), 100 kpc apertures:

| Position | Gas (10¹² M☉) | Stars (10¹² M☉) | κ̄ |
|---|---:|---:|---:|
| Main BCG | 5.5 | 0.54 | 0.36 ± 0.06 |
| Main plasma | 6.6 | 0.23 | 0.05 ± 0.06 (excess) |
| Sub BCG | 2.7 | 0.58 | 0.20 ± 0.05 |
| Sub plasma | 5.8 | 0.12 | 0.02 ± 0.06 (excess) |

Galaxy speeds are from Barrena et al. 2002: main 1249 (+109/−100) km/s from 71 galaxies;
subcluster 212 (+67/−52) km/s from 7.

**Model.**
* Two β-model gas clouds and two NFW-shaped stellar concentrations, fitted to the eight
  aperture masses. The bullet core has a 30 kpc floor, because the fit otherwise shrinks it
  to 5 kpc and plants a false lensing bump.
* The 3D field equation on a 192³ grid at 15 kpc resolution.
* κ with sources at z = 1 (Σ_crit = 2.84 × 10⁹ M☉/kpc²).
* The same two-circular-profile decomposition as Clowe et al.

**Prediction (A):** stars' speeds from each cluster's own gravity under our law, with the
subcluster's pre-collision gas restored at 1:8 of the main cluster.

| | Main BCG | Sub BCG | Main plasma | Sub plasma | Main peak | Sub peak |
|---|---:|---:|---:|---:|---|---|
| Observed | 0.36 ± 0.06 | 0.20 ± 0.05 | 0.05 ± 0.06 | 0.02 ± 0.06 | on galaxies | 43 kpc from BCG, toward gas |
| **Ours (A)** | 0.51 | 0.07 | 0.04 | 0.04 | **8 kpc from BCG** | **34 kpc from BCG, toward gas** |
| Ours, stars at M/L 1–1.5 | 0.31–0.41 | 0.07 | 0.05 | 0.04 | on galaxies | weak |
| Rounds 1–2 direction rule | 0.35 | 0.13 | **0.21** | 0.07 | **on the gas** | **on the gas** |
| Cold limit (no heat) | 0.14 | 0.07 | 0.06 | 0.04 | **on the gas** | **on the gas** |

* The law's own star speeds for the main cluster rise from 620 km/s at 20 kpc to 1,140 km/s
  at 400 kpc, against the measured 1,249 ± 100.
* **The pattern is reproduced:** lensing on the galaxies, nothing extra on the gas. Both
  changes are needed; each old rule puts the peaks on the gas.
* **The strengths are 2.5σ off in opposite directions.** Lighter stars, within Clowe's stated
  range, fix the main cluster but not the subcluster. The subcluster needs its stars to move
  at about 500 km/s or more. Our law gives 270–440 km/s; its 7 measured galaxies give
  212 ± 60. **This is the open item.**

### 10.5 The 72-collision stack

We used the Bullet as a template and varied the collision stage (gas lagging by 40–300 kpc)
and the mass ratio: 20 substructures in all.
* The lensing peak stays with the galaxies in every case, 4–22 kpc from them at 18 kpc
  resolution. The median fraction of the way to the gas is 0.06.
* Harvey et al. 2015 measure 5.8 ± 8.2 kpc (β = −0.04 ± 0.07).
* The drift back to the gas that sank round 2 is gone. The heat rides with the stars at
  every stage, with no memory needed.

### 10.6 Momentum: the companion carries it

**The round-2 action gave hot matter a reaction force** proportional to k. With stars as the
hot component (about 5% of the mass), that reaction pushes the stars outward by a median
32% of the gravity at 0.1 R500 and 73% at 0.3–0.5 R500. Cluster galaxies would then orbit as
if clusters held a third of their mass. Galaxy kinematics agree with lensing and X-ray
masses (for example Pizzuti et al. 2016), so that is excluded.

**The consistent alternative is how every streaming field behaves.** The companion carries
momentum, and an emitter recoils only if it emits unevenly. Gas and galaxies then feel the
same pull, so η = 1. The price is a small self-force on lopsided hot-and-cold systems, paid
from the companion's momentum flux `P/u = (a/2) Σ(1 + k)m`. For the Bullet it is 12% of the
mean pull, and 13% of that budget.

A dynamical, relativistic field theory for the companion is still to be written.

### 10.7 Everything else, re-run under round 3

* **KiDS early/late offset.** The early types' own stars (150–200 km/s, collisionless) give
  **+0.17 to +0.27 dex**, roughly constant beyond 100 kpc. Brouwer et al. 2021 measure 0.17
  (Sérsic split) and 0.27 (colour split). Hot haloes now count only as mass. They add
  +0.02–0.07 dex for 0.3–1 M*, so the prediction is that such haloes are modest.
* **SLACS**, standard bookkeeping:
  * stars need **1.05–1.35 × Salpeter** (round 2: 1.23–1.57);
  * the lensing-minus-kinematics gap is **−0.017 ± 0.024 dex**, so there is no slip.
* **Solar System.** Zero extra pull, from the release factor. The Sun is collisional plasma,
  so it has no heat term. Emission is ℓ = a u/2 = 6.5 × 10⁻⁶ W/kg, and the Sun's extra
  mass loss is **2.3 × 10⁻¹⁵ per year**, 6 times smaller than in round 2.
* **Guard.** The cold isolated mass takes a MOND form, `ν = 1 + e^(−y/λ)/√y`. The closest
  published function is RAR-exponential at 0.030 dex (simple 0.035, standard 0.084). Random
  speeds of 0–1000 km/s at the same Newtonian pull change the prediction by 0.93 dex, and
  SPARC points sit 0.022 dex beyond any local function of g_N. It is **not MOND, not Newton,
  and has no dark component or per-object parameters.**
* **MOND as the cold limit (§8.1)** is unchanged, with a₀ = 2ℓ/u and the new constants.

## 11. Open, and next (round 3)

1. **The subcluster's lensing strength.** It is 2.5σ low: its stars' speeds, its stellar mass,
   or the time since the collision.
2. **A field theory for the companion**, with momentum flux and a relativistic form, so
   lensing is derived rather than assumed equal to dynamics.
3. **Microphysics.** Derive u and g_d, and pin the companion wavelength and locking time
   that separate gas from stars.
4. **Wide binaries.** Near the Sun the Galaxy's pull is close to the release scale g_d, so the
   release factor is partly on. That makes our external-field boost smaller than MOND's, and
   the predicted wide-binary signal weaker. It needs a proper external-field calculation. The
   data are contested (Chae 2024; Banik et al. 2024).
5. **Checks others can make now:**
   * spectroscopic IMFs of the six SLACS lenses (1.05–1.35 × Salpeter);
   * KiDS offset flat with radius beyond 100 kpc;
   * galaxy and gas dynamics equal in relaxed clusters;
   * lensing always with the galaxies in mergers.

## 12. Round 4, 23 September 2026: the subcluster, and the companion's memory

Scripts: `code/bullet_v4.py`, `code/collisions_v4.py`, `code/stream_tidal_v4.py` and
`code/bullet_speeds_v4.py`. Their outputs are in `run-bullet-v4/`, `run-collisions-v4/`,
`run-stream-tidal-v4/` and `run-bullet-speeds-v4/`. `code/bullet_v3.py` gained three options:
`v_rel` (stream heat), `star_scale`, and `ghost` (a first version of the memory). Its defaults
reproduce round 3 exactly.

Round 3 left the Bullet's subcluster lensing at 0.07, against 0.20 ± 0.05. Three physical
effects were examined. **No new rule was added.**

### 12.1 Streams passing through each other (stream heat)

The heat weight is the mean-square speed of free-streaming matter about its *local* mean
motion. Where two galaxy streams overlap at relative speed V, that adds
`ρ_m ρ_s V² / (ρ_m + ρ_s)` to `ρ⟨|v − v̄|²⟩`. Streaming along one axis counts V², not 3V².

* **Galilean invariant, and gas-safe.** Colliding gas has no co-located streams.
* **Evaluated on the instantaneous field** (no memory), it raises the subcluster's κ from
  0.073 to 0.100, 0.109 and 0.162 at V = 2,600, 3,000 and 4,700 km/s.
* **But the heat sits where the two stellar densities are comparable,** between the clusters.
  So it smears the subcluster's peak away at 2,600–3,000 km/s, and at 4,700 km/s pulls it to
  45 kpc from the gas.
* **With the slow companion (§12.3), the heat that switched on in the collision has spread
  only about u × t ≈ 30–60 kpc.** Today it is negligible for the Bullet; it builds up over
  the next billion years.

### 12.2 Tidal shaking

As the subcluster crossed the main cluster's core, the main cluster's tides kicked its stars.
With our law's pull for the main cluster, and a 150 kpc miss distance at 2,600–4,700 km/s,
the kicks are 170–300 km/s at 50 kpc and 330–600 km/s at 100 kpc from the subcluster's centre.
A head-on pass diverges, because the pull stays finite at a cluster's centre.

**A fast pass gives every star at a given place the same kick.** That is ordered, converging
motion, not a spread of speeds, so it becomes heat only as the stars phase-mix, over roughly
an orbital time (~0.5–1 Gyr). It adds little to the Bullet today, and it predicts stronger
lensing around the smaller clump in older collisions.

### 12.3 The companion's memory: the decisive effect

Two consequences of the existing rules:

1. **The companion is slow.** At u = 197 km/s it takes about 500 Myr to travel 100 kpc. The
   field around each Bullet component today was therefore mostly emitted *before* the
   collision, about 150 Myr ago.
2. **The companion keeps its emitter's velocity (emission, not a medium).** This is required
   by Galilean invariance, because galaxies move at hundreds of km/s, well above u. If the
   companion moved at u through a fixed medium, a galaxy's gravity would depend on how fast it
   moves through space, yet all galaxies obey one rotation law.

**So the field around each cluster now is that of the settled cluster it was before the
collision, still moving with the pre-collision motion.** The collisionless galaxies kept that
motion, so the old field rides with them. Only inside a fresh sphere of radius u × t (30 kpc
for t = 150 Myr) has it been rebuilt around the stopped gas. The ordinary pull and the release
factor follow the matter where it is now.

**Model.** `bullet_v4.py`:
* **Main cluster before the collision:** today's gas and stars, centred on its galaxies.
* **Subcluster before the collision:** the bullet, its lost atmosphere (β-model, 150 kpc core)
  and its galaxies (today's core plus satellites), centred on its galaxies.
  * Baryons: 1/8 of the main's, from merger reconstructions (1:6–1:10).
  * Stars: 7% of the gas, typical of X-COP clusters.
  * Star speeds: from each pre-collision cluster's own gravity under our law.
* **Fresh spheres:** swapped in with a distance-masked vector kernel.
* **Nothing is fitted to the lensing.**

| Case | Main (0.36 ± 0.06) | Sub (0.20 ± 0.05) | Main gas (0.05 ± 0.06) | Sub gas (0.02 ± 0.06) | Peaks from galaxies (kpc) |
|---|---:|---:|---:|---:|---|
| A. No memory, pre-collision star speeds | 0.51 | 0.09 | 0.05 | 0.05 | 12 / 12 |
| **B. Memory: sub 1/8 of main, stars 7% of gas** | 0.52 | **0.175** | **0.04** | **0.04** | 12 / 17 |
| C. Memory, sub 1/10 | 0.53 | 0.13 | 0.03 | 0.04 | 10 / 14 |
| D. Memory, sub 1/6 | 0.51 | 0.25 | 0.05 | 0.04 | 18 / 20 |
| E / F. Memory, stars 5% / 10% of gas | 0.53 / 0.52 | 0.13 / 0.25 | 0.03 / 0.05 | 0.04 / 0.04 | 10 / 14; 17 / 19 |
| G. Memory, compact atmosphere (100 kpc core) | 0.52 | 0.19 | 0.04 | 0.03 | 12 / 17 |
| H. Memory, gas stopped 300 Myr ago (60 kpc sphere) | 0.52 | 0.17 | 0.04 | 0.04 | 13 / 17 |
| **I. Memory, star M/L 1.5 for both** | **0.43** | **0.19** | **0.04** | **0.04** | **18 / 20** |
| **J. Memory, star M/L 1.0 for both** | **0.33** | **0.21** | **0.05** | **0.04** | **32 / 19** |

Peaks are refined below the 15 kpc pixel. All of them sit on the gas side of their galaxies; the
observed subcluster peak is 43 kpc from its brightest galaxy, also on the gas side.

* **Faster stars alone do little.** Giving the subcluster's stars their pre-collision speeds
  without the memory (case A) lifts it only from 0.07 to 0.09.
* **The subcluster is solved by the memory.** It lifts the subcluster to 0.13–0.25 across standard
  pre-collision assumptions (0.175 at the central choice), against 0.20 ± 0.05.
* **The main cluster matches with lighter stars.** At mass-to-light 1–1.5 for both clusters,
  inside Clowe et al.'s 0.5–3 (their masses are upper limits because foreground galaxies were
  not removed), **all four measurements agree within about 1σ** (largest: +1.2σ for the main
  cluster at 1.5) and both peaks sit on the galaxies.
* **The subcluster's stars before the collision moved at 460–610 km/s** at 50–200 kpc under our
  law (520–690 km/s at 1:6; M/L 1.5). That is a prediction for a spectroscopic survey with more than the 7
  galaxies measured so far (212 +67/−52 km/s, a group picked out by its narrow velocity spread).
* **Independent support.** Barrena et al. (2002, A&A 386, 816), who measured those 7 galaxies,
  argued that the subcluster "is in fact the remnant core of a moderately massive cluster,
  stripped by the collision". Its X-ray temperature and luminosity correspond to a velocity
  dispersion of about 700 km/s and a pre-merger mass ratio of about 1:6. That is the cluster our
  memory model needs.

### 12.4 The 72-collision stack, with memory

Same 20 pieces as round 3 (the Bullet as a template; gas lagging its galaxies by 40–300 kpc;
subcluster mass ×1 and ×3), plus a reference with the gas left on its galaxies. The fresh
sphere is u × t with t = lag / 1,240 km/s (the Bullet's separation speed), so 6–48 kpc. Stars
at M/L 1.5. Peaks are refined below the 18 kpc pixel.

* **The lensing does not follow the gas.** Moving the gas 40–300 kpc from its galaxies moves
  the lensing peak by −2 to +10 kpc (median +3 kpc). As a fraction of the gas offset that is
  β = 0.03 (median; range −0.01 to 0.12). Straight-line fits over all six stages give
  −0.012 to +0.026.
* **Harvey et al. 2015 measure 5.8 ± 8.2 kpc (β = −0.04 ± 0.07).** We agree within 1σ.
* **Why it holds.** The fresh sphere grows at u = 197 km/s and never catches up with gas
  separating at about 1,000 km/s. This is exactly the change round 2 asked for (§8.6): its
  memory failed because the companion was then 874 km/s, so the stopped gas rebuilt its field
  within 0.3 Gyr. Round 3's refit lowered u for an unrelated reason (stars carry the heat), and
  that makes the memory safe.
* **A fixed offset from superposition, the same with the gas left in place.** In the summed map
  each peak sits 7–41 kpc toward the other cluster even at lag 0, because the other cluster's
  broad field tilts the map. Harvey et al. fitted each clump separately, which removes it, so
  the comparable number is the shift above. For the record, raw peak distances from the
  galaxies are 5–50 kpc (median 17), and round 3's pixel-level fraction was 0.06.
* **No lensing peak of its own on the main gas.** A weak separate bump on the subcluster's gas
  (from the gas's own mass) appears in 2 of 20 maps. Offsets across the collision axis are
  below 0.3 kpc.

### 12.5 A check on the main cluster's galaxy speeds

Round 3 compared our law's star speed at 400 kpc (1,140 km/s at M/L 2) with Barrena et al.'s
1,249 +109/−100 km/s from 71 main-cluster galaxies. The fair comparison is the line-of-sight
average over the region they sampled (their virial mass uses an aperture of about 1.5 Mpc):

| Star M/L | Inside 0.5 Mpc | Inside 1 Mpc | Inside 1.5 Mpc |
|---|---:|---:|---:|
| 1.0 | 900 km/s | 880 | 840 |
| 1.5 | 990 | 960 | 915 |
| 2.0 | 1,070 | 1,040 | 985 |

The 1.5 Mpc column is pulled down because our model's galaxies are cut off at 1.5 Mpc.

* **The lighter stars that match the lensing give galaxies about 20% slower than measured:**
  2.5–3σ. At the published M/L of 2 it is 1.8σ inside 0.5 Mpc.
* **Barrena et al. find the main cluster's dynamics undisturbed by the collision**, so we do
  not lean on merger-inflated speeds.
* **What would have to change:** the pull on the main cluster's galaxies at 0.3–1.5 Mpc would
  need to be about half as strong again, with the lensing inside 100 kpc unchanged. Our
  main-cluster model was built only from Clowe et al.'s 100 kpc apertures. The next steps
  are a full X-ray gas profile out to 2 Mpc, the same galaxy selection as Barrena et al., and
  a check on galaxy orbits being partly radial in cluster outskirts. The script is
  `code/bullet_speeds_v4.py`.

### 12.6 What round 4 does and does not change

* **Relaxed galaxies and clusters are unchanged.** A system moving steadily has its companion
  simply carried along. Round 3's fits (15.85 km/s; 0.227; u = 197 km/s) stand, and
  `bullet_v3.py`'s new options default off and reproduce round 3 exactly.
* **Formula check.** No new constant and no per-object number beyond the star mass-to-light
  ratio, which is kept inside the published range. The "old" companion is not invisible
  matter. It is fixed by the visible matter's own history and fades on a known timescale, and it
  is absent in settled systems.
* **Memory matters only where matter has recently changed its motion:** collisions, and to a
  lesser degree galaxies on tight orbits.
* **Predictions:**
  * Lensing in a collision reflects the clusters as they were before it, for roughly
    (distance)/u, which is 500 Myr per 100 kpc.
  * Lensing reappears on stopped gas only as a sphere growing at 197 km/s, about 200 kpc per
    Gyr.
  * Stream heat and tidal heat build up around the smaller clump over the following Gyr.

## 13. Open, and next (round 4)

1. **The main cluster's galaxy speeds (§12.5).** About 1,000 km/s from our law with the lighter
   stars, against 1,249 ± 100 measured. Needs a full X-ray gas profile to 2 Mpc and the same
   galaxy selection as Barrena et al.
2. **Test the collision predictions.** A larger spectroscopic sample of the subcluster (460–610
   km/s before the collision; 212 ± 60 from 7 galaxies now). Lensing in older collisions.
3. **A field theory for the companion**, now with its travel time built in (a retarded source),
   momentum flux, and a relativistic form.
4. **Microphysics** (u, g_d, the companion wavelength), **wide binaries**, **cosmology**:
   as in §11.

## 14. Round 5, 23 September 2026: the Bullet main cluster's galaxy speeds, and the lensing masses

Script: `code/bullet_main_v5.py`. Outputs are in `run-bullet-main-v5/` (2.9 Mpc box, 15 kpc
pixels) and `run-bullet-main-v5-wide/` (`--dx 25`, 4.8 Mpc box).
* Lensing strengths in 100 kpc apertures, gas residuals and peaks are quoted from the first.
* Masses inside 250 kpc are quoted from the second, which holds more of the line of sight.
  The narrow box gives 5–7% less.

Round 4 left the main cluster's galaxies about 20% too slow: 981 km/s (line of sight, averaged
over the 1.8 Mpc² Barrena et al. surveyed, equivalent radius 757 kpc), against 1,249 +109/−100.
Round 4's model was built only from Clowe et al.'s two 100 kpc apertures.

### 14.1 What the independent data say about the main cluster

* **Gas: not the problem.** The fitted β model agrees with:
  * the ROSAT β model of the whole system (Ota & Mitsuda 2004, as used by Paraficz et al. 2016):
    1.5 against 1.7 × 10¹⁴ inside 1 Mpc and 3.6 against 3.4 × 10¹⁴ inside 2 Mpc;
  * the gas share of the lensing mass inside 250 kpc (Paraficz et al.: 9 ± 3% of 2.5 × 10¹⁴).
* **Stars inside 340 kpc: consistent.** Barrena et al.'s R-band light, 1.0 × 10¹² L☉ in
  5.4 arcmin² (0.37 Mpc²) excluding the subcluster, is 2–3 × 10¹² M☉ for M/L_R 2–3. The model
  has 2.1–2.7 × 10¹² (M/L_I 1.5–2).
* **Stars beyond 340 kpc: never measured, and too few in the model.** In the 7 X-COP clusters
  with measured stellar profiles (Ghizzardi et al. 2021), the cumulative star/gas ratio is 0.074
  (median) at 0.5 R500 and 0.062 at R500 (range 0.035–0.073). Round 4's model has 0.045 and 0.035
  at M/L 2, and 0.034 and 0.026 at M/L 1.5. Its NFW was also cut off at 1.5 Mpc.

### 14.2 The galaxy speed

All cases use M/L 1.5 for the central stars unless stated. "Outer stars" is an extended
component (NFW, scale 800 kpc, cut at 3 Mpc).

| Case | Speed in Barrena et al.'s region | Star/gas at 0.5 / 1 R500 | Stars projected inside 340 kpc |
|---|---:|---:|---:|
| A. Round 4 (stars cut off at 1.5 Mpc) | 981 km/s (+2.7σ) | 0.034 / 0.026 | 2.1 × 10¹² |
| B. The same stars, not cut off | **1,150 (+1.0σ)** | 0.034 / 0.026 | 2.1 |
| C. Outskirts at the X-COP median ratio (fitted scale 300 kpc) | 1,310 (−0.6σ) | 0.080 / 0.058 | 5.1 (too much) |
| D. Outer stars from the speed: 8.2 × 10¹² | 1,249 | 0.053 / 0.043 | 3.2 |
| E / F. The same for 1,149 / 1,358 km/s: 4.9 / 12 × 10¹² | 1,149 / 1,358 | 0.045–0.061 / 0.036–0.051 | 2.7 / 3.7 |
| G. D with orbits partly radial (anisotropy 0.3) | 1,313 | as D | as D |
| J. Published M/L 2, outer stars from the speed: 6.2 × 10¹² | 1,249 | 0.059 / 0.048 | 3.6 |

* **Half of round 4's gap was a modelling artefact.** The model's galaxies stopped at 1.5 Mpc, so
  galaxies near that edge had almost no random motion, and the heat of the outskirts was
  missing. Letting the same profile continue gives 1,150 km/s.
* **The rest needs a normal amount of stars in the outskirts.** It takes a star/gas ratio of
  0.043–0.048 at R500: inside the X-COP range (0.035–0.073) and below its median. Orbits that
  are partly radial, common in cluster outskirts, would need fewer. The light inside 340 kpc
  stays close to Barrena et al.'s 2–3 × 10¹².
* **The speed is a sensitive measure of the stars in our law.** The heat term makes σ grow
  roughly as the square root of the hot stellar mass, so the galaxy speed is effectively a
  stellar census of the outskirts.

### 14.3 The lensing, re-examined: Clowe et al.'s kappas are floors

Clowe et al. (2006) state that their reconstruction "does, however, systematically underestimate
κ in the cores of massive clusters ... our measurements of κ in the peaks of the components are
only lower bounds." Rounds 3–4 treated 0.36 ± 0.06 and 0.20 ± 0.05 as measurements. They are
floors, so:
* round 4's "all four within about 1σ" should read "above both floors, with the gas residuals
  as measured";
* lowering the star masses to M/L 1–1.5 was never required.

The calibrated strengths are the projected masses inside 250 kpc of each brightest galaxy, from
strong (and weak) lensing:
* Bradač et al. 2006: main 2.8 ± 0.2, sub 2.3 ± 0.2 × 10¹⁴ M☉;
* Paraficz et al. 2016 (strong lensing, 14 multiple-image systems): main 2.5 ± 0.1, sub
  2.0 ± 0.2 × 10¹⁴ M☉.

| Case | κ main (floor 0.36) | κ sub (floor 0.20) | Gas residuals (0.05, 0.02 ± 0.06) | Main inside 250 kpc | Sub inside 250 kpc | Peaks from galaxies |
|---|---:|---:|---|---:|---:|---|
| A. Round 4 (M/L 1.5, sub 1:8) | 0.43 | 0.19 | 0.04, 0.04 | 1.66 | 0.92 | 19 / 20 kpc |
| D. Outer stars from the speed (M/L 1.5, sub 1:8) | 0.63 | 0.14 | 0.03, 0.04 | 2.43 | 1.05 | 11 / 19 |
| L. Published stars (M/L 2), outer stars, sub 1:8 | 0.67 | 0.14 | 0.03, 0.04 | 2.54 | 1.04 | 9 / 16 |
| L. The same, sub 1:4 | 0.64 | 0.36 | 0.06, 0.06 | 2.49 | 1.52 | 21 / 21 |
| **L. The same, sub 1:3** | **0.61** | **0.49** | **0.08, 0.10** | **2.42** | **1.93** | **31 / 17** |
| K. M/L 1.5, sub 1:2 | 0.44 | 0.72 | 0.06, 0.10 | 2.03 | 2.61 | 52 / 10 |

Masses are in 10¹⁴ M☉.

* **The main cluster now fits both its galaxy speed and its lensing mass.** They pull the same
  way: both needed more pull in the outskirts than round 4 had. With the published star masses
  and the outer stars the speed requires, the main cluster has 2.42–2.54 × 10¹⁴ inside 250 kpc
  (Paraficz et al. 2.5 ± 0.1; Bradač et al. 2.8 ± 0.2).
* **The subcluster's lensing mass says it was a bigger cluster than we assumed.** At 1:8 our law
  gives it 1.0 × 10¹⁴ inside 250 kpc, about half of the measured 2.0–2.3. It matches at about
  **1:3 in visible matter** (1.93 × 10¹⁴).
  * Its pre-collision galaxies then moved at 680, 810 and 900 km/s at 50, 100 and 200 kpc.
    Barrena et al. inferred about 700 km/s from its X-ray temperature and luminosity.
  * Dark-matter merger reconstructions use 1:6–1:10 in total mass; the lensing masses themselves
    are in the ratio 0.8 inside 250 kpc.
* **Costs of the bigger subcluster.**
  * The gas residuals rise to 0.08 and 0.10 (+0.5σ and +1.3σ).
  * The main peak moves to 31 kpc from its galaxies, because the subcluster's broad field tilts
    the summed map toward it (the superposition effect of §12.4).
* **Nothing is fitted to the lensing except the subcluster's size.** One number, the outer
  stars, is inferred from one measurement, the galaxy speed. Another, the subcluster's size, is
  inferred from the subcluster's lensing mass. Both become predictions.

### 14.4 Predictions from round 5

1. **The main cluster's outskirts** hold stars at about 5% of the gas mass at R500 (4.3–4.8%),
   within the X-COP range. Wide-field photometry out to 1.5 Mpc can check this.
2. **The subcluster was about a third of the main cluster before the collision.**
   * About 7 × 10¹² M☉ of stars came with it. They are now outside its compact core; tidal
     stripping spreads them around it, where most would be counted as main-cluster members.
   * About 1.1 × 10¹⁴ M☉ of its gas was stripped into the main cluster's gas.
   * Its original galaxies moved at about 800 km/s. A phase-space analysis of a larger
     spectroscopic sample can separate them, since they carry the subcluster's 616 km/s
     line-of-sight offset.
3. Unchanged from round 4: lensing stays with the galaxies; it returns around stopped gas at
   about 200 kpc per Gyr; older collisions show extra lensing around the smaller clump.

### 14.5 What round 5 does and does not change

* **The law and its three constants are unchanged.** Relaxed galaxies and clusters are
  unaffected.
* **Corrected:**
  * round 4's reading of Clowe et al.'s κ (floors, not measurements);
  * round 4's "subcluster solved" (only against the floor; against the lensing mass it needs a
    bigger pre-collision subcluster);
  * the round-4 idea that the stars had to be lighter.
* **Formula check (RULES.md §9):** no new formula. The inferred numbers are properties of
  visible matter, checkable by counting stars and galaxies. None is an invisible or free
  component.

## 15. Open, and next (round 5)

1. **Test the subcluster's size.** Look for the galaxies and stars that came with it, which
   should be about 7 × 10¹² M☉ moving with the subcluster. Also check whether a 1:3 collision
   reproduces the shock and X-ray morphology; dark-matter reconstructions use 1:6–1:10.
2. **A stellar census of the main cluster's outskirts** (predicted: star/gas about 0.045 at
   R500).
3. **A field theory for the companion**, with its travel time, momentum flux and a relativistic
   form.
4. **Microphysics, wide binaries, cosmology:** as in §13.

## 16. Round 6, 23 September 2026: three tests of round 5, and the wide-binary prediction

Scripts: `code/bullet_members_v6.py`, `code/bullet_light_v6.py` and `code/wide_binaries_v6.py`.
Outputs are in `run-bullet-members-v6/`, `run-bullet-light-v6/` and `run-wide-binaries-v6/`.
New data files in `data/`:
* `barrena2002_table1.json`: the 78 cluster members with positions and velocities, transcribed
  from Barrena et al. 2002, Table 1 (astro-ph/0202323);
* `ls_dr10_bullet.csv.gz`: 26,288 Legacy Survey DR10 galaxies with z < 22.5 and photometric
  redshifts, in a 1.6° × 0.9° box around the Bullet, from the NOIRLab Astro Data Lab. The query
  is in the script.

### 16.1 The subcluster's lost galaxies are not in the velocity sample

If the subcluster held a third of the main cluster before the collision (round 5), about
7 × 10¹² M☉ of its galaxies lie beyond its compact core but still move with it, at +616 km/s.

* **The core is recovered.** The 7 galaxies within 250 kpc of its brightest galaxy that are
  within 600 km/s of +616 are exactly Barrena et al.'s KMM group: mean +622, spread 205 km/s
  (theirs: +616 and 212).
* **A mixture test on the other 71 galaxies.** The model's projected star densities set each
  galaxy's prior odds of being lost. At 1:3 about 12–15 of them should be lost galaxies.
* **Result: every variant prefers none.** ln L(1:3) − ln L(1:8) is −3.9 to −0.7, where +0.2 to
  +0.7 is expected if 1:3 were true.
  * For a compact or cool lost population (scale 150–300 kpc, spread 600–850 km/s), the
    observed value falls below every simulated 1:3 sample: about 2.5–3σ.
  * For a population spread 600 kpc with a 1,200 km/s spread, it falls below 92%: about 1.4σ.
* **Forecast.** About 450–1,900 velocities would separate 1:3 from 1:8 at 3σ.

### 16.2 A star count from the Legacy Survey

Cluster galaxies are those with z < 21.5 and a photometric redshift within 0.04 (also tested:
0.06) of 0.296. The field level is measured at 4–6 Mpc. Stars are taken at M/L_I = 2 from the
z band.

**The main cluster's outskirts: confirmed.** Light inside R, eastern half doubled, relative to
all the light inside 250 kpc:

| R (kpc) | 500 | 750 | 1,000 | 1,500 | 2,000 | 2,500 | 3,000 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Data (window ±0.04) | 1.57 | 2.39 | 3.10 | 4.05 | 4.99 | 5.24 | 5.23 |
| Data (window ±0.06) | 1.76 | 2.52 | 3.19 | 4.61 | 5.74 | 6.58 | 6.40 |
| Round 4 model (stars cut at 1.5 Mpc) | 2.19 | 3.18 | 3.94 | 4.76 | 4.76 | 4.76 | 4.76 |
| Round 5 model (with the outer stars) | 2.26 | 3.35 | 4.24 | 5.39 | 5.78 | 6.05 | 6.19 |

* **The light keeps rising beyond 1.5 Mpc.** 23–28% of the total out to 3 Mpc lies beyond it.
  Round 5's model has 13% there; round 4's had none.
* **The total out to 3 Mpc is 85–103% of round 5's model**, which was inferred from the galaxy
  speed alone.
* **Inside 1.5 Mpc the observed light is more concentrated than either model.**
* **The absolute scale** from these magnitudes is about 1.8 times Clowe et al.'s inside 100 kpc,
  so we compare shapes.

**The subcluster: little light beyond its core.** We measured the excess light near the
subcluster, after removing the field and the main cluster's own light at the same distance
from the main cluster (taken from the eastern side). As a fraction of the main cluster's light
inside 250 kpc:

| Inside 250 kpc of the subcluster | Excess light ÷ main cluster's light inside 250 kpc |
|---|---:|
| Data (both windows) | 0.46 |
| Model: its core only | 0.32 |
| Model: core + lost galaxies still in place, 1:8 | 0.68 |
| Model: core + lost galaxies still in place, 1:3 | 1.55 |

Inside 500 kpc the excess is no larger than inside 250 kpc. **Neither the velocities nor the
light show the stars that 1:3 needs**, unless they have spread beyond 500 kpc, moved faster than
1,200 km/s, or are diffuse light too faint for a galaxy catalogue.

### 16.3 The subcluster's lensing mass: open again

* **In our law the subcluster's lensing comes mainly from its hot stars.** At 1:3, lowering its
  pre-collision stars from 7% to 3% or 1.5% of its gas drops its mass inside 250 kpc from 1.93
  to 1.23 or 1.05 × 10¹⁴. The 1:2 case with 3% stars gives 1.59. More gas alone does not replace
  the stars.
* **So the measured 2.0–2.3 × 10¹⁴ against our 1.0 at 1:8 is again the top open problem.**
* **What would have to change**, from the most to the least promising within our rules:
  1. **Heat gained in the passage, carried with the subcluster.** While its galaxies crossed the
     main cluster at 3,000–4,500 km/s, their random speed relative to the local mean was
     several times their own, so their companion was fed with a heat weight 10 times larger.
     That companion keeps their velocity (the memory), so it now surrounds the subcluster.
     Round 4 evaluated the stream heat only where the streams overlap *now*; the emission
     history has not been modelled.
  2. **Diffuse stars stripped from the subcluster's galaxies**, moving with it. A galaxy
     catalogue cannot see them; deep imaging of the diffuse light can.
  3. **A wider or faster lost population** (allowed at about 1.4σ by the velocities).

### 16.4 Wide binaries near the Sun: the prediction

* **Setup.** A binary is cold, so our law is h = g_N + e^(−|g_N|/g_d) √(a|g_N|) ĝ_N, with g_N
  the binary's pull plus the Galaxy's. The Galaxy's own Newtonian pull at the Sun is
  1.58 × 10⁻¹⁰ m/s², from a 230 km/s rotation at 8.2 kpc. That holds back half of the companion.
* **Method.** The extra ("as if") density has a closed form. Averaged over the binary's
  orientation, Gauss's law gives the mean inward pull as G(M + M_extra(<r))/r² exactly. With the
  Galaxy removed, the script reproduces the isolated law to 3 decimals.

Boost of the pull for a binary of 1.5 suns (velocities scale as its square root):

| Separation (thousand AU) | 1 | 2 | 3 | 5 | 7 | 10 | 15 | 30 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **Our law** | 1.000 | 1.000 | 1.003 | 1.079 | 1.172 | 1.189 | 1.191 | 1.192 |
| Our law, Galaxy's pull ±20% | | | | 1.07–1.08 | 1.14–1.21 | 1.14–1.26 | 1.14–1.26 | 1.14–1.27 |
| MOND (QUMOND, simple function) | 1.013 | 1.051 | 1.109 | 1.255 | 1.384 | 1.428 | 1.433 | 1.434 |
| Newton | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |

* **Our law: 19% more pull (about 9% faster orbits) beyond about 7,000 AU.** That is less than
  half of MOND's 43%, and it switches on later, past 3,000 AU where MOND already gives 11%.
  * Binaries of 1–2 suns give the same plateau.
  * Heat from the local disk (S = 0.1 g_e) adds 0.005.
* **The data are contested.** Chae (2023, 2024) reports about 1.4 times Newton's pull at the
  lowest accelerations. Banik et al. (2024) report Newton and exclude MOND. Our law agrees with
  neither as published, so it is a sharp test once the analyses converge.
* **If the Newtonian result holds**, the Galaxy's own heat-weighted flow would have to steer
  the extra pull near the Sun, so that a binary's weak internal field barely changes it. The
  direction rule (§10) already mixes in the heat flow; how much of the Milky Way's flow
  reaches the Sun has not been computed.

### 16.5 Predictions and next steps from round 6

1. **The main cluster's stars continue to 3 Mpc** at the level its galaxy speed requires.
   Confirmed by the star count.
2. **Around the subcluster**, look for diffuse light and for a wide, fast population of
   galaxies carrying +616 km/s. About 450–1,900 velocities would settle 1:3 against 1:8.
3. **Model the heat the subcluster's galaxies gained in the passage**, carried with them.
   This is the leading candidate for its missing lensing.
4. **Wide binaries:** 19% more pull beyond 7,000 AU (orbits 9% faster).

## 17. Round 7, 23 September 2026: the Milky Way and the standard tests, against the literature

The question: have we compared our law with the Milky Way, and with the stars and lenses the
field accepts as the standard tests? Until now the Milky Way entered only through the Sun and
the wide binaries. Round 7 runs the law, unchanged (a, g_d from SPARC, u from X-COP), on every
standard Galactic measurement we could verify, and on the standard lensing cases, side by side
with MOND (QUMOND, simple function, a0 = 1.2 × 10⁻¹⁰ m/s²), Newton with visible matter only,
and, where published, the dark-matter model.

Scripts: `code/mw_model.py` (solver), `code/milky_way_v7.py`, `code/mw_dwarfs_v7.py`,
`code/strong_field_v7.py`, `code/lensing_census_v7.py`, `code/transition_check_v7.py`.
Outputs in `run-milky-way-v7/`, `run-mw-dwarfs-v7/`, `run-strong-field-v7/`,
`run-lensing-census-v7/`, `run-transition-check-v7/`.

New data in `data/`:
* `mw_literature_v7.json`, `lensing_literature_v7.json`: every measured value used below, each
  with its paper, page and a verbatim snippet checked by script against the arXiv text;
* `mw_rotation_curves.json` (Eilers et al. 2019; Zhou et al. 2023; Ou et al. 2024; Jiao et
  al. 2023), `mw_dwarfs.json` (McConnachie 2012; Crater II; Antlia 2);
* `brouwer2021_kids/`: the KiDS-1000 lensing RAR data release (excess surface density per bin),
  and `mistele2024/table1_mrt.txt` (lensing circular velocities to 2.5 Mpc).

### 17.1 The solver

* **Visible matter.** McMillan (2017) best fit: an oblate bulge (8.9 × 10⁹ M☉), thin and thick
  stellar disks (3.5 and 1.0 × 10¹⁰), HI and H₂ disks with central holes (1.1 and 0.1 × 10¹⁰),
  plus a stellar halo of 1.4 × 10⁹ M☉ (Deason et al. 2019). Variants:
  * *local census*: the disks rescaled to the solar-neighbourhood census of McKee et al. 2015
    (stars and remnants 33.4, gas 13.7 M☉/pc²; McMillan has 45.8 and 12.2);
  * *BR13*: the short "maximal" stellar disk of Bovy & Rix 2013 (R_d = 2.15 kpc, 38 M☉/pc² at
    8 kpc);
  * *+ corona*: a hot gas corona of 1.3 × 10¹⁰ M☉ within 250 kpc (β model of Miller & Bregman
    2015; gas, so mass only).
* **g_N** on a 360 × 600 cylindrical grid out to 3 Mpc: Hankel transforms for the disks (exact
  for separable disks; against Freeman's razor-thin disk they agree to 10⁻⁴ beyond 12 kpc, and
  the 0.1–0.3% inside is the test disk's 5 pc thickness), the homoeoid formula for the bulge,
  multipoles beyond 60 kpc.
* **Heat** from the bulge (σ by the SPARC rule, 0.65 × its peak speed: 88 km/s) and the stellar
  halo (σ_r, σ_θ, σ_φ = 141, 75, 85 km/s). Disks are cold, as in the SPARC calibration.
* **The law in field-equation form**: ρ_ph = −∇·(h − g_N)/4πG on the grid, summed over rings
  (elliptic integrals) at target points on cell edges. The ring sum reproduces the analytic
  Newtonian pull to 0.05–0.08% in the plane and 0.03% at 1.1 kpc. The companion's extra pull
  stops beyond its reach, u × 13 Gyr = 2.6 Mpc.
* Masses inside radius r come from Gauss's law on spheres (the flux of the true pull equals
  the flux of h); escape speeds from the pull integrated outward along the plane.

### 17.2 Rotation curve (Gaia)

Circular speed (km/s), McMillan's visible matter:

| R (kpc) | 5 | 8.2 | 10 | 15 | 20 | 25 | 27.3 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Eilers et al. 2019 | 227–231 | 228.9 | 225.7 | 217.1 | 199.8 | 198.4 | |
| Zhou et al. 2023 | 225–234 | 233.2 | 230.3 | 219.6 | 206.7 | 200.6 (24 kpc) | |
| Ou et al. 2024 | | 232.8 | 229.4 | 218.3 | 203.0 | 191.5 | 173.0 |
| Jiao et al. 2023 | | | 223.3 (10.5) | 218.8 (15.5) | 202.5 (20.5) | 187.1 (24.5) | 175.7 (26.5) |
| **Our law** | 201.6 | **211.2** | 211.1 | 204.4 | 197.5 | 192.2 | 190.3 |
| Our law, BR13 disk | 216.2 | 217.4 | 214.0 | 203.2 | 195.9 | 191.2 | 189.4 |
| MOND | 215.4 | 222.8 | 221.3 | 213.0 | 206.4 | 202.2 | 200.7 |
| Newton | 186.2 | 177.4 | 168.1 | 142.8 | 124.3 | 111.3 | 106.5 |
| Dark matter (McMillan's halo) | 226.5 | 233.7 | 232.6 | 227.2 | 223.4 | 220.8 | 219.8 |

Mean offset (model − data), 5–10 kpc / 15–27 kpc, and rms:

| | Eilers | Zhou | Ou | Jiao |
|---|---|---|---|---|
| Our law | −19.5 / −5.7 (14.6) | −23.8 / −12.9 (18.3) | −21.3 / −8.5 (16.2) | −10.4 / −3.1 (11.4) |
| MOND | −7.4 / +3.3 (8.7) | −11.6 / −4.0 (8.6) | −9.5 / +0.5 (9.4) | +0.2 / +6.2 (10.1) |
| Dark matter (McMillan) | +3.6 / +19.9 (16.4) | −0.6 / +12.4 (8.8) | +1.4 / +17.2 (14.8) | +11.3 / +23.3 (21.7) |

* **Beyond 15 kpc our law follows the Gaia curves** (3–9 km/s low, 13 against Zhou), and it
  falls with radius as the newer Gaia DR3 curves do. McMillan's dark-matter model stays near
  220 km/s and is 17–23 km/s high there.
* **Near the Sun our law is 8% low**: 211 against 229–234 km/s.
* **That shortfall is not special to the Milky Way.** At the Sun's Newtonian pull
  (1.26 × 10⁻¹⁰ m/s²) the median SPARC point sits +0.033 dex above our law and −0.006 dex from
  MOND (`transition_check_v7.py`). The Milky Way shows the same 0.04 dex gap between the two
  laws. Our law, as fitted, releases the companion slightly too late in the transition.
* **What would have to change.** Either:
  * **a later switch-off**: g_d × 1.5 raises the Sun's speed to 217.6 km/s (× 2: 221; × 3:
    225). SPARC allows it: × 1.5 gives a typical miss of 15.80 km/s (now 15.85) and moves the
    median residual in the transition from +0.033 to +0.004. × 2 gives 15.93. Or
  * **25% more stellar disk** (55–60 M☉/pc² at the Sun: rms 8.3–8.5 km/s against Eilers). That
    is well above the star counts (33.4 ± 3, McKee; 38 ± 4, Bovy & Rix) and it would push the
    vertical pull too high (§17.3). So the fix belongs in the law, not the disk.

### 17.3 The vertical pull above the Sun

The pull towards the disk 1.1 kpc above the Sun, as a surface density K_z/2πG (M☉/pc²).
Measured: 68 ± 4 (Bovy & Rix 2013), 74 ± 6 (Holmberg & Flynn 2004).

| Visible matter | Ours | MOND | Newton | Dark matter |
|---|---:|---:|---:|---:|
| Local census (McKee 2015) | **74.1** | 83.5 | 49.0 | |
| McMillan 2017 | 85.0 | 96.3 | 60.3 | 74.0 |
| Bovy & Rix 2013 disk | 79.6 | 90.2 | 57.2 | |

* **With the directly counted local matter, our law gives 74**: on Holmberg & Flynn, 1.5σ above
  Bovy & Rix. MOND gives 84.
* The same matter gives 195 km/s at the Sun. Together, §§17.2–17.3 say the extra pull near the
  Sun must lean towards the Galactic Centre more than the Newtonian direction does. Our
  direction rule already tilts the pull towards the hot matter's flow (the bulge and stellar
  halo). Near the Sun that flow is only 9% of the Newtonian pull (g_hot = 1.2 × 10⁻¹¹ against
  g_N = 1.26 × 10⁻¹⁰ m/s²), so the tilt is small. This is the "direction-dependent gravity" door
  that round 7 leaves open.

### 17.4 Mass inside a radius, escape speed, the inner Galaxy

Mass inside r (10¹¹ M☉; our law and MOND, McMillan's matter):

| r (kpc) | Measured | Ours | MOND | Dark matter |
|---|---|---:|---:|---:|
| 20 | 1.91 +0.18/−0.17 (Posti & Helmi 2019) | 1.68 | 1.83 | 2.24 |
| 21.1 | 2.1 +0.4/−0.3 (Watkins et al. 2019) | 1.75 | 1.92 | 2.36 |
| 50 | 4.5 ± 0.4 (Correa Magnus & Vasiliev 2022); 5.4 +1.1/−0.8 (Vasiliev 2019) | 3.58 | 4.16 | 5.11 |
| 100 | 6.07 ± 0.29 ± 1.21 (Deason et al. 2021); 7.3 +0.9/−0.8 (CMV22); 8.5 +3.3/−2.0 (V19) | **6.57** | 7.99 | 8.54 |
| 200 | 11.0 +2.7/−2.2 (CMV22) | **12.5** | 15.6 | 12.8 |

* **At 100–200 kpc our law agrees; MOND runs high at 200 kpc (+1.8σ).**
* At 20–50 kpc our law is 10–20% low (1.2–2.3σ). That is the same transition shortfall as §17.2.

**Escape speed at the Sun (km/s).**
* **Measured:**
  * 533 +54/−41 (Piffl et al. 2014, to 3 R₃₄₀);
  * 528 +24/−25 (Deason et al. 2019, to 2 r₂₀₀);
  * 497 ± 8 (Koppelman & Helmi 2021; 552 after their bias correction);
  * 580 ± 63 (Monari et al. 2018);
  * 512 +94/−37 (Prudil et al. 2022);
  * 486 +10/−5 (Roche et al. 2024);
  * 445 +25/−8 (Necib & Lin 2022).
* **Ours:** 507 to 400 kpc, 522 to 540 kpc, 594 to the companion's reach.
* **MOND:** 546 and 564 (it has no finite escape speed without an external field).
* **Dark matter:** 560 and 570.
* **Newton, visible matter only:** 266.

**The inner Galaxy.** Star counts and bulge microlensing say visible matter supplies
0.88 ± 0.07 of the inner circular speed (Wegg, Gerhard & Portail 2016). The inner 2–3 kpc lie
above g_d, so our law is nearly Newtonian there: 0.97–0.99 (1.3–1.6σ). MOND gives 0.91–0.93.

### 17.5 The dwarf galaxies around the Milky Way

Method (`mw_dwarfs_v7.py`):
* each dwarf is a Plummer sphere with the measured half-light radius (circularised);
* stars weigh 2 × the Sun per unit V-band light, varied 1–3;
* the Galaxy's Newtonian pull, heat and hot flow at its distance come from the §17.1 model;
* the pull averaged over directions is exact by Gauss's law, and the speed spread follows from
  the isotropic Jeans equation, weighted by light.

MOND is computed identically. The isolated MOND limit reproduces σ⁴ = (4/81) G M a₀ to 0.1%.

| Dwarf | Measured σ (km/s) | Ours | Ours, isolated | MOND | Newton |
|---|---|---:|---:|---:|---:|
| Fornax | 11.7 ± 0.9 | **12.4** | 12.6 | 13.6 | 5.3 |
| Leo I | 9.2 ± 1.4 | **9.3** | 9.4 | 10.2 | 4.6 |
| Leo II | 6.6 ± 0.7 | 5.3 | 5.4 | 5.9 | 2.0 |
| Sculptor | 9.2 ± 1.4 | 6.7 | 7.2 | 7.5 | 2.9 |
| Carina | 6.6 ± 1.2 | 3.5 | 4.4 | 4.0 | 1.3 |
| Sextans | 7.9 ± 1.3 | 2.1 | 4.5 | 2.4 | 0.8 |
| Draco | 9.1 ± 1.2 | 2.8 | 4.1 | 3.2 | 1.2 |
| Ursa Minor | 9.5 ± 1.2 | 3.5 | 4.2 | 3.9 | 1.4 |
| Crater II | 2.7 ± 0.3 | 1.0 | 3.4 | 1.2 | 0.3 |
| Antlia 2 | 5.7 ± 1.1 | 1.1 | 4.2 | 1.3 | 0.4 |

* **Right for the four bright or distant dwarfs** (Fornax, Leo I within 1σ; Leo II, Sculptor
  within 2σ).
* **1.5–5 times too slow for the six faint ones near the Galaxy or very diffuse**, where the
  Galaxy's pull is as big as or bigger than their own. MOND, computed the same way, is only
  10–20% higher.
* Published MOND estimates for these dwarfs (e.g. 2.1 km/s for Crater II, McGaugh 2016) use the
  estimator σ² = a₀GM/(3 g_e r₁/₂), which gives about 1.7 times our Jeans value.
* **What would have to change: most of it is the Galaxy's hold.** Weakening how strongly the
  Galaxy's pull enters each dwarf's law, by a factor c:

  | Galaxy's hold | 1 | 0.3 | 0.1 | 0.03 | 0 |
  |---|---:|---:|---:|---:|---:|
  | χ² for 10 dwarfs | 135 | 105 | 82 | 65 | 60 |
  | Crater II / Antlia 2 (km/s) | 1.0 / 1.1 | 1.3 / 1.5 | 1.7 / 1.9 | 2.3 / 2.6 | 3.4 / 4.2 |

  With no hold, Crater II (3.4 against 2.7 ± 0.3) and Antlia 2 (4.2 against 5.7 ± 1.1) come
  close. Draco and Ursa Minor would still need about twice their counted stars; those two are
  also MOND's classic hard cases and are thought to be tidally disturbed.

### 17.6 The precision tests: planets, S2, pulsars, light bending, and Cassini

* **Planets.** The release factor is 10^(−7.6×10⁷) at Mercury and 10^(−12,590) at Neptune. Even
  at Sedna's far point (937 AU) the extra pull is 10⁻¹⁴ of Newton's.
* **S2 at the Galactic Centre.** The pull runs from 1.8 m/s² at closest approach to
  6.7 × 10⁻³ m/s² at the far end of the orbit, and the release factor is 10^(−10⁷) or smaller. So the precession is exactly GR's; GRAVITY measure
  1.10 ± 0.19 (2020) and 1.135 ± 0.110 (2024).
* **Binary pulsars.** The orbits are pure GR. Our only new effect, the companion's mass loss
  (2.3 × 10⁻¹⁵ per year), changes Ṗ_b by 0.4% (Hulse–Taylor) and 1.6% (Double Pulsar) of the
  measurement error.
* **Light bending at the Sun.** GR's 1.7516″ (γ = 1). Measured: γ − 1 = (2.1 ± 2.3) × 10⁻⁵
  (Cassini), γ = 0.99984 ± 0.00015 (VLBI).
* **Cassini's test of the Galaxy's field inside the Solar System (Q2): fails as the law
  stands.** In laws like ours and MOND, the Galaxy's pull reshapes the Sun's own extra pull
  3,000–20,000 AU out, and that reshaping reaches the planets as a quadrupole. Hees et al. 2014
  measured Q2 = (3 ± 3) × 10⁻²⁷ s⁻² from nine years of Cassini ranging.
  * Computed exactly from the Sun's "as if" density: **ours 2.4–3.1 × 10⁻²⁶**. That is 2.4
    with the Galaxy's heat, 2.8 for a 1.9 × 10⁻¹⁰ m/s² Galactic pull, 3.1 for 2.09 × 10⁻¹⁰.
    MOND gives 2.8–3.1 × 10⁻²⁶, inside the range Hees et al. quote (2.1 × 10⁻²⁷ to 4.1 × 10⁻²⁶).
    Both are excluded at about 8σ.
  * Half of it comes from beyond 4,300–4,600 AU. A sharper switch-off makes it larger, not
    smaller (n = 2: 3.2 × 10⁻²⁶).
  * A hot flow from the Galaxy steering the Sun's extra pull reduces it: to 9 × 10⁻²⁷ if that
    flow equalled the Galaxy's Newtonian pull, 4 × 10⁻²⁷ at three times. The Milky Way's
    actual hot flow is 0.09 of it.
  * **What would have to change: the release takes time.** If the companion is released only
    gradually as it travels, F → F·(1 − e^(−r/L)):

    | Release length L | 0 | 10⁴ AU | 3 × 10⁴ AU (0.15 pc, 720 yr at u) | 10⁵ AU (0.48 pc) | 3 × 10⁵ AU |
    |---|---:|---:|---:|---:|---:|
    | Q2 (10⁻²⁷ s⁻²) | 26–31 | 11–12 | **4.2–4.6** | **1.4–1.5** | 0.5 |
    | Wide-binary boost at 7,000 / 20,000 AU | 1.19–1.26 / 1.19–1.27 | 1.09–1.13 / 1.17–1.24 | 1.04–1.05 / 1.09–1.13 | 1.01–1.02 / 1.04–1.05 | 1.00 / 1.01–1.02 |

    * L ≳ 0.15 pc passes Cassini.
    * It leaves every galaxy, dwarf and cluster untouched, since all are more than a thousand
      times larger.
    * It makes wide binaries nearly Newtonian (1–5% extra pull), in line with Banik et al.
      2024.
    * The upper limit is set by the smallest dwarfs (tens of pc): L ≲ 10 pc.
    * Physically, this is "gravity that builds up over time": escape from attachment has a
      finite rate. It is a candidate, not yet adopted: the fits must be re-run with it.

### 17.7 Lensing: the standard examples

* **Solar light bending, Cassini γ:** GR exactly (§17.6).
* **Microlensing towards the bulge.** The Einstein radius of a 0.5 M☉ lens halfway there is
  2.9 AU, where the pull is 1.6 × 10⁶ g_d. So each event is standard, and the optical depths
  (OGLE-IV: τ₀ = 1.36 ± 0.04 × 10⁻⁶; MOA-II: 1.75 ± 0.04 × 10⁻⁶) count stars, as in Newton.
  Our inner Galaxy is nearly Newtonian (§17.4).
* **The Einstein Cross (Q2237+0305).** The ring lies 0.65 kpc from the lens centre, where the
  pull is 22.7 g_d. The companion's share there is 1.5 × 10⁻¹¹, so our law says the lensing
  mass inside the ring is the stars. Measured: dark matter ≤ 20% (van de Ven et al. 2010);
  7%, < 15% at 3σ (Trott et al. 2010).
* **SLACS (six lenses):** stars 1.05–1.35 × Salpeter; light and matter agree to
  −0.017 ± 0.024 dex (rounds 3 and 5).
* **The KiDS-1000 lensing RAR (Brouwer et al. 2021).** Around about 259,000 isolated galaxies.
  Each law goes through the paper's own steps (`lensing_census_v7.py`): pull → projected
  excess surface density → g_obs = 4GΔΣ at the radius where g_bar has the bin's value.
  * Spirals are cold (k = 0.1 for small bulges).
  * Ellipticals are hot (k = 3σ²/u², σ = 160 km/s).
  * Median log(observed/predicted) in the seven bins where the isolation is reliable
    (g_bar ≥ 10⁻¹³). χ² in brackets includes the paper's 0.1 dex conversion uncertainty:

  | Lenses | Ours | MOND (M16 fit) | Newton |
  |---|---|---|---|
  | All isolated | **+0.024 (5.4/7)** | +0.048 (6.2/7) | +1.08 |
  | Blue (spirals) | **−0.005 (6.6/7)** | −0.073 (6.5/7) | +0.85 |
  | Red (ellipticals) | **+0.021 (1.2/7)** | +0.113 (7.4/7) | +1.16 |
  | Sérsic n < 2 (disks) | +0.032 (3.0/7) | −0.058 (2.8/7) | |
  | Sérsic n > 2 (bulges) | **−0.023 (2.7/7)** | +0.066 (8.0/7) | |
  | GAMA (spectroscopic) | −0.032 (11.4/7) | −0.007 (12.6/7) | |

  * **Our law fits the absolute lensing of both spirals and ellipticals.** MOND's single curve
    runs 0.07 dex above the spirals and 0.11 below the ellipticals.
  * The early/late gap is 0.18 dex (ours), against 0.15 observed in the same bins (the paper's
    reliable-region values: 0.19 by colour, 0.14 by Sérsic index).
  * The earlier worry that our deep-regime strength, a ≈ 0.55 a₀, would under-predict galaxy
    lensing is not borne out: spirals sit where the cold law puts them.
* **Circular velocities from lensing to 2.5 Mpc (Mistele et al. 2024).**
  * Flat to about 1 Mpc, as in our law. Ours stays flat to the companion's reach (2.0–2.6 Mpc)
    and then falls: a prediction.
  * For ellipticals our speeds match: 184 / 221 / 270 km/s against 197 ± 12 / 211 ± 8 /
    266 ± 5.
  * For the lightest spirals we are 20–30% below their values: 106 and 149 against 137 ± 8 and
    180 ± 9. Mistele et al. re-derive masses and velocities from the same data with their own
    method, and in Brouwer et al.'s framing the spirals agree with us. So this is a check of
    the conversions as much as of the law.
* **Colliding clusters beyond the Bullet.** The published values are recorded in
  `data/lensing_literature_v7.json`, for modelling next:
  * Abell 520's disputed "dark core" (P3: 2.3–4.4 × 10¹³ M☉ in 150 kpc, depending on the team);
  * MACS J0025.4−1222 (two 2.5 × 10¹⁴ M☉ clumps, lensing offset from the gas);
  * El Gordo (M200 ≈ 2–3 × 10¹⁵);
  * Abell 1689 (Einstein radius 47″; where MOND fails badly).

  Round 4's memory rule predicts lensing around stopped gas inside a sphere growing at
  200 kpc per Gyr, which bears directly on Abell 520's P3.

### 17.8 The coverage table

| Test | Measured | Ours | MOND | Dark matter | Verdict for ours |
|---|---|---|---|---|---|
| SPARC rotation curves (149) | | 15.9 km/s | 16.1 | 7.5 (fitted) | pass |
| X-COP cluster masses (12) | | 25% | ×2.9 | 11% (fitted) | pass |
| MW rotation 15–27 kpc (4 Gaia curves) | 173–217 km/s | −3 to −13 km/s | −4 to +6 | +12 to +23 | **pass** |
| MW rotation at the Sun | 229–234 km/s | 211 | 223 | 234 | **8% low; fix: later switch-off (§17.2)** |
| Vertical pull at 1.1 kpc | 68 ± 4; 74 ± 6 | 74 (counted matter) | 84 | 74 | pass |
| Mass inside 100 / 200 kpc | 6.1–7.3 / 11.0 × 10¹¹ | 6.6 / 12.5 | 8.0 / 15.6 | 8.5 / 12.8 | **pass** |
| Mass inside 20–50 kpc | 1.9–4.5 × 10¹¹ | 10–20% low | 4–8% low | high | 1–2σ low |
| Escape speed at the Sun | 445–580 km/s | 507–522 | 546–564 | 560–570 | pass |
| Inner Galaxy: visible share of speed | 0.88 ± 0.07 | 0.97–0.99 | 0.91–0.93 | 0.88 (fitted) | pass at 1.5σ |
| Dwarfs: Fornax, Leo I, Leo II, Sculptor | 6.6–11.7 km/s | within 1–2σ | within 0–2σ | fitted | pass |
| Dwarfs: 6 faint or diffuse | 2.7–9.5 km/s | 1.5–5× low | 1.3–4× low | fitted | **fail; fix: weaker Galactic hold** |
| Planets, Sedna, S2, pulsars, light bending | GR | GR exactly | small effects | GR | **pass** |
| Cassini Q2 | (3 ± 3) × 10⁻²⁷ s⁻² | 2.4–3.1 × 10⁻²⁶ | 2.8–3.1 × 10⁻²⁶ | ≈ 0 | **fail; fix: gradual release, L ≳ 0.15 pc** |
| Wide binaries | disputed | 19% (9% at 20,000 AU with gradual release; §18.1) | 43% | 0 | open |
| Microlensing, Einstein Cross | stars | stars | stars | stars | pass |
| SLACS strong lenses | | light = matter | | | pass |
| KiDS lensing RAR, spirals and ellipticals | | −0.005, +0.021 dex | −0.073, +0.113 | tuned | **pass, better than MOND** |
| Bullet: lensing on galaxies; main half | | yes | no | yes | pass |
| Bullet: smaller half's lensing mass | 2.0–2.3 × 10¹⁴ | about half | no | fitted | open |
| 72 collisions: lensing with galaxies | | yes | no | yes | pass |
| MACS J0025.4−1222 (§18.2) | 2.5 / 2.6 × 10¹⁴; 835 km/s | 2.0 / 1.6; 770 | no | fitted | **pass** |
| Abell 520, six clumps and the dark core (§18.2) | 2.1–5.6 × 10¹³ | 5 of 6; P3 3.06 | no | a puzzle | **pass** |
| El Gordo (§18.2) | 15.8 × 10¹⁴ inside 1 Mpc | 11.3; 16.9 with twice the stars | no | fitted | ~ stars to measure |
| Abell 1689 | recorded | not yet run | | | next |

### 17.9 What round 7 changes, and what comes next

* **No change to the law.** All numbers above use the round-3 law and constants.
* **Two candidate amendments, each required by a specific measurement and each testable:**
  1. **Gradual release** (L ≳ 0.15 pc), required by Cassini. It turns the wide-binary
     prediction from 19% into 1–5%.
  2. **A later switch-off** (g_d × 1.5–2), suggested by the Sun's rotation speed and allowed by
     SPARC.

  Both must be refit together with SPARC, X-COP and the Bullet before either is adopted.
* **One open door:** how strongly the Galaxy holds a dwarf's companion. Partial loss of step
  between a dwarf's companion and the Galaxy's, which move relative to each other at 100–300
  km/s (comparable to u), is the natural candidate; it would also reduce Cassini's Q2.
* **Next:**
  * refit with the gradual release;
  * model the four recorded collisions and Abell 1689;
  * the Galaxy's hot flow near the Sun, which sets the direction of the extra pull there
    (§17.3).

## 18. Round 8, 23 September 2026: a regression suite, and three more colliding clusters

The request: build a full regression suite, so that everything can be re-tested whenever the law
changes, then add the cluster collisions. No change to the law in this round.

### 18.1 The regression suite

`regression/` (its README has the details):
* `law_config.py`: the law as one dictionary. Candidates are small JSON files in
  `candidates/` with the round-7 amendments as keys: `gd_scale`, `release_length_au`
  (gradual release), `external_hold`, and `refit` (a on SPARC, u on X-COP).
* `checks.py`: one grading rule for every test:
  * pass: within 2 standard errors, or inside the required range;
  * close: within 3;
  * fail: further out;
  * info: tracked, ungraded.

  Each check also carries a score, its distance in standard errors.
* `t_*.py`: the tests, grouped as machinery, galaxies, clusters, lensing, milky_way, dwarfs,
  precision and collisions. Each calls the code that produced the published numbers.
* `run_suite.py`: runs a tier (quick, about 1 minute; full, about 20 minutes) and compares with
  `baseline.json`. Every check is marked regressed / improved / worse / better / changed /
  known. The exit code is 1 on any regression.

**Validation.** On the round-3 law the suite reproduces every published number to the last
digit:
* SPARC 15.85 / 12.40 / 19.21 km/s (all / test / validation);
* X-COP 0.227, held-out 0.244;
* KiDS +0.024 / −0.005 / +0.021 dex;
* SLACS −0.017 dex;
* Milky Way 211.2 km/s, Σ(1.1 kpc) 74.1;
* dwarfs χ² 135;
* Cassini Q2 3.09 × 10⁻²⁶ s⁻²;
* Bullet κ 0.675 / 0.14, M(<250 kpc) 2.405 / 0.941 × 10¹⁴, outer stars 6.23 × 10¹²;
* collision stack β = 0.027.

**The baseline (round-3 law, full tier):** 89 checks: **61 pass, 8 close, 7 fail**, 13 tracked. Close: lensing.mistele_ltg, mw.v_sun, mw.mass_50, dwarfs.carina, abell520.m150_p6, abell520.sigma_peaks, elgordo.m500, elgordo.sigma_nw. Fail: dwarfs.draco, dwarfs.ursa_minor, dwarfs.sextans, dwarfs.crater_ii, dwarfs.antlia_2, precision.cassini_q2, bullet.m250_sub.

**Candidates (quick tier; full table in `regression/README.md`):**

| Candidate | Fixes | Breaks |
|---|---|---|
| gradual release, L = 30,000 AU | Cassini Q2 (3.1 × 10⁻²⁶ → 4.6 × 10⁻²⁷ s⁻²) | nothing |
| external hold 10% | Cassini (3.8 × 10⁻²⁷); Carina; dwarfs χ² 135 → 82 | wide binaries: 2.5× Newton at 20,000 AU |
| g_d × 1.25 | – (the Sun 211 → 215 km/s) | Σ(1.1 kpc) 74.1 → 76.4 |
| g_d × 1.5 | the Sun (217.6 km/s) | SLACS gap −0.057 dex; Σ(1.1 kpc) 78.1 |
| g_d × 1.5, a and u refitted | the Sun (216.8) | the same, and Mistele's spirals |
| all three, refitted | the Sun, Cassini, Carina | the same three, and wide binaries (1.76×, close) |

Findings:
* Gradual release is free.
* The external hold must be speed-dependent: weaker for a dwarf moving past the Galaxy's
  companion at 100–300 km/s, not for binary stars that move with it.
* A later switch-off trades the Sun's speed against Σ(1.1 kpc) and the SLACS gap. The disk's
  shape is the better lever: Bovy & Rix's shorter disk gives 217 km/s at the Sun (§17.2).

**A correction to round 7.** "1–5%" for wide binaries with gradual release holds at 7,000 AU
(3.9–5.3% for L = 30,000 AU) or for L ≥ 100,000 AU. At 20,000 AU with L = 30,000 AU it is
9–13% (`run-strong-field-v7/strong_field_v7.json`; the suite reproduces 1.093).

### 18.2 Three more colliding clusters

`code/collisions_v8.py` → `run-collisions-v8/`. It uses the Bullet machinery unchanged
(`bullet_v4.kappa_map_v4`: memory, ghosts, fresh sphere u t) and published inputs only. Every
input and its source is in the output JSON, with the assumptions not taken from a paper listed
as caveats. An independent read of the papers checked the inputs; its corrections are applied:
* Kim's 97 kpc cool-core offset;
* Clowe's unsmoothed light and masses;
* the 1.5 Mpc aperture for MACS J0025's galaxy speeds;
* timings.

**MACS J0025.4−1222** (z = 0.586; Bradač et al. 2008). Inputs:
* two NFW-shaped galaxy populations normalised to 2.7 and 1.9 × 10¹² inside 300 kpc;
* one β-model gas cloud through the published 3.6 × 10¹³ (500 kpc sphere) and 5.5 × 10¹³
  (projected), with r_c = 90 kpc;
* before the collision, each subcluster held half the gas about its galaxies;
* t = 0.5 Gyr (0.26–1.0 give the same numbers to 1%).

| | Measured | Ours |
|---|---|---|
| M(<300 kpc), SE | 2.5 (+1.0/−1.7) × 10¹⁴ | 1.96 (z = −0.3) |
| M(<300 kpc), NW | 2.6 (+0.5/−1.4) × 10¹⁴ | 1.60 (z = −0.7) |
| M(<500 kpc), about the gas | 6.2 (+1.2/−4.0) × 10¹⁴ | 4.11 |
| lensing peaks | on the galaxies (> 4σ from the gas) | 15 / 43 kpc from the galaxies, 356 / 127 from the gas |
| galaxy speed spread (1.5 Mpc) | 835 ± 59 km/s | 770 km/s |

Stars × 2 overshoots the speeds (1,017 km/s), so the published stars are right for our law here.

**Abell 520** (z = 0.201). Inputs:
* gas: one β-model about P3 (r_c = 356 kpc) fitted to Clowe et al.'s six gas columns;
* light: Clowe et al.'s unsmoothed F814W light at M/L = 2 (the Bullet's convention), NFW-shaped
  (scale 100 kpc);
* before the collision, each clump held gas in proportion to its light;
* t = 0.3 Gyr (Girardi et al.; 0.5 and 1.0 change the masses by ≤ 2%).

| M(<150 kpc), 10¹³ | P1 | P2 | P3 | P4 | P5 | P6 |
|---|---|---|---|---|---|---|
| Jee et al. 2014 | 2.10 ± 0.43 | 4.05 ± 0.28 | 3.35 ± 0.34 | 4.23 ± 0.28 | 2.93 ± 0.39 | – |
| Clowe et al. 2012 (ζ_c) | 2.81 ± 0.67 | 4.16 ± 0.67 | 2.84 ± 0.64 | 5.59 ± 0.68 | 3.17 ± 0.66 | 3.68 ± 0.68 |
| ours | 2.56 | 4.84 | 3.06 | 3.30 | 2.83 | 2.10 |
| ours, visible matter only | 0.38 | 0.61 | 0.76 | 0.49 | 0.60 | 0.54 |

The "dark core" P3 gets 3.06 × 10¹³ from 0.76 × 10¹³ of visible matter (mostly gas):
* the companion's boost of ordinary matter;
* plus the heat of the surrounding hot galaxies (the scalar sum S does not cancel);
* a lensing peak sits 25 kpc from P3.

The cluster inside 710 kpc: 4.33 × 10¹⁴ against 5.0 ± 0.55 (Mahdavi et al. 2007). Galaxy
speed spreads near P1/P2/P4/P5 are 432/587/526/438 km/s, against 811/749/579/668 (Girardi et al.
2008, 6–9 galaxies each; rms z 2.9). Twice the stars would match those speeds (550–760 km/s) but
overshoot the lensing (P3 4.4 × 10¹³), so Abell 520's lensing prefers its published light.

**El Gordo** (z = 0.870). Inputs:
* stars 7.5 and 5.6 × 10¹² (Menanteau et al. 2012, Chabrier SED fits);
* gas 2.1 × 10¹⁴ inside R500 in a β-model (r_c = 250 kpc) plus a 10¹³ cool core 97 kpc beyond
  the SE lensing peak;
* before the collision, 60/40 gas shares about the galaxies;
* t = 0.46 Gyr (outgoing; the returning 0.91 Gyr gives the same masses).

| | Measured | Ours | Ours, stars × 2 |
|---|---|---|---|
| M_2D(<0.5 / 1 / 1.5 Mpc) about the centre of mass | 6.70 / 15.78 / 22.57 (Kim et al. 2021 NFW sum, consistent with their aperture masses) | 4.61 / 11.34 / 17.04 | 7.05 / 16.86 / 24.80 |
| galaxy speed spread NW / SE (1 Mpc) | 1,290 ± 134 / 1,089 ± 200 | 944 / 839 | 1,187 / 1,047 |
| SE lensing peak from the cool core | about 97 kpc (12″), the peak nearer the centre | 159 kpc (peak 62 kpc from the galaxies, toward the centre) | 148 kpc |

* With the published stars, the lensing is 24–31% low and the speeds 23–27% low.
* The peak's distance from the cool core depends on where the bulk of the gas sits, which no
  paper maps (we assumed it), so the suite tracks it without grading. The order along the
  merger axis matches the observations: cool core, galaxies, lensing peak, centre.
* Twice the stars recovers both, within 10% and one standard error. That is inside the authors'
  stated factor-of-two uncertainty.
* At the X-COP median star/gas ratio (0.062 inside R500, against 0.042 for El Gordo as
  published) the stars would be 1.5 times the published value, halfway there.
* Twice the stars overshoots MACS J0025's galaxy speeds (1,017 against 835 ± 59 km/s) and
  Abell 520's lensing (P3 4.4 × 10¹³). So this is a statement about El Gordo's stars, not a
  rule for all clusters.

### 18.3 What round 8 changes, and what comes next

* **No change to the law.** The suite makes every future change testable against 89 checks in
  one command. RULES.md §12 makes that the standing practice.
* **Collisions:**
  * MACS J0025 passes everything with its published inputs.
  * Abell 520's "dark core" comes out of gas and heat alone, between the two teams' masses.
  * El Gordo passes if its stars are twice the colour-based estimate, a testable statement about
    its stars.
* **Next:**
  1. a speed-dependent external hold (dwarfs and Cassini without touching binaries);
  2. the Milky Way disk's scale length as the lever for the Sun's speed;
  3. El Gordo's and Abell 520's stellar masses from independent data (near-infrared,
     spectroscopy), and A520 with two main pre-collision subclusters instead of one group per
     clump;
  4. Abell 1689 and the elliptical gas shape of MACS J0025 (Riseley et al. 2017);
  5. the Bullet subcluster's crossing heat (§16.3) as a suite candidate.

## 19. Round 9, 23 September 2026: the Cassini fix adopted, what it opens, and the plan

The request: the Cassini fix was called free but not adopted. Say what has to change, check it
with the regression suite, and find out whether it can make everything better. Then give the
open items a proper plan.

### 19.1 The change: the companion's release takes time (adopted)

* **Before (round 3).** A companion leaving an emitter is released as soon as the pull drops
  below g_d: `extra = exp(−|g_N|/g_d) · √(a (|g_N| + S)) · direction`. Around the Sun that is
  beyond about 5,100 AU (where GM/r² = g_d).
* **Now.** The release also needs travel. Around each emitter, the released share builds up as
  `R(r) = 1 − exp(−r/L)`. It multiplies the extra pull: `extra → R · extra` (round 7, §17.5; the
  form tested there). L = 30,000 AU = 0.15 pc, about 720 years of travel at u. This is the
  smallest L that Cassini allows.
* **Where it acts.** Only within a fraction of a parsec of a star:
  * a wide binary's own pull;
  * the Sun's own companion, which the Galaxy's field distorts into Cassini's Q2.

  Anything larger is unchanged. At L = 1 pc, a dwarf galaxy's companion at its half-light radius
  is at most 0.54% unreleased (Ursa Minor; conservative scalar estimate). At 10 pc, at most 5.5%.
  Galaxies and clusters are thousands of times larger still.
* **Cost.** One new number. The law now has three constants fitted to galaxies and clusters
  (a, g_d, u) and a release length bounded by the Solar System: at least 0.15 pc (Cassini), at
  most a few parsecs (the dwarfs).

**The regression suite, full tier** (`regression/runs/`, new baseline `round9`):

| | Round 3 (baseline until now) | Round 9 (adopted) |
|---|---|---|
| pass / close / fail (76 graded) | 61 / 8 / 7 | **62 / 8 / 6** |
| Cassini Q2 | 3.09 × 10⁻²⁶ s⁻² (fail) | **4.59 × 10⁻²⁷ (pass**; measured (3 ± 3) × 10⁻²⁷) |
| Cassini Q2 with the Milky Way model's pull and heat | 2.39 × 10⁻²⁶ | 3.77 × 10⁻²⁷ |
| wide binaries, 20,000 AU | 1.192 | 1.093 (still between the two analyses) |
| every other check | | identical: 80 the same, 6 known fails, 0 regressed |

**One prediction changes.** Wide binaries: 4% more pull than Newton at 7,000 AU and 9% at
20,000 AU, instead of 19% (for L at its minimum; smaller for longer L).

The suite's comparison now also marks a number that moved while its grade and score did not
(**changed**; e.g. the wide binaries inside their passing range). It records the commit a run
came from, with `+changes` when the working tree differs from it. The pre-round-9 law remains
available as `--law round3`.

### 19.2 What the fix opens: the faint dwarf galaxies

`code/release_hold_v9.py` → `run-release-hold-v9/release_hold_v9.json`.

**Why the dwarfs are slow.** Inside each dwarf, the law's square root contains the total pull,
and that includes the Milky Way's:
* At the half-light radius the Galaxy's Newtonian pull is 1.9× the dwarf's own for Draco, 8.6×
  for Sextans, 54× for Crater II and 82× for Antlia 2.
* Under the square root, a large steady pull dilutes the dwarf's own contribution. This is the
  external field effect, familiar from MOND and stronger in ours.

Splitting it up (χ² over the ten dwarfs):
* the law as it stands: 135;
* without the Galaxy's heat: 133;
* without the Galaxy's pull: 80;
* without both: 60.

MOND with the same stars scores 119. So the Galaxy's pull is the cause.

**Why removing the hold was blocked, and why it no longer is.** The hold c scales the Galaxy's
pull and heat inside a separate system (a dwarf, the Sun, a wide binary). With c = 0 and no
release length, wide binaries would pull 2.97× Newton at 20,000 AU. The release length now
protects them:

| Release length | c = 1 | c = 0.3 | c = 0.1 | c = 0 |
|---|---|---|---|---|
| none | 1.19 | 1.72 | 2.48 | 2.97 |
| 0.15 pc = 30,000 AU (adopted) | **1.09** | 1.35 | 1.72 | 1.96 |
| 0.24 pc = 50,000 AU | 1.06 | 1.24 | 1.49 | 1.65 |
| 0.48 pc = 100,000 AU | 1.03 | 1.13 | 1.27 | **1.36** |
| 0.97 pc = 200,000 AU | 1.02 | 1.07 | 1.14 | **1.19** |
| 2.4 pc = 500,000 AU | 1.01 | 1.03 | 1.06 | 1.08 |

*The wide-binary boost at 20,000 AU; the two published analyses give 1.0 and 1.4–1.5. Cassini's
Q2 is ≤ 4.6 × 10⁻²⁷ in every row from 0.15 pc on, and 0 at c = 0.*

**The ten dwarfs as the hold weakens** (km/s; M/L_V = 2; grades as in the suite):

| Dwarf | Measured | c = 1 | 0.3 | 0.1 | 0.05 | 0 |
|---|---|---|---|---|---|---|
| Draco | 9.1 ± 1.2 | 2.83 | 3.56 | 3.93 | 4.03 | 4.14 |
| Ursa Minor | 9.5 ± 1.2 | 3.49 | 3.97 | 4.14 | 4.18 | 4.22 |
| Sculptor | 9.2 ± 1.4 | 6.72 ✓ | 7.08 ✓ | 7.19 ✓ | 7.22 ✓ | 7.24 ✓ |
| Sextans | 7.9 ± 1.3 | 2.14 | 2.79 | 3.57 | 3.98 | 4.49 ~ |
| Carina | 6.6 ± 1.2 | 3.53 ~ | 4.11 ~ | 4.32 ✓ | 4.37 ✓ | 4.43 ✓ |
| Fornax | 11.7 ± 0.9 | 12.36 ✓ | 12.52 ✓ | 12.56 ✓ | 12.57 ✓ | 12.58 ✓ |
| Leo II | 6.6 ± 0.7 | 5.29 ✓ | 5.35 ✓ | 5.37 ✓ | 5.37 ✓ | 5.37 ✓ |
| Leo I | 9.2 ± 1.4 | 9.35 ✓ | 9.36 ✓ | 9.37 ✓ | 9.37 ✓ | 9.37 ✓ |
| Crater II | 2.7 ± 0.3 | 1.01 | 1.34 | 1.73 | 2.05 ~ | 3.39 ~ |
| Antlia 2 | 5.71 ± 1.08 | 1.12 | 1.48 | 1.93 | 2.28 | 4.18 ✓ |
| **χ² (10 dwarfs)** | | 135 | 105 | 82 | 71 | **60** |

**The candidate `no_hold`** (c = 0 with L = 200,000 AU, about 1 pc; quick tier against the new
baseline): 41 pass, 5 close and 2 fail of the quick tier's 48 graded checks (the baseline's quick
subset is 39, 4, 5). The collisions read neither the hold nor the release, so in all it scores
**64 pass, 9 close, 3 fail**.
* Improved: Carina and Antlia 2 (to pass), Sextans and Crater II (to close).
* Draco and Ursa Minor are better but still fail.
* No check drops a grade. Fornax moves from 12.36 to 12.58 km/s against 11.7 ± 0.9, still a pass.
* Cassini's Q2 goes to 0 (a pass); wide binaries sit at 1.19.

**What is not solved.**
* **There is no derivation yet.** The law's heat rule already fixes how the companions of
  randomly moving galaxies combine: a scalar sum inside one square root. So "separate systems
  pull separately" is not available: it would change every cluster. The reason the Galaxy's
  pull weighs less inside a dwarf has to be something else, and must leave galaxies, clusters
  and the heat rule untouched.
* The data already constrain it. Crater II sits between: c = 0.05 gives 2.05 and c = 0 gives
  3.39, against 2.7 ± 0.3, so it wants a small hold. Antlia 2 wants none.
* **Draco and Ursa Minor** stay at 4.1–4.2 km/s against 9.1–9.5 whatever the hold. They are
  also the hardest cases for MOND. Both are among the closest of the ten (76–78 kpc), and a dwarf
  without dark matter is easily stirred by the Galaxy's tides.

### 19.3 The plan

Ordered by how directly each item can move the scoreboard. Each ends in a suite run. The owner's
ten longer-range proposals are evaluated in §19.5.
1. **Borrowed assumptions in the far-cluster inputs** (§19.4). Rebuild sizes, gas, stars and
   lensing for the Bullet, Abell 520, MACS J0025 and El Gordo with the project's static distance
   law and without Big-Bang age caps on stellar masses. Prefer aperture masses to NFW fits.
   Switch the suite's SLACS check to the project's distance convention.
2. **The dwarfs' hold** (5 fails, 1 close).
   * Find why the Galaxy's pull would weigh less inside a separate system, consistent with the
     heat rule (proposal 5).
   * Formulate it as a suite candidate. Grade it on the ten dwarfs, Cassini, the wide binaries
     and everything else.
   * Separately, test tides for Draco and Ursa Minor: their Jacobi radii under the law against
     their measured extents, and stretching or velocity gradients in Gaia and spectroscopic
     data.
3. **The Bullet Cluster's smaller half** (1 fail: 0.94 against 2.0–2.3 × 10¹⁴ inside 250 kpc).
   * Build §16.3's crossing heat, the shaking its galaxies got while crossing the main cluster,
     as a collisions candidate. Then build the full transport calculation (proposal 6).
   * Check whatever step 2 finds on it: the smaller half moves through the main cluster's
     companion at 4,500 km/s.
   * Redo it with the static distances.
4. **The Milky Way's matter** (2 close: the Sun's speed 211 against 229–234 km/s; the mass inside
   50 kpc, 3.6 against 4.5 ± 0.4 × 10¹¹).
   * Add the disk's scale length and the bar's mass as Milky Way model options in the suite.
     Bovy & Rix's (2013) shorter disk already gives 217 km/s, which would pass.
   * Run them against the pull above the disk, the 15–27 kpc curve and the masses together.
5. **Abell 1689**, not yet run: a massive relaxed cluster and strong lens, where MOND needs extra
   unseen mass. Model it with its published gas and stars (static distances), add it to the
   suite, and grade it against the lensing mass profile and the Einstein radius.
6. **Cluster stars** (4 close: El Gordo inside 500 kpc and its NW galaxy speeds; Abell 520's P6
   and its galaxy speeds).
   * Step 1 first.
   * Then infrared and spectroscopic stellar masses (proposal 9).
   * Abell 520 with two main subclusters before the crash instead of one group per clump.
7. **Lensing speeds of low-mass spirals** (1 close; rms z 3.0 against Mistele et al.). Our
   speeds agree with Brouwer et al.'s conversion of the same KiDS data. Grade the law against
   the published lensing profiles directly, with no conversion step and the static distances.
8. **Wide binaries** (data disputed). The prediction is now 4% at 7,000 AU and 9% at 20,000 AU,
   less for a longer release. With step 2's candidate it would be 19–36%. The binaries measure
   the release length, and they separate the two dwarf routes.
9. **Theory** (proposals 1, 3, 4, 8 and 10):
   * the field equation and a relativistic form with the release length;
   * what L is;
   * u and g_d from first principles.

   On L: one guess is that L is the companion's near zone, since a radio wave detaches from its
   antenna at about λ/2π. That gives λ ≈ 2πL ≈ 1 pc. Then postulate 3 needs cluster ions to be
   pitch-angle scattered along the magnetic field at least every ~1,000 years, far more often
   than Coulomb collisions (mean free path ~10 kpc). Plasma instabilities may do it; if not, L
   is an escape time unrelated to λ.
10. **Cosmology**: the microwave background and large-scale structure without expansion, still
   outside the law's tested scope.

### 19.4 Borrowed assumptions: an audit (prompted by the owner's round-9 instruction)

The owner asked that nothing inherited from dark matter, an expanding universe or a Big Bang
steer the work. Checking the inputs found four places where it had.

1. **Distances.** `collisions_v8.py` and `bullet_v3.py` convert angles and fluxes with flat ΛCDM
   (H0 = 70, Ωm = 0.3), as do the published masses they compare with. The project's own law is
   static and Euclidean, 1 + z = e^(αD) with α = 2.489 × 10⁻⁴ Mpc⁻¹, D_A = D and D_L = (1 + z) D.

   Static ÷ ΛCDM, at fixed observed angle and flux. Sizes scale as D_A, stars as D_L², X-ray gas
   as D_L D_A^1.5, and weak-lensing mass inside a fixed angle as D_l D_s/D_ls (typical source
   redshifts 1.0–1.3):

   | | z | kpc/″ ΛCDM | static | sizes | stars | gas | lensing |
   |---|---|---|---|---|---|---|---|
   | X-COP (typical) | 0.06 | 1.16 | 1.13 | 0.98 | 0.85 | 0.89 | – |
   | Abell 520 | 0.201 | 3.31 | 3.57 | 1.08 | 0.80 | 1.00 | 1.10 |
   | SLACS lens (typical) | 0.20 | 3.30 | 3.55 | 1.08 | 0.80 | 1.00 | 1.11 |
   | Bullet Cluster | 0.296 | 4.41 | 5.05 | 1.14 | 0.78 | 1.08 | 1.18 |
   | MACS J0025 | 0.586 | 6.61 | 8.98 | 1.36 | 0.73 | 1.36 | 1.39 |
   | El Gordo | 0.870 | 7.71 | 12.19 | 1.58 | 0.71 | 1.68 | 1.60 |

   At low z the difference is mostly the scale: cα = 74.6 km/s/Mpc against H0 = 70.
   * The SLACS lenses were computed both ways (`lenses_t35.py`). The light–matter gap is
     −0.012 ± 0.023 dex (static) against −0.017 ± 0.024 (ΛCDM), but the suite grades the ΛCDM
     version.
   * The far collisions need rebuilding. The direction of El Gordo's change is not obvious: its
     gas and lensing rise, its stars fall, and its scale grows.
2. **Stellar ages.** SED fits of distant galaxies usually cap the age at the ΛCDM age at the
   galaxy's redshift: 11.0 Gyr at z = 0.20, 10.1 at 0.30, 7.9 at 0.59 and 6.3 at 0.87 (13.5 today).
   Older populations have higher M/L. The cap bites hardest at El Gordo, the one cluster that
   "needs" twice its published stars. MACS J0025 does not. The colours of its galaxies will
   constrain their ages directly.
3. **Profile shapes.** El Gordo's comparison masses are Kim et al.'s NFW sums, consistent with
   their aperture masses. Aperture masses should be preferred wherever published.
4. **The reach.** u × 13 Gyr assumes emission since the oldest stars formed. With no Big Bang,
   matter may be older than its stars. No graded test reaches 2.6 Mpc, but prediction 13 (the
   outer lensing turnover) measures u times the emission time.

Collision timings taken from dark-matter simulations (El Gordo, MACS J0025) do not matter here:
the results change by ≤ 2% over 0.26–1.0 Gyr (§18.2).

### 19.5 The owner's ten proposals for completing the theory: evaluation

The owner proposed ten routes, each described as a hypothetical completed solution, not a result.
BLOG §9.1 evaluates each in plain language. In brief:

| # | Proposal | Status now | First concrete step | Depends on |
|---|---|---|---|---|
| 1 | Derive the full law from one microscopic matter–companion interaction; an independent simulator recovers it | nine postulates; Dicke toy and 3D field equation only. Postulate 4 (pull = amplitude) is the weakest link: a steady push from wave energy normally scales as A², giving 1/r² | a 1D wave-medium toy with bound and free states; test the interference (cross) term between a body's own companion and a passing one, which is linear in the passing amplitude | – |
| 2 | Locked, blind prediction recovered by independent teams | the suite re-checks known data; 18 predictions | tag the law and distance convention; commit a forecast file now (Gaia DR4 wide binaries: 1.04 at 7,000 AU, 1.09 at 20,000 AU); then lensing by σ-bins from KiDS-Legacy/Euclid inputs | 19.4 (distances) |
| 3 | Random motion gives a steady inward pull; collisions suppress it continuously | Dicke toy: scalar vs vector sums and the collisional suppression; no force measurement | add test bodies; σ → 2σ must give k → 4k and g_extra → 2 g_extra where heat dominates; vary the collision rate at fixed energy; check that catalogue splits change nothing | 1 |
| 4 | Light bending and matter motion from the same coupling | assumed (factor 2); SLACS −0.012 to −0.017 dex, KiDS | a light-like mode in the toy medium; the factor 2 must emerge | 1 |
| 5 | One mechanism for Cassini, dwarfs and binaries | L adopted; dwarfs need c ≈ 0 (reason unknown); Draco and UMi short even at c = 0 | a bound-to-free transition giving τ, and the overlap between companions in relative motion; the dwarf deficits must close, not shrink | 1 |
| 6 | Evolve the companion through mergers | memory imposed; Bullet subcluster 0.94 against 2.0–2.3 | a time-dependent transport solver started from observed gas and galaxies (static distances), checked against the steady memory rule; derive when lensing can peak on gas | 1, 3 |
| 7 | Measure u independently | u from X-COP only | the projected lensing of a growing fresh sphere (30 / 101 / 201 kpc after 0.15 / 0.5 / 1 Gyr); collisions timed by their shocks | 6 |
| 8 | Derive the constants (u² = κ/χ; a = 2ℓ/u; g_d = u²/2L_d) | fitted; clues a ≈ cα/10, g_d ≈ cα/3 with α the static redshift rate | after 1; also why L_d ≈ 2.8 kpc (barrier) and L = 0.15 pc (escape time) differ | 1 |
| 9 | Stellar masses from non-gravitational evidence | El Gordo ×2, SLACS 1.05–1.35 Salpeter, MACS ×1 (all with ΛCDM distances) | redo with static distances and no age cap (19.4); then NIR photometry and spectra | 19.4 |
| 10 | Derive and detect the energy cost | ℓ = 6.48 × 10⁻⁶ W/kg; Ṁ/M = −2.3 × 10⁻¹⁵ yr⁻¹ | energy and momentum balance in the toy. Already implied: the power cannot be thermal, since the Earth's 3.9 × 10¹⁹ W is 8 × 10⁵ times its internal heat flow and a 10⁻³ L☉ white dwarf's is 20 times its luminosity. The mass loss mimics Ġ/G = −2.3 × 10⁻¹⁵ yr⁻¹, about 30 times below lunar laser ranging's precision | 1 |

## 20. Round 10, 24 September 2026: first principles, part 1, and our own distances

The request: proceed to the next steps that solidify the theory and its first principles.
This round takes the first step of proposal 1 (§19.5), where "the companion pulls with its
amplitude" can come from. It tests what that step implies against the data, takes the far
collisions and the strong lenses out of the expanding universe's distances (§19.4), and locks
a first forecast (proposal 2).

### 20.1 The pull of a locked emitter: derived, and simulated

`code/first_principles_v10.py` → `run-first-principles-v10/first_principles_v10.json`.

**Model.** A scalar field ψ (the companion) obeys (1/u²)ψ_tt − ∇²ψ = Σ_i q_i(t) δ(x − x_i). The
only interaction is the local coupling q_i ψ(x_i), so a body feels F = q(t) ∇ψ_other. Each body
is a self-sustained emitter, q = q₀ sin(ωt + φ), whose phase locks to the wave passing it, with
offset Δ.

**Exact result** (retarded field of a point source Q; A = Q/4πR is the local amplitude of ψ):

  ⟨F_r⟩ = −(q₀ A / 2) [ k sin Δ + cos Δ / R ]

* **Δ = +π/2**, the emitter a quarter cycle ahead of the passing wave: a pull toward the source
  of exactly (q₀k/2) × the local amplitude, at every distance, near or far. This is postulate 4
  with a mechanism.
* **Δ = 0**, in step: a 1/R² pull, Newton-like. The companion gives no such term at quadrature.
* **Δ = −π/2**, a quarter cycle behind, which is what a passive absorbing response does: a push.

At Δ = +π/2 the emitter feeds the passing wave coherently at power P = F · v_phase: the pull is
the recoil of stimulated emission.

**The same from a 3D simulation of the wave equation alone.** The force law is not programmed
in; the emitter only reads the phase of the field it sits in.
* Grid: 192³, 16 cells per wavelength.
* Along the grid axis, pull ÷ (q₀k/2)A is 0.97–1.04 from kR = 1.6 to 19. Along the diagonal it
  is 0.94–0.98.
* A·R is constant to 2% out to R = 32 cells, and within 7% out to 64 cells, next to the absorbing
  layer.
* The quarter-behind case is pushed with the same size of force.

**Three conditions the microscopic model must meet** (derived here):
1. **An active medium.** Matter must feed the companion a quarter cycle ahead of the wave passing
   it, like a gain medium or a laser amplifier. Passive, absorbing matter would be pushed.
2. **Spectral locking, for universality.**
   * If all the emitting units in a body lock to one wave, they emit coherently. Power then grows
     as mass², so a galaxy's extra pull would scale as v² ∝ M, against v⁴ = G M a.
   * If they lock only weakly, the pull goes as the amplitude squared, which is Newton-like.
   * What works: each unit locks within its own narrow frequency band. Power then stays
     proportional to mass, and the pull stays proportional to the amplitude, with the same
     acceleration for every body.
3. **A slow phase.** The coherent power that pays for the pull is P = m g_extra v_phase per body.
   With v_phase = u, that would be 0.25 / 0.75 / 1.6 / 1.1 times the base feed ℓ at
   g_N = 10⁻¹² / 10⁻¹¹ / 10⁻¹⁰ / 3 × 10⁻¹⁰ m/s². It would amplify a galaxy's outgoing companion
   by order unity and spoil the tight v⁴ = G M a. So v_phase ≲ 0.05 u: the companion's crests
   must move much more slowly than its energy does.

A side result: the heat weight's velocity scale in the Dicke toy (§10) is Γ/k, the locking rate
over the wavenumber. Its travel speed is the group velocity. The law uses one u for both.
Measuring u from travel alone (proposal 7) therefore tests something new.

### 20.2 How many sources add up, and the clusters' verdict

For sources with independent phases (one frequency band), the pull on a locked emitter is, in
the Gaussian limit,

  F = (√π/2) (q₀k/2) Σ_j A_j² r̂_j / √(Σ_j A_j²)

Monte Carlo over 2,000 phase sets confirms it to 2–5% (one source: exact; opposite equal
sources: 0.0003 ± 0.0005 of one source's pull).

In the law's terms this is extra = √a (g_N + g_hot)/√(S_N + S), where S_N = G∫ρ/d² is the
unsigned sum of Newtonian pulls. It equals the round-3 law for a point mass and in the outskirts.
Inside extended systems it is weaker by the anisotropy factor √(|g_N|/S_N):

| Where | factor |
|---|---|
| uniform sphere at r/R = 0.25 / 0.5 / 0.75 / 0.99 | 0.29 / 0.43 / 0.57 / 0.79 |
| Milky Way plane at 4 / 8.2 / 12 / 20 / 30 kpc | 0.56 / 0.66 / 0.73 / 0.82 / 0.89 |
| cluster-like mix at 100 / 250 / 500 / 1,000 / 1,500 kpc (derived ÷ round 3) | 0.28 / 0.39 / 0.49 / 0.59 / 0.65 |

**Tested on data, it fails.**
* **X-COP**, with a and g_d held and u refitted: best u 75 km/s, rms 0.382. The mean residual
  runs from +0.45 at the innermost radius to −0.57 at the outermost. The round-3 rule, fitted the
  same way, gives u 197 km/s and rms 0.227, with residuals from +0.06 to −0.19.
* **The Milky Way** (local estimate): the Sun's speed would drop from 211 to 200 km/s.

So the companion inside galaxies and clusters does not add up like independent waves. What the
data select, round 3's rule, behaves like a single coherent flow:
* its strength is set by all the companion present, with no penalty for arriving from many sides;
* its direction is set by the net flow;
* for cold matter, the density follows the net current, the way a superfluid moves as one wave.

That is the target of the next toy: a coherent, condensate-like companion, tested the same way.
The analogy is ours, from the round-3 rule. Superfluid dark-matter theories are dark-matter
theories and are excluded.

### 20.3 The far collisions and the strong lenses in the project's own distances

`code/collisions_v10.py` → `run-collisions-v10/collisions_v10.json`. Every length of the round-8
models is rescaled by the size factor, gas by the X-ray factor and stars by the luminosity factor.
Measured lensing masses are multiplied by the lensing factor. Apertures are taken at the same
angles. The law's physics is unchanged. Effective source redshifts are 1.2 (MACS J0025), 1.0 (Abell 520) and 1.3 (El Gordo). These are
typical of their weak-lensing catalogues, and are to be replaced by the papers' own values. Only
the conversion of the measured lensing masses depends on them.

| | factors: size / stars / gas / lensing | Ours ÷ measured, lensing (round 8, LCDM) | Galaxy speeds, ours vs measured |
|---|---|---|---|
| MACS J0025, M(<300 kpc) SE / NW | 1.36 / 0.73 / 1.36 / 1.39 | 0.64 / 0.53 (0.78 / 0.61); z −0.5 / −0.9 | 707 vs 835 ± 59 km/s, z −2.2 (770) |
| Abell 520, six clumps inside 150 kpc | 1.08 / 0.80 / 1.00 / 1.10 | P1, P2 (now inside the range), P3, P5 pass; P4 z −2.1 (−1.4); P6 z −2.6 (−2.3) | – |
| Abell 520, inside 710 kpc | same | 0.77, z −2.1 (0.87, −1.2) | – |
| El Gordo inside 0.5 / 1 / 1.5 Mpc | 1.58 / 0.71 / 1.68 / 1.60 | 0.66 / 0.71 / 0.76 (0.69 / 0.72 / 0.75) | NW 917, SE 822 vs 1,290 / 1,089 (944 / 839) |
| El Gordo, twice the published stars | same | 0.92 / 0.94 / 0.98 (1.05 / 1.07 / 1.10) | NW 1,100, SE 974 (1,187 / 1,047) |

* **El Gordo:** the change of distance law is nearly a wash. Its gas and size rise as its stars
  fall. It still needs about twice the published stars, now 1.9 × 10¹³ M☉ in the static law.
* **MACS J0025 and Abell 520** do slightly worse, though mostly still passing. Their measured
  lensing rises by 10–39%; our lensing rises less, because their stars come out 20–27% lighter.
* **Not yet in the suite.** It still grades the round-8 (LCDM-geometry) collisions. Switching it
  waits on the papers' own source redshifts and on stellar masses without the age cap (§20.5).

**The six strong lenses** (the suite's check now uses the static convention):
* light = matter: −0.012 ± 0.023 dex (was −0.017 ± 0.024);
* the stars needed: 0.465 dex above Chabrier, that is 1.44–1.95 × Salpeter per lens (was
  1.05–1.35). Our distances give 1.48 times less starlight mass for the same light, and 1.11 times
  more lensing mass.

**What this means.** With our own distances the measured lensing masses rise by 10–60%, while the
stars inferred from the same light fall by 20–30%. In our law the stars carry the heat, so every
far test now leans on the stellar masses. The published ones still carry assumptions from the
expanding-universe timeline: the age cap of §19.4, where older stars weigh more for their light.
Redoing them is the next step.

### 20.4 A locked forecast: wide binaries for Gaia DR4 (proposal 2)

`code/forecast_wide_binaries_v10.py` → `forecasts/wide_binaries_gaia_dr4_v10.json`, committed
before the data. The law is the adopted round-9 form, and the Galaxy's pull at the Sun is as in
the suite. γ = pull ÷ Newton at 3D separation s, for solar-mass pairs:

| s (AU) | 3,000 | 5,000 | 7,000 | 10,000 | 20,000 | 30,000 | 50,000 |
|---|---|---|---|---|---|---|---|
| adopted law | 1.002 | 1.022 | 1.039 | 1.054 | 1.093 | 1.121 | 1.155 |
| candidate no_hold (not adopted) | 1.000 | 1.005 | 1.015 | 1.040 | 1.188 | 1.427 | 2.151 |

Pairs of 1.5 and 2 suns give the same values beyond 10,000 AU. The file states what would count
as support and as refutation, and carries a SHA-256 of the forecast
(94a7a470953cf63cb004e7ce39b0fb188720750206c7415520c013a09da22626).

### 20.5 Where this leaves the first principles

* **Derived now:**
  * the pull proportional to the companion's amplitude, from a local coupling, at every distance;
  * the direction along the net flow;
  * no energy exchange unless the pull is paid by coherent emission.
* **Required of the microscopic model, now explicit:** an active medium (a quarter cycle
  ahead), spectral locking (universality), and a slow phase (bookkeeping).
* **Learned from the data:** the companion adds up as one coherent flow, not as independent waves.
* **Next:**
  1. a coherent-flow (condensate-like) toy, graded on the same X-COP and Milky Way tests;
  2. stellar masses for the far clusters without the Big-Bang age cap;
  3. the Bullet Cluster, KiDS and Mistele's lensing in the project's distances;
  4. the release factor and release length from the bound-to-free transition (proposal 5).
