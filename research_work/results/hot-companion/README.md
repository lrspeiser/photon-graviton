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
