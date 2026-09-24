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

## 20. Round 10, 24 September 2026: first principles, and our own distances

The request: proceed to the next steps that solidify the theory and its first principles.
This round takes the first steps of proposal 1 (§19.5):
* part 1: where "the companion pulls with its amplitude" can come from (§20.1–20.2);
* part 2: how the companion of many pieces of matter must add up, with three candidate rules
  tested on the galaxies and the clusters (§20.5).

It also takes the far collisions and the strong lenses out of the expanding universe's distances
(§19.4), switches the suite to grade them that way (§20.3), asks how heavy the far clusters'
stars must then be (§20.6), and locks a first forecast (proposal 2, §20.4).

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
theories and are excluded. *(Part 2, §20.5: a condensate that conserves its quanta and their
momentum dilutes the pull exactly as independent waves do, so this analogy fails. What the data
require is sharper.)*

### 20.3 The far collisions and the strong lenses in the project's own distances

`code/collisions_v10.py` → `run-collisions-v10/collisions_v10.json`. Every length of the round-8
models is rescaled by the size factor, gas masses by the X-ray factor and star masses by the
luminosity factor, all at fixed angle and flux. Apertures are taken at the same angles. The law's
physics is unchanged. Each measured lensing mass is converted with its own paper's conventions,
from a check of each paper:

| Measurement | reference cosmology | background galaxies | lensing factor |
|---|---|---|---|
| MACS J0025 (Bradač et al. 2008) | flat, Ω_m 0.3, H₀ 70 | assumed z = 1.4 (one arc system at 2.38) | 1.38 |
| Abell 520 clumps (Jee et al. 2014; Clowe et al. 2012) | same | ⟨D_ls/D_s⟩ = 0.73 | 1.10 |
| Abell 520 inside 710 kpc (Mahdavi et al. 2007) | same | ⟨D_ls/D_s⟩ = 0.59 | 1.11 |
| El Gordo (Kim et al. 2021) | Planck 2015 (67.74, 0.3089) | ⟨D_ls/D_s⟩ = 0.254 (z_eff 1.31) | 1.55 |

El Gordo is now graded against Kim et al.'s aperture masses about the centre of mass (their
Fig. 11: 5.8, 14.9 and 20.1 × 10¹⁴ M☉ inside 0.5, 1 and 1.5 Mpc, ±9–14%; we use ±12%), not
against the sum of their two NFW fits. Star masses: El Gordo's come from colour fits that allowed
ages up to 7.0 Gyr (Menanteau et al. 2012); MACS J0025's from F814W light at a fixed M/L_K of 0.74;
Abell 520's from F814W light at M/L 2 (our choice in round 8, the Bullet's convention).

| | factors: size / stars / gas | lensing: ours ÷ measured (round 8, standard distances) | galaxy speeds, km/s (round 8) |
|---|---|---|---|
| MACS J0025, M(<300 kpc) SE / NW | 1.36 / 0.73 / 1.36 | 0.64 / 0.53, z −0.53 / −0.87 (0.78 / 0.61) | 707 vs 835 ± 59, z −2.17 (770) |
| Abell 520, six clumps inside 150 kpc | 1.08 / 0.80 / 1.00 | P1, P2, P3, P5 pass; P4 z −2.11 (−1.4); P6 z −2.64 (−2.3) | near P1, P2, P4, P5: rms z 3.18 (2.90) |
| Abell 520, inside 710 kpc | same | 0.76, z −2.17 (0.87, −1.2) | – |
| El Gordo inside 0.5 / 1 / 1.5 Mpc (apertures) | 1.58 / 0.71 / 1.68 | 0.75 / 0.75 / 0.86, z −2.12 / −2.08 / −1.19 (0.69 / 0.72 / 0.75 against the NFW sums) | NW 917, SE 822 vs 1,290 ± 134 / 1,089 ± 200 (944 / 839) |
| El Gordo, twice the published stars | same | 1.04 / 1.00 / 1.11 | NW 1,100, SE 974 |

**The suite now grades these collisions in the static law** (`regression/t_new_collisions.py`;
baseline saved with the full tier). Its tally moves from 62 pass, 8 close, 6 fail to **58 pass,
11 close, 7 fail**. The five grades that change are all far-cluster checks:
* from pass to close: MACS J0025's galaxy speeds, Abell 520's clump P4 and its mass inside 710 kpc,
  El Gordo's mass inside 1 Mpc;
* from close to fail: Abell 520's galaxy speeds.

All five move the same way, and for the same reason: the same light gives 20–29% less star mass
in our distances (§20.6).

**The six strong lenses** (the suite's check now uses the static convention):
* light = matter: −0.012 ± 0.023 dex (was −0.017 ± 0.024);
* the stars needed: 0.465 dex above Chabrier, that is 1.44–1.95 × Salpeter per lens (was
  1.05–1.35). Our distances give 1.48 times less starlight mass for the same light, and 1.11 times
  more lensing mass. *(Round 11 correction, §21.4: 1.25 times less starlight mass, the branch the
  lens code itself used; lensing masses 1.08–1.17 times more. No result changes.)*

**What this means.** With our own distances the measured lensing masses rise by 10–55%, while the
stars inferred from the same light fall by 20–30%. In our law the freely moving stars carry the
heat, so every far test now leans on the stellar masses. The published ones still carry
assumptions from the expanding-universe timeline: the age cap of §19.4, where older stars weigh
more for their light. §20.6 asks how much heavier they would have to be.

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

### 20.5 How the companion adds up: three rules, tested on the galaxies and the clusters (part 2)

`code/combination_rules_v10.py` → `run-combination-rules-v10/combination_rules_v10.json`.

**The local form of the pull.** For any companion field of one frequency, ψ = Im[Ψ(x) e^{iωt}], an
emitter locked a quarter cycle ahead of the field it sits in feels, exactly,

  ⟨F⟩ = (q₀/2) Im(Ψ* ∇Ψ) / |Ψ| = (q₀/2) |Ψ| ∇θ,    Ψ = |Ψ| e^{iθ}

That is, the local amplitude times the rate at which the local phase turns.
* Checked against the exact time average for four random configurations: they agree to 10⁻¹⁷ of
  the force.
* For one travelling wave |∇θ| = k, which is part 1's result.
* Where companion arrives from several directions at once, |Ψ|∇θ averages to the net current
  divided by the amplitude. Any companion that does not move with the net flow therefore dilutes
  the pull.
* Example: two coherent waves in opposite directions, amplitudes 1 and b. The mean pull is
  (1 − b)(2/π) K(4b/(1 + b)²) times one wave's: 0.95, 0.81, 0.53, 0.28 for b = 0.25, 0.5, 0.75,
  0.9. Round 3's rule would give √(1 − b²): 0.97, 0.87, 0.66, 0.44.
* A condensate that conserves its quanta and their momentum (a superfluid-like companion) has
  |Ψ|² = the total density and ∇θ ∝ the net current ÷ the density. It dilutes exactly as
  independent waves do (§20.2).

**Three rules.** In the law's terms, with S_N = G∫ρ/d² the plain total of every piece's Newtonian
pull, and S and g_hot the heat term's plain total and net flow:

| Rule | extra pull | what the companion does |
|---|---|---|
| 1. independent waves (§20.2) | √a \|g_N + g_hot\| / √(S_N + S) | waves pass through each other |
| 2. one stream carrying everything | √(a (S_N + S)) | always moves at full speed, all in one direction; nothing cancels |
| 3. round 3 (adopted) | √(a (\|g_N\| + S)) | ordered matter's opposing flows cancel; scrambled companion adds in full; one stream |

All three point along the net flow and use the same release factor.

**SPARC.** S_N at every measured radius, from each galaxy's surface densities:
* stars from the 3.6-μm profile at M/L 0.5, with an exponential thickness of 0.196 R_d^0.633 kpc
  (after Bershady et al. 2010, as in SPARC's own mass models);
* gas by inverting V_gas into thick annuli, by non-negative least squares with light smoothing,
  0.2 kpc thick;
* bulges spherical.

The reconstruction reproduces SPARC's own disk and gas pulls to medians of 0.990 and 0.999.
S_N/|g_N| has a median of 4.0 (16–84%: 2.2–8.3). With a and g_d refitted for each rule:

| Rule | typical miss at the round-9 constants | refitted a, g_d (m/s²) | typical miss, refitted (test, validation) |
|---|---|---|---|
| 3. round 3 | 15.85 km/s | 6.56e-11, 2.26e-10 | **15.85** (12.40, 19.21) |
| 1. independent waves | 24.73 | 4.58e-10, 8.69e-11 | 19.57 (16.21, 21.90) |
| 2. one stream | 30.82 | 8.71e-12, 7.66e-10 | 18.97 (15.76, 22.06) |

MOND (simple) scores 16.13. With disks twice as thick and 0.5-kpc gas the refitted misses are
15.85, 19.25 and 18.88: the verdict does not depend on the thickness.

**X-COP** (u refitted, a and g_d held): round 3 u 197.4 km/s, rms 0.227; independent waves 74.8,
0.382; one stream 221.9, 0.222. In these clusters the heat term is 5.6–26 times the net
Newtonian pull, so the treatment of cold matter hardly matters there; the treatment of the heat
does.

**The Milky Way** (local estimate as in §20.2, each rule with its own SPARC constants), circular
speed at 4 / 8.2 / 12 / 20 / 30 kpc:
* round 3: 198 / 211 / 207 / 195 / 185 km/s;
* independent waves: 188 / 202 / 216 / 237 / 249, rising outward;
* one stream: 205 / 206 / 190 / 162 / 143, falling steeply.

**Verdict: only round 3 passes both the galaxies and the clusters.**

| Rule | SPARC, typical miss (MOND 16.1 km/s) | X-COP, typical miss | Milky Way (local) |
|---|---|---|---|
| 1. independent waves | 19.6 km/s | 47% | too fast outside |
| 2. one stream | 19.0 km/s | 25% | too slow outside |
| 3. round 3 | **15.9 km/s** | **25%** | 211 km/s at the Sun, 185 at 30 kpc |

(Typical miss for X-COP: e^rms − 1, as elsewhere in the blog.) So the data, not the choice of
postulates, fix three properties of the companion:
1. **Ordered matter's opposing flows cancel.** Only the net flow survives (SPARC excludes rule 2).
2. **Scrambled companion adds its full strength** (X-COP excludes rule 1's dilution).
3. **Matter reads it as one stream whose phase turns at the full rate.** By the local form, any
   slower turning dilutes the pull. That happens for waves that pass through each other and for
   a superfluid-like companion.

A picture with these properties: waves that always travel at full speed and annihilate when they
meet head-on. Fronts in active (excitable) media behave this way, such as chemical waves, or waves
of activity in heart tissue. In those media, though, the medium sets the waves' strength. Here the
strength must still fall off away from its source. A model that does both is the next toy.

### 20.6 How heavy must the far clusters' stars be?

`code/collisions_star_sweep_v10.py` → `run-collisions-v10/star_sweep_v10.json`. The three
collisions rerun with every star mass multiplied by one common factor, on top of the static-law
conversion (z-scores; galaxy speeds in km/s):

| Star factor | MACS J0025 | Abell 520 | El Gordo |
|---|---|---|---|
| 1 | M300 −0.53 / −0.87; speeds 707 (z −2.17) | P4 −2.11, P6 −2.64, M710 −2.17; speeds rms z 3.18 | apertures −2.12 / −2.08 / −1.19; speeds 917 (−2.78), 822 (−1.33) |
| 1.4 | −0.27 / −0.67; 793 (−0.71) | P2 +1.45, P4 −1.09, P6 −2.28, M710 −1.09; speeds 2.74 | −1.18 / −1.27 / −0.38; 993 (−2.22), 884 (−1.02) |
| 1.8 | 0.00 / −0.46; 872 (+0.63) | P2 +2.98, P4 −0.05, P6 −1.90, M710 +0.05; speeds 2.37 | −0.20 / −0.43 / +0.48; 1,060 (−1.68), 945 (−0.72) |
| 2 | – | – | +0.30 / +0.01 / +0.93; 1,100 (−1.42), 974 (−0.58) |

* With a common factor of 1.4, the 17 graded checks of these three clusters would read 14 pass
  and 3 close (now 9 pass, 7 close, 1 fail), and the whole suite **63 pass, 7 close, 6 fail**.
  That is better than before the switch of distances.
* El Gordo prefers 1.6–2.0.
* Abell 520's clumps pull different ways. P2 is already heavy enough, while P6 and the galaxy
  speeds near P1 want more. That points at its clump-by-clump star masses.

**Is 1.4 plausible?** The factor 1.4 roughly undoes the static-law conversion of the star masses
(1/0.73, 1/0.80 and 1/0.71). The published masses carry assumptions that lower them:
* El Gordo's colour fits cap the stars' ages at the Big Bang's timeline (7.0 Gyr);
* MACS J0025's use one fixed light-to-mass ratio;
* Abell 520's use a ratio of 2 that is on the low side for old stars.

In standard stellar-population models (e.g. Bruzual & Charlot 2003), an old population's
mass-to-light ratio grows roughly as its age to the power 0.5 (near-infrared) to 0.85 (visible).
Stars 10–13 billion years old, like those of nearby cluster ellipticals, would weigh roughly
1.3–1.7 times more for the same light than stars 6–7 billion years old. So the factor is within
reach. It is not yet a result: the star masses have to be refitted from the photometry with no
age cap.

*Round 11: withdrawn. The calibrating clusters' stars had been counted in projection and u measured
in the standard distances; measured consistently, the far clusters need no common star factor
(§21.3–21.4). Distant cluster stars are observed to be younger, not older (§21.2).*

### 20.7 Where this leaves the first principles

* **Derived:**
  * the pull proportional to the companion's amplitude, from a local coupling, at every distance
    (part 1);
  * its local form: amplitude × the rate at which the phase turns (part 2);
  * the direction along the net flow;
  * no energy exchange unless the pull is paid by coherent emission.
* **Required of the microscopic model:**
  * from the model itself (part 1): an active medium (a quarter cycle ahead), spectral locking
    (universality), and a slow phase (bookkeeping);
  * from the data (part 2): ordered opposing flows cancel, scrambled companion adds in full, and
    matter reads one stream at the full phase rate.
* **Excluded by the data:**
  * independent waves;
  * a superfluid-like companion that conserves momentum;
  * one stream in which nothing cancels.
* **Next:**
  1. a toy of full-speed waves that annihilate head-on while their strength falls off from the
     source, graded on SPARC, X-COP and the Milky Way with the code of §20.5;
  2. star masses for the far clusters from their photometry with no age cap (§20.6);
  3. the Bullet Cluster, KiDS and Mistele's lensing in the project's distances;
  4. the release factor and release length from the bound-to-free transition (proposal 5).

## 21. Round 11, 24 September 2026: why orderly matter adds up like Newton, the stars audited, and everything in our own distances

The request: proceed with the next steps of round 10 (§20.7), which were:
* a model for the companion as one stream: full-speed waves that cancel head-on while weakening away
  from the source;
* star masses for the far clusters without the Big-Bang age cap;
* the Bullet Cluster, KiDS and X-COP in the project's own distances.

This round:
1. derives, from energy bookkeeping, why orderly matter's companion must add up like Newton's pulls,
   and measures how tightly the data demand it (§21.1). The dynamical toy of the waves is not built
   yet; it is the first next step (§21.7);
2. audits how every cluster's star masses were made (§21.2). Two corrections follow: X-COP's stellar
   profiles are projected masses, which rounds 1–10 used as spherical ones, and MACS J0025's are on a
   heavier (Salpeter) basis. Distant cluster stars turn out younger, not older;
3. moves every distance-dependent test into the static distance law, which recalibrates the
   companion's speed to u = 162.6 km/s (§21.3). These constants are now the adopted law;
4. follows what that does to the far collisions, the Bullet Cluster and the strong lenses (§21.4);
5. finds that galaxy lensing (KiDS) now measures two things: the distance law's scale and the
   companion's speed (§21.5);
6. compares the two static geometries in the repository (§21.6).

### 21.1 Why orderly matter's companion adds up like Newton's pulls

`code/companion_flow_v11.py` → `run-companion-flow-v11/companion_flow_v11.json`.

**The argument.**
* In a steady state the companion's energy flux J obeys div J = ℓρ: every watt emitted flows out, and
  nothing is lost on the way.
* By round 10's local form (§20.5), matter feels the companion's full pull only where the companion
  passes as one stream at the full phase rate. Its energy density is then |J|/u, and the pull is
  proportional to √|J|.
* If the flow also has no whirlpools (curl J = 0), J is unique: J = (ℓ/4πG)(−g_N), Gauss's geometry.
* So the stream's intensity is |g_N|, and the pull is ∝ √|g_N|. That is round 3's rule for ordered
  matter, obtained from one stream plus energy conservation.

**The alternative is not a consistent flow.** "One stream carrying everything" (round 10's rule 2,
intensity S_N along the net pull) would carry more power out of a closed surface than the matter inside
it emits:

| Where | power carried out ÷ power emitted inside |
|---|---|
| uniform sphere, at 0.25 / 0.5 / 1 / 2 radii | 11.8 / 5.5 / 1.5 / 1.06 |
| Plummer sphere, at 0.25 / 0.5 / 1 / 2 scale radii | 8.1 / 4.2 / 2.3 / 1.5 |
| Hernquist sphere, at 0.25 / 0.5 / 1 / 2 scale radii | 3.3 / 2.5 / 2.0 / 1.6 |
| the Milky Way model (spheres of 4 / 8.2 / 12 / 20 / 30 kpc) | 1.70 / 1.36 / 1.22 / 1.11 / 1.06 |

Independent waves (rule 1) conserve energy (their net flux is J), but they dilute the pull. Round 3 is
the one rule that is both lossless and undiluted.

**How tightly do the data demand it?** A three-parameter family contains all the rules:

  I = w|g_N| + (1 − w)S_N + hS + (1 − h)|g_hot|,  D = |g_N + g_hot| / (S_N + S),  extra = release · √(aI) · D^γ

Round 3 is (w, h, γ) = (1, 1, 0). The scans run through it, refitting a and g_d on SPARC and u on
X-COP at every point. X-COP is taken as the suite now takes it, in our distances with deprojected
stars (§21.2–21.3), and S_N comes from the surface densities of §20.5. The preferred value comes from
400 bootstrap resamplings:

| Dial | what it measures | SPARC typical miss along the scan | X-COP rms along the scan | preferred |
|---|---|---|---|---|
| w = 0 … 1 | how completely orderly flows cancel | 18.89 → 15.94 km/s | 0.212 → 0.221 (u 60 → 163 km/s; a 7.5 times smaller at w = 0) | SPARC: w = 1 (every resampling); X-COP leans to 0 |
| γ = 0 … 1 | how much the pull is diluted | 15.94 → 27.21 km/s | 0.221 → 0.505 (u 163 → 305 km/s) | SPARC 0.05 (0–0.1); X-COP 0 (every resampling) |
| h = 0 … 1 | how much the heat adds unsigned | 15.93 → 15.94 km/s | 0.288 → 0.221 (u 97 → 163 km/s) | X-COP: h = 1 (every resampling); SPARC flat |

So:
* **the galaxies demand complete cancellation** of orderly flows (w = 1; at w = 0.95 the miss is
  already 16.11 km/s). X-COP alone leans the other way by 0.009 in rms, but only through a trade:
  with SPARC's refit at w = 0, a is 7.5 times smaller and u falls to 60 km/s;
* **galaxies and clusters both demand an undiluted pull** (γ ≤ 0.1);
* **the clusters demand the heat unsigned** (h = 1 in every resampling). With the heat cancelling
  like arrows, u falls to 97 km/s and the rms rises to 0.288.

The first version of this scan, run before the X-COP correction of §21.2, gave the same verdicts,
with h ≥ 0.5 instead of h = 1.

### 21.2 How each cluster's star masses were made: an audit, and two corrections

A literature audit, with every number traced to its paper (`literature/star_mass_audit_v11.md`), put
every star mass the law uses on one basis. The law's u is calibrated on X-COP, whose stellar masses
assume a Chabrier IMF (Ghizzardi et al. 2021, arXiv:2007.01084, Sect. 4.1). Everything must therefore
be compared on that basis.

**Correction 1: X-COP's stellar profiles are projected.**
* The release's cumulative stellar profiles are masses "within a projected radius": a cylinder along
  the line of sight. The paper converts them to the R500 sphere by multiplying by 0.75, using a gNFW
  galaxy distribution (c = 0.72, α = 1.64; van der Burg et al. 2015).
* 0.75 × the release profile at R500 reproduces the paper's Table 3 for six of the seven clusters
  within 0.4–4% (A644: 12%).
* Rounds 1–10 (`run.load_xcop`) used the projected profiles as spherical masses. That overstated the
  stars by 1.3 at R500 and by up to about 2 in the middle. It also put the star/gas ratio at R500 at
  0.035–0.074, where the paper's spherical values give 0.023–0.049.
* `xcop_static_v11.deproject_xcop` fits each profile with two parts:
  * a BCG (Hernquist, a = 15 kpc);
  * satellites with that gNFW, R200 = R500/0.65, cut along the line of sight at 1.3 R200 so that
    their sphere/cylinder ratio at R500 is the paper's 0.75.

  It then replaces the profile by the fit's spherical mass. The five clusters without optical data
  take the median spherical star/gas ratio of the seven, as before.
* Check against Table 3 (M_star inside the R500 sphere, 10¹² M☉):

  | | A1795 | A85 | A644 | A2319 | ZW1215 | A2029 | A2142 |
  |---|---|---|---|---|---|---|---|
  | deprojected here | 3.17 | 2.27 | 4.19 | 5.30 | 3.17 | 6.82 | 7.02 |
  | Table 3 | 3.02 | 2.10 | 3.70 | 5.11 | 3.34 | 6.51 | 6.97 |

  The deprojected star/gas ratio at R500 is 0.027–0.056 (median 0.049).

**Correction 2: MACS J0025's star masses are Salpeter-based.**
* Bradač et al. 2008 convert F814W light to rest-frame K with a non-evolving elliptical template, and
  multiply by M/L_K = 0.74 ± 0.30 "following Drory et al. (2004)".
* That value is Drory et al.'s mean for massive galaxies at z = 0.5 (their Table 1), whose fits assume
  a Salpeter IMF.
* Bradač et al. do not state the IMF; this is inferred. We convert with the SLACS code's convention,
  × 10^−0.25 (`collisions_v10.chabrier_basis`).

**The rest are on the same basis:**
* El Gordo: Menanteau et al. 2012, Chabrier.
  * Their SED grid's ages run to 7.0 Gyr (the universe's age at z = 0.87 in their cosmology is
    6.6 Gyr), and they do not report fits piling up at that edge.
  * An independent estimate agrees within their stated factor of two: Hilton et al. 2013 give
    10.8 × 10¹² M☉ inside R500 on a Chabrier basis, against Menanteau's 13.1 × 10¹² inside r200.
* The Bullet: Clowe et al. 2006, M/L_I = 2 after Kauffmann et al. 2003, a Kroupa IMF close to Chabrier.
  They call it an upper limit, since no colour selection was made.
* Abell 520: no paper gives a stellar M/L for its galaxies. We keep M/L = 2 on F814W light (the
  Bullet's convention).

**Star masses without the age cap: distant cluster stars are younger (a correction to §20.6).**
Round 10 suggested the far clusters' stars could be older, and so heavier, than their published fits
allow. The observations say otherwise, and not through any distance law:
* at fixed velocity dispersion, cluster ellipticals' rest-frame U − V is 0.24 ± 0.02 mag bluer at
  z = 0.83 than in Coma (Holden et al. 2010);
* the median 4000 Å break of massive passive galaxies falls from 1.99 at z = 0.16 to 1.71 at
  z = 1.02 (Moresco et al. 2012);
* the fundamental plane gives d log(M/L_B)/dz = −0.555 ± 0.042 for cluster galaxies (van Dokkum & van
  der Marel 2007; Holden et al. 2010: −0.60; Saglia et al. 2010: −0.54). Its sizes and brightnesses
  depend on distances, but its trend agrees with the colours.

Distant cluster galaxies' stars are younger and lighter per unit light, so the age cap does not hide
extra star mass. With u calibrated consistently (§21.3), the far clusters no longer need heavier stars
(§21.4).

**The IMF of massive ellipticals, for context (and for §21.7):**
* about 1.7–1.9 × Chabrier (≈ Salpeter) at σ = 250–300 km/s from lensing with dynamics (Treu et al.
  2010, arXiv:0911.3392: 1.78; Auger et al. 2010, arXiv:1007.2409: 1.4–1.9);
* 0.95–1.02 × Salpeter from dynamics (Cappellari et al. 2013);
* about 1.85 × Chabrier from spectral features (Conroy & van Dokkum 2012, arXiv:1205.6473). This
  method needs neither dark matter nor any cosmology.

### 21.3 X-COP in our own distances, and the companion's speed recalibrated

`code/xcop_static_v11.py` → `run-xcop-static-v11/xcop_static_v11.json`.

The release gives every profile for flat LCDM (H0 = 70, Ω_m = 0.3). At the clusters' redshifts
(0.047–0.090), at fixed angle and flux, the static law gives:
* radii and hydrostatic masses × 0.970–1.000;
* gas × 0.886–0.917;
* stars × 0.841–0.859.

Temperatures and the stars' Jeans speeds are unchanged. u refitted with a and g_d held:

| X-COP sample | u (km/s) | rms ln(M_hydrostatic/M_predicted) |
|---|---|---|
| the release's units, projected stars (rounds 1–10) | 197.4 | 0.227 |
| our distances, projected stars | 179.6 | 0.227 |
| the release's units, deprojected stars | 179.4 | 0.221 |
| **our distances, deprojected stars** | **163.4** | **0.220** |

The stars carry the heat, with weight k = 3σ²/u². Fewer stars at the same speeds must be balanced by a
slower companion, and the fit gets slightly better.

A joint refit alternates three times: a and g_d on SPARC, where u enters only through the bulges'
heat, and u on X-COP. It gives:

| | a (m/s²) | g_d (m/s²) | u (km/s) | SPARC | X-COP rms |
|---|---|---|---|---|---|
| round 9 (the release's units, projected stars) | 6.561e-11 | 2.262e-10 | 197.4 | 15.85 km/s | 0.227 |
| **round 11 (our distances, deprojected stars)** | **6.547e-11** | **2.107e-10** | **162.6** | 15.94 km/s | 0.221 |

**Adopted.** `regression/law_config.load_law()` now returns these constants ('round11'); 'round9'
still loads the old ones. It is the same law, with its three constants measured in the project's own
distances and on correctly deprojected stars. Two consequences:
* the reach, u × 13 billion years, becomes 2.16 Mpc instead of 2.62;
* the fresh companion around stopped gas grows at 166 kpc per billion years instead of 202.

### 21.4 The far collisions, the Bullet Cluster and the strong lenses with the adopted constants

**The far collisions** (`collisions_v10`, static distances, published stars on the Chabrier basis):

| | round 10 (static distances; u = 197, fitted in the standard ones) | round 11 (u = 163, fitted in our distances) |
|---|---|---|
| MACS J0025: lensing SE / NW; speeds | z −0.53 / −0.87; 707 km/s (z −2.17) | z −0.65 / −0.97; 665 km/s (z −2.88) |
| Abell 520: P4, P6, 710 kpc; speeds | −2.11, −2.64, −2.17; rms z 3.18 | **−1.01**, −2.28, **−1.03**; rms z **2.70** |
| El Gordo: apertures 0.5 / 1 / 1.5 Mpc; speeds NW / SE | −2.12 / −2.08 / −1.19; 917 / 822 km/s | **−1.07 / −1.18** / −0.28; **1,002 / 892** km/s (z −2.15 / −0.98) |

Of the 17 graded checks of these clusters, 12 now pass, 4 are close and 1 fails (round 10: 9, 7, 1).
So round 10's "the far clusters need 1.4 times their stars" came mostly from mixing conventions. The far
clusters were in our distances and their stars on one basis, but the companion's speed had been
calibrated in the standard distances on projected stars.

**MACS J0025 is the exception, for two reasons.**
* Its stars, now converted to the Chabrier basis, are 1.78 times lighter than published. Its galaxy
  speeds fall to 665 km/s (z −2.88); with the published masses they would be 804 (z −0.53).
* Its NW lensing peak moves from the galaxies (59 kpc from them) to the gas (209 kpc; the galaxies
  lie 231 kpc from the gas). `code/macs_peak_scan_v11.py` → `run-collisions-v10/macs_peak_scan_v11.json`
  shows why, and how close it is:

  | time since closest approach | stars × 0.56 (Chabrier) | × 0.75 | × 1 (as published) |
  |---|---|---|---|
  | 0.26 Gyr | 73 kpc (pass) | 60 | 55 |
  | 0.35 Gyr | 71 (pass) | 59 | 54 |
  | 0.5 Gyr (the suite) | **209, at the gas (fail)** | 58 | 54 |

  Along the 231 kpc from the NW galaxies to the gas, the lensing map is a nearly flat ridge. At
  0.5 Gyr with Chabrier-basis stars, the map at the galaxies is 9% below its value at the gas. A younger collision
  (≤ 0.35 Gyr; Bradač et al. 2008: closest approach "a few 10⁸ years ago") or 1.33 times the stars
  returns the peak to the galaxies. The suite keeps grading at 0.5 Gyr, the value fixed in round 8, so
  this check now fails.

**The Bullet Cluster** (`code/bullet_static_v11.py` → `run-bullet-static-v11/bullet_static_v11.json`).
At z = 0.296, lengths are × 1.14, gas × 1.08, stars × 0.78, and lensing masses inside a fixed angle
× 1.18; κ is compared as measured. The suite now grades the Bullet this way:

| | round 9 (standard distances) | round 11 (our distances) | measured, converted |
|---|---|---|---|
| stars in the main cluster's outskirts (10¹² M☉) | 6.23 (target 3.9–6.7) | 3.32 | 3.0–5.2 |
| κ on the main / sub galaxies | 0.675 / 0.140 | 0.715 / **0.259** | ≥ 0.36 ± 0.06 / ≥ 0.20 ± 0.05 |
| extra κ on the main / sub gas | 0.033 / 0.042 | 0.046 / 0.048 | 0.05 ± 0.06 / 0.02 ± 0.06 |
| lensing peaks from their galaxies, main / sub | 9 / 16 kpc | 14 / 21 kpc | on the galaxies |
| mass inside 250 kpc (286 in ours), main (10¹⁴ M☉) | 2.40 (target 2.50–2.80) | 2.83 (z −0.61) | 2.94–3.29 |
| the same, sub | 0.94 (target 2.00–2.30; z −5.30) | 1.29 (z −4.53) | 2.35–2.70 |

The smaller half gains 37% but still has about half the measured mass inside 286 kpc.

**The six strong lenses**, with the adopted constants:
* light = matter −0.032 ± 0.023 dex;
* stars needed 0.445 dex above Chabrier, 1.37–1.85 × Salpeter.

*Correction to §20.3:* our distances make the lenses' stars 1.25 times lighter (1.23–1.27), not 1.48
times. The number quoted there was the energy-loss-only branch. The lens code itself used the branch
with arrival-rate stretching throughout, so no result changes. Their lensing masses are 1.08–1.17
times heavier.

### 21.5 Galaxy lensing in our own distances: it now measures the distance law's scale and the companion's speed

*Round 12 correction (§22.2–22.3): neither reading held. The early/late gap came from two assumed star
speeds, and with the constants refitted at each α the common level does not measure α.*

`code/kids_static_v11.py` → `run-kids-static-v11/kids_static_v11.json`. The KiDS lenses lie at mean
z = 0.25, with sources at an effective z = 0.75, and the paper's distances are flat LCDM
(Ω_m = 0.2793, H0 = 70). At z = 0.25 our distances give:
* g_bar × 0.64, which is (1 + z)⁻²: the static law's surface brightness dims as (1 + z)² instead of
  (1 + z)⁴;
* g_obs × 0.934, the ratio of critical densities.

Median log10(observed/predicted):

| | all | blue | red | disks (Sérsic n < 2) | bulges (n > 2) | GAMA | gap (obs 0.153 ± 0.04) |
|---|---|---|---|---|---|---|---|
| standard distances, round-9 constants | +0.024 | −0.005 | +0.021 | +0.032 | −0.023 | −0.032 | 0.178 |
| our distances, round-9 constants | +0.098 | +0.077 | +0.092 | +0.107 | +0.048 | +0.043 | 0.181 |
| **our distances, round-11 constants** | **+0.063** | **+0.077** | +0.040 | **+0.107** | −0.006 | +0.007 | **0.234** |

Mistele et al.'s lensing speeds, spirals: observed ÷ predicted 1.40, 1.32, 1.22, 1.17 from 50 to 300 kpc
(rms z 4.36); ellipticals 1.09, 0.98, 1.01 (0.95). Two separate things have moved.

**The level: the distance law's scale.** g_bar does not depend on α: the stars and the area both scale
as 1/α². The critical density, and with it every measured g_obs, scales as α. Holding the constants:

| α (1 + z = e^{αD}) | all | blue | red | disks | bulges | GAMA |
|---|---|---|---|---|---|---|
| × 0.8 | −0.034 | −0.020 | −0.057 | +0.010 | −0.103 | −0.090 |
| × 0.9 | +0.017 | +0.031 | −0.006 | +0.061 | −0.052 | −0.039 |
| × 1 (adopted) | +0.063 | +0.077 | +0.040 | +0.107 | −0.006 | +0.007 |
| × 1.1 | +0.104 | +0.119 | +0.081 | +0.149 | +0.035 | +0.049 |

The spirals' excess (+0.077) does not depend on u: their heat weight is small. The whole sample
centres at α ≈ 0.87 times the adopted value, an "H0-like" scale of 65 instead of 74.6 km/s/Mpc; the
spirals alone at 0.84 (63), the disks at about 0.78. The adopted α comes from the project's own
light-transport fit (`research_work/results/joint-light-forward`). SPARC's Hubble-flow distances,
about half its galaxies, would move with α too, which would change a. So the consistent test is a
joint one (§21.7). The geometry barely matters here (the "metric" rows of `alpha_scan` sit
0.004–0.014 dex lower).

**The gap between ellipticals and spirals: the companion's speed.** The red lenses' stars (160 km/s)
carry heat k = 3σ²/u², which is 2.9 at u = 162.6 and 2.0 at 197.4. The model's gap rises from 0.18 to
0.234 against the measured 0.153 ± 0.04. The gap does not depend on α. So the clusters (X-COP, with
correctly counted stars) want a slower companion than single ellipticals do. §21.7 lists the
candidates, first among them the ellipticals' IMF.

### 21.6 The two static geometries

`code/distance_variants_v11.py` → `run-distance-variants-v11/distance_variants_v11.json`. Both keep
1 + z = e^{αD} and D_L = (1 + z)D. They differ in how an angle becomes a size:
* "fixed": D_A = D, the fixed-material transport branch, used since round 10;
* "metric": D_A = D/(1 + z), the material-coasting geometry of
  `five_candidate_tests/prior/conformal_action_derivation_derivation.md`, eq. 13.

Each geometry refits its own constants (X-COP deprojected):

| | fixed (adopted) | metric |
|---|---|---|
| u; SPARC; X-COP rms | 162.6 km/s; 15.94; 0.221 | 161.2 km/s; 15.94; 0.221 |
| KiDS all / blue / disks; gap | +0.063 / +0.077 / +0.107; 0.234 | +0.055 / +0.063 / +0.100; 0.233 |
| Mistele, spirals / ellipticals (rms z) | 4.36 / 0.95 | 4.25 / 0.90 |
| SLACS: light = matter; stars needed | −0.032 ± 0.023; 1.37–1.85 × Salpeter | +0.000 ± 0.025; 1.24–1.61 × Salpeter |
| MACS J0025: lensing SE / NW; speeds | −0.65 / −0.97; z −2.88 | −0.72 / −1.09; z −3.43 |
| Abell 520: P4, P6, 710 kpc; speeds | −1.01, −2.28, −1.03; rms z 2.70 | −1.12, −2.47, −1.43; 2.77 |
| El Gordo: apertures 0.5 / 1 / 1.5 Mpc; speeds NW / SE | −1.07 / −1.18 / −0.28; z −2.15 / −0.98 | −2.33 / −2.52 / −1.89; z −2.75 / −1.36 |
| Bullet: κ sub; mass inside 286 kpc, main / sub | 0.259; z −0.61 / −4.53 | 0.117 (z −1.66); z −0.71 / −5.61 |

The galaxy lenses lean slightly (0.01 dex) towards "metric"; the clusters clearly prefer "fixed". The
fixed geometry stays.

### 21.7 The suite, where this leaves things, and next

**The suite with the adopted constants** (full tier, baseline saved): **57 pass, 9 close, 10 fail** (round 10's baseline:
58 pass, 11 close, 7 fail).
* **Improved (5):** Abell 520's clump P4 and its mass inside 710 kpc (close → pass), its galaxy
  speeds (fail → close); El Gordo's lensing inside 0.5 and 1 Mpc (close → pass).
* **Regressed (6):** KiDS all (pass → close), blue and disks (pass → fail), the early/late gap
  (pass → close); Mistele's spirals (close → fail); MACS J0025's NW peak (pass → fail).
* **Still failing, as before (6):** the five faint dwarfs and the Bullet's smaller half.
* Nearby tests move by at most a few percent: SPARC 15.85 → 15.94 km/s; X-COP 0.227 → 0.221
  (held out 0.244 → 0.236); the Sun's speed 211.2 → 210.7 km/s; Cassini's Q2 4.6 → 4.5 × 10⁻²⁷
  s⁻²; wide binaries at 20,000 AU 1.093 → 1.083 (at 7,000 AU 1.039 → 1.035, inside the locked
  forecast's window of 1.03–1.05 and 1.08–1.10); the collision stack's β 0.027 → 0.016.

**Where this leaves the first principles.**
* **Derived:** the ordered companion's rule, from energy conservation plus one stream. The "one stream
  of everything" rule is excluded not only by the data (§20.5) but by energy conservation.
* **Required:** the companion's flow has no whirlpools. A dynamical model still has to supply this.
* **Measured more cleanly:** the companion's speed, u = 162.6 km/s, now in the project's own distances
  and on correctly counted stars.

**Next:**
1. **The dynamical toy** of full-speed waves (not built this round). A concrete route: a companion that
   scatters off itself relaxes towards a flow without whirlpools, as heat does when it diffuses, and
   such a flow's steady flux follows Newton's field lines. The toy would show whether a local rule
   gives the curl-free condition of §21.1, graded with the code of §20.5.
2. **The distance law's scale, from galaxy lensing.** Refit α jointly on KiDS, SPARC's Hubble-flow
   distances and the light-transport fit that set it.
3. **The ellipticals' heat.** Clusters want u = 163 km/s, single ellipticals' lensing gap wants about
   200. Candidates:
   * the IMF: spectra say massive ellipticals carry 1.7–1.9 times Chabrier's star mass, which would
     change X-COP's calibration and KiDS's ellipticals together;
   * the speeds assumed for KiDS's red lenses (160 km/s).
4. **MACS J0025**: its collision age and star masses decide its NW peak.
5. **The Bullet's smaller half**, still about half its measured lensing mass.

## 22. Round 12, 24 September 2026: the next steps, one at a time

The request: proceed with round 11's next steps (§21.7), keeping the blog and `main` current as each one
lands. In the order done:
1. MACS J0025's collision dated from its own shock fronts (§22.1);
2. the lenses' heat measured instead of assumed (§22.2);
3. the distance law's scale, fitted jointly (§22.3);
4. a dynamical model of the companion (§22.4).

### 22.1 MACS J0025's age, from its shock fronts

`code/macs_timing_v12.py` → `run-collisions-v10/macs_timing_v12.json`.

Round 11 found that MACS J0025's NW lensing peak sits on its galaxies only if the collision is younger
than about 0.35–0.40 Gyr (§21.4). The suite used 0.5 Gyr, chosen in round 8 between two kinds of
estimate. This step measures the age with the collision's own clocks, both independent of our law:
* **Shock fronts.** Riseley et al. 2017 (A&A 597, A96; arXiv:1611.01273) found two radio relics, the
  usual tracers of merger shocks, NW and SE of the centre and perpendicular to the merger axis.
  * We measured the centroids of their 325 MHz contours (the 5, 7 and 9σ levels of their low-resolution
    map) from the centre they adopt: NW 23–25″, SE 49–51″ (150–158 and 315–326 kpc at their
    6.416 kpc/″).
  * The relics sit almost exactly where the galaxies are: 26″ (NW) and 56″ (SE) from the X-ray peak
    (Bradač et al. 2008). So the shocks have not yet outrun the galaxies: a young collision.
  * The NW relic's spectral index (α < −1.3) limits the shock to Mach < 1.87 (their eq. 1). With the
    gas's sound speed of about 1,300 km/s (Bradač et al. 2008), a shock launched at closest approach has
    taken at least d/(Mach c_s): 0.06–0.08 Gyr (NW) and 0.13–0.16 Gyr (SE), for Mach 1.87–1.5.
* **Separation over speed.** The galaxy concentrations are 540 kpc apart (Bradač et al.'s units), and
  the collision speed is about 2,000 km/s (their Sect. 3.2): 0.26 Gyr, their "a few 10⁸ years".
* In our distances sizes are 1.36 times larger and speeds are unchanged, so the clocks read
  **0.08–0.22 Gyr (shocks) and 0.36 Gyr (separation)**.
* The post-starburst clock (Ma et al. 2010: 0.5–1 Gyr since first core passage) dates when star
  formation was triggered and quenched, which can start before closest approach. It is listed, not used,
  for the time since the gas stopped.

**The suite now uses 0.3 Gyr** (`regression/t_new_collisions.py`), inside the dynamical range. The model
at nearby ages (round-11 constants, Chabrier-basis stars):

| time since closest approach | fresh companion around the gas | NW peak from its galaxies | grade |
|---|---|---|---|
| 0.3 Gyr (adopted) | 50 kpc | 72 kpc | pass |
| 0.4 Gyr | 67 kpc | 70 kpc | pass |
| 0.45 Gyr | 75 kpc | 207 kpc, at the gas | fail |

The lensing masses (1.94 and 1.72 × 10¹⁴) and the galaxy speeds (665 km/s) do not depend on the age.
So the NW peak's knife edge (0.40–0.45 Gyr) lies beyond every dynamical clock.

**The suite** (full tier, baseline saved): **58 pass, 9 close, 9 fail** (round 11: 57, 9, 10). The only
grade that moves is MACS J0025's NW peak, fail → pass.

### 22.2 The lenses' heat, measured

`code/lens_heat_sdss_v12.py` → `data/lens_heat_sdss_v12.json`; `code/kids_heat_v12.py` →
`run-kids-heat-v12/kids_heat_v12.json`.

Round 11 found KiDS's early/late gap too large with the new u (0.234 against 0.153 ± 0.04; §21.5) and
asked whether the ellipticals' heat or their IMF was to blame. The KiDS comparison had assumed one
number per class: 160 km/s for every star of a red or bulge-dominated lens, and k = 0.1 for every blue
or disk-dominated one. Both are now measured.

**The measurement.** SDSS DR17 spectra of 119,000 galaxies at 0.1 < z < 0.25 with MPA-JHU stellar
masses (Kroupa, close to Chabrier), split as Brouwer et al. 2021 split KiDS: observed u − r at 2.5
(red/blue), and fracDeV_r at 0.5 as a proxy for their Sérsic n = 2.
* The hot part of each galaxy is its de Vaucouleurs fraction, at the fibre dispersion corrected to one
  effective radius (σ_R/σ_e = (R/R_e)^−0.066, Cappellari et al. 2006). The rest is a disk at 30 km/s.
* Median σ_e of red galaxies: 125, 142, 161 and 186 km/s at log M* = 10.0, 10.5, 10.75 and 11.0. Blue
  galaxies of the same masses have bulges of 73–138 km/s, with de Vaucouleurs fractions of 0.1–0.5.
* Mean heat weight over log M* = 10.3–10.9 (the KiDS stacks' typical 10.6), with u = 162.6 km/s:

  | | red | blue | bulge-dominated | disk-dominated |
  |---|---|---|---|---|
  | measured (SDSS) | 2.14 | 0.55 | 2.15 | 0.30 |
  | assumed until round 11 | 2.90 | 0.10 | 2.90 | 0.10 |

  The red lenses are cooler than assumed. The blue lenses are warmer, because their bulges carry heat.

**The result** (static distances, round-11 constants; median log10(observed/predicted)):

| | all | blue | red | disks | bulges | GAMA | colour gap (obs 0.153) | Sérsic gap (obs 0.154) | Mistele spirals / ellipticals (rms z) |
|---|---|---|---|---|---|---|---|---|---|
| assumed heat (round 11) | +0.063 | +0.077 | +0.040 | +0.107 | −0.006 | +0.007 | 0.234 | 0.234 | 4.36 / 0.95 |
| **measured heat** | +0.070 | **+0.021** | +0.081 | +0.077 | +0.037 | +0.015 | **0.132** | **0.163** | 2.98 / 2.52 |

* **Both gaps now match** (0.132 and 0.163 against 0.153 and 0.154). The gap was never a problem of u.
  It came from two assumed numbers.
* What remains is a common level. Every sample sits 0.02–0.08 dex above the prediction, and Mistele's
  lensing speeds of both types exceed ours. §22.3 asks whether the distance law's scale sets it.

**The IMF of giant ellipticals, checked.** Spectra show bottom-heavy stars mainly in the centres of the
most massive ellipticals. Averaged over whole galaxies, the stellar M/L rises by at most 0.09 dex, and
only above 2 × 10¹¹ M☉ (Domínguez Sánchez et al. 2019, as used by Brouwer et al. 2021 for these
lenses). So:
* it does not reach the KiDS lenses (all below 10¹¹);
* in X-COP it touches only the BCGs (2–19% of the stars inside R500), moving u by at most 2%;
* it trims the SLACS lenses' star need from 0.445 to about 0.36 dex above Chabrier.

It is not the lever for u or KiDS.

**Adopted in the suite** (`regression/t_lensing.py`): each KiDS sample and each of Mistele's mass bins
takes its measured heat, and the Sérsic gap is graded too (a new check; 90 checks, 77 graded).

**The suite** (full tier, baseline saved): **59 pass, 10 close, 8 fail** of 77 graded (step 1: 58, 9, 9 of
76). Moved: KiDS blue fail → pass, the colour gap close → pass, Mistele's spirals fail → close, and the
new Sérsic gap passes; KiDS red pass → fail (+0.081) and Mistele's ellipticals pass → close (2.52), as
they join the common level.

### 22.3 The distance law's scale, fitted jointly

`code/sn_scale_v12.py` → `run-distance-scale-v12/sn_scale_v12.json`; `code/distance_scale_v12.py` →
`run-distance-scale-v12/distance_scale_v12.json` and (with `--adopt 0.95`) `adopted_v12.json`;
`code/kids_level_v12.py` → `run-distance-scale-v12/kids_level_v12.json`.

The static distance law, 1 + z = e^(αD), has one scale. Rounds 10–11 used α0 = 2.4890 × 10⁻⁴ Mpc⁻¹
(H0-like 74.6, the project's calibration on 164 nearby groups). Round 11 read KiDS's common level as a
measurement of α: 13% lower would centre it (§21.5). Here α is fitted jointly, with the constants
refitted at each α. SPARC's 81 Hubble-flow galaxies (f_D = 1, D = cz/H0 with H0 = 73 in Lelli et al.
2016) now go into the static law too: D × 73/(αc), radii ∝ D, component speeds ∝ √D, g_obs ∝ 1/D. The
other 68 keep their Cepheid, TRGB, Ursa Major and supernova distances. `run.load_sparc(alpha)` does
this and reproduces the old loader exactly at alpha = None.

**What each data set prefers:**

| Data | best α/α0 | range | H0-like |
|---|---|---|---|
| Pantheon+ supernovae, Cepheid-calibrated, all 1,365 at z > 0.023 | **0.955** | ± 0.013 | 71.2 |
| the same, z = 0.023–0.15 only (the SH0ES Hubble-flow range) | 0.967 | ± 0.014 | 72.2 |
| the same, z = 0.1–0.3 only (the lenses' depth) | 0.944 | ± 0.014 | 70.5 |
| SPARC's Hubble-flow galaxies, by the fit statistic (400 bootstraps) | **0.87** | 0.82–0.92 | 64.9 |
| the same, by rms speed | 0.94 | 0.90–0.97 | 69.4 |
| X-COP | flat: rms 0.225 → 0.220 over ×0.80–1.05 (u 183 → 162 km/s) | | |
| KiDS lensing level | flat: +0.052 → +0.072 over ×0.80–1.05 | | |

* The supernova fit reuses the project's Pantheon+ reduction (brightness-distance-consistency): the 77
  Cepheid-host rows calibrate M with the static law's two factors of S = e^(αD) at their geometric
  distances, and the far rows are scored by generalized least squares on the full STAT+SYS covariance.
  At α0 it reproduces the earlier score for z = 0.1–0.3 (χ² 468.0).
* The supernovae's best α drifts with depth, 0.967 → 0.944. That drift is the law's shape against the
  supernovae, not its scale. The bounded beam-area term found earlier (flux × 1/(1 + ηf),
  f = 1 − 1/(1 + z)) removes it: fitted jointly, ×0.979 at every depth with η = 0.44–0.48, and χ² 1210
  against 1311 for all 1,365.
* TRGB (Freedman 2021) and JAGB (Lee et al. 2024) calibrations sit 0.097 and 0.160 mag fainter than the
  Cepheids. They would lower α by ×0.956 and ×0.928.

Combined by inverse variance (supernovae and SPARC's fit statistic): ×0.949 ± 0.013. **Adopted: ×0.95,
α = 2.3645 × 10⁻⁴ Mpc⁻¹ (H0-like 70.9).** Refitted there (a and g_d on SPARC, u on X-COP, alternated to
convergence in four rounds):

| | a (m/s²) | g_d (m/s²) | λ | u (km/s) | reach | SPARC | X-COP rms |
|---|---|---|---|---|---|---|---|
| round 11 (α0; SPARC at its published distances) | 6.547 × 10⁻¹¹ | 2.107 × 10⁻¹⁰ | 3.218 | 162.6 | 2.16 Mpc | 15.94 km/s | 0.221 |
| α0 with SPARC in the static law | 6.743 × 10⁻¹¹ | 2.212 × 10⁻¹⁰ | 3.28 | 165.6 | 2.20 Mpc | 15.98 | 0.221 |
| **round 12 (×0.95, SPARC in the static law)** | **6.298 × 10⁻¹¹** | **2.027 × 10⁻¹⁰** | **3.219** | **169.4** | **2.25 Mpc** | **15.87** | **0.222** |

Round 11 had left SPARC's distances at H0 = 73, 2% off its own scale. Every data set is now in the same
distances.

**Correction of round 11 (§21.5).** Lowering α by 13% centres KiDS only if the constants are held
fixed. Refitted, a follows the Hubble-flow galaxies' distances, and the level moves by only 0.02 dex
over ×0.80–1.05. The level does not measure α.

**What the level does measure** (`kids_level_v12.py`, the round-12 law; median log10(observed/predicted)):

| one change at a time | all | blue | red | disks | bulges | GAMA |
|---|---|---|---|---|---|---|
| round 12 as adopted | +0.065 | +0.011 | +0.078 | +0.066 | +0.034 | +0.009 |
| the other geometry (D_A = D/(1 + z)) | +0.058 | −0.002 | +0.074 | +0.059 | +0.030 | +0.003 |
| brightness distance 5% longer | +0.042 | −0.013 | +0.056 | +0.044 | +0.012 | −0.013 |
| angular-size distance 5% longer | +0.066 | +0.014 | +0.079 | +0.068 | +0.034 | +0.011 |
| star masses +0.13 dex | −0.005 | −0.062 | +0.010 | −0.003 | −0.034 | −0.060 |
| circumgalactic gas 0.5 × the stars, within 100 kpc | +0.026 | −0.020 | +0.052 | +0.007 | −0.003 | −0.032 |
| circumgalactic gas 1 × the stars (Brouwer et al.'s nominal) | −0.002 | −0.053 | +0.018 | −0.031 | −0.021 | −0.068 |

* Geometry and angular sizes hardly matter. In the lensing regime the level follows the lenses' mass at
  fixed light, and through it their brightness distance. The supernovae fix that distance at the
  lenses' depth to ±1.4%, and at ×0.95 the law matches them there.
* So the level is mass the lenses have that KiDS's baryon count (stars and cold gas) leaves out, or a
  bias in the star masses. A circumgalactic gas halo of 0.5–1 times the stars closes it; that is the
  amount Brouwer et al. 2021 take as their nominal estimate (isothermal, within 100 kpc). So would star
  masses 0.13 dex higher, inside their ±0.2 dex systematic range.
* Hot X-ray haloes are seen mainly around early types, so gas would lower the red lenses more than the
  blue ones, which is the direction the samples need (red +0.078, blue +0.011).

**Adopted: the round-12 law** (`regression/law_config.py`: `load_law()` now returns 'round12'; 'round11'
keeps the previous constants). Two settings travel with the law: `alpha_per_Mpc`, which
`common.apply_distances` writes into `collisions_v10.ALPHA` before any test runs (X-COP, KiDS, Mistele,
the Bullet, the far collisions), and `sparc_distances = 'static'`. The strong lenses' stored static
distances (at α0) are rescaled by α0/α and their star masses by (α0/α)² (`lenses_t35.patched_readers`).

**The suite** (full tier, baseline saved from this run): **59 pass, 11 close, 7 fail** (step 2: 59, 10, 8).
* No grade got worse. Improved: KiDS disks (fail → close, +0.077 → +0.066).
* Better without a change of grade (14): SPARC's median residuals (0.030 → 0.027), KiDS all, blue,
  bulges, GAMA and the Sérsic gap, SLACS (−0.032 → −0.028 dex), the Milky Way's vertical pull and mass
  inside 200 kpc, Fornax, Abell 520 inside 710 kpc (4.92 → 5.27 × 10¹⁴ against 5.84 ± 0.64) and El Gordo
  inside 500 kpc and 1 Mpc (19.8 → 21.4 × 10¹⁴ against 24.3 ± 12%).
* Worse without a change of grade (7): six Milky Way numbers, because a is 4% lower (the Sun's speed
  210.7 → 209.2 km/s; mass inside 50 kpc 3.63 → 3.56 × 10¹¹), and the KiDS colour gap (0.132 → 0.126
  against 0.153).
* The locked wide-binary forecast (round 10) at the new constants: γ = 1.032 at 7,000 AU and 1.076 at
  20,000 AU, against the support windows 1.03–1.05 and 1.08–1.10 "within the analysis errors" and the
  refutation limits 1.02 and 1.2. It stands as locked. Cassini's Q2: 4.4 × 10⁻²⁷ s⁻².
* Other moves: the Bullet's smaller half 1.29 → 1.38 × 10¹⁴ (measured 2.47–2.85 in these distances;
  still failing), MACS J0025's NW peak 72 → 79 kpc from its galaxies (pass), the collision stack's β
  0.016 → 0.019.

**Next for the level:** weigh the circumgalactic gas of KiDS-like galaxies from X-ray stacks (early types)
and UV absorption (late types) and predict the level from it; check the gas's effect on SPARC's outermost
points (an isothermal halo of M* inside 100 kpc puts 0.3 M* inside 30 kpc, a flatter β-model less).

### 22.4 A dynamical toy of the companion: which local rule gives one stream and no whirlpools?

`code/companion_toy_v12.py` → `run-companion-toy-v12/companion_toy_v12.json` and (with `--convergence`)
`convergence_v12.json`.

Round 11 derived the ordered rule from three requirements on the companion's steady energy flow J, fed
at ℓρ (§21.1): no energy lost (div J = ℓρ), one stream (energy density |J|/u, so f = |J|/(nu) = 1), and
no whirlpools (curl J = 0). The first and third fix J = (ℓ/4πG)(−g_N) uniquely. The third was assumed.
This toy builds the companion from local rules for its quanta on a 2D grid ([−1, 1]², 128², nothing
entering at the edges) and grades each rule on all three. Rules:
* **free:** quanta fly straight from where they are made, in every direction (waves passing through each
  other);
* **scatter:** the same, plus isotropic scattering (σ = 10 per unit length: diffusion);
* **annihilate:** the same, plus head-on pairs annihilating (counter-streams removed);
* **align:** the same, plus each quantum turning towards the local mean flow (BGK relaxation to a von Mises
  distribution of width 9° about J, energy kept exactly in every cell: a Vicsek-type rule);
* **guided:** quanta move along the local Newtonian field line, away from the matter, at u (one stream
  advected along ĝ; the way Alfvén waves are guided along magnetic field lines).

The first four are solved as a kinetic equation on 64 directions (explicit upwind, to a steady state); the
guided rule as one density advected along e = J_N/|J_N|. Grades, inside |x|, |y| < 0.8, weighted by
|J_N|: the power leaving the box over the power fed in; f; the rms log ratio |J|/|J_N| and the mean
angle between them; and "whirl", the share of J's spatial variation that is rotation, |curl J| /
(|∂xJx| + |∂yJx| + |∂xJy| + |∂yJy|). The Newtonian pattern itself scores whirl 0.001 on this grid.

| sources | rule | power out ÷ in | one stream f | f < 0.9 on | flux vs J_N (rms dex) | direction (deg) | whirl |
|---|---|---|---|---|---|---|---|
| uniform disk | (J_N itself) | | 1 | | 0 | 0 | 0.001 |
| uniform disk | free | 1.000 | 0.873 | 19% | 0.008 | 0.4 | 0.022 |
| uniform disk | scatter | 1.000 | 0.123 | 100% | 0.055 | 2.1 | 0.008 |
| uniform disk | annihilate | 0.760 | 0.863 | 20% | 0.135 | 0.4 | 0.022 |
| uniform disk | align | 1.000 | 0.984 | 0% | 0.350 | 7.6 | 0.428 |
| uniform disk | **guided** | 1.000 | 1.000 | 0% | 0.008 | 0.0 | 0.012 |
| two equal blobs | (J_N itself) | | 1 | | 0 | 0 | 0.001 |
| two equal blobs | free | 1.000 | 0.842 | 36% | 0.015 | 0.9 | 0.048 |
| two equal blobs | scatter | 1.000 | 0.124 | 100% | 0.058 | 2.4 | 0.015 |
| two equal blobs | annihilate | 0.805 | 0.844 | 37% | 0.099 | 2.2 | 0.070 |
| two equal blobs | align | 1.000 | 0.983 | 1% | 0.689 | 21.2 | 0.579 |
| two equal blobs | **guided** | 1.000 | 1.000 | 0% | 0.020 | 0.0 | 0.063 |
| unequal pair 4:1 | (J_N itself) | | 1 | | 0 | 0 | 0.001 |
| unequal pair 4:1 | free | 1.000 | 0.873 | 27% | 0.016 | 1.0 | 0.049 |
| unequal pair 4:1 | scatter | 1.000 | 0.126 | 100% | 0.058 | 2.6 | 0.014 |
| unequal pair 4:1 | annihilate | 0.791 | 0.869 | 29% | 0.111 | 2.5 | 0.053 |
| unequal pair 4:1 | align | 1.000 | 0.983 | 1% | 0.749 | 25.2 | 0.520 |
| unequal pair 4:1 | **guided** | 1.000 | 1.000 | 0% | 0.014 | 0.0 | 0.026 |
| three blobs | (J_N itself) | | 1 | | 0 | 0 | 0.001 |
| three blobs | free | 1.000 | 0.828 | 38% | 0.014 | 1.2 | 0.048 |
| three blobs | scatter | 1.000 | 0.123 | 100% | 0.057 | 2.3 | 0.015 |
| three blobs | annihilate | 0.800 | 0.828 | 40% | 0.107 | 3.6 | 0.093 |
| three blobs | align | 1.000 | 0.984 | 1% | 0.721 | 32.4 | 0.547 |
| three blobs | **guided** | 1.000 | 1.000 | 0% | 0.018 | 0.0 | 0.051 |

**Resolution** (two equal blobs; N = 64 / 128 / 256): guided: flux vs J_N 0.038 / 0.020 / 0.011 dex, whirl 0.107 / 0.063 / 0.040 (first-order grid error); free: whirl 0.079 / 0.048 (N = 64 / 128); align: whirl 0.594 / 0.579 and 26° / 21° off (steady: the rule's own whirlpools).

**Reading.**
* **Free streaming** keeps every watt and carries exactly Newton's pattern (in 2D, like 3D, the ballistic
  flux from isotropic sources is Gauss's field), but where streams from different places cross, the energy
  is not one stream: f = 0.83–0.87, below 0.9 on 19–38% of the region. This is round 10's
  rule 1, which dilutes the pull.
* **Scattering** keeps every watt and its flux is curl-free (a gradient), close to Newton's pattern inside
  the box, but it is the opposite of one stream: f ≈ 0.12.
* **Head-on annihilation** loses 19–24% of the energy and still leaves oblique crossings (f 0.83–0.87).
* **Alignment** makes one stream (f 0.98) and keeps every watt, but the flow does not follow Newton's
  pattern: 0.35–0.75 dex and 8–32° off, with whirlpools that form on their own (whirl 0.43–0.58, steady
  with resolution), as flocks do in the Vicsek model. Merged streams carry energy where Newton's field is
  weak, the "one stream of everything" that round 11 excluded.
* **Guided streaming** meets all three: f = 1, every watt kept, and Newton's pattern to the grid's
  precision (0.008–0.02 dex; its small whirl reading falls as the grid is refined).

**Why the guided rule works, exactly.** Take a thin tube of Newtonian field lines. Along it, Gauss's law
says the field's flux through the tube, |g|A, grows by 4πG times the mass the tube passes through. If the
companion flows along the same lines, its flux through the tube, |J|A, grows by ℓ times the same mass,
and no companion crosses the tube's walls. Both start at zero where the tube starts (the tube leaves a
potential minimum with zero cross-section). So |J| = (ℓ/4πG)|g| along every line: J = (ℓ/4πG)(−g_N)
everywhere, with no whirlpools and no energy lost, and one stream by construction. The curl-free condition
is not needed as an input; it follows from guiding. What remains postulated is the guiding itself: the
companion is a wave whose energy runs along gravity's own field lines, as an Alfvén wave's energy runs
along the magnetic field whatever its wave vector. Two consequences to keep:
* **Heat.** A randomly moving source's contribution has no coherent direction and is not guided (round
  10's scrambled phases); it adds as a scalar, which the clusters demand (h = 1, §21.1). The guided stream
  is the ordered part; the scalar sum S is the scrambled part.
* **Memory.** An Alfvén wave is carried along by the plasma it lives in. The analogue here, a guided
  stream carried by the matter whose field guides it, would give the memory the collisions need (§3.10 of
  the blog). The next toy moves the sources.

### 22.5 Where round 12 leaves things, and next

* **Settled this round:** MACS J0025's NW peak (dated by its own clocks, §22.1); the KiDS early/late gap
  (the lenses' heat measured, §22.2); the distance scale (fitted jointly, ×0.95, with every data set in
  the static law, §22.3); the no-whirlpool condition of §21.1 (supplied by guiding, §22.4).
* **Corrected:** round 11's readings of galaxy lensing as a measurement of α and of u (§21.5).
* **Suite** (round-12 law, full tier): **59 pass, 11 close, 7 fail** (round 11: 57, 9, 10). The seven:
  five faint dwarfs, the Bullet's smaller half, KiDS's red lenses.
* **Next:**
  1. **The guided companion with moving sources:** does a companion carried with its matter keep the
     memory the collisions use (§3.10 of the blog), and what does the Bullet's smaller half get?
  2. **The lenses' circumgalactic gas:** weigh it (X-ray stacks for early types, UV absorption for late
     types), predict KiDS's level from it, and check the gas's effect on SPARC's outermost points.
  3. **The distance law's shape:** the supernovae's best α drifts with depth (0.967 → 0.944); the
     bounded beam-area term (η ≈ 0.45) removes the drift. Derive it, or find what else does.
  4. **MACS J0025's star masses** from infrared light; its galaxy speeds prefer the published ones.
  5. From round 11's list, still open: the faint dwarfs, the Sun's speed, Abell 1689.

## 23. Round 13, 24 September 2026: the pull and the scrambling in one experiment

The request: "Connect the local-force model to the coherence model in one experiment. Use the same matter
distribution and the same microscopic interaction rules. Change only whether the sources move orderly, move
randomly and freely, or undergo frequent direction-changing collisions. Let test bodies respond to the local wave
interaction—not to the finished gravity formula. The decisive outcome would be that freely moving random sources
generate a stronger time-averaged attraction, while collisions suppress that enhancement. In the heat-dominated
regime of the proposed formula, doubling the velocity dispersion should approximately double the extra pull; that
scaling should emerge rather than be programmed into the force." It is proposal 3 of §9.1 of the blog.

Mid-round, after part 1's result: "if our exploration of derivation works keep enhancing it, if not look at another
solution to explore, be optimistic that we can pull these together."

**Short answer.**
* **Part 1, as asked** (§23.1–23.6). The experiment is built and validated, and the coherence half works in it:
  collisions protect the cold state by the amount round 3's Dicke formula gives. With round 10's rule the decisive
  outcome does not appear: free random motion *halves* the pull on the test bodies, and collisions protect it. The
  reason is exact. The pull on a locked emitter is the power it feeds the passing wave divided by the wave's speed,
  and in a locked state all identical emitters, test bodies included, sit at the same phase to their own local wave.
  So a test body is pulled exactly when the cloud's sources feed one another, which makes the cloud a bright chorus,
  and random motion can only break a chorus up (checked at 16 offsets: no exception). With the opposite offset the
  cloud's *output* follows the heat term's pattern (dark when cold, released by free motion, kept dark by
  collisions), but its test bodies are pushed.
* **Part 2, another solution** (§23.7). One rule escapes: an emitter runs a quarter cycle *ahead* of a weak wave
  (feeding it) and a quarter cycle *behind* a strong one (absorbing it), sources and test bodies alike. Inside a
  dense cloud the waves are strong, so a cold cloud goes dark; far away they are weak, so test bodies feed and are
  pulled. Then in the densest cloud (87 times darker than
  independent emitters at rest) free random motion raises the pull on the test bodies ×1.32, 1.60, 1.83 and 1.99 as
  σ doubles from q = 1/32 to 1/4; frequent collisions raise it only ×1.07–1.39, and orderly rotation not at all
  (×0.90–1.01). The first doubling multiplies the extra pull by 1.9 and the released power by 3.3; later ones less,
  as the test bodies begin to lose step. The wave at the test bodies grows exactly as √(cold + released power): the
  law's √(|g_N| + S), unprogrammed.
* **The data** (§23.8): the heat weight's exponent is pinned by galaxy lensing, not by the clusters: doubling σ
  multiplies the extra pull by 1.8–2.0, as the law has it.
* The law is unchanged by this round (the suite still reads 59 pass, 11 close, 7 fail). What changes is the
  mechanism behind its heat term: the naive one is excluded, and a candidate that works in the experiment, at speeds
  below a quarter of the re-timing rate, takes its place.

### 23.1 The experiment

`code/coherent_force_v13.py` → `run-coherent-force-v13/coherent_force_v13.json` (135 runs, 870 s on 4 cores).

* **Matter:** N = 100 point emitters of the companion wave, uniform in a ball of radius Rb = 1 companion
  wavelength (24 per λ³). Each has a fixed strength, its own pitch (natural frequencies spread uniformly over
  ±w, w = 0.1Γ) and a phase that locks, at rate Γ = 0.1 (per unit time; the wave period is 1), to the phase of
  the wave around it plus a fixed offset δ:
  ```
  dθ_j/dt = dω_j + Γ sin(arg E_j + δ − θ_j),   E_j = Σ_{l≠j} G(R_jl) e^{iθ_l},   G(R) = e^{i k0 R} / (4πR)
  ```
  δ = −π/2 is round 10's rule, "a quarter cycle ahead": the offset that feeds the passing wave and is pulled.
* **Test bodies:** 32 emitters of the same kind, with the same rule and rate, held on a sphere of radius 6λ
  around the cloud. Each feels only the local wave: F = ½ Re(e^{−iθ_p} ∇E(x_p)), the time average of q∇ψ.
  Nothing of the law's formula enters.
* **Motion, the only thing that changes:** at rest; **ordered** (rigid rotation, line-of-sight speed spread σ
  in the equator); **free** (Maxwellian, 1D spread σ, bouncing off the ball's surface); **colliding** (the same
  Maxwellian, each velocity redrawn at rate ν = 20 k0σ, so a source changes direction about 80 times while crossing one
  wavelength). The speed is set by q = k0σ/Γ, the Doppler drift rate over the locking rate:
  q = 0.3, 1, 3, 10.
* **Reading q in the law's terms.** Round 3's toy gives a free source's phase variance as (k0σ/Γ)² = q², and
  equating it with the heat weight, k = 3σ²/u², gives u = √3 Γ/k0 (round 3's mapping). So **q ≈ √k**:
  the four speeds are heat weights k ≈ 0.09, 1, 9 and 100, from a spiral's bulge (k ≈ 0.1–0.5) through the
  KiDS red lenses (2.0) to a rich cluster (σ ≈ 1,000 km/s, k ≈ 105). The heat-dominated regime is q ≳ 1.
* **Measured:** the mean pull on the test bodies toward the cloud (also split into equator and poles); the
  cloud's radiated power relative to independent emitters, Σ_jl Re(a_j a_l*) sinc(k0R_jl)/N; the locking order
  ⟨cos(θ − target)⟩. Time averages over 550 wave periods after a 350-period start, three seeds each.
* **Validation:** one source and one test body locked a quarter cycle ahead: the pull equals (k0/2)|E| (round 10)
  to five digits at 0.3, 1 and 3 wavelengths (0.83333, 0.25000, 0.08333).
* **Run sets:** A, round 10's rule (δ = −π/2); B, the same with a pitch spread of ±10Γ, so the sources cannot
  lock (a "lamp"); D, the opposite offset, δ = +π/2 (the passing wave is absorbed); E, a cloud twice as big
  (N = 200, same density); C, the offset over half a turn, at rest.

### 23.2 Results

**Round 10's rule (A).** At rest the sources lock into a chorus that radiates 4.39 times what independent
emitters would, and the test bodies are pulled 1.99 times as hard as by independent emitters. Pull relative to
rest (radiated power relative to independent emitters in brackets):

| q = k0σ/Γ (k ≈ q²) | ordered | free | colliding |
|---|---|---|---|
| 0.3 (0.09) | 0.95 ± 0.04 (4.39) | 0.84 ± 0.02 (3.74) | 0.96 ± 0.01 (4.24) |
| 1 (1) | 0.64 ± 0.02 (4.35) | **0.49 ± 0.02** (2.35) | 0.94 ± 0.01 (4.06) |
| 3 (9) | 0.27 ± 0.04 (4.39) | **0.36 ± 0.05** (1.76) | 0.87 ± 0.02 (3.67) |
| 10 (100) | 0.25 ± 0.03 (4.39) | **0.47 ± 0.01** (1.92) | 0.64 ± 0.01 (2.63) |

* **Free random motion weakens the pull,** to about half, by breaking up the chorus (4.39 → 1.8–2.4).
* **Collisions protect it** (0.87–0.96 up to q = 3), because they protect the chorus. That is Dicke narrowing
  in the same experiment, with the same rules.
* **Orderly rotation leaves the chorus intact** (4.35–4.39) but its wave pattern sweeps past the test bodies:
  at q = 3 those at the equator keep 0.20 of their pull and those near the poles, who see the pattern turn about
  their line of sight, 0.43. Random motion is isotropic (0.35 and 0.37).
* **Twice the cloud (E, N = 200 at the same density):** at rest 1.91 times the independent pull (radiated 4.83);
  at q = 1, free 0.52 of rest (2.50), colliding 1.08 (4.84). The same as N = 100.

**The opposite offset (D, δ = +π/2).** At rest the cloud is *dark*: it radiates 0.193 of what independent
emitters would, and its test bodies are *pushed* (−0.27 of the independent pull). Radiated power (push relative to
rest in brackets):

| q (k ≈ q²) | ordered | free | colliding |
|---|---|---|---|
| 0.1 (0.01) | 0.182 (0.97) | 0.189 (1.03) | 0.173 (0.94) |
| 0.3 (0.09) | 0.189 (0.98) | 0.217 (0.96) | 0.187 (1.01) |
| 1 (1) | 0.190 (0.72) | **0.425** (0.93) | 0.206 (0.98) |
| 3 (9) | 0.192 (0.33) | **0.669** (0.56) | 0.266 (0.89) |
| 10 (100) | 0.187 (0.15) | **0.767** (0.39) | 0.442 (0.80) |

* The cloud's output does exactly what the heat term needs: dark when cold, **released by free random motion**
  (four times brighter at q = 10), **kept dark by collisions** (0.21 against 0.43 at q = 1; 0.27 against 0.67 at
  q = 3), and **untouched by orderly motion** (0.18–0.19).
* The released power grows as σ² while it is small: +0.024 at q = 0.3 and +0.232 at q = 1, ×9.7 for σ ×3.3
  (σ² would give ×11). §23.4 doubles σ step by step.
* The collisions' suppression matches round 3's formula: colliding phase variance k0²σ²/(νΓ) against free
  k0²σ²/Γ², a ratio Γ/ν = 1/(20q); at q = 1 the released power is 0.013 against 0.232 (ratio 0.056; formula 0.05).
* But the test bodies are pushed, and the push *shrinks* as the cloud brightens (0.39 of rest at q = 10 with four
  times the power): the released wave changes faster than they can follow (§23.5).

**Sources that cannot lock (B, pitch spread ±10Γ).** Radiated power 1.14 at rest and 1.12–1.17 with any motion;
the test bodies get only **0.034** of the pull a body following the wave perfectly would get, and motion changes it
by factors of 1.0–2.3 with large scatter (±0.1–0.6).

**The offset at rest (C):**

| offset δ | −π/2 | −π/4 | 0 | +π/4 | +π/2 | π |
|---|---|---|---|---|---|---|
| radiated ÷ independent | 4.70 | 4.66 | 3.48 | 1.26 | **0.198** | 3.89 |
| pull ÷ independent pull | +2.02 | +2.02 | +1.54 | +0.26 | **−0.28** | +1.67 |

The pull's sign follows the cloud's brightness, not the offset by itself: the only dark cloud pushes. §23.3 explains
why, and scans the offset over a full turn.

### 23.3 Why random motion weakens the pull here: two exact statements

**1. The pull is the recoil of the power a body feeds the wave.** For an emitter of strength a in a field whose
local complex amplitude is E, the time-averaged power it gives the field is P = (ω/2) Im(a* E), and the
time-averaged force is F = ½ Re(a* ∇E). Where the wave is locally travelling, ∇E = i k E, so
```
F = −(P / c) k̂
```
exactly: a body that feeds the passing wave (P > 0) is pushed backwards along it, toward the wave's source,
by the momentum the extra wave carries forward. That is round 10's pull, now as momentum bookkeeping. A body
that absorbs (P < 0) is pushed away. (A spherical wave adds a small term along ∇|E|, ½ cos(θ − arg E)|a|∇|E|,
smaller by 1/(k0 r) = 2.7% at the test bodies.) **Attraction requires feeding.**

**2. In a locked state, all identical emitters sit at the same phase to their own wave.** Locked at a common
frequency Ω, each emitter obeys Γ sin(arg E_j + δ − θ_j) = Ω − dω_j; with equal pitches the stable solution is
the same for all: arg E_j − θ_j = arcsin(Ω/Γ) − δ. Test bodies locked to the cloud obey the same equation. So
test bodies and sources all feed, or all absorb, together. And the power the sources feed one another is the
mutual part of the cloud's radiated power: Σ_j Im(a_j* E_j) = (k0/4π) Σ_{j≠l} sinc(k0R_jl) Re(a_j* a_l). So:
```
test bodies pulled  ⇔  sources feed one another  ⇔  cloud brighter than independent emitters
```
Random motion can only disturb a chorus; so wherever the test bodies are pulled, random motion weakens the pull,
and wherever it brightens the cloud, the test bodies are being pushed.

**Checked over a full turn** (`code/joint_checks_v13.py`, part "offsets": 16 offsets at rest, two seeds each; "feeding"
is the mean sin(arg E − θ), positive when an emitter feeds its local wave):

| offset δ/π | −1 | −7/8 | −3/4 | −5/8 | −1/2 | −3/8 | −1/4 | −1/8 | 0 | 1/8 | 1/4 | **3/8** | **1/2** | **5/8** | 3/4 | 7/8 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| radiated ÷ independent | 3.87 | 4.32 | 4.31 | 4.27 | 4.36 | 4.22 | 4.26 | 3.80 | 3.31 | 1.99 | 1.30 | **0.53** | **0.20** | **0.40** | 1.39 | 2.27 |
| sources feeding | +0.75 | +0.88 | +0.88 | +0.87 | +0.90 | +0.88 | +0.87 | +0.85 | +0.74 | +0.37 | +0.11 | **−0.35** | **−0.63** | **−0.40** | +0.11 | +0.39 |
| test bodies feeding | +0.81 | +0.90 | +0.90 | +0.89 | +0.90 | +0.88 | +0.88 | +0.86 | +0.75 | +0.43 | +0.19 | **−0.23** | **−0.56** | **−0.24** | +0.30 | +0.48 |
| pull ÷ independent pull | +1.63 | +1.86 | +1.85 | +1.87 | +1.95 | +1.90 | +1.90 | +1.75 | +1.46 | +0.68 | +0.27 | **−0.15** | **−0.26** | **−0.15** | +0.37 | +0.74 |

At every offset the test bodies feed when the sources feed, and are pulled exactly when the cloud outshines
independent emitters. The three dark clouds all push.

### 23.4 Doubling the speed spread

`code/joint_checks_v13.py` → `run-coherent-force-v13/joint_checks_v13.json` (part "doubling"): the absorbing offset,
σ doubled from q = 1/8 to 8, free and colliding (ν = 20 k0σ, as in the main run), two seeds, in the main cloud and in
one twice as small (N = 100 in Rb = 0.5: 191 per λ³). Power released by the motion, radiated − radiated at rest, in
units of independent emitters' power:

| q | 1/8 | 1/4 | 1/2 | 1 | 2 | 4 | 8 |
|---|---|---|---|---|---|---|---|
| main cloud (at rest 0.197), free | 0.00 | 0.00 | 0.076 | 0.222 | 0.387 | 0.536 | 0.593 |
| main cloud, colliding | −0.01 | −0.01 | −0.01 | 0.013 | 0.038 | 0.108 | 0.201 |
| dense cloud (at rest **0.028**), free | 0.028 | 0.059 | 0.104 | 0.182 | 0.266 | 0.384 | 0.427 |
| dense cloud, colliding | 0.007 | 0.019 | 0.059 | 0.071 | 0.094 | 0.161 | 0.206 |

* **The denser cloud is 35 times darker than independent emitters at rest,** so the heat-dominated regime
  (released ≫ held) starts at a slower speed: at q = 1/4 the released power is already twice what the cloud emits at
  rest.
* **The release grows with σ, steeply at first and then more slowly:** per doubling, ×2.9, 1.7, 1.4 and 1.1 in the
  main cloud from q = 1/2; ×2.1, 1.75, 1.76, 1.46, 1.44 and 1.1 in the dense one. In the law's terms (power ∝ k ∝ σ^p)
  that is p ≈ 1.5 falling to 0.5, flatter than the data's 1.75–2 (§23.8). Independent emitters' level is the ceiling.
* **Collisions hold it back:** 3–17 times less release in the main cloud, 1.5–4 times in the dense one, where
  ν = 20 k0σ is only 2.5–10 Γ at the slower speeds.
* The test bodies are pushed throughout (−0.26 → −0.09 in the main cloud, −0.11 → −0.03 in the dense one).

### 23.5 The second obstacle: slow bodies cannot use a scrambled wave

A test body follows the wave's phase at its locking rate Γ. A wave from sources whose Doppler drift k0σ is faster
than Γ changes faster than it can follow. Set B shows the extreme: a lamp-like wave, with intensity at the test
bodies equal to independent emitters', gives them 3.4% of the pull a perfect follower would get. In D at q = 10
the cloud radiates four times more than at rest, yet its test bodies feel 0.39 of the push. In the law's terms,
Γ/k0 = u/√3 = 98 km/s: **the heat-dominated regime (σ ≳ 100 km/s) is exactly where the scrambled companion
changes faster than a body locking at the companion's own rate can follow.** So the heat term cannot reach test
bodies through slow phase-locking, whatever its sign; it needs a receiver that responds to the wave's energy
rather than its phase, or one that is itself broadband.

### 23.6 Emitters whose strength can change: a first look

A natural escape from §23.3 is to let each emitter's strength vary, as a laser's does: then a cloud could hold its
output back without every emitter giving up feeding. `code/joint_checks_v13.py`, part "amplitude": Stuart–Landau
oscillators, dA_j/dt = (μ + i dω_j − (1 + i c)|A_j|²)A_j + K e^{iβ}E_j, μ = 1, with β = −π/2 + arctan c so that a
lone emitter in a weak wave still locks a quarter cycle ahead (§23.3's reduction: θ − arg E = β − arctan c),
same cloud, at rest (one seed, 300 wave periods; the lone emitters' lock is what varies with c):

| shear c | K = 0.1 | K = 0.3 | K = 1 |
|---|---|---|---|
| 0 | 5.51 (+0.12) | 7.52 (+0.11) | 14.99 (−0.22) |
| −1 | 5.35 (+0.20) | 7.42 (−0.01) | 14.17 (−0.53) |
| −3 | 5.48 (+0.28) | 6.72 (−0.11) | 12.07 (−0.77) |
| +3 | 4.94 (+0.51) | 5.16 (+0.28) | 7.19 (+0.30) |

(radiated power relative to independent emitters of the same rule at rest; pull relative to independent emitters in
brackets.) **Brighter, not darker:** 5–15 times independent emitters, rising with the coupling, as a gain-type medium
should. The test bodies lock poorly (0.17–0.45), because the cloud's own strong field shifts its frequency further than
a weak passing wave can pull a lone emitter. Free amplitudes do not escape §23.3.

### 23.7 Part 2, another solution: an offset that depends on the strength of the wave

§23.3 needs one phase for all locked emitters. It fails if an emitter's offset depends on something that differs
between the inside of a body and the space around it. The simplest such thing is the strength of the local wave:
inside a dense cloud the neighbours' waves overlap and are strong; far away the wave is weak. The rule, for sources
and test bodies alike:
```
δ(|E|) = −π/2 + π x⁴ / (1 + x⁴),   x = |E| / E_s
```
a quarter cycle **ahead** of a weak wave (feeding it) and a quarter cycle **behind** a strong one (absorbing it).
**A physical reading:** an amplifier whose gain saturates while its losses do not (every laser medium has some fixed
loss). A weak wave sees net gain and is fed; a strong one saturates the gain below the loss and is absorbed. E_s is
the strength at which saturated gain equals the loss.

`code/strength_offset_v13.py` → `run-coherent-force-v13/strength_offset_v13.json` (182 runs). `simulate` in
`code/coherent_force_v13.py` takes the rule as an option; with it off, the stored runs of §23.2 reproduce exactly.

**At rest, over the switch strength E_s** (two seeds; radiated ÷ independent emitters / pull ÷ independent pull;
sources' and test bodies' feeding in brackets):

| cloud | E_s = 0.03 | **0.1** | 0.3 | 1 | 3 |
|---|---|---|---|---|---|
| Rb = 1 (24 per λ³) | 0.178 / −0.20 (−0.67, −0.35) | **0.169 / +0.19** (−0.67, +0.48) | 0.281 / +0.38 | 1.23 / +0.82 | 3.29 / +1.59 |
| Rb = 0.5 (191 per λ³) | 0.027 / −0.01 | **0.030 / +0.11** (−0.57, +0.59) | 0.125 / +0.19 | 0.62 / +0.48 | 2.67 / +1.14 |
| Rb = 0.35 (557 per λ³) | 0.011 / +0.03 | **0.011 / +0.06** (−0.42, +0.50) | 0.030 / +0.09 | 0.41 / +0.24 | 1.60 / +0.66 |

* Too small an E_s puts the test bodies in the strong regime too (Rb = 1 at E_s = 0.03: all absorb, dark and pushed,
  as in §23.3); too large a one puts the sources in the weak regime (all feed: a bright chorus that pulls).
* **In between, the cloud is dark and its test bodies are pulled:** E_s = 0.1–0.3 at every density, and up to 1 in the
  denser clouds. The sources absorb (feeding −0.39 to −0.67) while the test bodies feed (+0.48 to +0.66). §23.3's
  statement no longer binds, because the two no longer sit at one phase.

**Moving** (E_s = 0.1; four seeds; collisions redraw each velocity 20 times per locking time or per Doppler time,
whichever is shorter, ν = 20Γ max(1, q)). Pull on the test bodies ÷ the same cloud at rest (± s.e.); radiated power ÷
at rest in brackets:

| densest cloud (Rb = 0.35; at rest radiated 0.0115 = 1/87, pull 0.057) | q = 1/32 | 1/16 | 1/8 | 1/4 | 1/2 | 1 |
|---|---|---|---|---|---|---|
| **free** | **1.32** ± 0.05 (1.44) | **1.60** ± 0.07 (2.47) | **1.83** ± 0.06 (4.50) | **1.99** ± 0.05 (6.86) | 1.91 ± 0.09 (9.18) | 1.66 ± 0.10 (14.1) |
| colliding | 1.07 ± 0.04 (1.20) | 1.07 ± 0.06 (1.12) | 1.18 ± 0.04 (1.31) | 1.39 ± 0.03 (1.85) | 1.91 ± 0.07 (4.62) | 1.99 ± 0.06 (8.05) |
| rotating | 1.01 ± 0.06 (1.00) | 0.98 ± 0.06 (1.00) | 0.96 ± 0.05 (1.00) | 0.90 ± 0.04 (1.00) | 0.81 ± 0.04 (1.00) | 0.55 ± 0.04 (1.00) |

| dense cloud (Rb = 0.5; at rest 0.029 = 1/34, pull 0.108) | q = 1/32 | 1/16 | 1/8 | 1/4 | 1/2 | 1 |
|---|---|---|---|---|---|---|
| **free** | 0.97 ± 0.04 (0.92) | **1.10** ± 0.02 (1.13) | **1.26** ± 0.04 (2.05) | **1.35** ± 0.03 (3.09) | 1.25 ± 0.04 (4.67) | 1.04 ± 0.06 (7.15) |
| colliding | 0.93 ± 0.04 (0.91) | 0.90 ± 0.01 (0.98) | 0.95 ± 0.03 (0.92) | 0.94 ± 0.01 (1.06) | 1.23 ± 0.03 (1.88) | 1.32 ± 0.05 (3.44) |
| rotating | 1.01 ± 0.04 (1.00) | 0.99 ± 0.02 (1.00) | 0.98 ± 0.03 (1.00) | 0.95 ± 0.03 (1.01) | 0.83 ± 0.04 (0.98) | 0.65 ± 0.02 (0.99) |

* **The decisive outcome, at speeds up to a quarter of the re-timing rate:** free random motion strengthens the
  pull (to ×1.99 in the densest cloud, ×1.35 in the other); frequent collisions suppress the enhancement (×1.07–1.39
  and ×0.90–0.95), because they keep the cloud dark (released power ×1.1–1.9 against ×1.4–6.9 free); orderly
  rotation gives none (×0.90–1.01). The same rules, only the motion changed, and test bodies that feel only the wave.
* **The law's square root, unprogrammed.** The wave at the test bodies (mean |E|, densest cloud, free) grows ×1.20,
  1.56, 2.12, 2.62, 3.04, 3.76; the square root of the radiated power's growth is 1.20, 1.57, 2.12, 2.62, 3.03, 3.75.
  Cold and released power add, and the pull follows the square root of the sum: √(|g_N| + S) in the law.
* **Doubling.** From q = 1/32 to 1/16 the released power grows ×3.3 (σ^1.7), the wave's extra strength at the test
  bodies ×2.8 and **the extra pull ×1.9**, close to the law's ×4 and ×2 and inside the data's ×3.4–4 and ×1.8–2.0
  (§23.8). The next doublings give less: released ×2.4 and 1.7, extra pull ×1.4 and 1.2.
* **Where it stops.** Beyond q ≈ 1/4 the free clouds keep brightening (×9–14) but the pull falls back: the test
  bodies lose step (their feeding falls from +0.50 to +0.22), as in §23.5. The colliding clouds, whose wave changes
  more slowly, catch up at q ≥ 1/2. Rotation sweeps the pattern past the test bodies (×0.55 at q = 1).
* **Reading u again.** Here heat dominates (released > held) from q ≈ 1/32–1/16, far below the re-timing rate,
  because the cold cloud is so dark. Round 3's mapping, u = √3 Γ/k0, assumed the release comes from phase variance
  alone; with a dark reservoir, u is set by the locking rate and the darkness together, and the darker cold matter is,
  the more of the heat-dominated regime the test bodies can follow.

**Why the rule may be more than a device.**
* It is the same shape as the law's release factor (§3.4 of the blog): **where the companion is strong it is held
  back, where it is weak it is free.** The law applies that to the Newtonian field, exp(−|g_N|/g_d), and the ordered
  companion's intensity is |g_N| (§21.1), so g_d plays the part of E_s². Whether one saturation gives both the
  release factor and the heat term is the first thing to derive.
* It puts the heat term's energy where §23.9 says it must be: a cold, dense body sits in its own strong wave and
  absorbs, holding its supply back; random motion scrambles that wave, more emitters find themselves in weak waves
  and feed, and the supply is tapped.

### 23.8 The heat weight's exponent, from the data

The request's scaling, "doubling σ doubles the extra pull", is the law's heat weight k = 3σ²/u², since the extra
pull goes as √k where heat dominates. `code/heat_exponent_v13.py` → `run-coherent-force-v13/heat_exponent_v13.json`
asks the data what power of σ they want: k = 3(σ/u)^p, u refitted on X-COP's 12 clusters (round-12 law, static
distances) for each p, then the KiDS early/late lensing gap (the lenses' star speeds measured, §22.2) predicted with
that u. ⟨fσ^p⟩ is approximated as ⟨f⟩(⟨fσ²⟩/⟨f⟩)^{p/2} in each 0.1-dex bin.

| p | u (km/s) | X-COP rms (ln) | k red / blue lenses | gap, colour split (0.153 ± 0.04) | gap, Sérsic split (0.154) | KiDS level, all lenses |
|---|---|---|---|---|---|---|
| 1 | 32.8 | 0.2231 | 11.7 / 5.2 | 0.149 | 0.201 | −0.221 |
| 1.25 | 63.0 | 0.2221 | 7.5 / 2.8 | 0.161 | 0.216 | −0.132 |
| 1.5 | 97.7 | 0.2215 | 4.8 / 1.6 | 0.160 | 0.210 | −0.055 |
| 1.75 | 133.7 | 0.2213 | 3.1 / 0.9 | 0.148 | 0.188 | **+0.011** |
| **2 (the law)** | **169.4** | 0.2215 | 2.0 / 0.5 | 0.126 | 0.155 | +0.065 |
| 2.25 | 203.9 | 0.2221 | 1.3 / 0.3 | 0.101 | 0.121 | +0.107 |
| 2.5 | 236.6 | 0.2230 | 0.8 / 0.2 | 0.077 | 0.089 | +0.139 |
| 3 | 296.3 | 0.2260 | 0.3 / 0.06 | 0.039 | 0.043 | +0.179 |

* **The clusters alone cannot tell** (rms 0.221–0.226 for p = 1–3): their galaxies span too narrow a range of
  speeds (typically 800 km/s), and u absorbs the change.
* **The lenses, about six times slower (red lenses' stars ≈ 140 km/s against the clusters' galaxies' ≈ 800), pin it:** within 1σ of the colour-split gap, p = 1–2; of the Sérsic-split gap,
  p = 1.75–2.25; of both, **p = 1.75–2**. Doubling σ multiplies the extra pull by 2^(p/2) = **1.8–2.0**. The data
  demand the request's scaling, so whatever the mechanism, it must produce it.
* **A side lead:** at p = 1.75 the KiDS level, all lenses' common 16% excess over the law (+0.065 dex, §22.3),
  falls to +0.011 (3%) with the clusters fitted as well (u = 134 km/s). Whether the rest of the suite accepts
  p = 1.75 is untested (the heat weight enters the galaxies' bulges, every cluster and every collision).

### 23.9 Where the heat term's extra power must come from

If the heat term's extra power, k times ℓ per kilogram, came out of the sources' motion, it would drain their
kinetic energy per kilogram, (3/2)σ², in
```
t = (3/2)σ² / (k ℓ) = u² / (2ℓ) = u / a = 1.694 × 10⁵ m/s ÷ 6.30 × 10⁻¹¹ m/s² = 2.69 × 10¹⁵ s = 85 million years
```
(a = 2ℓ/u), whatever the temperature. Clusters and ellipticals are far older and still hot, so the extra power
comes from matter's own supply, the one that feeds the companion at ℓ in the first place. Motion only unlocks it.
That is what set D's dark cloud does: its emitters hold their output back when cold, and random motion releases
it. Collisions keep it held.

### 23.10 What this means, and next

* **Derived and confirmed in one experiment:** round 10's pull, now as momentum bookkeeping (F = −(P/c) k̂: the
  recoil of the power fed), and round 3's protection of coherence by collisions, with the pull measured and the
  suppression as round 3's formula gives. "Gas does not count as hot" survives the join.
* **Excluded:** the heat term's mechanism as the blog stated it ("scrambled contributions don't cancel, so they pull
  harder"), for any one fixed offset and for free amplitudes. Scrambling weakens a chorus, and a slow body cannot use
  a scrambled wave. The heat term itself stands on the data (clusters, the elliptical/spiral lensing gap, the
  collisions), and §23.8 shows the data demand its scaling.
* **Found:** one rule, the same for every emitter, that puts the heat term's pattern into the pull itself: an offset
  set by the wave's strength, ahead of weak waves and behind strong ones. Cold dense clouds go dark and still pull;
  free random motion strengthens the pull (to ×2), frequent collisions suppress the gain, orderly rotation gives
  none; the pull follows √(cold + released power), the law's √(|g_N| + S); the first doubling of σ gives ×1.9 in the
  extra pull. The extra energy comes from matter's own supply, as §23.9 requires.
* **Not yet:**
  * the gain holds while the motion is slower than about a quarter of the re-timing rate; beyond that the test bodies
    lose step (§23.5). The darker the cold state, the more of the heat-dominated regime lies below that limit;
  * the release flattens with σ (released power ×3.3, 2.4, 1.7 per doubling; the data want ×3.4–4 sustained over a
    factor of ten in σ, from the lenses to the clusters);
  * the rule is a hypothesis with a physical reading, not yet a derivation.
* **The test bench stays:** any candidate rule can be dropped into `simulate` in `code/coherent_force_v13.py` and
  graded by the same runs.
* **Next (round 14):**
  1. **Derive the strength-dependent offset** from a saturable emitter (gain that saturates, loss that does not), and
     test whether the same saturation gives the law's release factor, exp(−|g_N|/g_d).
  2. **Bodies that keep in step:** a pull drawn from the companion's energy flow (round 12's guided stream) rather
     than its phase, or a cold state dark enough that the whole heat-dominated regime lies at slow speeds; then u
     mapped onto the locking rate and the darkness together.
  3. **The release's steepness:** why the toy's release flattens (σ^1.7 → σ^0.8), and what keeps it near σ² over a
     factor of ten.
  4. **The heat exponent p = 1.75 on the full suite** (§23.8): it removes KiDS's common level.
  5. Round 12's list (§22.5): the guided companion with moving sources, the lenses' gas, the distance law's shape,
     MACS J0025's star masses, the faint dwarfs, the Sun's speed, Abell 1689.

## 24. Round 14, 24 September 2026: the next steps, one at a time

The request: "proceed with the next steps, keep the blog and main current": round 13's list (§23.10). In the
order done:
1. the heat weight's exponent on the full suite (§24.1);
2. where "feed quiet waves, absorb loud ones" comes from, and whether it is the release factor (§24.2);
3. how the release grows with the speed spread, and when test bodies lose step (§24.3).

### 24.1 The heat exponent on the full suite

Round 13 (§23.8) found that the KiDS early/late gap allows k = 3(σ/u)^p with p = 1.75–2 and that p = 1.75 removes
KiDS's common level with the clusters unchanged. Here the whole suite judges it.

* **One switch.** Every heat weight in the code now goes through `law.heat_weight` (from σ) or `law.k_from_sig2`
  (from a mean-square speed, with its prefactor: 3 for a 1D dispersion, 1 for a 3D mean square, 3 − 2β for an
  anisotropic radial one), with the exponent in `law.HEAT_P`. Nine modules were rewired (`bullet_v3`, `bullet_v4`,
  `bullet_main_v5`, `mw_model`, `milky_way_v7`, `mw_dwarfs_v7`, `kids_static_v11`, `lens_heat_sdss_v12`, `law`). At
  p = 2 the arithmetic is each script's own: the quick tier reproduces the baseline to the last digit (54 same, 0
  changed). For p ≠ 2 an effective 1D dispersion σ_eff² = (prefactor) σ²/3 enters as 3(σ_eff/u)^p; the lenses'
  bulge and disk parts are weighted separately, as in §23.8. The law's `heat_exponent` sets it
  (`regression/law_config.py`, `common.apply_distances`).
* **A slip in the suite, fixed.** Its refit of u (`common.refit_constants`) used the round-3 X-COP sample
  (projected stars, published distances) while the clusters test grades the static, deprojected one that the
  adopted constants were fitted on. For the adopted law it returned u = 192 km/s instead of 169.4. Laws based on
  round 11 or 12 now refit u on the graded sample; older candidates keep the old one, so they still reproduce.
  The control (`heat_p2_refit`: the adopted law through the same refit) now returns a = 6.298 × 10⁻¹¹ m/s² and
  u = 169.44 km/s and reproduces the baseline exactly (83 checks the same).
* **The result** (full tier; `regression/runs/heat_p175-full/`):

| law | a (m/s²) | u (km/s) | pass / close / fail |
|---|---|---|---|
| adopted (p = 2) | 6.298 × 10⁻¹¹ | 169.4 | 59 / 11 / 7 |
| p = 1.75, a and u refitted | 6.181 × 10⁻¹¹ | 132.2 | **60 / 10 / 7** |

| check | p = 2 | p = 1.75 |
|---|---|---|
| KiDS, all lenses: level (0 ± 0.025 dex) | 0.065 (close) | **0.012 (pass)** |
| KiDS, red lenses | 0.078 (fail) | **0.018 (pass)** |
| KiDS, disk-dominated lenses | 0.066 (close) | **0.038 (pass)** |
| KiDS early/late gaps (0.153, 0.154) | 0.126, 0.155 | 0.149, 0.189 |
| Mistele et al.: ellipticals' lensing speeds, 50–300 kpc (rms z) | 2.6 (close) | **5.2 (fail)** |
| Mistele et al.: spirals | 2.9 (close) | 2.6 (close) |
| SLACS: lensing minus kinematic star mass (0 ± 0.023 dex) | −0.028 (pass) | −0.054 (close) |
| SPARC's 25 bulge-dominated galaxies (MOND 30.35 km/s) | 29.4 (pass) | 30.5 (close) |
| X-COP rms | 0.2215 | 0.2213 |
| Abell 520: P2 / P6 | pass / close | close / pass |

* **Reading.** p = 1.75 trades one tension for another. Brouwer et al.'s KiDS relation wants a little more heat
  at the red lenses' ≈ 140 km/s; Mistele et al.'s circular speeds from the same survey, and SLACS's massive
  lenses (≈ 250 km/s), want a little less. Not adopted: the law keeps p = 2, and the exponent is now a switch for
  when the two KiDS analyses are reconciled. Round 13's toy release grew as σ^1.7 at its onset (§23.7), so the
  microscopic model does not choose between the two either.

### 24.2 Where "feed quiet waves, absorb loud ones" comes from

`code/singer_v14.py` → `run-coherent-force-v13/singer_v14.json` (297 runs); `code/release_shape_v14.py` →
`run-coherent-force-v13/release_shape_v14.json`. `simulate` in `code/coherent_force_v13.py` gains options for both
emitters below (and `stop_at`, to stop the motion mid-run); with them off, every stored run reproduces exactly.

**A. A singer that also absorbs.** Each emitter keeps round 10's locked oscillation (a quarter cycle ahead of its
wave) and gains a passive part driven by the wave, i c0 E_j, which takes energy in proportion to |E|²; the sources
are solved each step as s = (1 − i c0 G)⁻¹ a (c0 ≤ 4π/k0 keeps the passive part passive). A lone one in a wave of
strength E is pulled by (k/2)(|E| − c0|E|²): it feeds quiet waves and absorbs loud ones, with no switch put in.

| at rest (two seeds): radiated ÷ independent / pull ÷ independent pull | c0 = 0.25 | 0.5 | 1 | 2 |
|---|---|---|---|---|
| Rb = 1 | 1.92 / +1.25 | 1.11 / +0.90 | 0.47 / +0.56 | 0.28 / +0.23 |
| Rb = 0.5 | 1.66 / +1.16 | 0.59 / +0.70 | **0.17 / +0.37** | 0.22 / +0.05 |
| Rb = 0.35 | 1.13 / +1.01 | 0.34 / +0.54 | **0.021 / +0.075** | 0.27 / −0.01 |

* **Dark and pulled, derived.** The singers all feed (they keep round 10's timing), and the passive parts swallow the
  cold chorus: at c0 = 1 the densest cloud radiates 1/47 of what independent singers would, and its test bodies are
  pulled. §23.3's statement is escaped by a second channel, not by changing the singers' timing.
* **Moving** (densest cloud, c0 = 1, three seeds; pull ÷ rest, radiated ÷ rest in brackets): free ×1.14, 1.59, 1.83,
  3.62, 3.98, 3.83 (q = 1/32 … 1; radiated ×1.3–4.6); colliding ×1.08, 1.21, 1.07, 1.47, 2.16, 3.84; rotating ×1.00 →
  0.60. **But the gain is un-jamming, not release:** at rest the dense cloud's singers are only half locked (locking
  0.55, unchanged over 900 locking times), and free motion at q ≥ 1/4 lets them lock (0.90–0.99), which makes the
  cloud five times brighter. It has memory: moved at q = 1/2 and then stopped, 2 of 4 clouds stayed bright and locked
  (radiated 0.10–0.12, locking 0.997) and 2 went back (0.021–0.022).
* **With resonant passive parts** (line width γ_a = Γ; narrower ones cannot follow the singers, and the cold cloud is
  bright: 0.9 at γ_a = 0.03, 7 at 0.01), the cold cloud is dark and well locked (0.110; locking 0.93), and motion
  releases nothing: ×1.04–1.07 up to q = 1, then dimmer (free q = 8: ×0.73). A separate absorber keeps absorbing
  whatever reaches it, so motion cannot free what it holds.

**B. A singer with a power balance.** Give each emitter a fixed supply P_s (the law's ℓ per kilogram) and a loss to
the wave around it, c_L|E|² (the same absorption as in A). In a steady state it must feed its wave the difference,
|E| sin(lead) = P_s − c_L|E|², so its timing is set as a synchronous machine's load angle is: a generator, ahead of a
quiet wave; a motor, behind a loud one (the switch at |E| = √(P_s/c_L)). Round 13's rule, now from two physical
numbers: δ(|E|) = −arcsin(clip(P_s/|E| − c_L|E|, −1, 1)).

| at rest (two seeds): radiated ÷ independent / pull ÷ independent pull | P_s 0.05, E_s 0.1 | P_s 0.1, E_s 0.1 | P_s 0.05, E_s 0.2 | P_s 0.2, E_s 0.3 |
|---|---|---|---|---|
| Rb = 1 | 0.179 / +0.13 | 0.170 / +0.20 | 0.514 / +0.47 | 0.294 / +0.38 |
| Rb = 0.5 | 0.032 / +0.10 | 0.030 / +0.12 | 0.218 / +0.26 | 0.137 / +0.20 |
| Rb = 0.35 | 0.012 / +0.06 | **0.013 / +0.07** | 0.113 / +0.16 | 0.070 / +0.13 |

All twelve are dark and pulled (sources feeding −0.24 to −0.67, test bodies +0.35 to +0.64).

| moving (P_s 0.1, E_s 0.1; four seeds): pull ÷ rest (radiated ÷ rest) | q = 1/32 | 1/16 | 1/8 | 1/4 | 1/2 | 1 |
|---|---|---|---|---|---|---|
| densest (Rb 0.35), **free** | **1.30** (1.5) | **1.51** (2.2) | **1.81** (4.6) | 1.71 (6.0) | 1.65 (8.6) | 1.69 (13.4) |
| densest, colliding | 0.98 (1.1) | 1.01 (1.1) | 1.08 (1.2) | 1.31 (1.6) | 1.68 (4.4) | 1.75 (7.4) |
| densest, rotating | 0.98 | 0.95 | 0.93 | 0.88 | 0.80 | 0.52 |
| dense (Rb 0.5), free | 0.94 (1.0) | 1.05 (1.2) | 1.22 (2.1) | 1.20 (3.2) | 1.11 (4.4) | 0.91 (7.0) |
| dense, colliding | 0.96 | 0.95 | 0.92 | 0.96 (1.1) | 1.18 (1.8) | 1.18 (3.4) |

* **Round 13's pattern, from the power balance:** free random motion strengthens the pull (×1.3–1.8 in the densest
  cloud up to q = 1/8), frequent collisions hold the gain back (×0.98–1.08 there), rotation gives none. The first
  doublings multiply the released power by 2.5 and 3.0 and the extra pull by 1.7 and 1.6.
* **No memory:** moved at q = 1/4 and stopped, all four clouds return to their cold values (radiated 0.009–0.015,
  pull 0.050–0.071 of independent, against 0.0125 and 0.063 at rest). The gain is a steady response to motion.
* **Limits, as in round 13:** the pull stops growing beyond q ≈ 1/8 while the release goes on (×13 at q = 1): the
  test bodies lose step. The less dense cloud gains less (×1.2). §24.3 looks at both.

**The hold-back shape, on the galaxies.** A lone emitter of either kind is held back in loud waves: A by
1 − |E|/E_s, B by (P_s − c_L|E|²)/|E| (the full round-10 pull in quiet waves). The ordered companion's intensity
is |E|² ∝ |g_N| (§21.1), so each gives a release factor R(|g_N|) to set against the law's exp(−|g_N|/g_d). The 149
SPARC galaxies judge them (a and the scale refitted; u and the heat held at the round-12 law):

| R(g), x = √(g/g_s) | constants | halves the pull at g_N (m/s²) | typical miss (km/s) | fit statistic |
|---|---|---|---|---|
| exp(−g/g_d) (the law) | a, g_d | 1.41 × 10⁻¹⁰ | 15.87 | 0.21111 |
| A: 1 − x, stopped at 0 | a, g_s | 1.35 × 10⁻¹⁰ | **15.86** | 0.21169 |
| B: clip(A (1/x − x), 0, 1) | a, g_s, A | 1.23 × 10⁻¹⁰ | 16.08 | **0.21103** |
| A with its push (1 − x) | a, g_s | 2.64 × 10⁻¹⁰ | 16.19 | 0.21475 |
| an absorber that saturates, 1 − x/(1 + x²) | a, g_s | 6.95 × 10⁻¹⁰ | 16.78 | 0.21709 |

**The galaxies cannot tell the law's exponential from the absorbers' hold** (15.86 and 16.08 km/s against 15.87; the
fit statistic, which the fits minimise, is lowest for B). Holds that push in strong fields, or let the pull return in
them, do worse. So the release factor, "strong gravity holds the companion back", may be the same loss that sets the
emitters' timing. What stops the hold at zero rather than letting it push is the next thing to derive.

**Reading.** "Feed quiet waves, absorb loud ones" follows from a fixed supply and a loss that grows with the
wave's intensity; with it, the heat pattern of round 13 appears with no memory, and the same loss gives a hold-back
the galaxies accept in place of the exponential release factor. A passive absorber alone gives the darkness but not
the release: the heat term needs the cold hush to be arranged by the emitters' own timing.

### 24.3 How the release grows, and when test bodies lose step

`code/release_depth_v14.py` → `run-coherent-force-v13/release_depth_v14.json` (153 runs). The power-balance singer
(P_s = 0.1, c_L = 10) in three clouds of 100 of increasing density, free and colliding, q = 1/64 … 2, three seeds.
Power released by the motion over the power held at rest ("released/held"; the analogue of the heat weight k):

| free motion | at rest: radiated ÷ independent | q = 1/32 | 1/16 | 1/8 | 1/4 | 1/2 | 1 | 2 |
|---|---|---|---|---|---|---|---|---|
| Rb 0.5 (191 per λ³) | 1/33 | −0.05 | 0.08 | 1.07 | 2.0 | 3.3 | 5.8 | 9.8 |
| Rb 0.35 (557 per λ³) | 1/76 | 0.39 | 1.09 | 3.2 | 4.9 | 7.3 | 11.9 | 19.2 |
| Rb 0.25 (1,528 per λ³) | 1/107 | 0.78 | 3.3 | 5.1 | 7.7 | 9.3 | 14.3 | 23.8 |

| pull ÷ rest (test bodies' feeding) | at rest | 1/32 | 1/16 | 1/8 | 1/4 | 1/2 | 1 | 2 |
|---|---|---|---|---|---|---|---|---|
| Rb 0.5, free | (0.60) | 0.93 | 1.03 | **1.21** (0.51) | 1.17 | 1.10 | 0.89 (0.23) | 0.44 (0.10) |
| Rb 0.35, free | (0.52) | 1.22 | 1.43 | 1.71 (0.44) | 1.68 | 1.63 | **1.73** (0.25) | 0.82 (0.10) |
| Rb 0.25, free | (0.47) | 1.52 | **1.95** (0.45) | 1.80 | 1.68 | 1.68 | 1.39 (0.18) | 0.76 (0.08) |
| Rb 0.35, colliding | | 0.98 | 0.94 | 1.05 | 1.32 | 1.62 | 1.69 | 1.54 |
| Rb 0.25, colliding | | 0.93 | 1.05 | 1.24 | 1.65 | 1.80 | 1.55 | 1.46 |

* **The square-root law holds at every speed:** the wave at the test bodies grows as the square root of the
  radiated power in every row to 1–2% (e.g. ×4.47 against √20.2 = 4.49): cold and released power add.
* **Darkness sets where heat takes over.** Released exceeds held from q* ≈ 1/8, 1/16 and 1/24 in clouds 33, 76 and
  107 times darker than independent singers: q* ≈ 4 × (radiated ÷ independent at rest). Collisions delay the onset by
  a factor of 4–5 in speed (Rb 0.35: q ≈ 1/3 instead of 1/16; Rb 0.25: 1/6 instead of 1/24), as round 3's formula
  has them.
* **Test bodies keep step up to about the re-timing speed.** Their feeding falls from ≈ 0.5 at rest to ≈ 0.25 at
  q = 1, but the wave grows faster, so the pull stays raised (×1.4–1.7 in the two densest clouds at q = 1); at q = 2
  the feeding collapses (0.08–0.10) and the pull with it. The range of speeds over which heat dominates and bodies
  keep step, q* to 1, widens with darkness: ×8 (Rb 0.5), ×16 (0.35), ×24 (0.25). Round 3's reading, u = √3 Γ/k0,
  put the heat-dominated regime (k > 1) exactly where bodies lose step; with a dark reservoir the two separate, and
  the darker cold matter is, the wider the gap.
* **The steepness does not match yet.** After its onset (×3–6 per doubling, σ^1.6–2.6) the release grows by ×1.5–1.7
  per doubling (σ^0.6–0.8) over four doublings in every cloud. The law, and the lensing data (§23.8), want about ×4
  per doubling (σ²) sustained over a factor of ten in σ. In the toy the reservoir empties too gently: scrambled, most
  emitters still sit in strong waves and keep absorbing, and only a slowly growing share flips to feeding. What makes
  the share grow as σ² (a larger cloud, a loss that rises faster than |E|², emitters of several strengths) is the
  next question.

### 24.4 Where round 14 leaves things, and next

* **Settled this round:** the heat exponent on the full suite (p = 1.75 trades KiDS's level for the massive
  ellipticals' lensing; p = 2 kept); a slip in the suite's refit of u (fixed); the origin of round 13's timing rule
  (a power balance: fixed supply against a loss that grows with the wave's intensity), with no memory.
* **Found:** the singers' hold-back in loud waves fits the 149 galaxies as well as the law's exponential release
  factor (15.86 and 16.08 against 15.87 km/s), so the release factor and the heat mechanism may be one loss. Darker cold
  matter moves the start of the heat-dominated regime to lower speeds (q* ≈ 4 × darkness), and test bodies keep step
  up to about the re-timing speed, so a deep enough hold covers the whole heat range.
* **Not yet:** the release grows as σ^0.6–0.8 after its onset, against σ² in the law and the lensing data; why the
  hold stops at zero rather than pushing; where the absorbed energy goes (it must return to the emitters' supply, or
  cold matter would warm).
* **Suite** unchanged: 59 pass, 11 close, 7 fail (the law is unchanged; the exponent is a switch, `heat_exponent`).
* **Next:**
  1. **The release's steepness:** larger clouds (a deeper hold), a loss rising faster than |E|², and emitters of
     several strengths, in the same test bench; the target is released/held ∝ σ² over a factor of ten.
  2. **The release factor from the loss:** replace exp(−|g_N|/g_d) by the singers' hold in the whole suite (Cassini,
     the Milky Way, the dwarfs, lensing), and derive why it stops at zero.
  3. **The two KiDS analyses** (Brouwer et al.'s relation against Mistele et al.'s circular speeds): what in their
     methods makes one want more heat for ellipticals and the other less.
  4. Round 12's list: the guided companion with moving sources (the collisions' memory; the Bullet's smaller half),
     the lenses' gas, the distance law's shape, MACS J0025's star masses, the faint dwarfs, the Sun's speed, Abell 1689.

## 25. Round 15, 24 September 2026: an independent calculation, joined to the local force

The request: "We need to try and lock down a solid first principles proof. Explore this one next", with an
independent calculation, "Motion-opened radiation from ordinary matter's internal oscillations" (kept as received
in `independent-r15/first_principles_test/`). Its idea: each piece of matter holds a **quiet internal oscillation**
D, which carries its internal energy and barely radiates, and three **radiating** ones B_m. Random relative motion
shifts their frequencies apart by δ_m = χ w_m, and that mixes D into B linearly:
```
dD/dt = −γ₀ D − i Σ_m δ_m B_m        dB_m/dt = −γ B_m − i δ_m D
```
So the extra radiation is 3χ²σ²/(γ₀γ) times the quiet leak: the law's heat weight k = 3σ²/u², with
u² = γ₀γ/χ², from linear mixing followed by quadratic energy. Collisions (δ redrawn at rate ν) leave γ/(γ + ν)
of it. Its README asks for the next step: "use emitted waves from these evolving internal modes to drive test
emitters, derive forces from a common local interaction, and measure force, phase, energy and momentum
separately. Do not convert sqrt(radiated power) into a force."

In the order done:
1. the calculation reproduced (§25.1);
2. does the wave alone give its postulated coupling? (§25.2);
3. which test bodies a passing wave pulls, from their own equations (§25.3);
4. the integration test (§25.4);
5. what it costs, what it says about the constants, and next (§25.5).

### 25.1 The independent calculation, reproduced

`independent-r15/REPRODUCTION.md`. Rerun here (numba 0.67 instead of 0.65), four of its result files come out
byte-identical and the fifth (a control in a different wrapper) number for number.

* **Finite store (its part 3; γ = 1, γ₀ = 10⁻⁶, 512 pieces):** the extra radiation per doubling of the detuning
  is ×4.00, 3.99, 3.98, 3.92, 3.69, 3.01. It is σ² until the stores start to empty. Direction changes at
  ν = 1, 10 and 20 γ leave 50.7%, 9.2% and 4.8% of it, against γ/(γ + ν) = 50%, 9.1% and 4.8%. Stopped, the extra
  output dies away; only the spent fuel is gone.
* **A gain-reservoir oscillator (its part 1)** gives more to a weak incoming wave than it takes, and takes more
  from a strong one, with the switch at an incoming amplitude of 3.1225. There is no switch in its equations.
* **Its control (its part 2):** a fixed pump with no internal loss radiates exactly the pump at every speed. A
  fixed supply cannot shine brighter when hot; a finite store can.

### 25.2 Does the wave alone give the coupling?

`code/wave_dark_v15.py` → `run-reservoir-force-v15/wave_dark_v15.json` (176 runs, 4 minutes). Nothing is
postulated here. There are N = 40 ordinary oscillators in a ball of radius λ/4. Each has energy |a_j|² and
radiation rate 1, and they talk only through the scalar wave:
```
da_j/dt = −a_j/2 + (i/2) Σ_l C_jl a_l,   C_jl = cos(kR)/(k√(R² + ε²)) + i sin(kR)/(kR)
```
The soft core ε = 0.02λ acts on the near-field part only, so d(a†a)/dt = −(radiated power) stays exact. The
radiating part depends only on the distances between oscillators. A "dissipative" variant drops the near-field
part. Motion is measured as q = kσ in units of the radiation rate.

| extra radiation ÷ radiation at rest | q = 0.002 | 0.008 | 0.032 | 0.128 | colliding, q = 0.016: ν = 1, 10 | rigid rotation, uniform boost |
|---|---|---|---|---|---|---|
| over 200 time units, after settling 400 at rest (full) | 0.56 | 1.17 | 2.01 | 2.15 | 1.33, 0.72 (free: 1.74) | < 10⁻¹² |
| the same, dissipative | 3.5 | 13.4 | 31.8 | 40.3 | 5.3, 0.93 (free: 22.9) | < 10⁻¹² |

| steady leak rate P/E (second half of 3000) | at rest | q = 0.002 | 0.008 | 0.032 | 0.128 | ×per doubling of q |
|---|---|---|---|---|---|---|
| full | 7.5 × 10⁻⁴ | 3.9 × 10⁻³ | 6.8 × 10⁻³ | 1.64 × 10⁻² | 4.05 × 10⁻² | 1.29–1.59 |
| dissipative | 4.5 × 10⁻⁵ | 8.3 × 10⁻⁴ | 2.2 × 10⁻³ | 7.0 × 10⁻³ | 1.91 × 10⁻² | 1.53–1.82 |

* **Derived from the wave alone:**
  * Rigid rotation and uniform boosts release exactly nothing: the difference from rest is below 10⁻¹², because
    the wave's coupling depends only on the distances between emitters. The calculation's coupling, "the rate of
    change of the distance between particles", is the only kind of motion that can matter.
  * Free random motion opens the cloud's dark (subradiant) states.
  * Collisions hold that back, roughly as 1/(1 + ν/γ_eff). At q = 0.032 the extra leak with ν = 1, 10 and 100 is
    38%, 7.3% and 0.4% of the free one (dissipative; γ_eff ≈ 0.6), or 76%, 31% and 8% (full). This is Dicke
    narrowing.
* **Not derived: the σ².** The steady leak grows by ×1.3–1.8 per doubling of the speed (σ^0.4–0.9). That is
  the same shortfall as round 14's toy (×1.5–1.7 per doubling, §24.3). The spectra show why. The cloud's
  collective modes decay at rates spread over every decade:
  * full: the quietest decays at 1.2 × 10⁻⁴ to 1.5 × 10⁻³, and every decade from there to 10 holds 1–17 modes;
  * dissipative: every decade from 10⁻¹⁰ to 10 holds 1–6 modes.

  There is no gap between quiet and loud. Motion moves the state down this continuous ladder, and the leak settles
  where the ladder's crowding balances the scrambling: a fractional power of σ.
* **Reading.** σ² needs one quiet store separated by a wide gap from the radiating channels, γ₀ ≪ δ²/γ ≪ γ (in
  the calculation, 10⁻⁶ against 1). A cloud's collective quietness has no such gap. So the quiet store must be a
  property of each piece of matter, an internal oscillation, and not of how the pieces are arranged. This also
  explains round 14's steepness gap: its darkness was collective.

### 25.3 Which test bodies a passing wave pulls, from their own equations

`code/receivers_v15.py` → `run-reservoir-force-v15/receivers_v15.json`. One test body sits 6 wavelengths from a
steady source. It feels only the local wave:
* pull = −½Re(a* ∇E)·r̂;
* power it feeds the wave = (ω/2) Im(a* E);
* lead = sin(arg E − arg a): +1 means a quarter cycle ahead (feeding), −1 a quarter cycle behind (absorbing).

The self-sustained bodies start at 8 random phases, and the table gives their zero-detuning values. With a
detuning of half the locking range, each keeps its sign (lead ±0.78 to ±0.87 in the two weaker waves).

| body (its own equations) | pull at \|E\| = 0.0133 / 0.133 / 1.33 | lead | grows as |
|---|---|---|---|
| ordinary oscillator, driven (da/dt = −a/2 + iE) | −1.1 × 10⁻³ / −0.11 / −11 | −1 | intensity (pushed) |
| inverted, below its own threshold (stimulated emission: da/dt = −a/2 − iE) | +1.1 × 10⁻³ / +0.11 / +11 | +1 | intensity |
| rounds 10–14: fixed amplitude, held a quarter cycle ahead (the rule, imposed) | +0.042 / +0.42 / +4.2 | +1 | amplitude |
| **the calculation's gain-reservoir oscillator**, complex amplitude, physical coupling | **−0.071 / −0.76 / −13.3** | **−1** | amplitude (pushed) |
| **an inverted self-sustained emitter** (pumped ensemble kept oscillating by its own collective emission; Bloch equations; the wave's torque carries the inversion) | **+0.0138 / +0.129 / +0.537** | **+1** | **amplitude** (pulled) |
| the same, with the wave's torque given the ordinary sign (control) | −0.0138 / −0.134 / −1.24 | −1 | amplitude (pushed) |

* **Pull equals fed power over the wave speed in every case** (1.0000). Attraction requires feeding (§23.3).
* **The sign comes from the body's own oscillation.**
  * An ordinary self-sustained oscillator, the calculation's included, locks a quarter cycle *behind* a passing
    wave. It takes energy from the wave and is pushed, in proportion to the wave's amplitude. Its "gives more to
    weak waves" is its own free emission, which goes out evenly in all directions and carries no net momentum.
  * An **inverted** self-sustained emitter (one that holds energy it is ready to give, like the atoms of a laser)
    locks a quarter cycle *ahead*, by itself. It feeds the wave, and is pulled in proportion to the wave's
    amplitude: ×9.3 for ×10 in |E| in weak waves.
* **So round 10's rule is derived for inverted, self-sustained matter:** the quarter-cycle lead and the pull ∝
  amplitude. The Bloch equations are standard (the steady-state "superradiant laser"; Meiser et al., PRL 102,
  163601, 2009; Bohnet et al., Nature 484, 78, 2012). A body that only amplifies, below its own threshold, is
  pulled in proportion to the intensity. That falls off as 1/r², like Newton's pull, not the law's amplitude term.
* **In strong waves** (the `strong` part) the inverted emitter's pull stops growing. At |E| = 0.13, 1.3, 4.0, 13
  and 133 it is 0.129, 0.537, 0.602, 0.612 and 0.613: the wave drains the inversion (from 0.29 to 5 × 10⁻⁶), and the body
  feeds exactly as fast as its pump refills it. The lead stays +1. That is a ceiling, not the switch-off the Solar
  System needs; the release factor stays separate (§24.2).

### 25.4 The integration test

`code/reservoir_force_v15.py` → `run-reservoir-force-v15/reservoir_force_v15.json` (90 runs, 8 minutes).

**Sources.** 100 pieces of matter are fixed in a ball of radius λ (round 13's geometry). Each has its own D and
B_m, integrated exactly as the calculation writes them (γ = 1, γ₀ = 10⁻⁶), a random phase of D, and its own
detunings:
* at rest;
* free: fixed, Gaussian, rms q per component;
* colliding: redrawn at rate ν.

D's quiet leak radiates as a monopole and each B_m as a dipole along axis m. The strengths are set so that
the power each radiates into the field equals what its modes lose. That was checked through a sphere: 0.6000
against 0.6, 2.0000003 against 2.

**Test bodies.** 64 test bodies sit on a sphere of radius 6λ. Each feels only the local wave. Three kinds share
every wave:
* round 10's body, locked a quarter cycle ahead at rate Γ = 0.1, or Γ = 10;
* the inverted self-sustained emitter of §25.3 (coupling b = 50), whose timing nobody sets;
* an amplifier.

Nothing about the law enters.

Five arrangements for each case. Values are relative to rest, per unit of fuel left, to remove the fuel spent:

| q (rms detuning ÷ γ) | k = 3q²γ/γ₀ | released ÷ cold leak (1 + k) | wave intensity at the bodies (its √) | pull, round-10 body | pull, inverted body | pull, amplifier | fuel left at the end |
|---|---|---|---|---|---|---|---|
| 0.00025 | 0.19 | 1.19 (1.19) | ×1.27 (1.13) | ×1.13 ± 0.02 | ×1.15 ± 0.02 | ×1.27 | 1.00 |
| 0.0005 | 0.75 | 1.77 (1.75) | ×2.01 (1.41) | ×1.43 ± 0.06 | ×1.50 ± 0.06 | ×2.01 | 1.00 |
| 0.001 | 3 | 4.08 (4) | ×4.93 (2.20) | ×2.26 ± 0.14 | ×2.41 ± 0.16 | ×4.93 | 0.99 |
| 0.002 | 12 | 13.2 (13) | ×16.3 (4.00) | ×4.09 ± 0.30 | ×4.38 ± 0.34 | ×16.3 | 0.97 |
| 0.004 | 48 | 47.8 (49) | ×59.4 (7.62) | ×7.77 ± 0.59 | ×8.24 ± 0.66 | ×59.4 | 0.89 |
| 0.008 | 192 | 162 (193) | ×202 (14.0) | ×14.3 ± 1.1 | ×14.8 ± 1.2 | ×202 | 0.66 |
| 0.016 | 768 | 418 (769) | ×509 (22.3) | ×22.6 ± 1.7 | ×23.1 ± 1.8 | ×508 | 0.29 |

* **The heat term's full form comes out of the equations.**
  * The released power is 1 + k times the cold leak to within 3% up to k ≈ 50. Beyond that the fastest pieces
    burn their fuel first.
  * The round-10 bodies' pull follows the square root of the wave's intensity at their place to 1–3%. That is
    √(cold + released), the law's √(|g_N| + S_hot) with S_hot = k|g_N| for a point mass.
  * Nothing converts power into force: the pull is ½Re(a* ∇E), summed over bodies that feel only the local wave.
* **Doubling σ doubles the extra pull in the heat-dominated regime,** the decisive outcome asked for in round 13.
  The extra pull's ratio per doubling is 3.55, 2.95, 2.48, 2.20 and 1.97 (± 0.03–0.36), against the law's 3.60,
  3.10, 2.61, 2.30 and 2.15. The last step, 1.63, is where the fuel runs low.
* **The inverted bodies do the same with no timing put in.** Their lead is 0.77–0.93 in the faint wave at rest
  and 0.999–1.000 once motion brightens it. Their pull is ×1.15 … ×23.1.
* **The amplifier follows the intensity instead** (×1.27 … ×508): an intensity law, not the law's.
* **Collisions,** at q = 0.002 (k = 12):

| motion | extra radiation ÷ free | γ/(γ + ν) | slow round-10 bodies (Γ = 0.1): extra pull ÷ free, lead | fast ones (Γ = 10) | inverted bodies |
|---|---|---|---|---|---|
| ν = 0.3γ | 0.768 | 0.769 | 0.112, 0.33 | | 0.034, 0.24 |
| ν = γ | 0.502 | 0.500 | 0.039, 0.33 | 0.59, 0.99 | 0.015, 0.27 |
| ν = 3γ | 0.253 | 0.250 | 0.020, 0.41 | | 0.005, 0.35 |
| ν = 10γ | 0.094 | 0.091 | 0.010, 0.57 | 0.16, 0.96 | −0.002, 0.49 |
| ν = 30γ | 0.039 | 0.032 | 0.006, 0.71 | | −0.002, 0.61 |

  Collisions hold the heat back twice.
  * At the source, the extra radiation falls as γ/(γ + ν), to 1%.
  * At the test bodies, the colliding pieces' waves flicker. Bodies that re-time faster than the flicker keep
    step and are pulled by the square root of what is radiated (0.59 and 0.16, against 0.64 and 0.18 expected).
    Slower ones, and the inverted bodies with their weak coupling, lose step and gain almost nothing. The law's
    rule, "colliding matter adds no heat", holds whichever way the bodies re-time.
* **Other checks.**
  * *Moving in bulk as well,* the pieces' waves re-timed by their own Doppler shifts: extra pull ×1.02 (round-10
    bodies) and ×0.97 (inverted) of the fixed case.
  * *No memory:* moved at q = 0.004 until t = 900 and then stopped, the pull returns to 0.954 of the rest value.
    That is √(fuel left) = √0.915 = 0.957: only the spent fuel is missing.
* **Energy and momentum, measured separately.**
  * The stores drain exactly by what they radiate: at rest 0.9976 is left after 1200, which is e^(−2γ₀T).
  * The flux of the pieces' waves through a sphere is 0.60–1.40 times what their modes lose, depending on the
    arrangement (mean 1.02 over 85 runs). The pieces' waves interfere with each other, and that is not fed back
    into their modes. A fully consistent version would let each piece feel the others' waves.
  * Pull ÷ (fed power / wave speed) = 0.989–1.002 over every run, for every kind of body. It is below 1 because
    the sources are spread over the ball, so the pull is the radial part of the momentum fed.

### 25.5 What it costs, what it says about the constants, and next

* **Being pulled costs power: F = P/v_phase.** *Corrected in round 16 (§26.1):* this paragraph first used the
  travel speed u. A wave carries momentum k/ω = 1/v_phase per unit energy, so the phase speed sets the cost, as
  round 10 had it (§20.1, condition 3).
  * Pulled at the law's a = 6.30 × 10⁻¹¹ m/s², a kilogram feeds a·v_phase = 2(v_phase/u)ℓ. The galaxies allow
    v_phase up to about u/2 (§26.1), so at most ℓ, the companion power the kilogram emits when cold. That is at
    most 1.9 × 10⁻¹⁵ of its mass per year.
  * A light-speed phase would cost a·c = 3,500ℓ. The Sun, pulled at about 5 × 10⁻¹¹ m/s² by the Galaxy's
    companion, would then feed 78 solar luminosities into it and lose 5 × 10⁻¹² of its mass a year. Planetary
    ephemerides limit the change of the Sun's GM to about 10⁻¹³ per year (e.g. Pitjeva & Pitjev, MNRAS 432, 3431,
    2013).
  * **So energy and momentum require slow crests,** v_phase ≲ u/2 ≈ 85 km/s.
* **The constants.**
  * If a piece's internal energy is its rest energy, the quiet leak is γ₀ = ℓ/(2c²) = 3.0 × 10⁻²³ s⁻¹ (ℓ = au/2
    = 5.3 × 10⁻⁶ W/kg). Then u² = γ₀γ/χ² fixes γ/χ² = u²/γ₀ = 9.7 × 10³² m²/s.
  * The law's two rules bound the radiating modes' lifetime 1/γ:
    * stars carry heat, so their direction changes (orbits, ~200 million years at the Sun) must be slower than
      γ;
    * cluster gas does not, so its ions' gyration (~100 s in a microgauss field) must be faster.

    So 1/γ lies between about 100 seconds and 200 million years.
  * If the frequency shift is a Doppler shift of the companion itself (χ = its wavenumber), the companion's
    wavelength lies between about 60 parsecs and 500 megaparsecs.
* **Where the proof chain stands:**

| link | status after round 15 |
|---|---|
| random motion releases companion power ∝ σ², collisions hold it back as γ/(γ + ν), rotation and bulk motion release none | derived in the reduced model, given a quiet internal store with a gap and a frequency shift ∝ the rate of change of distances (§25.1, §25.4). The wave alone gives the symmetries and the collision rule, not the gap (§25.2) |
| released and cold power add, and reach test bodies as waves | derived, energy checked (§25.4) |
| the pull is the amplitude, √(cold + released), a quarter cycle ahead | derived for inverted, self-sustained test bodies from their Bloch equations (§25.3). Every run has pull = fed power / wave speed |
| doubling σ doubles the extra pull where heat dominates | measured: ×2.48, 2.20, 1.97 per doubling (§25.4) |
| the energy bill | a·v_phase per kilogram; ≤ ℓ with v_phase ≲ u/2, which the galaxies allow; a light-speed phase is excluded by the planets (§25.5, corrected in §26.1) |
| u from microscopic rates | u² = γ₀γ/χ²: one relation between two unknown rates, not yet a number |
| the release factor (strong fields) | still separate: the inverted body saturates rather than switching off |

* **Suite** unchanged: 59 pass, 11 close, 7 fail (the law is unchanged).
* **Next:**
  1. **A physical quiet store with a gap:** an internal two-oscillator structure whose frequencies are pulled apart
     by relative motion. The motional mixing of a metastable atomic state is the laboratory example. Derive χ
     and γ, and so u.
  2. **One kind of matter:** every piece both a source (a quiet store opened by motion) and a receiver (inverted
     and self-sustained), in one cloud with mutual waves fed back. That removes the ±40% interference noise.
  3. **The release factor from the receivers:** whether a body's inversion, drained by its own neighbours' waves in
     strong fields, can switch its pull off rather than cap it.
  4. Round 14's list: the two KiDS analyses; round 12's list.

## 26. Round 16, 24 September 2026: the next steps, one at a time

The request: "proceed with the next steps, keep the blog and main current": round 15's list (§25.5). In the
order done:
1. the energy bill of being pulled, fed back into the companion, on the galaxies and clusters (§26.1); this also
   corrects §25.5;
2. what opens the quiet store (§26.2);
3. one kind of matter, with every wave fed back (§26.3).

### 26.1 The energy bill of being pulled, fed back into the companion

`code/feeding_feedback_v16.py` and `code/feeding_checks_v16.py` → `run-feeding-feedback-v16/`.

**A correction first.** §25.5 priced the pull with the companion's travel speed u. A wave carries momentum
k/ω = 1/v_phase per unit energy, so a body pulled with the companion's extra acceleration g_c feeds
P/m = g_c v_phase into it. That is round 10's condition 3 (§20.1), and §25.5 is corrected above.

**The fed power joins the stream.** It goes into the wave passing the body, which flows on outward. So the stream's
power through radius r is ℓM(r), each kilogram's own feed, plus v_phase ∫₀ʳ g_c dM, the feeding. With β =
v_phase/u and a = 2ℓ/u, the ordered companion's intensity becomes
```
I(r) = |g_N(r)| (1 + 2β W(r)/(a M(r))),   W(r) = ∫₀ʳ g_c dM,   g_c = exp(−|g_N|/g_d) √(a (I + S))
```
solved from the centre out (spherical bookkeeping, as for |g_N| = GM(r)/r²). β = 0 is the law. Round 10
estimated that β ≈ 1 would spoil the tight v⁴ = G M a and asked for β ≲ 0.05. Here the data measure it.

| β = v_phase/u | a (10⁻¹¹ m/s²) | g_d (10⁻¹⁰ m/s²) | 149 galaxies: typical miss (km/s) | fit statistic | outer residual against W/(aM) |
|---|---|---|---|---|---|
| 0 (the law) | 6.30 | 2.03 | 15.87 | 0.2111 | +0.074 ± 0.062 |
| 0.05 | 6.04 | 2.03 | 15.83 | 0.2103 | |
| 0.2 | 5.30 | 1.97 | 15.76 | 0.2082 | |
| 0.5 | 3.99 | 1.71 | 15.79 | 0.2051 | −0.023 ± 0.034 |
| 1 | 2.32 | 1.29 | 16.35 | 0.2054 | −0.054 ± 0.019 |

a and g_d are refitted at each β, with u held. The residual slope is log₁₀(v²_obs/v²_pred) at each galaxy's last
point against W/(aM) there, with 400 bootstrap resamples of the galaxies.

* **Up to β ≈ 0.5 the galaxies are as happy as with the law** (15.76–15.83 against 15.87 km/s). A smaller a
  absorbs most of the extra intensity: a·(1 + 2β⟨W/aM⟩) stays close to 6.3 × 10⁻¹¹.
* **β ≈ 1 is ruled out.**
  * The outer residuals trend with W/(aM) at 2.8σ, the other way from β = 0.
  * Fitted on the 89 training galaxies alone, β runs to 0.95, but the 60 held-out galaxies then do worse than with the
    law: 19.88 against 19.10 km/s (validation) and 13.80 against 12.40 (test).
  * The all-galaxy fit statistic's 3% gain at β = 0.5–1 does not survive out of sample.
* **The clusters accept any of these** with u lowered: X-COP's rms is 0.222 at β = 0 (u = 170 km/s), 0.223 at
  β = 0.5 (u = 142) and 0.228 at β = 1 (u = 113).
* **Reading.** The crests must move at no more than about half the companion's travel speed, v_phase ≲ 85 km/s.
  * Round 10's bound (≲ 0.05u) was ten times too strict: a refitted a takes up most of the feedback.
  * At the law's a the bill is then at most ℓ per kilogram, and at most 1.9 × 10⁻¹⁵ of its mass per year.
  * A light-speed phase is out twice over: the planets (§25.5), and the galaxies.
  * The law is unchanged (β = 0 fits as well as β = 0.5). A small surface-density trend in the outer residuals
    (+0.055 ± 0.025 at β = 0, +0.032 at β = 0.5) is noted for later.
