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
(rms z 4.36); ellipticals 1.09, 0.98, 1.01 (0.95). Two separate things have moved. *(Correction, 25 September: the
four ratios are bins of stellar mass, log M* 10.10–11.29, each averaged over 50–300 kpc, not four radii; §30.)*

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
* **Reading** *(narrowed in round 16, after an independent audit)*. σ² needs one quiet store separated by a wide gap
  from the radiating channels, γ₀ ≪ δ²/γ ≪ γ (in the calculation, 10⁻⁶ against 1).
  * The tested cloud's collective quietness does not provide that separation.
  * An internal quiet store is a concrete candidate that would. This one construction does not rule out every
    collective arrangement.
  * The same lack of a gap is a likely reason for round 14's steepness shortfall, whose darkness was collective.

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

* **A conditional demonstration of the heat term's form.** *(Round 16: narrowed after an independent audit; round 15
  first wrote "the heat term's full form comes out of the equations", which was premature.)*
  * The mixing δ = χw is an input here, not derived, and the pieces sit at fixed points.
  * Given it, the released power is 1 + k times the cold leak to within 3% up to k ≈ 50. Beyond that the fastest
    pieces burn their fuel first.
  * The round-10 bodies (their lead imposed) are pulled by the square root of the wave's intensity at their place,
    to 1–3%. That is √(cold + released), the law's √(|g_N| + S_hot) with S_hot = k|g_N| for a point mass.
  * **The inverted bodies, whose lead comes out of their own equations,** are the first-principles result. Their lead
    is 0.77–0.93 in the faint wave at rest and 0.999–1.000 once motion brightens it. Their pull, ×1.15, 1.50, 2.41,
    4.38, 8.24 and 14.8, follows √(1 + k) = 1.09, 1.32, 2.00, 3.61, 7.00 and 13.9 in trend, 10–20% above it. Finite
    geometry, interference and the receivers' own response all enter, so "exact" is not established.
  * Nothing converts power into force: the pull is ½Re(a* ∇E), summed over bodies that feel only the local wave.
  * A pull ∝ √(intensity) outside one cloud is not yet the full law. The spatial integral, the pull's direction in an
    asymmetric system, the potential equation and light bending are still to be derived.
* **Doubling σ doubles the extra pull where heat dominates,** with the imposed-lead bodies: the extra pull's ratio
  per doubling is 3.55, 2.95, 2.48, 2.20 and 1.97 (± 0.03–0.36), against the law's 3.60, 3.10, 2.61, 2.30 and 2.15.
  The last step, 1.63, is where the fuel runs low. The inverted bodies show the same trend.
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
    into their modes.
  * **This is the biggest gap** (round 16, from an independent audit):
    * The audit reproduced the mismatch (seed 2: 0.6039; seed 3: 1.229; unchanged from 400 to 8,000 points on the
      sphere).
    * It showed why with two emitters a tenth of a wavelength apart: in phase they radiate 1.94 times their assigned
      loss, and in opposite phase 0.065 times.
    * It tested a repair, damping taken from the same coupling that radiates: ż = −iHz − ½W†Wz, with W the map from
      internal amplitudes to outgoing waves. That closes the energy budget exactly.
    * In its reduced test (20 sources, no test bodies, the mixing still assumed), the extra radiation per doubling of
      the detuning fell from 4.13, 4.05, 4.00, 3.91 to 3.80, 3.28, 2.02, 1.64.
    * So the σ² scaling has to be re-established with the feedback included (§26.3). Averaging over arrangements is
      not enough: each cloud needs a budget that closes for its own state.
  * Pull ÷ (fed power / wave speed) = 0.989–1.002 over every run, for every kind of body. It is below 1 because
    the sources are spread over the ball, so the pull is the radial part of the momentum fed.

### 25.5 What it costs, what it says about the constants, and next

* **Being pulled costs power: F = P/v_phase.** *Corrected in round 16 (§26.1):* this paragraph first used the
  travel speed u. A wave carries momentum k/ω = 1/v_phase per unit energy, so the phase speed sets the cost, as
  round 10 had it (§20.1, condition 3). This toy's wave has the same phase and travel speeds; whether the
  companion's crests move more slowly than its energy is still to be established.
  * Round 15's "exactly 2ℓ at v = u" followed from the law's own definition a = 2ℓ/u. It was an identity, not a
    check.
  * Pulled at the law's a = 6.30 × 10⁻¹¹ m/s², a kilogram feeds a·v_phase = 2(v_phase/u)ℓ. The galaxies allow
    v_phase up to about u/2 (§26.1), so at most ℓ, the companion power the kilogram emits when cold. That is at
    most 1.9 × 10⁻¹⁵ of its mass per year.
  * A light-speed phase would cost a·c = 3,500ℓ, about 1,770 times the cost at 169 km/s. The Sun, pulled at about
    5 × 10⁻¹¹ m/s² by the Galaxy's companion, would then feed 78 solar luminosities into it and lose 5 × 10⁻¹² of
    its mass a year.
  * Planetary ephemerides limit the change of the Sun's GM to about 10⁻¹³ per year (e.g. Pitjeva & Pitjev, MNRAS
    432, 3431, 2013). Taken as a simple budget, that bounds the crests' speed at about 5,700 km/s. Applying it needs
    two more steps: how this power loss changes the Sun's mass and GM, and the Sun's other losses.
  * **So energy and momentum require slow crests.** The planets allow up to about 5,700 km/s; the galaxies about
    85 km/s (§26.1). Neither measures u = 169 km/s; both are compatible with it.
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
| random motion releases companion power ∝ σ², collisions hold it back as γ/(γ + ν), rotation and bulk motion release none | conditional: shown in the reduced model with independent radiation channels, given a quiet store with a gap and the mixing δ = χw as inputs (§25.1, §25.4). The wave alone gives the symmetries and the collision rule, not the gap (§25.2). With the shared wave's feedback the scaling changes (independent audit); to be re-established (§26.3) |
| released and cold power add, and reach test bodies as waves | shown, but the shared wave's energy budget does not close (flux 0.60–1.40 × the modes' loss) |
| the pull is the amplitude, √(cold + released), a quarter cycle ahead | an attractive local response demonstrated: inverted, self-sustained test bodies choose the lead themselves (§25.3), pulled ∝ \|E\| in weak waves; 10–20% above √(1 + k) in the integration. Every run has pull = fed power / wave speed |
| doubling σ doubles the extra pull where heat dominates | measured with imposed-lead bodies: ×2.48, 2.20, 1.97 per doubling (§25.4); to be redone with the feedback |
| the energy bill | a·v_phase per kilogram; ≤ ℓ with v_phase ≲ u/2, which the galaxies allow; a light-speed phase is excluded by the planets (§25.5, corrected in §26.1) |
| u from microscopic rates | u² = γ₀γ/χ²: one relation between two unknown rates, not yet a number |
| the release factor (strong fields) | still separate: the inverted body saturates rather than switching off |

* **Suite** unchanged: 59 pass, 11 close, 7 fail (the law is unchanged).
* **Next** *(reordered in round 16, as an independent audit recommends)*:
  1. **One kind of matter, with the feedback from the start:** every piece both a source (a quiet store opened by
     motion) and a receiver (inverted and self-sustained, powered by the same finite store), one local coupling
     giving its emission, its loss, its response and its force, in one cloud with every wave fed back so the energy
     budget closes. Then cold, free and colliding motion again. No strong-field switch in that test.
  2. **A physical quiet store with a gap:** an internal two-oscillator structure whose frequencies are pulled apart
     by relative motion. The motional mixing of a metastable atomic state is the laboratory example. Derive χ
     and γ, and so u.
  3. **The release factor from the receivers:** whether a body's inversion, drained by its own neighbours' waves in
     strong fields, can switch its pull off rather than cap it.
  4. Round 14's list: the two KiDS analyses; round 12's list.

## 26. Round 16, 24 September 2026: the next steps, one at a time

The request: "proceed with the next steps, keep the blog and main current": round 15's list (§25.5). In the
order done:
1. the energy bill of being pulled, fed back into the companion, on the galaxies and clusters (§26.1); this also
   corrects §25.5;
2. what opens the quiet store, and the crossing heat in the Bullet Cluster (§26.2);
3. the shared wave with its energy budget closed, as an independent audit of round 15 asked (§26.3);
4. one kind of matter, every piece both sender and receiver, with every wave fed back (§26.4).

The audit also narrowed several round-15 claims; §25 is corrected in place and marked.

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

### 26.2 What opens the quiet store: the companion flowing past

`code/stream_store_v16.py` → `run-stream-store-v16/`; `code/crossing_heat_v16.py` → `run-crossing-heat-v16/`.

**The hypothesis.** The independent calculation (§25.1) mixes a piece's quiet store into its radiators by δ = χw,
with w "a relative velocity" left unspecified, and it needs a separate quiet leak γ₀. Here w is the piece's velocity
relative to the companion flowing through it, w = v − V_flow. The companion streams away from the matter that emits it
at u and keeps its sources' mean motion (§12, the memory). There is no quiet leak of its own (γ₀ = 0). The equations
are the calculation's own, integrated as it writes them (γ = 1; 8 seeds × 64 pieces per point).

| case | measured | expected |
|---|---|---|
| at rest in its own flow: the leak per unit store | 2.0002 × 10⁻⁴ | χ²u²·2/γ = 2 × 10⁻⁴ |
| free random motion, σ/u = 0.25, 0.5, 1, 2, 4, 8: released ÷ cold | 1.16, 1.70, 3.88, 12.7, 48.6, 198 (± 1–5%) | 1 + 3σ²/u² = 1.19, 1.75, 4, 13, 49, 193 |
| σ = 2u, collision rate ν = 1, 3, 10, 30 (units of γ): extra | 5.90, 2.96, 1.07, 0.383 | 12γ/(γ + ν) = 6, 3, 1.09, 0.387 |
| a steady drift along the flow, v_r/u = −2, −1, −0.5, +0.5, +1, +2 | 9.01, 4.00, 2.25, 0.250, 0.000, 1.00 | (1 − v_r/u)² |
| moving across the flow at v/u = 0.5, 1, 2: with memory / without | 1.000, 1.000, 1.000 / 1.25, 2.00, 5.00 | 1 / 1 + v²/u² |

* **The leak comes from the flow.** A piece at rest in its own outflowing companion sees |w| = u, so its cold leak is
  χ²u²/γ. The calculation's relation u² = γ₀γ/χ² then holds identically, with u the companion's travel speed. The law
  has used one u for both jobs, the travel speed and the heat scale in k = 3σ²/u², since round 3; here the two are the
  same speed by construction.
* **The factor 3 comes out too:** three random velocity components against the flow's one.
* **Collisions** scramble only the random part, as the law's rule for gas needs.
* **Ordered rotation carries no heat only if the flow co-moves with its source** (the memory). Without it, a disk
  turning at 200 km/s through a flow fixed in space would count as hot (k ≈ 1.4). So this reading requires the
  memory that the colliding clusters already need. Differential rotation (shear) is not yet tested.
* **A new prediction: drift.** Matter falling inward against the flow leaks more, (1 + |v_r|/u)²; matter flowing out
  at u leaks nothing. A steady drift survives collisions, so even infalling gas would take this heat.
* **Status: conditional.** The quiet store with its gap and the linear mixing are still the calculation's postulates
  (§25.5). What this adds is what the relative velocity is, and that it makes u one speed.

**The crossing heat.** If the store opens to the velocity relative to the flow, matter moving through another
system's flow is heated. A piece sees its own system's flow (share 1 − p) and the other's (share p, moving outward
from the other system at u and carrying its velocity), so
```
k = 3σ²/u² + p (v² − 2uv cos θ)/u²
```
with v its speed relative to the other system and θ the angle between its motion and the other flow (receding
cos θ = +1, approaching −1). Gas collides, so it takes none.

The Bullet Cluster (round 11's static-distance model, `code/bullet_static_v11.py`) is rerun with this one term added.
* p is the other system's companion intensity (|g_N| + S, each system's own spherical model before the collision) over
  the sum, at each of four shells of stars.
* The history is a straight pass: receding at v_out since pericentre, approaching at v_in before it, from where the
  other system's flow ends (u × 13 Gyr).
* The heat emitted a time t ago sits a distance ut from the stars now (the memory), so it enters through a
  distance-dependent kernel.
* Grid 192², 15 kpc. The lensing targets span Bradač et al.'s and Paraficz et al.'s masses within 250 kpc (§14), in our
  distances.

| case | M(<250 kpc), smaller half (10¹⁴ M☉; target 2.47–2.85) | main half (target 3.09–3.46) | κ at the galaxies, smaller / main | crossing heat ÷ all heat (mid-plane) |
|---|---|---|---|---|
| the law, no crossing heat | 1.385 | 2.969 | 0.276 / 0.717 | 0 |
| v_out 3,900 km/s, v_in 3,000, impact parameter 150 kpc | 1.603 | 3.119 | 0.489 / 0.766 | 0.22 |
| slower: 2,700 and 2,000 | 1.556 | 3.073 | 0.443 / 0.750 | 0.16 |
| faster: 4,700 and 4,000 | 1.641 | 3.160 | 0.521 / 0.780 | 0.26 |
| 3,900 and 3,000, head-on | 1.608 | 3.130 | 0.496 / 0.772 | 0.23 |
| 3,900 and 3,000, impact parameter 300 kpc | 1.590 | 3.106 | 0.472 / 0.760 | 0.21 |

* **The crossing heat lifts the smaller half by 12–18%** (1.56–1.64 against 1.39), and moves the main half into its
  measured range (3.07–3.16 against 3.09–3.46) in all but the slowest crossing.
* The gas residuals come out 0.050 and 0.018–0.024, against the measured 0.05 and 0.02 (the law alone: 0.050 and
  0.052), and the lensing peaks sit on the galaxies (6–15 kpc from them).
* **Not enough alone.** The smaller half still has 63–66% of the lower end of its measured mass; the term closes about
  a fifth of the gap. The speeds and the impact parameter matter little: the heat is dominated by v² ≫ u², and the
  share p where the stars sit.
* Next for this term: MACS J0025, Abell 520 and El Gordo, where it adds heat to galaxies that crossed; and the dwarfs
  moving through the Milky Way's flow.

### 26.3 The shared wave with its energy budget closed

`code/shared_wave_v16.py` → `run-shared-wave-v16/`.

An independent audit of round 15 found the combined experiment's biggest gap (§25.4): each piece lost energy at its own
rate while the pieces' waves were added coherently. It proposed the repair, damping taken from the same coupling that
radiates: ż = −iHz − ½W†Wz, so that d‖z‖²/dt = −P_out exactly. In its reduced test (20 pieces, dissipative coupling
only, the mixing still assumed) the extra radiation per doubling of the detuning fell from 4.13, 4.05, 4.00, 3.91 to
3.80, 3.28, 2.02, 1.64. Its conclusion: the heat law has to be re-established with the feedback included.

**The complete coupling.** Eliminating the wave between pieces with monopole (quiet store D, rate γ₀) and dipole
(radiators B, rate γ) channels gives one complex symmetric matrix M, built from the exact Green's function e^{ikR}/(4πR)
and its derivatives (near and far field). Its imaginary part Γ = W†W is the radiated power; its real part J is the
reactive near field, which moves energy between pieces but radiates none. A passive piece obeys ż = −iHz + (i/2)Mz.
Checks (`validate()`):
* z†Γz equals the flux through a distant sphere for 1, 2 and 12 pieces in random states (7.01166 against 7.01162;
  29.7195 against 29.7195; 164.391 against 164.390);
* the forces on all pieces balance the momentum radiated, to 10⁻⁵–10⁻⁴ (the sphere's quadrature), including the recoil
  of a piece whose monopole and dipole are coherent;
* two monopoles a tenth of a wavelength apart radiate 1.9355 and 0.0645 times their own-rate loss in and out of phase,
  the audit's numbers.

**Three versions of the damping, same clouds.**
* *independent*: each piece's own rate (round 15);
* *dissipative*: −½Γ, the audit's reduced repair;
* *shared*: the whole coupling, −½Γ + (i/2)J.

The clouds are 20 pieces in a ball of radius 0.5, 1 or 3 wavelengths (at least 0.15 apart), and round 15's 100 pieces
in radius 1. The pieces' own frequencies are identical, or spread with rms 0.3γ. γ₀ = 10⁻⁶, γ = 1, and the rms
detuning q runs 0.00025 … 0.016 (k = 0.19 … 768 for an isolated piece). Stores start full at random phases. The
evolution is the exact matrix exponential; 5 arrangements each.

Extra radiation per unit stored energy, (P/E)_q/(P/E)_0 − 1, ratio per doubling of q (means over arrangements; cold
leak relative to an isolated piece's 2γ₀):

| cloud, damping | cold leak | t = 50: q 0.00025 → … → 0.016 | t = 1000 |
|---|---|---|---|
| any cloud, independent | 1 | 4.00, 4.00, 4.00, 3.99, 3.97, 3.87 | 4.00, 3.99, 3.95, 3.83, 3.39, 2.46 |
| 20 in radius 1, dissipative | 0.03–0.05 | 3.98, 3.99, 3.98, 3.95, 3.81, 3.27 | **3.88, 3.53, 2.45, 1.42**, 1.36, 1.35 |
| 20 in radius 1, shared | 0.91 | 3.50, 3.71, 3.84, 3.91, 3.90, 3.75 | 3.50, 3.70, 3.77, 3.63, 3.01, 2.18 |
| 20 in radius 0.5, shared | 0.91 | 3.28, 3.56, 3.75, 3.85, 3.87, 3.73 | 3.28, 3.54, 3.68, 3.57, 2.94, 1.88 |
| 100 in radius 1, shared | 0.82 | 9.05, 5.11, 4.43, 4.17, 4.00, 3.70 | 9.10, 5.08, 4.30, 3.71, 2.81, 2.12 |
| 20 in radius 3, shared | 1.00 | 4.08, 4.04, 4.02, 4.00, 3.96, 3.84 | 4.08, 4.03, 3.96, 3.80, 3.30, 2.48 |
| 20 in radius 1, frequencies spread, shared | 0.92 | 3.99, 4.00, 3.99, 3.99, 3.94, 3.79 | 5.07, 4.40, 4.12, 3.82, 3.18, 2.23 |

* **The budget closes for each cloud's own state.** Stored energy lost equals the radiated energy to (0.25–1.4) × 10⁻⁵
  over t = 3000 (trapezoid quadrature), and the flux through a distant sphere equals z†Γz to 10⁻⁵.
* **The audit's result is reproduced by the dissipative half alone** (bold, t = 1000; its 3.80, 3.28, 2.02, 1.64). The
  mechanism: in a dense cloud the radiators, driven by their neighbours' quiet monopoles, re-radiate against them. That
  screens the cold leak down to 5% of an isolated piece's or less, and the doubling ratio falls once the fastest stores
  have emptied.
* **The complete coupling removes most of the screening.** The reactive near field is the real part of the same
  Green's function whose imaginary part is Γ; causality ties the two. With it, the cold leak stays at 82–100% of an
  isolated piece's, and the σ² survives: 3.3–3.9 per doubling in dense clouds, 4.0 in dilute ones, until the fastest
  stores empty. Their emptying lowers the late ratios for independent damping just as much.
* **A dense cloud's heat coefficient is renormalized.** At q = 0.002 (k = 12 for an isolated piece) the extra is
  11.7, 9.9, 6.9 and 4.3 for the clouds of radius 3, 1, 0.5 and the 100-piece cloud, against 11.8–12.6 with
  independent damping. In the law that reads as a somewhat larger u for tightly packed matter, a correction to size
  rather than a failure.
* **A term linear in the detuning** appears at the smallest q: one piece's quiet monopole interferes with another's
  radiators. It depends on the stores' phases, averages away over them and shrinks as 1/√N. It is why some ratios at
  the smallest q depart from 4 (9.05 for the 100-piece cloud).
* **Collisions,** q = 0.002, the detunings redrawn at ν = 1 and 10 (units of γ), share of the free extra kept at
  t = 200: 0.55 and 0.091 (radius 3), 0.48 and 0.095 (radius 1), against γ/(γ + ν) = 0.50 and 0.091. The densest
  cloud (radius 0.5) keeps more, 0.98 and 0.26, with large scatter.
* **Reading.** With the energy fed back, the motion-opened σ² and the collision rule hold in dilute and moderately
  dense clouds and survive, renormalized, in dense ones. The mixing δ = χw and the quiet store's gap are still inputs.

### 26.4 One kind of matter: every piece both sender and receiver

`code/one_matter_v16.py` → `run-one-matter-v16/`.

The audit's recommended construction: "Keep the quiet and bright internal modes. Give every piece the same
finite-energy reservoir and the same local coupling to the wave. Derive its outgoing radiation, its energy loss, its
response to incoming radiation and its mechanical force from that coupling. Include mutual feedback from the
beginning." Round 15 had two descriptions of matter: sources with a finite quiet store, and receivers with an applied
pump. Here there is one, and sources and receivers differ only in where they sit and how they move.

**The piece.** Every piece holds:
* **a quiet oscillation s**, the coherence of an inverted ensemble. It holds energy ready to give, which is the kind
  of body round 15 found locks ahead of a passing wave. It radiates as a monopole at rate γ₀.
* **its inversion w**, with energy (1 + w)/2, and **a finite store n** that refills it at W₀n/n₀, debited exactly.
* **an internal collective gain G** that keeps the quiet oscillation going (its output leaves by an internal channel,
  tallied). This is what makes the piece self-sustained, as in a superradiant laser.
* **three passive radiators B**, dipoles that radiate at rate γ = 1 and also lose energy inside at γ_i. Their share
  of energy that goes into the companion is f = γ/(γ + γ_i).
* **the motion's mixing** of s into the B's, δ = χw. This is the independent calculation's postulate, still an input.

**The coupling.** Every piece, source or receiver, talks to every other through the single complex symmetric matrix M
of §26.3 (the exact Green's function, near and far field). The equations are
```
ds/dt = −(i/2) w (M z + K z)_s + (G/2) w s − λ s        (−w: inverted; an ensemble with w = −1 responds like a passive oscillator)
dB/dt = +(i/2)   (M z + K z)_B − γ_i B
dw/dt = W (1 − w) − g∥ (1 + w) − 2G|s|² − 2 Im(s* (M z + K z)_s),     dn/dt = −W (1 − w)/2
```
with K the mixing block [[0, −2δᵀ], [−2δ, 0]]. Each channel loses exactly the power it feeds the wave, so the pieces'
total loss is z†Im(M)z: the interference between pieces is inside every piece's own budget. The force on a piece is
(1/ω) Re Σ z_j†∂M_jl z_l, from the same M. Nothing converts power into force, no piece's timing is imposed, and there
is no strong-field switch.

**Set-up.** 48 source pieces in a ball of radius 3 wavelengths (at least 0.15 apart), and 8 receivers of the same
matter on a sphere of radius 6, at rest. The constants are γ₀ = 0.01, G = 4, W₀ = 1.6, g∥ = 0.01, γ_i = 4 (so
f = 0.2), and n₀ = 10⁵, so the stores lose at most 4% in a run. The sources are at rest, move freely (δ fixed,
rms q per component), or collide (δ redrawn at rate ν). Each run lasts 8,000 time units and is measured over the
second half. There are three arrangements, and the isolated-piece heat weight is k = 3q²γ/(γ₀(γ + γ_i)²).

**First, an exact recoil rule** (`one_matter_v16.py --recoil` → `recoil_check.json`: one receiver 6 wavelengths
from a monopole source, γ₀ = 0.02, its quiet oscillation held a quarter cycle ahead of the source's wave, its radiators
in their steady response):

| radiators' share to the companion f | quiet channel, pulled by the source's wave | recoil from its own radiators | radiators | net, lead +1 |
|---|---|---|---|---|
| 1 (γ_i = 0) | +1.050 × 10⁻⁴ | −1.050 × 10⁻⁴ | −1.126 × 10⁻⁴ | −1.126 × 10⁻⁴ (pushed) |
| 1/2 (γ_i = 1) | +1.050 × 10⁻⁴ | −0.525 × 10⁻⁴ | −0.563 × 10⁻⁴ | −0.038 × 10⁻⁴ |
| 0.1 (γ_i = 9) | +1.050 × 10⁻⁴ | −0.105 × 10⁻⁴ | −0.113 × 10⁻⁴ | +0.833 × 10⁻⁴ (pulled) |

* The quiet channel's pull is exactly the power it feeds the wave, divided by the wave's speed.
* A piece's own radiators, driven by the same wave, radiate coherently with its quiet oscillation. Their interference
  takes back a fraction f of the pull on the quiet channel and f again on the radiators: **net (1 − 2f) × fed/c**,
  minus the radiators' ordinary radiation pressure, which grows as the wave's intensity.
* Why, in one line: a plane wave A e^{ik·x} drives the radiators to B = −(ωkC₁/(4(γ + γ_i))) n̂A. The coherent
  monopole-dipole recoil, −(C₀C₁k³/(24π)) Im(s*B), is then exactly −f times the quiet channel's pull, −(C₀k/2) n̂ Im(s*A)
  (using C₁² = 48πγ/(ωk³)), and the same again on the radiators' side.
* The rule holds at 6 and at 20 wavelengths. A full run with a cold source (`--set probe`: γ_i = 9, γ₀ = 0.02, round
  15's receiver constants) gave a net pull of 4.89 × 10⁻⁴ against 0.8 × 6.18 × 10⁻⁴ = 4.94 × 10⁻⁴ from the fed power.
* **So one kind of matter can be pulled only if its radiators send less than half of their energy into the
  companion** and keep the rest inside. With f = 1 a piece that feeds the wave is pushed. This is a new, sharp
  condition on what matter must be like, and it comes from the equations, not from a choice.

**The decisive test: cold, warm and colliding sources, one matter throughout.** Means over three arrangements (the
spread between arrangements in brackets where it matters). "Wave" is the sources' wave at the receivers, "from the
sources" the pull on the receivers' quiet oscillation from that wave, "keeping step" its lead sin(arg E − arg s) on it:

| sources | k (isolated piece) | glow ÷ rest | wave ÷ rest (√ of the glow) | pull from the sources ÷ rest | keeping step | net pull on a receiver | radiators | sources' rhythm spread |
|---|---|---|---|---|---|---|---|---|
| at rest | 0 | 1 | 1 | 1 | 0.55 (± 0.08) | **+1.21 × 10⁻⁴** | −1.22 × 10⁻⁴ | 4.7 × 10⁻⁴ |
| free, q = 0.20 | 0.48 | 1.21 | 1.01 (1.10) | 0.56 | 0.26 | +0.43 × 10⁻⁴ | −0.97 × 10⁻⁴ | 7.7 × 10⁻⁴ |
| free, q = 0.41 | 2.0 | 2.11 | 1.31 (1.45) | 0.71 | 0.27 | +0.25 × 10⁻⁴ | −1.52 × 10⁻⁴ | 1.2 × 10⁻³ |
| free, q = 0.82 | 8.1 | 5.00 | 2.24 (2.24) | 1.06 | 0.27 | −1.04 × 10⁻⁴ | −3.47 × 10⁻⁴ | 2.6 × 10⁻³ |
| free, q = 1.15 | 15.9 | 8.07 | 2.81 (2.84) | 0.91 | 0.16 | −2.92 × 10⁻⁴ | −5.08 × 10⁻⁴ | 3.9 × 10⁻³ |
| colliding, ν = 5 | 8.1 | 3.92 | 1.91 (1.98) | 0.87 | 0.25 | −0.20 × 10⁻⁴ | −2.27 × 10⁻⁴ | 3.9 × 10⁻⁴ |
| colliding, ν = 50 | 8.1 | 1.62 | 1.29 (1.27) | 1.06 | 0.44 | **+1.07 × 10⁻⁴** | −1.48 × 10⁻⁴ | 4.3 × 10⁻⁴ |

* **The books balance.** Energy (stores, ensembles, radiators, internal output and radiated power against the start)
  closes to 10⁻¹² with free motion and 4 × 10⁻⁷ with collisions. The forces balance the momentum carried away to the
  accuracy of the far-field quadrature (10⁻³–6 × 10⁻², 2,000 directions around a 12-wavelength arrangement).
* **Same matter, pulled.** A cold source pulls distant pieces of the same matter: they fall a quarter beat ahead by
  themselves (0.55 on the sources' wave over the measured stretch; 0.92 once settled, follow-up 3) and are pulled net,
  +1.21 × 10⁻⁴, after their radiators' push. The pull on the
  quiet channel is the fed power over the wave's speed, less the recoil rule's share.
* **The source's glow follows the heat rule.** Glow ×1.21, 2.11, 5.00 and 8.07 at k = 0.5, 2, 8 and 16. That is
  somewhat below 1 + k, because the inverted mixing draws on the store's throughput.
* **The wave carries the law's square root.** Its amplitude at the receivers is ×1.01, 1.31, 2.24 and 2.81, against √ of
  the glow, 1.10, 1.45, 2.24 and 2.84.
* **But the receivers fall out of step.** Their lead drops from 0.55 to 0.16–0.27, so the pull on them stays flat (×0.56
  … 1.06). Meanwhile their radiators' push grows with the glow (−1.2 to −5.1 × 10⁻⁴). The net pull turns into a push at
  k ≥ 8. **The decisive test is not passed** in this construction.
* **Why: warming puts the source's own pieces out of tune.** The spread of the sources' rhythms grows from 4.7 × 10⁻⁴
  at rest to 3.9 × 10⁻³ at k = 16.
  * A warm piece's quiet oscillation drives its radiators, and the radiators feel their neighbours. That shifts each
    piece's rhythm by an amount that depends on its own motion (∝ |δ|²).
  * The pieces of a cold source settle into a common beat. A warm source's pieces can't, and a wave with a smeared beat
    can't be followed from a distance.
  * The receivers' re-timing rate, (w/2)|E|/|s|, against the sources' rhythm spread: 1.6 at rest (keeping step 0.55),
    then 0.98, 0.82, 0.58 and 0.43 as the source warms (0.26, 0.27, 0.27, 0.16). Warming strengthens the wave, which
    speeds up the re-timing as √ of the glow, but spreads the rhythms faster.
  * Frequent collisions (ν = 50) keep the rhythms together (4.3 × 10⁻⁴), and the receivers keep step again (0.44) and
    are pulled as by the cold source. At ν = 5 the rhythms stay together too, but the radiators' part of the wave
    flickers at the collision rate, far faster than the receivers can re-time (0.25).
* **What this settles and what it leaves.**
  * Settled: in one kind of matter with every wave fed back, energy closes. A piece pulled by feeding a wave loses
    (1 − 2f) of it to its own radiators. And the extra glow of warm matter reaches distant matter as the square root of
    its power, as the law needs.
  * Open: distant matter has to keep step with a warm source. In this construction warmth scrambles the source's own
    beat faster than it strengthens its hold on distant matter.

**Follow-up 1: radiators that whisper** (f = 0.05: γ_i = 19; G = 40 and W₀ = 16 so the extra drain stays small;
dt = 0.01; two arrangements):

| sources | k | glow ÷ rest | wave ÷ rest (√ glow) | pull from the sources ÷ rest | keeping step | net pull | radiators | rhythm spread |
|---|---|---|---|---|---|---|---|---|
| at rest | 0 | 1 | 1 | 1 | 0.40 | **+1.64 × 10⁻⁴** | −0.22 × 10⁻⁴ | 3.5 × 10⁻⁴ |
| free, q = 1.63 | 2.0 | 1.83 | 1.18 (1.35) | 0.56 | 0.19 | **+1.00 × 10⁻⁴** | −0.25 × 10⁻⁴ | 9.3 × 10⁻⁴ |
| free, q = 3.27 | 8.0 | 4.73 | 2.40 (2.17) | 1.05 | 0.21 | **+1.25 × 10⁻⁴** | −0.73 × 10⁻⁴ | 2.2 × 10⁻³ |
| colliding, ν = 200 | 8.0 | 1.74 | 1.38 (1.32) | 0.84 | 0.28 | **+1.31 × 10⁻⁴** | −0.25 × 10⁻⁴ | 3.9 × 10⁻⁴ |

* **The push is gone as a problem.** The radiators' push is 4–7 times smaller, and the receivers are pulled net at
  every heat.
* **Keeping step is not fixed.** The sources' rhythms spread just as much with warmth (3.5 × 10⁻⁴ → 2.2 × 10⁻³), the
  receivers' lead falls (0.40 → 0.19–0.21), and the pull still doesn't grow. So the detuning does not travel through
  the radiators' share. A likely route, still to be tested: each warm piece tugs its neighbours' quiet oscillations
  through its stirred radiators. At a given heat weight that coupling scales as γ₀√(k/3), whatever f is.

**Follow-up 2: the cause confirmed.** A cold source, no heat at all, with its pieces' own rhythms spread by hand (rms
Δ₀), against the main run's rest (two arrangements):

| cold source, rhythms spread by hand | glow ÷ rest | wave ÷ rest | sources' rhythm spread | keeping step | pull from the sources ÷ rest |
|---|---|---|---|---|---|
| Δ₀ = 0 (the main run) | 1 | 1 | 4.7 × 10⁻⁴ | 0.55 | 1 |
| Δ₀ = 0.001 | 0.96 | 0.99 | 7.0 × 10⁻⁴ | 0.28 | 0.46 |
| Δ₀ = 0.0025 | 0.83 | 1.02 | 2.4 × 10⁻³ | 0.12 | 0.21 |
| Δ₀ = 0.004 | 0.79 | 1.03 | 3.9 × 10⁻³ | 0.19 | 0.39 |

* Being out of tune, by itself, loses the receivers. It does so at the same spreads as warmth does: 0.28 at 7.0 × 10⁻⁴
  (warm, k = 0.5: 0.26 at 7.7 × 10⁻⁴), and 0.12–0.19 at 2.4–3.9 × 10⁻³ (warm, k = 8–16: 0.27–0.16, a little better
  because the warm wave is stronger).
* **So the one obstacle is now pinned down:** warm, collisionless matter's pieces must stay in tune closely enough
  (rhythm spread well below the receivers' re-timing rate, (w/2)|E|/|s|) for distant matter to follow their wave.

**Follow-up 3: given time** (`code/one_matter_rhythm_v16.py`: the main run's matter and constants, arrangement 1,
24,000 time units instead of 8,000, reported every 1,000):

| source | sources' rhythm spread | the wave's jitter at the receivers | receivers keeping step |
|---|---|---|---|
| at rest | 1.2 × 10⁻³ at t = 1,000, settling to 3–4 × 10⁻⁴ from t ≈ 4,000 | 1–3 × 10⁻⁴ once settled | −0.35 → 0.54 (t = 5,000) → **0.82–0.92** (t = 15,000–24,000) |
| free, k = 2 | 1.3–2.3 × 10⁻³ throughout | 1–16 × 10⁻³ | wanders between −0.24 and +0.45 (mean 0.10 over the last 12,000) |

* **A cold source settles into one beat, and distant pieces of the same matter lock onto it almost perfectly**, a
  quarter cycle ahead, by themselves (0.91–0.92 over the last 4,000). The main runs' 0.55 was measured while this was
  still settling.
* **A warm source never settles.** Its pieces stay out of tune, its wave jitters ten times more, and the receivers
  never lock.

### 26.5 Where round 16 leaves the proof

| link | status after round 16 |
|---|---|
| the energy bill of being pulled | a·v_phase per kilogram; the galaxies allow v_phase ≲ u/2, so at most ℓ (§26.1) |
| random motion releases companion power ∝ σ², collisions hold it back, rotation none | shown with the shared wave's energy balanced, 3.3–4.0 per doubling (§26.3), given the quiet store's gap and the mixing δ = χw. With δ measured against the flowing companion, u is one speed for both of the law's jobs (§26.2, conditional) |
| a piece pulled by feeding a wave | exact: the fed power over the wave's speed on its quiet channel, less 2f to its own radiators' recoil (§26.4) |
| released and cold power reach distant matter as waves, as √ of the glow | shown in one kind of matter: the wave's amplitude at the receivers ×2.24 for glow ×5.00 (§26.4) |
| distant matter keeps step by itself | yes with a cold source (0.92 once settled) and with frequently colliding ones (0.44); **not with a warm collisionless source** (0.16–0.27), whose own pieces fall out of tune. A cold source put out of tune by hand loses the receivers the same way (§26.4) |
| a piece's radiators | must send less than half their energy into the companion (the recoil rule); at f = 0.05 the push is small and the pull net at every heat (§26.4) |
| the complete energy-conserving law | not yet: the pull on distant matter does not grow with the source's heat in this construction |
| u from microscopic rates | the flow picture makes u the travel speed by construction; χ and γ still unknown |
| the release factor (strong fields) | still separate; no strong-field switch was added |

* **Suite** unchanged: 59 pass, 11 close, 7 fail (the law is unchanged).
* **Next:**
  1. **Keeping warm matter in tune.** The obstacle is now one number: the spread of a warm source's rhythms against
     distant matter's re-timing rate. Test what the tugging between warm pieces depends on (the density, the size of
     the source, whether a large source locks itself to one beat), and whether a matter rule that lets the radiators
     shift only the store's strength, not its rhythm, keeps the beat. Then add real motion's Doppler shifts to the same
     budget.
  2. **A physical quiet store with a gap:** derive χ and γ.
  3. **The crossing heat** in MACS J0025, Abell 520, El Gordo and the dwarfs.
  4. Round 14's list: the two KiDS analyses; round 12's list.

## 27. Round 17, 24 September 2026: keeping warm matter in tune

Round 16 (§26.4) ended on one open link: in one kind of matter, a warm source's pieces fall out of tune and distant
matter cannot keep step with its wave, so the pull does not grow with the source's heat. An independent review set this
round's objective: "derive a mechanism that lets motion increase the companion's output without destroying the rhythm
needed to produce attraction", in three steps:
* **A**, find exactly which term of the energy-balanced equations shifts a warm piece's rhythm, as an equation with its
  terms measured;
* **B**, a structural protection, the same inside every piece, not phase resets; its candidate was a pair of internal
  radiating responses at ±Δ whose pulls cancel;
* **C**, require the net pull, not just the wave, to follow: R(σ, r) = (F(σ, r)/F(0, r))·√(I(0, r)/I(σ, r)) ≈ 1 at
  increasing distance, for cold, free, colliding and stopped sources, with the microscopic constants fixed, the power in
  proportion to mass, and the same acceleration for any receiving matter. Only then the astrophysical fits.

### 27.1 Step A: the rhythm budget, and the term that does the damage

`code/rhythm_budget_v17.py` → `run-rhythm-budget-v17/rhythm_budget_v17.json`.

**The reduction.** Each piece's radiators relax at γ + γᵢ = 5, a thousand times faster than the rhythms move, so they can
be eliminated exactly in linear response. With every piece at its steady inversion w₀ and amplitude, the rhythms obey
```
dφ_j/dt = ω_j − (w₀/2) Σ_{l≠j} |C_jl| cos(φ_l − φ_j + arg C_jl),      ω_j = −(w₀/2) Re C_jj,
C = M_ss − (i/2) Xᵀ Y⁻¹ X,     X = M_Bs − 2D,     Y = (i/2) M_BB − γᵢ
```
(D holds each piece's mixing δ_j). ω_j is a piece's own rhythm offset: zero for an isolated piece, nonzero in a cloud,
where its radiators' response picks up its neighbours' waves. C_jl is how piece l tugs piece j. The imaginary part of
C is a "gradient" coupling: it acts like the ordinary pull toward a common beat, and a system coupled only that way
always settles. The real part is not: it shifts a piece's rhythm by an amount that depends on its neighbours' phases,
which can keep a system from ever settling.

**The budget is exact.** In the full model (arrangement 1), each piece's rhythm equals the sum of three terms to
1.5 × 10⁻⁷ (correlation 1.0000):

| source | rhythm spread | tug of the others' quiet oscillations | of their radiators | of its own radiators, through its mixing |
|---|---|---|---|---|
| at rest | 4.4 × 10⁻⁴ | 4.8 × 10⁻⁴ | 3.6 × 10⁻⁴ | 0 |
| k = 2 | 1.2 × 10⁻³ | 8.6 × 10⁻⁴ | 1.1 × 10⁻³ | 1.4 × 10⁻³ |
| k = 8 | 2.6 × 10⁻³ | 6.1 × 10⁻⁴ | 1.7 × 10⁻³ | 2.4 × 10⁻³ |

**The reduced model reproduces the full one** (three arrangements): rhythm spread 3.3 × 10⁻⁴, 1.1 × 10⁻³ and 3.0 × 10⁻³
at k = 0, 2 and 8 (full model, §26.4: 4.7 × 10⁻⁴, 1.2 × 10⁻³, 2.6 × 10⁻³), receivers keeping step 0.63, 0.31, 0.17
(full: 0.55, 0.27, 0.27).

**The equation, term by term.** X is linear in the mixing and Y does not depend on it, so C = C₀ + C₁ + C₂ exactly: the
cold coupling (radiators passive), a part linear in δ (one piece's quiet oscillation reaching another's through one set
of stirred radiators) and a part quadratic in δ (through two). Means over three arrangements, with k = 12q² for these
constants:

| heat weight k | 0.5 | 2 | 8 | 16 | fitted |
|---|---|---|---|---|---|
| own offsets: cold part | 2.0 × 10⁻⁴ | 2.0 × 10⁻⁴ | 2.0 × 10⁻⁴ | 2.0 × 10⁻⁴ | 2.0 × 10⁻⁴ |
| own offsets: linear part | 2.3 × 10⁻⁴ | 4.7 × 10⁻⁴ | 9.3 × 10⁻⁴ | 1.3 × 10⁻³ | 3.3 × 10⁻⁴ √k |
| own offsets: quadratic part | 0.9 × 10⁻⁴ | 3.9 × 10⁻⁴ | 1.6 × 10⁻³ | 3.1 × 10⁻³ | 1.9 × 10⁻⁴ k |
| tugs on a piece (rms sum): cold | 1.9 × 10⁻³ | 1.9 × 10⁻³ | 1.9 × 10⁻³ | 1.9 × 10⁻³ | 1.9 × 10⁻³ |
| tugs: linear | 1.9 × 10⁻³ | 3.8 × 10⁻³ | 7.6 × 10⁻³ | 1.1 × 10⁻² | 2.7 × 10⁻³ √k |
| tugs: quadratic | 0.9 × 10⁻³ | 3.8 × 10⁻³ | 1.5 × 10⁻² | 3.0 × 10⁻² | 1.9 × 10⁻³ k |

So the warm source's pieces are tugged by their neighbours' stirred radiators up to 16 times harder than the cold
coupling that lets a source settle, and about half of every heat-induced tug is of the non-gradient kind (gradient share
0.49 linear, 0.53 quadratic).

**Which part does the damage** (`code/rhythm_protect_v17.py`, diagnostics; three arrangements; the receivers keep step on
the sources' part of their actual drive):

| round 16's matter, parts of the heat switched off by hand | keeping step, k = 2 | k = 8 |
|---|---|---|
| everything (the model as it is) | 0.34 | 0.18 |
| **heat tugs: only their gradient part** | **0.75** | **0.88** |
| heat tugs: only their non-gradient part | 0.15 | 0.00 |
| no heat tugs among the source's pieces | 0.66 | 0.63 |
| no heat offsets (tugs kept) | 0.36 | 0.24 |
| for comparison: a cold source | 0.63 | 0.63 |

* **The term that puts warm matter out of tune is the non-gradient part of the tugs between warm pieces**, the part that
  shifts a piece's rhythm by an amount set by its neighbours' phases.
* **With it removed, warmth helps.** Keeping only the gradient part of the heat tugs, the warm source settles (rhythm
  spread 3.9 × 10⁻⁴ at k = 8) and distant matter keeps step at 0.88, better than with a cold source, because the
  warm source's wave is stronger.
* The own offsets matter much less: removing them alone barely helps (0.36, 0.24).

### 27.2 Step B: what an internal structure can and cannot protect

`code/paired_channel_v17.py` → `paired_channel_v17.json`; `code/rhythm_protect_v17.py` → `rhythm_protect_v17.json`,
`rhythm_protect_compact_v17.json` (all in `run-rhythm-budget-v17/`).

**The review's pair, reproduced for one piece.** A quiet mode D mixed by the motion (strength g) into two radiating
families at +Δ and −Δ, each radiating at γ through its own channel:

| check | result | review |
|---|---|---|
| added decay per doubling of g (0.0025 → 0.02), Δ = γ | ×4.00000, 4.00000, 4.00000 | ×4.00008, 4.00030, 4.00120 |
| largest rhythm shift | 4.3 × 10⁻¹⁹ | < 1.1 × 10⁻¹⁸ |
| opposite 1% deviations of the two couplings: shift ÷ added decay | 0.0200 | about 2% |
| collisions at ν = γ: extra glow kept, pair / one resonant family | 0.800 / 0.500 | 80% |

* At Δ = γ the first correction to the g² law cancels as well, so we get ×4 exactly. The review's small excesses are what
  a single resonant family gives, 4(1 + 3g₁²/γ²); either way the added decay grows as g².
* The pair keeps more of its glow under collisions than a single family (0.8 against 0.5). Our law needs collisions to
  hold the heat back, so this is a cost to watch.

**Why the pair's protection is lost in a shared wave.** The review's cancellation protects the rhythm against a common
shift of the two families' own frequencies: to first order the sensitivity is proportional to Σ± 1/(γ ± iΔ)², which
vanishes at Δ = γ. In one kind of matter both families radiate into the same wave, and what the surroundings change is
the channel they share (the sum B₊ + B₋), not each family's frequency. The first-order sensitivity is then proportional
to χ², with χ = Σ± 1/(γᵢ ± iΔ) the channel's own response, which is real and never zero:

| Δ/γ | 0 | 0.5 | 1 | 2 |
|---|---|---|---|---|
| rhythm shift per unit reactive change of the surroundings, per unit added decay: separate channels | 1.00 | 0.60 | **0.00** | −0.60 |
| the same, both families in one shared channel (f = 0.2) | 1.67 | 1.32 | 0.82 | 0.32 |

**A general limit.** Twelve internal structures, all energy-consistent and identical in every piece, compared at equal
glow in the reduced model (three arrangements):

| structure | mixing q at k = 8 | keeping step, k = 0 | k = 2 | k = 8 | k = 16 | rhythm spread, k = 8 | own offsets, k = 8 | tugs, k = 8 |
|---|---|---|---|---|---|---|---|---|
| single (round 16) | 0.82 | 0.63 | 0.34 | 0.18 | 0.17 | 2.7 × 10⁻³ | 1.9 × 10⁻³ | 1.65 × 10⁻² |
| single, odd parity | 0.82 | 0.63 | 0.34 | **0.33** | 0.23 | 2.2 × 10⁻³ | 1.6 × 10⁻³ | 1.70 × 10⁻² |
| paired, Δ = γ + γᵢ (the review's) | 1.41 | 0.71 | 0.33 | 0.18 | 0.16 | 2.6 × 10⁻³ | 1.6 × 10⁻³ | 1.66 × 10⁻² |
| paired, Δ = γ + γᵢ, odd | 1.41 | 0.71 | 0.39 | 0.26 | 0.17 | 2.2 × 10⁻³ | 1.3 × 10⁻³ | 1.71 × 10⁻² |
| paired, Δ = (γ + γᵢ)/2 | 0.87 | 0.53 | 0.35 | 0.18 | 0.15 | 3.3 × 10⁻³ | 2.5 × 10⁻³ | 1.64 × 10⁻² |
| paired, Δ = γ + γᵢ, 1% asymmetry | 1.41 | 0.71 | 0.28 | 0.09 | 0.07 | 3.2 × 10⁻³ | 3.8 × 10⁻³ | 1.66 × 10⁻² |
| paired, Δ = 4(γ + γᵢ) | 12.2 | 0.83 | 0.40 | 0.25 | 0.19 | 2.3 × 10⁻³ | 1.9 × 10⁻⁴ | 1.71 × 10⁻² |
| **paired, Δ = 4(γ + γᵢ), odd** | 12.2 | **0.83** | 0.39 | **0.41** | **0.26** | 1.6 × 10⁻³ | 1.5 × 10⁻⁴ | 1.76 × 10⁻² |
| bright + dark partner, Δ = γ + γᵢ | 1.63 | 0.49 | 0.33 | 0.23 | 0.17 | 3.1 × 10⁻³ | 1.5 × 10⁻³ | 1.69 × 10⁻² |
| bright + dark partner, odd | 1.63 | 0.49 | 0.44 | 0.34 | 0.22 | 2.0 × 10⁻³ | 1.4 × 10⁻³ | 1.73 × 10⁻² |
| valve (damped chamber, gR = 20, gB = 5) | 4.08 | 0.68 | 0.44 | 0.34 | 0.23 | 1.8 × 10⁻³ | 1.4 × 10⁻³ | 1.69 × 10⁻² |
| valve, odd | 4.08 | 0.68 | **0.67** | 0.21 | 0.20 | 7.5 × 10⁻⁴ | 1.3 × 10⁻³ | 1.70 × 10⁻² |

* **At equal glow the heat tugs have the same size in every structure**, 1.64–1.76 × 10⁻² at k = 8. This is
  reciprocity: a piece that sends its motion-opened glow into the shared wave receives its neighbours' glow through the
  same channel, with the same strength. And the Onsager–Casimir relations forbid a mixing linear in the velocity from
  being one-way (sending without receiving) at every speed. So no internal structure can remove the heat tugs. It can
  only change their character, and a piece's own offsets.
* **The velocity's parity, a correction from first principles.** A velocity changes sign when time runs backwards. A
  coupling proportional to it, between two internal modes that do not, must then be imaginary and antisymmetric (like
  the Coriolis force, or the drag of a flowing medium on a sound wave), not real and symmetric as round 16 and the
  independent calculation had it. Both are energy-conserving, and for one isolated piece they give the same glow. In a
  cloud the odd form makes the linear part of the coupling antisymmetric, so each piece's own linear offset vanishes
  identically (10⁻¹⁹ against 3.3 × 10⁻⁴ √k). Distant matter's keeping step at k = 8 goes from 0.18 to 0.33.
* **The review's pair at a large detuning makes each piece a poor scatterer.** Families at ±4(γ + γᵢ) answer an outside
  wave weakly (their responses nearly cancel) while the motion still drives them hard. The pieces' own offsets fall
  tenfold (1.5 × 10⁻⁴ at k = 8), a cold source settles better (keeping step 0.83 against 0.63), and with the odd parity
  this is the best structure at k = 8 and 16 (0.41, 0.26).
* **The valve with the odd parity keeps a mildly warm source fully in tune.** At k = 2, rhythm spread 3.0 × 10⁻⁴ (cold:
  3.3 × 10⁻⁴) and keeping step 0.67 (cold: 0.68).
* **A second requirement, found with the valve.** At k = 8 the valve-odd source does settle (arrangement 2: spread
  1 × 10⁻⁴ from t = 10,000 on), but all together at a rhythm shifted by +3.1 × 10⁻³, while the cold receivers stay near
  +1.0 × 10⁻³ and cannot follow. A source's beat has to stay where cold matter's is, or distant matter, whose hold on the
  beat weakens with distance, eventually cannot follow at all.
* **A source about a wavelength across is worse**, even cold: keeping step 0.02–0.04 at radius 0.6 wavelengths
  (corrected 25 September: the label read "wavelengths/2π"; the models use wavelength 1), where the pieces' near
  fields dominate every coupling (`rhythm_protect_compact_v17.json`). This ball is also 125 times denser than the
  main runs, so size and density are not separated here.

**The velocity's parity in the full model** (`code/one_matter_v17.py --set parity` → `run-one-matter-v17/parity.json`;
round 16's set-up, three arrangements, 8,000 time units; `code/one_matter_summary_v17.py`):

| source | glow ÷ rest | keeping step: round 16's mixing / the velocity's parity | pull of the sources' wave on the receivers' quiet channel ÷ rest (√ of the glow) | net pull on a receiver |
|---|---|---|---|---|
| at rest | 1 | 0.55 / 0.55 | 1 | +1.21 × 10⁻⁴ / +1.21 × 10⁻⁴ |
| free, k = 2 | 2.07 | 0.38 / **0.46** | 1.02 / **1.14** (1.44) | +0.52 / +0.73 × 10⁻⁴ |
| free, k = 8 | 5.00 | 0.23 / **0.33** | 0.85 / **1.22** (2.24) | −1.11 / −0.48 × 10⁻⁴ |
| free, k = 16 | 8.10 | 0.16 / **0.24** | 0.70 / **1.13** (2.85) | −2.95 / −2.43 × 10⁻⁴ |
| colliding, k = 8, ν = 5 | 3.9 | 0.25 / 0.16 | 0.87 / 0.62 | −0.19 / −0.70 × 10⁻⁴ |
| colliding, k = 8, ν = 50 | 1.6 | 0.44 / 0.48 | 1.06 / 1.09 | +1.07 / +1.13 × 10⁻⁴ |

* The full model confirms the reduced one: with the velocity's parity distant matter keeps step better at every heat,
  and the pull of the warm source's wave on it now grows with the heat (×1.14, 1.22 and 1.13 of the cold value), where
  round 16's form made it shrink (×0.85 and 0.70 at k = 8 and 16). It still grows only about half as fast as the
  square root of the glow that our law needs. Energy closes to 10⁻¹² (free) and 4 × 10⁻⁷ (collisions).
* The net pull still turns negative at k ≥ 8. That is the receivers' own radiators, pushed by the brighter wave in
  proportion to its intensity (−1.2 × 10⁻⁴ at rest, −3.1 at k = 8, −5.0 at k = 16): round 16's second problem, which
  radiators that mostly ring inside (f = 0.05) remove (§26.4, follow-up 1).

### 27.3 A wave that only travels outward: a hypothesis test

`code/one_way_v17.py` → `run-rhythm-budget-v17/one_way_v17.json`, `one_way_mass_v17.json`, `one_way_order_v17.json`.

Step A says what does the damage: the two-way exchange between warm pieces (each piece's rhythm pushed around by
neighbours whose own rhythms it pushes back). Step B says that no structure inside a piece can remove it, because at a
given glow reciprocity fixes its size. What would remove it is a wave that cannot come back. The companion streams
outward from the matter that emits it (§26.2). If its crests are carried outward by that stream faster than they can
move against it (the review's "one medium for crest and transport speeds"), a piece can hear only matter nearer the
source's centre, and never its own echo.

**The test.** The reduced model of §27.1, with the coupling made one-way before the reduction: the part carrying a wave
from piece l to piece j is kept when l is nearer the source's centre than j, and multiplied by a leak ε otherwise (ε = 1
is round 16's two-way wave). Each piece's own self-coupling is kept. Receivers of the same matter sit on four shells,
radius 6, 9, 13.5 and 20 (four on each). In the reduced model a locked receiver's pull is its lead times the wave's
amplitude, so R(σ, r) = lead(σ, r)/lead(0, r). Three arrangements, 16,000 time units (the second half measured):

| wave | source | sources' rhythm spread | keeping step at r = 6 | 9 | 13.5 | 20 | R(σ, r) at 6, 9, 13.5, 20 |
|---|---|---|---|---|---|---|---|
| two-way (ε = 1), round 16's matter | at rest | 3.1 × 10⁻⁴ | 0.52 | 0.29 | 0.28 | 0.48 | 1 |
| | k = 2 | 1.5 × 10⁻³ | 0.31 | 0.25 | 0.34 | 0.11 | 0.60, 0.86, 1.23, 0.24 |
| | k = 8 | 4.5 × 10⁻³ | 0.19 | 0.23 | 0.10 | 0.09 | 0.36, 0.79, 0.35, 0.19 |
| | k = 16 | 8.8 × 10⁻³ | 0.12 | 0.06 | −0.01 | 0.06 | 0.24, 0.20, −0.02, 0.12 |
| **outward only (ε = 0)** | at rest | 1.5 × 10⁻⁴ | 0.88 | 0.79 | 0.59 | 0.76 | 1 |
| | k = 2 | 7.7 × 10⁻⁵ | 0.97 | 0.73 | 0.83 | 0.58 | 1.11, 0.92, 1.41, 0.77 |
| | k = 8 | 6.3 × 10⁻⁵ | **0.98** | **0.96** | **0.85** | **0.81** | **1.11, 1.21, 1.45, 1.07** |
| | k = 16 | 6.0 × 10⁻⁶ | **1.00** | **0.97** | **0.84** | **0.77** | **1.14, 1.22, 1.43, 1.01** |
| outward only, odd parity | k = 8 | 1.1 × 10⁻⁵ | 0.97 | 0.92 | 0.89 | 0.82 | 1.11, 1.17, 1.51, 1.07 |
| | k = 16 | 3.0 × 10⁻⁵ | 0.87 | 0.94 | 0.94 | 0.97 | 0.99, 1.20, 1.60, 1.27 |

* **With a wave that only travels outward, warmth no longer spoils the beat.** The warm source's pieces share one
  rhythm (spread 6 × 10⁻⁵ at k = 8, 6 × 10⁻⁶ at k = 16, against 4.5 × 10⁻³ and 8.8 × 10⁻³ two-way), its common beat
  stays where cold matter's is (within 3 × 10⁻⁵), and distant matter keeps step at every distance tested, better the
  warmer the source. **R(σ, r) is 1.0 to 1.5 at all four distances,** where the two-way wave gives 0.1–0.4 once the
  source is hot.
* **Why:** no piece hears its own echo, so there are no own offsets (exactly zero); the innermost matter sets the beat
  and everything farther out locks onto it a quarter cycle ahead, with nothing coming back to disturb it.
* **What matters is that the wave has no way back, not the direction.** With the sources put in a random one-way order
  instead (each hears only those before it; `one_way_order_v17.json`), a warm source again keeps one beat (spread
  8 × 10⁻⁶ at k = 8) and distant matter keeps step at 1.00, 0.99, 0.81 and 0.84 (k = 8) and 0.99, 0.98, 0.85, 0.93
  (k = 16). Any coupling without loops protects the beat; an outflow supplies the natural order, by distance from the
  centre.
* **How strict it has to be:**

| inward leak ε | 1 | 0.5 | 0.3 | 0.1 | 0 | 0, but two-way inside a core of half the radius |
|---|---|---|---|---|---|---|
| keeping step at r = 6, 9, 13.5, 20, k = 8, odd parity | 0.23, 0.14, 0.15, 0.12 | 0.28, 0.16, 0.17, 0.04 | 0.25, 0.22, 0.18, 0.10 | 0.84, 0.77, 0.61, 0.32 | 0.97, 0.92, 0.89, 0.82 | 0.66, 0.73, 0.43, 0.34 |

  The protection needs the wave to be nearly one-way (an inward leak of about 10% or less) throughout the source; a
  two-way core of half the radius loses much of it.
* **The power stays in proportion to mass.** A chain of pieces each locking ahead of the passing wave could build it up
  coherently, so that the power grew as the square of the mass. It does not: at the same density, the wave's intensity
  at the outer shells per source piece is about the same for 24, 48 and 96 pieces (at r = 13.5: 9.8, 5.2, 5.6 × 10⁻⁸
  cold; 3.1, 3.2, 3.7 × 10⁻⁷ at k = 8), and keeping step improves with size (0.92–1.00 at k = 8 with 96 pieces).
* **A consequence worth noticing.** With a wave that only travels outward, a piece feels only the matter nearer the
  centre than itself. For a round source that is exactly the rule that Newton's gravity obeys (only the mass inside a
  radius pulls there), and it is the form our law already takes, through g_N of the enclosed mass.
* **What this is and is not.** It is a test of a hypothesis about the companion, made by hand in the reduced model:
  the one-way coupling is imposed, not derived, and the reduced model tracks rhythms, not energy or forces. Making it
  real means deriving the companion as a flowing medium (the crests' speed against its outflow speed, and whether the
  wave's energy is exchanged with the flow) and then building the one-way coupling into the full, energy-balanced model.

### 27.4 Step C in the full model: the one-way wave with every watt booked

`code/one_matter_v17.py --set oneway_distance` → `run-one-matter-v17/oneway_distance.json`; summary by
`code/one_matter_summary_v17.py` → `oneway_distance_summary.json`.

**The test.** Round 16's full model (§26.4: 48 source pieces of one kind of matter, every piece coupled to every other
through the complete wave, near field and far field), with §27.3's one-way coupling imposed on the whole coupling, quiet
oscillations and radiators alike: the wave from piece l reaches piece j only when l is nearer the centre (ε = 0), and
each piece keeps its own self-coupling. What the pieces then give the wave, Im(z†Mz), is shared between the outgoing
wave and the flow that carries it, and the recoil of the upstream pieces is taken up by the flow; both are booked as "to
the wave and its flow", and the books balance to 2 × 10⁻¹² (3 × 10⁻⁷ with collisions). Receivers of the same matter sit
on four shells, r = 6, 9, 13.5 and 20 (four on each). Sources at rest, free at k = 2 and 8, colliding (k = 8, ν = 50),
and warm (k = 8) but stopped at t = 4,000. Two arrangements each, 16,000 time units, the second half measured. Two kinds
of matter: round 16's single family ("single"), and the same with the velocity's time-reversal parity of §27.2
("single, odd"). R_net(σ, r) = (F(σ, r)/F(0, r))·√(I(0, r)/I(σ, r)) with F the net pull on the receivers and I the
source wave's intensity at them; "from the quiet channel" uses only the pull of the source wave on the receivers' quiet
oscillations.

**One-way wave, round 16's matter (single)**

| source | glow ÷ rest | keeping step (all receivers) | net pull at r = 6, 9, 13.5, 20 (10⁻⁴) | R_net(σ, r) at 6, 9, 13.5, 20 | R from the quiet channel |
|---|---|---|---|---|---|
| at rest | 1.00 | 0.63 | +1.61, +1.14, +0.41, +0.53 | 1 | 1 |
| free, k = 2 | 1.98 | 0.78 | +2.17, +1.64, +1.28, +0.54 | 0.85, 1.05, 1.83, 0.77 | 1.06, 1.25, 2.00, 0.75 |
| free, k = 8 | 4.75 | 0.93 | +2.36, +2.12, +1.75, +1.07 | 0.68, 1.01, 1.65, 0.89 | 1.09, 1.24, 1.96, 1.25 |
| colliding, k = 8, ν = 50 | 1.62 | 0.46 | +1.53, +0.99, +0.38, +0.50 | 0.72, 0.66, 0.65, 0.66 | 0.81, 0.72, 0.78, 0.71 |
| k = 8, stopped at t = 4,000 | 1.01 | 0.62 | +1.85, +1.06, +0.63, +0.11 | 1.06, 0.88, 1.70, 0.24 | 1.08, 1.01, 1.46, 0.42 |

**One-way wave, with the velocity's parity (single, odd)**

| source | glow ÷ rest | keeping step (all receivers) | net pull at r = 6, 9, 13.5, 20 (10⁻⁴) | R_net(σ, r) at 6, 9, 13.5, 20 | R from the quiet channel |
|---|---|---|---|---|---|
| at rest | 1.00 | 0.63 | +1.61, +1.14, +0.41, +0.53 | 1 | 1 |
| free, k = 2 | 2.05 | 0.87 | +2.01, +1.74, +1.40, +0.99 | 0.90, 1.05, 1.62, 0.97 | 1.03, 1.16, 2.04, 1.08 |
| free, k = 8 | 4.78 | 0.79 | +1.50, +1.80, +1.73, +0.88 | 0.59, 0.86, 1.67, 0.64 | 0.93, 1.08, 2.12, 0.90 |
| colliding, k = 8, ν = 50 | 1.61 | 0.47 | +1.45, +1.06, +0.37, +0.49 | 0.71, 0.71, 0.64, 0.66 | 0.80, 0.77, 0.74, 0.72 |
| k = 8, stopped at t = 4,000 | 1.00 | 0.70 | +1.67, +1.15, +0.57, +0.76 | 0.99, 1.00, 1.11, 1.02 | 1.00, 1.03, 1.19, 1.13 |

Keeping step by shell at k = 8: 0.98, 0.98, 0.84, 0.93 (single) and 0.94, 0.86, 0.78, 0.60 (odd), against 0.84, 0.68,
0.32, 0.68 at rest.

* **Warm matter now pulls distant matter harder than cold matter, at every distance.** Single matter at k = 8:
  +2.36, +2.12, +1.75 and +1.07 × 10⁻⁴ at r = 6, 9, 13.5 and 20, against +1.61, +1.14, +0.41 and +0.53 at rest, with
  k = 2 in between. Averaged over the receivers the net pull is ×1.98 of cold at k = 8 (×1.53 at k = 2), and the
  source wave's pull on the receivers' quiet channel ×2.69 (×2.14 for odd matter), where the two-way wave gave ×0.85
  with round 16's mixing form and ×1.22 with the velocity's parity (§27.2).
* **The pull follows the square root of the wave's extra strength.** R_net averages 1.06 over the four shells (single,
  k = 8; 1.12 at k = 2) and 0.94 (odd, k = 8; 1.14 at k = 2): what the law needs is 1. Shell by shell it scatters
  (0.59–1.83), and the scatter follows the cold reference: at r = 13.5 the cold source's receivers keep step at only
  0.32, so the warm-to-cold ratio there is large (1.65). From the quiet channel alone R is 1.26–1.39 on average; the
  difference is the receivers' own radiators, which push back in proportion to the wave's intensity (they take 40% of
  the quiet channel's pull at rest and 52% at k = 8, with this matter's f = 0.2).
* **Heat now helps the beat.** The warm source's pieces share one rhythm (spread 4.6 × 10⁻⁵ at k = 8, against
  1.7 × 10⁻⁴ at rest and 2.6 × 10⁻³ with the two-way wave, §27.1), and distant matter keeps step at 0.93 on average
  (0.63 at rest).
* **Colliding sources pull like cold ones.** At ν = 50 the glow falls from ×4.75 to ×1.62 of rest and the net pull is
  +1.53, +0.99, +0.38 and +0.50: the cold pull, within 5–13%. The little extra glow that survives the collisions adds
  no pull (R_net ≈ 0.67), which is the law's rule for gas (k = 0).
* **A warm source brought to rest returns to the cold pull.** Stopped at t = 4,000, its glow returns to ×1.01 and its
  pull to +1.85, +1.06, +0.63 and +0.11 (single) or +1.67, +1.15, +0.57 and +0.76 (odd), the cold values within the
  scatter except at the outermost shell of single matter.
* **The velocity's parity is not needed here:** with the one-way wave, round 16's mixing form does as well as the
  velocity's (R_net 1.06 against 0.94 at k = 8). The parity still matters for the two-way part of any real wave.

**Matter whose radiators mostly ring inside** (`--set oneway_lowf` → `oneway_lowf.json`, `oneway_lowf_summary.json`):
round 16's follow-up matter (§26.4: the radiators send a twentieth of their energy into the companion, f = 0.05), with
the velocity's parity, at the same four distances, two arrangements each.

**One-way wave, radiators that mostly ring inside (f = 0.05, with the velocity's parity)**

| source | glow ÷ rest | keeping step (all receivers) | net pull at r = 6, 9, 13.5, 20 (10⁻⁴) | R_net(σ, r) at 6, 9, 13.5, 20 | R from the quiet channel |
|---|---|---|---|---|---|
| at rest | 1.00 | 0.65 | +2.14, +1.28, +0.49, +0.50 | 1 | 1 |
| free, k = 2 | 2.01 | 0.80 | +2.64, +1.77, +1.46, +0.73 | 0.94, 1.02, 1.52, 1.01 | 0.96, 1.04, 1.80, 0.96 |
| free, k = 8 | 4.76 | 0.73 | +2.85, +1.84, +1.98, +0.77 | 0.81, 0.85, 1.67, 0.67 | 0.85, 0.92, 1.92, 0.82 |

* The receivers' own radiators now take 8% of the quiet channel's pull at rest and 11% at k = 8, as the rule 2f = 10%
  says (§26.4), against 40–52% at f = 0.2.
* R_net averages 1.12 at k = 2 and 1.00 at k = 8, and the net pull is ×1.50 and ×1.69 of cold. The sources' rhythm
  spread at k = 8 is 1.9 × 10⁻⁵ (1.8 × 10⁻⁴ at rest); the books balance to 5 × 10⁻¹².

**The same distances with the two-way wave** (`--set distance` → `distance.json`, `distance_summary.json`; the
baseline for R_net: the same matter, arrangements and times, round 16's two-way coupling):

**Two-way wave, round 16's matter (single)**

| source | glow ÷ rest | keeping step (all receivers) | net pull at r = 6, 9, 13.5, 20 (10⁻⁴) | R_net(σ, r) at 6, 9, 13.5, 20 | R from the quiet channel |
|---|---|---|---|---|---|
| at rest | 1.00 | 0.30 | +0.75, +0.71, +0.38, +0.22 | 1 | 1 |
| free, k = 2 | 1.95 | 0.29 | +0.57, +0.87, +0.34, +0.00 | 0.57, 1.02, 0.73, 0.01 | 0.89, 1.15, 0.72, 1.04 |
| free, k = 8 | 4.69 | 0.11 | -0.89, -0.55, -0.03, -0.09 | -0.53, -0.42, -0.04, -0.24 | 0.49, 0.33, 0.30, 0.27 |

**Two-way wave, with the velocity's parity (single, odd)**

| source | glow ÷ rest | keeping step (all receivers) | net pull at r = 6, 9, 13.5, 20 (10⁻⁴) | R_net(σ, r) at 6, 9, 13.5, 20 | R from the quiet channel |
|---|---|---|---|---|---|
| at rest | 1.00 | 0.30 | +0.75, +0.71, +0.38, +0.22 | 1 | 1 |
| free, k = 2 | 1.92 | 0.28 | +0.42, +0.41, +0.35, +0.32 | 0.37, 0.48, 0.73, 1.07 | 0.79, 0.94, 0.67, 2.70 |
| free, k = 8 | 4.65 | 0.16 | -0.64, -0.23, -0.03, -0.02 | -0.39, -0.18, -0.04, -0.05 | 0.61, 0.55, 0.31, 1.33 |

* **With the two-way wave, a warm source pushes the nearer receivers instead of pulling them:** −0.89 × 10⁻⁴ at r = 6
  and −0.55 at r = 9 (k = 8, single), where the same source with the one-way wave pulls +2.36 and +2.12. Distant matter
  keeps step at only 0.11–0.16 (0.30 at rest); the receivers' own radiators, pushed in proportion to the wave's
  intensity, then outweigh their quiet channel (−1.30 against +0.91 × 10⁻⁴ averaged, k = 8). The sources' rhythm spread
  is 4.1 × 10⁻³ at k = 8, ninety times the one-way wave's.
* **The cold pull is weaker too:** +0.75, +0.71, +0.38 and +0.22 × 10⁻⁴, against +1.61, +1.14, +0.41 and +0.53 with the
  one-way wave, because distant matter keeps step less well (0.30 against 0.63).
* So in the full model the one-way wave turns a warm source's push into a pull, and doubles the cold pull at the inner
  shells.

**Mass, with the two-way wave** (`--set mass` → `mass.json`, `mass_summary.json`: sources of 24, 48 and 96 pieces at
the same density, the velocity's parity, receivers at r = 12 and 18, two arrangements, 8,000 time units):

| source pieces | glow per piece at rest | glow at k = 8 ÷ rest | keeping step, rest / k = 8 | net pull at r = 12, 18, rest (10⁻⁴) | at k = 8 |
|---|---|---|---|---|---|
| 24 | 2.74 × 10⁻³ | 5.24 | 0.30 / 0.05 | +0.20, +0.17 | -0.27, +0.18 |
| 48 | 2.78 × 10⁻³ | 4.75 | 0.17 / 0.36 | +0.21, -0.04 | +0.42, +0.55 |
| 96 | 2.88 × 10⁻³ | 4.72 | 0.40 / 0.09 | +0.69, +0.49 | -0.77, -0.11 |

* **The glow is in proportion to mass in the full model:** 2.74, 2.78 and 2.88 × 10⁻³ per piece at rest for 24, 48 and
  96 pieces (within 5%), and heat multiplies it by about the same factor at every size (×5.24, 4.75 and 4.72 at k = 8).
* With the two-way wave the pull is erratic at every size (keeping step 0.05–0.40), as §27.1's damage predicts.

**Mass, with the one-way wave** (`--set oneway_mass` → `oneway_mass.json`, `oneway_mass_summary.json`: sources of
24, 48 and 96 pieces at the same density, the f = 0.05 matter with the velocity's parity, receivers at r = 12 and 18,
two arrangements, 8,000 time units):

| source pieces | glow per piece at rest | glow at k = 8 ÷ rest | keeping step, rest / k = 8 | net pull at r = 12, 18, rest (10⁻⁴) | at k = 8 | R_net at r = 12, 18 |
|---|---|---|---|---|---|---|
| 24 | 2.19 × 10⁻³ | 5.94 | 0.62 / 0.63 | +1.07, +0.30 | +1.89, +1.00 | 0.88, 1.40 |
| 48 | 2.30 × 10⁻³ | 5.22 | 0.46 / 0.71 | +1.03, +0.69 | +3.11, +1.65 | 1.42, 1.16 |
| 96 | 2.37 × 10⁻³ | 5.46 | 0.61 / 0.88 | +1.13, +1.21 | +3.32, +2.79 | 1.51, 1.06 |

* **The glow is in proportion to mass** (2.19, 2.30 and 2.37 × 10⁻³ per piece at rest, within 8%), with about the same
  heat multiplier at every size (×5.2–5.9 at k = 8).
* **The pull grows as the square root of the mass, as the law needs.** Four times the mass (24 → 96 pieces) raises the
  net pull averaged over the receivers ×1.71 at rest and ×2.11 at k = 8, where the law's √(G M a)/r gives ×2. Shell by
  shell the doublings scatter (×0.97–2.28), following how well each shell keeps step (0.42–0.94) in 8,000 time units.
* **At every size the warm source pulls harder in proportion to the square root of its extra wave:** R_net averages
  1.14, 1.29 and 1.28 for 24, 48 and 96 pieces.

**What this establishes, and what not.** In the full, energy-booked model, with one kind of matter for sources and
receivers, a companion wave that only travels outward lets motion raise the output and the pull together, with the net
pull growing as the square root of the wave's extra strength (R ≈ 1) at four distances and as the square root of the
source's mass; collisions and stopping switch the extra pull off, and with a two-way wave the same warm source pushes. The one-way coupling is still imposed by hand: whether the companion's wave is one-way is the
question the next round has to answer from the companion as a flowing medium.

### 27.5 Where round 17 leaves the proof

| link | status after round 17 |
|---|---|
| which term puts warm matter out of tune | **found:** the non-gradient part of the two-way tugs between warm pieces, through their stirred radiators. With it switched off, distant matter keeps step at 0.88 with a k = 8 source (cold: 0.63) (§27.1) |
| a protection inside each piece | **none can remove the tugs:** at equal glow reciprocity fixes their size (twelve structures, ±4%), and a coupling linear in the velocity cannot be one-way. The review's pair protects only when its halves radiate through separate channels, not in a shared wave. What helps: the velocity's time-reversal parity (a first-principles correction) and a widely detuned pair (a poor scatterer): 0.18 → 0.41 at k = 8 (§27.2) |
| the motion's mixing | **corrected:** odd under time reversal, as a velocity's coupling must be. Full model, k = 8: keeping step 0.23 → 0.33, the warm wave's pull on the receivers' quiet channel ×0.85 → ×1.22 of cold (§27.2) |
| a companion wave that only travels outward | **removes the damage in the reduced model:** one beat, R(σ, r) = 1.0–1.5 at four distances, power ∝ mass; needs an inward leak ≲ 10% throughout the source (§27.3). **In the full, energy-booked model it does the same:** warm matter pulls distant matter harder than cold matter at all four distances (net pull ×1.98 of cold at k = 8), R_net = 1.06 on average (0.94 with the velocity's parity), and colliding and stopped sources pull like cold ones (§27.4). The one-way coupling is imposed, not yet derived |
| the frequency requirement | a warm source's common beat must stay where cold matter's is, or far receivers cannot follow (§27.2); the one-way wave meets it (within 3 × 10⁻⁵) |
| the net pull | the receivers' radiators are pushed in proportion to the wave's intensity; radiators that mostly ring inside (f ≪ 1/2) are needed (round 16); with the one-way wave and radiators that mostly ring inside (f = 0.05) they take only 8–11% of the pull, and R_net = 1.00 on average at k = 8 (§27.4) |
| power in proportion to mass | **holds in the full model:** the glow per piece is the same within 8% for 24, 48 and 96 pieces at the same density, with either wave, and heat multiplies it by about the same factor at every size. With the one-way wave the pull at a fixed distance grows as the square root of the mass: ×1.71 at rest and ×2.11 at k = 8 for four times the mass, where the law gives ×2 (§27.4) |
| every kind of matter falls alike | not tested: the model does not yet say what a piece's inertia is. In it, a locked receiver's pull is (1 − 2f) × its fed power over the wave's speed, so equal acceleration needs every kind of matter to share f (or f ≪ 1) and an inertia in proportion to its quiet amplitude — requirements for the microscopic model |
| the law and the suite | unchanged (59 pass, 11 close, 7 fail); no astrophysical fit was changed this round, as the review asked |

**Next:**
1. **The companion as a flowing medium.** Derive the companion's wave in its own outflow: the crests' speed against the
   stream (round 16's galaxy bound, crests at no more than half the travel speed, would put it between u/2 and u if the
   crests are the stream's slow wave), whether the companion inside a source streams as one medium, and the wave's
   exchange of energy and momentum with the stream. That decides whether the wave is one-way, and it gives the
   coupling to build into the full model in place of the imposed one.
2. **Then Step C with that wave:** R(σ, r) at more distances, cold, free, colliding and stopped sources, mass scaling
   (physical mass and refinement separately), receivers of different kinds, with the energy and momentum the stream
   takes booked explicitly.
3. **What inertia is** for a piece of this matter, so that "every kind of matter falls alike" can be tested.
4. A physical quiet store with a gap (χ and γ, and u from them); the release factor from drained receivers; the
   Doppler shifts of real motion.
5. Only then the astrophysical fits.

## 28. Round 18, 25 September 2026: the review's checklist, and the companion as a medium

An independent review of the fresh write-up and of round 17 (25 September 2026) set out ten items between this work
and a paper, in dependency order: (1) derive the outward-only wave from a local medium, with the medium's own energy
and momentum; (2) recover the law's numbers, not only its trends, including a hot-shell benchmark; (3) define inertia
and show why different matter falls alike; (4) screening, the external field and light's coupling, derived or stated
as assumptions (with Cassini's 2026 value); (5) collisions as a causal calculation; (6) a statistical inference in
place of the scoreboard; (7) the remaining discrepancies explained or delimited; (8) a literature search and one
bounded claim; (9) a frozen, reproducible release; (10) a manuscript written to the strength of the evidence. Its
shortest path: derive the medium and its one-way response first, then force scaling and inertia with it, then freeze
the law and validate it statistically. This round does the write-up items at once and starts item 1.

### 28.1 The write-up matched to the evidence

`BLOG.md` and the page: the subtitle says "a few stated assumptions, worked out step by step" instead of "from first
principles"; a new §3.10 lists every piece of the law as assumed, derived (and where), borrowed or open; the scoreboard
and the comparison table carry a statement that they are typical misses, not a statistical comparison, with the
uncertain inputs uncounted; the regression suite is described as an engineering tool (a pass is within two standard
deviations of one measurement; the checks are not independent); "Einstein's, exactly" becomes "no measurable extra
pull, with Einstein's theory assumed in strong fields"; strong-lens agreement is stated with light's response assumed;
"no structure inside a piece can remove it" becomes "none of the twelve we tried, and reciprocity suggests why for
structures that send and receive through one channel"; "exactly what the law requires" becomes the measured values
with their spread (R_net 1.06 on average, 0.6–1.8 by distance; ×1.7 and ×2.1 for four times the mass, two
arrangements); the collision rule is conditional (the lensing follows the companion's net flow; it can sit on dense
gas among converging hot galaxies, as in Abell 520's clump, and returns around long-stopped gas); the "explains why"
row is gone; §9 is the review's checklist.

### 28.2 Cassini, 2026

`code/cassini_2026_v18.py` → `run-cassini-2026-v18/cassini_2026_v18.json`. Park, Hees, Famaey, Desmond & Durakovic
(Phys. Rev. D, 28 July 2026; arXiv:2602.17884) re-estimated the Galaxy's quadrupole in the Sun's field with the DE440
data, simultaneously with the other ephemeris parameters: Q₂ = (1.6 ± 1.8) × 10⁻²⁷ s⁻², 40% tighter than Hees et al.
2014's (3 ± 3) × 10⁻²⁷, which set the release length L in round 9. The adopted law (round-12 constants, the Galaxy's
Newtonian pull at the Sun 1.638 × 10⁻¹⁰ m/s², the suite's formula):

| release length L | Q₂ (10⁻²⁷ s⁻²) | from the 2026 value | wide binaries: extra pull at 7,000 / 20,000 AU |
|---|---:|---:|---|
| 0 (released at once) | 29.4 | +15.4σ | 15.4% / 15.6% |
| 0.10 pc | 6.13 | +2.5σ | 4.4% / 9.7% |
| **0.145 pc (adopted)** | **4.38** | **+1.54σ** | **3.2% / 7.6%** |
| 0.20 pc | 3.26 | +0.9σ | 2.4% / 6.0% |
| 0.30 pc | 2.22 | +0.3σ | 1.6% / 4.3% |
| 0.50 pc | 1.35 | −0.1σ | 1.0% / 2.7% |
| 1.0 pc | 0.69 | −0.5σ | 0.5% / 1.4% |

The adopted law sits 1.54σ above the new value (it was 0.46σ from the old one). A release length of at least 0.19 pc
meets 1σ, at least 0.12 pc 2σ; the wide-binary forecast falls with it. The suite's check now uses the 2026 value (still
a pass under its two-sigma rule, flagged "worse"). L is not refitted here: it is one of the pieces to be derived
(review item 4), and the law is to be frozen only after item 1.

### 28.3 Why the model's glow grows more slowly than 1 + k

`code/heat_factor_v18.py` → `run-heat-factor-v18/heat_factor_v18.json`. The heat weight k is defined for one piece in
linear response (the mixing q is chosen so that one piece's released glow is k times its cold glow: k = 12 q² for round
16's matter, exactly linear). In the full model at k = 8 a source of 48 pieces glows ×4.75 of its cold output, not ×9.
The same matter, 48 pieces in a dense ball (radius 3 wavelengths, mean spacing about 1.3 wavelengths) and in a
dilute one (radius 30, mean spacing about 13 wavelengths; both labels corrected 25 September, the models use
wavelength 1), and one isolated piece, 4,000 time units:

| glow ÷ cold at k = 2, 8, 16 (1 + k = 3, 9, 17) | round 16's matter (f = 0.2) | f = 0.05 matter |
|---|---|---|
| dilute ball: pieces nearly independent | 2.60, 6.97, 11.8 (76–87%) | 2.62, 7.23, 12.9 (76–87%) |
| dense ball | 1.81, 4.54, 7.57 (45–60%) | 1.82, 4.40, 7.63 (45–61%) |
| one isolated piece (a single velocity draw) | 1.54, 3.14, 5.19 | 1.55, 3.17, 5.30 |
| the pieces' inversion at k = 16 ÷ cold (dilute) | 0.76 | 0.88 |

* **Two causes.** Nearly independent pieces reach 76–87% of 1 + k: the extra glow draws on each piece's supply, and its
  inversion falls (by 14% at k = 8 and 24% at k = 16 for round 16's matter), while its quiet amplitude hardly changes
  (−3%). In a dense ball the shared wave holds back a further part (45–61% of 1 + k), round 16's "renormalization"
  of dense clouds (§26.3).
* A single isolated piece depends on its one random velocity (its |δ|² is one draw of three squared normal numbers), so
  it is not the ensemble value.
* **In the law** this reads as an effective u larger than the microscopic one, by 1/√0.8 in dilute matter and 1/√0.5
  in dense; the fitted u absorbs a common factor, but a factor that depends on packing would make the heat term differ
  between ellipticals and clusters. Which way real matter goes is part of review item 2.

### 28.4 The companion as a flowing medium, exactly in one dimension

`code/flowing_medium_v18.py` → `run-flowing-medium-v18/flowing_medium_v18.json` (units c = ω = 1). Round 17 made the
wave one-way by deleting the inward couplings and booked the energy "to the wave and its flow". Here the wave gets its
own local equations, in two versions, and everything the review asked for follows from them.

**A. Waves carried by the stream.** The companion streams at U and its waves move at c relative to it:
```
(∂_t + U ∂_x)² ψ − c² ∂_x² ψ = Σ_j Q_j(t) δ(x − x_j)
```
(the one-dimensional analogue-acoustics equation; its three-dimensional form, for an irrotational stream of density
ρ₀, is Unruh's). The pieces enter through L_int = Σ Q_j ψ(x_j): piece j feeds the medium the power Q_j ∂_tψ(x_j) and
feels the force Q_j ∂_xψ(x_j). The steady medium conserves the pseudo-energy H = π²/2 − Uπψ_x + c²ψ_x²/2 (π = ψ_t + Uψ_x)
and the pseudo-momentum.

* **Speeds, with their frames.** A disturbance moves at U ± c in the lab (c relative to the stream). For U > c both
  move downstream: nothing travels upstream. In the time-harmonic problem the two waves have lab wavenumbers
  k_f = ω/(U + c) and k_s = ω/(U − c); their crests move at U + c and U − c, and so does their energy (no dispersion).
  In the stream's frame the slow wave's frequency is negative (ω' = −ωc/(U − c)), so it carries negative lab
  pseudo-energy and pseudo-momentum pointing upstream.
* **The Green's function is exact:** for U > c, g(x) = (i/2ωc)(e^{ik_f x} − e^{ik_s x}) downstream and 0 upstream; for
  U < c, (i/2ωc) e^{ik_± x} on either side.
* **One-way is derived:** for U > c no piece hears any piece downstream of it. In three dimensions, for a stream that
  moves at u everywhere along the field lines of g_N (the guided stream of rounds 12–13, whose density falls as |g_N|
  toward a source's centre while its speed stays u), every characteristic moves outward (dr/dt ≥ u − c > 0), so the
  domain of influence of every point lies strictly outward, right into the centre. A pressure-driven wind would not do:
  its speed grows from zero at the centre and is below its own wave speed inside the source (the Chevalier–Clegg
  solution), so waves travel both ways exactly where the warm pieces are.
* **A lone piece in the stream:**

| U/c | one-way | pseudo-energy fed, per \|q\|² | push downstream, per \|q\|² | power in the stream's frame | pull on a locked receiver, per \|q_R\|\|φ\| | power fed ÷ pull |
|---|---|---:|---:|---:|---:|---|
| 0 | no | 0.250 | 0 | 0.250 | 0.500 | c |
| 0.5 | no | 0.250 | 0.167 | 0.333 | 0.333 | 3.0 U |
| 1.1 | yes | 0 | 1.190 | 1.310 | 2.619 | 0.17 U |
| 1.2 | yes | 0 | 0.568 | 0.682 | 1.364 | 0.31 U |
| √2 | yes | 0 | 0.250 | 0.354 | 0.707 | 0.50 U |
| 2 | yes | 0 | 0.083 | 0.167 | 0.333 | 0.75 U |

* **The energy bill.** A receiver that leads the local wave by a quarter cycle feeds P = (ω/2)|q_R||φ| and is pulled by
  F = (k_f + k_s)|q_R||φ|/4 toward the source, exactly (the two waves' beating cancels in the ratio). Both waves pull.
  So P/F = 2ω/(k_f + k_s) = (U² − c²)/U, not a crest speed. The galaxies' bound on the feeding feedback (round 16,
  P/F ≤ u/2) then needs c ≥ U/√2; with one-way (c < U) the window is **0.71 ≤ c/U < 1**. (The blog's earlier hint, "c
  between U/2 and U", used the slow wave's crest speed alone.)
* **But the stream pushes every emitter downstream.** A lone piece in a supersonic stream feeds no net pseudo-energy
  (its two waves' pseudo-energies cancel) but feels a downstream push |q|²/(4(U² − c²)), exactly its power in the
  stream's frame divided by U: wave drag, like wind on a sail. The push is independent of any other source, while the
  pull on a receiver falls with the incoming wave: their ratio is |q_R|/(2(U/c)|φ|). Per kilogram, with the emission
  equal to the cold glow ℓ, the push is the stream-frame power over u: (U/c)/((U/c)² − 1) × ℓ/u in one dimension,
  at least 0.71 a in the window. Three dimensions change the factor, not its order (a source moving through its medium
  radiates at least its static power, so the push is at least of order ℓ/u = a/2): larger than the law's own pull
  √(a g_N) through the outskirts of galaxies. **As it stands this medium is ruled out.** (Subsonic streams push too,
  and are two-way.)

**B. Waves the stream absorbs when they move against it.** Waves travel at c through the matter's frame; the stream
absorbs the part that moves against it, amplitude e^{−κ} per unit length travelled inward:
```
(∂_t − c ∂_x) u₊ = S − c κ_L(x) u₊,    (∂_t + c ∂_x) u₋ = S − c κ_R(x) u₋,    ψ_t = (u₊ + u₋)/2,  ψ_x = (u₊ − u₋)/2c
```
(the undamped equations are ψ_tt − c²ψ_xx = S; an outflow from x = 0 damps left-movers at x > 0 and right-movers at
x < 0). Emission is symmetric, so a lone piece feels **no push**; the pull and the energy bill are the static ones (P/F =
c, so the galaxies need c ≤ u/2); a wave from x = 2 reaches x = 1 at 0.37 of the outward coupling for κ = 1, 0.10 for
κ = 2.3 and 0.007 for κ = 5. The absorbed energy and momentum go to the stream.

**Checks against a direct simulation of the local equations** (grid 0.01, RK4, pieces as prescribed oscillators; the
receiver leads the source's wave by a quarter cycle; force and power as time averages of Q ∂_xψ and Q ∂_tψ):

| case | force on source / receiver: simulated | exact | power fed, source / receiver: simulated | exact |
|---|---|---|---|---|
| A, U = 0.5 | +0.570 / −0.005 | +0.574 / 0.000 | 0.455 / 0.499 | 0.454 / 0.500 |
| A, U = 1.3 | +0.360 / −0.504 | +0.362 / −0.499 | 0.001 / 0.460 | 0 / 0.457 |
| B, κ = 0 | −0.177 / −0.251 | −0.179 / −0.250 | 0.073 / 0.502 | 0.071 / 0.500 |
| B, κ = 2 | −0.008 / −0.259 | 0 / −0.250 | 0.242 / 0.494 | 0.250 / 0.500 |

In A at U = 0.5 the receiver's pull and its drag happen to cancel; at U = 1.3 the source, which nothing upstream
reaches, feels only its drag (+0.36).

### 28.5 The absorbing stream in the full model

`code/absorbing_stream_v18.py` (sets `distance`, `strong`) → `run-absorbing-stream-v18/distance.json`, `strong.json`
and their summaries; the option `absorb` in `code/one_matter_v17.py`; `code/absorb_reduced_v18.py` →
`reduced_scan.json`.

**The medium in the full model.** Round 16's model of one kind of matter (48 source pieces in a ball of radius 3,
receivers on shells at r = 6, 9, 13.5 and 20, every piece coupled to every other through the complete wave), with the
stream of §28.4 B: a wave travelling in direction n through a point where the stream moves along e = x/|x| loses
amplitude at the rate κ max(0, −n·e). Along the straight path from piece l to piece j this integrates exactly to
T_jl = exp(−κ L_jl), L_jl = r_l minus the smallest radius the path reaches before it turns outward or ends: κ is the
attenuation per unit of distance travelled inward. Every block of the coupling between two pieces is multiplied by
T_jl, the force carries its derivative with respect to the receiving piece's position, and each piece's own terms are
unchanged (emission stays symmetric: no drag). The energy the pieces give the wave, Im(z†Mz), is split into what
reaches infinity (every piece's far-field pattern, each direction attenuated by its own inward travel along the ray)
and what the stream absorbs. Round 16's matter (f = 0.2), two arrangements, 16,000 time units, the second half measured.

| wave | source | net pull at r = 6, 9, 13.5, 20 (10⁻⁴) | keeping step | sources' rhythm spread | R_net at 6, 9, 13.5, 20 | absorbed |
|---|---|---|---|---|---|---|
| two-way (§27.4) | cold | +0.75, +0.71, +0.38, +0.22 | 0.30 | 3.3 × 10⁻⁴ | 1 | – |
| | k = 8 | −0.89, −0.55, −0.03, −0.09 | 0.11 | 4.1 × 10⁻³ | pushed | – |
| absorbing, κ = 2.3 | cold | +1.54, +0.95, +0.62, +0.37 | 0.64 | 3.1 × 10⁻⁴ | 1 | 33% |
| | k = 8 | +0.68, +1.51, +0.72, +0.60 | 0.50 | 3.3 × 10⁻³ | 0.21, 0.87, 0.51, 0.88 | 38% |
| | colliding, ν = 50 | +1.29, +0.77, +0.52, +0.35 | 0.46 | 3.3 × 10⁻⁴ | 0.74, 0.72, 0.75, 0.79 | 33% |
| absorbing, κ = 5 | cold | +1.36, +0.76, +0.56, +0.47 | 0.72 | 2.5 × 10⁻⁴ | 1 | 37% |
| | k = 2 | +1.42, +1.33, +1.04, +0.27 | 0.70 | 9.0 × 10⁻⁴ | 0.86, 0.95, 1.00, 0.56 | 43% |
| | k = 8 | +1.72, +1.38, +1.23, +0.60 | 0.67 | 3.0 × 10⁻³ | 0.54, 0.81, 0.82, 0.75 | 44% |
| | colliding, ν = 50 | +1.22, +0.72, +0.55, +0.44 | 0.51 | 2.5 × 10⁻⁴ | 0.72, 0.70, 0.74, 0.76 | 37% |
| absorbing, κ = 20 (not passive: see below) | cold | +1.25, +0.94, +0.44, +0.50 | 0.75 | 1.0 × 10⁻⁴ | 1 | 42% |
| | k = 2 | +1.02, +1.14, +0.81, +0.44 | 0.65 | 7.2 × 10⁻⁴ | 0.66, 0.97, 1.36, 0.71 | 47% |
| | k = 8 | +1.63, +1.60, +1.26, +0.56 | 0.70 | 1.9 × 10⁻³ | 0.59, 0.89, 1.25, 0.73 | 48% |
| | colliding, ν = 50 | +1.15, +0.89, +0.29, +0.49 | 0.52 | 1.1 × 10⁻⁴ | 0.71, 0.79, 0.55, 0.77 | 42% |
| imposed radius mask (§27.4) | cold | +1.61, +1.14, +0.41, +0.53 | 0.63 | 1.7 × 10⁻⁴ | 1 | – |
| | k = 8 | +2.36, +2.12, +1.75, +1.07 | 0.93 | 4.6 × 10⁻⁵ | 0.68, 1.01, 1.65, 0.89 | – |

* **No push, and the pull grows with heat.** With the absorbing stream a warm source pulls distant matter at every
  distance, where the two-way wave pushed it; the gain over cold comes closer to the law's square root as the
  absorption strengthens: R_net averages 0.62 at κ = 2.3 and 0.73 at κ = 5 (k = 8; 0.81 and 0.84 at k = 2), and 0.87
  (0.93) at κ = 20, where, however, the ray form of the medium is no longer passive (below). The imposed mask gave
  1.06.
* **Colliding sources pull like cold ones** at every κ (R_net 0.7–0.8: their small surviving extra glow adds no pull).
* **Every watt booked.** The pieces' own budgets close to 2 × 10⁻¹² (3 × 10⁻⁷ with collisions); the stream absorbs
  33–48% of what they give the wave (the part emitted inward), and absorbs more from warm sources.
* **Passivity.** The dissipative part of a coupling, (M − M†)/2i, must have no negative eigenvalue, or some state of the
  pieces could draw energy from the medium. For the same arrangement: the static two-way coupling, smallest eigenvalue
  +9 × 10⁻⁶ of a largest 5.59 (passive); round 17's imposed mask −0.87 of 4.90 (not passive: the hand-made coupling
  could create energy); the absorbing stream +3 × 10⁻⁴, −8 × 10⁻⁴ and −8 × 10⁻³ of about 4.8 at κ = 1, 2.3 and 5
  (passive to within the ray approximation's accuracy), but −0.09 at κ = 10 and −0.46 at κ = 20: an absorption length (1/κ = 0.05) far
  below the wavelength (1; corrected 25 September from "2π") and the pieces' spacing (about 1.3) is outside the ray form's range, so the κ = 20 run is an
  extrapolation of the approximation, not a physical medium (κ = 10 is borderline); strong absorption needs a full
  wave treatment.
* **But a warm source keeps no single beat** (rhythm spread 2–3 × 10⁻³, against 4.6 × 10⁻⁵ for the imposed mask). The
  reduced model (three arrangements, the same geometry) shows why:

| medium | sources' rhythm spread, k = 8 | keeping step at r = 6, 9, 13.5, 20: k = 8 (cold) | pairs still two-way |
|---|---:|---|---:|
| two-way | 4.5 × 10⁻³ | 0.19, 0.23, 0.10, 0.09 (0.52, 0.29, 0.28, 0.48) | 100% |
| imposed radius mask | 6.3 × 10⁻⁵ | 0.98, 0.96, 0.85, 0.81 (0.88, 0.79, 0.59, 0.76) | 0 |
| inward absorption, κ = 2.3 | 4.1 × 10⁻³ | 0.51, 0.53, 0.42, 0.39 (0.46, 0.45, 0.27, 0.25) | 36% |
| κ = 5 | 3.6 × 10⁻³ | 0.74, 0.63, 0.60, 0.53 (0.75, 0.59, 0.43, 0.42) | 17% |
| κ = 10 (ray form not passive) | 3.0 × 10⁻³ | 0.81, 0.70, 0.81, 0.45 (0.82, 0.59, 0.47, 0.36) | 8% |
| κ = 20 (not passive) | 2.3 × 10⁻³ | 0.91, 0.88, 0.83, 0.61 (0.92, 0.82, 0.60, 0.52) | 3% |
| inward and sideways absorption ("kinetic"), κ = 5 | 1.5 × 10⁻³ | 0.59, 0.33, 0.55, 0.21 (0.74, 0.17, 0.45, −0.02) | 15% |
| κ = 10 | 8.2 × 10⁻⁴ | 0.58, 0.46, 0.21, 0.15 (0.45, 0.32, 0.46, 0.16) | 6% |

  ("Pairs still two-way": the sources' coupling, summed over pairs, that the weaker direction keeps.) Stronger inward
  absorption lets distant matter keep step, better with a warm source than a cold one, but even with 97% of the
  two-way coupling gone the source keeps no single beat. **The reason is what a physical one-way medium lets a piece
  hear.** In the absorbing stream a wave that must pass the centre to get from one side of the source to the other is
  absorbed, so for strong absorption a point hears only the matter inside the ball whose diameter joins it to the
  centre, and a distant receiver mainly the near half of a source. Round 17's rule let every piece hear every piece
  nearer the centre, so the innermost pieces set one beat for all; in the absorbing stream the two sides follow
  their own inner pieces. A stream that also absorbs sideways waves tightens the beat but cuts distant matter off from
  most of the source.
* **What this changes.** Two of round 17's results were properties of the imposed rule rather than of any medium found
  so far: one beat for a warm source, and "a piece feels exactly the matter nearer the centre, Newton's rule for a
  round source". With the physical media the heard region is smaller, which is the review's hot-shell benchmark in
  general form: the law's |g_N| and S count all matter, a one-way medium does not. The pull's growth with heat does
  survive in a physical, drag-free, nearly passive medium: 73–84% of the square root the law needs where the medium's
  ray form is passive (κ ≤ 5). *(Correction, 25 September: R_net divides by the model's own intensity gain, not by
  √(1 + k). The model's glow grows only ×4.75 at k = 8, so against the law's heat gain the raw pull ratios at k = 8
  (κ = 5: 1.26, 1.82, 2.20, 1.28 by shell) are about 1.6 against the law's 3, about half; §30.)*

### 28.6 Where round 18 leaves the proof, against the review's checklist

| review item | status after round 18 |
|---|---|
| 1. the one-way wave from a local medium | **started.** Waves carried by the stream: one-way derived (exactly in 1D; by characteristics in 3D for a stream at u everywhere, into the centre), energy bill window 0.71 ≤ c/u < 1, but wave drag of order the emitted power over u, larger than the law's pull in galaxies' outskirts: ruled out as it stands. Waves the stream absorbs when they move against it: one-way with no drag, passive to within the ray approximation, and in the full 3D model removes the push, lets colliding sources pull like cold ones and makes the warm source's pull grow with heat at 73–84% of the square root (κ ≤ 5, where its ray form is passive); but a warm source keeps no single beat, because each piece hears only the inner matter on its own side (§28.5). Still to do: the absorption from the stream's own microphysics (what sets κ), the wave's momentum given to the stream followed as the stream's evolution, and where the streams of two bodies meet |
| 2. the law's numbers from the models | the heat factor explained (76–87% of 1 + k for nearly independent pieces, 45–61% in dense balls; §28.3); the hot-shell benchmark stated: in the absorbing medium a system inside a shell of hot matter hears nothing from the shell, while the law's S counts it, so the medium predicts a one-way S whose effect on the cluster fits is the next quantitative test; the mass exponent with its uncertainty still to be fitted |
| 3. inertia and universal free fall | open |
| 4. screening, light, Cassini | Cassini updated to the 2026 value (1.54σ; L ≥ 0.19 pc for 1σ); the release factor and light's response still assumed, now stated as such |
| 5. collisions as a calculation | open |
| 6. statistical inference | open; the write-up no longer presents the suite's tally or the parameter count as evidence of probability |
| 7. remaining discrepancies | listed in BLOG §9, each with the measurement that would settle it |
| 8. novelty and scope | literature search owed; the review's suggested first claim ("motion-enhanced attraction in an active streaming medium") adopted as the working scope |
| 9. reproducible release | open (figures drawn with earlier fits are marked) |
| 10. the manuscript | the write-up's claims matched to the evidence (§28.1) |

**Next:**
1. **Strong absorption as a wave problem,** not in the ray form: the absorbing stream's full wave equation (a
   direction-dependent loss in the medium), its Green's function between pieces closer than a wavelength, and whether
   a passive medium with strong absorption reaches the law's square root and a single beat.
2. **What sets κ:** the absorption of counter-moving waves from the stream's own degrees of freedom (a two-component kinetic model: wave quanta scattering off the streaming companion, with the rate set by the relative speed), energy and momentum followed in both.
3. **The one-way law it implies,** tested before any fit is changed: S and the pull's direction with the absorption factor along each path (the hot-shell benchmark), first against the X-COP cluster profiles and the collision maps, with every other constant fixed.
4. **Mass scaling and distance with the absorbing medium:** the exponent with its uncertainty over a wider range of mass, physical mass and resolution varied separately.
5. **Inertia**, then the statistical programme (review items 3, 6, 9).

## 29. Round 19, 25 September 2026: the medium tested against the data, and the Milky Way refitted

A second review (25 September 2026) set the order of work: finish the medium and the law it produces (the absorbing
medium in full wave form: what absorbs, energy and momentum, convergence, passivity, no self-push; then the law's
dependence on mass, distance, dispersion, density and shape, with uncertainties, the hot-shell benchmark and the
absorbed power, tested on clusters before anything is refitted); improve the joint Milky Way and Solar-System
analysis alongside (vary the Galaxy's matter within independent uncertainties, fit the radial and vertical pull
together, add the January 2026 Cepheid rotation curve, constrain the release length jointly, keep the original
wide-binary forecast and label revisions as amendments); complete the novelty review now; then test a frozen
prediction. This round does each of these once.

### 29.1 What is borrowed, what may be new

`NOVELTY.md` maps every central claim to the closest published work found (checked on 25 September 2026, abstracts
and summaries only): the field equation's architecture is Milgrom's QUMOND (borrowed); the heat weight has the form of
Tolman and Whittaker's active mass ρ + 3p/c² with the companion's u for c (known form; our differences are the speed
and that colliding gas does not count); collisions switching the heat off is Dicke narrowing (borrowed physics); a
streaming medium mediating attraction is the Fatio–Le Sage class with its drag and heating objections (known class);
in-phase emitters attracting is Bjerknes' (known class); gain bodies pulled toward the light they amplify is published
(Mizrahi & Fainman 2010; Gao et al. 2017); one-way coupling through a reservoir is cascaded-systems theory and
Metelmann & Clerk's reservoir engineering (known class). Not found in this search: an above-threshold, self-sustained
emitter pulled in proportion to the wave's height; lensing that follows the galaxies in colliding clusters because
hot galaxies carry the companion; a velocity-dispersion-weighted source of lensing at fixed visible mass. The document
ends with a candidate novelty statement and the searches still owed.

### 29.2 The absorbing stream as an exact wave medium

`code/full_wave_v19.py` → `run-full-wave-v19/full_wave_v19.json`; `code/self_force_sparc_v19.py` →
`run-full-wave-v19/self_force_sparc_v19.json`.

Round 18 built the stream that absorbs counter-moving waves in its ray form. Here it is solved exactly, for a uniform
stream (the geometry a receiver sees far from a source), in Fourier space: a medium adds −iΠ(q) to the wave equation,
G(q) = 1/(q² − k² − iΠ(q)), and it is passive (no arrangement of emitters can draw energy from it) exactly when
Π(q) ≥ 0 for every wavevector q. Units: wavelength 1, k = 2π.

**What absorbs.** Absorbers moving along the stream at the wave's own speed see a wave of wavevector q at the frequency
ω − u q·e = k − q_z (the calculation sets u = c = 1; see the correction below the Reading). A wave moving with the stream is at zero frequency for them, and every passive absorber is
transparent at zero frequency: the medium is one-way for propagating waves with no tuning (the Doppler effect does
it), with the attenuation (κ/2)(1 − cos θ) per unit length for a wave at angle θ to the stream (a cardioid). Two kinds:

* **Point-like absorbers** (Π = κ(k − q_z) for every q) act on the emitters' near fields too, where components move
  more slowly than the stream (q_z > k, Π < 0): the stream pumps them and drags every emitter downstream. The drag,
  in units of the emitted power over the wave speed, grows without bound with the emitter's sharpness (the integral
  cut at 4, 8, 16, 32 k): +0.06, +0.13, +0.27, +0.54 at κ = 0.5; +0.25, +0.59, +1.2, +2.5 at κ = 2.3. **Ruled out.**
* **Absorbers about a wavelength across**, which cannot respond to finer structure (a form factor e^{−(|q| − k)²/2w²},
  direction read from q̂): Π ≥ 0 everywhere, **passive at every strength, by construction.** The smallest eigenvalue
  of 24 emitters' dissipative matrix stays at the free medium's +0.016 (+0.016 to +0.040) at κ = 0.5–20.

| medium (w = k/2) | κ | backward transmission at 1, 2, 4 λ (ray form) | sideways at 2 λ (ray) | with the stream at 2 λ | near field changed at 0.15 λ | self-force / (P/c) | power fed |
|---|---|---|---|---|---|---|---|
| cardioid (the Doppler absorbers) | 0.5 | 0.619, 0.372, 0.135 (0.607, 0.368, 0.135) | 0.611 (0.607) | 0.996 | 7% | −0.015 | ×0.94 |
| | 2.3 | 0.108, 0.009, 0.000 (0.100, 0.010, 0.000) | 0.098 (0.100) | 0.973 | 23% | −0.044 | ×0.81 |
| | 5 | 0.053, 0.004, 0.000 (0.007, 0.000, 0.000) | 0.008 (0.007) | 0.921 | 37% | −0.063 | ×0.70 |
| | 20 | 0.022, 0.007, 0.001 (0, 0, 0) | 0.009 (0) | 0.620 | 63% | −0.088 | ×0.50 |
| inward only (round 18's rule) | 0.5 | 0.611, 0.367, 0.133 (0.607, 0.368, 0.135) | 0.924 (1) | 0.997 | 4% | −0.015 | ×0.97 |
| | 2.3 | 0.091, 0.033, 0.013 (0.100, 0.010, 0.000) | 0.744 (1) | 0.986 | 15% | −0.048 | ×0.90 |
| | 5 | 0.171, 0.059, 0.030 (0.007, 0, 0) | 0.618 (1) | 0.971 | 24% | −0.075 | ×0.85 |
| | 20 | 0.279, 0.224, 0.135 (0, 0, 0) | 0.523 (1) | 0.913 | 42% | −0.127 | ×0.75 |

(With w = k/4, larger absorbers, the same pattern with more leakage and about twice the self-force: backward
transmission 0.225 and 0.073 at κ = 2.3 for the inward rule; self-force −0.03 to −0.22.)

* **The ray form is right where absorption is weak, and wrong where it is strong.** At κ = 0.5 (an absorption length
  of two wavelengths) the exact medium matches the ray form to 1–5%; the cardioid matches it to 10% up to κ = 2.3.
  Beyond, a passive medium cannot be more one-way than its absorbers can resolve: the backward coupling levels off at
  a few per cent (cardioid) to 10–28% (inward rule, whose sharp edge at 90° diffracts) at 1–2 wavelengths, where the
  ray form gives 0.7% and less. Round 18's strong-absorption runs (κ = 5, 20) overstated the one-way property.
* **The near field** is kept to within 4–7% at a sixth of a wavelength at κ = 0.5, but not at strong absorption.
* **Energy and momentum.** An emitter feeds the absorbing medium 3–30% less power than empty space (its surroundings
  absorb part of its own field). What the stream absorbs, with the momentum of the absorbed waves, goes to the stream.
* **A self-pull, not a push.** A lone emitter is pulled against the stream (toward a source's centre) by 1.5% of its
  emitted power over the wave speed at κ = 0.5, 4–9% at κ = 2.3–20 (twice that with larger absorbers). Per kilogram
  this is a constant acceleration η a (a = 2ℓ/u) where the companion is released, η ≈ (0.015–0.09) × u/2c. **SPARC
  allows it:** with a refitted, η = 0.01, 0.03, 0.05, 0.1 give 15.89, 15.95, 16.03, 16.24 km/s against 15.87 (the
  held-out galaxies 18.98, 18.93, 18.92, 18.96 against 19.01), a absorbing most of it (6.3 → 5.9, 5.3, 4.8, 4.0 × 10⁻¹¹).

**Reading.** A physically consistent absorbing stream exists: Doppler absorbers carried at the wave's speed, at least
about a wavelength across; it is passive at every strength, keeps the near field and pulls emitters only slightly.
It is accurately described by round 18's ray form at weak absorption (an absorption length of a wavelength or more),
and it cannot be strongly one-way over distances shorter than a few wavelengths. The next section finds that the data
allow only weak absorption anyway.

*Correction (25 September, §30).* The calculation takes the absorbers to move at the wave's own speed (u = c). But the
energy bill of a stream that absorbs counter-moving waves (§28.4B) needs the companion's stream to move at least twice
as fast as the wave's crests (c ≤ u/2). Absorbers carried by that stream at u would see forward waves within
cos θ > c/u (within 60° of the stream at c = u/2) at negative frequency and amplify them, the pumping that rules out
point absorbers above. So the passive, one-way medium found here exists only if the absorbers are a separate component
moving at the crests' speed, not the companion's stream itself. What such absorbers are is open.

### 29.3 What the medium hears, tested on the clusters and the full suite (the hot-shell benchmark)

`code/hot_shell_v19.py` → `run-hot-shell-v19/hot_shell_v19.json`; the suite runs in `run-hot-shell-v19/suite/`
(candidates `regression/candidates/stream_k3*.json`, `stream_k10_refit`, `stream_k30_refit`, `one_way_refit`); the
hearing rules in `code/law.py` (`HOT_GEOMETRY`, `STREAM_KAPPA`, `stream_weights`).

**A geometric result first.** For a round source, a stream that absorbs every inward-travelling wave lets a receiver
at radius R hear, of a shell of radius s < R, exactly the cap cos θ ≥ s/R, and that cap carries **exactly half** of
the shell's 1/d² sum (artanh(x)/2x against artanh(x)/x, x = s/R) and **exactly half** of its net flux (1/2 against
Gauss's 1). So for round sources the absorbing stream only halves what is heard of all inner matter, cold and hot
alike, which rescales a; the one real difference from the law is that shells outside the receiver are not heard.

**How much of the law's heat term comes from outside.** In the 12 X-COP clusters, the median share of S at 0.1, 0.2,
0.3, 0.5, 0.7 and 1 R500 that comes from hot galaxies farther out than the receiver: 69, 63, 53, 46, 37, 28%.

| hearing (hot glow; cold glow for the stream) | X-COP rms, adopted constants | u refitted | X-COP rms | mean miss by radius (0.1 … 1 R500) |
|---|---|---|---|---|
| two-way (the law) | 0.222 | 169.4 | 0.222 | +0.10 +0.12 +0.08 +0.02 −0.09 −0.24 |
| one-way, inner shells only | 0.408 | 115.9 | 0.281 | +0.30 +0.23 +0.09 −0.03 −0.20 −0.40 |
| one-way net flux (Gauss) | 0.486 | 101.4 | 0.288 | +0.31 +0.24 +0.10 −0.03 −0.21 −0.42 |
| absorbing stream, 1/Mpc | 0.226 | 155.5 | 0.215 | +0.10 +0.13 +0.08 +0.02 −0.09 −0.24 |
| 3/Mpc (333 kpc) | 0.258 | 140.9 | 0.210 | +0.10 +0.13 +0.08 +0.02 −0.10 −0.25 |
| 10/Mpc (100 kpc) | 0.340 | 121.8 | 0.211 | +0.12 +0.15 +0.09 +0.02 −0.11 −0.27 |
| 30/Mpc (33 kpc) | 0.429 | 106.9 | 0.221 | +0.15 +0.17 +0.09 +0.01 −0.13 −0.30 |
| 100/Mpc (10 kpc) | 0.512 | 95.6 | 0.237 | +0.19 +0.19 +0.09 +0.00 −0.15 −0.33 |

(The stream hears the near side of outer shells within an absorption length, which a strict one-way rule does not:
that is why it keeps the clusters' shape where the strict rule doubles the radial trend.)

**Without refitting**, the stream at 3/Mpc changes only the clusters in the whole suite (35 pass, 8 close, 6 fail
against 36/7/6; the X-COP radial trend 0.239 → 0.285). **With a refitted on SPARC and u on X-COP:**

| candidate | a (10⁻¹¹), u | pass / close / fail | SPARC bulges (km/s) | X-COP rms, trend | KiDS all, red | lensing speeds of ellipticals (rms z) | SLACS light vs matter (dex) |
|---|---|---|---|---|---|---|---|
| the law | 6.30, 169.4 | 36 / 7 / 6 | 29.4 | 0.222, 0.239 | 0.065, 0.078 | 2.55 | −0.028 |
| stream 3/Mpc | 6.21, 139.7 | 35 / 8 / 6 | 30.4 | 0.210, 0.251 | 0.024, 0.029 | 5.15 | −0.060 |
| stream 10/Mpc | 6.12, 119.9 | 32 / 10 / 7 | 31.7 | 0.211, 0.271 | −0.012, −0.014 | 7.67 | −0.094 |
| stream 30/Mpc | 6.02, 104.4 | 30 / 9 / 10 | 33.4 | 0.221, 0.300 | −0.048, −0.056 | 10.2 | −0.131 |
| one-way (inner shells) | 6.08, 113.8 | 32 / 6 / 11 | 32.3 | 0.281, 0.399 | −0.026, −0.030 | 8.60 | −0.083 |

**Reading.** The data need the hot glow to reach inward. A strictly one-way medium fails (11 failures after
refitting); an absorbing stream is tolerated only with an absorption length of about 300 kpc or more (3/Mpc: the
clusters and KiDS improve, the ellipticals' lensing speeds and SLACS worsen, because the lower u the clusters then ask
for makes hot stars louder in every galaxy; the same trade as round 14's heat exponent 1.75, §24.1). So on the scale of galaxies and clusters the medium can be at most weakly
one-way: an absorption length longer than a galaxy and comparable to a cluster's radius. That is also where the exact
medium of §29.2 agrees with its ray form. It also means the strong one-way absorption the models used to keep a warm
source in step (κ R ≈ 7–15 across the source) is excluded at the scale of clusters: the warm source must keep its
beat some other way, or be smaller than the medium's wavelength (§29.4 tests sources of about a wavelength).

### 29.4 The law the models produce: mass, distance, heat and density, with their spread

`code/emergent_law_v19.py` → `run-emergent-law-v19/emergent_law_v19.json` (84 runs).

Round 16's one kind of matter (structure 'single'), sources of 24, 48 and 96 pieces at round 16's density (radius
2.38, 3, 3.78), and 48 pieces squeezed or spread (radius 1, 2 and 4.5; radius 1 is about a wavelength), cold, warm
(k = 2, 8) and colliding (k = 8, ν = 50); 12 receivers of the same matter at each of r = 9 and 18; three
arrangements each. Two media: the two-way wave, and the stream that absorbs inward waves at κ = 1 per wavelength in its
ray form, the strongest the clusters allow (§29.3: a source three wavelengths across is then three absorption lengths
across) and where the ray form is exact to 5% (§29.2). Mean net pull on the receivers (10⁻⁵; the spread is the
standard error over the three arrangements):

| source | medium | r = 9 | r = 18 |
|---|---|---|---|
| cold, 24 / 48 / 96 pieces | two-way | +3.6 ± 1.3 / +7.1 ± 1.2 / +13.5 ± 2.3 | +1.8 ± 1.5 / +3.5 ± 1.1 / +5.3 ± 0.8 |
| | absorbing, κ = 1 | +3.2 ± 0.6 / +7.5 ± 0.4 / +11.2 ± 1.9 | +1.5 ± 1.1 / +3.6 ± 1.1 / +3.4 ± 0.6 |
| warm k = 2, 48 pieces | two-way / κ = 1 | +6.5 ± 2.7 / +7.1 ± 1.0 | +2.5 ± 2.3 / +3.3 ± 1.1 |
| warm k = 8, 24 / 48 / 96 pieces | two-way | −0.6 / −0.3 / −9.8 ± 2.4 | +0.6 / +1.0 / −3.8 ± 1.7 |
| | κ = 1 | −0.6 / −4.2 ± 1.7 / −3.6 ± 4.2 | +1.1 / −1.5 ± 0.8 / +2.5 ± 2.2 |
| colliding k = 8, 48 pieces | two-way / κ = 1 | +4.8 ± 1.7 / +6.6 ± 0.3 | +3.4 ± 1.2 / +3.5 ± 1.3 |
| cold, 48 pieces, radius 1 / 2 / 3 / 4.5 | two-way | +0.5 / −1.6 / +7.1 / +8.9 | +1.4 / +0.3 / +3.5 / +2.5 |
| | κ = 1 | −0.6 / +3.9 / +7.5 / +6.3 | +1.5 / +1.4 / +3.6 / +2.8 |

* **Cold matter gives about the law's √M/r.** Fitting pull = A M^p r^(−q) over the cold sources, with the arrangements
  resampled: p = 0.58 ± 0.13 and q = 1.11 ± 0.22 in the absorbing stream, 0.80 ± 0.29 and 0.93 ± 0.20 two-way (the
  law: 0.5 and 1). The pull halves from r = 9 to 18 (×2.1 for 48 pieces), and four times the mass gives ×2.3 at
  r = 18 but ×3.5 at r = 9, where the largest source is only 2.4 of its radii away: far receivers see the square root,
  near ones closer to a straight proportion.
* **The heat gain is not there.** Against cold sources of the same size and arrangement, warm sources at k = 2 pull
  0.90 ± 0.37 (two-way) and 0.97 ± 0.16 (κ = 1) as hard, where the law needs √3 = 1.73; at k = 8 they push nearer
  receivers (gain 0.47 ± 0.74 and −0.52 ± 0.15; the law 3). Colliding sources pull 0.78 ± 0.11 and 0.90 ± 0.05 as hard
  as cold ones (the law: 1). So in the medium the data allow, the models reproduce the cold law and the collision rule
  but not the law's heat term; the absorption that made warm sources pull harder in rounds 17–18 (κR ≈ 7–15 across the
  source, or a rule imposed by hand) is what the clusters exclude. *(Correction, 25 September, §30: the quoted errors
  treat two distances × three arrangements as six independent values, which understates them. The arrangements
  scatter widely: with the two-way wave at k = 2 one arrangement gives 2.21 and 1.79 at r = 9 and 18, the others
  0.85/0.29 and 0.26/0.01; at k = 8 one gives 1.70 and 3.47 while the other two push. The fair reading is "not
  reproduced on average, with very large scatter"; three arrangements cannot settle the heat gain in the two-way
  medium. At κ = 1 the k = 8 sources push in every arrangement.)*
* **Size matters in the models, and should not.** At fixed mass, sources of radius 1–2 wavelengths pull 2–15 times less
  than sources of radius 3–4.5 (and sometimes push): their pieces, within a wavelength of each other, settle into
  collective states that radiate weakly. Real galaxies of one mass rotate alike whatever their size (the baryonic
  Tully–Fisher relation has no size term), so either real sources are always many companion wavelengths across (a
  wavelength well below a kiloparsec, which with the clusters' absorption length makes the absorption per wavelength
  tiny) or the models' collective states are wrong for real matter.
* **Energy** is booked to 10⁻¹² (10⁻⁷ with collisions) in every run; the warm sources' rhythm spread is 3–7 × 10⁻³
  against 3–6 × 10⁻⁴ cold, in both media.

**Reading.** Review item 2 asked for the law the medium produces. The cold part comes out close to the law's
square-root-of-mass, inverse-distance form (within its uncertainties at the far receivers), and colliding matter pulls
like cold matter, as the law says; the heat term does not, in the only media the data allow. The data's case for the
heat term (the clusters' masses, the ellipticals' lensing, the collisions) stands on its own; its mechanism is open
again, and it now has two hard constraints: the hot glow must reach inward (§29.3), and whatever keeps a warm source in
step must work without strong one-way absorption.

### 29.5 The Milky Way fitted to independent constraints; Cassini and wide binaries from it

`code/mw_joint_v19.py` → `run-mw-joint-v19/mw_joint_v19.json` (the bulge's heat by the SPARC rule, as in every galaxy)
and `mw_joint_v19_bulge120.json` (the bulge's measured dispersion, 120 km/s); the literature values (checked against
the papers) are listed in the script's header.

**The data and the priors.** The rotation curve of Feng, Huang, Zhang & Liu 2026 (MNRAS 546, stag011; arXiv:2512.21780):
903 classical Cepheids from Gaia DR3, 6.6–17.6 kpc in 12 bins, with R0 = 8.275 kpc and the Sun's azimuthal speed
250.2 km/s (the Sgr A* proper motion); statistical errors 0.9–2.2 km/s, to which a 3 km/s floor is added for the dip
and bump no axisymmetric model has, and a common scale (1 ± 0.02, the Sun's own speed) marginalised in closed form.
The curves already in the suite (Eilers 2019, Zhou 2023, Ou 2024) are fitted separately, each in its own convention
(R0 = 8.122–8.178 kpc; all four tie the Sun's speed to R0 through the same Sgr A* motion). The vertical pull:
K_z(1.1 kpc)/2πG = 70 ± 5 M☉/pc² (Bland-Hawthorn & Gerhard 2016's consensus of 67–74). The matter, varied within
independent measurements: stars at the Sun 33.4 ± 3 M☉/pc² (McKee, Parravano & Hollenbach 2015), the thick disk's
share 0.12 ± 0.04 and the scale lengths 2.6 ± 0.5 (thin) and 2.0 ± 0.6 kpc (thick; Bland-Hawthorn & Gerhard 2016),
gas 13.7 ± 1.6 M☉/pc² (McKee et al.), the bulge (1.2 ± 0.35) × 10¹⁰ M☉ (Licquia & Newman 2015's 0.91 to
Bland-Hawthorn & Gerhard's 1.4–1.7). Scale lengths on a grid (thin 2.1, 2.6, 3.1; thick 2.0, 3.0 kpc), heights fixed
at 0.3 and 0.9 kpc, normalisations fitted continuously. The law's constants are held (SPARC and X-COP).

| curve (points) | law | best χ² (curve part) | common scale | stars Σ* (33.4 ± 3) | bulge (10¹⁰) | K_z(1.1) |
|---|---|---|---|---|---|---|
| Feng 2026 (12) | ours | 35.0 (23.2) | 0.949 | 34.6 | 2.28 | 77.1 |
| | Newton, same matter | 393 (335) | 0.730 | 49.1 | 1.64 | 73.6 |
| | ours, bulge at its measured 120 km/s | 38.1 (27.6) | 0.935 | 35.4 | 2.14 | 78.1 |
| Eilers 2019 (38) | ours | 33.9 (31.5) | 0.966 | 33.8 | 1.30 | 75.3 |
| Zhou 2023 (34) | ours | 35.4 (32.6) | 0.947 | 34.6 | 1.19 | 76.4 |
| Ou 2024 (35) | ours | 61.7 (57.4) | 0.956 | 35.2 | 1.15 | 76.5 |

(The curve part includes the scale's own prior term; with Feng's curve the best thin disk is 2.6 kpc, with the others
2.1 kpc, close to Bovy & Rix's mass-weighted 2.15.)

* **The shape fits; the level is 3–6% low.** With the matter at its measured values the law follows every curve's
  shape (after the common scale, the Cepheid residuals are within ±6 km/s, χ² ≈ 17 for 12 points with the 3 km/s
  floor), but each measured curve must be scaled down by 3.4–6.5% to meet it: the law's speed is 3–6% below the
  Galaxy's. Newton with the same matter fails by far (27%).
* **This is the offset the law has in the galaxies too.** In SPARC, at the Milky Way's acceleration (log g_N = −10 to
  −9.75), the median log(g_obs/g_law) is +0.026 ± 0.015 dex, and +0.022 to +0.038 dex in every bin from −11.25 to
  −9.25: the law runs about 6% low in acceleration, 3% in speed, for the typical galaxy (the fit minimises a mean
  square, not the median). The Milky Way is 1–3% lower still in speed: not an outlier. The round-7 shortfall ("209
  against 229–234 km/s", 20–25% in acceleration) was mostly McMillan's matter model, itself fitted together with a
  dark halo, whose stars at the Sun (45.8 M☉/pc²) exceed the local census by 4σ and whose disk is long.
* **What the Cepheid curve still asks.** Its high inner points (243 km/s at 6.6 kpc) want a heavier bulge than the
  prior (2.1–2.3 × 10¹⁰, 3σ; with the bulge's measured dispersion as heavy), and the vertical pull comes out 1.4σ high
  (77–78 against 70 ± 5). The other three curves are met with the bulge at its prior.

**Cassini and wide binaries from the fitted Galaxy** (`regression/t_precision.sun_in_galaxy` with the fitted pull and
heat at the Sun: g_N = 1.37 × 10⁻¹⁰, S = 9.7 × 10⁻¹¹ m/s² for the best fit):

| Galactic field at the Sun | Q₂ at L = 0.15 pc (10⁻²⁷ s⁻²) | from (1.6 ± 1.8) | L for 1σ, 2σ (pc) | wide binaries at 7,000 / 20,000 AU |
|---|---|---|---|---|
| the suite's (230 km/s inverted, no heat) | 4.38 | +1.54σ | 0.19, 0.12 | +3.2%, +7.6% |
| fitted Galaxy, bulge by the SPARC rule | 1.96 (1.47–2.30 over the six disk shapes) | +0.20σ (−0.07 to +0.39) | 0.087, 0.049 | +2.4%, +5.7% |
| fitted Galaxy, bulge at 120 km/s | 2.37 | +0.43σ | 0.10, 0.06 | +2.8%, +6.8% |

With the Galaxy fitted to its own measurements, Cassini's 2026 value is met at the adopted release length (0.2–0.4σ;
it was 1.5σ with the suite's rough field). **Wide binaries: the forecast locked in round 10 is unchanged**
(`forecasts/wide_binaries_gaia_dr4_v10.json`: 1.039 at 7,000 AU, 1.093 at 20,000 AU, SHA-256 94a7a470…); the values
above are **amendments** (the round-12 constants and the fitted Galactic field), labelled as such: 1.024–1.028 and
1.057–1.068.

### 29.6 One process for screening and release (a supplied derivation), carried to the Solar System

`screening-blockers-v19/` (supplied by the user on 25 September 2026; `PROVENANCE.md`); `code/blocker_release_v19.py` →
`run-blocker-release-v19/blocker_release_v19.json`.

**The proposal.** A coupling element of the companion holds n temporary "blockers", created at rate η g and each
removed at rate γ, and couples only when it holds none. The stationary population is Poisson, so the open fraction is
exp(−ηg/γ) = exp(−g/g_d) with g_d = γ/η: the law's screening, emerging from the chance of holding no blocker. If a newly
emitted element inherits the local equilibrium population plus one blocker made by the emission, its open fraction in
a constant field is exp(−g/g_d)(1 − e^{−r/L}) with L = u/γ: the round-9 release, from the same lifetime (868 years for
L = 0.15 pc; g_d L = u/η). In a changing field the load m obeys dm/dt = ηg − γm and the open fraction is (1 − s)e^{−m}:
the logarithm of the open fraction relaxes, not the fraction itself, which a strong-to-weak field step distinguishes
from a plain relaxation (half release after 2.7 lifetimes against 0.7 for g/g_d = 10 → 0), while small disturbances
cannot. **Checked here:** the package's master equation and event simulation reproduce (`PROVENANCE.md`); in a
constant field our implementation gives back 1 − e^{−r/L} to 10⁻⁹.

**What it does to the Sun's companion.** The companion is launched in the Sun's own field, g/g_d ≈ 1.4 × 10¹² at the
surface, and moves outward at u while the field falls much faster than blockers clear (for r ≪ L). It carries a load
that clears only after about ln(load) lifetimes, so the release becomes a sharp switch far beyond L:

| launch preparation (Galactic field from §29.5) | half release at L = 0.15 pc | Q₂ (10⁻²⁷) | binaries 7,000 / 20,000 AU | lifetime for the adopted law's Q₂ → half release, binaries |
|---|---|---|---|---|
| adopted scalar release (1 − e^{−r/L}) | 20,800 AU | 2.02 (+0.23σ) | +2.5%, +5.8% | – |
| equilibrium at the surface + one blocker | 850,000 AU (4.1 pc) | ≈ 0 | 0, 0 | 6 yr → 5,800 AU, +11.7%, +12.0% |
| empty at the surface + one blocker | 380,000 AU | ≈ 0 | 0, 0 | 9 yr → 5,900 AU, +11.2%, +12.0% |
| K = 10 sites, full at launch | 79,000 AU | ≈ 0 | 0, 0 | 58 yr → 7,400 AU, +5.3%, +11.9% |
| K = 100 sites, full at launch | 149,000 AU | ≈ 0 | 0, 0 | 29 yr → 6,400 AU, +7.7%, +12.0% |

**Reading.** The derivation is a real step: one lifetime gives both of the law's factors in constant fields, and a
changing field distinguishes it from a plain relaxation. For the Solar System it changes the logic: with a load made
at launch, the release is a switch at L·ln(load), so Cassini sets a *lower* bound on the lifetime (about 6–60 years
for the adopted Q₂, depending on the launch), not a value to be tuned, and at the derivation's own calibration
(868 years) the Sun's companion stays closed out to 0.4–4 pc: no Q₂ at all (the measurement is 0.9σ from zero) and
wide binaries exactly Newtonian. Wide binaries then decide: a switch shows as a **step** in the extra pull with
separation (nothing inside the switch radius, the full effect beyond), unlike the adopted law's gradual rise; with the
same Cassini standing the step gives +12% at 20,000 AU, twice the adopted law. Galaxies are unaffected (every star's
launch load clears within a few parsecs). It adds no constant (η, γ replace g_d, L), the blockers are not yet found in
the matter model, and whether the force follows the open fraction or its square root is open; it is recorded as a
separately named candidate, not adopted.

### 29.7 A frozen prediction: lensing against the stars' own speeds

`code/frozen_prediction_v19.py` → `run-frozen-prediction-v19/frozen_prediction_v19.json` and `SHA256`
(01219170346b60123a0dc5bf979cf847976ad9e7cbd178efca52c10a8b0fd093), written before this project looked at any lensing
split by velocity dispersion. *(Note, 25 September: the fingerprint is the SHA-256 of the prediction text as the
script wrote it, without the file's final newline; `sha256sum` of the file itself gives
c0b1056771313364fae91530f51adf921236d64d9738e523fad63cd8b201c427. The content is unchanged.)*

At fixed stellar mass and isolation, the law's heat term makes a galaxy whose stars move faster lens more: far out,
√(a G M (1 + k)), k = 3σ²/u². For isolated bulge-dominated lenses (the whole stellar mass at the central dispersion
σ_e) of log M* = 10.6 at z = 0.25 (the KiDS-1000 stacks' typical lens), computed with the adopted law and the suite's
own forward model at the survey's g_bar bins, the lensing acceleration relative to σ_e = 200 km/s, median over the
reliable bins below 10⁻¹² m/s²:

| σ_e (km/s) | 100 | 150 | 200 | 250 | 300 |
|---|---|---|---|---|---|
| k = 3σ²/u² | 1.05 | 2.35 | 4.18 | 6.53 | 9.41 |
| log g_obs − log g_obs(200) | −0.183 | −0.087 | 0 | +0.077 | +0.144 |

So 250 against 150 km/s: +0.164 dex (×1.46; the deep-regime √((1 + k₁)/(1 + k₂)) gives 0.176). MOND predicts no
dependence at fixed visible mass; in dark-matter models any dependence runs through the halo mass that goes with σ at
fixed stellar mass, a separate, measurable relation. The test: split isolated lenses of one stellar-mass bin by their
measured σ_e (SDSS/GAMA spectra), stack each group's excess surface density, and compare the ratios with this table.

### 29.8 Where round 19 leaves the checklist

| review item | status after round 19 |
|---|---|
| 1. the one-way medium in full wave form | **done for a uniform stream** (§29.2): what absorbs (Doppler absorbers at the wave's speed, at least a wavelength across; point absorbers drag and are ruled out), passive at every strength, the near field kept at weak absorption, the self-force (−0.015 to −0.13 P/c, allowed by SPARC up to ≈ 0.05 a), the ray form's range of validity. Still to do: the radial geometry exactly, and what sets κ |
| 2. the law the medium produces | **the hot-shell benchmark run on data** (§29.3): the hot glow must reach inward; a stream is allowed only with absorption lengths ≳ 300 kpc. **Exponents** (§29.4): cold matter gives about the law's √M/r (p = 0.58 ± 0.13, q = 1.11 ± 0.22 in the absorbing stream, 0.80 ± 0.29 and 0.93 ± 0.20 two-way; the law 0.5 and 1) and colliding matter pulls like cold matter (0.78–0.90); the heat gain is absent in both media (0.90–0.97 at k = 2 against the law's 1.73; pushes at k = 8), and compact sources pull 2–15 times less at the same mass, where the law has no size term |
| 3. inertia and universal free fall | open |
| 4. screening, light, Cassini | Cassini with the fitted Galaxy: 0.2–0.4σ at the adopted L (§29.5); a candidate process for screening and release from one lifetime, with a distinct Solar-System signature (§29.6); light's response still assumed |
| 5. collisions as a calculation | open |
| 6. statistical inference | a frozen prediction written (§29.7); the likelihood programme open |
| 7. remaining discrepancies | the Milky Way's shortfall reduced to the law's own 3–6% level with independent matter (§29.5) |
| 8. novelty and scope | `NOVELTY.md` (§29.1) |
| 9. reproducible release | open (the suite's candidate runs and every new script's outputs are in the repository) |
| 10. the manuscript | the write-up updated to the evidence (blog, page) |

**Next, after round 19:**
1. **Keep a warm source in step without strong one-way absorption.** The data allow at most weak absorption (§29.3), and
   there the models give no heat gain (§29.4). Candidates: the reservoir-engineered (Metelmann–Clerk) one-way coupling
   applied inside a source only, with the far field two-way; receivers whose locking survives a detuned source; and
   the models' collective states in sources about a wavelength across, which also make size matter where it should
   not.
2. **The radial medium exactly** (partial waves), with κ set by the stream's density, and the self-force's scaling.
3. **The blocker candidate in the matter model:** look for a metastable blocking excitation, measure its formation
   against the field and its lifetime separately, and derive whether the force follows the open fraction or its square
   root; then wide binaries as the test (a step against a gradual rise).
4. **The law's 3% at the Milky Way's pull:** the SPARC median offset (+0.026 dex at g_N ≈ 10⁻¹⁰) with a robust fit
   statistic; then the Milky Way again.
5. **The frozen prediction:** σ-split lensing from KiDS × GAMA or SDSS × HSC.
6. **Inertia;** the likelihood programme with the frozen law; a reproducible release.

## 30. A step back, 25 September 2026: what did not go our way, and which doors are still open

The owner asked for a list of the results we would prefer were different, what was tried and why it was not
enough, and whether the derivation work closed doors or boxed itself in. The full answer is
[STEP-BACK-AUDIT.md](../../../STEP-BACK-AUDIT.md), from eight read-throughs of the whole record (§1–§29, THEORY, the
roadmap, the novelty audit, the rules and the earlier photon–companion era), with the key numbers re-checked against
the result files. It changes nothing in the law, its constants or the locked forecasts. In short:

* **The law is not cornered.** Its standing misses (five faint dwarfs, the Bullet's smaller half, KiDS's red lenses;
  heavy stars in strong lenses; the Milky Way's 3–6%; the clusters' radial trend) each have untried remedies, and
  several trace to our own conventions: the static distance law's exact form (in the survey's own distances KiDS's
  level is +0.013, red lenses +0.029), star-mass conventions, and cluster stars' speeds taken from the X-ray-measured
  gravity they are used to predict (with the law's own gravity the X-COP miss rises from 0.227 to 0.329, §10.2).
* **The derivation work is in a corner of its own making, with a marked way out.** Since round 10 every matter model
  makes the pull by rhythm-locking (§20.1). Every setback with warm matter since then follows from needing a shared
  rhythm: dilution by independent sources (§20.2), scrambled waves (§23.5), warm sources out of tune (§26.4, §27),
  the one-way wave (§27.3–§28.5), no heat gain in the medium the data allow (§29.4). §23.5 and §23.10 named the way
  out, "a receiver that responds to the wave's energy rather than its phase" and "a pull drawn from the companion's
  energy flow", and it was never tried. The data-selected rule (strength from all the companion present, direction
  from its net flow, warm contributions adding without cancelling) is what such a pull would give.
* **Most model exclusions are narrow.** The force experiments never had moving pieces (heat was a frozen random
  mixing, no Doppler shifts; round 15's glow-only `wave_dark_v15.py` did move its pieces, §31.1), never went beyond k = 16 (cluster galaxies have about 25–150), never had more than about 100
  pieces or sources more than 4.5 wavelengths across, and always used one scalar channel with senders and receivers
  identical. The strength-dependent rhythm rule of §23.7, which gave the heat term's pattern, was dropped with the
  move to the store model, not refuted. The reciprocity argument of §27.2 covers linear, passive, time-invariant,
  single-channel structures only; the two best structures were never run in the full model.
* **Owner's decisions:** the QUMOND-form field equation against RULES.md §1; the dwarf fix `no_hold` shelved for
  lack of a derivation while L was adopted without one; the distance law's form, and writing the no-expanding-universe
  rule into RULES.md.

**Corrections made with the audit** (none changes the law, a constant, a forecast or a figure):
1. §29.2: the passive one-way medium needs absorbers moving at the wave's speed; the companion's stream must be at
   least twice as fast (§28.4B), so the absorbers cannot be the stream itself (note after §29.2's Reading).
2. §28.5: "73–84% of the square root the law needs" is relative to the model's own glow; against the law's heat gain
   the raw ratio at k = 8 is about 1.6 against 3 (note in place; blog and page corrected).
3. §29.7: the fingerprint is over the text without the file's final newline (note in place).
4. §29.4: the heat-gain errors are understated and the arrangements scatter widely (note in place).
5. §27.2 and §28.3, §28.5: unit labels (wavelength 1, not 2π; the dense and dilute balls' spacings are about 1.3 and
   13 wavelengths).
6. §21.5: Mistele et al.'s four ratios are stellar-mass bins, not radii.
7. §27.4: the one-way full-model runs (`run-one-matter-v17/oneway_*.json`) booked energy but ran with momentum
   bookkeeping off, so momentum balance was not checked there.

## 31. Round 20, 25 September 2026: a second mechanism track, and the cluster heat without X-ray input

The owner supplied a review of the step-back audit (§30) with a new calculation ([energy-shift-v20/](energy-shift-v20/),
reproduced here to 1 part in 10⁹; its provenance in `energy-shift-v20/PROVENANCE.md`). Its recommendations:
* Make the energy-shift approach the next independent mechanism track: "gravity comes from the energy and stress of a
  coupled matter–medium state", not "gravity requires synchronized emitters".
* Keep the adopted law and the phase-locked models as comparison cases.
* Make the immediate goal a nonlinear, finite-population response.
* In parallel, run the strongest data-side test: the cluster heat input without the X-ray gravity.

This round takes the first steps on each. Nothing in the adopted law, its constants or the locked forecasts changes.

### 31.1 The audit narrowed

Four statements in `STEP-BACK-AUDIT.md` were broader than the evidence, and are now corrected there (its §8), in §30
and in the blog and page:
1. **"Pulled only if it feeds the wave"** holds for the travelling-wave recoil mechanism of §20.1, not for every force.
   Conservative forces between bodies that share a field (Casimir) and optical pulling by redirected momentum are
   different.
2. **The 85-million-year drain (§23.9)** applies to continuous-emission mechanisms paid for by motion, not to every
   static attraction.
3. **"No moving matter".** Round 15's `wave_dark_v15.py` moves its pieces (free, colliding, rotating, boosted) and
   recomputes their couplings, but measures only the released glow. The missing experiment is hot motion with the
   forces, the wave and the recoil evolved together.
4. **"No energy-sensitive receiver".** Round 15's `receivers_v15.py` has passive and below-threshold receivers whose
   forces follow the intensity (§25.3). What is missing is a nonlinear, many-direction response.

### 31.2 The supplied check: attraction through correlations, at any partial filling

**Setup.** Eight identical two-level pieces share one medium, exchanging excitations through
J(r) = −C e^{−κr}/r (C = 0.04/4π, κ = 1.354; the kernel comes from the owner's earlier structured-medium derivation,
which was not supplied here). Every excitation number is diagonalised exactly.

| excited pieces | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| rms internal force ÷ one excitation | 0 | 1 | 1.72 | 2.15 | 2.29 | 2.15 | 1.72 | 1 | 0 |

**Results:**
* In every partially filled sector, all 28 pair contributions attract.
* Every piece's own mean oscillation is exactly zero: the force comes from pair correlations, not from a rhythm on
  any piece.
* The force equals the energy's derivative to 6 × 10⁻¹³.
* With every piece excited, the exchange force vanishes: more stored energy is not automatically more attraction.

**Two algebraic limits.**
* A receiver whose energy depends only on the local intensity I = KM/r² feels F = 2I f′(I)/r, so f = I^q gives
  M^q/r^(2q+1): never √M/r.
* A fixed-weight positive mixture of Yukawa forces never falls more slowly than 1/r².

### 31.3 A stiffening medium, energized: only the range changes

`code/anharmonic_medium_v20.py` → `run-anharmonic-medium-v20/` (`mc.json`, `hot.json`, `critical.json`).

**The test.** This is the review's item 3. The gapped, positive-energy medium on a lattice gets a local quartic term,
E = Σ[(m₀²/2)X² + (β/4)X⁴] + ½Σ(X_n − X_m)² with m₀² = 0.1, so the bare range is 3.2 spacings. It is energized to a
temperature T; the classical statistics are exact (Metropolis, 32³, 40,000 sweeps).

**How the results are read.** The static response χ(r) = ⟨X₀X_r⟩/T is the exchange kernel: two sources coupled to
the medium feel a force proportional to χ′(r).

| β | T | ⟨X²⟩ (MC; self-consistent estimate) | range (MC; estimate; bare 3.18) | force exponent at r = 3, 6 |
|---|---|---|---|---|
| 0 (control) | 1 | 0.2267; 0.2267 | 3.18; 3.18 | 2.76, 3.34 |
| 0.3 | 0.03 / 0.3 / 3 | 0.0068 / 0.0659 / 0.568 | 3.13 / 2.64 / 1.35 | 2.79, 3.58 / 2.93, 3.84 / 3.74, 5.78 |
| 1 | 0.03 / 0.3 / 1 / 3 | 0.0067 / 0.0625 / 0.187 / 0.473 | 2.76 / 1.91 / 1.24 / 0.93 | 2.83, 3.56 / 3.23, 4.49 / 3.81, 5.51 / 4.64, 6.21 |

* **Uniform energizing only shortens the range**, as the self-consistent harmonic estimate m_eff² = m₀² + 3β⟨X²⟩
  predicts (ranges agree to within 7%). The force falls more steeply as the excitation grows.
* **Energized regions** do the same, in the local form of the same estimate (96³ static Green's function). Tried: a
  hot ball of radius 6, heating falling as 1/r, and heating falling as 1/r².
  * The smallest force exponent anywhere is 2.9–3.2, against 2.2 in the cold medium.
  * The response does depend on where the medium is energized, but always toward a shorter reach.
* **With the gap closed** (m₀ = 0; the exact radial solution of −∇²φ + βφ³ = source):
  * weak sources give Newton's 1/r²;
  * strong sources saturate. The source-strength exponent falls from 1 to 0, passing 1/2 only in a narrow crossover,
    and the distance exponent is 2.1–2.8, tending to 2 with logarithmic corrections.
  * The law's √M/r never appears.

**Reading.** This is the review's stop rule: the nonlinearity only changed a Yukawa length. The general reason: in
equilibrium, a medium with a gap has static correlations that die off exponentially, and a positive quartic term only
widens the gap. Closing the gap gives at best Newton's 1/r², with the source's strength saturating. So the long reach
and the square root cannot come from a settled medium. They must come from a medium held out of equilibrium, and
energy flowing continuously outward from matter is the companion's defining property. That is the next construction.

### 31.4 A probe outside a source: the pull's mass and distance laws in the linear medium

`code/finite_population_probe_v20.py` → `run-finite-population-v20/finite_population_probe_v20.json`.

**Setup.** The supplied Hamiltonian and constants, for a source cluster of Ns = 3–13 pieces at unit density (its
"mass") and one probe piece at distance D = 2.5–6 from the cluster's centre:
* six random clusters each;
* the lowest state of each excitation sector, as in the supplied check (how it forms is not shown);
* the pull on the probe is −⟨∂H/∂R⟩ (checked against the energy's finite difference to 5 × 10⁻¹³).

| filling | distance exponent, D = 2.5 → 6 (Ns = 9 and 13) | mass exponent (median pull, D = 3 and 6) |
|---|---|---|
| exactly half the pieces excited | 4.9, 5.6, 6.2, 7.2, 8.5: exactly the kernel's own Yukawa slope, 2 + (κD)²/(1 + κD) | 1.13, 1.11 (Ns 3–13, odd) |
| one excitation | 9.7, 11.0, 12.4, 14.3, 17.0 (twice the kernel's decay: the probe is correlated only at second order) | 1.0 for Ns 3–9, then levels off (0.8 overall) |
| two excitations | the same as one | 0.34, 0.39 (Ns 4–13) |

**Results:**
* All 1,260 cases attract.
* The cluster's arrangement moves the pull by up to 4.7× at one excitation and 1.6× at half filling.

**Reading.** In the linear medium the pull is at best Newton-like in mass, never falls more slowly than the medium's
own Yukawa kernel, and its mass law depends on the filling. That confirms the review: the next model must change the
medium's collective spatial response itself, not its gap or its occupation alone.

### 31.5 What the cold law asks of any static, conservative medium (for comparison, as the review asks)

**Read the companion as a static field Φ:**
* its energy density |∇Φ|²/8πG is the energy density of the outflow;
* that outflow is fed at ℓ per kilogram and carried outward at u along ∇Φ.

Conservation of that energy flux, ∇·(u |∇Φ| ∇Φ/8πG) = ℓρ, is ∇·(|∇Φ| ∇Φ) = 4πG a ρ with a = 2ℓ/u. That is exactly the
deep-regime field equation of Bekenstein & Milgrom's AQUAL (1984), with their a₀ replaced by a; for a point mass it
gives |∇Φ| = √(G M a)/r. (The law's guided stream, §22.4, gives the algebraic QUMOND form instead.)

**What follows:**
* Any conservative "energy and stress" construction that reproduces the cold law has one of the MOND family's field
  equations at long range, as the galaxy data force.
* What can be new lies elsewhere: the microscopic origin, the heat term with its collision rule, and the release near
  stars.
* This bears on the owner's field-equation decision (audit §4.3). Deriving our own field equation from the companion's
  energy balance leads to the AQUAL form, which is also named in RULES.md §1.

### 31.6 The cluster heat without the X-ray input

`code/xcop_selfconsistent_v20.py` → `run-xcop-selfconsistent-v20/xcop_selfconsistent_v20.json`.

**Setup.** The suite takes the stars' random speeds from the Jeans equation in the measured (hydrostatic) gravity, the
quantity being predicted (§10.2; audit §1.6). Here they come from the law's own gravity, iterated to a fixed point
with round 3's solver, on today's inputs: the round-12 law, static distances, deprojected stars.

| stars' speeds from | u (km/s) | typical miss (rms ln) | mean miss at 0.1 … 1 R500 |
|---|---|---|---|
| the X-ray gravity (the suite) | 169.4 | 25% (0.222) | +0.10 +0.12 +0.08 +0.02 −0.09 −0.24 |
| the law's own gravity | 169.4 | 40% (0.335) | +0.04 +0.03 −0.04 −0.11 −0.24 −0.40 |
| the law's own gravity, u refitted | 186.2 | 36% (0.306) | +0.18 +0.17 +0.09 +0.01 −0.12 −0.28 |

**The comparison, at 0.1 … 1 R500 (medians):**
* The law's own gravity over the measured one: 0.83, 0.83, 0.94, 1.07, 1.26, 1.51.
* The stars' speeds, law's own over X-ray: 0.94, 1.00, 1.05, 1.15, 1.22, 1.25.

**Reading.** The cluster fit leans on the X-ray input. With the law's own gravity, the outskirts are over-predicted by
up to 1.5×. That is well beyond the non-thermal pressure X-COP's own analysis finds, about 6% at R500 and 10% at R200
(Eckert et al. 2019; that estimate uses the standard cosmology's baryon fraction). Three possible causes:
* the heat term is too strong in cluster outskirts;
* the outer galaxies' orbits are radially biased: the law's weight is 3 − 2β for an anisotropic spread (§24.1), and
  β ≈ 0.3–0.5 in cluster outskirts would lower it by 20–33%;
* the outer galaxies are still falling in.

The decisive data are independent velocity-dispersion profiles of the cluster galaxies and weak-lensing masses.

### 31.7 The no-hold rule on today's law (a registered comparison, not adopted)

`regression/candidates/no_hold_r12.json` → `run-no-hold-v20/` (the suite report).

**Setup.** Round 9's rule (the Galaxy's pull does not hold back a separate system's companion; release over 1 pc) on
the round-12 law. It is held to the same standard as the release length, since neither is derived.

**Results:**
* **Tally:** 61 pass / 12 close / 4 fail, against 59 / 11 / 7.
* **Improved:** Carina and Antlia 2 pass; Sextans and Crater II are close; the ten dwarfs' χ² falls from 137.7 to 60.4.
* **Unchanged:** Draco and Ursa Minor stay at 4.1–4.2 km/s against 9.1 and 9.5.
* **Changed predictions:** Cassini's Q₂ goes to zero, and wide binaries to +18% at 20,000 AU (against 7.6%).

It stays a registered comparison until the rule has a physical basis; the blocker picture (§29.6) is one candidate.

### 31.8 Where round 20 leaves the tracks

* **The rhythm track** (§20.1–§29): kept as the active-recoil control.
* **The energy-shift track:**
  * attraction through correlations is robust (§31.2, §31.4);
  * equilibrium media, stiffening or at their gap's closing, cannot give the law's reach or its square root (§31.3);
  * the linear medium's pull is at best Newton-like in mass with a Yukawa range (§31.4).
  
  **Next:** a medium carrying the companion's steady outward energy flux, a non-equilibrium steady state. Derive its
  energy density, energy flux, stress and correlations around a source, and the energy shift and force on a probe
  (§31.5 gives the macroscopic target), with every watt booked. Then add real motion with its recoil.
* **Data:**
  * the cluster heat needs independent galaxy dispersions and weak-lensing masses (§31.6);
  * the frozen lensing test (§29.7) stands;
  * `no_hold_r12` and alternative static distance laws are registered comparisons.


## 32. Round 21, 25 September 2026: three proposals checked (a faster glow from collisions, a distance law with its own time stretch, a relativistic action)

The owner supplied feedback with three proposals:
1. **A separate spreading speed for the heat made in collisions.** Keep u = 169.4 km/s for the heat weight and for the
   cold companion's memory, but let the heat a collision generates spread at its own speed v_h. The feedback reports
   that v_h ≈ 600–625 km/s puts both halves of the Bullet Cluster inside their lensing masses (smaller half
   2.50 × 10¹⁴, main 3.44 × 10¹⁴ suns inside 250 kpc).
2. **A distance law read as a path correction.** Round 12's bounded beam-area term (flux × 1/(1 + η b), b = z/(1 + z),
   §22.3) becomes a correction to the distance itself, D* = (ln(1 + z)/α) √(1 + η b), with D_A = D*/(1 + z) and
   D_L = (1 + z) D*, which keeps distance duality. b = z/(1 + z) solves db/dy = 1 − b with y = ln(1 + z), and η = 1/2
   may be allowed. Redshift and the supernovae's time stretch would come together from a temporal-dilation operation,
   a_out(ω) = √s a_in(sω), s = 1 + z.
3. **A relativistic action.** The companion's occupation I sets the size of the extra pull and its flow direction s^μ
   the direction, through a constraint on a response field χ; matter and light see one metric with Φ_N + χ in both of
   its parts (no slip).

Each was checked against the data and the law's own bookkeeping. Nothing in the adopted law, its constants or the
locked forecasts changes.

### 32.1 The Bullet with a faster glow: reproduced, and what it costs in energy

`code/hot_mode_speed_v21.py` → `run-hot-mode-v21/` (`hot_mode_speed_coarse.json`, `_fine.json`, `_power.json`).

**Setup.** Round 16's crossing-heat model (§26.2, `code/crossing_heat_v16.py`), unchanged except for where the heat
sits. There, the heat a star picked up a time t ago now sits a distance d = u t from it; here d = v_h t. The history is
round 16's: a straight pass, receding at 3,900 km/s since pericentre and approaching at 3,000 before it, impact
parameter 150 kpc; the round-12 law in the static distances. Two versions:
* **as proposed**: only the timing changes;
* **energy booked**: a glow emitted with the same power but spreading v_h/u times faster is v_h/u times thinner
  (energy density = power per area ÷ speed), so its heat weight carries u/v_h. This is the bookkeeping behind the
  law's own a = 2ℓ/u (§31.5). "As proposed" is the same as "energy booked" with the collision putting v_h/u times more
  power into the fast glow.

Masses inside 250 kpc in 10¹⁴ suns (targets 2.47–2.85 for the smaller half, 3.09–3.46 for the main one); gas
residuals measured 0.05 ± 0.06 (main) and 0.02 ± 0.06 (smaller half). Grid 22.5 kpc; the 600 km/s rows repeated at
15 kpc (in brackets).

| v_h (km/s) | as proposed: smaller half | main | peaks from galaxies (kpc) | gas residuals | energy booked: smaller half | main | gas residuals |
|---|---|---|---|---|---|---|---|
| no crossing heat | 1.372 (1.385) | 2.953 (2.969) | 16, 26 | 0.050, 0.048 | | | |
| 169.4 (round 16) | 1.586 | 3.101 | 15, 8 | 0.050, 0.018 | 1.586 | 3.101 | 0.050, 0.018 |
| 300 | 1.863 | 3.204 | 15, 8 | 0.055, 0.009 | 1.644 | 3.101 | 0.053, 0.009 |
| 450 | 2.241 | 3.303 | 17, 8 | 0.059, 0.069 | 1.711 | 3.100 | 0.055, 0.024 |
| **600** | **2.536 (2.559)** | **3.372 (3.392)** | 20, 9 (19, 8) | 0.051, 0.100 (0.052, 0.101) | 1.745 (1.759) | 3.096 (3.113) | 0.056, 0.036 |
| 800 | 2.833 | 3.427 | 24, 10 | 0.029, 0.120 | 1.755 | 3.089 | 0.057, 0.043 |
| 1200 | 3.195 | 3.527 | 37, 13 | 0.016, 0.138 | 1.720 | 3.078 | 0.059, 0.048 |

**How much power the fast glow needs** (energy booked, with the collision's power into the fast glow multiplied by m;
smaller half / main):

| v_h (km/s) | m = 2.5 | m = 3 | m = v_h/u (as proposed) |
|---|---|---|---|
| 600 | 2.256 / 3.274 | 2.396 / 3.324 | 2.536 / 3.372 (m = 3.54) |
| 800 | 2.273 / 3.256 | 2.414 / 3.302 | 2.833 / 3.427 (m = 4.72) |

The smaller half reaches the bottom of its range, 2.47, at m ≈ 3.2–3.3 at either speed: the requirement is set by the
power more than by the speed.

**Findings:**
* **Reproduced.** At 600 km/s both halves land inside their measured masses: 2.54 and 3.37 (2.56 and 3.39 on the finer
  grid), against the feedback's 2.50 and 3.44. Both stay in range for v_h ≈ 570–820 km/s. The lensing peaks stay on
  the galaxies (19–20 and 8–9 kpc). The price is lensing over the smaller half's gas: 0.10 against 0.02 ± 0.06 measured
  (1.3σ), growing with v_h (0.14 at 1,200 km/s).
* **Why it works.** Since the crossing (213 Myr ago in this history) the heat has spread to u t ≈ 37 kpc at 169 km/s,
  but to 130 kpc at 600 km/s, filling the 250 kpc aperture with the heat made when the two clusters were closest.
* **With the energy booked, most of the gain goes.** The smaller half reaches at most 1.76 (71% of the lower end of
  its range), against round 16's 1.59; the main half stays at 3.08–3.11. So the proposal works only if the collision
  puts at least about 3.3 times more power into the fast glow than the law's heat rule gives (the table above), or if
  this glow's pull follows its power per area rather than its energy density. That is a definite requirement for the
  microscopic model, as the feedback says for the ratio v_h/u itself.

### 32.2 The same rule in MACS J0025 and El Gordo

`code/hot_mode_collisions_v21.py` → `run-hot-mode-v21/hot_mode_collisions_v21.json`. Round 16's crossing heat had only
been applied to the Bullet; the suite's far collisions have none. Here it is added to the two far collisions with a
clean two-body history (Abell 520's clumps have no agreed history): straight passes with the suite's times since
pericentre (MACS J0025 0.3 Gyr, from its shock fronts, §22.1; El Gordo 0.46 Gyr, outgoing, Ng et al. 2015), receding at
the galaxies' separation over that time, approaching at 2,000 km/s (MACS J0025, Bradač et al. 2008) and 2,400 km/s (El
Gordo's pericentre speed), impact parameter 150 kpc (not measured). Grid and apertures as in the suite (the "none"
rows reproduce it).

In the static distances the histories are: MACS J0025 now 772 kpc apart, receding at 2,420 km/s; El Gordo 1,246 kpc
apart, receding at 2,590 km/s (faster than its pericentre speed: the static distances make it larger than in the
papers' convention, a caveat on the straight-pass history).

| cluster: lensing masses (10¹⁴ suns) | no crossing heat (the suite) | round 16's rule (v_h = u) | v_h = 600, as proposed | v_h = 600, energy booked |
|---|---|---|---|---|
| MACS J0025 SE, inside 300 kpc (3.64, +1.46/−2.48) | 2.05 (z −0.64) | 2.24 (−0.57) | 2.58 (−0.43) | 2.23 (−0.57) |
| MACS J0025 NW, inside 300 kpc (3.79, +0.73/−2.04) | 1.85 (−0.95) | 2.03 (−0.86) | 2.50 (−0.63) | 2.06 (−0.85) |
| MACS J0025 lensing peaks from their galaxies (kpc; gas 530 and 244 kpc away) | 24, 79 | 13, 18 | 18, 25 | 21, 44 |
| El Gordo inside 500 kpc (9.45 ± 12%) | 8.42 (−0.91) | 8.86 (−0.52) | 9.75 (+0.27) | 8.84 (−0.54) |
| El Gordo inside 1,000 kpc (24.3 ± 12%) | 21.4 (−0.98) | 22.1 (−0.75) | 23.7 (−0.19) | 22.1 (−0.74) |
| crossing heat ÷ all heat (mid-plane), MACS J0025 / El Gordo | 0 | 0.26 / 0.15 | 0.53 / 0.37 | 0.24 / 0.14 |

* **The heat of crossing helps both clusters, already at the law's own speed.** Both move toward their measured
  masses, and MACS J0025's NW lensing peak, 79 kpc from its galaxies in the suite, moves to 18 kpc.
* **The faster glow as proposed passes here too:** El Gordo lands on its measured masses (z +0.27 and −0.19) and
  MACS J0025 rises further (−0.43 and −0.63), with the peaks on the galaxies. Nothing overshoots.
* **With the energy booked,** the results are those of round 16's rule, as in the Bullet.
* The 72-collision stack (Harvey et al. 2015) has no individual histories and is not rerun here; in all three clusters
  the extra lensing rides with the galaxies, which is what the stack measures.

### 32.3 The distance law: η = 1/2 fits, the photons must be conserved, and the stretch needs change over time

`code/distance_eta_v21.py` → `run-distance-eta-v21/distance_eta_v21.json`: round 12's Pantheon+ reduction
(`code/sn_scale_v12.py`: Cepheid-calibrated magnitudes, the full STAT+SYS covariance), 1,365 supernovae at z > 0.023.

| | all 1,365 | z < 0.15 (490) | z > 0.15 (875) |
|---|---|---|---|
| η = 0 (the adopted law): Δχ² | +100.9 | +7.4 | +27.5 |
| best η (Δχ² ≤ 1 range) | **0.449** (0.41–0.49) | 0.56 (0.35–0.77) | 0.38 (0.31–0.46) |
| η = 1/2: Δχ² | **+1.09** | +0.08 | +2.15 |
| the scale at η = 1/2, as H₀ (km/s/Mpc) | 73.2 | 73.0 | 73.2 |
| pulse energy conserved (D_L = (1 + z)^½ D*), η refitted: Δχ² | **+59.9** (η = 2.4) | −0.2 | +12.1 |
| the power of (1 + z) left free: power, η, Δχ² | 0.89, 0.80, −2.5 | – | 0.80, 1.26, −3.5 |
| for scale only, flat ΛCDM fitted the same way: its χ² minus ours | −2.1 (Ω_m 0.33) | +0.2 | −2.4 |

**Findings:**
* **η = 1/2 is allowed** (Δχ² = 1.1 against the best, η = 0.45), and b = z/(1 + z) is exactly the solution of
  db/dy = 1 − b from b = 0 (checked to 4 × 10⁻¹²): a state relaxing at the redshift's own rate. η = 1/2 is worth
  deriving rather than fitting.
* **It fits as well as the standard model's distances:** χ² 1210.4 against 1208.3, with two fitted numbers each (the
  standard model is shown only as a yardstick; it is excluded as an answer). To second order,
  α D* = z − (1 − η) z²/2: the expanding-universe form with deceleration parameter q₀ = −η. The standard fit has
  q₀ = −0.51; η = 1/2 is q₀ = −1/2; the adopted law (η = 0) is the coasting form, q₀ = 0. The supernovae measure the
  same bend in both readings and cannot choose between them: the path factor is our law's version of what the standard
  reading calls cosmic acceleration.
* **The operation must conserve the photons, not the energy.** a_out(ω) = √s a_in(sω) keeps each pulse's energy and
  multiplies the number of photons by 1 + z, so the flux is dimmed by the stretching alone: D_L = (1 + z)^½ D*. The
  data reject that by Δχ² = 60 (η would have to be 2.4). What fits is a_out(ω) = a_in(sω): each photon arrives with
  1/(1 + z) of its energy, stretched by 1 + z, so the medium must take up the fraction z/(1 + z) of the light's
  energy. (Left free, the power of (1 + z) comes out 0.89 with η 0.80: Δχ² −2.5 for one more number.)
* **The stretch needs something that changes in time.** If nothing along the way changes with time, pulses sent one
  second apart arrive one second apart, and a linear medium cannot change light's frequency at all; a gravitational
  redshift between two points is a blueshift the other way. The measured stretch is exact (b = 1.003 ± 0.011 in
  (1 + z)^b, 1,504 supernovae; White et al. 2024), so the travel time from a fixed galaxy must grow by z seconds every
  second: the optical length of the path grows at α c ≈ 7.5 × 10⁻¹¹ per year. The laboratory "time lens" the feedback
  cites works exactly this way: an electronic modulator changes the medium while the pulse passes. If that growth also
  happened inside laboratories, optical cavities would drift against atomic clocks by 2.4 × 10⁻¹⁸ per second; the
  earlier-era comparison (`nature_tests/report.md`) found that 9–47 times the drifts measured in four silicon cavities
  (Lee et al. 2025), a tension rather than an exclusion, since each instrument's own drift is not independently known.
  So the changing property has to be confined to the space between galaxies (for example, a property of the
  intergalactic companion), or be shared by rulers and clocks, in which case it cannot be told apart from expansion.

**The suite** (full tier, registered comparisons, a and u refitted; reports in `run-distance-eta-v21/suite/`):

| candidate | pass / close / fail | what moves |
|---|---|---|
| round 12, adopted (fixed geometry D_A = D; η = 0; H₀-like 70.9) | 59 / 11 / 7 | |
| `dist_fixed_eta_r12`: η = 1/2 with its scale (H₀-like 73.2), fixed geometry | **62 / 8 / 7** | KiDS's all, red and disk lenses pass (+0.035, +0.048, +0.036 dex, from +0.065, +0.078, +0.066); El Gordo's NW galaxy speeds pass; lensing masses rise in MACS J0025 (2.43 and 2.17 × 10¹⁴), El Gordo (1.09 and 1.07 of its two aperture masses) and the Bullet's smaller half (0.61 of the lower end of its range, from 0.56); Mistele's ellipticals 2.55 → 3.71 (fail); Abell 520 P2 → close |
| `dist_metric_eta_r12`: the proposal (η = 1/2, H₀-like 73.2, D_A = D*/(1 + z)) | **62 / 8 / 7** | KiDS's all, red and disk lenses pass (+0.027, +0.042, +0.029 dex); Mistele's ellipticals 4.03 (fail); the Bullet's smaller half 0.48 of the lower end of its range; El Gordo 0.89 and 0.86 of its aperture masses |
| `dist_metric_r12`: the metric geometry alone (η = 0, H₀-like 70.9) | 57 / 13 / 7 | MACS J0025's galaxy speeds fail (631 km/s against 835 ± 59); El Gordo falls to 0.73 and 0.70 of its aperture masses (close); the Bullet's smaller half 0.45; KiDS's red lenses fail → close |

* **The path factor helps; the metric geometry alone does not.** With η = 1/2 (and the scale its fit gives), the
  lensing level of 259,000 KiDS galaxies, 16% above the law in our distances until now (§30), falls to 6–8%, inside
  the test's tolerance. The geometry D_A = D/(1 + z) shrinks the far clusters (El Gordo at z = 0.87, MACS J0025 at
  0.59) and costs them mass and galaxy speed. The price of η = 1/2 in both geometries is Mistele et al.'s ellipticals
  (lensing speeds 50–300 kpc), which fail.
* All three stay registered comparisons (`regression/candidates/dist_*.json`); none is adopted. SPARC's Hubble-flow
  distances take the new scale but not the path factor (under 1% at their redshifts), and SLACS keeps its own
  conversion.

### 32.4 The relativistic action: a good skeleton, and what has to change

`code/action_checks_v21.py` → `run-action-checks-v21/action_checks_v21.json` (weak field; a Milky-Way-sized mass,
6 × 10¹⁰ suns, where a number is needed).

**What it gets right.** The split, with the occupation I setting the size of the pull and the flow s^μ its
direction, is what the clusters' hot shells asked for (§29), and deriving everything from one action is the right
way to settle where momentum goes (§10.6).

**The checks:**
1. **The sign.** With s pointing outward, the constraint as written, s·∇χ = −f√(aI), makes the extra pull repel: at
   10 kpc a Milky-Way-sized mass would hold a circular speed of 105 km/s instead of Newton's 161, and beyond about
   14 kpc the net pull points outward. With s·∇χ = +f√(aI): 202 km/s at 10 kpc, 168 at 40. A one-sign fix.
2. **The constraint field λ is harmless as mass, but it acts instantly.** Varying χ, with matter coupled to Φ_N + χ,
   gives div(λ s) = −8πGρ, so λ = −2 G M(r)/r² around a round source. Its energy between 1 kpc and 1 Mpc is
   1.2 × 10⁻⁶ of the source's mass (of order (v_f/c)² ln(r₂/r₁)): no hidden mass. But χ and λ are fixed by
   integrating along each flow line through space, so a change in I anywhere along a line resets χ along all of it at
   once, which needs a preferred frame. Imposing the constraint along the companion's world lines (moving at u)
   instead of its direction in space would make changes travel at u: the law's memory (§12, §26.2).
3. **I and s must be the companion's own fields.** If I contains the heat term S computed from the matter's positions
   and speeds, varying the matter adds a force on warm matter of the same form as round 2's reaction (§8.2): 40–77% as
   strong (0.65–0.75 where the heat matches or exceeds the cold pull, as in clusters). Round 3 excluded round 2's push:
   32–73% of gravity on cluster stars (§10.6). With I and s as independent fields that matter feeds by emitting evenly
   in all directions, the reaction lands on the companion, as §10.6 concluded.
4. **Light and gravitational waves.** A single field χ that enters matter's metric as an overall rescaling bends no
   light (light's paths ignore such a rescaling). Getting χ into the spatial part with the same sign needs a unit
   time-like vector field: the construction of TeVeS and AeST, which RULES.md excludes. And gravitational waves must
   feel χ as light does. From the Milky Way's χ alone (the law's constants, the companion's 2.25 Mpc reach), light from
   GW170817 would have arrived about 3.6 years after the waves if the waves did not feel it; they arrived within 1.7 s
   (Boran et al. 2018 make the same argument against dark-matter emulators). So χ has to live in the metric itself:
   the metric's own field equation must change, not only matter's coupling to it.
5. **√(aI) is inserted.** Round 20 (§31.5) found where it can come from: the companion's energy balance, with energy
   density |∇Φ|²/8πG carried at u and fed at ℓ, gives the deep equation. An action built from that balance would
   derive the square root instead of assuming it.

### 32.5 Where round 21 leaves things

* **The Bullet's smaller half has a candidate fix with a definite physical requirement.** A faster glow from the
  collision, v_h ≈ 570–820 km/s, puts both halves in range and keeps the lensing peaks on the galaxies, and the same
  heat of crossing moves MACS J0025 and El Gordo toward their measured masses (§32.2). It needs the collision to put
  at least about 3.3 times as much power into that glow as the law's heat rule gives (§32.1), or the glow's pull to
  follow its power per area rather than its energy density. Next: derive v_h/u and that power from the medium's
  dispersion (the feedback's own programme), with every watt booked; then the 72-collision stack with the heat of
  crossing.
* **The distance law has a strong candidate refinement.** The path factor with η = 1/2 has a derivable form (a state
  relaxing at the redshift's own rate), fits the supernovae as well as the standard model, and brings KiDS's lensing
  level into the suite's tolerance. Its conditions: the photons are conserved (the intergalactic medium takes up
  z/(1 + z) of the light's energy), and something in intergalactic space changes slowly with time (7.5 × 10⁻¹¹ per
  year) without reaching laboratories. Next: derive η = 1/2 and the photon-conserving stretch from one intergalactic
  process, and say what in the intergalactic companion changes; test the geometry (fixed or metric: angular sizes
  differ by 1 + z) with standard rulers and surface brightness.
* **The relativistic action keeps its skeleton with four changes:** the sign; the constraint along the companion's
  world lines; I and s as the companion's own fields, fed by emission even in all directions; and χ in the metric that
  gravitational waves feel. The last is the hardest and decides whether a relativistic version can exist outside the
  MOND family's constructions.
* **Unchanged:** the adopted law, its constants, the locked forecasts and the frozen lensing test (§29.7).

## 33. Round 22, 26 September 2026: the fast glow from collisions, taken as far as the model allows

The request: "proceed with the next step", the first item of round 21's list (§32.5): derive the fast glow's speed and
power from the medium, with every watt booked, then test it on the collision stack. Four steps:
1. **Could the glow's speed be a matter of frame?** If the glow a collision makes is left partly behind in the other
   system's frame, like a boat's wake, it would spread away from the stars at hundreds of km/s with no new speed and no
   extra power (§33.1).
2. **The energy bill** of the heat of crossing, and of the fast glow's extra power (§33.2).
3. **What an evenly spreading fast glow requires of the medium** (§33.3).
4. **The collision stack:** which way each version moves the lensing, the quantity the 72-collision stack measures
   (§33.4).

Nothing in the adopted law, its constants or the locked forecasts changes.

### 33.1 A glow left behind (a wake) does not work

`code/crossing_frame_v22.py` → `run-crossing-frame-v22/` (`frame_a.json` … `frame_d.json`).

**Setup.** Round 16's heat of crossing and history, with one change: the glow a star emitted a time t ago keeps only a
fraction μ of the star's motion relative to the other system. It still travels at u in its own frame, so the sphere it
fills has radius u t, but the sphere's centre has fallen behind the star by (1 − μ) × the path the star has travelled
relative to the other system since. μ = 1 is round 16; μ = 0 is a glow that stays where it was emitted. Each time slice
carries exactly round 16's energy (Gaussian shells of width 11 kpc on the 22.5 kpc grid); only where it sits changes.
Seen from the star, the glow drifts backwards at up to (1 − μ) v + u: 560–950 km/s for μ = 0.9–0.8. The check at μ = 1
against round 16's own kernel: 1.603 / 3.104 against 1.586 / 3.101.

Masses inside 250 kpc in 10¹⁴ suns (targets 2.47–2.85 for the smaller half, 3.09–3.46 for the main one); peak offsets
along the collision axis, positive toward each system's own gas (the suite allows a quarter of the galaxy–gas
separation: 63 and 58 kpc); gas residuals measured 0.05 ± 0.06 and 0.02 ± 0.06. Grid 22.5 kpc.

| μ (share of the star's motion the glow keeps) | glow drifts back from the star at up to (km/s) | smaller half | main | peaks toward the gas: main, smaller (kpc) | gas residuals |
|---|---|---|---|---|---|
| no heat of crossing | – | 1.372 | 2.953 | +16, +25 | 0.050, 0.048 |
| 1 (round 16) | 169 | 1.603 | 3.104 | +14, +7 | 0.051, 0.013 |
| 0.9 | 560 | 1.675 | 3.133 | +27, **+93** | 0.054, **0.154** |
| 0.85 | 760 | 1.631 | 3.142 | +21, **+125** | 0.056, **0.174** |
| 0.8 | 950 | 1.590 | 3.143 | +20, **+189** | 0.064, **0.160** |
| 0.7 | 1,340 | 1.534 | 3.151 | +29, +45 | 0.080, 0.125 |
| 0.5 | 2,120 | 1.505 | 3.147 | +22, +35 | 0.078, 0.089 |
| 0 (a wake) | 4,070 | 1.486 | 3.082 | +17, +29 | 0.066, 0.064 |
| for comparison: the even fast glow of round 21 (600 km/s, as proposed) | 600 | 2.536 | 3.372 | +19, +9 | 0.051, 0.100 |

**Findings:**
* **A glow left behind never lifts the smaller half.** Its mass stays at 1.49–1.68, at most 68% of the lower end of
  its range, against 2.54 for the evenly spreading glow at the same speed. The energy is the same; it is simply
  deposited behind the stars instead of around them.
* **It drags the lensing onto the gas.** For μ = 0.8–0.9, where the drift matches round 21's speeds, the smaller half's
  lensing peak moves 93–189 kpc toward its gas (the suite allows 58), and the lensing over that gas rises to 0.15–0.17
  (2.2–2.6σ above the measured 0.02 ± 0.06). A glow left further behind (μ ≤ 0.7) moves the peak less, because the
  glow ends up beyond the gas, near the other cluster, but it adds nothing to the smaller half either.
* **So the fast glow cannot be a matter of frame.** It has to spread evenly around the stars that made it: a
  genuinely faster mode of the companion's medium, with its own power (§33.2, §33.3).

### 33.2 The energy bill

`code/crossing_energy_v22.py` → `run-crossing-frame-v22/crossing_energy_v22.json`.

The cold companion carries ℓ = a u / 2 = 5.3 × 10⁻⁶ W per kilogram of matter; warm matter carries (1 + k) ℓ. Integrated over
round 16's crossing history (the other system's flow reaches the stars for the last 950 Myr), the heat of crossing
costs, per kilogram of stars:

| stars | largest k | ∫ k dt (Myr) | energy (J/kg) | peak power (W/kg) |
|---|---|---|---|---|
| smaller half, inner 60 kpc | 214 | 76,700 | 1.3 × 10¹³ | 1.1 × 10⁻³ |
| smaller half, 60–180 kpc | 313 | 129,000 | 2.2 × 10¹³ | 1.7 × 10⁻³ |
| main cluster, inner 60 kpc | 89 | 18,000 | 3.0 × 10¹² | 4.7 × 10⁻⁴ |

**For scale,** per kilogram of the smaller half's inner stars:
* the kinetic energy of the relative motion (centre-of-mass frame) is 7.3 × 10¹² J/kg at 3,900 km/s: **round 16's glow
  already carries 1.8 times that**, and the fast glow's 3.3 times the power would carry 5.8 times it. The glow cannot be
  paid for by slowing the collision down (it would stop the smaller cluster); it has to come from energy stored inside
  matter, as round 16's quiet store assumes;
* the fast glow's total, 4.3 × 10¹³ J/kg, is 4.7 × 10⁻⁴ of the stars' rest energy, and a quarter (0.24) of what the
  law's settled heat already asks of a cluster galaxy's stars over 10 Gyr (k = 3σ²/u² at σ = 1,000 km/s:
  1.8 × 10¹⁴ J/kg). So the extra power is large, but in line with what the heat term already requires of warm matter;
* at its peak, round 16's glow is 5.9 times a Sun-like star's light output per kilogram (19 times for the fast glow).
  It is dark (companion, not light), but the store that feeds it must hold at least this energy: a quantitative
  target for the store's model, which has been assumed since round 16.

### 33.3 What an evenly spreading fast glow requires of the medium

`code/fast_glow_medium_v22.py` → `run-crossing-frame-v22/fast_glow_medium_v22.json`.

**Where a hot glow would live.** Round 10 (§20.1, condition 3) requires the companion's crests to move much more
slowly than its energy, v_phase ≲ 0.05 u. A star crossing another system's flow at speed w therefore meets that flow's
crests at the Doppler-shifted frequency ω′ = ω₀ |1 − w_r/v_phase|: the collision drives each star hundreds to
thousands of times faster than the companion's own frequency ω₀. That is the natural home of a distinct "hot" glow,
and its speed is the medium's group velocity at ω′.

| relative speed | crest speed | ω′/ω₀ | needed rise of the group velocity, as a power of frequency | single power law ω ∝ kⁿ with v_g/v_phase = 20: v_g rises by | two-regime medium: curvature needed |
|---|---|---|---|---|---|
| 3,900 km/s | 0.05 u | 459 | ω^0.21 | ×338 | β = 0.12 u/k₀ |
| 3,900 km/s | 0.01 u | 2,301 | ω^0.16 | ×1,563 | 0.12 u/k₀ |
| 3,000 km/s | 0.05 u | 353 | ω^0.22 | ×263 | 0.16 u/k₀ |

* **A single power-law medium cannot do it.** Such a medium has v_g/v_phase = n at every frequency, so slow crests at
  ω₀ (n ≥ 20) make the group velocity rise several hundredfold by ω′: the hot glow would race away at 45,000–265,000
  km/s and, with its energy booked, be far too thin to matter. The 3.5-fold rise needs n ≈ 1.2–1.3, which contradicts
  the slow crests.
* **A two-regime medium can.** A frequency that crosses zero at a finite wavenumber k₀ and rises as ω = u q + β q²
  above it (q = k − k₀) has slow crests just above k₀ (the carrier at q ≈ 0.05 k₀) and a group velocity u + 2βq that
  reaches 3.5 u at ω′ for β ≈ 0.12–0.16 u/k₀.
* So "v_h/u ≈ 3.5" becomes a concrete statement about the medium's dispersion, which the microscopic model has to
  produce, together with the power (§33.2). It is not yet derived.

### 33.4 The collision stack: which way the heat moves the lensing

`code/crossing_offsets_v22.py` → `run-crossing-frame-v22/crossing_offsets_v22.json`, with the frame runs of §33.1.

The stack's own family (regression/t_collisions.py) holds the galaxies fixed and moves only the gas, so a heat of
crossing with one history adds the same lensing to every member and cannot change its β. What changes β is heat that
moves toward the gas as time goes on. Both the gas lag and any such drift grow with the time since the crossing, so on
the Bullet the change in β is estimated as (peak offset with the heat − without) / (galaxy–gas separation):

| version | peaks toward the gas: main, smaller (kpc) | change in β: main, smaller, mean |
|---|---|---|
| no heat of crossing | +15.9, +25.1 | 0 |
| round 16's rule (u) | +14.3, +8.1 | −0.006, −0.073, **−0.04** |
| even fast glow, 600 km/s, as proposed | +19.2, +8.6 | +0.013, −0.071, **−0.03** |
| even fast glow, 600 km/s, energy booked | +17.2, +15.7 | +0.005, −0.040, **−0.02** |
| glow left behind, μ = 0.9 / 0.85 / 0.8 | +27 / +21 / +20, +93 / +125 / +189 | **+0.17 / +0.22 / +0.36** |
| glow left behind, μ = 0.5 / 0 | +22 / +17, +35 / +29 | +0.03 / +0.01 |

* **The evenly spreading glow passes the stack.** The suite's β without the heat of crossing is +0.02; the even glow
  moves it by −0.02 to −0.04, to between 0.00 and −0.02, against Harvey et al.'s −0.04 ± 0.07. Heat centred on the
  galaxies pulls the smaller half's peak onto them (from 25 to 8 kpc).
* **The glow left behind fails it** where it matches round 21's speeds: β +0.19 to +0.38, 3–6σ from the measurement.
* This is an estimate from the Bullet's stage, not a rerun of the twelve-member family, which would need collision
  histories for every member.

### 33.5 Where round 22 leaves things

* **The fast glow is now a sharply defined target.** It must spread evenly around the stars that made it (a glow left
  behind drags the lensing onto the gas and adds no mass, §33.1), at about 3.5 times the settled companion's speed,
  with about 3.3 times the power the heat rule gives (§32.1). That power is 5.8 times the collision's energy of motion,
  so it must come from matter's internal store; the store must hold at least 4 × 10¹³ J/kg for the Bullet's stars,
  a quarter of what the settled heat of a cluster galaxy already costs over 10 Gyr (§33.2).
* **In the medium,** the collision drives each star at the frequency at which it crosses the other system's slow
  crests, 350–2,300 times the companion's own. An even glow at 600 km/s needs the group velocity to rise about 3.5-fold
  up to there: impossible for a single power-law medium with slow crests (which would make the glow hundreds of times
  faster and far too thin), possible for a two-regime one (§33.3).
* **It passes the stack** (β between 0.00 and −0.02, measured −0.04 ± 0.07), MACS J0025 and El Gordo (§32.2).
* **Next:** a microscopic medium with that two-regime dispersion, driven at the crest-crossing frequency, to compute the
  hot glow's speed and power from first principles with every watt booked; and a store model with the capacity §33.2
  requires. The rest of round 21's list stands (§32.5).
* **Unchanged:** the adopted law, its constants, the locked forecasts and the frozen lensing test (§29.7).

## 34. Round 23, 26 September 2026: when the two clusters' companions interfere

The owner asked: "I don't understand the problem. Can the situation be different if there are two different sources
and they create an interference pattern?"

**The problem, in plain terms.** Lensing weighs each half of the Bullet Cluster. Inside 250 kpc of the smaller half it
finds 2.47–2.85 × 10¹⁴ suns (in our distances). The law gives 1.37 from the visible matter and its settled glow, and
1.59 once the glow its stars made while crossing the big cluster is added (round 16). That extra glow stays close to
the stars: it spreads at the companion's 169 km/s, so in the ~210 Myr since the crossing it has reached only ~37 kpc,
while the lensing counts everything out to 250 kpc. A glow spreading at ~600 km/s would fill the circle and fit
(round 21), but a faster glow with the same energy is thinner; it fits only if the collision puts about 3.3 times more
energy into it (round 21; round 22 showed it must spread evenly, not be left behind).

### 34.1 What interference can and cannot do

* **A pattern of bright and dark bands adds no energy.** It only moves it around: the bright bands are paid for by the
  dark ones.
* **For the pull, a pattern that is averaged over is slightly worse than none.** The law's pull follows the
  companion's height (amplitude). Two equal waves with a random relative phase average to 4/π = 1.27 times one wave's
  height, against √2 = 1.41 when they add without interfering: 10% less (5% at a height ratio of 0.5, 2% at 0.3).
  And the clusters pass each other at 3,000–3,900 km/s, so any pattern sweeps across every star hundreds of times per
  companion cycle (the crest-crossing rate, §33.3): for the pull, it is averaged over.
* **What interference adds is flicker.** At a star sitting in both companions, the combined intensity rises and
  falls by 2√(I_other I_own) around the mean I_other + I_own, at the crest-crossing rate. As a fraction of the mean
  that is m = 2√(I_other I_own)/(I_other + I_own). Round 16's heat of crossing used the other companion's share,
  p = I_other/(I_other + I_own), as if the two added without interfering. The flicker is much larger when the other
  companion is faint, because interference grows with the fainter wave's height, not its intensity (the principle of
  heterodyne detection): a companion with 1% of the intensity makes the total flicker by 20%. For the Bullet's smaller
  half, m/p = 2.3 at closest approach (p = 0.44), 3.4 now (p = 0.26) and 6.6 during the approach (p = 0.085).
* **The flicker is the shaking that drives the heat of crossing.** Round 22 found that a crossing star is driven at
  the rate it crosses the other companion's crests; the flicker is that drive. So the heat of crossing is weighted
  here by the flicker's depth m instead of the share p. Nothing else changes, and no constant is refitted.

### 34.2 The Bullet Cluster with the heat of crossing driven by the flicker

`code/beat_heat_v23.py` → `run-beat-heat-v23/bullet_*.json` (grid 22.5 kpc; energy booked whenever the glow spreads
faster than u).

| heat of crossing | smaller half (target 2.47–2.85) | main (3.09–3.46) | peaks from galaxies (kpc) | toward gas | gas residuals (0.05 ± 0.06, 0.02 ± 0.06) |
|---|---|---|---|---|---|
| share p, spreading at u (round 16) | 1.586 | 3.101 | 15, 8 | +14, +8 | 0.050, 0.018 |
| flicker m, spreading at u | 1.695 | 3.194 | 8, 5 | +8, +4 | 0.050, −0.017 |
| flicker m, 300 km/s | 1.874 | 3.203 | 9, 6 | +9, +6 | 0.050, 0.016 |
| flicker m, 600 km/s | 1.995 | 3.190 | 11, 10 | +10, +10 | 0.050, 0.073 |
| flicker m, 600 km/s, power × 1.5 | 2.216 | 3.300 | 9, 8 | +9, +8 | 0.048, 0.088 |
| flicker m, 600 km/s, power × 2 | 2.401 | 3.404 | 8, 7 | +8, +6 | 0.046, 0.096 |
| for comparison: share p, 600 km/s, energy booked (round 21) | 1.745 | 3.096 | 18, 16 | +17, +16 | 0.056, 0.036 |

**Findings:**
* **Interference helps, on its own terms.** With no extra energy assumed, the smaller half rises from 1.59 to 1.70,
  and to 2.00 with the faster glow (against 1.75 without interference): about a third of the remaining gap.
* **It cuts the extra power the fast glow needs from about 3.3 times to about 2.2 times** (2.47 is reached at
  × 2.2 by interpolation; the main cluster is then at the top of its range, 3.44).
* **The lensing sits on the galaxies** (5–11 kpc; the suite allows 58–63) and the gas stays clean, because the
  flicker is strongest right around the stars.

### 34.3 MACS J0025 and El Gordo

`code/beat_heat_far_v23.py` → `run-beat-heat-v23/far_*.json` (the models and histories of §32.2).

| cluster: lensing masses (10¹⁴ suns) | share p, at u (round 16's rule) | flicker m, at u | flicker m, 600 km/s (energy booked) |
|---|---|---|---|
| MACS J0025 SE, inside 300 kpc (3.64, +1.46/−2.48) | 2.24 (z −0.57) | 2.23 (−0.57) | 2.23 (−0.57) |
| MACS J0025 NW, inside 300 kpc (3.79, +0.73/−2.04) | 2.03 (−0.86) | 2.00 (−0.88) | 2.02 (−0.87) |
| MACS J0025 lensing peaks from their galaxies (kpc) | 13, 18 | 4, 9 | 12, 23 |
| El Gordo inside 500 kpc (9.45 ± 12%) | 8.86 (−0.52) | 8.68 (−0.68) | 8.61 (−0.74) |
| El Gordo inside 1,000 kpc (24.3 ± 12%) | 22.1 (−0.75) | 21.6 (−0.91) | 21.7 (−0.90) |

* **In the two near-equal collisions the flicker changes little.** Both stay within their errors, slightly lower than
  with the share rule. Where one companion dominates (each system's outer stars sit mostly in the other's glow,
  p > 0.8), the flicker is shallow (m < p), which offsets the gain for the inner stars. MACS J0025's lensing peaks move
  onto their galaxies (4 and 9 kpc, from 13 and 18).
* **So the flicker helps most in unequal collisions**, where a small system sits in a big one's companion: that is the
  Bullet's smaller half. This is a testable pattern: in collisions of unequal clusters, the smaller one should carry
  relatively more extra lensing than the same rule gives near-equal pairs.

### 34.4 Where round 23 leaves things

* **Interference is part of the answer.** Driving the heat of crossing by the two companions' flicker instead of
  their shares supplies about a third of the Bullet's missing lensing on its own terms, keeps the lensing on the
  galaxies, leaves the near-equal collisions within their errors, and cuts the extra power the fast glow needs from
  3.3 times to about 2.2 times.
* **What remains:** that factor of about 2.2 in power, from matter's internal store (§33.2), and the faster spreading
  from the medium's two-regime dispersion (§33.3). The flicker rule itself (heat ∝ m) is motivated, not derived: the
  store's response to a modulated companion is the next calculation for the matter model.
* **A registered comparison,** not adopted: the suite's Bullet has no heat of crossing, so none of rounds 16–23 enters
  its tally.
* **Unchanged:** the adopted law, its constants, the locked forecasts and the frozen lensing test (§29.7).
