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
