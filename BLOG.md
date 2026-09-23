# Gravity that streams

*Hot-companion gravity: one law for rotating galaxies, bending light, galaxy clusters and
colliding clusters, with no dark matter and no expanding universe.*

**Rewritten from scratch on 23 September 2026 (rev 12).**
* Every number below is computed from public data by a script in this repository, named
  where the number appears (§11).
* The earlier notebook (revisions 1–11), with all its retracted and retired claims left
  visible, is archived in
  [research_work/blog-archive/BLOG-notebook-rev11.md](research_work/blog-archive/BLOG-notebook-rev11.md).
* The technical record is in
  [research_work/results/hot-companion/README.md](research_work/results/hot-companion/README.md).

---

## In one minute

**The puzzle.** Galaxies and clusters of galaxies pull harder than the matter we can see
should allow:
* Stars at a galaxy's edge orbit too fast.
* Light passing a galaxy bends too much.
* Clusters hold their hot gas with 4 to 27 times the pull their visible matter provides.

When two clusters collide, the extra pull follows the galaxies, not the gas, even though the
gas outweighs the galaxies ten to fifty times. The usual answer is invisible "dark matter".

**Our answer, in plain words.**
1. **Every piece of ordinary matter feeds a faint, streaming companion to its gravity.**
2. **Orderly matter feeds it in step.** Contributions from opposite sides cancel, just as
   Newton's pulls do. Far from a galaxy, what survives falls off slowly enough to keep the
   outer stars orbiting fast.
3. **Matter whose parts move randomly and freely, like stars in elliptical galaxies or
   galaxies in clusters, feeds it out of step.** Nothing cancels, so it pulls much harder.
   That is why clusters need so much extra pull.
4. **Gas does not count as "hot" in this sense.** Its particles collide so often that the
   companion sees them as sitting still. That is why, in a collision, the extra pull stays
   with the galaxies and leaves the gas behind.
5. **Strong gravity holds the companion back,** so the Solar System feels nothing.

**The scoreboard.**

| Test | Ours | MOND | Dark matter |
|---|---|---|---|
| 149 galaxies' rotation speeds (typical miss) | **15.9 km/s** | 16.1 km/s | 7.5 km/s |
| 12 galaxy clusters' masses (typical miss) | **25%** | 2.9× | 11% |
| Bullet Cluster: lensing on the galaxies, not the gas | **yes** | no | yes |
| 72 colliding clusters: lensing stays with the galaxies | **yes** | no | yes |
| Ellipticals bend light 0.17–0.27 dex more than spirals | **0.17–0.27** | no difference | yes, with tuned haloes |
| Solar System | **no extra pull** | small extra pull | no extra pull |
| Adjustable numbers | **3 in total** | 1 | 2 per galaxy or cluster (~320) |

**What is still open:**
* the strength of the lensing around the Bullet Cluster's smaller half;
* a full relativistic version of the law.

Both are in §8, with what could close them.

---

## 1. The problem, as plainly as we can put it

Four measurements have to be explained together:

1. **Rotation.** In a spiral galaxy, stars far from the centre should orbit slowly, the way
   Neptune orbits the Sun more slowly than Mercury. They don't; the speeds stay flat. Across
   149 well-measured galaxies (the SPARC catalogue, 3,150 speed measurements), Newton's law
   applied to the visible stars and gas misses the measured speeds by 45.6 km/s on average.
2. **Light bending (lensing).** Galaxies and clusters bend the light of objects behind them
   more than their visible matter can.
3. **Clusters.** In 12 well-measured galaxy clusters (the X-COP sample), the pull needed to
   hold the hot gas in place is up to 27 times what the visible matter supplies near the
   centre, and 4 times at the edge. Newton misses these masses by a factor of about 9.
4. **Collisions.** In the Bullet Cluster, two clusters passed through each other about 150
   million years ago. Their galaxies sailed through. Their gas, which outweighs the galaxies
   ten to fifty times, collided and was left behind in the middle. The lensing, which traces
   the pull, sits on the galaxies. A stack of 72 such collisions says the same.

There are three standard responses:
* **Dark matter.** Add invisible mass, tuned galaxy by galaxy, until the sums work. It fits
  well, but it needs about two adjustable numbers per object, and the particle has never been
  detected.
* **MOND.** Change the law of gravity at low accelerations. It predicts galaxy rotation well
  with one constant, but it fails in clusters by a factor of about 3, and it cannot put the
  collision's lensing on the galaxies.
* **Newton with visible matter only.** It fails everywhere outside the Solar System.

## 2. The ground rules we work under

The owner of this project has set four rules (recorded in [RULES.md](RULES.md)):
* no dark matter;
* no MOND or anything derived from it as an input;
* no plain Newtonian gravity as the answer;
* no expanding universe.

Every candidate formula goes through an automatic check,
[research_work/tools/formula_guard.py](research_work/tools/formula_guard.py). It asks three
questions:
* **Is it secretly MOND?** For an isolated cold mass, does it reduce to a function of the
  Newtonian pull alone? If so, how close is it to each published MOND formula? Across real
  galaxies, how much of its prediction depends on anything else?
* **Is it secretly Newton?** Is its boost just a constant, which is Newton with rescaled masses?
* **Is it secretly dark matter?** Does it need extra, adjustable, invisible mass?

We also work in a fixed order: first principles, then a toy model, then real data. When a
test fails, the first question is what would have to change to make it pass.

## 3. The idea, built from first principles

### 3.1 Matter feeds a companion to its gravity

Suppose every kilogram of ordinary matter steadily feeds a tiny amount of energy, ℓ watts
per kilogram, into a companion to its gravity. The energy streams outward at a speed u. This
is the project's founding intuition, that gravity is produced and carried rather than simply
attached. It is the only new ingredient; everything else follows from how such a stream
adds up.

### 3.2 Orderly matter adds up the way Newton does

Energy flowing away from a source spreads over a sphere, so its flux falls as 1/distance².
That is the same geometry as Newton's pull. Add up the flows from many orderly sources
arrow by arrow, as vectors, and the total flow points exactly along Newton's pull:

```
F_coherent = (ℓ / 4πG) · g_N
```

Here g_N is the ordinary Newtonian pull of all the visible matter. This is Gauss's geometry
and holds exactly; our script checks it to 2 parts in 10¹⁵. Opposite contributions cancel,
just as Newton's pulls do.

### 3.3 The companion pulls with its strength

Give the companion an energy density A²/8πG, and let its pull equal its amplitude A. This
has the same form as the usual energy density of a gravitational field, g²/8πG. Energy balance through a sphere
around a mass M then says that what streams out equals what is fed in:

```
4πr² · u · A²/8πG = ℓ M   ⟹   A = √(G M a) / r,   with   a = 2ℓ / u
```

This pull falls as 1/r, not 1/r². A star orbiting at distance r then has
v² = r·A = √(G M a), which is **the same at every radius: a flat rotation curve.** It also
gives the observed rule **v⁴ = G M a**, which links a galaxy's mass to its rotation speed.

The constant a is not a new constant of nature here. It is twice the power per kilogram
divided by the companion's speed. With the fitted values, ℓ = 6.5 × 10⁻⁶ watts per
kilogram. Over 13.8 billion years that is three thousandths of one percent of a kilogram's
rest energy.

### 3.4 Strong gravity holds the companion back

Near the Sun the pull is billions of times stronger than in a galaxy's outskirts. If the
companion pulled there, planetary orbits would show it; they don't. So the companion must
be held back, "attached", where ordinary gravity is strong, and released where it is weak:

```
released fraction = exp( −|g_N| / g_d )
```

Here g_d = 2.26 × 10⁻¹⁰ m/s², measured on the galaxies. At every planet the released
fraction is zero to the precision of a computer. The same factor also improves the galaxy
fits.

### 3.5 Random motion scrambles the companion, and scrambled contributions don't cancel

A moving source Doppler-shifts whatever it emits. Take a population whose parts move
randomly, with a spread of speeds σ. Each part shifts its companion differently, so the
contributions drift out of step. Out-of-step contributions add as plain totals, not as
arrows, so they no longer cancel. We call the plain total the "scalar sum":

```
S = G ∫ k ρ / d²,   heat weight  k = 3σ² / u²
```

The heat weight k comes from the same Doppler effect. A moving source's received power goes
up with the cube of its Doppler factor. Averaged over random directions, that gives an extra
3σ²/u².

We checked the scrambling numerically. With emitters all in step, the time-averaged strength
is the arrow sum: 160.50, against a predicted 160.49. With emitters moving randomly at
0.3 u, it is the plain total: 210.55, against 209.20.

Why this matters is geometry. Inside a big, spread-out cloud, Newton's pulls from all sides
largely cancel, but the plain total does not. In a cluster, that is exactly the extra pull
needed: large near the centre and smaller toward the edge.

### 3.6 Collisions switch the scrambling off

**This is the idea that solved colliding clusters.** Laboratory physics has a well-known
exception to Doppler scrambling. If an emitter changes direction many times before it has
moved one wavelength, its Doppler shifts average away. Its emission then looks like that of
an emitter at rest. The effect is called Dicke narrowing (R. H. Dicke, 1953). Its cousin,
the Mössbauer effect, is why atoms locked in a crystal emit perfectly sharp lines.

Gas particles collide, and in clusters they also spiral around magnetic fields, constantly.
Stars never collide; galaxies pass through each other. So two things at the same
temperature behave oppositely:

| | Moves freely? | Companion |
|---|---|---|
| Stars in an elliptical galaxy; galaxies in a cluster | yes | scrambled: **hot** |
| Gas and plasma, including the Sun's interior | no, constant collisions | in step: **cold** |

Our simulation shows it:
* An emitter that collides 100 times per companion wavelength stays in step about 70 times
  longer than a free one (102 time units against 1.4).
* A population of colliding emitters adds up almost like cold matter: 1.2 times the arrow
  sum, where free movers at the same temperature reach 11 times it.

So **only freely moving matter, stars and galaxies, carries the heat term.** Gas still counts
in full as ordinary mass.

### 3.7 The companion pulls along its net flow

Out-of-step waves don't add their strengths as arrows, but their energy *flows* still add
as arrows. So the net flow of the companion points along Newton's pull plus the heat-weighted
flow from the free-moving matter. The extra pull goes that way:

```
direction of the extra pull = (g_N + g_hot) / (|g_N| + |g_hot|)
g_hot = the Newtonian pull of the free-moving matter, weighted by its heat k
```

In a round, settled galaxy or cluster, both flows point to the centre, and nothing changes.
In a collision, the gas dominates Newton's pull but the galaxies dominate the heat-weighted
flow. The pull then bends toward the galaxies, and so does the lensing.

Where the two flows point in opposite directions, they partly cancel and the extra pull fades
smoothly to zero. An earlier version used the flow's direction alone, which flipped abruptly
there and would have planted thin sheets of *negative* lensing. The version above does not.

### 3.8 The whole law

```
g_N  = Newton's pull of all ordinary matter (stars, gas, everything)
S    = G ∫ k ρ_free / d²                   the scrambled total of the free-moving matter
g_hot = G ∫ k ρ_free (x′ − x) / |x′ − x|³   its heat-weighted flow
k    = 3σ²/u² for stars and galaxies,  0 for gas and plasma

extra pull  = exp(−|g_N|/g_d) · √( a (|g_N| + S) )   along   (g_N + g_hot)/(|g_N| + |g_hot|)
h = g_N + extra pull ;   ∇²Φ = −∇·h ;   the pull on matter and on light is  g = −∇Φ
```

The last line turns the rule into a proper field that conserves energy for anything moving
through it. Mathematically it has the same shape as Milgrom's QUMOND field equation, which
we use for its structure, not its physics (§10).

**Three constants, all universal, none set per object:**

| Constant | Value | Measured on |
|---|---|---|
| a | 6.56 × 10⁻¹¹ m/s² | 149 galaxies |
| g_d (release level) | 2.26 × 10⁻¹⁰ m/s² | 149 galaxies |
| u (companion speed) | 197 km/s | 12 clusters |

**For round, settled objects the law is simply**

```
g = g_N + exp(−|g_N|/g_d) · √( a (|g_N| + S) )
```

That is what we fit to galaxies and clusters.

### 3.9 Momentum: the companion carries it

A field that streams can carry momentum, as light does. A lamp does not recoil when its light
later pushes on something far away.

We tried the alternative, making ordinary matter alone balance every push and pull. It fails:
* the stars, only about 5% of a cluster's mass, would have to absorb the reaction to steering
  the gas;
* they would be pushed outward with a third to three quarters of the force of gravity;
* cluster galaxies would then orbit as if clusters held a fraction of their mass, and they
  don't.

So we take the streaming option: gas and galaxies feel the same pull, as observed.

The cost is a small self-push on lopsided systems. For the Bullet Cluster it is 12% of the
average pull. That is 13% of the momentum the companion's stream can carry, so the books
balance.

## 4. MOND, as the cold limit of our law, not an input

Switch off the heat (S = 0), as for a cold disk galaxy, and the law becomes

```
g = g_N · ν(g_N/a),   ν(y) = 1 + e^(−y/λ) / √y,   λ = g_d/a = 3.45
```

That is a MOND-type formula. **But it came out of the companion mechanism; it was not put
in:**

| MOND's ingredient | In MOND | Here |
|---|---|---|
| Its constant a₀ | a new constant of nature | 2 × power per kilogram ÷ companion speed |
| The switch between Newton and the deep regime | chosen by hand | set by the release factor. Not identical to any published choice; the nearest (McGaugh, Lelli & Schombert's exponential) is 0.030 dex (7%) away |
| Why galaxies obey it | assumed | disks are cold |
| Why clusters don't | needs extra matter | cluster galaxies are hot |
| Why collisions don't | needs extra matter | galaxies carry the heat, gas doesn't |

This is how Newton's law comes out of Einstein's in weak gravity: a special case of a larger
law. The formula check still flags the cold limit as MOND-like, as it should. Every result we
claim as ours rests on where the law is not MOND.

## 5. Tests against real data

### 5.1 Rotation curves: 149 galaxies, 3,150 measured speeds

Average miss per galaxy (km/s; lower is better). The galaxies were split in advance into
training, validation and test sets.

| | All | Training | Validation | Test | Adjustable numbers |
|---|---:|---:|---:|---:|---|
| **Ours** | **15.85** | 15.95 | **19.21** | **12.40** | 3 in total |
| MOND | 16.13 | 15.85 | 19.80 | 13.52 | 1 |
| Newton, visible matter only | 45.58 | 45.94 | 48.64 | 41.68 | 0 |
| Dark-matter halo fitted to each galaxy | 7.52 | | | | 298 |

* Ours beats MOND overall and on the held-back galaxies; MOND is slightly better on the
  training set.
* **Galaxies with big bulges** were the risk. Their bulges are hot, and the new heat rule
  makes hot stars count much more. Those 25 galaxies came out better: 29.2 km/s against
  MOND's 30.4.
* Dark matter fits more tightly because it has a hundred times more adjustable numbers.

Script: `code/run_v3.py`.

### 5.2 Galaxy clusters: 12 clusters with measured gas, temperatures and stars

The quantity is the typical miss in the mass needed to hold the gas, at six radii from the
core to the edge.

| | Typical miss | Adjustable numbers |
|---|---:|---|
| **Ours** | **25%** (rms 0.227 in ln M) | same 3 constants |
| Ours, clusters held out of the fit (924 splits) | 28% (0.244) | |
| Ours, stars' speeds computed by our law from visible matter alone | 39% (0.329) | |
| MOND | ×2.9 (1.062) | 1 |
| MOND with its constant refitted on clusters | ×1.4 (0.338); needs a constant 9.7× its galaxy value | 1 |
| Newton | ×9.1 (2.211) | 0 |
| Dark matter (NFW fit per cluster) | 11% (0.101) | 24 |

**Why the stars can do it.** The stars are concentrated in the middle, where clusters need
the most extra pull:

| Radius (fraction of the cluster's size, R500) | 0.02 | 0.05 | 0.1 | 1 |
|---|---:|---:|---:|---:|
| Mass in stars ÷ mass in gas | 2 to 11 | about 1 | 0.3 to 0.5 | 0.02 to 0.07 |

Their random speeds, roughly 300 to 1,200 km/s, make each kilogram count dozens of times
over.

The companion speed comes out at u = 197 km/s. Across the 924 ways of splitting the clusters
in half, it stays between 178 and 218 km/s in 90% of them.

Script: `code/run_v3.py`, with stellar mass profiles from the X-COP release (Ghizzardi et al.
2021).

### 5.3 Colliding clusters

**The Bullet Cluster, on its published numbers.** Clowe et al. (2006) measured the gas, the
stars and the lensing strength κ inside 100-kiloparsec circles at four places:

| Place | Gas (10¹² suns) | Stars (10¹² suns) | Measured lensing κ |
|---|---:|---:|---:|
| Main cluster's galaxies | 5.5 | 0.54 | 0.36 ± 0.06 |
| Main cluster's gas | 6.6 | 0.23 | 0.05 ± 0.06 (extra) |
| Subcluster's galaxies | 2.7 | 0.58 | 0.20 ± 0.05 |
| Subcluster's gas | 5.8 | 0.12 | 0.02 ± 0.06 (extra) |

**How we modelled it.**
* We built the gas and the stars to match those eight masses.
* We solved our field equation in 3D and projected the result exactly as a lensing map is
  made.
* We applied the same analysis the observers used. Nothing about the lensing was fitted.
* The stars' speeds were computed by our own law from the visible matter.

| | Main galaxies | Subcluster galaxies | Main gas | Subcluster gas | Where the peaks are |
|---|---:|---:|---:|---:|---|
| Measured | 0.36 ± 0.06 | 0.20 ± 0.05 | 0.05 ± 0.06 | 0.02 ± 0.06 | on the galaxies |
| **Our law** | 0.51 | 0.07 | 0.04 | 0.04 | **8 kpc and 34 kpc from the two brightest galaxies** |
| Our law, but with the old direction rule (pull along ordinary gravity only) | 0.35 | 0.13 | **0.21** | 0.07 | **on the gas** |
| Our law, but with the gas counted as hot (round-1 rule) | 0.46 | 0.21 | **0.34** | 0.08 | **on the gas** |
| Our law with no heat at all (MOND-like) | 0.14 | 0.07 | 0.06 | 0.04 | **on the gas** |

**What this shows:**
* **The pattern is reproduced.** Both lensing peaks sit on the galaxies. The observed
  subcluster peak sits 43 kpc from its brightest galaxy, on the gas side; ours sits 34 kpc
  away, on the same side. The gas regions carry almost no extra lensing, as measured.
  **This is the pattern usually called impossible without dark matter.**
* **Both new ideas are needed.** Drop either one and the lensing lands on the gas.
* **Our law predicts the galaxies' speeds, and they check out.** For the main cluster it
  gives 620 km/s near the centre, rising to 1,140 km/s at 400 kpc. The measured value is
  1,249 ± 100 km/s from 71 galaxies.
* **The strengths are not yet right.**
  * The main cluster comes out 40% too strong. Lighter stars fix it: the published star masses
    are upper limits with a mass-to-light ratio uncertain between 0.5 and 3, and at 1–1.5 we
    get 0.31–0.41.
  * The subcluster comes out too weak, 0.07 against 0.20 ± 0.05. It needs its stars to have
    moved faster than our law gives, about 500 km/s against 270–440. Its 7 measured galaxies
    give 212 ± 60 km/s.
  * This is the main open item (§8).

**72 colliding clusters.** Harvey et al. (2015) stacked 72 pieces of colliding clusters. The
lensing sits within 5.8 ± 8.2 kpc of the galaxies.
* We ran 20 pieces of simulated collision under our law, at different stages and mass ratios.
* The lensing stayed with the galaxies every time, 4–22 kpc from them. Our map resolution is
  18 kpc.
* It never drifted to the gas. That drift was the failure of our previous attempt, which is
  now gone.

Scripts: `code/bullet_v3.py`, `code/collisions_v3.py`.

### 5.4 Lensing by galaxies

**Ellipticals versus spirals (KiDS survey).** Brouwer et al. (2021) found that elliptical
galaxies bend light more than spirals with the same stars. The mean difference is 0.17 dex
(splitting by shape) and 0.27 dex (splitting by colour), at more than 5.7σ. A law that treats
all matter alike cannot produce that.
* In ours, an elliptical's stars move randomly at 150–200 km/s and never collide, so they are
  hot. A spiral's stars circle in step.
* **Our prediction: 0.17 to 0.27 dex, roughly constant beyond 100 kpc.**
* No hidden gas is needed. If ellipticals do have large hot gas haloes, the gas adds only its
  mass: 0.02–0.07 dex more for haloes of 0.3–1 times the stars. So we predict those haloes
  are modest.

Script: `code/kids_v3.py`.

**Six strong lenses (SLACS).** Each lens has a measured ring size and measured star speeds.
With standard distances and published star masses:
* **The stars need to be 1.05–1.35 times heavier than the "Salpeter" standard assumption.**
  Spectra of giant ellipticals already suggest star populations this heavy. Measuring these
  six lenses' star populations directly is the clean test.
* The mass needed to bend the light and the mass needed to move the stars agree to
  −0.017 ± 0.024 dex. **Light and matter feel the same pull.**

Script: `code/lenses_t35.py --constants ../run-v3/results.json`.

### 5.5 The Solar System and the Sun

* **Planets.** At every planet, the release factor switches the extra pull off entirely.
* **The Sun.** Its interior is colliding plasma, so it has no heat term. It feeds the
  companion at ℓ = 6.5 × 10⁻⁶ W/kg: 1.3 × 10²⁵ W, or 3.4% of its light. That means an
  extra mass loss of 2.3 × 10⁻¹⁵ of its mass per year, below what planetary tracking can see.

### 5.6 The formula check

* **Not MOND.** For the same Newtonian pull, random speeds of 0–1,000 km/s change the
  prediction by a factor of 8.6 (0.93 dex). Across the 3,150 galaxy measurements, 0.022 dex
  of the prediction depends on more than the local Newtonian pull.
* **Not Newton.** The boost varies from object to object.
* **Not dark matter.** Nothing invisible is added, and no number is set per object.
* **The cold limit is MOND-like**, and the check reports it (§4).

## 6. How this compares

| | Ours | MOND | Dark matter |
|---|---|---|---|
| Rotation speeds, 149 galaxies | **15.9 km/s** | 16.1 km/s | 7.5 km/s (298 numbers) |
| Cluster masses, 12 clusters | **25%** | ×2.9 | 11% (24 numbers) |
| Bullet Cluster: lensing on galaxies | **yes** | no | yes |
| Bullet Cluster: lensing strengths | main within the star-mass range; subcluster 2.5σ low | no | yes (fitted) |
| Collisions: lensing stays with galaxies | **yes** | no | yes |
| Ellipticals lens more than spirals | **yes, 0.17–0.27 dex, from their stars** | no | yes, via tuned haloes |
| Strong lenses: light and stars agree | **yes** | | yes |
| Solar System | **silent** | small effect | silent |
| Explains *why* | disks are cold; cluster galaxies are hot; gas collides | no | no |
| Adjustable numbers | **3** | 1 | ~320 |

Dark matter fits individual objects more tightly because it is tuned object by object. Ours
fits everything with three shared numbers, and says why each kind of system behaves as it
does.

## 7. Predictions anyone can check

1. **Heavy stars in the six SLACS lenses:** 1.05–1.35 times Salpeter, measurable from their
   spectra. Dark-matter models expect about 1.0.
2. **The elliptical/spiral lensing gap is flat beyond 100 kpc.** A hot-gas-halo explanation
   makes it grow with radius.
3. **In relaxed clusters, galaxy orbits and gas agree.** Both feel the same pull.
4. **In every collision, lensing stays with the galaxies**, at every stage. Where collisionless
   streams overlap at high speed, the lensing should be boosted.
5. **The extra pull tracks random motion.** At equal visible mass, systems whose stars move
   randomly and freely pull harder than those whose stars circle in step.
6. **The Sun's extra mass loss** is 2.3 × 10⁻¹⁵ per year.

## 8. What is still open, and why we are optimistic

1. **The subcluster's lensing strength (2.5σ low).** Its stars need to have moved at about
   500 km/s; our law gives 270–440 km/s from what remains of it. Stars keep the speeds they
   had before a collision strips the gas. Better estimates would settle it either way:
   * the subcluster's mass before the collision;
   * its star masses (the published ones are upper limits);
   * a speed measured from more than 7 galaxies.

   There is also a physical effect we have not yet modelled. Streams of galaxies passing
   through each other at 3,000 km/s look "hot" to each other, which would add pull exactly
   there.
2. **A full field theory for the companion.** It needs its momentum and a relativistic form,
   so that lensing is derived rather than assumed equal to the pull on matter. The data
   support that equality: the six lenses agree to 0.017 dex.
3. **The microscopic numbers.** We want to derive u and g_d, and pin down the companion's
   wavelength. That wavelength must lie between the gas's collision length and the stars'
   orbits. For cluster gas this requires particles to be scattered, mainly by magnetic
   fields, at least every 0.1–1 million years. That is thought likely but is not measured.
4. **Wide binary stars.** Near the Sun, the Galaxy's pull is close to the release level. Our
   law therefore predicts a weaker signal than MOND does. The current data are disputed
   between groups; we will compute our exact prediction next.
5. **Cosmology.** The cosmic microwave background and the growth of large-scale structure are
   outside this law's scope as tested so far. Under the project's no-expansion rule they need
   their own treatment.

Each of these is a concrete calculation or measurement, not a wall.

## 9. How we got here

* **Revisions 1–9.** A "streaming gravity" law plus a light-versus-matter "slip". The law
  turned out to be MOND's, so under the new rules it was retired (revision 10). The slip later
  proved an artefact of measuring against MOND.
* **Revision 10 (round 1).** Hot-companion gravity, with the gas as the hot component. It fit
  galaxies and clusters (u = 874 km/s), but put the Bullet's lensing on the gas.
* **Revision 11 (round 2).**
  * MOND derived as the cold limit.
  * A field equation, and an action to conserve momentum.
  * Hot gas haloes to explain the ellipticals.
  * The Bullet partly solved by "companion memory", but that failed the 72-collision stack.
* **Revision 12 (round 3, this page).** The collision rule and the flow direction. Stars now
  carry the heat (u = 197 km/s), and the Bullet's pattern is reproduced.

**Superseded along the way, kept on the record:**
* round 2's hot-gas-halo explanation of the ellipticals (now it is their stars);
* round 2's reaction force on hot matter (now the companion carries momentum);
* round 1's companion speed (874 km/s, now 197);
* round 2's solar mass-loss figure (1.4 × 10⁻¹⁴ per year, now 2.3 × 10⁻¹⁵);
* earlier retractions:
  * a cluster claim (revision 5, retracted in revision 7);
  * a misattributed group-lensing figure (revision 8).

## 10. What is borrowed and what is ours

**Borrowed, and credited:**
* Newton's and Einstein's gravity in strong fields.
* Gauss's flux geometry.
* Dicke narrowing (Dicke 1953) and the Mössbauer effect, as known physics that we apply to
  the companion.
* The mathematical form of Milgrom's QUMOND field equation.
* The data:
  * SPARC (Lelli, McGaugh & Schombert);
  * X-COP (Eckert, Ettori, Ghirardini, Ghizzardi and collaborators);
  * Clowe et al. 2006;
  * Barrena et al. 2002;
  * Harvey et al. 2015;
  * KiDS-1000 (Brouwer et al. 2021);
  * SLACS (Bolton, Auger, Treu and collaborators).

**Ours, as far as we have found:**
* the companion mechanism;
* heat as loss of step;
* collisions switching it off;
* the pull along the companion's net flow;
* MOND's constant and switch derived from it.

A full literature search is still owed before any claim of priority.

## 11. Reproduce it

All scripts are in `research_work/results/hot-companion/code/`. Each writes to a fresh
output folder.

```
python dicke_toy.py      --output-dir ../run-dicke          # collisions switch scrambling off
python run_v3.py         --output-dir ../run-v3             # galaxies, clusters, constants, guard
python bullet_v3.py      --output-dir ../run-bullet-v3      # the Bullet Cluster
python collisions_v3.py  --output-dir ../run-collisions-v3  # 20 simulated collision pieces
python kids_v3.py        --output-dir ../run-kids-v3        # ellipticals vs spirals
python lenses_t35.py     --output-dir ../run-lenses-v3 --constants ../run-v3/results.json
python derive_mond.py    --output-dir ../run-derive         # MOND as the cold limit
python field_equation.py --output-dir ../run-field          # the 3D field equation
```
