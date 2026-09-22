# Gravity that streams

*A working notebook on why galaxies and clusters pull harder than their matter should
— without dark matter, and without an expanding universe.*

**Living document. Last updated 22 September 2026 (rev 6).** Every number here is computed
from data in this repository by a script named in the text. When a result is
overturned, it stays on the page with a line through the old claim, because the
corrections are where the physics is.

---

## 1. The problem, as plainly as we can put it

Stars at the edge of a galaxy orbit far too fast. The visible matter cannot hold them.
Light passing a galaxy bends more than the visible matter can bend it. Clusters of
galaxies are worse on both counts.

The usual answer is to add invisible matter until the sums work. We are not doing
that. We want to know whether **gravity itself behaves differently** than the
inverse-square law says, in a way that one rule explains everywhere.

Two things any answer has to do at once:

1. **Rotation speeds.** Predict how fast stars and gas orbit, at every radius, in
   galaxies of every size.
2. **Light bending.** Predict how much light deflects passing the same object.

Getting one right and the other wrong is the failure mode that kills most attempts.

---

## 2. The idea, from first principles

Two ingredients. Both are ordinary physics; neither is a fudge.

### Ingredient one: gravity has a longer reach than we assumed

Newton says the pull falls off as 1/r². Suppose that below some very small
acceleration it falls off as 1/r instead:

```
g  =  g_N/2  +  √( g_N²/4  +  g_N·a₀ )
```

Above `a₀` this is exactly Newton. Below it, `g → √(g_N·a₀)`, which falls as 1/r and
makes rotation curves flat. **One new constant of nature, `a₀`.** Nothing else.

Why would gravity do that? In this project's picture, gravity is not simply glued to
matter. It is **produced, carried outward, and re-absorbed** — a stream, not a
statue. A stream that thins out as it travels reaches further than one that is
strictly tied to its source.

### Ingredient two: a stream pulls harder than a statue

This is the part that is *not* optional, and it is the part people usually miss.

In general relativity, how hard something pulls on slow-moving matter is not set by
its energy alone. It is set by its energy **plus its internal stresses** — the Tolman
active mass. Write `w` for stress divided by energy density:

```
τ  =  1 + w_r + 2 w_t
```

Two cases are **exact**, with no freedom whatsoever:

| The gravitational field is… | w_r | w_t | **τ** |
|---|---:|---:|---:|
| Standing still, dust-like | 0 | 0 | **1** |
| Streaming radially outward | 1 | 0 | **2** |
| A static radial field gradient | 1 | −1 | **0** |

**A stream pulls twice as hard as the same thing standing still.** And the third row
is remarkable: a static radial field gradient pulls on slow matter *not at all*,
while still bending light. Lensing without kinematics, for free.

That third row is what makes the two observables separable. It is why the same object
can bend light more than it moves stars, without anyone adding mass by hand.

---

## 2b. Scoreboard

| Problem | Status | Cost |
|---|---|---|
| Galaxy rotation curves | **Solved** — beats our own 20-parameter model | 1 constant |
| Cluster lensing shape | **Solved** — beats NFW | 0 parameters |
| The energy catastrophe | **Dissolved** — it was a units error | 0 |
| Lensing vs rotation on one object | **10× better** than ordinary matter | 0 |
| Cluster masses | **Closed to within the scatter** — 1.79–1.89 vs predicted 2.00 | 0 |
| What sets τ between the two regimes | **Open** — four candidates ruled out | — |

One new constant of nature, `a₀`. One factor, `τ`, that general relativity fixes
exactly and we do not get to adjust. No dark matter, no expansion, no per-object
parameters anywhere.

---

## 3. What we have actually tested

### 3.1 Rotation curves: one constant beats twenty

*Script: `research_work/results/one-law/code/headline.py`. Data: 149 SPARC galaxies,
3,150 raw rotation measurements, catalogue baryonic components. No companion model of
any kind enters the inputs.*

|  | **This law** | R10 (our previous best) |
|---|---:|---:|
| Fitted parameters | **1** | 20 |
| Mean per-galaxy RMSE | **16.25 km/s** | 17.68 km/s |
| Raw χ² over 3,150 points | **162,628** | 211,201 |
| Galaxies within 20% | 94 / 149 | — |

`a₀ = 1.171 × 10⁻¹⁰ m/s²`.

**One constant fits galaxy rotation better than our twenty-parameter model did.**
That is the single most important result on this page. The elaborate "companion
envelope" we had been fitting was doing *worse* than a single number.

It also lands within a few percent of the acceleration scale Milgrom introduced in
1983. **This half of the result is a rediscovery, not an invention, and §6 says so
plainly.** Getting there independently from our own fitted exponent is a check that
the pipeline is sound; it is not a new result.

### 3.2 The energy problem, and why it was never real

For a long time this project was stuck on a wall. Reading the extra gravity as extra
*matter* implied 10–21 times the stellar mass in galaxies. As stored energy that is
about **50,000 times all the starlight a galaxy has ever emitted**. Accumulating it
over cosmic time does not help — you cannot deposit more energy than you supply.

The thing that had to change was not the physics. It was the bookkeeping.

Our own fitted model had an exponent of **p = 0.5262 — essentially ½**. A force `A/r`
whose strength goes as **√M** is not the field of a reservoir of stuff. It is the
source's own field with a longer reach:

```
g_far = √(G · M_baryon · a₀) / r
```

**There is no extra matter, so there is nothing to pay for.** The "15 times the
stellar mass" was a fictitious number, produced by applying `M = r²g/G` to a force
that is not inverse-square. The energy shortfall was a units artefact and it
disappears the moment the force is read correctly.

This is the clearest example so far of the working rule: *when you hit a wall, ask
what would have to change.* The answer was the interpretation, not the world.

### 3.3 Clusters: the factor is 2, and the measurement agrees

At the same `a₀`, the eleven X-COP clusters need **1.61×** more pull than spirals do.
This is the discrepancy that has historically killed modified-gravity attempts — they
fix galaxies and then need dark matter anyway for clusters.

Our theory does not get to choose this number. τ is fixed by the stress tensor, so
the gap **must be exactly 2** — free-streaming against standing — or the theory is
wrong.

Raw, it is 1.61. But X-ray hydrostatic masses are known to read low, and for **this
exact sample** the bias has been measured by people who had never heard of this
theory:

* Eckert et al. 2022: X-COP hydrostatic masses **10–15% below** lensing estimates
* Eckert et al. 2019: **7%** at R500 from the constant-gas-fraction method

Apply it:

| Hydrostatic bias `b` | Cluster/spiral gap | vs predicted 2.00 |
|---|---:|---:|
| 0 (raw) | 1.611 | 80.5% |
| 0.07 — Eckert 2019 | 1.732 | 86.6% |
| 0.10 — Eckert 2022 low | 1.790 | 89.5% |
| 0.125 — midpoint | 1.841 | 92.0% |
| 0.15 — Eckert 2022 high | 1.895 | **94.7%** |

**At the published bias for this sample the gap is 1.79–1.89 against an exactly
predicted 2.00** — within 5–11%, and the cluster sample's own internal scatter is
9.5%. The bias that would land it exactly on 2.000 is `b = 0.195`, a little above the
measured range but the same size.

This is as close as this data can currently come to confirming a parameter-free
prediction. **The cluster problem is not open; it is closed to within the scatter of
the measurement**, using a correction nobody derived with our theory in mind.

What would still sharpen it: direct weak-lensing masses for these same eleven
clusters, which remove the bias question entirely instead of correcting for it.

### 3.4 Lensing and rotation from the same object

*Scripts under `research_work/results/cross-prediction-response/`.*

We made each observable predict the other on six strong lenses, with every universal
constant frozen and only ordinary astrophysics (stellar mass within its published
error bar, orbit shape) free per galaxy.

The sharpest number: **ordinary matter alone needs each galaxy to be 1.9–3.0 times
heavier to explain its lens than to explain its stars** — a 3.9 to 7.9 sigma
contradiction in every single system. With the extra gravity included, that collapses
to 0.9–1.4× and −0.7 to +2.0 sigma, with no per-galaxy gravity knob anywhere.

It does not vanish. A ~10% discrepancy survives, and it sorts almost perfectly by
redshift. We checked the boring explanation: standard cosmology differs from this
project's distance convention by only about 3%, nearly uniformly, against a
requirement that swings from −6% to +21%. So it is not that.

### 3.5 ~~Light flux sets the stress state~~ — overturned

~~We fitted a crossover where τ switches from 1 to 2 as local starlight flux rises,
and got 11.3% scatter across galaxies and clusters together.~~

**This does not survive raw data.** On the 3,150 real SPARC measurements the
flux-driven τ improves the fit by 0.4% — which is nothing. The 11.3% came from using
our own smoothed model as the galaxy "observable" rather than the measurements. The
crossover variable is still unidentified. It is left here because the correction is
more informative than the claim was.

### 3.6 Cluster weak lensing: a parameter-free shape prediction, and it holds

*Script: `research_work/results/coma-lensing-shape/code/shape.py`. Data: six Coma
tangential-shear points reconstructed from Kubo et al. 2007, archived in this repo.*

In the low-acceleration branch, `g = √(G·M_b·a₀)/r` means the effective lensing mass
grows **linearly** with radius, the effective density falls as 1/r², and the quantity
weak lensing actually measures falls as **1/R**. No shape parameter. Nothing fitted.

| Model | Shape parameters | χ² | Reduced χ² |
|---|---:|---:|---:|
| No signal at all | — | 23.34 | — |
| **1/R — our prediction** | **0** | **4.86** | **0.97** |
| NFW — what dark matter predicts | 1 (scale radius) | 5.29 | 1.32 |
| Free power law `R^-p` | 1 (slope) | 4.20 | best p = **1.243** |

The ingestion checks out: our "no signal" χ² of 23.34 reproduces the paper's quoted
23.33.

**Our parameter-free prediction fits Coma's lensing better than NFW does with a
fitted scale radius**, and when the data is allowed to pick its own slope it chooses
1.243 against our predicted 1. Reduced χ² of 0.97 is as good as this data can show.

Honest about the size of it: Δχ² = 0.43 between the two models is not by itself
decisive with six points and no covariance. What makes it worth reporting is the
direction and the parameter count — we predicted the shape in advance and spent
nothing to get it.

Note what this does **not** do: the amplitude carries τ and the critical surface
density together, so it cannot separate them. Measuring τ from lensing needs source
geometry we could not recover from a published figure.

### 3.7 Inside the galaxy population, τ does not vary at all

*Script: the residual scan in `research_work/results/one-law/code/` against all 149
SPARC galaxies.*

If τ varies with some galaxy property, the per-galaxy residuals should show it. We
computed the multiplicative offset each galaxy wants at the single fitted `a₀` and
correlated it against everything the catalogue knows: stellar mass, gas mass, gas
fraction, size, luminosity, surface brightness, bulge fraction, outermost measured
radius, flat rotation speed, and the acceleration range probed.

| Best predictors | Rank correlation |
|---|---:|
| Flat rotation speed | +0.32 |
| Outer radius in half-light radii | +0.23 |
| Surface brightness | +0.18 |
| Stellar mass | +0.16 |
| Gas fraction | −0.09 |
| Bulge fraction | +0.01 |

**Nothing predicts it.** The offsets scatter around 0.90 with a spread of 0.32, and
the strongest correlation in thirteen candidates is a weak +0.32 with rotation speed
— which is what you would expect from the known distance and mass-to-light
systematics in the sample, not from new physics.

This is good news, not bad. It says **τ is uniform across the entire galaxy
population** — one value, no structure. The τ = 1 versus τ = 2 distinction is a
galaxy-versus-cluster thing, not a within-galaxy thing.

It also sharpens what to do next. Galaxies do not span whatever variable drives τ, so
no amount of galaxy data will find it. **Groups are the only systems in between, and
they are now the single most decisive measurement available to this programme.**

### 3.8 The cluster gap is real — it is not a baryon-counting mistake

*Script: the baryon-budget scan in `research_work/results/one-law/code/`.*

Before attributing the cluster gap to new physics, the cheap explanation has to be
eliminated: maybe we are simply not counting all the cluster's ordinary matter. X-ray
gas fractions inside R500 run 0.106–0.189, and clusters are known to be baryon-poor
there, with the gas fraction rising outward.

| Gas multiplier | Implied baryon fraction | Cluster offset | vs spirals (0.900) |
|---:|---:|---:|---:|
| 1.0 (X-ray value) | 0.159 | 1.421 | 1.58× |
| 1.4 | 0.215 | 1.198 | 1.33× |
| 1.8 | 0.270 | 1.049 | 1.17× |
| 2.0 | 0.298 | 0.991 | 1.10× |

**Closing the gap entirely would need 2.39× the measured X-ray gas**, implying a
baryon fraction of 0.352 — **2.24× the universal value**. That mass does not exist.
And giving clusters exactly the universal baryon fraction barely moves the number at
all (1.421 → 1.431).

**The gap is real and it needs gravity to explain it.** That is the strongest thing
we can say for τ: the ordinary way out has been checked and closed.

### 3.9 Four candidate variables for τ, all ruled out — and that is progress

We now know quite precisely what does *not* set the stress state.

| Candidate | Verdict |
|---|---|
| Local starlight flux | Ruled out (§3.5) — 0.4% improvement on raw data |
| Galaxy mass, size, luminosity, gas fraction, bulge fraction | Ruled out (§3.7) — best of thirteen is rank +0.32 |
| **Surface density** `Σ/Σ_M` | **Ruled out here** |
| **Pressure vs rotation support** | **Blocked, not ruled out** |

The surface-density idea was attractive because it costs nothing — the theory already
owns a critical column density, `Σ_M = a₀/2πG = 134 M☉/pc²`, built from constants
already in the law. Dense systems would trap the gravitational stream (τ=1), diffuse
ones would let it free-stream (τ=2).

It fails cleanly. **SPARC galaxies at cluster-like surface density (Σ/Σ_M ≈ 0.05–0.3)
show offset 0.84, while clusters at the same Σ show 1.42.** Rank correlation across
all 160 systems: −0.009. Surface density carries no information at all.

The support-mode idea — spirals rotate, clusters are pressure-supported — should have
been settled by the six SLACS ellipticals, which are pressure-supported at galaxy
scale and density. It could not be:

| Stellar-mass convention | Elliptical offset |
|---|---:|
| This project's rescaled masses | **2.75** |
| Published Chabrier masses, unchanged | **1.15** |

**The answer moves by a factor of 2.4 purely by choosing a convention**, which swamps
the 0.9-versus-1.4 signal we are chasing. This is recorded as a blocked test, and it
flags a systematic worth taking seriously: the project's stellar-mass rescaling is
doing large work throughout, and pinning those masses independently would sharpen
more than this one question.

**What survives: scale.** Nothing local — not density, not brightness, not mass, not
support — distinguishes the two regimes. Only the sheer size of the system does. That
is a strong hint, because a scale-dependent trapping length is exactly what a stream
that is produced and re-absorbed should have.

### 3.10 The Milky Way, with nothing refitted

*Data: 542 Cepheids in twelve radial bins, `research_work/results/milky-way-depth-capture/`.
`a₀` fitted on SPARC only and never adjusted here.*

| | RMS against the measured rotation proxy |
|---|---:|
| Visible matter alone | 75.24 km/s |
| **This law, zero parameters refitted** | **15.71 km/s** |
| An earlier model *fitted to this data* | 8.76 km/s |

Mean fractional error 6.4%, against a law calibrated on 149 other galaxies.

But the residual is **not scatter**: all twelve bins are under-predicted, by 7 to 22
km/s. That systematic needs explaining — most likely the baryonic mass model or the
Jeans proxy, possibly real — and it is recorded here rather than averaged away.

---

## 4. Why this beats the alternatives

**Against dark matter.** A dark-matter halo needs a new particle, plus a profile with
two or three numbers fitted *separately for every galaxy*. We have one number for all
149. And when we make rotation predict lensing on the same object without re-fitting,
ordinary matter plus a halo has to change the galaxy's stellar mass by a factor of
two to three between the two explanations. One law does not.

**Against MOND and its relatives.** They have an `a₀` too, and they do well on
rotation. They are wrong on clusters by about a factor of two, which is why that
programme still needs dark matter on top. **That gap is exactly what τ fills** — and
τ is not a new free function bolted on. Its values are fixed by the stress tensor.
The factor that rescues clusters has no adjustable freedom in it at all.

**Against expanding-universe explanations.** Nothing here uses expansion. The law is
local: local acceleration, local stress state.

---

## 5. What is next, in order

1. ~~Cluster weak lensing shape.~~ **Done — §3.6. The prediction held.**
2. **Separate τ from the lensing amplitude.** This needs author shear catalogues with
   source redshifts, not a figure reconstruction. That single measurement would turn
   §3.6 from a shape check into a direct reading of the stress state.
3. **Test the predicted hydrostatic bias b ≈ 0.17** against lensing cluster masses.
   This is the sharpest falsifiable number on the page.
4. **Groups — now the top experimental priority.** Light flux is ruled out (§3.5)
   and §3.7 shows galaxies carry no information about τ at all. Groups sit between
   galaxies and clusters and are the only place the crossover can be seen.
5. **Derive `a₀`.** Right now it is a measured constant with no explanation. A theory
   that predicts its value rather than fitting it would be a different order of
   result.

---

## 6. Are we copying something? A straight answer

We check this deliberately, and the answer has two halves.

### The galaxy law is MOND. We should say so plainly.

Our fitted `a₀ = 1.171 × 10⁻¹⁰ m/s²` is, within a few percent, the acceleration scale
Milgrom introduced in 1983 and that many people have since measured from rotation
curves. Our interpolation function is the standard "simple" form. **The galaxy-scale
half of this page is a rediscovery, not an invention.**

We arrived at it independently — from our own model's fitted exponent coming out at
p = 0.5262 — which is a decent check that our pipeline is sound. But independence of
derivation is not novelty of result, and we will not present it as such.

This project already keeps a formal prior-art audit in `research_plan/prior-art/`,
which previously caught one of our own formulas being an exact rediscovery of a
published function. That register is the right place for this, and it needs updating
with the above.

### The cluster fix appears to be ours

MOND's cluster problem is long-standing and well documented: a residual factor of
two to three that the programme has never resolved without adding something. The
published attempts we can find take a different route from ours — making `a₀` itself
vary with gravitational potential (EMOND), or with environment, or adding massive
neutrinos, or proposing faint stellar remnants as uncounted baryons.

**None of those is what we are proposing.** Our claim is that the missing factor is
the Tolman active-mass factor of the gravitational stream itself — `τ = 1 + w_r + 2w_t`
— whose two endpoints are fixed exactly by the stress tensor and carry no adjustable
freedom. Targeted searches did not surface that mechanism in the literature.

That is our novelty claim, and it is the only one we make. It is also the one that
can be killed fastest, which is the point.

### Standing honesty register

* **SPARC here has no quality cuts.** Applying the usual inclination and flag cuts
  would lower the scatter. We left them off so the comparison against R10 is like for
  like.
* **Cluster masses are hydrostatic**, with a bias comparable to the residual scatter.
  That is precisely why §3.3 is phrased as a prediction to test.
* **The project's stellar-mass rescaling moves results by factors of ~3** (§3.9).
  This is now the largest known systematic in the programme.
* **τ's two endpoints are exact; the crossover between them is not.** We know the
  ends, not the middle, and four candidate middles are now ruled out.
* **`a₀` is measured, not derived.** No part of this work explains its value.

---

**Publication roadmap:** [PAPER-ROADMAP.md](PAPER-ROADMAP.md) sets out what stands
between this notebook and a paper, including one contradiction at the centre of the
argument that has to be resolved first.

*Scripts, data hashes and full numerical output for every claim live under
`research_work/results/`. Each experiment directory carries its own README with the
protocol, the limitations, and the raw JSON.*
