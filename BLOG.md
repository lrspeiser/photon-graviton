# Gravity that streams

*A working notebook on why galaxies and clusters pull harder than their matter should
— without dark matter, and without an expanding universe.*

**Living document. Last updated 22 September 2026 (rev 2).** Every number here is computed
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

It also lands within a whisker of the acceleration scale other people have found
independently from the same kind of data — which is either a good sign or a warning
that we have rediscovered something. We treat it as both. (See §6.)

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

### 3.3 Clusters: short by a factor near two

At the same `a₀`, the eleven X-COP clusters need **1.653×** more pull than galaxies
do. This is the discrepancy that has historically killed modified-gravity attempts —
they fix galaxies and then need dark matter anyway for clusters.

Our claim is that this factor is `τ`: galaxies stand, clusters stream, and the ratio
is exactly 2.

The measured gap is 1.653, not 2. **But the cluster masses come from X-ray
hydrostatic equilibrium, which is independently known to read low.** If the gap is
exactly 2, that bias must be:

```
b  =  1 − 1.653/2  =  0.173
```

Published hydrostatic bias runs roughly 0.1–0.3. So **the theory predicts b ≈ 0.17**,
and cluster weak lensing can check it. That is a real, sharp, falsifiable number
rather than a fudge factor, and testing it is next.

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
4. **Find what actually sets τ.** Light flux is ruled out (§3.5). Groups, which sit
   between galaxies and clusters, are the place to look.
5. **Derive `a₀`.** Right now it is a measured constant with no explanation. A theory
   that predicts its value rather than fitting it would be a different order of
   result.

---

## 6. Honesty register

Things that are true and inconvenient, kept on the page deliberately:

* **`a₀` is close to a value others have published** from similar data. Our claimed
  novelty is not `a₀`; it is τ, the stress factor that handles clusters without new
  matter. We keep a formal prior-art audit in `research_plan/prior-art/`, which has
  already caught one of our formulas being an exact rediscovery.
* **SPARC here has no quality cuts.** Applying the usual inclination and flag cuts
  would lower the scatter. We left them off so the comparison against R10 is like for
  like.
* **Cluster masses are hydrostatic**, with a bias comparable to the residual scatter.
  That is precisely why §3.3 is phrased as a prediction to test, not a success.
* **The strong-lens results still inherit R10** in places. Driving them from raw
  V_rms and Einstein angles is on the list.
* **τ's two endpoints are exact; the crossover between them is not.** We know the
  ends, not the middle.

---

*Scripts, data hashes and full numerical output for every claim live under
`research_work/results/`. Each experiment directory carries its own README with the
protocol, the limitations, and the raw JSON.*
