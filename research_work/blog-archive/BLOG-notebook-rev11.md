# Gravity that streams

*A working notebook on why galaxies and clusters pull harder than their matter should
— without dark matter, and without an expanding universe.*

**Living document. Last updated 23 September 2026 (rev 11).** Every number here is computed
from data in this repository by a script named in the text. When a result is
overturned, it stays on the page with a line through the old claim, because the
corrections are where the physics is.

> **Rev 11: MOND is now an output, not an input.** Starting from our companion mechanism
> and letting the matter cool, the law becomes MOND. MOND's constant turns out to be the
> companion's emitted power per kilogram divided by its speed, and MOND's interpolating
> function is fixed rather than chosen (§3.15).
>
> The rev-10 next steps were also run:
> * a field equation that conserves momentum;
> * hot haloes that close the elliptical/spiral lensing gap;
> * lenses that no longer need "slip";
> * the Bullet Cluster, only partly solved and now the open problem.

> **Rev 10: a new rule, and a new law that is ours.** The owner has ruled out MOND and
> anything derived from it, plain Newtonian gravity, and dark matter, and asked us to check
> every formula against all three. The rule is recorded in [RULES.md](RULES.md). The galaxy
> law this notebook used until now *is* MOND, so it is retired below and kept visible.
> Its replacement, **hot-companion gravity** (§3.14), comes from this project's own ideas:
> * **Galaxies:** beats MOND, 15.93 against 16.13 km/s.
> * **Clusters:** fits them about five times better than MOND, error 0.22 against 1.06.
> * **Solar System:** stays silent.
> * **Formula check:** it is not MOND, because heat changes what it predicts.

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

> **Retired in rev 10.** The formula above is MOND with Milgrom's "simple" function. Under
> the owner's rule it cannot be our answer. It stays here as the record, and its
> replacement is §3.14.

Why would gravity do that? In this project's picture, gravity is not simply glued to
matter. It is **produced, carried outward, and re-absorbed** — a stream, not a
statue. A stream that thins out as it travels reaches further than one that is
strictly tied to its source.

### Ingredient two: a stream pulls harder than a statue

> **Correction, rev 7.** The Tolman argument below was this project's original
> motivation and is kept for the record, but it does not survive scrutiny as stated:
> it needs the extra gravity to carry a stress-energy tensor, and that reading fails on
> energy by four orders of magnitude. The object that survives is the **gravitational
> slip** η, which gives the same lensing-versus-dynamics split with no energy cost. See
> §3.11 and [THEORY.md](THEORY.md).

> **Rev 9: the intuition comes back, in its correct form.** Stress really is what separates
> lensing from dynamics; that part was right. The radial component of Einstein's
> equations makes it exact (§3.13): inside any radius, lensing and dynamics can differ
> only if the extra gravity carries a radial stress there. What was wrong was the
> carrier. No substance holds the stress: it comes from how matter couples to the shape of
> spacetime.
>
> The measured slips sit *between* the rows of the table below:
> * clusters at τ = 1/η ≈ 0.70–0.79;
> * the under-bent lenses at 0.65–0.74.
>
> Both are a radial *tension* that bends light more than it pulls stars.

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

## 2b. Scoreboard (rev 10: hot-companion law)

| Problem | Status | Cost |
|---|---|---|
| Galaxy rotation curves, 149 galaxies | **Beats MOND on every split**: 15.93 against 16.13 km/s | 3 universal constants |
| Cluster masses, 12 X-COP clusters | **Fitted to about 20%**: error 0.223 against MOND's 1.062. Held-out cross-check 0.234 | the same 3 |
| Solar System | **Passes**: zero extra pull at every planet | — |
| Cluster lensing shape (Coma) | **Consistent**: 1/R at 1.7–14 Mpc | — |
| Ellipticals lens harder than spirals (KiDS-1000) | **Right direction**, a tenth to a third of the size | — |
| Strong lenses (SLACS) | **Open**: need heavy stellar populations | — |
| Not MOND, not Newton, not dark matter | **Passes the formula check** | — |
| MOND derived as the cold limit (rev 11) | **Done**: a₀ = 2 × power per kg ÷ speed | — |
| Ellipticals vs spirals (KiDS), with hot haloes (rev 11) | **Reproduced**: 0.19–0.22 dex with about 1 stellar mass of hot gas | — |
| Lensing vs star motions in 6 lenses (rev 11) | **No slip needed** (gap +0.04 ± 0.03 dex) | — |
| Momentum conservation (rev 11) | **Passes**, once written as an action | — |
| Merging clusters (Bullet, 72 collisions) | **Open**: lensing moves toward the galaxies, but not durably enough | — |

Three universal constants, all fitted to data, with nothing fitted per object:
* `a` = 6.6 × 10⁻¹¹ m/s²
* the release level, 2.8 × 10⁻¹⁰ m/s²
* the companion speed, 874 km/s

<details><summary>The rev 9 scoreboard, built on the retired MOND law</summary>

| Problem | Status | Cost |
|---|---|---|
| Galaxy rotation curves | **Solved** — beats our own 20-parameter model | 1 constant |
| Cluster lensing shape | **Solved** — beats NFW | 0 parameters |
| The energy catastrophe | **Dissolved** — it was a units error | 0 |
| Lensing vs rotation on one object | **10× better** than ordinary matter | 0 |
| ~~Cluster masses~~ | ~~Closed to within the scatter~~ — **retracted, see §3.3** | — |
| Cluster dynamics | **Open** — the classic deficit, measured at 1.6–1.9× | — |
| **Gravitational slip η ≠ 1** | **The novel claim**: measured on 6 lenses and 11 clusters | 0 |
| Gravitational-wave speed (GW170817) | **Passes by construction**, via the one slot that survives (§3.13) | 0 |
| Solar system | **Passes**: slip moves γ by 4 × 10⁻¹², bound 2 × 10⁻⁵ | 0 |

One new constant, `a₀` — Milgrom's, attributed. One novel measured quantity, the
slip `η`. No dark matter, no per-object parameters. The theory and its field
equations are in [THEORY.md](THEORY.md).

---

</details>

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

### 3.3 ~~Clusters: the factor is 2, and the measurement agrees~~ — retracted

~~At the published hydrostatic bias the gap is 1.79–1.89 against an exactly predicted
2.00. The cluster problem is not open; it is closed to within the scatter.~~

**This does not survive the theory work in [THEORY.md](THEORY.md).** The argument
used τ, the Tolman active mass of a stress-energy tensor, while elsewhere arguing that
the extra gravity is not a substance at all. Tested numerically, the substance reading
fails by four orders of magnitude on energy, so the only viable frame has no stress
tensor for τ to belong to.

What survives is narrower and cleaner. τ had been doing two different jobs. In the
lensing channel it is the **gravitational slip** η — a standard, well-defined object
that costs no energy (§3.11). In the *dynamical* channel, which is what X-ray
hydrostatic masses measure, slip cannot act at all: it is a ratio between two
potentials, not a change in either one.

And the bias correction runs the wrong way for the rescue. With no τ, the
cluster/spiral gap should be 1.00. Raw it is 1.611; correcting for the published
hydrostatic bias moves it to 1.895 — *away* from 1.

**So the cluster dynamical deficit is open.** It is the classic MOND cluster problem,
now measured precisely at **1.6–1.9×**, and this framework does not currently explain
it. That the number sits near 2 is, on present analysis, a coincidence. The measured
facts in the table stand; the interpretation does not.

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

### 3.11 What the extra gravity actually is — and what is new

> **Retired in rev 10.** The second equation below is Milgrom's QUMOND, and the slip values
> were measured against it. Under the owner's rule neither is our answer; see §3.14.

*Full derivation and verification: [THEORY.md](THEORY.md).*

The theory now has field equations. With baryons as the only source:

```
∇²Φ_N = 4πG ρ
∇²Φ   = ∇ · [ ν(|∇Φ_N|/a₀) ∇Φ_N ]        what stars feel
Ψ     = Φ_N + η (Φ − Φ_N)                 lensing partner
```

The second equation reproduces the fitted law to **4.4 × 10⁻¹⁶** — machine precision.
It is Milgrom's quasi-linear MOND, and we say so.

The third is ours. Light is bent by (Φ + Ψ)/2, so the extra gravity enters lensing at
weight **(1 + η)/2**. Standard relativistic MOND is built to give η = 1 everywhere, so
lensing tracks dynamics. We measure it is not 1:

| η | Where |
|---:|---|
| 0.82 – 0.88 | three over-bent SLACS lenses |
| 1.36 – 1.54 | three under-bent SLACS lenses |
| 1.27 – 1.42 | eleven X-COP clusters, read from the published hydrostatic bias (§3.12) |

**This is the claim that survives, and it is a good one.** It has no internal
contradiction, costs no energy, is not a rediscovery, and predicts something standard
theory does not: *lensing and dynamics disagree, by an amount you can measure object by
object.*

It also survives GW170817. That event killed TeVeS and most relativistic MOND theories,
but a class with gravitational waves on the light cone exists — which constrains the
slip to live in the scalar sector, and tells us where to build. §3.13 builds there.

### 3.12 Two independent routes to the same slip

*Scripts: `research_work/results/cluster-slip/` and `research_work/results/groups-slip/`.*

**A correction first.** Revision 7 said galaxy groups show weak-lensing masses about
50% above dynamical ones. The paper cited says no such thing — the only "50%" in it is
the fraction of galaxies that live in groups. The figure came from a search-engine
summary, not the source. It has been removed everywhere. The paper's real conclusion
is that group velocity dispersions read *low*.

**Clusters.** X-ray hydrostatic masses measure what matter feels; weak lensing
measures what light feels. By the slip equation their ratio is
`1 + f(η − 1)/2`, where f is the non-Newtonian share of the mass — which the measured
gas fractions fix at 0.79–0.87. So the hydrostatic bias everyone already publishes is
a slip reading:

| Published bias | Implied η |
|---:|---:|
| 0.10 | 1.25 – 1.28 |
| 0.125 | 1.33 – 1.36 |
| 0.15 | 1.40 – 1.45 |

**The X-COP bias implies η = 1.27–1.42. The under-bent strong lenses gave 1.36–1.54
by a completely different route** — X-ray gas on megaparsec scales against stellar
orbits on kiloparsec scales. They overlap.

The field explains that bias with non-thermal pressure, and in one comparison the two
readings are degenerate. They separate cleanly, though: pressure should concentrate the
bias in disturbed, merging clusters; slip should track the non-Newtonian share and
ignore dynamical state. **Splitting a lensing-and-X-ray sample into relaxed and
disturbed halves at matched f is now the sharpest test in the programme**, and the
catalogues are public.

**Groups.** Three SL2S groups carry real weak-lensing masses. At face value all three
point to η above 1 (3.1–4.7), and one excludes η = 1 at P = 0.05 — but a 30%
dispersion underestimate, which the source paper says is present, erases all of it.
Suggestive, not established.

### 3.13 Where the slip lives: a relativistic home that survives GW170817

> **Retired in rev 10.** AeST is a relativistic MOND theory, so this "home" is ruled out by
> the owner's rule. The identity that slip equals radial stress is ordinary general
> relativity and still stands. The hot-companion law will need its own relativistic form.

*Script: `research_work/results/relativistic-slip/`. Full working in
[THEORY.md](THEORY.md). Every identity is checked symbolically. The sources were read in
the papers themselves this time, not in summaries.*

**Slip is directional gravity, exactly.** Take any theory at all. The radial part of
Einstein's equations says that inside a radius R, the mass that moves stars and the mass
read from the shape of space differ by exactly

```
M_dynamics − M_space  =  4π R³ × (radial stress at R) / c²
```

So lensing and dynamics can disagree *only* if the extra gravity pushes or pulls
differently along the radius than across it. Light bent more than stars are pulled
(η > 1) means a radial *tension*. This is the "gravity behaves differently in different
directions" idea, and it turns out not to be optional: it is what slip *is*.

In a toy galaxy with η = 1.34, the tension is 13% of the extra gravity's energy density at 1 kpc
and 25% far out. The sideways stress falls to zero there, so far out it points purely
along the radius. The equation itself is textbook, and Faber & Visser used it in 2006 to
propose weighing dark matter's pressure. We use it with no dark matter at all.

**Four places the extra gravity could live, and one survives.**

| Where it lives | Light feels it at | GW170817 | Our data |
|---|---:|---|---|
| A uniform stretch of spacetime (1984's relativistic AQUAL) | 0% | passes | **excluded**: we see 91–127% |
| A separate metric for light and matter (TeVeS) | 100% | **fails** | — |
| One shared metric (AeST, 2021) | 100% | passes | close, misses by ±0.4 |
| **AeST + a small uniform stretch felt only by matter** | **1/(1+κ)** | **passes** | **fits** |

*Light feels it at* is relative to the pull the same extra gravity puts on stars.

Why the stretch is invisible to light: light doesn't care about the overall scale of
spacetime, only its shape. So light and gravitational waves still travel together, at
the same speed and with the same delay. That is exactly what GW170817 checks. The slip
hides where only slow-moving matter can see it.

**What the numbers say:**
* **Clusters:** the matter-only part pushes outward by 12–17%, so gas and galaxies feel less
  of the extra gravity than light does.
* **Three under-bent lenses:** it pushes outward by 15–21%.
* **Three over-bent lenses:** it pulls inward by 6–10%.

**Checks it passes today:**
* GW170817 speed and delay: passes by construction.
* The solar system: slip moves Cassini's γ by at most 3.7 × 10⁻¹², against a bound of 2.3 × 10⁻⁵.
* **The one published direct cluster slip measurement** (Pizzuti et al., MACS J1206, galaxy
  orbits against lensing, immune to gas pressure) gives 1.01 +0.31/−0.28. We predict
  1.21–1.37. That is consistent at about 1σ.

**What it does not do yet: set κ.** One value cannot give a push in clusters and a pull
in three lenses. There is a promising lead, though. The KiDS-1000 survey (Brouwer et al. 2021) found that
early- and late-type galaxies of the same stellar mass bend light differently at ≥6σ,
which the authors note no universal law of gravity can produce. A coupling that depends
on galaxy type can. That test is public and next in line.

### 3.14 A new rule, and a law that is ours: hot-companion gravity

*Full working: [research_work/results/hot-companion/](research_work/results/hot-companion/README.md).
The rule: [RULES.md](RULES.md).*

**The rule.** No MOND or anything derived from it, no plain Newton, no dark matter, and
every formula checked against all three. That check is now a program,
`research_work/tools/formula_guard.py`. Pointed at MOND itself, it correctly flags it.

**What the data force.** Disk galaxies obey an observed law, the radial acceleration
relation, so tightly that we could not beat it with anything else. We tried three of our
own ideas on 149 galaxies:

| Idea | Mean velocity error (lower is better) |
|---|---:|
| A companion pool, set by radius and the galaxy's total mass | 20.35 km/s |
| A companion stream, where every piece of matter adds and nothing cancels | 18.7 – 19.4 km/s |
| MOND (for reference) | 16.13 km/s |

So in cold, isolated disk galaxies any working law has to look like MOND, just as any law
has to look like Newton in the Solar System. Our law has to be different *everywhere else*.

**The clue came from clusters.** The extra pull a cluster needs is 27 times its ordinary
mass near the centre but only 4 times at the edge. Cluster gas is hot. Hot gas lying farther
out adds to a companion field without cancelling, and that effect is 5–12 times larger near
the centre and about 2 times larger at the edge. It is the same steep shape.

**The law, in plain words.** Every piece of ordinary matter feeds a companion field.
* **Cold, orderly matter** (disk galaxies) feeds it in step. Pulls from opposite sides
  cancel, just like Newton's.
* **Hot, randomly moving matter** (cluster gas, galaxies whizzing through clusters,
  ellipticals) feeds it out of step. Nothing cancels, and it feeds harder the faster it
  moves, measured against the companion's own speed of **874 km/s**.
* **The companion pulls with the square root of its strength.**
* **Strong gravity holds it back.** It only acts where gravity is weak, which is why the
  Solar System can't see it.

```
g  =  g_N  +  exp(−g_N / g_d) · √( a · ( g_N + S_hot ) )
S_hot  =  G ∫ (3σ²/u²) ρ / d²     ← hot matter, counted without cancelling
```

This is the project's own companion idea: energy (heat) turning into gravity, spreading
without cancelling, attaching and detaching. It is not a correction bolted onto MOND.

**How it does on real data:**

| | Hot companion (ours) | MOND | Newton | Dark matter (NFW) |
|---|---:|---:|---:|---:|
| 149 galaxies, velocity error | **15.93 km/s** | 16.13 | 45.58 | 7.52 |
| 12 clusters, error in log mass | **0.223** | 1.062 | 2.211 | 0.101 |
| Adjustable numbers | **3 in total** | 1 | 0 | 298 + 24 |

* **Held-out test.** Fit the speed on six clusters and predict the other six: the error is
  0.234, across all 924 ways of splitting them. The speed stays at 874 km/s.
* **MOND's best effort on clusters.** Refitting its own constant on the clusters gets it
  to 0.338, but only with a constant nearly ten times its galaxy value. So MOND cannot fit
  both with one number, and ours does.
* **Dark matter still fits each object more closely.** But it uses two free numbers for
  every galaxy and every cluster. Ours uses three numbers for everything.

**Where it is still open, and hopeful.**
* **Ellipticals versus spirals (KiDS-1000).** Ellipticals bend light at least 0.2 dex
  more than spirals of the same stellar mass, which the survey's authors say no
  modified-gravity theory can produce. Ours leans the right way, because ellipticals are
  hot, but gets only a tenth to a third of the effect. Hot gas haloes around ellipticals
  may carry the rest; that is the next calculation.
* **The six SLACS lenses.** They need star populations 1.5–2.2 times heavier than the
  standard assumption. That is a sharp, testable prediction: spectroscopic studies of stellar
  populations can check it. The known stellar-mass convention problem (T3.5) sits here too.

**The formula check, in one line:** it matches a MOND-shaped function only for perfectly
cold, isolated matter, and even there it is not identical to any published one. Heat the
matter and the prediction moves by up to 0.34 dex at the same Newtonian pull. So it is not
MOND, not Newton, and not dark matter.

### 3.15 MOND, derived: what our law becomes when matter is cold

*Scripts and numbers: [research_work/results/hot-companion/](research_work/results/hot-companion/README.md), §8.*

You asked whether MOND could be *derived* from our cluster formulas rather than assumed.
It can, in five steps:

1. **Every kilogram of matter gives off companion energy,** a tiny 2.9 × 10⁻⁵ watts. It
   streams away at 874 km/s.
2. **Cold, orderly matter gives it off in step,** so the flows add like arrows. Their sum
   has exactly the shape of Newton's gravity. That is a geometric fact, checked to 15 digits.
3. **Heat scrambles it.** Randomly moving matter Doppler-shifts each emitter, so hot matter
   adds without cancelling. The toy confirms it: cold emitters give the arrow sum, and hot
   emitters give the plain sum.
4. **The companion's energy grows as the square of its strength,** and its pull equals its
   strength. Energy is conserved as it spreads, so its pull falls as 1/distance and scales as
   √mass. That is MOND's deep law: flat rotation curves, with v⁴ = G × M × a₀.
5. **Strong gravity holds it back,** which gives the normal Newtonian limit.

MOND takes its constant a₀ and its interpolating function as given. Here both come out:

* **a₀ = 2 × (power per kg) ÷ (companion speed).**
* **The interpolating function is fixed.** The closest published version is 0.031 dex away.

So **galaxies obey MOND because disks are cold, and clusters defy it because their gas is
hot.** MOND becomes a special case of our law, the way Newton is a special case of
Einstein's.

**One number to check in our own backyard.** The Sun, hot inside, would pour 0.2 of its
luminosity into the companion. It would lose an extra 1.4 × 10⁻¹⁴ of its mass per year,
about 15% more than sunlight and solar wind. Planetary tracking is just reaching that
precision.

**The next steps, run:**

| Step | Result |
|---|---|
| **A field equation** | Built and solved in 3D, and conservative. As first written, the heat term broke Newton's third law: a hot-plus-cold pair pushed itself along at 25% of the pair force. Rewritten as an action, momentum is conserved (the residual shrinks to 1% as the box grows). Hot cluster gas feels a small extra inward pull, about 5% at the edge. |
| **Ellipticals vs spirals (KiDS)** | Give ellipticals a hot halo of about one stellar mass at 0.6–1 keV. Our law then yields 0.19–0.22 dex extra lensing beyond 300 kpc, matching the observed ≥ 0.2. MOND, which counts the halo's mass but not its heat, gets 0.14. |
| **Lens bookkeeping (T3.5)** | Standard distances and published star masses remove about 0.13–0.16 dex of the old lens tension. Our law then needs stars 1.2–1.6× heavier than Salpeter. Dark-matter models of the same lenses need about 1.0. Spectroscopic studies of massive ellipticals report super-Salpeter stars, so there is a clean test available. |
| **Slip, re-derived** | Lensing and star motions now ask for the same stellar mass within +0.04 ± 0.03 dex. The earlier "slip" was an artefact of the MOND baseline. |
| **Bullet Cluster** | If companions keep the velocity of the matter that emitted them, lensing moves off the gas toward the galaxies: 42–75% of the way, depending on time since the collision. But in 72 stacked collisions, lensing sits within 6 ± 8 kpc of the stars (Harvey et al. 2015). Our gas rebuilds its own companion too fast to match that. **This is the open problem:** the gas must rebuild far more slowly. |

---

## 4. Why this beats the alternatives

> **Rev 10.** Under the new law the comparison is sharper:
> * **Galaxies:** it beats MOND.
> * **Clusters:** it fits them, where MOND fails by a factor of about 2.6.
> * **Dark matter:** it does all this with three numbers for everything, where dark matter
>   uses two per object.
>
> The paragraphs below describe the retired MOND-based law and are kept for the record.

**Against dark matter.** A dark-matter halo needs a new particle, plus a profile with
two or three numbers fitted *separately for every galaxy*. We have one number for all
149. And when we make rotation predict lensing on the same object without re-fitting,
ordinary matter plus a halo has to change the galaxy's stellar mass by a factor of
two to three between the two explanations. One law does not.

**Against MOND and its relatives.** On dynamics we *are* MOND. Our equation for what
stars feel is Milgrom's, and we say so. ~~They are wrong on clusters by about a factor
of two, and that gap is exactly what τ fills.~~ **Retracted (§3.3):** we do not
currently solve the cluster deficit either.

Where we differ is lensing. Relativistic MOND theories are built so that light and
matter respond to the same potential — slip η = 1. We measure η ≠ 1 on individual
galaxies, with every constant frozen. If that holds, it is a clean, falsifiable
departure from the whole MOND family, and it is something dark matter does not predict
either: a halo bends light and moves stars through one potential. And the one
relativistic MOND theory that survives GW170817 only needs a small addition to carry
it (§3.13).

**Against expanding-universe explanations.** Nothing here uses expansion; every test
is local, below z = 0.3. The roadmap recommends leaving the position out of any paper,
because no result here requires it.

---

## 5. What is next, in order

*Rev 11:*

1. **Mergers.** Find what makes the gas rebuild its companion slowly after a collision.
   Then test against the 72-collision stack.
2. **Derive the companion speed (874 km/s) and the release level.** Both are still measured.
3. **Spectroscopic star counts** for the six SLACS lenses: our law says 1.2–1.6× Salpeter,
   dark matter says about 1.0.
4. **The Sun's extra mass loss** of 1.4 × 10⁻¹⁴ per year, against planetary ranging.
5. **The KiDS offset against radius:** our law says it grows outward.

*Rev 10, for the hot-companion law:*

1. **Hot gas haloes around ellipticals**, to see whether the heat term closes the KiDS
   elliptical/spiral gap.
2. **A field equation** for the law, to settle how whole stars and planets move and to
   compute lensing without assuming light and matter feel the same pull.
3. **Fix the stellar-mass convention (T3.5)**, then test the heavy-star prediction on the
   SLACS lenses.
4. **The Bullet Cluster**, using how long the companion takes to spread after a merger shock.
5. **Re-derive the slip η** against the new law.

*The rev 9 list, for the retired law:*

*Rev 9. This list replaces the rev 6 list, which was built around τ. The old list is kept
below and struck through.*

1. **Slip from galaxy orbits in a stack of clusters.** Galaxy orbits ignore gas pressure, so
   this measures η directly. One cluster has already been done by others (§3.13); the
   CLASH-VLT sample has more.
2. **Relaxed versus disturbed clusters at matched f** (§3.12). This separates slip from gas
   pressure using public catalogues.
3. **Lensing by galaxy type (KiDS-1000).** This is the most direct handle on what sets κ.
4. **Fix the stellar-mass convention** (roadmap T3.5), then re-read the three over-bent
   lenses.
5. **Derive `a₀`.** It is still a measured constant with no explanation.

<details><summary>The rev 6 list, superseded</summary>

1. ~~Cluster weak lensing shape.~~ Done: §3.6. The prediction held.
2. ~~Separate τ from the lensing amplitude.~~ τ became slip; items 1–3 above do this job.
3. ~~Test the predicted hydrostatic bias b ≈ 0.17.~~ The τ = 2 cluster claim behind this
   number is retracted (§3.3).
4. ~~Groups.~~ Done: §3.12, suggestive but inconclusive.

</details>

---

## 6. Are we copying something? A straight answer

We check this deliberately, and the answer has two halves.

> **Rev 10.** Because of the paragraph below, the galaxy law is now retired as our answer
> ([RULES.md](RULES.md)). The hot-companion law's closest relatives we know of are:
> * **EMOND**, which raises MOND's constant in deep potentials;
> * **superfluid dark matter**, whose behaviour switches with the dark matter's temperature.
>
> Ours changes no constant and adds no dark matter. It uses the *ordinary* matter's random
> motions, counted without cancelling. We have not found it published; a full check is
> roadmap T3.3.

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

### ~~The cluster fix appears to be ours~~ — superseded in rev 9, see below

MOND's cluster problem is long-standing and well documented: a residual factor of
two to three that the programme has never resolved without adding something. The
published attempts we can find take a different route from ours — making `a₀` itself
vary with gravitational potential (EMOND), or with environment, or adding massive
neutrinos, or proposing faint stellar remnants as uncounted baryons.

**None of those is what we are proposing.** Our claim is that the missing factor is
the Tolman active-mass factor of the gravitational stream itself — `τ = 1 + w_r + 2w_t`
— whose two endpoints are fixed exactly by the stress tensor and carry no adjustable
freedom. Targeted searches did not surface that mechanism in the literature.

~~That is our novelty claim, and it is the only one we make.~~ The Tolman cluster
mechanism was retracted in rev 8 (§3.3), so this novelty claim went with it.

### What is ours now, exactly

* **Ours:** the measurement that the extra gravity's slip is **not 1**, object by
  object, with every constant frozen. That covers six lenses, plus the reading of the published
  cluster hydrostatic bias as a slip.
* **Ours, pending a full literature check:** the relativistic home. Both ingredients are
  published: AeST (Skordis & Złośnik 2021) and a conformally coupled scalar (Bekenstein &
  Milgrom 1984). We have not found the combination proposed as a source of slip.
* **Not ours:**
  * the dynamics equation (Milgrom's QUMOND);
  * the stress identity (the radial Einstein equation, used by Faber & Visser 2006);
  * the idea of measuring slip in clusters (Pizzuti et al. 2016 did it for one cluster).

### Standing honesty register

* **SPARC here has no quality cuts.** Applying the usual inclination and flag cuts
  would lower the scatter. We left them off so the comparison against R10 is like for
  like.
* **Cluster masses are hydrostatic**, with a bias comparable to the residual scatter.
  That is precisely why §3.3 is phrased as a prediction to test.
* **The project's stellar-mass rescaling moves results by factors of ~3** (§3.9).
  This is now the largest known systematic in the programme.
* ~~**τ's two endpoints are exact; the crossover between them is not.**~~ τ became
  the slip η (rev 8). **κ is measured, not derived, and one κ does not yet fit every
  system.**
* **`a₀` is measured, not derived.** No part of this work explains its value.

---

**Publication roadmap:** [PAPER-ROADMAP.md](PAPER-ROADMAP.md) sets out what stands
between this notebook and a paper, including one contradiction at the centre of the
argument that has to be resolved first.

*Scripts, data hashes and full numerical output for every claim live under
`research_work/results/`. Each experiment directory carries its own README with the
protocol, the limitations, and the raw JSON.*
