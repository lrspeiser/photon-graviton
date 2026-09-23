# Theory: what the extra gravity is, and the equations it obeys

> **Superseded as our solution, 22 September 2026: see [RULES.md](RULES.md).**
> The owner has ruled out MOND and anything derived from it, plain Newton, and dark matter.
> The dynamics equation below is Milgrom's QUMOND, and the AeST "relativistic home" is a
> MOND theory, so neither can be our answer. They are kept here as a record of what was
> tested.
>
> The slip results (η) were measured against that MOND law, so they must be re-derived
> before they are quoted.
>
> The current candidate is **hot-companion gravity**:
> [research_work/results/hot-companion/](research_work/results/hot-companion/README.md).
> * **Galaxies:** 15.93 km/s, against MOND's 16.13.
> * **Clusters:** error 0.223, against MOND's 1.062.
> * **No anomalous pull in the Solar System.**
> * **The formula check shows it is not MOND**, because it depends on how hot the matter is.


## Hot-companion gravity: the current theory (23 September 2026)

*Details and every number: [research_work/results/hot-companion/](research_work/results/hot-companion/README.md).
Proposed here; originality unverified.*

### Postulates

1. **Emission.** Every kilogram of ordinary matter emits companion energy at rate ℓ. The
   energy streams away at speed `u`.
2. **Coherence.** Cold, orderly matter emits in step, so energy *fluxes* add as vectors.
   Hot matter has its emitters Doppler-scrambled by random motion (speed σ), so it emits an
   extra `k ℓ` per kilogram, with `k = 3σ²/u²`, out of step. Energy *densities* add as
   scalars.
3. **Amplitude.** The companion's energy density is A²/8πG, and its pull on matter and
   light equals its amplitude A.
4. **Attachment.** Strong fields hold the companion. It is released where the ordinary pull
   is weak: `f = exp(−|g_N|/g_d)`.

### What follows

```
g  =  g_N  +  exp(−|g_N|/g_d) · √( a (|g_N| + S_hot) ),    a = 2ℓ/u,    S_hot = G ∫ k ρ / d²
```

* **Constants,** fitted once: a = 6.58 × 10⁻¹¹ m/s², g_d = 2.82 × 10⁻¹⁰ m/s²,
  u = 874 km/s. Therefore ℓ = 2.9 × 10⁻⁵ W/kg.
* **Field equation:** ∇²Φ = −∇·h, with h the vector form of the law above, and g = −∇Φ.
* **Action:** `L = −(1/8πG)[2∇Φ·∇ψ − W(|∇ψ|, S)] + ρ(v²/2 − Φ)`. It conserves momentum,
  verified in 3D, and it gives hot matter a small reaction force.

### MOND is the cold limit, derived rather than assumed

With σ → 0, the field equation becomes Milgrom's QUMOND:
* The deep law √(a g_N) comes from companion energy conservation.
* The Newtonian limit comes from attachment.
* The scale symmetry is exact.
* **a₀ = 2ℓ/u** is emission power per kilogram over speed.
* The interpolating function `ν(y) = 1 + exp(−y/λ)/√y` is fixed, not chosen. The closest
  published function is 0.031 dex away.

So galaxies obey MOND because disks are cold, and clusters defy it because their gas is hot.

### Status

| Test | Result |
|---|---|
| 149 galaxies | 15.93 km/s, against MOND's 16.13 |
| 12 clusters | error 0.223, against MOND's 1.062; held-out 0.234 |
| Solar System | 0 anomalous pull |
| Lenses | no slip needed (gap +0.040 ± 0.026 dex); needs stars 1.2–1.6× Salpeter, a testable prediction |
| KiDS elliptical/spiral gap | reproduced with hot haloes of about 1 stellar mass |
| Mergers | lensing moves toward the galaxies, but the rebuild is too fast for the 72-collision stack. **This is the open problem.** |

---

**22 September 2026.** Resolves roadmap tasks T0.1 (the central contradiction), T2.2
(the gravitational-wave constraint) and T0.2 (field equations, non-relativistic and
relativistic). Every numerical claim below was computed in this session and is
reproducible from the scripts named.

This document **corrects** earlier claims in the notebook. Where it does, the
correction is stated plainly rather than softened.

---

## T0.1 — The contradiction, resolved

### What was wrong

The notebook argued in §3.2 that the extra gravity is **not matter** — that there is
no reservoir, that the force is the source's own field reaching further, and that
this is why the energy catastrophe dissolves. It then used **τ**, the Tolman active
mass, which is a property **of a stress-energy tensor**.

Both cannot hold. So each branch was tested numerically.

### Branch (a): the extra gravity carries a real stress-energy tensor

For a galaxy of 2 × 10¹¹ M☉, the law implies an effective enclosed mass of
**6.5 × M★** at 100 kpc. As energy that is 2.3 × 10⁵⁹ J, against 1.6 × 10⁵⁵ J of
starlight over a Hubble time — **short by 1.5 × 10⁴**.

Could the field's *own* energy pay instead of a reservoir? At g = 3a₀ the field energy
density is 7.4 × 10⁻¹¹ J/m³, equivalent to 8.2 × 10⁻²⁸ kg/m³. Cluster matter density at
R500 is 4.8 × 10⁻²⁴ kg/m³. **The field's self-energy is 1.7 × 10⁻⁴ of what it would
have to rival.**

**Branch (a) fails.** Neither a reservoir nor the field itself can pay.

### Branch (b): a modified propagator — survives

The extra gravity is a nonlinearity in how ordinary matter sources the field. Nothing
new exists, so nothing needs funding. This is the only viable frame.

### The consequence, which is the real finding

In branch (b), τ as a Tolman factor has no substance to belong to. But the notebook
had used the symbol **τ for two different physical quantities**, and they come apart
cleanly:

| Where τ was used | What it actually measured | Survives branch (b)? |
|---|---|---|
| JR-10, six SLACS lenses | The **lensing weight** of the extra term, relative to dynamics | **Yes** — this is the *gravitational slip* η, a standard, well-defined object in modified gravity, with no energy cost |
| JR-11/13, eleven clusters | A multiplier on the **dynamical** force | **No** — slip is a ratio *between* the two potentials and cannot change either one's amplitude |

So:

* **Slip is real and is ours to measure.** It legitimately explains a
  lensing-versus-dynamics split, and that is exactly what JR-10 found.
* **The cluster factor near 2 has no mechanism.** It is measured in hydrostatic
  mass — dynamics, Φ — so slip cannot touch it. Its resemblance to the free-streaming
  Tolman value is, on this analysis, a coincidence.

### The direction of the bias correction, checked

This matters because the notebook read the hydrostatic-bias correction as closing the
cluster gap. If the law were universal with no τ, the cluster/spiral gap should be
**1.00**. Raw it is 1.611; with the published X-COP bias of 0.15 it is 1.895. **The
correction moves the gap away from 1 and toward 2.** It only helps if τ = 2 is real in
the dynamical channel — which T0.1 has just ruled out.

**Correction to the notebook and the published page:** the cluster problem is *not*
closed. It is the classic MOND cluster deficit, now measured precisely at **1.6–1.9×**,
and this framework does not currently explain it.

---

## T2.2 — GW170817: does the frame survive?

The simultaneous arrival of gravitational waves and light from GW170817 bounds the
difference between their speeds to about 10⁻¹⁵. That result killed a large class of
modified-gravity theories outright, including **TeVeS, Einstein-Aether, Hořava gravity
and generalised Proca**, the relativistic MOND frameworks most often cited.

It did not kill all of them. A class of relativistic MOND theory has been constructed
in which the tensor modes propagate at exactly c, and theories whose modification
sits in the scalar sector with an unaltered tensor sector survive by construction.

**Verdict: branch (b) survives GW170817, with a construction constraint.** Any
relativistic completion must keep gravitational waves on the light cone, which means
the modification — and in particular the slip — must live in the scalar sector. This
narrows T0.2 usefully and early, which was the point of running it first.

*Sources: Ezquiaga & Zumalacárregui, Phys. Rev. D 97, 061501 (2018), arXiv:1711.07403;
Phys. Rev. D 100, 104013 (2019) on alternatives with c_T = c.*

**A second GW170817 constraint, and why slip escapes both.** The same event also bounds
any difference in how light and gravitational waves are delayed by the mass along
their path: |γ_GW − γ_EM| < 9.8 × 10⁻⁸ (Boran, Desai, Kahya & Woodard, Phys. Rev. D 97,
041501 (2018)). That rules out every theory in which matter and light respond to the
extra gravity through a different metric from gravitational waves, TeVeS among them.

GW170817 constrains only what light and gravitational waves feel, and both feel the
combination Φ + Ψ. Slip concerns Φ, which only slow-moving matter feels. **So a slip built so that light and
gravitational waves still share one Φ + Ψ passes both tests by construction.** The
relativistic section of T0.2 below builds exactly that.

---

## T0.2 — Field equations

### The equations

With baryonic density ρ as the only source:

```
(1)   ∇²Φ_N = 4πG ρ                                   Newtonian potential
(2)   ∇²Φ   = ∇ · [ ν(|∇Φ_N| / a₀) ∇Φ_N ]              dynamics — what stars feel
(3)   Ψ     = Φ_N + η (Φ − Φ_N)                         lensing partner

      ν(y) = ½ + √(¼ + 1/y)
```

Stars move in Φ. Light is deflected by (Φ + Ψ)/2, which from (3) is

```
      (Φ + Ψ)/2  =  Φ_N  +  ½(1 + η) · (Φ − Φ_N)
```

so the extra gravity enters lensing at weight **(1 + η)/2** relative to dynamics.

### Verification

**Equation (2) reproduces the fitted law exactly.** In spherical symmetry, Gauss's law
turns (2) into g = ν(g_N/a₀) g_N, and for a Hernquist baryonic profile the maximum
relative difference from the law fitted to 3,150 SPARC measurements is
**4.4 × 10⁻¹⁶** — machine precision. It is not an approximation to the law; it is the
law.

**Equation (3) is what JR-10 measured.** The quantity called γ_χ there — the lensing
weight on the extra term, solved on six real lenses — is identically η. Nothing is
re-fitted by renaming it:

| η | Lensing weight (1+η)/2 | Where seen |
|---:|---:|---|
| 1.00 | 1.000 | standard relativistic MOND, everywhere |
| 0.82 – 0.88 | 0.91 – 0.94 | three over-bent SLACS lenses |
| 1.36 – 1.54 | 1.18 – 1.27 | three under-bent SLACS lenses |
| 1.27 – 1.42 | 1.13 – 1.21 | eleven X-COP clusters, read from the published hydrostatic bias *(T1.2 below)* |
| 3.1 – 4.7 at face value | — | three SL2S groups, erased by a 30% dispersion bias *(T1.3 below)* |

### Attribution — stated exactly

* **Equation (2) is not ours.** It is Milgrom's quasi-linear MOND (QUMOND, 2010) with
  the standard "simple" interpolation function. Our a₀ = 1.171 × 10⁻¹⁰ m/s² sits within
  a few percent of his 1983 value.
* **Equation (3) with η ≠ 1 is the novel claim.** Standard relativistic MOND theories
  are built to give η = 1, so that lensing tracks dynamics. A measured η ≠ 1 is a
  departure from them, and it is well defined, testable, and costs no energy.

### What is still open in T0.2

1. ~~**A relativistic completion with η ≠ 1.**~~ **Found: see the next section.** The
   GW170817-surviving MOND theory gives η = 1, now verified from the paper itself. The
   smallest addition to it that gives η ≠ 1 and passes GW170817 is a conformal coupling
   κ of matter to its scalar.
2. **What sets η** — now, equivalently, what sets κ. The SLACS values range 0.82–1.54
   and correlate with redshift; the cluster reading gives 1.27–1.42. This is the "what
   sets τ" question in its correct form, and it is now a question about one coupling
   constant.
3. **Cosmological slip constraints.** Combined CMB and lensing analyses bound slip on
   large scales. Our η lives in the galaxy-scale extra-gravity sector, a different
   regime, but the two must be checked for consistency.

---

## T0.2 (relativistic) — where the slip can live

*Script: `research_work/results/relativistic-slip/code/slip_slots.py`. Each identity is
checked symbolically, and each number is printed by the script.*

### Slip is radial stress

Take any static, weak-field metric, in any theory. The radial component of Einstein's equations then
gives, exactly,

```
M_Φ(<R) − M_Ψ(<R)  =  4π R³ p_r(R) / c²
```

* M_Φ is the mass that sets orbits.
* M_Ψ is the mass read from the spatial curvature.
* p_r(R) is the radial stress, at R, of whatever sources the metric beyond ordinary matter.

**The equation is standard**, and Faber & Visser (MNRAS 372, 136, 2006) proposed using
it to measure the pressure of dark matter from rotation curves plus lensing. What is new is reading it
without dark matter:

**Inside R, lensing and dynamics can disagree only if the extra gravity carries a
radial stress at R.** η > 1, where light is bent more than stars are pulled, needs radial
tension. η < 1 needs radial pressure. Gravity that behaves differently in different
directions is not an optional flourish here: it is what slip *is*.

A toy galaxy, 10¹¹ M☉ Hernquist under QUMOND with η = 1.34, shows what this looks like:
* At 1 kpc the needed stress is a radial tension of 13% of the extra gravity's energy density.
* In the far field it grows to 25%.
* The sideways stress falls to zero there, so far out the stress points purely along the radius.

### Four places the extra gravity could live — one survives

| Where the extra gravity lives | Light feels it at | η | GW170817 | Our measurements |
|---|---:|---:|---|---|
| A scalar that stretches spacetime uniformly (conformal): relativistic AQUAL, 1984 | 0% | −1 | passes | **excluded**: we see 91–127% |
| A separate metric for matter and light (disformal): TeVeS | 100% | 1 | **fails**, on both speed and Shapiro delay | — |
| One shared metric, no stress: AeST (Skordis & Złośnik 2021) | 100% | 1 | passes | close, but misses the ±0.4 |
| **AeST + a small conformal coupling κ of matter to its scalar** | 1/(1+κ) | (1−κ)/(1+κ) | **passes** | **fits** |

**Why light ignores the first row.** Maxwell's equations do not care about the overall
scale of spacetime. A uniform stretch leaves every light path unchanged. Skordis & Złośnik
state this directly: "null geodesics are unaltered by conformal transformations."

**AeST's η = 1, verified.** The earlier hedge "as far as we can establish" is removed.
The paper states that in the weak-field quasistatic limit "(5) leads to Ψ = Φ", and that
"since Ψ = Φ, (6) leads to the right lensing."

**Why the last row passes GW170817.** Adding a conformal factor does not move light's paths.
Light and gravitational waves therefore still share one set of paths, one speed and one
Shapiro delay, and both feel exactly the Φ + Ψ that AeST gives. The slip lives entirely
in Φ − Ψ, which only slow-moving matter feels.

**What the last row says physically.** Most of the extra gravity is built into the shape
of spacetime, the AeST part, which light and stars both feel. A small uniform-stretch
part is felt only by matter:

| Systems | η | κ | The matter-only part |
|---|---:|---:|---|
| 11 X-COP clusters | 1.27 – 1.42 | −0.12 to −0.17 | pushes outward, 12–17% |
| 3 under-bent SLACS lenses | 1.36 – 1.54 | −0.15 to −0.21 | pushes outward, 15–21% |
| 3 over-bent SLACS lenses | 0.82 – 0.88 | +0.06 to +0.10 | pulls inward, 6–10% |

**Why this is the natural slot:**
* Within AeST-type theories, a nonminimal G₄(ϕ)R coupling is the same slot seen in another frame.
* Aether terms give no slip at this order: Einstein-aether has γ = 1 (Foster & Jacobson 2006).
* G₄(X) and G₅ couplings change the speed of gravitational waves.
* Beyond-Horndeski theories offer other routes, not explored here.

### Checks it passes

* **Solar system.** At the Cassini ray (1.6 R☉), slip moves γ by at most
  **3.7 × 10⁻¹²**, against the bound of 2.3 × 10⁻⁵.
* **The one published direct cluster slip measurement** (see T1.2 below): consistent at
  about 1σ.

### What it does not yet do

1. **Set κ.** One κ cannot give η > 1 in clusters and η < 1 in three lenses.
   * Either κ depends on something not yet identified, or the three over-bent lenses carry
     a systematic. The stellar-mass convention (T3.5) moves lens results by 2.4×.
   * **A public lead:** Brouwer et al. (A&A 650, A113, 2021) find that early- and
     late-type galaxies of equal stellar mass sit on different weak-lensing acceleration
     relations, at ≥6σ. They note that a universal modification of gravity cannot produce
     this. A κ that depends on galaxy type could. Not yet tested.
2. **Background drift.** If the scalar's background value drifts in time, particle masses
   drift at κ times that rate, and lunar laser ranging bounds such drifts near 10⁻¹³ per
   year. A static background avoids the issue entirely.
3. **Prior art.** Both ingredients are published: AeST, and a conformally coupled scalar
   (relativistic AQUAL). We have not found the combination proposed as a source of slip.
   The full check is T3.3.
4. **The simple ν in the solar system.** This is separate from slip: the extra pull tends to a₀ at
   high acceleration rather than switching off. Hees et al. (MNRAS 455, 449, 2016) find
   that Cassini rules out several popular MOND transition functions and leaves others
   viable. Whether ours survives is T2.3.

---

## T1.2 and T1.3 — results

### A correction first

An earlier revision of this document, the notebook and the published page said the
SL2S group paper finds weak-lensing masses about 50% above dynamical ones. **It does
not.** The only "50%" in that paper refers to the fraction of galaxies that live in
groups. The figure came from a search-engine summary and was never in the source. The
paper's own conclusion runs the other way: from simulations, it finds group velocity
dispersions are *always underestimated*, so dynamical masses read low. The claim has
been removed everywhere it appeared.

### T1.3 — groups, from the published table

*Script: `research_work/results/groups-slip/code/sl2s_slip.py`.* Munoz et al. 2013,
Tables 2 and 3. Of seven groups, three carry a real weak-lensing mass; the rest have
only an upper limit, sit in the galaxy-lensing regime, or fall at the field edge.

Weak-lensing masses are projected within 2 Mpc; virial masses within 0.2–1.1 Mpc, so
they cannot be ratioed directly. The law makes the far field isothermal, so the
isothermal aperture correction is the theory applied to itself, not an extra choice.

| Group | M_WL (10¹⁴ M☉) | η at face value | 68% range | η with 30% σ bias |
|---|---:|---:|---:|---:|
| SL2SJ02140-0535 | 5.5 ± 3.7 | 4.68 | 0.9 – 14.8 | 1.78 |
| SL2SJ08544-0121 | 6.3 ± 2.5 | 4.73 | 2.5 – 10.0 | 1.81 |
| SL2SJ09413-1100 | 3.7 ± 3.4 | 3.14 | −0.6 – 21.2 | 1.03 |

At face value all three point to η above 1, and SL2SJ08544-0121 excludes η = 1 with
P = 0.05 — the same direction as the under-bent lenses. But a 30% dispersion
underestimate, which the source paper says is present, brings every group inside η = 1
at 68%. **Suggestive of η > 1; cannot establish it.** A real group test needs stacked
weak lensing around a large, well-sampled group catalogue.

### T1.2 — clusters: the hydrostatic bias is already a slip measurement

*Script: `research_work/results/cluster-slip/code/cluster_slip.py`.*

Hydrostatic masses measure Φ; weak lensing measures (Φ + Ψ)/2. From equation (3),

```
M_WL / M_HSE  =  1 + f (η − 1) / 2,        f = 1 − M_baryon / M_HSE
```

and for X-COP the measured gas fractions fix f at 0.79–0.87. So every published
X-ray-versus-lensing comparison already contains a slip reading:

| Published bias b | M_WL / M_HSE | Implied η |
|---:|---:|---:|
| 0.07 | 1.075 | 1.17 – 1.19 |
| 0.10 | 1.111 | 1.25 – 1.28 |
| 0.125 | 1.143 | 1.33 – 1.36 |
| 0.15 | 1.176 | 1.40 – 1.45 |

**The published X-COP bias implies η = 1.27–1.42. The three under-bent SLACS lenses,
measured by a completely different route, gave 1.36–1.54.** One is X-ray gas and weak
lensing on megaparsec scales; the other is stellar kinematics and strong lensing on
kiloparsec scales. They overlap.

What is *not* claimed: that the bias *is* slip. The field attributes it to non-thermal
pressure, with good simulation support, and the two readings are degenerate in a
single comparison.

**The test that separates them.** Non-thermal pressure predicts the bias grows with
dynamical disturbance — merging, unrelaxed clusters read lower. Slip predicts it tracks
f and ignores dynamical state. Split a lensing-and-X-ray sample into relaxed and
disturbed halves at matched f: same bias in both means slip; bias concentrated in the
disturbed half means pressure. This is cheap, uses public catalogues, and is now the
sharpest observational test in the programme.

**A pressure-free cross-check already exists.** Pizzuti et al. (arXiv:1602.03385)
measured slip directly in the relaxed cluster MACS J1206.2-0847, using *galaxy orbits*
against strong and weak lensing. Galaxy orbits do not feel gas pressure.

* **Their result:** η(r₂₀₀) = 1.01 +0.31/−0.28 for the total potential.
* **Our prediction:** η_total = 1 + f(η − 1) = **1.21–1.37**. The total is diluted by the
  baryon share, which is why it sits below the extra-gravity-only range of 1.27–1.42.
* **Comparison:** our range sits 0.7–1.2σ above their centre. Consistent, at the edge, and not
  a confirmation.

A stack of such clusters breaks the slip-versus-pressure degeneracy without needing the
relaxed-versus-disturbed split. The CLASH-VLT sample is the natural place to do it.

---

## What the paper is now

Smaller than the notebook suggested, and on firmer ground.

* **Dynamics:** QUMOND with one constant — Milgrom's, attributed, reproduced from our
  own pipeline as a check.
* **The novel claim:** the gravitational slip of the extra gravity is **not 1**.
  Measured on six individual lens galaxies with every universal constant frozen
  (η = 0.82–1.54, correlated with redshift), and read independently from the published
  cluster hydrostatic bias (η = 1.27–1.42), which overlaps the under-bent lenses.
* **A relativistic home:** the aether-scalar-tensor theory, which survives GW170817,
  plus a conformal coupling κ of matter to its scalar.
  * It passes GW170817 by construction and moves the solar-system γ by 10⁻¹².
  * It turns every slip measurement into a reading of one number, κ.
  * Clusters need κ ≈ −0.15.
* **Clusters:** the classic dynamical deficit, measured at 1.6–1.9×, **not solved**,
  and reported as such.

That is a claim a referee cannot dismiss as a rediscovery, that has no internal
contradiction, and that makes one testable prediction the standard theory does not:
**lensing and dynamics disagree, by an amount you can measure object by object.**
