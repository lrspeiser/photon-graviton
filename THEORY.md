# Theory: what the extra gravity is, and the equations it obeys

**22 September 2026.** Resolves roadmap tasks T0.1 (the central contradiction), T2.2
(the gravitational-wave constraint) and T0.2 (field equations). Every numerical claim
below was computed in this session and is reproducible from the commands shown.

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

1. **A relativistic completion with η ≠ 1.** The GW170817-surviving MOND class gives
   η = 1 in the quasi-static limit, as far as we can establish. Getting η ≠ 1 needs
   anisotropic stress in the scalar sector — a specific construction task, not a
   conceptual obstacle, but not done.
2. **What sets η.** The SLACS values range 0.82–1.54 and correlate with redshift; the
   cluster reading gives 1.27–1.42. This is the "what sets τ" question in its correct
   form, and it is now a question about a well-defined object rather than a borrowed one.
3. **Cosmological slip constraints.** Combined CMB and lensing analyses bound slip on
   large scales. Our η lives in the galaxy-scale extra-gravity sector, a different
   regime, but the two must be checked for consistency.

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

---

## What the paper is now

Smaller than the notebook suggested, and on firmer ground.

* **Dynamics:** QUMOND with one constant — Milgrom's, attributed, reproduced from our
  own pipeline as a check.
* **The novel claim:** the gravitational slip of the extra gravity is **not 1**.
  Measured on six individual lens galaxies with every universal constant frozen
  (η = 0.82–1.54, correlated with redshift), and read independently from the published
  cluster hydrostatic bias (η = 1.27–1.42), which overlaps the under-bent lenses.
* **Clusters:** the classic dynamical deficit, measured at 1.6–1.9×, **not solved**,
  and reported as such.

That is a claim a referee cannot dismiss as a rediscovery, that has no internal
contradiction, and that makes one testable prediction the standard theory does not:
**lensing and dynamics disagree, by an amount you can measure object by object.**
