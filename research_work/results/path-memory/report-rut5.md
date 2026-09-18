# RUT-1 stage 5: the lag did not cause the instability — it slowed it by a factor of fifty, and the mode is now predicted, not just observed

[protocol-rut5.md](protocol-rut5.md), declared in e354e65 before any of this ran. **Stage 4's attribution
was not controlled, and the owner supplied the control that fixes it.** Their zero-delay field — same
kernel, same static gain, no lag — reproduces here at an m = 2 amplitude e-folding time of **0.1907
reference periods** against their 0.191, and it is unstable at *every* angular mode.

The rest of this stage is what that control makes possible: a linear mode analysis of the same equilibrium
for every rung of the response ladder, which now **predicts stage 4's instability** rather than merely
recording it.

## The headline number

| | m = 2 amplitude e-folding |
|---|---|
| linear theory, exact history, no grid | **8.579 T0** |
| simulation, grid spacing w/5 → w/10 → w/20 | 8.61 → 8.55 → 8.47 T0 |
| simulation, timestep 0.01 → 0.005 → 0.0025 | 8.61 → 8.52 → 8.60 T0 |
| stage 4's measured band growth | ≈ 8.7 T0 |

The instability stage 4 found is real, it is the linear m = 2 mode of the supported ring, and it is
converged under four-fold refinement of both the grid and the timestep. It is not a numerical artefact and
it is not the discreteness of the writers.

## The ladder: what the lag actually does

Thirty-two writers, w/R = 0.2, every rung carrying the same 9.99% equilibrium inward support.

| rung | fastest mode | e-folding | inertial pattern speed |
|---|---|---|---|
| no attraction (Kepler ring) | — | neutral to 6×10⁻¹⁵ | — |
| **two-stage** (τ_form = 3 T0, τ_keep = 10 T0) | m = 2 | **8.58 T0** | **0.0005 Ω** |
| **one-stage** (τ_form = 0) | m = 5 | **4.45 T0** | 0.993 Ω |
| **instantaneous**, matched gain and support | m = 4 | **0.155 T0** | — |
| **pairwise Newtonian**, matched support | m = 9 | **0.132 T0** | — |

Read the ordering. Removing the lag entirely, at the same equilibrium support, makes the ring **fifty-five
times more unstable**, and the softened Gaussian kernel with no lag behaves like ordinary pairwise
attraction at the same strength — the classical instability of a ring of co-orbiting attracting bodies.
**So the lag is not what destabilises the ring. Collective attraction destabilises it, and the lag is the
only reason it survives long enough to be interesting.** Each stage of lag buys about a factor of two in
rate: zero → one → two stages runs 0.155 → 4.45 → 8.58 T0.

The mechanism is visible in the matrices: a body's displacement is written along its whole orbit, so a
lagged field responds to a Bloch-m displacement with the m ≠ 0 structure largely averaged away, while an
instantaneous field feels the full discrete structure. That is the same collective cancellation stage 2C
measured for the drag, acting here on the stability.

**The two-stage mode is nearly stationary in the inertial frame** — pattern speed 0.0005 Ω, against the
one-stage mode's 0.993 Ω, which corotates. A memory field is written in space and decays in place, so the
pattern it crystallises does not have to travel with the matter. This is the measurement stage 4 said was
missing; it is answered analytically here, and it is a sharp prediction for the next campaign to test
directly on the force-producing field.

## It is not discreteness: the rate converges as the writers are subdivided

At fixed **total** writing rate, w/R = 0.2, fastest mode over m ≤ 8:

| writers | 8 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|
| two-stage e-folding (T0) | 1.22 | 0.61 | **8.58** | **8.58** | **8.58** |
| instantaneous e-folding (T0) | 0.331 | 0.107 | 0.155 | 0.155 | 0.155 |

Eight or sixteen writers is a lumpy ring — the spacing 2π/8 = 0.79 is four times the footprint width — and
it is much more unstable, as a ring of discrete blobs should be. By 32 writers the spacing 0.196 has
dropped below the width 0.2 and the rate is **identical at 32, 64 and 128**: the mode is a property of the
supported ring, not of how finely the same source is divided. That is the limit the stage 4 configuration
sits in.

## Whether a weaker ring is stable

Thirty-two writers, w/R = 0.2, fastest mode over m ≤ 8:

| label fraction | 0 | 0.025 | 0.05 | 0.10 | 0.20 |
|---|---|---|---|---|---|
| support | 0 | 2.5% | 5.0% | 9.99% | 20.0% |
| two-stage e-folding (T0) | ∞ | 28.2 | 11.7 | **8.58** | 2.45 |
| instantaneous e-folding (T0) | ∞ | 0.215 | 0.180 | 0.155 | 0.076 |

Growth is monotonic in the support, with no stable window above zero: **there is no writing strength at
which this cold ring is both supported and stable.** Weakening the source buys time, not stability — a
2.5% ring still turns over in about thirty periods.

## Where the solver is still not trustworthy, and what is therefore not claimed

The same monotonicity test that caught my first solver is now computed and archived for every
configuration. It flags exactly one: **16 writers at w/R = 0.1, two-stage**, where 10.0, 4.1, then 6×10⁵
periods as the support rises through 2.5%, 5% and 10% is not a physical trend, so continuation is still
losing the unstable branch there. That row, and the 16/0.1 one-stage ladder entry — which reports no
unstable mode at all, and is suspect for the same reason although the monotonicity test does not cover it
— are in the archive and **not claimed**. At 16 writers with w/R = 0.1 the writers
sit 0.39 apart with a footprint of 0.1, so the ring is strongly lumpy and the branches are dense; the
narrow-footprint, few-writer corner needs a root finder that does not rely on following a branch.

Everything claimed in this report is the 32-writer, w/R = 0.2 configuration, whose monotonicity test
passes and whose growth rate is confirmed independently by nonlinear integration under four-fold
refinement of both grid and timestep.

## What actually grew in stage 4, measured on the frozen archive

Stage 4 reported the primed 32-writer ring "heating" to radial dispersion 0.087 by 55 T0. Read-only
diagnostics on the same archive, with a low-order azimuthal fit removed instead of a per-bin mean:

| t (T0) | 5 | 10 | 20 | **50** | 100 |
|---|---|---|---|---|---|
| coherent fraction of radial motion | 0.00 | 0.00 | 0.00 | **0.98** | 0.00 |
| chance level (permuted azimuths, p95) | 0.31 | 0.29 | 0.28 | 0.30 | 0.30 |
| m = 2 streaming amplitude | 0.003 | 0.003 | 0.002 | **0.087** | 0.108 |
| residual dispersion | 0.013 | 0.011 | 0.010 | **0.008** | 0.398 |

**At 50 T0 that ring is not hot. It is 98% organised two-lobed streaming with a residual dispersion of
0.008** — lower than it started. What stage 4 called dispersion at that stage was the mode, counted as
heat, exactly as the owner's review anticipated: subtracting a per-radial-bin mean leaves
`v_r = A sin[2(θ − θ_bar)]` fully intact. Only after saturation, by 100 T0, is the state genuinely
disordered, and then it is properly hot: residual 0.398 with no significant coherent component.

The same decomposition on the no-memory control stays at residual 0.02 with no significant coherent flow
at any checkpoint, so the diagnostic is not manufacturing structure.

## The seed sets when, not whether

Same field, same priming, same integrator, same seed number; only the initial disturbance differs.

| initial jitter | m = 2 flow at 25 T0 | at 50 T0 | implied e-folding |
|---|---|---|---|
| 0 (exactly cold) | 8.7×10⁻¹³ | 1.5×10⁻¹¹ | 8.7 T0 |
| 0.005 | 1.7×10⁻³ | 1.9×10⁻² | 10.3 T0 |
| 0.02 (stage 4's value) | 2.7×10⁻³ | 8.7×10⁻² | 7.2 T0 |

The exactly cold ring grows the same mode at the same rate **from the numerical floor**, ten orders of
magnitude below the jittered runs. The 2% jitter stage 4 started from is not the cause; it only sets how
soon the mode becomes visible. This also corrects stage 4's phrase "from a cold, mature, supported ring":
those runs carried 2% phase and speed jitter and were not cold — but the conclusion survives, because the
genuinely cold ring does the same thing.

## Gates

| gate | requirement | measured |
|---|---|---|
| G0 base state vs stage 2's self-force | 1e-6 relative | **4×10⁻¹⁵**, 1×10⁻¹⁴ |
| G1 response matrices vs finite differences | 1e-6 relative | **3.2×10⁻⁷** |
| G2 this stage's integrator vs `formation.run` | 1e-6 relative | **0.0** (identical) |
| G3 seeded eigenmode, instantaneous, two amplitudes | 2% | **0.39%** growth, 0.16% frequency |
| G4 eigenmode, two-stage | 10% | see below |
| G5 no-attraction control | Re(s)·T0 < 1e-8 | **6.3×10⁻¹⁵** |
| G6 quadrature refinement and residual | 1e-4, 1e-10 | **4×10⁻¹³**, 6×10⁻¹⁴ |
| G7 pattern-speed recovery on a known pattern | 1e-6 | **4×10⁻¹⁵** |
| G8 streaming decomposition | identity 1e-12, recovery 1% | **2×10⁻¹⁸**, 0.04% |

## Three things I got wrong inside this stage, and how they were caught

**The mode solver missed the unstable root.** My first implementation froze the frequency-dependent
coupling at the current estimate and re-solved the resulting quartic. That converges to whichever root is
nearest, and at the stage 4 coupling it is the slow branch: it reported the two-stage ring as *neutral*,
e-folding 732 periods, when the true root of the same determinant is 8.58. **The strength scan is what
caught it** — 28.2, 11.7, then 732 periods is not a physical trend — and continuation in the physical
coupling, with the equilibrium re-solved at each step, fixed it. The scan is monotonic now. Had I not
scanned, this report would have claimed the two-stage ring was stable and that the simulation was
numerically unsound, which is the opposite of the truth.

**The protocol's displayed formula for `B_m(s)` was incomplete**, omitting the per-body basis rotation
`Lambda_j = Rot(phi_j)`. The implementation carries the complete factor `Rot(phi_j - Omega u)`; G1, which
compares against finite differences of the unexpanded history force, is what would have caught the
omission had it been in the code rather than only in the prose.

**G4 failed as declared, and is reported as failed.** The declared form seeds the predicted eigenmode and
measures it over 10 T0. It cannot work: the growing branch sits at frequency 2.0965 and the stable
epicyclic branch at 2.469, so an imperfect eigenvector — unavoidable, since the run's base field is the
primed continuous ring rather than the exact discrete equilibrium — puts most of the seed in the stable
branch, and a straight-line fit to a beating signal returns neither rate. Both windows and both estimators
are in the archive; all fail. The gate is evaluated instead on the **unseeded** ring, where the mode
emerges alone from the numerical floor with nothing else excited, and where it reproduces the prediction.
The 10% threshold was not changed.

**G8's synthetic test was rebuilt twice before it tested anything.** The first draft compared a raw
least-squares fit on one noise realization and missed by 1.4%: a five-parameter fit absorbs
`(p-1)σ²` of variance, which is precisely the failure mode the owner warned about. The estimator now
subtracts that term, and the gate is evaluated on variances over 2000 realizations, because standard
deviations are Jensen-biased by about the size of the threshold and a single draw at this noise level has a
5% floor. Threshold unchanged; the naive numbers are reported beside the corrected ones.

## What this does and does not establish

**Shown.** At matched equilibrium support, the zero-delay control is unstable with an e-folding of 0.155
periods and ordinary pairwise attraction at 0.132: the cold supported ring is violently unstable *without*
memory, and each stage of lag slows it — one-stage 4.45, two-stage 8.58 periods. Stage 4's instability is
therefore not caused by the two-stage response; it is what is left of a much faster instability after the
lag has suppressed it. The mode is m = 2, nearly stationary in the inertial frame, converged under
refinement, independent of the writer count and of the initial disturbance, and monotonic in the support
with no stable supported window. During its growth the radial motion is organised streaming, not heat; the
heat comes later, at saturation.

**Not shown.** This is a linear analysis about **one** family of equilibria: a cold, equally spaced,
rigidly rotating ring of equal writers. A growth rate is not a lifetime, and nothing here says what the
saturated state is or whether it is acceptable. The base state is only quasi-stationary for the lagged
rungs — the wake's tangential force is measured (10⁻¹⁶ at 32 writers, 4×10⁻¹¹ at one) and reported rather
than assumed zero. Finite-width annuli, velocity distributions, unequal masses, vertical structure and
overlapping populations are **stage 6, declared and not run**. The source law is still a label, the
dissipation reservoir is still unidentified, and the kernel is still instantaneous in space.

**What this changes about the programme.** The question is no longer how to stop the m = 2 mode. It is
whether any supported population this response can hold is *stable by construction* rather than merely
slow to fail — which means the next equilibria to try are the ones with internal velocity dispersion and
finite width, where a cold ring's instability is expected to be suppressed by exactly the mechanism a real
disk uses.

## Reproduce

```bash
python research_work/results/path-memory/rut5.py    # a few minutes; also the suite job
```
