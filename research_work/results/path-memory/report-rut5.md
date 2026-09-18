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
converged under four-fold refinement of both the grid and the timestep. It is not a numerical artefact,
and its rate does not change when the same total source is divided among 32, 64 or 128 writers.

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
times more unstable** in this configuration, and ordinary pairwise attraction at matched support is
unstable too, at 0.132 T0. **So the lag is not what destabilises the ring. Collective attraction
destabilises it, and the lag is the only reason it survives long enough to be interesting.** Each stage of
lag buys about a factor of two in rate here: zero → one → two stages runs 0.155 → 4.45 → 8.58 T0. These
are ratios for this configuration, not universal factors.

**What the Newtonian rung does and does not say.** It establishes that a cold ring with ordinary
attraction of the same strength is also violently unstable, which is enough to reject the claim that
memory is *necessary* for the instability. It is not an identification of the two: its fastest mode is
m = 9 against the Gaussian control's m = 4, so the branches differ in structure even where the timescales
are close. Similar timescales are not the same instability.

The mechanism is visible in the matrices: a body's displacement is written along its whole orbit, so a
lagged field responds to a Bloch-m displacement with the m ≠ 0 structure largely averaged away, while an
instantaneous field feels the full discrete structure. That is the same collective cancellation stage 2C
measured for the drag, acting here on the stability.

**The two-stage mode's orientation is nearly stationary in the inertial frame** — pattern speed
0.0005 Ω, against the one-stage mode's 0.993 Ω, which corotates. A memory field is written in space and
decays in place, so the pattern it crystallises does not have to travel with the matter.

**That is a nearly stationary orientation with a growing amplitude, not a stationary field, and the
distinction matters.** At a fixed point in space the m = 2 harmonic varies as

    s_spatial = s − i m Ω = 0.018552 − 0.001027 i,

so the orientation turns slowly but the amplitude changes substantially on the same timescale. It would be
wrong to conclude that the mode therefore sees the zero-frequency response: its growth rate enters the
response as much as its frequency does, which is exactly why the mode condition has to be solved with the
full complex H(s) rather than evaluated at a guessed real frequency. This answers the measurement stage 4
said was missing, and it remains a prediction for the next campaign to test directly on the
force-producing field, whose own phase history has still never been recorded.

## The large-scale rate converges as the writers are subdivided

At fixed **total** writing rate, w/R = 0.2, fastest mode over m ≤ 8:

| writers | 8 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|
| two-stage e-folding (T0) | 1.22 | 0.61 | **8.58** | **8.58** | **8.58** |
| instantaneous e-folding (T0) | 0.331 | 0.107 | 0.155 | 0.155 | 0.155 |

Eight or sixteen writers is a lumpy ring — the spacing 2π/8 = 0.79 is four times the footprint width — and
it is much more unstable, as a ring of discrete blobs should be. By 32 writers the spacing 0.196 has
dropped below the width 0.2 and the rate is **identical at 32, 64 and 128**: this mode is a property of the
supported ring rather than of how finely the same source is divided, and the stage 4 configuration sits in
that converged limit.

That is a statement about *this* large-scale mode in *this* family. It does not establish that particle
discreteness is irrelevant to every mode or to the nonlinear outcome — the search here covers m ≤ N/2 in
the ladder and m ≤ 8 in the scans, and a converged linear rate says nothing about how a discrete source
behaves once the mode saturates.

## Whether a weaker ring is stable

Thirty-two writers, w/R = 0.2, fastest mode over m ≤ 8:

| label fraction | 0 | 0.025 | 0.05 | 0.10 | 0.20 |
|---|---|---|---|---|---|
| support | 0 | 2.5% | 5.0% | 9.99% | 20.0% |
| two-stage e-folding (T0) | ∞ | 28.2 | 11.7 | **8.58** | 2.45 |
| instantaneous e-folding (T0) | ∞ | 0.215 | 0.180 | 0.155 | 0.076 |

Growth is monotonic across the four sampled positive strengths, with no stable case among them:
**none of the sampled supported strengths gives a stable cold ring.** Weakening the source buys time, not
stability — a 2.5% ring still turns over in about thirty periods. Four discrete couplings, with modes
searched through m = 8, are not a proof that no stable interval exists between or beyond them; what they
show is that the trend runs the wrong way for finding one.

## The root count certifies nothing, and is itself wrong

The spectrum records a contour count beside the roots found by continuation. **It does not establish
completeness, and in its shipped form it is unreliable.** The owner's review caught this, and the archive
already contained the evidence:

* For the **no-attraction control**, whose determinant is exactly `s²(s² + Ω²)` with no root in the open
  right half-plane, the counter returns **one unstable root for every angular mode**. The double root at
  the origin sits on the contour's left edge, 10⁻⁴ Ω away, and the phase winding there is under-resolved.
  The false count is a systematic +1.
* The same +1 appears throughout the two-stage and one-stage spectra, so `all_counts_match` is false for
  configurations whose dominant root is perfectly well determined — including the m = 2 root this report
  rests on, whose residual is 6×10⁻¹⁴.
* One two-stage mode carries a Newton residual of **4.8×10⁻⁶**, far outside the 10⁻¹⁰ that G6 verified for
  m = 1 to 4. That mode's root is not converged.

So a count mismatch currently cannot distinguish "continuation missed a root" from "the counter
miscounted". **Nothing in this report claims a complete spectrum or the absence of unstable modes**; every
claim rests on roots that were located, with their residuals, and on the nonlinear agreement. Establishing
that one mode grows is a much weaker requirement than establishing that none does, and the second is what
stage 6 will need. The repair — a contour-integral eigensolver that locates all roots in a declared
region, validated first against spectra that are known exactly, with resolution and contour varied — is
declared in stage 6 and is a precondition for any claim that a configuration is stable.

Monotonicity is also only a warning. It caught my first solver, and it still flags the 16-writer case
below, but no monotonicity theorem for this coupled nonlinear eigenproblem is established anywhere here,
so passing it cannot certify that every relevant root was found.

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
| G4 eigenmode, two-stage | 10% | **failed as declared**; see below |

**G4's status, stated as three separate facts:**

| | outcome |
|---|---|
| the originally declared seeded test | **failed** (both windows, both estimators, archived) |
| the unseeded verification that replaced it | **passed**: 0.26% in rate, 0.002% in frequency |
| a full-state seeded test, with the field perturbations initialised | **not yet done** |

Keeping the tolerance unchanged does not make a changed initial condition and a changed estimator into the
experiment that was declared, and the headline should not blur them.
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

**G4 failed as declared, and is reported as failed — and my diagnosis of *why* was incomplete.** I
attributed it to the nearby stable branch at frequency 2.469 and to the approximate base field, and both
matter. But the owner's review identified a more direct cause: **`seed_mode` perturbs positions and
velocities only.** The state of this model is four things — positions, velocities, the excitation field E
and the force-producing field C — and a perturbation of the first two is not an eigenmode of all four, so
it necessarily launches transients on other branches. Calling the seeded test intrinsically
ill-conditioned was premature; it has a setup defect that can be fixed, by building the mode's own source
history into E and C before release. That is declared in stage 6 and not yet done. Both windows and both estimators
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
lag has suppressed it. The identified mode is m = 2, its orientation nearly stationary in the inertial
frame while its amplitude grows, converged under refinement, unchanged from 32 to 128 writers, and
independent of the initial disturbance, which sets only when it becomes visible. Growth rises with support
across every sampled strength. During its growth the radial motion is organised streaming, not heat; the
heat comes later, at saturation.

**Not shown.** **No spectrum here is certified complete**, for the reasons above: one mode was located and
verified, not all modes excluded. This is a linear analysis about **one** family of equilibria: a cold,
equally spaced, rigidly rotating ring of equal writers. A growth rate is not a lifetime, and nothing here
says what the saturated state is or whether it is acceptable. The base state is only quasi-stationary for the lagged
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
