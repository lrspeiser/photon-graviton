# RUT-1 stage 8: the linear modes of the constructed populations, predicted before they were measured

Protocol `protocol-rut8.md` (20a6401), pushed before any stage 8 run. Predictions archive
`rut8-predictions.json`, pushed (aaf0bfe) **before any new population was simulated**. Amendment
`protocol-rut8-amendment-1.md` (cf908ee), written after a declared gate failed and declared before its
replacement was run. Library `warm_modes.py`, units of work `rut8_tasks.py`, driver `rut8.py`, measurement
archive `rut8-results.json`, series `rut8-series/`, suite job `rut8_checks.py`. Stage 7's code builds, checks,
draws and evolves every population here, by import and unchanged.

This is item 6 of the owner's handoff: the linear mode calculation around the populations of stage 7, which
stage 7 said was the only thing that could explain its result, locate the threshold, or exclude slower modes.

## The three statuses, kept separate

| status | result |
|---|---|
| **reproduction** | **passes.** `rut8_checks.py` re-hashes code, both protocol files and every series file in both archives; recomputes from scratch gate L1 on all eight populations and its control, L5's control, the sub-rectangle holding each unstable root, a member of the family scan, a ring-limit annulus with its exact-ring target, and a new population; re-derives both ring-limit gates, the threshold fit and every part M reading from the committed files; and replays a 2 T0 prefix of a live and frozen pair. Zero differences. |
| **numerical verification** | **FAILED, and recorded as failed: the ring-limit gate, twice.** L3 failed as declared on one of its four criteria — a ceiling on the order of convergence, 2.507 against 2.5. Amendment 1, written after that failure and before its replacement was run, declared L3b; **L3b failed too**, on its final bound — 0.555% against 0.5%. The amendment said that if L3b failed nothing further would be declared in its place without the owner, and nothing has been. Both misses have one cause, found afterwards by an exploratory check that is not a gate: the declared rule for the kernel nodes gives a very narrow annulus **three nodes**, and the basis is under-resolved there. It does not touch the populations, whose 9 to 26 nodes are converged to 5×10⁻⁸ (gate L4). The other four gates of part L and all six of part M pass, and every negative control is rejected. `rut8_checks.py` exits non-zero, as the protocol says it must. |
| **scientific outcome** | **reported as it fell, and under the protocol not yet claimable: a failed numerical gate blocks the conclusions that depend on it, and whether the ring-limit failure does is the owner's ruling, not mine.** As it fell: the calculation accounts for stage 7 and locates the threshold; of three blind predictions, two are confirmed by the declared reading and one is unresolved. Every qualitative prediction held: which populations are unstable, that only m = 2 is, that the pattern is nearly stationary in the inertial frame, that growth weakens with temperature, and that the population beyond the threshold is quiet. The growth *rate* of the best-measured new population is 22% below the prediction at five standard errors, and I do not know why. |

## The calculation

A perturbation δC_m(r) e^{imθ + st} of the memory field is a potential perturbation −δC for the bodies. In
action-angle variables of the unperturbed potential Φ₀ = −GM/r − C₀(r), harmonic by harmonic in the radial
angle,

    δf_l(J) = F_l(J) ψ_l(J) / (ν_l(J) − i s),    ν_l = l Ω_r + m Ω_θ,    F_l = ν_l ∂f₀/∂E + m ∂f₀/∂L|_E,

and the memory field answers the density locally in the Laplace variable, δC = α H(s) K ∗ δΣ with
H(s) = τ_keep/[(1 + sτ_keep)(1 + sτ_form)]. Expanding δC_m in the footprint's own azimuthal harmonics
k_m(r_b, ·) = exp[−(r−r_b)²/2w²] I_m^e(r r_b/w²) and collocating on the same nodes gives T(s)c = 0 with

    T(s) = P + α H(s) M(s),    M_ab(s) = (2π)² Σ_l ∫ dJ_r dL  G_al G_bl F_l / (ν_l − i s),

G_al being the orbit transform of k_m(r_a, ·). Overlapping Gaussians are nearly dependent, so T is posed in
the orthonormal basis of P's well-represented subspace — the same roots and a determinant of order one. At
the declared node spacing that is 9 to 26 functions; when the nodes are refined the basis stops growing (13
for `A_cold` at a spacing of 0.05, 14 at 0.025) and the root does not move. Roots in Re s > 0 are located by stage 7's contour method and
counted independently by the argument principle. **A damped mode is invisible to this search**, and the floor
is a growth rate of 0.006 (an e-folding of 26.5 T0), where the orbit quadrature was measured to stop
resolving the resonant denominators: no root above the floor is not stability.

## Part L. The calculation, verified

| gate | result | the control it had to reject |
|---|---|---|
| L1 the orbit library | 12,288 orbits per population: mass and four density moments agree with the solver's density to 1.2×10⁻⁵ on all eight populations (declared 5×10⁻⁵); for the most nearly circular orbit of every L, Ω_r matches κ to 3.3×10⁻⁶ and Ω_θ matches Ω_c to 1.1×10⁻⁷ (declared 10⁻³) | with the memory field removed from the potential the moments are off by 223%: rejected |
| L2 the response matrix, in the time domain | 200,000 bodies under an imposed growing perturbation, +ε and −ε on identical bodies: the kernel-weighted harmonic of the difference equals −M_ab(s)ε to **1.30%** of the largest element for `A_cold` at s = 0.05 − 0.05i and **1.95%** for `A_warm` at 0.03 − 0.10i, real and imaginary parts at every node (declared 3%). In the prototype the same comparison was 2.1% at 100,000 bodies and 0.58% at 400,000: it is sampling noise that remains | keeping only the l = 0 harmonic is off by 99% and 117%: rejected |
| **L3 the ring limit** | **FAILED as declared.** Differences from the exact-history ring 23.35%, 4.67%, 0.63%: falling at every step ✓, last under 2% ✓, control rejected ✓, **order 2.507 against a declared ceiling of 2.5 ✗** | the ring at R = 1 instead of the annulus's radius is off by 72%: rejected |
| **L3b the ring limit, amendment 1** | **FAILED as amended.** Four annuli: 23.35%, 4.67%, 0.63%, **0.555%** — falling at every step ✓, order 1.86 ≥ 1.5 ✓, control rejected ✓, **last under 0.5% ✗** | the ring at R = 1 is off by 70%: rejected |
| L4 convergence | doubling the orbit quadrature, doubling the angle nodes, raising the harmonics to 12, refining the kernel nodes to 0.075 and widening their range move the m = 2 root of `A_cold` and of B16 by at most **4.8×10⁻⁸** (declared 10⁻⁴) | with only l = 0 the root disappears: rejected |
| L5 completeness | 864 sub-rectangles: the winding number of the full determinant equals the number of located roots in every one; singular gap ≥ 3.9×10⁹ wherever a root is located; relative residual ≤ 2×10⁻¹⁶; the two partitions return identical roots | the 2×2 determinant formula applied to the n×n problem returns a winding number of 0 around `B_cold`'s root: rejected |

### The ring limit, and what failed

Three narrow cold annuli at L₀ = 1, each solved for a 9.99% mean support, against `ring_modes.Ring` — the
exact-history theory of stage 5, certified by stage 7's contour solver — at the annulus's **own** mean radius
and total writing rate:

| dL, dE | radial width | annulus, m = 2, as declared | exact ring at the same radius and rate | difference in growth rate (the gates use its size) |
|---|---|---|---|---|
| 0.04, 0.004 | 0.02793 | 0.044634 − 0.012339i | 0.036185 − 0.006049i | +23.35% |
| 0.02, 0.001 | 0.01335 | 0.034482 − 0.005182i | 0.032944 − 0.004512i | +4.67% |
| 0.01, 0.00025 | 0.00660 | 0.032423 − 0.004323i | 0.032222 − 0.004214i | +0.63% |
| 0.005, 0.0000625 (L3b only) | 0.00329 | 0.031868 − 0.004125i | 0.032046 − 0.004144i | **−0.555%** |

**Both gates failed, and the reason is the same.** I found it afterwards, with an exploratory check that is
not a gate and changes no recorded outcome.

The declared kernel nodes run from 0.1 inside the smallest pericentre to 0.1 outside the largest apocentre, 0.1
apart. For the populations that is 9 to 26 nodes. For an annulus 0.007 wide it is **three**, and three Gaussians
0.1 apart do not represent a field centred between two of them. On the third annulus, with everything else as
declared:

| kernel basis | nodes | basis size | m = 2 root | against the ring |
|---|---|---|---|---|
| **as declared: spacing 0.1, margin 0.1** | 3 | 3 | 0.032423 − 0.004323i | +0.626% |
| as declared, but the 64 × 32 orbit quadrature | 3 | 3 | 0.032423 − 0.004323i | +0.626% |
| margin 0.3 | 7 | 7 | 0.032554 − 0.004350i | +1.030% |
| spacing 0.05 | 6 | 6 | 0.032554 − 0.004350i | +1.032% |
| spacing 0.05, margin 0.3 | 14 | 11 | 0.032554 − 0.004350i | +1.032% |
| spacing 0.025, margin 0.5 | 43 | 16 | 0.032554 − 0.004350i | +1.032% |

Every adequate basis gives the same root to six figures, and it is not the declared one. The orbit quadrature,
which I blamed in the amendment, is not the cause: it changes nothing. **So the third difference is 1.03%, not
0.63%**; the artificially small 0.63% is what pushed L3's fitted order over its ceiling, and the fourth
annulus's −0.555% — the only negative difference in the sequence — is the same defect again. On the fourth annulus the same
table reads −0.555% with the declared three nodes, at either orbit quadrature, and **+0.177% to +0.180%** with
5, 7 or 13 nodes — 0.032104 − 0.004172i from every adequate basis.

Two things follow, and I can claim only one of them. What I can claim: the amendment's diagnosis of L3 was
wrong, and said so with more confidence than it had earned. What I cannot: that the ring limit "really passes".
With an adequate basis the prototype's sequence was 23.35%, 4.70%, 1.03%, 0.18%, which would have met L3 and
L3b alike — but that is a prototype and an exploratory table, not a declared gate run on its declared
configuration, and two declared gates failed. A properly specified third version needs a node rule that
cannot leave a narrow annulus with three nodes, and a convergence check made *on the annuli themselves*. It is
the owner's to ask for.

**What the failure does and does not reach.** Gate L4 refined the kernel nodes and widened their range on
`A_cold` and on B16 and moved their roots by 5×10⁻⁸; every other population has as many nodes or more. So the
populations' predictions do not depend on the defect, as far as L4 tested.
But L4 checked convergence on the populations and not on the annuli of L3, which is where it was needed, and
that is my omission.

**It also answers a question stage 7 left open.** Stage 7 found the cold annulus growing twice as fast as the
cold ring's 8.58 T0 and I wrote that I did not know why. The ring theory evaluated where the annulus actually
sits — R ≈ 0.91, where a unit of angular momentum puts it under 10% extra support, with the annulus's total
writing rate — gives 0.032 to 0.036, an e-folding of 4.4 to 5.0 T0, not 8.58. Most of the difference was the
radius; finite width and temperature account for the rest.

## Part P. The predictions

**In the declared strip — growth rate at least 0.006, pattern speeds from −1/m to 2 + 1/2m, m = 1 to 4 — every
population at or below dE = 0.018 has exactly one unstable root, at m = 2, and no other population has any.**

| population | dE | σ_r | m = 2 root | e-folding | pattern speed | m = 1, 3, 4 |
|---|---|---|---|---|---|---|
| A_cold | 0.010 | 0.0940 | 0.050876 − 0.052004i | 3.13 T0 | +0.0260 | none |
| B_cold | 0.010 | 0.0940 | 0.041315 − 0.077451i | 3.85 T0 | +0.0387 | none |
| **B13** (new) | 0.013 | 0.1077 | 0.032160 − 0.087854i | 4.95 T0 | +0.0439 | none |
| **B16** (new) | 0.016 | 0.1201 | 0.020226 − 0.097944i | 7.87 T0 | +0.0490 | none |
| B_mid | 0.020 | 0.1348 | none above the floor | — | — | none |
| **B24** (new) | 0.024 | 0.1481 | none above the floor | — | — | none |
| A_warm, B_warm | 0.030 | 0.166 | none above the floor | — | — | none |

The pattern speeds are 2% to 4% of the orbital frequency: **the unstable mode of a warm annulus is the same
inertially stationary two-lobed pattern as the cold ring's**, and the calculation says so from first
principles.

**The threshold.** Through family B — mass, α, footprint and both response times fixed — the m = 2 rate falls
smoothly with temperature:

| dE | 0.010 | 0.011 | 0.012 | 0.013 | 0.014 | 0.015 | 0.016 | 0.017 | 0.018 | 0.019 – 0.024 |
|---|---|---|---|---|---|---|---|---|---|---|
| growth rate | 0.0413 | 0.0386 | 0.0355 | 0.0322 | 0.0285 | 0.0245 | 0.0202 | 0.0156 | 0.0106 | none above 0.006 |

A quadratic through the five roots below 0.03 reaches the axis at **dE\* = 0.01993** (0.01992 to 0.01994 under
dropping each point in turn), σ_r ≈ 0.134. **Stage 7's `B_mid`, at dE = 0.020, sits at marginal stability.**
That is consistent with what stage 7 saw and could not explain — its radial spread rose 10.9% at 64 bodies,
5.6% at 256 and 0.9% at 1,024. A population at the edge of instability answers shot noise with a large,
slowly decaying m = 2 disturbance, and the heating that causes falls with N without being two-body
relaxation: the case stage 7's discreteness reading could not tell apart. It is an interpretation that the
calculation makes available, not something this stage tested.

**Retrospective, and labelled so.** Against stage 7's exploratory growth-phase fits, which I had seen before
computing these: `A_cold` predicted 0.0509 and +0.0260, measured 0.0480 ± 0.0071 and +0.0258 ± 0.0009 at 1,024
bodies; `B_cold` predicted 0.0413 and +0.0387, measured 0.0370 ± 0.0055 and +0.0448 ± 0.0056.

## Part M. The three populations that had never been simulated

Built with stage 7's construction at the family's α; V1, V2, V3, V5 and V6 (four realizations of the declared
draw and of each of its three negative controls, 48 runs of 16,384 bodies) all pass and every control is
rejected; 72 runs at 256 and 1,024 bodies, live and frozen from byte-identical bodies, all completed. The
horizon of each was sized from its prediction, and the growth-phase rule — stage 7's, fixed before these
populations existed — was applied as declared.

| population | N | fit window | measured rate | predicted | measured pattern speed | predicted | reading |
|---|---|---|---|---|---|---|---|
| B13 | 256 | 6 – 9.75 T0 | 0.0365 ± 0.0068 | 0.0322 | +0.0451 ± 0.0023 | +0.0439 | |
| B13 | 1,024 | 6 – 14.75 T0 | **0.0252 ± 0.0013** | **0.0322** | +0.0445 ± 0.0019 | +0.0439 | **unresolved** |
| B16 | 256 | 6 – 8.25 T0 | 0.0346 ± 0.0171 | 0.0202 | +0.0458 ± 0.0109 | +0.0490 | |
| B16 | 1,024 | 6 – 11.75 T0 | 0.0270 ± 0.0071 | 0.0202 | +0.0552 ± 0.0111 | +0.0490 | **confirmed** |
| B24 | 256, 1,024 | no growth phase | — | none | — | — | **confirmed** |

**What held.** Both populations predicted unstable grow an m = 2 pattern at both body counts, and the one
predicted to hold no root does not: its |C₂| peaks during the shot-noise imprint and decays. The pattern
speeds are right — 0.0445 ± 0.0019 against 0.0439 for B13. The ordering is right, in the rates and in what
the instability does to the population (D at 1,024 bodies, percent):

| quantity | B13 (dE 0.013) | B16 (0.016) | B24 (0.024) | for comparison, stage 7: B_cold (0.010), B_mid (0.020) |
|---|---|---|---|---|
| radial spread | +43.7 ± 0.5 | +22.5 ± 2.2 | +0.17 ± 0.49 | +73.2 ± 3.4, +0.87 ± 0.21 |
| residual radial dispersion | +47.8 ± 0.8 | +21.6 ± 2.1 | +0.70 ± 0.42 | +72.4 ± 4.4, +0.19 ± 0.10 |
| mean angular momentum | −1.31 ± 0.03 | −0.73 ± 0.08 | −0.02 ± 0.00 | −2.01 ± 0.05, −0.021 ± 0.005 |
| mean radius | −0.43 ± 0.05 | −0.39 ± 0.08 | +0.01 ± 0.07 | −0.63 ± 0.04, +0.06 ± 0.03 |

Five temperatures, physics fixed: the damage falls monotonically with temperature and vanishes between
dE = 0.016 and 0.020, where the calculation puts the threshold.

**What did not.** B13 at 1,024 bodies is the best-measured cell in the stage — four realizations at 0.0242,
0.0270, 0.0220 and 0.0274, a clean exponential whose local slope is 0.025 from 8 to 18 T0 — and it is **22%
below the prediction, at 5.4 standard errors**. By the declared reading that is neither confirmed (it is not
within three standard errors or 15%) nor contradicted (it is not beyond 30%): **unresolved**. B16's "confirmed"
is a weak one and should be read as such: its rate is 33% *above* the prediction and passes only because its
standard error is large, the four realizations running from 0.010 to 0.045.

**I do not know why B13 grows more slowly than predicted, and I have not tested any explanation.** The
candidates I can name: the simulated background is not the smooth equilibrium the theory perturbs — shot noise
heats it before and while the mode grows, and the rate falls by 0.0035 for every 0.001 in dE, so a 15% rise
in dE would do it; the field's axisymmetric part is live in the simulation and frozen in the theory; and the
saturation amplitude is only about eight times the shot-noise level at 1,024 bodies, so the fit window holds
less than two e-foldings and is never far from either noise or saturation. The stage 7 populations show the
same thing in both directions: over 8 to 16 T0 at 1,024 bodies `A_cold`'s local slope is 0.0455, 11% below
its prediction, and `B_cold`'s is 0.0473, 14% above, wandering from 0.003 to 0.057 on the way. None of this
is a finding.

**What this says about the measurement, not the theory.** Shot noise is a poor way to start a slow mode. With
the seed at 2×10⁻⁴ and saturation near 10⁻³, a mode with an e-folding of 8 T0 is visible as an exponential for
about one e-folding. The clean test is the one stage 6 built for the ring and the owner's item 6 names:
**seed the predicted eigenvector in the full state — bodies, excitation and field — at an amplitude well above
the noise and well below saturation**, where stage 6 measured the ring's rate to 0.014%. `warm_modes.py`
already produces the eigenvector.

## What I got wrong on the way

* **Both ring-limit gates, above.** L3's ceiling was a check sized on a configuration other than the one it was
  run on, for the third time in two stages. Then I diagnosed its failure wrongly, in an amendment, before
  testing the diagnosis — I blamed the orbit quadrature, which changes nothing, when the kernel basis had three
  nodes — and declared a replacement that failed for the reason I had not found. The convergence gate, L4, was
  run on the populations and not on the narrow annuli. A number in a protocol has to be checked on the exact
  configuration it will be run on, and a diagnosis has to be tested before it is published.
* **`spectrum.located()` takes the contour determinant with the 2×2 formula whatever the size of T.** I wrote
  it in stage 7, described it as general, and never gave it anything but 2×2 matrices. On this stage's n×n
  problem it returned a winding number of zero around a root. Stage 7 is unaffected and `spectrum.py` is
  pinned by its archive, so it is not edited; this stage carries its own routine, and the old formula is L5's
  negative control.
* **The first kernel basis was nearly singular.** Posed directly on overlapping Gaussians, refining the kernel
  nodes from 0.1 to 0.075 raised the condition number from 2×10⁵ to 5×10⁸ and lost the root. In the orthonormal
  basis the root is identical to seven figures from a spacing of 0.1 down to 0.025.
* **I launched the campaign by accident.** A toy-protocol smoke run in the scratchpad picked up the real
  protocol, which by then was installed in the repository and is found first. It ran the declared computation
  with byte-identical code, after the protocol was pushed; the official run from the repository reused its
  finished tasks — the cache is keyed on the code hashes, which matched — and reproduced it identically. No
  result depends on the accident, but it was one.

## What this stage does not establish

* **Below the floor.** No root above a growth rate of 0.006 is not stability. The cold ring's spectrum has
  unstable roots with growth times of hundreds to 10⁵ periods at every m, and this discretization cannot see
  their counterparts in an annulus. A damped mode is invisible to any search of the right half-plane, so
  `B_mid` is "marginal" only from the unstable side.
* **The growth rate to better than about 20%.** The pattern speeds, the selection of m = 2, the threshold's
  location and the ordering with temperature are confirmed by simulation; the rate is confirmed at that level
  and not beyond it.
* **What sets the threshold.** One footprint, one pair of response times, one angular-momentum width, one
  radius, one coupling. dE\* = 0.0199 is a number for that family, not a criterion.
* Writing proportional to mass is still a label with a declared mass; the dissipation reservoir is
  unidentified; the kernel is instantaneous in space; and the bodies' lost angular momentum still has no
  ledger. No observational comparison is made.

## Next, in the owner's order

The seeded, full-state eigenmode test of these predictions — which would settle B13's 22% — and with it the
same calculation pushed below the present floor, which needs a better treatment of the resonant denominators
than a tensor Gauss-Legendre grid; then formation from an empty field, and the reciprocal energy and
angular-momentum accounting through growth and saturation. Each in its own protocol file, with every number
in it checked on the configuration it will be run on.
