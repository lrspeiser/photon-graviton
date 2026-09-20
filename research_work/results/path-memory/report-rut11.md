# RUT-1 stage 11: the quiet region

Protocol `protocol-rut11.md`, pushed with the code before any stage 11 run. Units of work `rut11_tasks.py` and
`resonant_response.py`, driver `rut11.py`, archive `rut11-results.json`, series `rut11-series/`, suite job
`rut11_checks.py`. Stage 9's response, stage 10's quiet sample, stage 7's simulator and population builder and
the owner's sampler are used by import, unchanged.

This is the third of the four pieces of work the owner set in the review of 3a80fec — *establish the quiet
region's scope* — and the last clause of the milestone.

## The three statuses, kept separate

| status | result |
|---|---|
| **reproduction** | see `rut11_checks.py` in the suite: the rule recomputed from scratch, one disturbed run prepared and replayed bit for bit, and every gate, reading and fit re-derived |
| **numerical verification** | **FAILED**, 8 gates of 10 |
| **scientific outcome** | archived as it fell; it never sets the exit status |

## The two gates that failed, and why

Both are faults in what I declared, not in what the rule does, and by the owner's rule on 3a80fec neither is
reworded: they stay failed here.

**A1, the rule on a known integral.** The protocol declares the pair (192×64, 384×128) as the rule and sets its
tolerance at 10⁻⁶, quoting prototype errors of "9×10⁻¹⁰ to 5×10⁻⁸". Those numbers are the *finer* pair,
(384×128, 768×256): I read the wrong line of my own prototype and attributed one configuration's accuracy to
another. The declared pair reaches 1.64e-06 at the lowest growth rate, over its threshold. The rule is
unaffected — the finer pair is at 3.48e-08, the negative control is rejected at 6, and every
contour in parts L and T carried its own error monitor, which passed — but the threshold was set from a
measurement of something else, which is exactly the mistake the prototype rule exists to prevent.

**L1, the search is sound.** Of 540 rectangles,
2 counted a root by the argument principle and declined to locate
it — `box:B_dE0.019:m1:s5` at a Beyn rank gap of 937, and `box:B_mid:m1:s5` at a Beyn rank gap of 973 against a threshold of 1000. That threshold is `singular_gap`, which this stage inherited from stage 9's block
instead of declaring in its own: a second declaration fault. In each case the other partition located the same
root cleanly, B_dE0.019 at a gap of 11079 and B_mid at a gap of 3902, so no reading below depends on a rectangle that failed.

## Part A. The rule

Stage 8's floor of 0.006 is not a tolerance anyone chose. Its quadrature replaces the continuum of orbital
frequencies by 12,288 of them, so the discretised response has a pole at s = −iν for each, all on the imaginary
axis about 1.2×10⁻⁴ apart. This stage integrates each cell of a refined grid in closed form instead, so the pole
is integrated rather than sampled.

Against a reference integral of the same shape as the response's own, good to fourteen digits:

| growth rate Re s | the rule (the declared pair) | the finer pair | stage 9's quadrature |
|---|---|---|---|
| 0.03 | 1.46e-08 | 9.11e-10 | 1.15e-10 |
| 0.006 | 5.19e-08 | 5.85e-10 | 7.04e-04 |
| 0.001 | 6.93e-07 | 4.87e-08 | 0.061 |
| 1.0e-04 | 1.47e-06 | 3.85e-08 | 1.9 |
| 1.0e-06 | 1.64e-06 | 3.48e-08 | 6 |

Where stage 9's quadrature *is* valid, the two rules must agree, and they do:

| population | stage 9's root | the rule's root | relative difference | under the finer pair | declared error |
|---|---|---|---|---|---|
| A_cold | 0.0508781-0.0520359i | 0.0508781-0.0520359i | 1.14e-08 | 1.07e-08 | 1.76e-04 |
| B13 | 0.0321618-0.0878664i | 0.0321618-0.0878664i | 2.66e-08 | 2.50e-08 | 1.91e-04 |
| B16 | 0.0202286-0.0979527i | 0.0202286-0.0979527i | 3.79e-08 | 3.69e-08 | 1.67e-04 |
| B_cold | 0.041317-0.0774698i | 0.041317-0.0774698i | 1.64e-08 | 1.54e-08 | 1.84e-04 |
| B_dE0.017 | 0.0155846-0.101321i | 0.0155846-0.101321i | 4.50e-08 | 4.20e-08 | 1.83e-04 |
| B_dE0.018 | 0.0105855-0.104756i | 0.0105855-0.104756i | 2.38e-07 | 4.84e-08 | 1.94e-04 |

Used in the low box, where its poles lie, stage 9's own rule breaks: on the declared control rectangle it
returns a winding number of 3 while locating 1,
and its root differs from the analytic rule's by 0.011.

The error the rule declares about itself, over every contour in parts L and T: worst
**0.0013** against a budget of 0.002, on 548 contours.

## Part L. The low box

From Re s = 5.0e-04 — an e-folding time of
318 reference periods — to stage 8's
floor of 0.006, over the strip stage 8 declared, m = 1, 2, 3, 4,
in two partitions offset by half a rectangle.

| population | rectangles | m | roots located | worst declared error | nearest pole |
|---|---|---|---|---|---|
| A_warm | 108 | 1, 2, 3, 4 | 2 | 0.0013 | 0.3 |
| B24 | 108 | 1, 2, 3, 4 | 2 | 6.76e-04 | 0.31 |
| B_dE0.019 | 108 | 1, 2, 3, 4 | 3 | 4.48e-04 | 0.3 |
| B_mid | 108 | 1, 2, 3, 4 | 1 | 5.28e-04 | 0.31 |
| B_warm | 108 | 1, 2, 3, 4 | 2 | 9.65e-04 | 0.31 |

The search can see a root that is there: in **B_dE0.019**, which stage 9 reported no root in, having
searched only above 0.006, both partitions locate **0.00522803-0.108323i** and
**0.00522803-0.108323i**, agreeing to 9.65e-16 — an e-folding time of
30.4 periods.

## Part M. The mode domain

Where the spectral norm of αH(s)M(s) stays below one, T = 1 + αHM cannot be singular and no mode exists at that
m **at any resolution**. The maximum over the declared grid of the whole strip, from the new floor to
Re s = 0.5:

| population | m = 1 | m = 2 | m = 3 | m = 4 | m = 5 | m = 6 | m = 8 | m = 10 | m = 12 | m = 16 |
|---|---|---|---|---|---|---|---|---|---|---|
| A_warm | 1.04 | 1.63 | 0.324 | 0.204 | 0.123 | 0.0708 | 0.0355 | 0.0163 | 0.00669 | 8.130e-04 |
| B13 | 1.02 | 7.01 | 0.41 | 0.307 | 0.114 | 0.078 | 0.0368 | 0.016 | 0.00607 | 5.827e-04 |
| B24 | 1.03 | 3.07 | 0.308 | 0.267 | 0.108 | 0.0774 | 0.0373 | 0.0166 | 0.00657 | 7.139e-04 |
| B_mid | 1.03 | 3.96 | 0.344 | 0.281 | 0.111 | 0.0778 | 0.0371 | 0.0164 | 0.00638 | 6.614e-04 |
| B_warm | 1.04 | 2.2 | 0.269 | 0.244 | 0.108 | 0.0765 | 0.0376 | 0.017 | 0.0069 | 8.054e-04 |

Refining the grid moved the largest of these by 0. The control —
B13 at m = 2, which holds a root — reaches
7.01, so the diagnostic discriminates.

## Part T. The family's boundary, bracketed

| dE | root above the floor | root | e-folding, periods |
|---|---|---|---|
| 0.019 | yes | 0.00522803-0.108323i | 30.4 |
| 0.0193 | yes | 0.00355265-0.10943i | 44.8 |
| 0.0195 | yes | 0.00241915-0.11018i | 65.8 |
| 0.0197 | yes | 0.00127371-0.11094i | 125 |
| 0.0198 | yes | 0.000696988-0.111327i | 228 |
| 0.0199 | no | — | — |
| 0.0201 | no | — | — |
| 0.0203 | no | — | — |

The warmest member with a root this rule resolves is **dE = 0.0198**, growing at 6.970e-04 — an e-folding time of 228 periods. The coolest with none above the floor is **dE = 0.0199**. The boundary of the family is therefore **bracketed between 0.0198 and 0.0199**, by resolution rather than by extrapolation. Nothing is said about growth slower than the floor, and a damped mode is invisible to this search: the cooler member is not called stable.

Stage 8's extrapolated crossings, reported beside it and not replaced by it: linear
0.02044, quadratic 0.01993,
cubic 0.01991.

## Part D. The controlled disturbance

A quiet sample, folded so it carries no source at the harmonic being disturbed, and into its field a declared
kernel pattern of relative amplitude 1.0e-05 at m = 2. The bodies are
not moved: there is no mode to prepare them in.

| population | realizations | fitted rate over the window | amplitude at the horizon / at the start | peak / start |
|---|---|---|---|---|
| A_warm | 4 | -0.0179 ± 8.8e-05 | 0.018 | 1 |
| B24 | 5 | -0.0084 ± 8.1e-04 | 0.096 | 1 |
| B_dE0.019 | 0 | — | — | — |
| B_mid | 4 | -4.564e-05 ± 0.005 | 0.64 | 1.2 |
| B_warm | 4 | -0.0149 ± 8.8e-04 | 0.036 | 1 |

The apparatus would see growth: B13, kicked at
1.0e-07 and read over its growth phase, grows at
**0.0344 ± 0.001** against stage 9's recalculated root of
0.032162 — 0.069 relative, over 4
realizations.

### What the wrong preparations and the refinements do

| variant | amplitude | fitted rate | difference from the base run | amplitude at the horizon |
|---|---|---|---|---|
| coarse_grid | 1.0e-05 | -0.0105 | -0.0035 | 0.07 |
| fine_grid | 1.0e-05 | -0.00695 | 4.54e-06 | 0.11 |
| half_step | 1.0e-05 | -0.00695 | 1.01e-06 | 0.11 |
| no_excitation | 1.0e-05 | -0.00448 | 0.0025 | 0.1 |
| ordinary_draw | 1.0e-05 | 0.00229 | 0.0092 | 12 |

The base run it is compared with: rate -0.00695, amplitude at the horizon 0.11 of the start.

## The declared reading

**Every population searched holds a mode at m = 1 that stage 8 could not have seen.** Its growth rate runs from 8.33e-04 to 9.80e-04 — e-folding times of 162 down to 191 reference periods — against a floor of 5.0e-04, seven times lower than the one that hid it, and its pattern speed is within 10⁻⁵ of zero: a lopsided figure that barely turns. The rate **rises** with the population’s internal velocity width, the opposite of the m = 2 branch, so warming the population does not remove it. At m = 2 nothing is found above the floor in any population warmer than the bracket below, and at m = 3, 4, 5, 6, 8, 10, 12, 16 the operator norm forbids a mode at any resolution. So the quiet region of this family is quiet at m ≥ 2 and is **not** quiet at m = 1, and that is the answer to whether the warmer population is genuinely stable, weakly unstable or simply untested: weakly unstable, on a timescale no previous search could reach.

- **A_warm holds a root above the new floor**: 0.00098013-1.17781e-05i. The operator norm bounds m = 3, 4, 5, 6, 8, 10, 12, 16 below 0.32, where no mode can exist at any resolution. A declared disturbance of 1.0e-05 at m = 2, over 4 realizations, decays: the fitted rate over periods 10 to 40 is -0.0179 ± 8.8e-05 and the amplitude at the horizon is 0.018 of its value at the start.
- **B24 holds a root above the new floor**: 0.000880706-7.54749e-06i. The operator norm bounds m = 3, 4, 5, 6, 8, 10, 12, 16 below 0.31, where no mode can exist at any resolution. A declared disturbance of 1.0e-05 at m = 2, over 5 realizations, decays: the fitted rate over periods 10 to 40 is -0.0084 ± 8.1e-04 and the amplitude at the horizon is 0.096 of its value at the start.
- **B_dE0.019 holds a root above the new floor**: 0.00083263-6.68133e-06i, 0.00522803-0.108323i.
- **B_mid holds a root above the new floor**: 0.000841727-6.8099e-06i. The operator norm bounds m = 3, 4, 5, 6, 8, 10, 12, 16 below 0.34, where no mode can exist at any resolution. A declared disturbance of 1.0e-05 at m = 2, over 4 realizations, neither grows nor decays measurably: the fitted rate over periods 10 to 40 is -4.564e-05 ± 0.005 and the amplitude at the horizon is 0.64 of its value at the start.
- **B_warm holds a root above the new floor**: 0.000947606-9.32826e-06i. The operator norm bounds m = 3, 4, 5, 6, 8, 10, 12, 16 below 0.27, where no mode can exist at any resolution. A declared disturbance of 1.0e-05 at m = 2, over 4 realizations, decays: the fitted rate over periods 10 to 40 is -0.0149 ± 8.8e-04 and the amplitude at the horizon is 0.036 of its value at the start.

## What I got wrong on the way

Two faults surfaced in a dry run started after the protocol was pushed, before any result existed, and were
corrected in 90bdd4d with a paragraph in the protocol saying so. Stage 9's `family_name` rounds a width to
three decimals, so six of this stage's eight ladder rungs would have been built as the wrong population,
silently; and the growth control, given the same amplitude as the quiet populations, would have been fitted
through its own saturation. The dry run had finished one task, the reference-integral check, and its output was
discarded.


## A post-hoc check on the m = 1 roots

**Computed after the campaign, after seeing its result, and not a gate.** Part L's rectangles do not carry the
refinement check that part T's ladder does, so the m = 1 roots came out of the search without one. Because
they are the stage's central finding, they were put through the finer pair and the three representation
variants afterwards. Saying so is the point: this was chosen once the answer was known, it cannot change the
stage's status, and it is reported as relative movement of the root from the declared representation.

| population | m = 1 growth rate | the finer pair (2, 4) | kernel spacing w/8, margin 2.5w, shifted | orbits 288 × 96, 256 angles | harmonics to l = 12 |
|---|---|---|---|---|---|
| A_warm | 0.00098012956 | 2.1e-06 | 1.1e-13 | 2.6e-06 | 1.0e-05 |
| B24 | 0.000880705536 | 8.9e-07 | 5.1e-13 | 6.1e-07 | 1.2e-06 |
| B_dE0.019 | 0.000832630294 | 1.6e-06 | 2.5e-12 | 9.9e-07 | 1.1e-07 |
| B_mid | 0.000841727498 | 1.3e-06 | 7.2e-14 | 9.9e-07 | 2.1e-07 |
| B_warm | 0.000947605647 | 7.0e-07 | 1.1e-13 | 2.0e-07 | 6.7e-06 |

The largest movement of any root under any of these is 1.0e-05. That is evidence the roots are real and resolved,
gathered after the fact; it is not a declared verification, and a stage that wants one should declare it.

## What this stage does not establish

Damped modes: the rule is used only for Re s > 0, and a decaying mode needs a continuation across the cut on
Re s = 0 that is not attempted. Anything outside the declared strip in Im s. Azimuthal numbers above the
declared set, except through the operator norm, which is only as good as its grid. Growth slower than
5.0e-04, which the band-edge prototype shows is where this rule's convergence begins to
degrade. The sharply truncated distribution function is still the model. And the disturbance is one declared
shape at one declared amplitude: another shape may excite something this one does not.

## Next, in the owner's order

The fourth piece of work: formation from an empty field with the reciprocal energy and angular-momentum budget
carried through the transition. Before that, two things this stage points at: the analytic rule's closed forms
are explicit analytic functions of s, so continuing them across the cut on Re s = 0 would reach the damped
modes and turn the decay this stage measures in the live system into a prediction; and the operator bound is
evaluated on a grid, where a Lipschitz argument would make it a statement rather than a scan. Then the owner's
prior-art programme in `research_plan/prior-art/next-comparisons.md`, beginning with the published
response-matrix and contour-solver benchmarks.
