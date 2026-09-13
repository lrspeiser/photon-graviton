# Depth-dependent capture: results

The original exact-third reference remains selected. Across 149 galaxies and 3150 radii, all three families choose zero depth dependence under the declared training logarithmic loss. This rejects these bounded replacements under that criterion, not every possible depth law. All data were previously exposed; the old splits are transfer checks, not new blind evidence.

## Hypothesis and provenance

Let B=-Phi be positive well depth, zero at infinity, and S=B/(B+(150 km/s)^2). Test opacity multipliers 1+beta*S (boost), 1-beta+beta*S (shallow suppression), and 1+beta*(2*S-1) (both). Beta is shared and scanned over 0, 0.25, 0.5, 1. The saturating function, Newtonian potential and attenuation integrals are known mathematics. Their use as a companion-capture law is a proposed phenomenological modification, not an established interaction or proven novelty. The original one-third retention exponent and source normalization stay fixed.

The same multiplier changes both absorption at a location and interception along incoming rays. The deposited potential feeds back into capture until a stationary fixed point converges. Ordinary-matter midplane depth is integrated from the saved rotation baseline and extended spherically, with a finite-mass outer tail. This is not a full disk/bulge potential, physical formation history, energy-supply measurement or stability proof. Stronger capture may increase total stored inventory; it is not mass-preserving redistribution.

## Same-resolution comparison

RMSE is the square root of the mean of per-galaxy mean squared speed errors; log RMS weights proportional errors. Selection uses training log RMS. Every trial is reported.

| Rule | beta | Train RMSE | Validation RMSE | Test RMSE | Train log RMS |
|---|---:|---:|---:|---:|---:|
| Reference | 0 | 29.013 | 32.480 | 23.575 | 0.138475 |
| boost | 0.25 | 30.516 | 34.409 | 24.793 | 0.140347 |
| boost | 0.5 | 31.884 | 36.113 | 25.888 | 0.142333 |
| boost | 1 | 34.169 | 38.881 | 27.772 | 0.146309 |
| gate | 0.25 | 28.959 | 32.251 | 23.295 | 0.141132 |
| gate | 0.5 | 29.081 | 32.105 | 23.118 | 0.146712 |
| gate | 1 | 31.478 | 33.206 | 25.013 | 0.199314 |
| tilt | 0.25 | 30.429 | 34.185 | 24.476 | 0.142486 |
| tilt | 0.5 | 31.787 | 35.737 | 25.259 | 0.147997 |
| tilt | 1 | 34.681 | 38.631 | 26.889 | 0.179770 |

Boosting deeper wells increases errors. Weak shallow suppression offers a small km/s improvement but worsens training and validation proportional errors. Picking the lowest test RMSE would change the selection rule after observing the outcome.

An exploratory grouping by maximum ordinary-matter predicted speed explains the tradeoff: the weak gate retains about 83%, 91% and 93% of reference deposited mass in the <80, 80-160 and >=160 km/s groups. The weakest group loses the most inventory and its average proportional error grows. These group boundaries are post hoc diagnostics, not fitted or independently validated physics.

## Refinement and radial behavior

Because every family chose zero, the weak gate was refined only as a diagnostic of its metric-dependent gain; it does not replace the selected reference.

| Fine calculation | Train RMSE | Validation RMSE | Test RMSE |
|---|---:|---:|---:|
| Reference | 29.022 | 32.491 | 23.587 |
| Gate beta=0.25 | 28.968 | 32.261 | 23.307 |

Point-weighted mean speed errors for R/Rd <1, 1-3, >=3:
- reference: inner -8.549 km/s, middle -13.067 km/s, outer +6.623 km/s.
- gate_0.25: inner -8.554 km/s, middle -13.436 km/s, outer +4.404 km/s.

All coarse branches converge. Fine gate solutions also converge from initial density factors 1 and 0.1; maximum speed difference is 2.38704e-06 km/s. Maximum coarse/fine gate speed change is 0.120233 km/s. Fine reference differs from archived reference by at most 0.0393568 km/s. Analytic ray relative error is 2.53e-13; an independent analytic spherical-potential profile agrees within 0.000243 relative error. Numerical convergence is not dynamical stability.

## Interpretation and next step

Well depth can enter the math, but deepening the inner well also changes the gravity felt farther out. Increasing capture is not guaranteed to correct opposite inner and outer residuals. This test does not establish a shared improvement or solve lensing. Keep the unchanged reference. A useful next bounded hypothesis would redistribute a fixed captured inventory toward intermediate radii while limiting central accumulation, with its capture, support and energy rules specified before fitting. Do not add separate correction factors for each galaxy. Direction-dependent depth feedback and uncertainty in the unmeasured outer potential remain untested.

Reproduce: run depth-capture.py, depth-capture.py --refine, then depth-capture-report.py. Source hashes and all predicted curves are retained in the result JSON files.
