# Combined stellar gradients and radial orbits at calibrated lens angles

13 September 2026. A conditional nuisance-model experiment on exposed data; no new companion law is adopted.

## Outcome

Combining stellar mass-to-light gradients with radial orbital anisotropy improves the motion fit beyond either separate extension, but does not supply a viable joint solution. Inner chi-squared decreases from 106.88/105.62 for the expanded radial-only configuration to 82.40/81.67 (Chabrier/Salpeter). The covariance-conditioned outer residual-square sum decreases from 115.71/114.32 to 103.20/101.32, still substantially above the free-mass motion control at 55.01/53.42.

The newly fitted orbital profiles also fail a necessary density-slope/anisotropy condition in three galaxies under both population assumptions. These failures occur in the extrapolated central light profiles. They rule out the specified separable augmented-density completion with these exact profiles, not all possible nonspherical, nonseparable or modified-central-profile constructions. Even the cases passing this check have no constructed positive distribution function or demonstrated stability.

## Model and provenance

The fixed exact-third companion law remains unchanged. Re is the fixed projected half-light radius. Use Upsilon(r)=Upsilon_out[1+h/(1+(r/Re)^2)] and beta(r)=beta0+(beta_infinity-beta0)r^2/(r^2+Re^2). The first is the earlier chosen phenomenological stellar gradient; the second is a specialization of the established [Baes and Van Hese anisotropy family](https://arxiv.org/abs/0705.4109). Combining them is a project diagnostic, not a first-principles derivation of companion physics. These parameters change the stellar mass distribution and orbit moments, not the capture exponent or shared companion amplitude.

For every trial gradient, the catalogue lens angle fixes stellar mass. The angle is consumed calibration, not a prediction. Fit three nuisance parameters to inner stellar-motion bins, using their covariance; the outer bin is evaluated through the same conditional residual as preceding experiments. Keep h in [-0.8,9], beta0 in [-2,0.45], beta_infinity in [-2,0.95]. The expanded outer bound came from prior adaptive work on these observations. No physical measurement establishes any of these bounds. Positive density follows from h>=-0.8, since the gradient factor is at least 0.2.

This experiment uses 20 starts per population/system, including both prior separate-fit solutions. Multiple local minima occur; the lowest successful result is reported. Between 18 and 20 starts succeed per case. This is not proof of a global optimum. Added freedom is expected to improve the inner fit and is not evidence for the radiation-derived origin of the fitted gravity.

## Results

Each population row is a different assumption applied to the same six galaxies, not an independent sample.

| Population | Inner chi-squared, radial / combined | Outer residual-square sum, radial / combined | Gradient / central-beta / outer-beta boundary counts |
|---|---:|---:|---:|
| Chabrier | 106.885 / 82.403 | 115.715 / 103.198 | 2 / 3 / 2 |
| Salpeter | 105.623 / 81.667 | 114.316 / 101.316 | 2 / 2 / 2 |

| Galaxy | Population | h | beta0 | beta_infinity | Inner chi-squared | Outer standardized residual | Slope check |
|---|---|---:|---:|---:|---:|---:|---|
| J0037-0942 | Chabrier | 9.0000 | -0.0534 | 0.9500 | 22.096 | 3.160 | No sampled violation |
| J0037-0942 | Salpeter | 9.0000 | -0.0366 | 0.9500 | 22.061 | 2.965 | No sampled violation |
| J1112+0826 | Chabrier | -0.6700 | 0.3438 | 0.7286 | 12.185 | 5.143 | No sampled violation |
| J1112+0826 | Salpeter | -0.6701 | 0.3486 | 0.7355 | 12.147 | 5.140 | No sampled violation |
| J1204+0358 | Chabrier | 9.0000 | -1.1839 | 0.5599 | 3.088 | 5.274 | No sampled violation |
| J1204+0358 | Salpeter | 9.0000 | -1.1572 | 0.5745 | 3.103 | 5.268 | No sampled violation |
| J1402+6321 | Chabrier | -0.5964 | 0.4500 | 0.9500 | 40.283 | 5.366 | Violation |
| J1402+6321 | Salpeter | -0.5894 | 0.4500 | 0.9500 | 39.621 | 5.316 | Violation |
| J1621+3931 | Chabrier | -0.6499 | 0.4500 | 0.5336 | 2.451 | -1.685 | Violation |
| J1621+3931 | Salpeter | -0.6420 | 0.4480 | 0.5326 | 2.449 | -1.709 | Violation |
| J1630+4520 | Chabrier | -0.5717 | 0.4500 | 0.6404 | 2.299 | 2.704 | Violation |
| J1630+4520 | Salpeter | -0.5660 | 0.4500 | 0.6425 | 2.286 | 2.679 | Violation |

J0037-0942 and J1204+0358 reach h=9, a central-to-outer mass-to-light ratio of ten. Four other systems prefer negative gradients, with central mass-to-light ratios around one-third to two-fifths of the outer value. These are fit choices, not measurements of stellar populations. A universal inward-increasing gradient is not supported by these best fits. The physical stellar-population plausibility remains untested.

## Necessary orbital check

Apply gamma=-d ln nu/d ln r>=2 beta on the same 4097 log radii, 1e-6 to 100 Re, where nu is the fixed luminosity tracer density. For spherical separable augmented density and beta0<=1/2, this is a known necessary condition for a nonnegative distribution function: [Van Hese, Baes and Dejonghe](https://arxiv.org/abs/1010.4301). The earlier radial-only result cannot be transferred to changed parameters.

The minimum margin is -0.07228. Failing cases and the outermost violating sampled radius are:

| Galaxy | Population | Minimum gamma-2 beta | Outermost violating sampled r/Re |
|---|---|---:|---:|
| J1402+6321 | Chabrier | -0.072276 | 1.5820447e-05 |
| J1402+6321 | Salpeter | -0.072276 | 1.5820447e-05 |
| J1621+3931 | Chabrier | -0.057799 | 7.7387689e-06 |
| J1621+3931 | Salpeter | -0.053739 | 6.8847888e-06 |
| J1630+4520 | Chabrier | -0.072276 | 1.5820447e-05 |
| J1630+4520 | Salpeter | -0.072276 | 1.5820447e-05 |

All violations reach the inner grid boundary at 1e-6 Re; this calculation does not measure the actual light or orbits there. Central-profile uncertainty is therefore material. At Re the condition passes in every case. Allowing a different unresolved central profile could change the check, but it must be specified and tested, not silently substituted. No constrained refit or new positive distribution function has yet been constructed.

## Verification and limits

The h=0 limit reproduces previous radial-only inner chi-squared to within 5.0e-13; equal endpoint anisotropies reproduce the gradient-only result to within 9.1e-13. The imposed lens equation closes to fractional error below 1e-10. Deprojection order doubling and the independent log-density derivative checks retain the numerical tolerances documented in orbit-slope-check-report.md. These verify calculation consistency, not physical validity. This experiment does not introduce a new spatial grid or establish numerical convergence of an arbitrary distribution-function inversion.

The improved inner score does not justify replacing the companion law. Remaining failures divide into: poor conditional outer motions; nuisance-parameter boundary dependence; central orbital consistency in the specified construction; and missing independent evidence for extreme stellar gradients. The next distinguishing calculation is a fit constrained by orbital consistency, followed by an actual distribution-function construction if residuals warrant it. Alternative central light profiles or nonspherical orbit models must be treated as distinct model changes with their own evidence.

## Reproduction

Run `python research_work/results/companion-extensions/gradient-orbits.py`, then `python research_work/results/companion-extensions/orbit-slope-check.py --combined`. The latter writes a distinct gradient-orbits-slope-results.json and preserves the earlier radial-only results. Input hashes, all predicted and measured motion bins, individual scores, optimizer outcomes and boundaries are recorded. The pre-execution choices are in gradient-orbits-protocol.md.
