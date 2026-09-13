# Fitting all stellar-motion bins with lens-calibrated mass and constrained orbits

13 September 2026. Descriptive all-bin diagnostic, not a held-out prediction. Exact-third companions remain unchanged.

## Outcome

Fitting all 40 motion bins across the six galaxies reduces total covariance-weighted motion chi-squared from 193.62/190.91 to 161.81/159.83 (Chabrier/Salpeter). This is an improvement of about 16%, but significant residual structure remains within the tested model. The improvement exchanges a worse inner fit for a better outer fit: inner totals rise from 87.58/86.79 to 104.25/103.00, while conditional outer residual-square sums fall from 106.03/104.12 to 57.57/56.84.

J1402+6321 contributes approximately 47% of the remaining total under either population assumption. The all-bin fit retains an outer standardized residual around 5.3 for that system. This is a conditional residual based on the released covariance and model assumptions, not a model-independent exclusion significance. Multiple numerical starts do not prove a global impossibility result.

The formerly reserved outer bin has now been used in optimization. Its lower residual is not an independent predictive success. The result tests how well the family can describe all available bins under the same lens calibration and necessary orbital constraints.

## Fixed physics, fitted choices and provenance

Retain the adjusted-reference exact-third companion density and galaxy-trained amplitude. Keep conditional nonexpanding geometry, light profiles, PSF, projected half-light radius Re and covariance unchanged. At every trial stellar gradient, the catalogue lens angle fixes stellar mass. No independent lensing prediction is made here, and the companion capture law is not refitted.

The stellar parameters remain h, beta0 and beta_infinity in the previous gradient and radial anisotropy equations. The gradient is a project-selected phenomenological mass-to-light prescription; the orbital family is known stellar dynamics. Enforce the same positive-gradient factor, endpoint parameter bounds, central asymptotic condition beta0<=0.375 and sampled gamma>=2 beta condition as in constrained-gradient-orbits-report.md. The latter remains necessary for the specified separable augmented-density completion, not sufficient for a positive distribution function or stability.

The only new optimization choice is the objective: include the outer measurement in the full covariance-weighted motion residual. With residual vector e=y-model and covariance C,

    chi_squared_all = e^T C^-1 e
                    = chi_squared_inner + (conditional outer residual)^2.

This is the standard block-covariance/Schur-complement identity, not new companion physics. It provides a direct check of consistent treatment of the correlations. The total must use the conditional outer residual, not add the raw outer residual independently to the inner score.

The experiment uses 22 starts per system/population, including the earlier constrained inner-only optimum. All 22 converge in every case; local minima and objectives remain recorded. Three stellar nuisance parameters per galaxy are optimized, with stellar mass set by lens calibration rather than freely optimized. No statistical model-selection ranking is inferred from the descriptive totals or nominal parameter count.

## Aggregate results

The same six galaxies and 40 bins are used under each population proxy. The two proxies are alternate assumptions, not independent observations.

| Population | Previous inner-only fit: all-bin chi-squared | New all-bin chi-squared | New inner contribution | New conditional outer contribution |
|---|---:|---:|---:|---:|
| Chabrier | 193.617 | 161.812 | 104.246 | 57.566 |
| Salpeter | 190.914 | 159.834 | 102.996 | 56.838 |

The historical free-mass, constant-anisotropy motion fit has smaller inner residuals but does not satisfy the same lens calibration. It is useful context, not a matched all-bin model ranking. This diagnostic neither improves the galaxy rotation law nor establishes an advantage over MOND or dark-halo models.

## Per-system results

| Galaxy | Population | h | beta0 | beta_infinity | All-bin chi-squared | Inner chi-squared | Conditional outer residual |
|---|---|---:|---:|---:|---:|---:|---:|
| J0037-0942 | Chabrier | 0.5665 | 0.0957 | 0.9500 | 30.155 | 25.887 | 2.066 |
| J0037-0942 | Salpeter | 0.7456 | 0.0922 | 0.9500 | 29.292 | 25.338 | 1.989 |
| J1112+0826 | Chabrier | 0.0721 | -0.2258 | -2.0000 | 25.996 | 14.122 | 3.446 |
| J1112+0826 | Salpeter | 0.0971 | -0.2305 | -2.0000 | 26.038 | 14.164 | 3.446 |
| J1204+0358 | Chabrier | 9.0000 | -1.6372 | 0.8443 | 20.136 | 10.668 | 3.077 |
| J1204+0358 | Salpeter | 9.0000 | -1.6005 | 0.8498 | 20.157 | 10.636 | 3.086 |
| J1402+6321 | Chabrier | -0.5061 | 0.3750 | 0.9500 | 75.860 | 47.203 | 5.353 |
| J1402+6321 | Salpeter | -0.4960 | 0.3750 | 0.9500 | 74.587 | 46.442 | 5.305 |
| J1621+3931 | Chabrier | -0.1247 | -0.0191 | 0.3498 | 3.447 | 3.129 | -0.564 |
| J1621+3931 | Salpeter | -0.0754 | -0.0299 | 0.3465 | 3.466 | 3.143 | -0.569 |
| J1630+4520 | Chabrier | 0.8609 | -0.0239 | -2.0000 | 6.217 | 3.236 | 1.727 |
| J1630+4520 | Salpeter | 0.9445 | -0.0278 | -2.0000 | 6.293 | 3.273 | 1.738 |

J1112+0826 and J1630+4520 move from positive outer anisotropy in the inner-only fit to beta_infinity=-2, the tangential endpoint bound. Their fitted stellar gradients also change sign. This sensitivity shows that the previously extrapolated orbital structure was not robustly constrained by the inner bins alone. The physical interpretation is not that real stars changed their orbits; different data constraints select different model parameters.

J0037-0942 and J1402+6321 retain the radial endpoint beta_infinity=0.95, while J1204+0358 retains h=9. Four galaxies thus reach an outer anisotropy bound, and one reaches the gradient bound, under both population assumptions. J1402+6321 also reaches the central slope limit. No extreme parameter is treated as an independently measured property.

J1621+3931 is comparatively well described in this diagnostic. This does not justify selecting it alone as evidence for the theory: all six targets and both population assumptions remain in the assessment.

## Verification and limitations

The full covariance score equals the sum of inner score and squared conditional outer residual to within 2.67e-14. Evaluating the previous inner-only parameters reproduces their all-bin score to within 1e-5, and every optimized total is no worse than that feasible starting point. The earlier gradient-only and radial-only inner-score reproduction checks remain active. All stellar masses and predicted moments are positive; imposed lens angles close to fractional error below 1e-10.

Every fitted profile satisfies the central limit and the refined 8193-point slope check; the minimum finite-radius margin is 0.0777244. Positivity and stability of the full distribution function still do not follow. These are numerical consistency checks in the inherited spherical Jeans approximation, not proof that any model fits all physical observables. The calculation does not include a new raw lens-image reconstruction, revised light profile, evolving source distribution or companion microphysics.

## Interpretation and next distinguishing work

The bad outer-bin behavior was partly an extrapolation problem: including it moves the solution and reduces its residual. It was not solely an extrapolation problem, because the all-bin optimum still retains a substantial mismatch, especially in J1402+6321. The present family has not produced a demonstrated joint solution.

The remaining possibilities include restrictive stellar-gradient/orbit forms, fixed transition scale, spherical geometry, central light-profile uncertainty, or a problem in the adopted gravitational/optical response. These are different changes, not interchangeable explanations. A bounded next test can vary the fixed orbit transition scale while retaining the lens-calibrated mass, necessary orbital condition and explicit complexity cost, or obtain independent stellar-population constraints to prevent gradients from merely absorbing discrepancies. No such change is assumed to succeed here.

The underlying photon-conversion, energy supply, retention and propagation problems remain open. Improving stellar nuisance fits alone cannot establish the radiation-derived origin of the extra gravity.

## Reproduction

Run `python research_work/results/companion-extensions/gradient-orbits.py --slope-constrained --all-motion-bins`. The script requires the constraint flag for this diagnostic and writes all-motion-orbits-results.json without replacing either earlier fit. The output records all predicted and observed bins, total and decomposed scores, necessary-condition margins, bounds, optimizer results and input hashes. Choices were written before execution in all-motion-orbits-protocol.md.
