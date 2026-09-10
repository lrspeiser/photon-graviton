# Reserved validation result for the frozen lens pilot

The frozen empirical companion pilot performs worse than the ordinary-matter benchmark on all seven eligible reserved validation lenses, judged by absolute Einstein-angle residual under the primary assumptions. The companion model overpredicts all seven published SIE angles. This is a negative validation result for the tested approximation, not a reason to silently retune its force law or drop lenses.

## Primary result

The primary case was fixed before the validation calculation: updated I-band effective radii, Gaussian seeing FWHM 1.5 arcsec, companion source boundary 20 effective radii, spherical isotropic Hernquist tracers, the archived static geometry, and equal temporal/spatial metric potentials. All empirical coefficients and the mass-from-dispersion procedure are unchanged. Each validation galaxy's measured dispersion is a predictor used to infer its mass; its Einstein angle is the withheld target, not a mass-fitting input.

| Model | Training RMS, 33 systems | Validation RMS, 7 systems | Validation median predicted / SIE angle |
|---|---:|---:|---:|
| Ordinary-matter benchmark | 0.1840 arcsec | 0.1712 arcsec | 1.0614 |
| Frozen empirical companion | 0.2083 arcsec | 0.2505 arcsec | 1.1221 |

The companion model has a larger absolute residual than the benchmark for each of the seven systems. It is median 12.2% high in angle. The ordinary-matter benchmark overpredicts six of seven. No calibrated probability, significance or full-likelihood model preference is claimed.

Why this matters in plain language: the same rule that overshot the training galaxies also overshoots galaxies whose prediction errors had not been used to build or adjust the pilot. Its difficulty is therefore not confined to the original training examples. But the pilot's spherical geometry, unconstrained stellar mass normalization and orbital assumptions may themselves be inadequate; the result does not isolate photon conversion as the cause.

## Sample and exposure

The nine original validation systems were filtered using the same criteria as training: early-type classification and available spectroscopic dispersion. J0157-0056 lacks a dispersion; J1103+5322 is unclassified. The other seven are all retained:

J0029-0055, J0216-0813, J0912+0029, J0946+1006, J1420+6019, J1531-0105 and J1627-0053.

The public catalog, basic fields and sample summaries had been inspected earlier. This was a reserved prediction-residual test, not pristine blindness. The validation results are now exposed and must be treated as development evidence for any future revisions. The separate test predictions have **not** been evaluated. The command currently permits training and validation roles only.

## Prespecified sensitivity cases

| Assumed seeing FWHM | Boundary / effective radius | Companion validation RMS | Median predicted / SIE angle |
|---|---:|---:|---:|
| 0 arcsec | 5 | 0.2281 arcsec | 1.1027 |
| 0 arcsec | 20 | 0.2428 arcsec | 1.1139 |
| 0 arcsec | 100 | 0.2475 arcsec | 1.1174 |
| 1.5 arcsec | 5 | 0.2355 arcsec | 1.1108 |
| 1.5 arcsec | 20 | 0.2505 arcsec | 1.1221 |
| 1.5 arcsec | 100 | 0.2552 arcsec | 1.1256 |

The zero-seeing ordinary-matter benchmark has RMS 0.1656 arcsec. None of the listed companion sensitivity cases improves upon the corresponding ordinary-matter benchmark in aggregate RMS. These are assumed seeing and cutoff values, not per-system measurements or fitted capture radii. We do not promote the numerically best sensitivity case after inspecting validation.

## Verification and limits

The protocol file was written before validation prediction execution; its hash is recorded in the result. The verification script confirms unchanged shared parameters, exact eligible-system membership, no residual-based exclusions, and correct validation roles. Doubling grid/quadrature resolution changes any validation prediction by less than 0.000068 arcsec. Original training files and results remain preserved.

The formulas are the same known Hernquist, Jeans and weak-lensing equations plus the project's previously fitted empirical power law and explicit geometric postulates described in report.md. No new physical law or novelty claim is introduced here. The updated photometry is independently sourced but its earlier age/mass-to-light assumptions have not been fully refitted; the mass-normalization problem identified in lens-photometric-audit remains.

In particular, this is not an imaging likelihood. It compares a spherical predicted angle with an intermediate-axis SIE summary. It lacks shape, anisotropy, per-object seeing, environmental and lens-model covariance treatment. Dispersion error propagation alone does not cover those uncertainties, and fitting a mass from dispersion is not an independent stellar-mass measurement. A good ordinary-matter score therefore does not prove that its required masses are photometrically allowed.

## Decision

Record the present frozen pilot as underperforming on validation. Preserve the test systems while constructing a physically specified revision, a consistent stellar-population likelihood and a more realistic shared dynamical/lensing treatment. Do not adjust each lens independently until it agrees or claim the empirical rotation fit has supplied the missing photon-production, capture or metric equations. The timing, energy, transport, lensing and broad observational requirements remain active.

Reproduce with run.py --role validation --updated-profile --refined and the corresponding coarser run, then verify_validation.py. Full per-system outputs, all 56 system/configuration combinations, protocol and numerical checks are archived alongside this report. Source provenance and earlier formulas are linked from the existing training and data-readiness reports.
