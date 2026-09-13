# Expanded synthetic timing calibration results

All 160 frozen base cases and eight selected numerical refinements completed. These are artificial fluxes on exposed DES cadence patterns, not measurements of astronomical time dilation.

| Shape | SNR | True b | Mean bias | Monte Carlo SE | Accepted nominal-95% truth membership | Exact 95% interval for acceptance | Invalid fits |
| --- | ---: | ---: | ---: | ---: | --- | --- | ---: |
| split_gaussian | 5 | 0 | +0.0489 | 0.0350 | 13/20 | [0.408, 0.846] | 7 |
| split_gaussian | 5 | 1 | +0.0159 | 0.0354 | 16/20 | [0.563, 0.943] | 4 |
| split_gaussian | 20 | 0 | +0.0079 | 0.0215 | 20/20 | [0.832, 1.000] | 0 |
| split_gaussian | 20 | 1 | +0.0056 | 0.0213 | 20/20 | [0.832, 1.000] | 0 |
| split_gaussian_shoulder | 5 | 0 | +0.0121 | 0.0317 | 16/20 | [0.563, 0.943] | 4 |
| split_gaussian_shoulder | 5 | 1 | +0.0066 | 0.0338 | 18/20 | [0.683, 0.988] | 1 |
| split_gaussian_shoulder | 20 | 0 | +0.0019 | 0.0216 | 20/20 | [0.832, 1.000] | 0 |
| split_gaussian_shoulder | 20 | 1 | +0.0108 | 0.0225 | 19/20 | [0.751, 0.999] | 0 |

Screening thresholds: absolute mean bias <=0.1, at least 85% nominal-95 coverage, and zero invalid fits in each cell. These were declared before execution. They are not precision validation of 95% coverage. Per-cell screening results: [False, False, True, True, False, False, True, True].

## Selected numerical refinements

| Refined case | Change in b | Mean event log-likelihood change | Pass |
| --- | ---: | ---: | --- |
| split_gaussian-snr5-b0-seed917-refined | 0.00006 | 0.00423 | True |
| split_gaussian-snr5-b1-seed911-refined | 0.00105 | 0.00393 | True |
| split_gaussian-snr20-b0-seed920-refined | 0.00891 | 0.08768 | True |
| split_gaussian-snr20-b1-seed919-refined | 0.00373 | 0.06405 | True |
| split_gaussian_shoulder-snr5-b0-seed911-refined | 0.00252 | 0.00456 | True |
| split_gaussian_shoulder-snr5-b1-seed915-refined | 0.00054 | 0.00640 | True |
| split_gaussian_shoulder-snr20-b0-seed919-refined | 0.00311 | 0.08632 | True |
| split_gaussian_shoulder-snr20-b1-seed919-refined | 0.00501 | 0.07029 | True |

Selected refinements retain both estimates and are not substituted selectively into the base coverage counts. Refining the largest estimation error does not guarantee that every other case has converged numerically.

## Interpretation

The exact binomial intervals show how uncertain these coverage fractions remain with 20 seeds per cell. Seeds are paired across cells; aggregate counts cannot be treated as independent. Under the frozen protocol, invalid fits count as unsuccessful interval acceptance and are explicitly listed. This operational rate is not the fraction of otherwise valid intervals excluding truth: all 16 invalid cases hit the scatter floor, and LR-only membership is reported separately in the base diagnosis. No case is reclassified. The bias Monte Carlo standard error describes variation over these synthetic draws, not total astrophysical uncertainty.

All eight selected refinements pass the preset numerical gate; the largest exponent change is below 0.009. This does not repair the 16 base scatter-boundary failures. A separate boundary-capable estimator and full recalibration are the next step, preserving this frozen baseline. Source evolution, real filter transmission, spectroscopic selection, amplitude-prior sensitivity and physically calibrated brightness remain outside this batch. A numerical or screening pass does not establish that photon-companion conversion explains observed redshift, timing or brightness. All six scientific goals remain open.
