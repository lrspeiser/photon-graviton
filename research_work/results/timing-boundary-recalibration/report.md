# Continuous-scatter timing recalibration

All 160 exposed artificial samples were retained. This is a development comparison, not a test on untouched observations and not evidence of photon conversion. All declared source hashes and event-array hashes were verified.

| Shape | SNR | Truth b | Mean bias | Monte Carlo SE | Old accepted /20 | Revised accepted /20 | Revised exact 95% interval | Old/new invalid | Screen |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| split_gaussian | 5 | 0 | 0.04932 | 0.03509 | 13 | 20 | [0.832, 1.000] | 7/0 | True |
| split_gaussian | 5 | 1 | 0.01630 | 0.03528 | 16 | 20 | [0.832, 1.000] | 4/0 | True |
| split_gaussian | 20 | 0 | 0.00794 | 0.02153 | 20 | 20 | [0.832, 1.000] | 0/0 | True |
| split_gaussian | 20 | 1 | 0.00565 | 0.02129 | 20 | 20 | [0.832, 1.000] | 0/0 | True |
| split_gaussian_shoulder | 5 | 0 | 0.01240 | 0.03167 | 16 | 20 | [0.832, 1.000] | 4/0 | True |
| split_gaussian_shoulder | 5 | 1 | 0.00654 | 0.03383 | 18 | 19 | [0.751, 0.999] | 1/0 | True |
| split_gaussian_shoulder | 20 | 0 | 0.00187 | 0.02165 | 20 | 20 | [0.832, 1.000] | 0/0 | True |
| split_gaussian_shoulder | 20 | 1 | 0.01081 | 0.02250 | 19 | 19 | [0.751, 0.999] | 0/0 | True |

Acceptance means both valid fit and nominal likelihood-ratio interval membership at the injected exponent. Invalid fits count as unsuccessful acceptance, not necessarily a valid interval falsely excluding truth. The nominal chi-square threshold is a known asymptotic statistical approximation; this experiment checks its behavior rather than assuming it is exact. The normal width population and continuous interpolation are statistical tools, not new physical formulas.

The screening limits are unchanged: absolute mean bias at most 0.1, accepted fraction at least 0.85, and no invalid fits in each cell. Twenty trials per cell leave wide uncertainty even if screens pass. Seeds are paired across cells; do not pool them as independent trials. Original results are preserved, not selectively replaced.

Zero scatter is now permitted. Finite-support normalization, source priors and bounds still matter. Optimization start differences, width-grid sensitivity, independent simulation seeds, intrinsic source evolution, actual filters and calibrated brightness remain required. These results cannot establish a propagation law, absolute energy transfer or gravitational deposition. All six research goals remain open.
