# Outer stellar motions with mass and orbit uncertainty

The empirical extra-force model still misses three of six outer measurements at the level of its central 95% conditional posterior predictive intervals. The baryonic baseline misses four. Neither count is a calibrated rejection probability for the physical hypothesis. These are previously exposed training galaxies, not the reserved independent prediction test.

## Observed values and conditional predictions

All entries are stellar root-mean-square line-of-sight velocities in km/s, not circular rotation speeds. Intervals include measurement covariance and mass/orbit parameter uncertainty, with uniform log mass and uniform orbit anisotropy as the reference prior.

| Galaxy | Observed outer bin | Baryonic 95% interval | Extra-force 95% interval |
|---|---:|---:|---:|
| J0037-0942 | 242.91 | 225.82 to 237.76 | 233.75 to 245.54 |
| J1112+0826 | 255.30 | 234.69 to 245.52 | 237.65 to 248.46 |
| J1204+0358 | 219.37 | 212.89 to 225.77 | 216.37 to 229.70 |
| J1402+6321 | 255.57 | 237.84 to 247.16 | 241.72 to 251.03 |
| J1621+3931 | 225.72 | 214.26 to 232.14 | 219.76 to 237.93 |
| J1630+4520 | 240.80 | 216.34 to 232.17 | 222.47 to 238.21 |

Using uniform mass instead changes any predictive mean by at most 0.0664 km/s. It leaves the inside/outside interval classifications unchanged. J1112+0826, J1402+6321 and J1630+4520 remain above the extra-force intervals.

## Calculation and formula provenance

The Jeans projection, conditional multivariate Gaussian, Bayes integration and trapezoidal quadrature are established mathematics. The extra-force profile is a frozen empirical prescription from the earlier galaxy fit; it is not derived from companion production, transport or capture.

Let t denote the inner bins and h the outer bin. Only the inner measurements enter the mass-and-anisotropy posterior. At each parameter pair:

\[\mu_{h|t}=\mu_h+C_{ht}C_{tt}^{-1}(y_t-\mu_t),\qquad
V_{h|t}=C_{hh}-C_{ht}C_{tt}^{-1}C_{th}.\]

The final prediction integrates this Gaussian over the inner-data posterior. Its variance includes both conditional measurement variance and variation of the conditional mean across that posterior. The outer observation is used only to evaluate the prediction, never to weight the parameter grid.

Mass bounds are 1e7 to 1e14 solar masses; constant anisotropy bounds are -2 to 0.45. Both priors are uniform in anisotropy. Light profiles, geometry, seeing, force coefficients, cutoff and covariance remain fixed.

## Numerical evidence and limits

All 24 model/prior calculations pass the declared coarse/fine grid checks: 801 by 321 versus 1601 by 641 mass/anisotropy nodes, changes no larger than 0.005 in observed predictive CDF and 0.1 km/s in mean or standard deviation. This checks parameter integration, not physical accuracy. Saved source hashes are checked when this report is regenerated.

Maximum posterior probability within 0.05 of an anisotropy boundary is 1.18%. Finite parameter bounds remain an assumption.

The six objects were selected by available released light components and the release use flag, as documented in the preceding audits. J1538 lacks the required component profile; J0330 has release use flag zero. Their absence must not be interpreted as successful predictions.

These intervals omit light-profile, seeing, geometry, force-law, covariance estimation and selection uncertainty. The two mass priors are a limited sensitivity check, not all possible prior choices. No population coverage claim or hypothesis p-value follows from six already-exposed systems.

## Consequence for the research goal

Mass and constant-orbit uncertainty alone does not remove the radial-profile mismatches under the current fixed assumptions. Further adjustment must have independent physical or observational support, rather than being selected solely to match these six outer bins. Lensing must then be recomputed from the same gravitational source. Companion-source derivation and the reserved prediction test remain unfinished.

Reproduce with `python research_work/results/slacs-outer-predictive/run.py`, then `python research_work/results/slacs-outer-predictive/summarize.py`.
