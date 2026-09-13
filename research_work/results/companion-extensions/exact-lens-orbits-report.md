# Can constant orbital anisotropy accommodate the lens-required stellar masses?

13 September 2026. Conditional exact-angle compatibility experiment.

## Result

No tested branch recovers the previous stellar-motion fit after the stellar mass is forced to match the catalogue lens angle. Across the six systems, total inner-stellar chi-squared increases from about 41-46 to 585-625. The sums of squared conditional outer-bin residuals also increase, from about 43-55 to 120-136. The same three systems identified in the bending-budget audit dominate the inner-fit cost: J0037-0942, J1204+0358 and J1402+6321.

Constant anisotropy is refitted and all recovered optima remain inside the existing bounds. Thus this particular adjustment does not reconcile the data within the assumed mass profiles and geometry. It does not exclude radially varying orbital structure, stellar-population gradients, changes in companion spatial shape or revised optical geometry.

Zero lens-angle error in this experiment is imposed calibration, not a successful prediction. The catalogue angles are treated as exact targets solely to quantify compatibility cost; no lens-angle uncertainty likelihood or raw-image fit is supplied, so the chi-squared changes are not translated into exclusion probabilities.

## Calculation

At the catalogue impact parameter, write the fixed-profile effective deflection as

\[
\widehat\alpha=\frac{M_\star}{10^{11}M_\odot}\widehat\alpha_{\star,11}+\widehat\alpha_c.
\]

With the same geometry used previously, the exact required mass is

\[
M_{\star,E}=10^{11}M_\odot\frac{\theta_E D_s/D_{ls}-\widehat\alpha_c}{\widehat\alpha_{\star,11}}.
\]

This is ordinary algebra using the known linear dependence of the adopted weak-field deflection on mass. It is an inverse constraint, not a new physical law or a source-energy derivation. All required masses are positive and within the existing 1e7-1e14 Msun numerical bounds for this sample.

Freeze each branch's companion density, galaxy-trained shared amplitude, capture scale, opacity and exact-third response. Fix M_star to the value above, then fit only constant beta=1-sigma_t^2/(2 sigma_r^2) to the inner Vrms bins using their released covariance. This beta convention uses sigma_t squared as the sum of the two tangential dispersions. Retain the prior -2<=beta<=0.45 bounds and four starting values. The outermost bin stays out of the fit and is evaluated conditionally on inner residuals as before.

The stellar mass and lens angle remain tied to the same fixed light-profile shape and conditional distance mapping. None of the capture or optical parameters is adjusted on the lenses. Stellar population proxies and the two nuisance assumptions are not promoted to measured independent masses.

## Aggregate results

Each row contains six systems. The free fit adjusts stellar mass and constant beta using inner motions; the exact-angle case fixes mass using the lens constraint and adjusts beta using those same inner motions. Lower residual sums are better, but these summaries are not combined into an undocumented joint statistic.

| Profile | Population | Free inner chi-squared | Exact-angle inner chi-squared | Free outer residual square sum | Exact-angle outer residual square sum |
|---|---|---:|---:|---:|---:|
| Original reference | Chabrier | 42.16 | 623.50 | 44.32 | 122.56 |
| Original reference | Salpeter | 41.55 | 625.37 | 42.70 | 120.07 |
| Adjusted reference | Chabrier | 45.81 | 611.48 | 55.01 | 136.23 |
| Adjusted reference | Salpeter | 45.27 | 613.61 | 53.42 | 134.17 |
| Local filling | Chabrier | 41.63 | 586.25 | 50.41 | 128.44 |
| Local filling | Salpeter | 40.86 | 585.35 | 48.69 | 125.91 |
| Full recycling | Chabrier | 42.89 | 586.81 | 54.92 | 132.83 |
| Full recycling | Salpeter | 42.33 | 586.65 | 53.58 | 130.98 |

All eight branch/population cases have the same three systems with individual inner chi-squared increases above 50. That threshold is a descriptive grouping, not a significance criterion. The other systems have smaller costs and must not be described as exhibiting the same severe discrepancy. Complete masses, mass ratios, betas, angle residuals and per-system motion costs are retained in the machine-readable results.

## What the result means

The bending-budget diagnosis was not merely an artifact of holding beta at its previous best fit. Refitting constant beta does not allow the lens-required stellar mass to retain the previous quality of the inner-motion fit. New shared capture amplitudes from the local and recycling branches do not change that conclusion.

This identifies an incompatibility within a specific combination of assumptions: fixed projected light components, constant stellar mass-to-light ratio, spherical dynamics with one anisotropy, the specified companion profile, and the inherited conditional optical geometry. It does not identify which assumption is physically wrong. The released covariance, lens-model summaries and unmodeled structural uncertainties also matter; the present large numerical costs cannot be interpreted as a model-independent rejection of companions.

The next distinguishing models should change a specified piece of that combination, such as a physically constrained radial stellar mass-to-light profile or orbital distribution, and use both motions and lensing to assess it. Refitting independent values for every lens without a source or population rationale would add flexibility without deriving our gravity mechanism. Likewise, changing the lens geometry requires checking the propagation and distance assumptions jointly rather than inserting target-specific scale factors.

## Verification and reproduction

Run `python research_work/results/isotropic-galaxy-transfer/lensing.py --capacity-branch=BRANCH --capacity-exact-lens` for original, reference, local and recycling, then run `python research_work/results/companion-extensions/exact-lens-orbits.py`. New filenames preserve the unconstrained outputs. See the [protocol](exact-lens-orbits-protocol.md).

The verifier checks all 48 branch/population/system cases: identical observed arrays, geometry and population mapping; unchanged shared amplitudes; positive bounded masses; successful optimizer starts; no beta boundary; exact-angle residual below 1e-7; and constrained inner minima no better than unconstrained minima beyond numerical tolerance. For the original and local branches it also compares the inferred mass ratios with the independent earlier bending decomposition on the same grids. Input/output and engine hashes are recorded in `exact-lens-orbits-results.json`. The stored scope explicitly marks the lens angle as consumed calibration.

This calculation tests compatibility within the retained constant-anisotropy model. It does not prove global optimality for arbitrary orbital distributions, dynamic stability, independent distances, lens-image agreement or the complete photon-companion theory. The reference is not replaced and no physical goal is declared complete.
