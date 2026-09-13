# Combined gradient/orbit fits constrained by a necessary orbital condition

13 September 2026. Conditional refit on the same exposed lens systems. No change to exact-third companions or the lens geometry.

## Outcome

Most of the unconstrained combined fit's improvement survives imposing the necessary density-slope/anisotropy condition. Inner chi-squared rises from 82.40/81.67 to 87.58/86.79 (Chabrier/Salpeter), approximately 6.3% in both cases. The constrained fit is still better than radial-only fits at 106.88/105.62. Outer residual-square sums rise from 103.20/101.32 to 106.03/104.12, still roughly twice the free-stellar-mass motion control at 55.01/53.42.

All twelve configurations satisfy the imposed central asymptotic bound and the independently refined finite-radius check. Three galaxies under both populations reach the central bound exactly. These results remove the previously identified necessary-condition violations; they do not establish a positive distribution function, stable orbits, or an acceptable joint motion/lensing model.

## Constraint and mathematical provenance

The retained orbit family is beta(r)=beta0+(beta_infinity-beta0) f(r), with f(r)=r^2/(r^2+Re^2). The stellar gradient is Upsilon(r)=Upsilon_out[1+h/(1+(r/Re)^2)]. These are the known orbital family and project-selected phenomenological stellar gradient from the preceding experiment, not new companion interactions.

For a spherical separable augmented-density completion with beta0<=1/2, gamma=-d ln nu/d ln r>=2 beta is a known necessary condition for a nonnegative distribution function. The theorem and scope are documented by [Van Hese, Baes and Dejonghe](https://arxiv.org/abs/1010.4301). It is not a universal sufficient condition for arbitrary stellar orbit models. Nu is the luminosity tracer, not total gravitating density.

Rearranging this condition on the radial grid gives a cap at each beta_infinity:

    beta0 <= min_r [(gamma(r)/2 - beta_infinity*f(r))/(1-f(r))].

This is algebra, not an additional physical law. We also enforce the existing beta0<=0.45 bound and the analytical central limit of the adopted light model. For a positive Sersic component with n>1, Abel deprojection has leading central behavior nu proportional to r^(1/n-1). This follows by substituting projected radius R=r*cosh(u) in the derivative of the central Sersic surface profile; the resulting integral has a finite leading coefficient for n>1. In a positive mixture the largest n supplies the most divergent central cusp. An n=1 logarithmic term is subdominant.

All six adopted light models contain n=4 as their largest index. Consequently gamma(0)=1-1/4=0.75 and beta0<=0.375. This is an implication of the adopted extrapolated light profiles, not a measured central anisotropy. The gradient is finite and positive at zero and does not change that tracer cusp. No unresolved light profile is silently changed to rescue the fit.

A unit-interval numerical coordinate maps beta0 from -2 to the tightest cap. This enforces the inequalities without a penalty term or a new fitted physical parameter. The finite grid spans 1e-6 to 100 Re with 4097 log points; the final check doubles that to 8193 points. Satisfying the limit at zero plus a finite grid is still not an all-radii mathematical proof, nor a distribution-function inversion in the actual potential.

## Fitting and comparison

All capture constants, companion normalization, light profiles, Re, PSF, geometry and covariance are unchanged. At each gradient the catalogue lens angle fixes stellar mass. This consumes lensing calibration; it is not a fresh lens prediction. The three stellar nuisance parameters are fitted only to inner motion bins. The outer bin is evaluated using its covariance-conditioned residual and is not used for selection.

There are 21 optimizer starts per case: the preceding 20 starts plus the unconstrained combined optimum mapped into the permitted region. Twenty or twenty-one converge per case. Multiple local minima are retained in the record; the minimum successful objective is selected. This does not prove a global optimum. The feasible radial-only solutions remain a comparison control and each new inner score is no worse than that control.

| Population | Free-mass inner | Radial-only inner | Unconstrained combined inner | Constrained combined inner | Constrained outer sum |
|---|---:|---:|---:|---:|---:|
| Chabrier | 45.812 | 106.885 | 82.403 | 87.583 | 106.034 |
| Salpeter | 45.268 | 105.623 | 81.667 | 86.789 | 104.125 |

| Galaxy | Population | h | beta0 | beta_infinity | Inner chi-squared | Outer standardized residual | Central bound active |
|---|---|---:|---:|---:|---:|---:|---|
| J0037-0942 | Chabrier | 9.0000 | -0.0534 | 0.9500 | 22.096 | 3.160 | No |
| J0037-0942 | Salpeter | 9.0000 | -0.0366 | 0.9500 | 22.061 | 2.965 | No |
| J1112+0826 | Chabrier | -0.6700 | 0.3438 | 0.7286 | 12.185 | 5.143 | No |
| J1112+0826 | Salpeter | -0.6701 | 0.3486 | 0.7355 | 12.147 | 5.140 | No |
| J1204+0358 | Chabrier | 9.0000 | -1.1839 | 0.5599 | 3.088 | 5.274 | No |
| J1204+0358 | Salpeter | 9.0000 | -1.1572 | 0.5745 | 3.103 | 5.268 | No |
| J1402+6321 | Chabrier | -0.4595 | 0.3750 | 0.9500 | 45.196 | 5.730 | Yes |
| J1402+6321 | Salpeter | -0.4483 | 0.3750 | 0.9500 | 44.475 | 5.677 | Yes |
| J1621+3931 | Chabrier | -0.5923 | 0.3750 | 0.5112 | 2.473 | -1.502 | Yes |
| J1621+3931 | Salpeter | -0.5842 | 0.3750 | 0.5110 | 2.469 | -1.531 | Yes |
| J1630+4520 | Chabrier | -0.4770 | 0.3750 | 0.4930 | 2.544 | 2.587 | Yes |
| J1630+4520 | Salpeter | -0.4695 | 0.3750 | 0.4965 | 2.534 | 2.567 | Yes |

The first three systems retain effectively the same solutions as the unconstrained fit. J1402+6321, J1621+3931 and J1630+4520 move to beta0=0.375. The larger cost is dominated by J1402+6321. J1112+0826, J1204+0358 and J1402+6321 still have conditional outer residuals around five to six. These depend on the data covariance and stipulated model; they are not model-independent exclusion significances.

J0037-0942 and J1204+0358 still reach h=9, while other galaxies prefer negative gradients. The accompanying stellar-population interpretation remains unverified. J0037-0942 and J1402+6321 still reach beta_infinity=0.95. Enforcing one necessary condition does not resolve these boundary dependences.

## Verification

The minimum refined finite-radius margin gamma-2 beta is 0.0777244; the minimum central-limit margin is zero, attained by the six central-bound cases. Both populations are checked on all 8193 radii with the same order-512 deprojection used in the original verified slope evaluation. The old radial-only and constant-anisotropy gradient limits continue to reproduce their stored inner chi-squared to better than 1e-4 (actual differences below 1e-12). Imposed lens angles close below 1e-10 fractional error. These checks validate the implementation within the adopted approximation, not the existence or stability of real stellar distributions.

The common Abel slope routine was extracted into orbit_density.py without changing its mathematics. Re-running the preceding unconstrained combined-fit check reproduced its JSON byte-for-byte. The constrained output is separate and does not overwrite historical fits.

## What this changes

The earlier orbital-condition failure was repairable with a modest fitting cost in this family. It was not the main cause of the remaining joint mismatch. The next physical assessment still requires a nonnegative distribution function in the fitted potential, and ultimately stability. Even a successful construction would leave large outer residuals and extreme stellar-gradient assumptions to address. Alternative central light profiles, nonspherical dynamics and independently constrained stellar populations remain distinct avenues; none is assumed to fix the discrepancy in this report.

No result here derives photon conversion, absolute energy supply or the one-third exponent. Those tracks remain open. This experiment isolates a stellar-dynamical limitation so it is not confused with evidence for or against the underlying radiation mechanism.

## Reproduction

Run `python research_work/results/companion-extensions/gradient-orbits.py --slope-constrained`. It writes constrained-gradient-orbits-results.json, recording every observed and predicted motion bin, optimizer outcome, parameter bound, central and refined-grid margins, and input hashes. Choices were recorded before execution in constrained-gradient-orbits-protocol.md. The historical unconstrained run remains available without the flag.
