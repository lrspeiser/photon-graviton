# Radially varying orbital anisotropy at the lens-required stellar masses

13 September 2026. Known orbital model used as a conditional nuisance-model test.

## Outcome

Radial variation in anisotropy reduces the stellar-motion cost of imposing the catalogue lens angles, but does not establish a joint solution. With both endpoint anisotropies restricted to the inherited -2..0.45 range, inner chi-squared totals fall from 611.48/613.61 to 447.76/453.84 (Chabrier/Salpeter). The three previously difficult systems all reach the outer endpoint bound.

A follow-up allowing their outer endpoint to reach 0.95 reduces the combined totals further to 106.88/105.62. Thus the earlier 0.45 restriction materially limited the fit. However, all three reach the new limit too, and combined conditional outer residual sums remain 115.71/114.32, compared with 55.01/53.42 for the free constant-anisotropy motion fits. Strongly radial outer orbits have not been shown to constitute a positive distribution function or a stable stellar system. No orbital prescription is adopted as a physical solution.

## Formula provenance

Use the smooth anisotropy profile

\[
\beta(r)=\beta_0+(\beta_\infty-\beta_0)\frac{r^2}{r^2+R_e^2}.
\]

This is a fixed-sharpness case of established generalized anisotropy profiles, not a new companion formula. See [Baes and Van Hese, Dynamical models with a general anisotropy profile (2007)](https://arxiv.org/abs/0705.4109). The present choice fixes the transition at the observed projected half-light radius; that equality is a modeling assumption, not an empirical orbital measurement.

For radial pressure p=nu sigma_r^2, the standard spherical Jeans equation is

\[
\frac{dp}{dr}+\frac{2\beta(r)}r p=-\nu(r)g(r).
\]

Its integrating factor, up to an irrelevant constant, is

\[
F(r)=\left(\frac r a\right)^{2\beta_0}
\left(\frac{r^2+R_e^2}{a^2+R_e^2}\right)^{\beta_\infty-\beta_0},\qquad
p(r)=\frac1{F(r)}\int_r^\infty\nu(s)g(s)F(s)ds.
\]

These are established equation-solving steps applied to the chosen profile. Projection uses the existing PSF/annulus kernels with W-beta(r)T inside the radial integral. A variable beta must not simply be inserted as one constant outside that integral. When beta0=betainfinity, the existing constant-beta solution is recovered.

## Controlled setup

Keep the adjusted-reference companion profile, its galaxy-trained amplitude, constant stellar mass-to-light ratio, light deprojection, conditional geometry, and exact-lens-required stellar mass fixed. Do not combine the previous stellar M/L gradient with this orbital change. Fit beta0 and betainfinity using only inner Vrms bins and the same covariance. The outer bin remains a covariance-conditioned diagnostic outside that fit.

The lens angle has already fixed stellar mass and is consumed calibration. It is unchanged by this orbital fit because the mass profiles are held fixed. This adds one orbital nuisance degree of freedom relative to constant beta, with no fitted transition scale or new capture parameters. Existing lens targets and partitions are already exposed, so reduced residuals are not independent validation of the theory.

## Results and follow-up

| Case | Inner chi-squared: Chabrier / Salpeter | Outer residual square sum: Chabrier / Salpeter |
|---|---:|---:|
| Free constant-beta motion fits | 45.81 / 45.27 | 55.01 / 53.42 |
| Exact-lens mass, constant beta | 611.48 / 613.61 | 136.23 / 134.17 |
| Exact-lens mass, radial beta, original bounds | 447.76 / 453.84 | 103.09 / 101.47 |
| Outer-bound follow-up combined with remaining fits | 106.88 / 105.62 | 115.71 / 114.32 |

The follow-up is adaptive: after three systems hit betainfinity=0.45, those systems alone are refitted with -2<=betainfinity<=0.95, retaining the previous central bounds. The other three retain their initial interior fits. The combined row describes that tested configuration, not a proof of the global optimum over all enlarged-bound fits or all anisotropy profiles. No bound is asserted to be an observationally measured limit.

For the Chabrier cases at the expanded boundary:

| System | beta0 | betainfinity | Inner chi-squared | Conditional outer standardized residual |
|---|---:|---:|---:|---:|
| J0037-0942 | 0.2066 | 0.9500 | 27.90 | 1.68 |
| J1204+0358 | -0.0460 | 0.9500 | 4.92 | 6.98 |
| J1402+6321 | 0.2447 | 0.9500 | 54.21 | 6.54 |

The Salpeter counterparts exhibit the same endpoint pattern and are recorded in full. The outer residuals show why a much better inner fit is not enough. These standardized residuals depend on the released covariance and model assumptions; they are not model-independent exclusion levels.

Betainfinity=0.95 describes a highly radial asymptotic velocity ellipsoid under the conventional beta definition. It is not the anisotropy at every observed radius: beta(Re) is the mean of the two endpoints. Positivity of the projected second moments is necessary but insufficient for a physical phase-space distribution. Formation, orbital stability, nonspherical structure and independent kinematic constraints remain open.

## Interpretation

This result changes the diagnosis: a restrictive constant/weakly radial orbital model accounts for a substantial part of the earlier exact-lens tension. It does not justify declaring the companion model successful. The fit still depends on a new boundary, outer motions remain discrepant, and the same conditional optical mapping and stellar mass profile are retained.

The next physical discriminator is whether the demanded orbit structure can be supported by the tracer density and total potential, and whether a constrained orbital model can match both inner and outer data. Changing transition scales, combining gradients or altering geometry would add further choices and must be assessed explicitly. These orbital variations are ordinary stellar dynamics; they do not derive photon conversion, capture states or the one-third retention mechanism.

## Verification and reproduction

Run `python research_work/results/companion-extensions/radial-orbits.py`; see the [protocol](radial-orbits-protocol.md). The solver uses direct Jeans integration and 17 starting points for the initial endpoint fits, then multiple starts for the three-system outer-bound follow-up under each population proxy. Constant-limit predictions reproduce the stored exact-lens constant-beta velocities to 1.14e-13 km/s and are additionally checked against the original solver at three fixed beta values. Nested initial fits do not worsen the constant-beta inner optimum. Predicted squared velocities remain positive and finite. The mass and lens constraints are inherited without alteration.

Input hashes, complete per-system predictions, bounds and residuals are saved in `radial-orbits-results.json`. This verifies the variable-beta implementation on the inherited grid, not a full distribution-function construction, numerical refinement of every orbit model or dynamical stability. No source data or reference parameters are changed.
