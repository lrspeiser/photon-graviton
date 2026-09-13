# Necessary phase-space consistency check for the fitted orbital profiles

13 September 2026. No new fit or change to the companion reference.

## Outcome

No sampled violation of gamma(r)>=2 beta(r) occurs in the 12 original-bound radial-orbit fits or the six expanded-outer-bound fits. The smallest margin is 0.32422 in the first set and 0.32546 in the second, over 1e-6<=r/Re<=100. The strongly radial fits are therefore not rejected by this necessary test for the specified distribution-function construction.

This is not a proof that a positive distribution function exists, a stability result, or an improvement to the observed-motion residuals. In particular, the outer stellar discrepancies recorded in the radial-orbit experiment remain unchanged.

## Established theorem and applicable class

For a spherical stellar system whose augmented density separates into a function of relative potential and a function of radius, with central anisotropy beta0<=1/2, positivity implies the global density-slope/anisotropy condition gamma>=2 beta. This result is established in [Van Hese, Baes and Dejonghe, On the universality of the global slope--anisotropy inequality](https://arxiv.org/abs/1010.4301). It is not a universal criterion for arbitrary distribution functions, and it is not sufficient for positivity.

Here gamma=-d ln nu/d ln r refers to the tracer density, not the total companion-plus-stellar gravitating density. The checked beta0 values satisfy the theorem's central-anisotropy condition. We consider separable augmented density as a possible completion of the previously fitted Jeans moments; separability was not established by those fits. The theorem is known stellar dynamics, not a new photon-companion law.

## Relation to the prior orbital model

For beta(r)=beta0+(betainfinity-beta0)r^2/(r^2+Re^2), the prior Jeans integrating factor F obeys d ln F/d ln r=2 beta. A candidate separable construction takes its radial factor proportional to 1/F and its potential-dependent factor along the physical profile proportional to nu F. Thus

\[
\frac{d\ln(\nu F)}{d\ln r}=-\gamma+2\beta.
\]

For a monotonically decreasing relative potential, a nonnegative gamma-2 beta corresponds to the required first-derivative sign of that potential-dependent factor. This algebra explains what is being checked; higher conditions and the actual distribution-function inversion remain unperformed. The total gravitational potential enters those further steps even though this necessary inequality only uses tracer slope and anisotropy.

## Density-slope calculation

Use the same fixed multi-component Sersic surface-light models as the lens/Jeans work. For one component, write R=r cosh(u) in the Abel integral. Its deprojected density is the integral of -Sigma'(R)/pi over u. Differentiating the integrand gives the density-weighted slope factor

\[
1-\frac1n+\frac{b_n}{n}\left(\frac R{R_e}\right)^{1/n}.
\]

Sum the density and slope numerators across components, then divide to obtain gamma(r). This is differentiation of the known Abel/Sersic formula, not a new empirical slope law. Component radii and amplitudes are retained; the effective half-light radius only sets the dimensionless coordinate. The same exponentially suppressed integration cutoff used by the preceding deprojection is retained. No dark-matter or expanding-universe model is used to supply the tracer slope.

## Executed checks

| Orbital set | Population/system cases | Sampled violations | Minimum gamma-2 beta |
|---|---:|---:|---:|
| Original endpoint bounds | 12 | 0 | 0.32422 |
| Expanded outer-bound cases | 6 | 0 | 0.32546 |

The radial grid has 4097 logarithmically spaced points from 1e-6 to 100 Re. These are checks on extrapolated fitted light profiles as well as measured-aperture scales; they are not new observations over that full range. A finite grid is not an all-radii mathematical proof.

Comparing deprojection quadrature orders 256 and 512 changes gamma by at most 2.14e-14. Independently differentiating the resulting log-density on the radial grid agrees with the integrated slope to within 7.27e-5 away from grid edges. The positive margins exceed these numerical differences substantially. These checks validate the slope evaluation for the adopted profiles, not their astrophysical accuracy.

## Consequence and remaining work

The large outer-anisotropy fits cannot be dismissed using this necessary condition. They remain candidates for a more demanding distribution-function calculation. They also remain boundary-dependent fits with poor conditional outer motions, so passing this check is not evidence that they describe real systems.

A subsequent completion needs a nonnegative phase-space distribution reproducing the tracer density in the specified potential, followed by dynamical/support assessment. Nonseparable alternatives are not excluded here. Radial-orbit instability, source formation, nonspherical structure and the conditional lens geometry remain outside this calculation. No global research goal or physical viability claim is closed by this test.

## Reproduction

Run `python research_work/results/companion-extensions/orbit-slope-check.py`. It records hashes of the orbital results and light-profile inputs, the theorem reference, per-case minimum margins and radii, and numerical derivative/quadrature checks in `orbit-slope-check-results.json`. It modifies no data, fitted orbit, mass profile or lens prediction.

Related: [radial-orbit fit and its limitations](radial-orbits-report.md).
