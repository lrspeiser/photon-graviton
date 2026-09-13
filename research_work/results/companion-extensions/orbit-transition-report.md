# Free stellar orbital transition radius with fixed companion gravity

13 September 2026. Four-parameter stellar nuisance fit to exposed motion data and calibrated lens angles. No companion-law change.

## Outcome

Freeing the orbital transition radius improves the all-motion chi-squared by approximately 30%: 161.81/159.83 becomes 113.56/111.27 (Chabrier/Salpeter). The previously fixed half-light transition was materially restrictive. The best radii lie between 0.303 and 1.333 projected half-light radii, away from the imposed 0.1 and 10 limits.

This is not a joint physical solution. J1402+6321 retains all-motion chi-squared 64.97/63.24, accounting for approximately 57% of the remaining total. Its conditional outer residual is 5.47/5.40, slightly worse than in the fixed-radius fit. It also reaches the lowest allowed stellar gradient. Several other fitted stellar parameters remain on boundaries. All twelve configurations pass the imposed necessary orbital checks, but a positive distribution function and stability remain unestablished.

The model fits every motion bin and uses the catalogue lens angle to set stellar mass. The improvement is descriptive fit flexibility, not independent prediction or evidence that photon energy causes the gravitational component.

## Formula and provenance

Replace the fixed stellar orbit transition radius Re with ra:

    beta(r) = beta0 + (beta_infinity-beta0) r^2/(r^2+ra^2).

This remains the known smooth anisotropy family used previously; see [Baes and Van Hese](https://arxiv.org/abs/0705.4109). The new choice is to fit ra rather than equate it to the observed projected half-light radius. The added scale is a stellar orbital nuisance parameter, not a new photon-companion constant.

The radial Jeans integrating factor is changed consistently:

    F(r) = (r/a)^(2 beta0)
           * [(r^2+ra^2)/(a^2+ra^2)]^(beta_infinity-beta0).

Here a is the existing numerical reference length and cancels from physical predictions; ra is the orbital transition scale. The factor follows by integrating d ln F/d ln r=2 beta, a standard linear-equation operation. The projection kernel uses the same beta(r). Only changing the projection, without changing the radial pressure integral, would be inconsistent.

The stellar mass-to-light gradient still transitions at Re. The companion density, exact-third response, shared amplitude, capture scale, stellar light profile, PSF, covariance and conditional optical geometry remain unchanged. Lens bending calibrates total stellar mass for each gradient, just as in the preceding diagnostic.

## Constraints and optimization

Fit log(ra/Re) with 0.1<=ra/Re<=10; retain h in [-0.8,9], beta0 in [-2,0.45], beta_infinity in [-2,0.95]. The tighter central cusp bound beta0<=0.375 remains active. The gamma>=2 beta cap is recalculated at each trial ra on 4097 log radii and checked on 8193 final radii. For the specified spherical separable augmented-density class this is a known necessary condition, not a sufficient or universal distribution-function test. The central bound assumes the adopted light profiles continue to zero radius.

There are now four fitted stellar nuisance parameters per galaxy, versus three in the fixed-radius diagnostic. Forty motion bins across six galaxies are used under each population proxy; the population proxies are alternatives, not independent observations. No formal model-comparison probability is assigned using naive parameter counts, since boundaries, nonlinear fits, conditional geometry and prior data exposure matter.

Use 27 optimizer starts per case: the 22 prior starts at ra=Re, plus the prior fixed-radius best fit at ra/Re=0.1, 1/3, 1, 3 and 10. Twenty-six or twenty-seven succeed in each case. The lowest successful value is recorded; this does not prove a global optimum over this family or more general orbit models.

## Results

| Population | Fixed-radius total chi-squared | Free-radius total | Inner contribution | Conditional outer contribution |
|---|---:|---:|---:|---:|
| Chabrier | 161.812 | 113.565 | 61.696 | 51.868 |
| Salpeter | 159.834 | 111.274 | 60.435 | 50.839 |

| Galaxy | Population | ra/Re | h | beta0 | beta_infinity | Total chi-squared | Conditional outer residual |
|---|---|---:|---:|---:|---:|---:|---:|
| J0037-0942 | Chabrier | 0.5305 | -0.6573 | 0.3750 | 0.9500 | 10.444 | 1.759 |
| J0037-0942 | Salpeter | 0.5302 | -0.6483 | 0.3750 | 0.9500 | 9.739 | 1.658 |
| J1112+0826 | Chabrier | 0.4238 | 0.2607 | 0.3750 | -2.0000 | 18.668 | 3.214 |
| J1112+0826 | Salpeter | 0.4221 | 0.2932 | 0.3750 | -2.0000 | 18.700 | 3.215 |
| J1204+0358 | Chabrier | 1.3326 | 9.0000 | -1.3236 | 0.9500 | 14.333 | 2.825 |
| J1204+0358 | Salpeter | 1.3281 | 9.0000 | -1.2934 | 0.9500 | 14.424 | 2.835 |
| J1402+6321 | Chabrier | 0.3071 | -0.8000 | 0.3750 | 0.8635 | 64.974 | 5.465 |
| J1402+6321 | Salpeter | 0.3034 | -0.8000 | 0.3750 | 0.8636 | 63.243 | 5.395 |
| J1621+3931 | Chabrier | 0.7205 | -0.2547 | 0.0399 | 0.3373 | 3.409 | -0.545 |
| J1621+3931 | Salpeter | 0.7188 | -0.2192 | 0.0304 | 0.3335 | 3.426 | -0.549 |
| J1630+4520 | Chabrier | 0.4924 | 1.8061 | 0.3443 | -2.0000 | 1.737 | 0.544 |
| J1630+4520 | Salpeter | 0.4883 | 1.9967 | 0.3477 | -2.0000 | 1.742 | 0.550 |

The aggregate improvement is dominated by inner motion: the inner contribution drops from 104.25/103.00 to 61.70/60.44, while the conditional outer sum drops from 57.57/56.84 to 51.87/50.84. The already-consumed outer bins must not be called withheld predictions.

J0037-0942 improves markedly, while J1630+4520 reaches a low descriptive residual. J1621+3931 changes little. These cases do not justify omitting the poor fits: J1112+0826, J1204+0358 and especially J1402+6321 retain discrepancies. The listed standardized residuals are conditional on the model and released covariance, not model-independent significance measures.

J1204+0358 retains h=9; J1402+6321 reaches h=-0.8, implying a central stellar mass-to-light ratio only 0.2 times the outer value. Their physical population interpretation is not independently verified. J0037-0942 and J1204+0358 reach the radial outer endpoint 0.95; J1112+0826 and J1630+4520 reach the tangential endpoint -2. Three systems reach the central cusp bound under both populations. No transition-radius boundary is active.

## Verification

At ra=Re the new implementation reproduces the stored all-motion objective exactly at the recorded precision in all twelve cases. Every final score is no worse than that feasible reference. The older radial-only and gradient-only inner-limit checks remain active. The block-covariance decomposition of total chi-squared into inner and conditional-outer contributions agrees to within 2.85e-14. The imposed lens equation closes below 1e-10 fractional error, and masses and predicted moments are positive.

All refined slope margins are positive, with a minimum of 0.0777244; the central asymptotic margin reaches zero in the boundary cases. The known necessary check does not construct a distribution function, establish dynamical stability, or validate the extrapolated central light profile. The inherited numerical integration grid and smooth orbital law are retained; no raw-image lens reconstruction or new observational calibration is performed.

## Consequence for the theory

A fixed orbit transition at Re was an avoidable nuisance-model restriction. Relaxing it is useful evidence about the stellar modeling, but it does not explain the remaining lens-calibrated motion profile. The dominant system still requires a better explanation under the adopted assumptions. Adding more independently fitted stellar scales could continue lowering residuals while weakening predictive content; that alone would not derive the companion mechanism.

Further options should be distinguished by evidence: independently constrained stellar gradients, nonspherical orbit structure, a positive distribution-function construction in the fitted potential, or a change in the gravitational/optical response. Any comparison with MOND or dark-halo models would need the same data and stellar freedoms. The photon-transfer, absolute supply, retention, event-duration and emission tracks remain unresolved and are not established by these orbital fits.

## Reproduction

Run `python research_work/results/companion-extensions/gradient-orbits.py --slope-constrained --all-motion-bins --free-orbit-radius`. Results are written to orbit-transition-results.json, preserving earlier outputs. The file records every predicted and observed motion bin, radius and other parameters, objective values from successful starts, failure messages, bound checks, slope margins, reproduction errors and input hashes. Pre-execution choices are in orbit-transition-protocol.md.
