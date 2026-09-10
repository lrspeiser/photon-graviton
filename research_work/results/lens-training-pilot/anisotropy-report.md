# How much of the lens discrepancy depends on stellar orbits?

Favoring radial stellar motions reduces the companion pilot's training angle overprediction, while favoring tangential motions increases it. Under the same orbital assumption for both models, neither tested alternative makes the companion pilot outperform the ordinary-matter benchmark in training RMS. The radial case brings their discrepancies much closer. This identifies an important modeling uncertainty; it does not establish a successful revised theory.

These are post-validation development calculations on the original 33 training systems only. The earlier isotropic validation result is preserved. Nonzero-anisotropy execution is explicitly restricted to training by the command. No test predictions were evaluated and no validation result is relabeled as passing.

## What changed

**Known stellar-dynamics definition, not a new gravity parameter:**

\[
\beta=1-\frac{\sigma_\theta^2+\sigma_\phi^2}{2\sigma_r^2}.
\]

Beta is zero for equal radial and tangential dispersions, positive for radial preference, and negative for tangential preference. We test the common constants −0.3 and +0.3 against the archived zero case. These are illustrative alternatives, not measured values, an inferred population distribution, or a confidence interval. Each value applies equally to every training galaxy and to both gravity models.

Changing beta changes how a gravitational field translates into the aperture dispersion. The mass is therefore re-inferred from each galaxy's measured dispersion under that orbital assumption, then the Einstein angle is predicted. The angle never enters the mass root. A, p, a*, distance geometry, light-profile sizes, seeing cases and companion boundary cases remain fixed.

## Training outcomes

Representative cases use the already stated 1.5-arcsec Gaussian seeing and 20-effective-radius companion source boundary.

| Shared orbital beta | Ordinary-matter angle RMS | Companion angle RMS | Companion median predicted / SIE angle |
|---|---:|---:|---:|
| −0.3, tangential preference | 0.1872 arcsec | 0.2294 arcsec | 1.1319 |
| 0, isotropic | 0.1840 arcsec | 0.2083 arcsec | 1.0973 |
| +0.3, radial preference | 0.1904 arcsec | 0.1942 arcsec | 1.0417 |

The radial case reduces the median companion overprediction from 9.7% to 4.2%; it does not erase the object-to-object discrepancies. The ordinary-matter model's median predicted/SIE ratio becomes 0.9789 in that case. We do not choose a different beta separately for each lens or for the two models to manufacture a preferred comparison.

The normalization-only stellar-population comparison also shifts. With event stretching S=1+z and the fixed population assumptions from the photometric audit, the companion pilot's median required ordinary mass divided by the population estimate is:

| Beta | Chabrier IMF | Salpeter IMF |
|---|---:|---:|
| −0.3 | 2.52 | 1.42 |
| 0 | 2.37 | 1.35 |
| +0.3 | 2.20 | 1.25 |

There are 32 mass-matched training galaxies. These ratios are not full stellar-population posteriors or significance tests. The radial alternative helps but does not remove the mass-normalization issue under those assumptions. Source-population priors, photometric errors and the two IMF alternatives remain as previously documented.

## Generalized operator

**Known constant-anisotropy spherical Jeans equation and solution:**

\[
\frac{d(j\sigma_r^2)}{dr}+\frac{2\beta}{r}j\sigma_r^2=-jg,
\qquad
j\sigma_r^2=r^{-2\beta}\int_r^\infty j(s)g(s)s^{2\beta}\,ds.
\]

The solution assumes the stated vanishing outer-pressure boundary. At a viewing angle theta relative to the radial direction, the line-of-sight variance is sigma_r²(1−beta sin²theta). With aperture/seeing inclusion probability P(r,theta), define

\[
W(r)=\int_0^1P(r,\mu)\,d\mu,\qquad
W_\beta(r)=\int_0^1P(r,\mu)[1-\beta(1-\mu^2)]\,d\mu,
\]

where mu=cos(theta). **Known projection and integration interchange** then give

\[
\sigma_{ap}^2=
\frac{\int_0^\infty j(s)g(s)s^{2\beta}K_\beta(s)\,ds}
{\int_0^\infty j(r)r^2W(r)\,dr},
\qquad K_\beta(s)=\int_0^s r^{2-2\beta}W_\beta(r)\,dr.
\]

For a sharp circular aperture with mu_min=sqrt(max(0,1−R_ap²/r²)), W=1−mu_min and

\[
W_\beta=W-\frac{\beta}{3}W^2(\mu_{min}+2).
\]

With Gaussian seeing the orientation integral is evaluated directly using the existing aperture-inclusion kernel. The luminosity denominator uses W, not W_beta; replacing it by W_beta would incorrectly change the observed light normalization.

The beta=0 branch retains the prior isotropic operator. All equations here are established moment/projection mathematics. Beta describes stellar orbits, not the production or gravity coupling of companion waves.

## Physical and numerical checks

The Hernquist tracer's central cusp slope is one. The tested beta values satisfy the applicable necessary central cusp condition gamma>=2beta, discussed in [An & Evans, A Cusp Slope–Central Anisotropy Theorem](https://arxiv.org/abs/astro-ph/0511686). This is only a necessary condition; it is not a proof that every proposed potential and beta combination has a nonnegative distribution function or stable supported population. Jeans moments alone cannot supply that proof.

Independent integration of Jeans pressure followed by aperture projection agrees with the swapped-integral calculation to better than 1.6×10⁻⁷ relative in the verification cases. The angular projection identity is checked separately. A complete spherical aperture recovers the known Hernquist self-gravitating result sigma_los²=GM/(18a) for all three beta values, as required by spherical averaging and the virial relation. Refining the radial and angular numerical grids changes any new predicted lens angle by less than 0.000018 arcsec. These checks verify the operator, not the astrophysical assumptions.

## Decision and reproducibility

Do not treat the reduced radial-orbit bias as a measured preference for beta=0.3 or as a fix to the failed validation pilot. Before promoting an orbital revision, obtain independent kinematic/profile constraints or a documented population model for beta, verify physically admissible distribution functions, and use the same treatment under both gravity hypotheses. The current lensing summaries cannot independently determine every mass, orbit shape, light-loss factor and companion response.

The broader requirements remain: a derived common propagation/arrival-time law, energy accounting, a capture and support mechanism, a field predicting motions and lensing together, and tests across the other observations. Refining stellar nuisance assumptions cannot fill those missing physical equations.

Reproduce the two new runs using run.py --updated-profile --refined --beta 0.3 and --beta -0.3, repeat without --refined for numerical comparisons, then run verify_anisotropy.py. All eight existing seeing/boundary configurations are preserved for each beta, totaling 528 new training configurations at each numerical resolution. The mass-ratio summary reuses the archived fixed-population normalization estimates; no stellar-population refit is claimed.
