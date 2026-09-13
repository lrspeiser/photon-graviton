# Lensing uncertainty from the same stellar-motion mass posterior

Only inner stellar-motion measurements determine the mass and constant orbit-anisotropy posterior. The outer stellar bin and the published lens angle are excluded from its likelihood. The frozen force and light profiles then map that mass uncertainty into a lens-angle interval.

These six galaxies are exposed training systems. The intervals below are central 95% conditional credible intervals for the model angle, not measurement-inclusive prediction intervals or validated coverage.

| Galaxy | Published SIE angle | Baryonic angle interval | Extra-force angle interval |
|---|---:|---:|---:|
| J0037-0942 | 1.530 | 1.463 to 1.526 | 1.530 to 1.608 |
| J1112+0826 | 1.490 | 1.139 to 1.210 | 1.194 to 1.280 |
| J1204+0358 | 1.310 | 1.268 to 1.343 | 1.326 to 1.419 |
| J1402+6321 | 1.350 | 1.352 to 1.405 | 1.375 to 1.438 |
| J1621+3931 | 1.290 | 1.080 to 1.173 | 1.070 to 1.193 |
| J1630+4520 | 1.780 | 1.327 to 1.403 | 1.394 to 1.487 |

All angles are in arcseconds. Reference prior: uniform log mass and uniform anisotropy. A catalog SIE angle is a fitted nonspherical lens summary; it is not an exact measured spherical Einstein radius.

The catalog summary falls outside the displayed model interval in 5/6 baryonic cases and 5/6 extra-force cases. This descriptive count is not a hypothesis rejection probability. It omits lens-summary uncertainty and the structural assumptions below.

Changing to a uniform mass prior changes the median model angle by at most 0.00120 arcseconds.

## Method and formula provenance

Bayesian marginalization, Jeans projection and the weak-field lens integral are established mathematics. The additional force remains a frozen empirical prescription, not a photon-derived deposit source.

For mass M and orbit anisotropy beta, the Gaussian inner-bin likelihood uses the full released covariance. Its posterior is integrated over beta to obtain p(M | inner motions). The same ComponentModel.angle(M) used in the earlier fixed-mass lens calculations supplies the angle. No lensing normalization is fitted. Mass quantiles at 2.5, 16, 50, 84 and 97.5 percent are mapped through that increasing lens-angle function; their ordering is checked at both numerical resolutions.

The prescribed extra force scales with M to the fixed positive power p. Both force contributions increase with M at fixed radius. The adopted single Einstein-root branch gives the increasing mapping used here.

## Numerical verification

The initial 801/1601 mass-node and 321/641 anisotropy-node comparison failed 15 of 24 interval checks. Its results remain saved. The repeated 3201/6401 mass-node and 641/1281 anisotropy-node comparison passes all 24 unchanged gates: at most 0.005 in log-mass quantiles and 0.005 arcseconds in angle quantiles. The largest refined angle change is 0.000681 arcseconds.

These checks establish integration stability at the stated tolerance, not physical or observational agreement. Both runs retain dependency hashes, checked by this report generator.

## Remaining assumptions and consequence

Light components, their spherical deprojection, constant mass-to-light ratio, seeing, covariance, static redshift-derived geometry, equality of the two metric potentials, force coefficients and outer cutoff are fixed. Their uncertainties are absent. The same omissions affect the inner stellar inference and can correlate its errors with lensing.

A genuine joint test needs image/shear measurements with uncertainties and a consistently modeled stellar distribution. These intervals cannot substitute for that likelihood. They expose which discrepancies persist when mass and orbit uncertainty is included, without allowing a separate lensing mass fit. The companion production, capture, and transport equations still have to supply the gravitational source.

Reproduce: run `run.py`, `refine.py`, then `summarize.py` from this folder or by repository-relative script paths. No reserved outcomes opened.
