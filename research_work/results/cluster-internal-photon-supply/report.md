# Internal stellar photons: conversion before cluster escape

## Outcome

Efficient companion capture does not mean efficient conversion of stellar light inside a cluster. We now connect photon energy to companions and deposits along the same finite path, instead of assuming a pre-existing companion bath. The test preserves every unit as escaping photons, escaping companions or retained deposits.

For a stipulated R=1 Mpc sphere and the archived path-loss coefficient alpha=0.0002488993265191759 per Mpc, uniform internal emission converts only 0.018665% of its photon energy before the boundary. Even perfect companion capture cannot exceed that fraction in this branch. Strong finite capture kappa R=10 retains 0.016362% of emitted photon energy, or 87.66% of the energy converted before exit. These are different denominators from the previous external-bath capture fractions.

The radius is an illustrative receiver size, not an adopted universe size or a measured cluster radius. The archived alpha is an exposed-data empirical coefficient, not a first-principles constant; see ../joint-galaxy-audit/report.md. The code also evaluates alpha R=0.01 and 1 as declared sensitivity cases, not observationally acceptable new fits.

## Formula provenance and assumptions

The project postulates a uniform distribution of isotropic internal photon emission, constant conversion alpha, constant companion absorption kappa, common straight paths at c, and permanent local storage. There is no external inflow, photon recycling, focusing or evolving gravity. The source pulse is followed until all its unbound energy leaves the sphere, without imposing a cosmic age. This is a restricted transport calculation, not a full dynamical energy-momentum closure.

Known linear transfer equations, applied to the proposed sectors along distance s:

    dP/ds = -alpha P
    dC/ds = alpha P - kappa C
    dD/ds = kappa C
    (P,C,D)(0) = (1,0,0).

P, C and D are fractions of initial photon energy. D is cumulative energy deposited along the path, not a moving intensity. Established matrix-exponential mathematics gives

    P(s) = exp(-alpha s)
    C(s) = alpha [exp(-alpha s)-exp(-kappa s)]/(kappa-alpha)
    D(s) = 1-P(s)-C(s).

For kappa=alpha the finite limit is C(s)=alpha s exp(-alpha s). These are conditional solutions of stated postulates, not claimed novel physical laws. The density, stresses and support of the deposit are not established by this ledger.

For uniform emission in a sphere, known geometric integration gives the forward exit-length probability density

    p(L) = 3/(4R) [1-L^2/(4R^2)], 0 <= L <= 2R.

To derive it, fix the ray direction by spherical symmetry. At cylindrical radius b, the source is uniformly distributed along a chord; integrating the volume element 2 pi b db dL over chords longer than L gives p(L). Its mean is 3R/4. This differs from the chord distribution for rays entering from outside. Averaging P,C,D over p(L) therefore counts internal emission once.

As kappa tends to infinity, D tends to 1-exp(-alpha L). The mean deposited fraction is bounded above by the mean conversion fraction, approximately 3 alpha R/4 when alpha R is small. More stars increase emitted energy, but do not change this fraction when the assumed geometry and rates remain fixed.

## Executed values

All fractions in this table are relative to emitted photon energy after following the complete pulse.

| alpha R | kappa R | Escaping photons | Escaping companions | Deposited |
|---|---|---:|---:|---:|
| 0.00024889933 | 0.1 | 0.999813350 | 0.000177096 | 0.000009554 |
| 0.00024889933 | 1 | 0.999813350 | 0.000117649 | 0.000069000 |
| 0.00024889933 | 10 | 0.999813350 | 0.000023028 | 0.000163621 |
| 0.01 | 0.1 | 0.992539834 | 0.007077883 | 0.000382283 |
| 0.01 | 1 | 0.992539834 | 0.004699875 | 0.002760291 |
| 0.01 | 10 | 0.992539834 | 0.000918834 | 0.006541332 |
| 1 | 0.1 | 0.527252194 | 0.446207435 | 0.026540371 |
| 1 | 1 | 0.527252194 | 0.284759506 | 0.187988301 |
| 1 | 10 | 0.527252194 | 0.050291910 | 0.422455896 |

## Verification and next physical step

All nine scenarios conserve energy. The closed forms agree with an independent matrix exponential within 2.1e-15. Exit-length integration agrees with direct source-volume/direction integration within 8.9e-11; doubling spatial/angular resolution changes fractions by at most 1.4e-9. The equal-rate limit is explicitly exercised. No astronomical holdout was used.

This does not establish an overall cluster energy shortage: external photons can convert over much longer paths and arrive as companions, while total stellar emission histories and any initial reservoir remain uncalculated. It also does not show that retained energy reproduces a lens map. The next supply calculation must include external sources and their history with causal travel and intervening absorption, then spatial capture and gravitational response. Escaped internal energy cannot also be assigned to this cluster unless a return path is modeled. Age is not fixed; additional time increases supply only insofar as additional energy is actually emitted.

All six physical research objectives remain open. Run `python research_work/results/cluster-internal-photon-supply/run.py`.
