# Modest tangential support is insufficient for three fitted lens profiles

The previous support calculation found that five of six morphology-mapped lens profiles cannot have isotropic bound-particle support. This run tests a fixed alternative without changing any deposited density, fitted ordinary stellar mass, or predicted lens angle: constant deposit anisotropy beta_d=-1/2.

## Physical meaning and known formula

Beta=1-(sigma_theta^2+sigma_phi^2)/(2 sigma_r^2). At beta=-1/2, summed tangential variance is three times radial variance. This describes more sideways than radial motion, with no required net rotation. It concerns deposited particles; it is not the stellar anisotropy previously fitted to spectroscopy.

Use the known constant-anisotropy distribution family f(E,L)=L F(E), where L is specific angular-momentum magnitude and E is relative binding energy. Direct velocity-space integration gives

    rho(r)/r=2pi^2 integral_0^Psi F(E)(Psi-E)dE,
    F(E)=1/(2pi^2) d^2[rho/r]/dPsi^2.

These are established distribution-function mathematics, not a new capture or gravity law. The general consistency framework is discussed by [An et al.](https://arxiv.org/abs/1202.0004). Applying this family is an explicit support hypothesis, not a derivation of how incoming companions acquire their orbits.

For h=rho/r, g=-dPsi/dr, l=dln h/dln r and m=dln g/dln r, the sign of F is the sign of

    Q=l^2-l+dl/dln r-l*m,
    d^2h/dPsi^2=h Q/(r^2 g^2).

The total g includes both the unchanged deposits and the ordinary mass obtained from the same deprojected Sersic components. A negative Q implies negative F and rejects this particular bound-particle equilibrium. It is not merely a poor least-squares fit.

## Results at fixed fitted galaxy profiles

| Galaxy | Minimum Q on tested range | Radius of minimum, kpc | Negative region resolved? |
|---|---:|---:|---|
| J0037-0942 | -0.60293 | 19.39 | Yes |
| J1112+0826 | 0.18726 | 11.02 | No |
| J1204+0358 | 0.64568 | 4.41 | No |
| J1402+6321 | -0.24891 | 15.93 | Yes |
| J1621+3931 | -0.19526 | 16.13 | Yes |
| J1630+4520 | 0.00682 | 13.20 | No |

For the three failing cases, negative regions extend approximately 12.81-24.90, 12.30-19.10, and 12.78-19.01 kpc, respectively. They persist when radial sampling is doubled; the largest Q change is below 0.000775, much smaller than those negative minima. The tested range is r/a_capture from 1e-3 to 1e3, away from integration endpoints. These are conditional model radii, not measured particle trajectories.

The other three only pass the sampled inversion-sign test. In particular the small positive minimum for J1630+4520 is conditional on the fixed input mass/light geometry. No full uncertainty propagation, global positive-distribution proof, dynamical stability or formation calculation is provided.

## Consequence

A single modest tangential bias does not support all the fitted profiles. More tangential or radially varying distributions, redistribution after capture, or a physically specified field stress remain possible, but none has been supplied by this test. Selecting a more extreme anisotropy independently for every galaxy would introduce additional freedom requiring its own physical explanation and prediction test.

Rotation and lensing outputs remain the earlier conditional results because their mass profiles are unchanged. The failure affects whether this particular particle population can realize those profiles. It neither establishes the companion theory nor rules out all bound states in the hypothetical universe. All six goals remain open.

Reproduce: `python research_work/results/isotropic-galaxy-transfer/tangential.py`. tangential-results.json records both resolutions, negative intervals, mass-integration checks and all six outcomes. The morphology and geometry assumptions of lensing-report.md remain in force.
