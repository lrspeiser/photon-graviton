# The Milky Way bar changes the nonlinear gravity source before any orbit fit

**The circularly averaged approximation is close at the tested locations, but it removes position-dependent structure relevant to the bulge test.** Applying the frozen nonlinear response before averaging the bar changes the mean response vector by up to 1.183% on 32 rings. Its azimuthally averaged effective Poisson density changes by up to 1.853% on 20 derivative probes. The response vector itself varies around the inner rings by as much as 21.057% RMS relative to the response of the mean field.

These are calculations in a specified ordinary-matter model, not measured deviations of stars, predicted changes in their speeds, or evidence for photon-produced gravity. They identify a concrete approximation to address before the full bulge orbit comparison. No training or held-out stellar velocities were loaded, and no fitted parameter was changed.

## Shared gravity rule and provenance

The [conservative field completion](../conservative-field-completion/report.md) uses a known QUMOND-style Poisson construction with the project's previously fitted response coefficient. For ordinary-matter acceleration g_b, define

`B = A*(|g_b|/a_star)^(p-1)*g_b`,

`div(g_extra)=div(B)`, `g_extra=-grad(Phi_extra)`.

Equivalently, `Laplacian(Phi_extra)=-div(B)`. This separates the conservative acceleration from the algebraic vector B. **B is not generally the acceleration to apply to stellar orbits.** Solving the Poisson equation with appropriate boundary conditions is still required.

**Provenance:** the mathematical construction is established in [Milgrom's QUMOND formulation](https://arxiv.org/abs/0911.5464). The power-law response and fixed A, p and a_star here are the project's empirical fit, not an independently derived photon/companion interaction or a new gravity theorem. They remain unchanged from the joint-galaxy training result: A=0.2422960666, p=0.4624587420 and a_star=7.2496089687e-10 m/s^2. The new work is a numerical audit of the averaging approximation.

Because B depends nonlinearly on g_b, taking the average first is not generally equivalent to taking it afterward:

`mean_phi[B(g_b)] != B(mean_phi[g_b])`.

This is a known property of nonlinear functions, not a novel physical law. All averages here use cylindrical radial, tangential and vertical components. Averaging Cartesian vectors around a ring would incorrectly erase its radial field.

## Ordinary matter and test locations

We use the existing conservative order-24 disk-potential spline, bar spherical harmonics through order 64, nuclear components through order 16, and the same central-mass approximation from the [rotating-bar calculation](../rotating-bar-orbits/report.md). No dark halo or synthetic companion ring/cap/shell is added. Published mass-model assumptions and uncertainties remain; these are not independent mass measurements.

This keeps the two averaging orders on exactly the same baryonic numerical field. It is not an exact rerun of the separate axisymmetric completion's higher-order disk representation. The field is evaluated in the bar frame; rotating the bar rigidly leaves each complete-ring average unchanged. A comparison with actual stars must also specify the Sun-to-bar orientation.

The 32 rings combine radii 0.5, 1, 2, 3, 5, 8, 12 and 20 kpc with heights 0, 0.1, 0.8 and 1.5 kpc. Each is sampled at 128 and then 256 equally spaced azimuths. Negative heights are symmetry-related in this adopted potential; this does not establish north/south symmetry of the observed stellar populations.

## Results relevant to the plane and bulge-height comparison

The mean difference is the norm of the two mean-vector prescriptions' difference divided by the norm of B(mean g_b). Azimuthal RMS is the vector RMS about mean B, using the same denominator. Neither column is a stellar-velocity percentage.

| Radius (kpc) | Height (kpc) | Mean response difference | Azimuthal response RMS |
|---:|---:|---:|---:|
| 1 | 0.1 | 0.817% | 17.259% |
| 1 | 0.8 | 0.419% | 12.216% |
| 3 | 0.1 | 1.167% | 20.478% |
| 3 | 0.8 | 0.655% | 13.352% |
| 8 | 0.1 | 0.0118% | 1.743% |
| 8 | 0.8 | 0.0107% | 1.646% |

The larger response variation near the inner plane is a prediction of this prescribed bar plus nonlinear response, before solving for the added potential. It must not be presented as observed evidence that photons are intercepted there, that stars speed up by those percentages, or that the deposited field has the same variation. Spatial differentiation and the nonlocal Poisson solution can change its magnitude and shape.

## Check the actual averaged Poisson source

To go beyond the vector comparison, use the known cylindrical divergence identity. The full-ring average of the azimuthal derivative vanishes, leaving

`mean_phi[div(B)] = (1/R)*d_R[R*mean_phi(B_R)] + d_z[mean_phi(B_z)]`.

The effective source density is `rho_extra=-div(B)/(4*pi*G)`. It is a mathematical density that generates the proposed additional potential, not a measured deposited-energy density or proof of a physical particle population.

We differentiate both averaging orders at 20 points: R=0.5,1,3,8,20 kpc and the four heights. Central-difference steps are 0.002 and 0.001 kpc, with 128/256 azimuth refinement. The largest fractional source difference occurs at R=3 kpc, z=0.8 kpc:

- Response before averaging: 3.1438851e7 solar masses/kpc^3 effective density.
- Averaging before response: 3.0856389e7 solar masses/kpc^3.
- Difference: 1.853%, normalized by the larger absolute density.

This is a local source comparison. A fractional source change is not a bound on the resulting gravitational-force error: its full spatial distribution and boundary conditions are needed. The calculation does not solve the three-dimensional potential or determine a statistically acceptable modeling tolerance from data.

## Numerical controls and practical consequence

Doubling azimuth samples changes the mean-difference/RMS statistics by at most 2.09e-8 in fractional units. Replacing bar order 64 with 40 changes them by at most 3.19e-5 and 3.88e-4 respectively; that finite comparison is not a continuum guarantee. The effective-density derivative-step change is at most 0.0681%, and its azimuth refinement change is at most 2.67e-9 fraction. These numerical errors remain separate from astrophysical mass-model uncertainty and from the uncomputed potential.

An exact m=0 bar projection reproduces the averaged ordinary-matter field at three independent controls. Its response commutes with averaging. A deliberately linear response also commutes in all ring controls. These checks distinguish the nonlinear bar effect from a coordinate-average mistake. The existing galpy setup emits its previously documented zero-radius helper and unavailable C-extension warnings; all sampled fields and reported diagnostics are finite, and the calculation uses the Python implementation and existing disk spline.

The next field calculation should retain the bar-dependent source and solve one conservative three-dimensional potential before fitting the position-conditioned stellar velocity distributions. The existing axisymmetric result remains a declared approximation for circularly averaged diagnostics. This audit supplies no reason to retune a separate speed multiplier for the plane, bulge top or bulge bottom. In particular, the previously explored synthetic deposition rings are not a substitute for the frozen empirical response tested here.

The companion interpretation still needs its independent source/capture law. Even a successful full-bar gravity fit would not establish photon conversion, observable redshift/event timing, lossless transport or the deferred total photon-supply requirement.

## Reproduction

Run `python research_work/results/bar-response-averaging/run.py`. `results.json` contains all ring values, density-source comparisons, numerical sensitivities, axisymmetric/linear controls, response parameters and hashes of the exact field caches and code. No holdout outcomes were opened. The full stellar orbit likelihood remains unfinished.
