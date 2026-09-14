# Coma inverse calculation: finite positive shapes fit, formation remains unresolved

13 September 2026. The six saved Kubo weak-shear bins admit positive finite total-density profiles. A simple supported-profile candidate fits about as closely as the old NFW fit. No shared energy-conserving companion capture law has yet been demonstrated: the inverse profiles include all gravitating matter, and absolute companion supply is not established.

## What was fitted

Use all six previously reconstructed vector-figure measurements and plotted diagonal errors, including the negative outer shear. Work in radius ratios R/R0; no cosmological physical distance or published halo mass is imported. The matched ordinary gas/star profile and calibrated source geometry needed for companion-only inversion are absent from this saved comparison. Treating the entire fitted density as companion deposits would double-count ordinary matter.

Nonnegative sums of 25 known spherical Plummer profiles provide smooth, finite, positive density. For each component, with dimensionless central convergence w and scale a:

    kappa(R) = w a^4/(R^2+a^2)^2
    gamma_t(R) = w a^2 R^2/(R^2+a^2)^2
    M/(pi Sigma_crit R0^2) = w a^2

These are established projection formulas, not project inventions. Direct density integration independently verifies the surface-density expression. Constrain total kappa<=0.05 at measured radii. Fit weak shear as in the inherited Coma diagnostic; this is not a full source-weighted reduced-shear likelihood. The largest reduced-versus-weak correction among retained profiles is 0.186 plotted standard errors, small for this exploratory exercise but not zero.

## Fit results

| Profile | Six-bin diagonal chi-squared |
|---|---:|
| Earlier fitted NFW | 3.855 |
| Single finite Plummer, fitted scale and amplitude | 3.727 |
| Positive mixture, outer scale sensitivity 10 Rlast | 3.558 |
| Positive mixture, outer scale sensitivity 30 Rlast | 3.566 |
| Positive mixture, outer scale sensitivity 100 Rlast | 3.572 |

The single-component fit was added as an exploratory follow-up after the mixture envelopes. It uses two locally fitted parameters, not frozen galaxy constants. The mixtures have many more degrees of freedom; their smaller chi-squared is not model-selection evidence. Minor score differences between extent grids reflect a changed basis spacing, not a physical preference for one boundary. These results show compatible positive shapes, not superiority over NFW or a companion formation prediction.

## How much outer material is permitted?

Linear optimization finds minimum and maximum mass while requiring every shear prediction to lie within two plotted errors, positive weights and kappa<=0.05. This box of feasible predictions is not a statistical confidence region. The maximum component scale is an imposed sensitivity, not a measured collection boundary; Plummer tails extend beyond it.

| Largest component scale / last measured radius | Minimum total mass units | Maximum total mass units | Feasible enclosed mass at last radius |
|---|---:|---:|---:|
| 10 | 0.00914 | 318.58 | 0.00914 to 0.46275 |
| 30 | 0.00914 | 2865.66 | 0.00914 to 0.46157 |
| 100 | 0.00914 | 31838.64 | 0.00914 to 0.46166 |

Mass units are pi Sigma_crit R0^2, not solar masses. The enclosed extrema are separate optimizations, not necessarily achieved by the total-mass extrema. A very broad component resembles an almost uniform projected sheet across the measured region: it can add substantial mass with little shear. Therefore these six measurements do not provide a useful upper bound on total outer inventory without external constraints. This is freedom in an inverse fit, not evidence that an enormous companion supply exists. Positive distributions can fit the negative outer point within its plotted uncertainty; no negative density is introduced to match it exactly.

## A concrete equilibrium candidate

A known pressure/support law provides a finite-density profile rather than simply placing material arbitrarily:

    grad P = -rho grad Phi
    P = K rho^(6/5)
    Laplacian Phi = 4 pi G rho

In the isolated self-gravitating limit, these admit the Plummer (n=5 polytropic) solution:

    rho(r) = 3M/(4 pi a^3) [1+r^2/a^2]^(-5/2)
    Phi(r) = -GM/sqrt(r^2+a^2)
    P(r)/rho(r) = GM/[6 sqrt(r^2+a^2)]

The Plummer profile is established stellar-dynamics mathematics, originating with [Plummer (1911)](https://academic.oup.com/mnras/article/71/5/460/972150); polytropic equilibrium is also discussed in [Chavanis (2003)](https://www.aanda.org/articles/aa/full/2003/13/aa2869/node5.html). Applying an effective support pressure to captured companions is hypothetical. Neither the exponent nor the profile is claimed as our discovery. Direct numerical differentiation checks pressure/gravity balance, and integration checks finite unit mass in the G=M=a=1 control.

In plain language, gravity draws deposits inward while their motions or stresses resist unlimited concentration. This offers a possible reason for an extended reservoir. But stationary balance alone does not show that capture populates those motions, that the state is stable, or that one universal K works across systems. A fluid and collisionless-particle interpretation have different stability conditions. No stability claim is made.

With actual cluster gas and stars, the equation must instead use Phi_b+Phi_c. The simple isolated Plummer solution no longer follows automatically. Mass M and scale a were fitted here; their implied K is not a galaxy-calibrated constant. This optional support branch does not replace the one-third retention law, the MOND-guided branch or redshift equations.

## What can be concluded and what cannot

There is no sign or finite-mass obstruction to fitting the saved Coma shear shape. The observational inverse problem is highly nonunique. A known equilibrium family supplies a concrete support candidate with a competitive fitted shape, but not yet a common transport/capture prediction.

The next physical calculation must use matched gas/stars and source geometry, subtract their gravitational contribution consistently, then solve supported deposits with a common K and an explicit captured inventory. It must account for settling/support energy and outgoing radiation rather than manufacture missing mass. Only after those inputs are fixed can another cluster be a meaningful transfer test. The Coma photometry catalog alone does not provide its bolometric external companion bath or exposure history. No galaxy parameter or redshift rate was changed here.

Reproduce with coma-inverse.py. All profiles, extrema, checks and source hashes are saved in coma-inverse-results.json. The manuscript supplement is updated; PDF v1.5 predates this diagnostic.
