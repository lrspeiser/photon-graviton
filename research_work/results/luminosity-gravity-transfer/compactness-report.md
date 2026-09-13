# Radiation-intensity proxy: mixed transfer results

The previous turn made concrete progress by fitting and publishing luminosity-dependent galaxy predictions. This follow-up tests a fixed size dependence against the same actual rotation data, rather than further data-resolution work.

## Formula

    X = (L_3.6 / 10^10 L_sun) / (R_disk / 3 kpc)^2
    g_extra = A a_star (g_b/a_star)^p X^q
    v_pred^2 = R (g_b + g_extra)

The known inverse-square dilution of point-source radiation motivates X as a characteristic intensity proxy. An extended disk is not a point source: this is not the local incident radiation integral, and it does not calculate companion production or deposition. The power-law modification is an empirical postulate for this project using known mathematics, not a first-principles derivation or an established originality claim.

Fitting only the 89 historical training galaxies gives A=0.16829707, p=0.37372303, q=0.10272639, retaining a_star=7.24960897e-10 m/s^2. The geometric inverse-square exponent is fixed, not independently fitted. At fixed ordinary acceleration and luminosity, doubling disk scale length reduces the predicted extra acceleration by 2^(-2q), approximately 0.867. This is a conditional comparison, not an operation that physically changes a galaxy while leaving all other properties fixed.

## Frozen predictions on actual rotation curves

| Group | Baseline RMS km/s | Luminosity only RMS km/s | Luminosity/area RMS km/s |
|---|---:|---:|---:|
| 89 training galaxies | 20.005 | 19.660 | 19.213 |
| 29 validation galaxies | 27.377 | 26.344 | 27.358 |
| 31 test galaxies | 17.200 | 16.613 | 16.109 |

| Group | Baseline log10 RMS | Luminosity only log10 RMS | Luminosity/area log10 RMS |
|---|---:|---:|---:|
| Training | 0.109349 | 0.107804 | 0.107825 |
| Validation | 0.097429 | 0.096080 | 0.096990 |
| Test | 0.079141 | 0.075108 | 0.076562 |

Every galaxy receives equal weight. The declared fitting objective is squared log10 velocity residual; absolute km/s error is complementary. The intensity proxy modestly improves both groups relative to the original baseline. Against luminosity alone it is worse in both groups on the declared fractional-error metric, while absolute error improves only in the test group. Consequently it does not justify replacing luminosity-only as the empirical candidate on these results. Both remain exploratory: repeated use of these partitions is not blind validation and model selection against them would further consume their independence.

The optimizer converged away from bounds. Identical-input q=0 predictions reproduce the previous baseline within 1e-8 km/s; input hashes match. compactness-predictions.json records all 3,150 measured and predicted speeds, while compactness-results.json records each galaxy's scores. There is no uncertainty-marginalized likelihood or claim that remaining residuals match measurement errors.

## Consequence for the mechanism

A positive fitted intensity dependence is compatible with more radiation contributing more extra gravity, but compatibility is not identification. Stellar mass and this intensity proxy share photometry, and fixed mass/light assumptions may generate the dependence. The rule lacks external illumination, spectral energy, age-independent radiation histories, capture probabilities and a spatial deposit calculation. The redshift scale remains absorbable into free A. Neither this result nor the preceding luminosity test derives photon conversion.

The useful constraint is that a more physical source/capture law must outperform these simple proxies with frozen parameters and predict the radial profile, rather than treating the small improvement as permission to assign arbitrary deposits. The motion/lensing, time/brightness and cluster requirements remain unresolved. All six goals stay open.

Reproduce: `python research_work/results/luminosity-gravity-transfer/compactness.py`.
