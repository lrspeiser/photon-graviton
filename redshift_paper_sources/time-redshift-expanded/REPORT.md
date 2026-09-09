# Time replacing expansion: expanded holdout test

**Result: on 25 withheld galaxy-group representatives, the hypothetical time model and a specified expanding-universe reference make essentially indistinguishable redshift predictions. Neither is selected by this test.** This is a real-data calculation using galaxies with distances measured from stellar populations, not a sample of individually resolved stars. It does not test replacing dark matter.

## The time number

The fitted static propagation model has a logarithmic time-stretching rate of **7.7314966e-05 per million light-years of geometric path**, or **0.007732% added stretching for each million light-years**. A one-second gap would grow by about 77.32 microseconds over that path. Applying the wave-spacing rule to a whole signal assumes the same stretching affects its envelope.

Let R be geometric path length in millions of light-years. Then:

    arrival_interval / emitted_interval = 1 + z = exp(7.7314966e-05 R)

The coefficient alpha is 0.000252167692 per Mpc, corresponding to c alpha = 75.5980 km/s/Mpc. A training sky-block bootstrap gives a conditional 95% range of 73.053 to 78.656 km/s/Mpc for that slope. This is not a complete systematic uncertainty; spatial correlations and common calibration remain.

| Geometric path | A one-second emitted interval becomes |
|---|---:|
| 1 million light-years | 1.000077318 seconds |
| 10 million light-years | 1.000773449 seconds |
| 100 million light-years | 1.007761462 seconds |

These are cumulative arrival-time factors, not measurements of clocks sitting in empty space. A fixed spatial clock-rate gradient by itself does not establish a cumulative propagation law. The hypothetical rule is stipulated and still needs a physical mechanism.

## Sample and separation

Cosmicflows-4 table 2 contains a method-specific surface-brightness-fluctuation (SBF) distance modulus. SBF uses the statistics of unresolved stellar light to estimate distance. It supplies distances that are not simply calculated by converting the target redshift through a Hubble law. However, published SBF calibration and photometric corrections still carry assumptions and were not reconstructed under the hypothetical theory.

We selected SBF distance estimates from 10 to 150 Mpc with distance-modulus uncertainty at most 0.30 mag. All 68 galaxy groups from the earlier pilot were excluded, including its training and held-out groups. There were 213 qualifying galaxy measurements, reduced to 164 group representatives by retaining the smallest SBF-error measurement in each dominant-PGC group (ties broken by PGC). No cut used a target redshift.

The sky was divided into 24 right-ascension bins and four equal sin(declination) bins. A fixed hash of the tile assigned training, validation, or test status. There were no split retries. The result was:

| Subset | Group representatives | Occupied sky tiles |
|---|---:|---:|
| Training | 104 | 36 |
| Validation | 35 | 16 |
| Final test | 25 | 13 |

All fits used only training labels. Validation was diagnostic, not a model-selection or retuning step. Parameters and predictions were saved and hashed before scoring either held-out subset. Catalog preparation necessarily extracts labels into separate files; this is a computational holdout, not independently blinded observations. Shared catalog calibration and neighboring structures can cross tile boundaries.

## Prediction comparison

All errors below compare predicted redshift to the catalog CMB-frame systemic redshift, using cz as the velocity-unit representation. The number in km/s is a redshift error multiplied by c, not a claim that all shifts are physical recession velocities.

| Model | Validation RMSE (km/s) | Final test RMSE (km/s) | Final test MAE (km/s) |
|---|---:|---:|---:|
| Static time model with photon-conserving brightness relation | 435.927 | 412.441 | 336.927 |
| Flat expanding reference, fixed matter fraction 0.3 | 435.900 | 411.934 | 336.252 |
| Earlier exponential distance-as-path approximation | 436.301 | 414.627 | 339.437 |
| Static time model with earlier pilot rate unchanged | 444.872 | 427.645 | 355.102 |

The expanding reference wins numerically by 0.507 km/s in final-test RMSE, about one tenth of one percent of the overall error. A paired 5,000-draw sky-tile bootstrap for time-minus-expansion RMSE spans -0.814 to 1.948 km/s. It includes zero. Mean negative-log-likelihood differences likewise do not distinguish the models. The expanding reference has slightly lower mean NLL on the final test; the time model has slightly lower mean NLL on validation. There is no convincing winner.

The nominal conditional 95% prediction bands cover 22 of 25 final-test objects (88%). These bands omit parameter uncertainty and common systematics and use an approximate motion model; they should not be treated as validated uncertainty intervals.

The earlier pilot static coefficient was c alpha = 69.5192 km/s/Mpc. The new fit is 8.74% higher. Applying that earlier rate unchanged gives a final-test RMSE of 427.64 km/s. Therefore a universal time coefficient has not been established. Sample differences, calibration, motion and model misspecification remain possible explanations for the shift.

## How brightness is handled without feeding in the answer

Our main time model assumes static Euclidean geometry, no opacity or photon loss, and equal stretching of wave spacing and signal duration. Photon energies and arrival rates each decrease by 1/(1+z), so the luminosity distance satisfies D_L = r(1+z). Meanwhile 1+z = exp(alpha r). This predicts redshift from the brightness-based distance alone:

    z_predicted = exp(W(alpha D_L)) - 1

W is the Lambert W function, the inverse of w exp(w). No target redshift is used on the right-hand side. The code uses distance in Mpc and alpha in inverse Mpc. For geometric path predictions in million light-years, use the earlier coefficient.

This brightness mapping is a stated model assumption. Published SBF distances and corrections were not reprocessed from raw imaging with an alternative transport model, and correlated calibration is retained. This prevents interpreting the exercise as a clean empirical proof of static geometry.

The earlier approximation uses z = exp(alpha D_L)-1 and treats luminosity distance as path length. It is included only as a diagnostic to connect to the first pilot; it is not the primary physical comparison.

The expansion reference is flat FLRW with Omega_m=0.3 and Omega_Lambda=0.7 fixed, negligible radiation, and a fitted H0:

    D_L = (c/H0)(1+z) integral_0^z [0.3(1+u)^3+0.7]^(-1/2) du

The code numerically inverts this relation to predict redshift from the same catalog distance. It fits H0 = 75.8848 km/s/Mpc. Both primary models have one fitted scale parameter. This is a specified expansion reference, not a global fit of every cosmological parameter or all cosmological observations.

## Statistical assumptions and checks

Each model minimizes a Gaussian negative log likelihood in catalog CMB-frame cz. Variance is a fixed 300 km/s motion floor squared plus the first-order propagated SBF distance error squared. No catalog-derived peculiar velocity or cosmological curvature-adjusted velocity is used. The representative galaxy's individual systemic velocity is used; group internal motion is not removed. The model does not rebuild correlated flow covariance, source gravitational shifts, survey selection, or common distance-scale uncertainty. A first-order propagated error model is an approximation rather than a full latent-distance analysis.

Training-only profile intervals and 1,000 sky-tile bootstrap refits are supplied. Paired final-test bootstrap draws keep predictions and measurements together within sky tiles. Tiles are an approximate blocking strategy; independence between all tiles is not guaranteed.

Verification checked unique groups, complete exclusion of earlier groups, no tile split across subsets, preservation of frozen prediction and parameter hashes, the static brightness identity, and agreement between FLRW interpolation and independent numerical quadrature. No held-out result was used to change the model, filter, or split.

## What has and has not been tested

This test asks whether a single fitted static time-propagation rate predicts additional nearby redshifts about as well as a specified expansion distance law, conditional on the published distance measurements. The answer is yes within this limited comparison, but it does not distinguish the causes. At low redshift, both relations approach the same linear distance-redshift law.

No density, void-path, gravitational potential or acceleration data were supplied. Therefore this run cannot test whether stretching is suppressed inside galaxies or enhanced near galaxy outskirts. It cannot show that actual stars orbit differently. There are no rotation curves, lensing measurements or dynamical equations in this analysis, so it does not test replacing dark matter. A universe with the same full observational behavior would need to reproduce those additional data, not merely this redshift relation.

The next decisive redshift test needs environments measured independently enough to predict different shifts at similar distances, or distance data with sufficient range and a self-consistent reduction to distinguish the two distance laws. The current final test is now exposed; it must not be reused as an untouched holdout after developing new formulas.

## Sources

- Tully et al., Cosmicflows-4: https://arxiv.org/abs/2209.11238
- Catalog and field definitions: https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/ApJ/944/94?format=html&tex=true
- Downloaded data: https://cdsarc.cds.unistra.fr/ftp/J/ApJ/944/94/table2.dat.gz
- SBF method/application: https://arxiv.org/abs/2101.02221
- Cosmological distance definitions: https://arxiv.org/abs/astro-ph/9905116

## Reproduce

Requires Python, NumPy and SciPy. In the expanded-test directory, run these stages in order:

    python run_test.py prepare
    python run_test.py train
    python run_test.py predictions
    python run_test.py verify
    python run_test.py score
    python make_report.py

The archive includes the full source catalog, its documentation, prior-pilot exclusions and parameters, protocol, fitting code, all predictor inputs, separately stored labels, frozen predictions, scores, and verification output. No web fetch is needed to rerun it. A fresh rerun reproduces the computational split, not a newly blinded scientific test.

## Every final-test prediction

The distances below are the published luminosity-distance estimates, not geometric path lengths. No final-test objects have been omitted.

| PGC | Catalog distance (Mpc) | Time-model predicted z | Expansion predicted z | Observed z, CMB frame |
|---|---:|---:|---:|---:|
| 27765 | 15.517 | 0.0039052 | 0.0039158 | 0.0049167 |
| 53201 | 19.534 | 0.0049139 | 0.0049259 | 0.0056606 |
| 69342 | 21.817 | 0.0054866 | 0.0054991 | 0.0027352 |
| 5766 | 22.019 | 0.0055372 | 0.0055498 | 0.0052636 |
| 26939 | 22.121 | 0.0055627 | 0.0055753 | 0.0057406 |
| 25915 | 22.532 | 0.0056658 | 0.0056785 | 0.0049768 |
| 25955 | 23.594 | 0.0059321 | 0.0059449 | 0.0076186 |
| 68870 | 24.706 | 0.0062108 | 0.0062237 | 0.0049501 |
| 51275 | 25.870 | 0.0065025 | 0.0065156 | 0.0057807 |
| 20047 | 26.230 | 0.0065927 | 0.0066057 | 0.0078588 |
| 51787 | 26.965 | 0.0067768 | 0.0067899 | 0.0063777 |
| 6983 | 27.593 | 0.0069341 | 0.0069472 | 0.0046265 |
| 19476 | 32.569 | 0.0081794 | 0.0081921 | 0.0091530 |
| 46115 | 32.719 | 0.0082170 | 0.0082296 | 0.0103338 |
| 71730 | 38.976 | 0.0097809 | 0.0097918 | 0.0095433 |
| 45657 | 39.337 | 0.0098709 | 0.0098817 | 0.0108608 |
| 33871 | 47.076 | 0.0118017 | 0.0118085 | 0.0111644 |
| 31586 | 53.064 | 0.0132931 | 0.0132955 | 0.0111310 |
| 4224 | 59.924 | 0.0149989 | 0.0149949 | 0.0166248 |
| 6962 | 60.145 | 0.0150539 | 0.0150496 | 0.0153006 |
| 66537 | 60.478 | 0.0151367 | 0.0151320 | 0.0146968 |
| 5037 | 63.096 | 0.0157867 | 0.0157792 | 0.0127688 |
| 50558 | 69.759 | 0.0174398 | 0.0174239 | 0.0169284 |
| 5324 | 73.858 | 0.0184554 | 0.0184337 | 0.0174154 |
| 68434 | 77.804 | 0.0194320 | 0.0194042 | 0.0186262 |
