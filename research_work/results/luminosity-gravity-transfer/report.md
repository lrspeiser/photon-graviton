# Luminosity term modestly improves galaxy transfer

We fit 89 galaxies and froze the common parameters before evaluating 29 validation and 31 test galaxies. All partitions have been used in earlier research: these results are exploratory transfer, not untouched confirmation. The dataset contains 3,150 measured rotation points. Every observation and prediction is in predictions.json; individual galaxy scores, parameters and input hashes are in results.json.

## Formula and meaning

Candidate empirical formula:

    g_extra = A a_star (g_b/a_star)^p (L_3.6 / 10^10 L_sun)^q
    v_pred = sqrt[R (g_b + g_extra)]

Use consistent SI units for the latter expression. g_b includes the catalog gas and fixed stellar mass/light contributions. a_star=7.24960897e-10 m/s^2 retains the previous redshift normalization. The fitted candidate has A=0.18221887, p=0.38118825, q=0.06007764. At fixed g_b, ten times the luminosity multiplies extra acceleration by approximately 1.148, not ten. The old q=0 model refits to A=0.24229607, p=0.46245874.

**Originality:** circular-speed kinematics and power-law regression are known mathematics. The luminosity factor is an empirical candidate for this project, not a first-principles derivation or a verified new formula in the literature. With A free, the scale a_star can be absorbed into A: this fit does not independently establish a redshift/gravity coupling.

## Actual versus predicted

Each galaxy receives equal weight. Parameters minimize mean squared log10 speed residual on training only; km/s RMS is a complementary descriptive score.

| Group | Galaxies | Old RMS km/s | Luminosity RMS km/s | Old log10 RMS | Luminosity log10 RMS |
|---|---:|---:|---:|---:|---:|
| Training | 89 | 20.0054 | 19.6602 | 0.109349 | 0.107804 |
| Validation | 29 | 27.3769 | 26.3443 | 0.097429 | 0.096080 |
| Test | 31 | 17.1995 | 16.6130 | 0.079141 | 0.075108 |

The frozen luminosity dependence improves both evaluation groups modestly. Both fits converge without hitting parameter boundaries. The q=0 reproduction agrees with the earlier audit within 0.001 km/s in all groups. No galaxy-specific parameters, changes in cuts, or evaluation-driven retuning were used.

## What follows physically

This is a candidate to retain for comparison, not an adopted mechanism. Catalog 3.6-micron luminosity is a stellar population proxy, not total photon power or historical external companion exposure. The same photometry enters ordinary mass, so correlated errors and population mass/light differences can produce an apparent luminosity dependence. An external companion bath need not correlate with a receiver's own luminosity. The result therefore neither identifies conversion nor tests all external-supply models.

Published distances are stipulated facts for this exercise; several were originally Hubble-flow estimates. No expansion law is used in the candidate, but those distances cannot count as independent evidence against expansion. Measurement, inclination, distance and population uncertainties are not jointly marginalized; these scores are not confidence levels or proof of agreement within observational errors.

The next physical step is to replace this luminosity proxy with a specified incident companion energy and capture prediction, and freeze that rule before another galaxy evaluation. It must predict spatial deposition and motion/lensing together, not merely gain another freely fitted normalization. No claim of a complete theory or resolution of the timing/brightness requirements follows. All six research goals remain open.

Reproduce: `python research_work/results/luminosity-gravity-transfer/run.py`.
