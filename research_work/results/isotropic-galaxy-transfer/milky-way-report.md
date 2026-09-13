# Milky Way transfer: rotation improves, joint agreement remains incomplete

The previous turn established improved SPARC transfer for an attenuated companion bath. This run keeps those fitted capture parameters unchanged and predicts Milky Way rotation and vertical pull from the same spherical deposits. Nothing is fitted to the Milky Way outcomes.

## Inputs and formula status

Use R_disk=2.6 kpc, the thin-disk star-count scale reported by [Juric et al.](https://arxiv.org/abs/astro-ph/0510520). Identifying this tracer scale with the 3.6-micron light scale used in SPARC is a conditional assumption. The star-count result includes population/binary corrections; this is not a precision determination of the Milky Way's total infrared scale.

The capture coefficient and retention assumptions are postulates documented in derivation.md. The additional radial and vertical forces use known spherical gravity:

    v_pred^2(R)=v_b^2(R)+G M_d(<R)/R
    |Kz_pred(R,z)|=|Kz_b(R,z)|+G M_d(<r)|z|/r^3,
    r=sqrt(R^2+z^2).

No separate vertical normalization is allowed. The direction of incoming companions does not determine the direction of the force: the deposited mass distribution does.

Both previously specified ordinary-matter baselines I and II are preserved, with their disk/gas/central-structure uncertainty. The comparison includes 38 [Eilers et al. circular-speed summaries](https://arxiv.org/abs/1810.09466). These are kinematically inferred speeds with equilibrium/tracer assumptions, not direct accelerometer readings. The 43 [Bovy-Rix vertical summaries](https://arxiv.org/abs/1309.0809) are provisional: their inference used potential/distribution-function families including halo assumptions. They are not adopted as independent dark-matter-free facts or counted as goal-4 completion. No halo source is added in our prediction.

## Frozen rotation predictions

RMS speed error in km/s across all 38 bins:

| Ordinary baseline | Ordinary matter only | Empirical transfer | Transparent capture | Intercepted capture |
|---|---:|---:|---:|---:|
| I | 52.575 | 12.488 | 18.107 | 14.035 |
| II | 62.340 | 7.544 | 26.173 | 21.518 |

Interception improves on transparent capture for both baselines, but the empirical relation still has lower overall error. Intercepted predictions have mean speed bias -11.02 km/s (I) and -19.06 km/s (II), so residual underprediction remains.

For the 18 outer bins, interception RMS is 10.813 km/s (I) and 15.381 km/s (II), compared with transparent 18.419 and 24.894. The 20 inner bins have intercepted RMS 16.403 and 25.823. These are descriptive residuals, not uncertainty-calibrated acceptance scores. All observations have been examined previously; transfer here means frozen parameters, not blind validation.

## Same deposits, vertical comparison

RMS discrepancy in |Kz|/(2piG), expressed in Msun/pc^2, across the 43 provisional summaries:

| Ordinary baseline | Ordinary matter only | Empirical transfer | Transparent capture | Intercepted capture |
|---|---:|---:|---:|---:|
| I | 15.579 | 23.813 | 17.006 | 16.917 |
| II | 19.036 | 16.560 | 16.507 | 16.849 |

The vertical result is mixed: interception slightly improves on transparent capture under I and worsens it under II. Do not pool these numbers with speed errors or treat the rows as independent likelihood measurements. Their inherited modeling assumptions prevent interpreting the table as independent confirmation or rejection of companion gravity.

## Assessment

The capture/interception rule produces useful rotation improvements with no Milky Way refit, but it does not yet reproduce the observed motions adequately or demonstrate a common three-dimensional gravitational solution. It cannot describe differences around a barred bulge with a spherical deposit profile and these outer-disk summaries alone.

The next spatial model must retain this rotation constraint while predicting an independently inferred vertical/stellar distribution or lensing observable. Changing ordinary matter, disk scale or exposure to improve a particular outcome would be a new declared fit, not continued frozen transfer. Actual photon supply, capture microphysics, stable retention, collisions and redshift/timing/brightness remain unresolved. All six goals remain open.

Reproduce: `python research_work/results/isotropic-galaxy-transfer/milky-way.py`. milky-way-results.json records input hashes and inner/outer/vertical scores; milky-way-predictions.json includes every observation, prediction and residual under both baselines and all four models. The integration radius 30 kpc only covers requested force evaluations; opacity rays still extend to infinity and it is not a physical deposit cutoff.
