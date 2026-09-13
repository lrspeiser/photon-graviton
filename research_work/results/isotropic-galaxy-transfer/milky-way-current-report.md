# Current model predictions for Milky Way circular speeds

No parameters were fitted to Milky Way speeds. All capture, migration and exchange parameters are transferred from previous external-galaxy fits. The exact one-third law remains fixed. Both archived ordinary-matter baselines and every predeclared input sensitivity are retained.

## Data and scope

The 38 cached [Eilers et al.](https://arxiv.org/abs/1810.09466) circular-speed summaries cover about 5-25 kpc and are inferred from stellar kinematics using an axisymmetric Jeans analysis. They are not individual star velocities or direct acceleration readings. The source reports 2-5% systematic uncertainty near the Sun; the tabulated small formal errors do not capture every modeling uncertainty. Its separate dark-halo inference is not used here. The inner barred Galaxy is outside this comparison.

Rd=2.6 kpc follows the earlier [Juric et al. star-count proxy](https://arxiv.org/abs/astro-ph/0510520). It is not a measured total 3.6-micron light scale. Needed luminosity is approximated by archived stellar mass / 0.5 Msun per Lsun. This nominal mass-to-light assumption is a conditional proxy, not a new photometric measurement. Predeclared sensitivity uses Rd=2.08,2.6,3.12 kpc and luminosity factors 0.5,1,2, without selecting a winning case.

## Formula provenance

Known gravity: v_c^2(R)=v_b^2(R)+G M_dep(<R)/R. The spherical deposited density and exact-third retention are the existing proposed capture model. Original ordinary disk/gas forces remain as archived. No dark halo is added.

The 0.37% partial redistribution is trained on SPARC rotation curves; the much larger halo-trained partial fraction comes from normalized halo-profile fits. They are different calibrations. Current binding and released-binding parameters come from the six-lens retained-geometry feedback fits. None was selected using this Milky Way result.

Feedback now uses a spherical angular average of the actual archived stellar Miyamoto-Nagai/Plummer potential plus the spherical monopole of its thin gas disks. It updates the compact and extended companion field. This is an explicitly spherical approximation for migration in a flattened galaxy; it does not infer three-dimensional bulge or vertical orbits. Gas and stellar baseline uncertainty is retained as two alternatives, not marginalized into a confidence interval.

## Fiducial predictions

RMS discrepancy across all 38 circular-speed bins, km/s:

| Model | Ordinary baseline I | Ordinary baseline II |
|---|---:|---:|
| ordinary_matter | 52.57 | 62.34 |
| exact_third | 6.76 | 10.41 |
| rotation_trained_partial | 6.70 | 10.02 |
| halo_trained_partial | 35.57 | 28.56 |
| feedback_binding | 15.08 | 8.32 |
| feedback_released_binding | 12.03 | 5.75 |

These are descriptive speed residuals, not chi-square significance or acceptance thresholds. The data were examined in earlier work, so frozen transfer is not a blind observation test.

## Selected actual-versus-predicted rows

Baseline I, fiducial inputs. Each selected radius is a real catalog row, not an interpolated measurement.

| R kpc | Inferred circular speed | Ordinary matter | Exact third | Rotation-trained partial | Feedback binding |
|---|---:|---:|---:|---:|---:|
| 5.27 | 226.83 | 197.32 | 218.97 | 220.01 | 256.92 |
| 8.19 | 228.86 | 188.43 | 226.60 | 227.24 | 248.24 |
| 12.25 | 222.23 | 170.34 | 223.67 | 224.05 | 234.64 |
| 15.22 | 217.07 | 158.70 | 217.79 | 218.06 | 224.96 |
| 20.27 | 199.84 | 142.65 | 205.76 | 205.94 | 209.66 |
| 24.82 | 198.42 | 131.15 | 194.87 | 195.00 | 197.35 |

## Outcome

The fixed exact-third model substantially improves the Milky Way circular-speed curve relative to both ordinary baselines: RMS 6.76/10.41 km/s versus 52.57/62.34. The SPARC rotation-trained small redistribution changes this to 6.70/10.02. These are conditional predictions using the declared disk-scale and luminosity proxies, not fits to these speeds.

The recent halo-trained large partial fraction overconcentrates the Milky Way and gives 35.57/28.56 km/s. The newer local feedback adjusts its mobile fraction to this galaxy instead of imposing the same approximately 36% split: binding gives about 14.2%/13.5% mobile, and released binding 9.6%/8.1%, for the full modeled inventory. Released-binding feedback improves baseline II to 5.75 km/s but worsens baseline I to 12.03. Therefore no variant dominates both ordinary-matter choices.

Near the Sun at the actual 8.19 kpc row, inferred speed is 228.86 km/s. Under baseline I, ordinary matter predicts 188.43, exact third 226.60, and the rotation-trained partial variant 227.24. This is a useful radial prediction; it neither establishes individual stellar orbits nor proves the photon-to-gravity mechanism.

The scale/luminosity sensitivity is material. For baseline I, exact-third RMS ranges roughly 6.6-18.5 km/s across the declared grid, while some feedback variants worsen substantially at the larger disk scale. These ranges are exploratory input sensitivity, not confidence intervals; the fiducial choice is retained, and no best-case proxy is selected after looking at the data. Further progress should independently constrain those light/matter inputs and test another kinematic sample, rather than adjusting them to force the curve to match.

## Sensitivity and checks

Every 18 baseline/scale/luminosity scenario is in the JSON with all 38 predictions, residuals and inner/outer scores. No parameter was adjusted to minimize these errors. Predicting individual speeds additionally requires orbital phases and velocity distributions; the radial curve cannot establish that.

Doubling radial and angular resolution changes a fiducial predicted speed by at most 0.0003286 km/s. Analytic capture inventory checks pass. Feedback conserves deposited inventory and converges from opposite initial fractions. Attenuation paths extend to infinity; the outer radial integration is a numerical approximation, not a universe-size limit.

Absolute photon supply, formation work, support, source histories and traveling radiation stresses remain unclosed. The previous provisional vertical summaries are not promoted to independent constraints here. This run tests actual circular-speed predictions, without claiming a complete Milky Way model or completion of the six broader goals.
