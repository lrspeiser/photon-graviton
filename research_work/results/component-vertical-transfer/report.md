# Vertical gravity after the Cepheid amplitude calibration

**The same calibrated source predicts a stronger inner vertical pull than the previous empirical field.** Across 90 inner off-plane probes, total vertical gravity is 1.036-1.290 times the old model, with median 1.207. This is an unrefitted model-to-model prediction, not a measurement or observational exclusion.

## Definition and provenance

The amplitude k=1.719397 is taken unchanged from the exposed Cepheid rotation-proxy calibration. Known linear Poisson gravity gives a_new=a_b+k a_component; the previous target gives a_reference=a_b+a_empirical. At z>0 the restoring vertical pull is -a_z. Ratios use the total force including ordinary matter, not only the added component. All source hashes are verified, and no vertical parameter is adjusted.

The stored source model is reflection-symmetric, so the corresponding below-plane restoring force has opposite sign and equal magnitude. That symmetry is imposed by the model, not newly established from observations. The exact plane is excluded from vertical ratios because the symmetric vertical force vanishes there.

| Probe subset | Points | Minimum ratio | Median ratio | Maximum ratio | Maximum vertical refinement / target |
|---|---:|---:|---:|---:|---:|
| inner | 90 | 1.0360 | 1.2068 | 1.2905 | 1.373% |
| all_offplane | 200 | 0.5185 | 1.1992 | 1.3158 | 1.373% |
| solar_R8 | 25 | 1.0360 | 1.1776 | 1.2563 | 1.373% |

The inner subset has R<=8 kpc and 0<z<=1 kpc. The all-offplane grid extends to R=20 kpc and z=4 kpc. This spans very different geometries; the global range is not an uncertainty interval.

## Example at R=8 kpc on the bar axis

| Height (kpc) | Previous total vertical pull | New total vertical pull | Ratio |
|---|---:|---:|---:|
| 0.1 | 754.628 | 787.264 | 1.0432 |
| 0.5 | 2064.036 | 2365.544 | 1.1461 |
| 1 | 2705.805 | 3306.831 | 1.2221 |
| 2 | 3238.567 | 4036.623 | 1.2464 |
| 4 | 3395.138 | 3976.355 | 1.1712 |

Accelerations above are (km/s)^2/kpc. The 1.37% maximum vertical refinement difference is larger than the earlier full-vector normalized difference because near-plane vertical forces have a smaller denominator. It is retained rather than replacing it with the more favorable vector check.

## Consequence for the observation test

Improving radial motion by increasing deposited-source strength also changes the predicted vertical motions. We cannot claim a successful common gravitational explanation by fitting rotation alone. The previous empirical field is not an authoritative measurement of vertical gravity; its disagreement with the new model does not establish which is closer to nature.

The APOGEE/Gaia bulge analysis already showed proper-motion/distance and association sensitivities. A valid comparison must predict velocity distributions with tracer density, cross moments, survey selection and measurement covariance. Raw vertical velocity dispersion cannot simply be equated to vertical force or multiplied by this ratio. These quantities are therefore predictions to carry into the pending three-dimensional stellar analysis, not a completed vertical-data test. Physical photon supply, capture and storage remain unclosed. All six objectives remain open.
