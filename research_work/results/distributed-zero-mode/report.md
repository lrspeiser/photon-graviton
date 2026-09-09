# Distributed-source first case: conditional rolling from finite radiation

A uniform radiation bath in a finite periodic, fixed three-dimensional volume drives the tested scalar field from rest toward persistent rolling. Its rate follows from the shared Hamiltonian and radiation supply, rather than being prescribed from redshift observations. The initial radiation loses energy and the field receives it; there is no continuing source in this test.

| Initial radiation density | Final n at t=80 | Final n_dot | Derived limiting n_dot |
|---|---:|---:|---:|
| 0 | 1 | 0 | 0 |
| 0.003 | 5.25210 | 0.0696965 | 0.0774597 |
| 0.006 | 7.60803 | 0.102092 | 0.109545 |

Units are dimensionless, K=g=c0=1. The rate tends to |g| sqrt(2 U0/K)/c0 in general units. This is a derived conditional rate for an optional mechanism, not an astronomical calibration or originality claim.

## Verification

The independently integrated background agrees with the exact implicit solution within 1.8e-9 reference-time units. Maximum relative energy error is below 2.8e-11; tightening numerical tolerance changes state values by less than 4.6e-9. The zero-radiation control remains stationary. Across 36 probe checks (three backgrounds, three launch times, two canonical momenta, two launch-offset resolutions), the independently differenced arrival map agrees with the carrier stretch within the frozen tolerances. Full values are in results.json. Reproduce with `python research_work/results/distributed-zero-mode/run.py`.

## Physical interpretation and limitations

Radiation can provide the initial push; field kinetic energy then supports motion without further energy creation. This result differs from the earlier localized source because homogeneous periodic geometry retains a zero spatial mode. It does not establish that ordinary discrete sources produce a uniform cosmic field. Spatial wavelength stays fixed while the reference frequency falls; observable spectral shifts require the unresolved atomic-clock response. n grows without bound and photon speed relative to the reference coordinates falls, which is a physical-completion question, not evidence of observational viability.

This step does not complete the distributed three-dimensional source calculation requested by the plan. Nonuniform transport, open-boundary behavior, potential/restoring terms, source history, gravity, matter standards and astronomical normalization remain untested. The immediate next step is to check whether rolling survives nonuniform finite sources without using a periodic zero mode as an unexplained cosmic assumption. The fresh-catalog audit remains in progress in parallel as a research workstream, not an independent running process.
