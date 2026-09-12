# Angular trajectory geometry audit

The current coarse volume cells are not a sufficiently accurate representation of the moving source. An independent comparison with actual source8 trajectories shows that source6 affine interpolation develops large position errors. This is a numerical representation problem; it does not reject the companion hypothesis or prove that the underlying trajectories are chaotic.

## What was compared

Each source8 launch ray is located in a source6 spherical-hull triangle. Positive affine weights reproduce the ray direction, and the same weights interpolate the three parent trajectories. The prediction is compared with the actual cached source8 trajectory. Errors use normalized source8 capture weights and include all four symmetry partners. Both sets of trajectories use the same ordinary-bar field and physical launch rule. This is a numerical reference comparison, not a fresh astronomical holdout.

All ray reconstructions agree within1e-14 and weights are nonnegative to the declared numerical tolerance and sum to one. At age zero the error is the initial sphere-versus-chord approximation. Later errors include the curved direction-to-position map; simply subdividing the time steps cannot fix this angular approximation.

| Launch radius (kpc) | Age (kpc/(km/s)) | RMS position error (kpc) | RMS / launch radius (%) | Source weight beyond 2% radius (%) |
|---:|---:|---:|---:|---:|
| 1 | 0 | 0.03841 | 3.84 | 100.00 |
| 1 | 0.005 | 0.01788 | 1.79 | 27.22 |
| 1 | 0.01 | 0.02594 | 2.59 | 71.45 |
| 1 | 0.025 | 0.04773 | 4.77 | 84.03 |
| 1 | 0.05 | 0.14744 | 14.74 | 100.00 |
| 1 | 0.1 | 0.35914 | 35.91 | 95.55 |
| 1 | 0.25 | 0.66657 | 66.66 | 100.00 |
| 3 | 0 | 0.11416 | 3.81 | 100.00 |
| 3 | 0.005 | 0.10096 | 3.37 | 100.00 |
| 3 | 0.01 | 0.08455 | 2.82 | 88.08 |
| 3 | 0.025 | 0.08461 | 2.82 | 69.33 |
| 3 | 0.05 | 0.35866 | 11.96 | 100.00 |
| 3 | 0.1 | 0.83708 | 27.90 | 100.00 |
| 3 | 0.25 | 2.10602 | 70.20 | 100.00 |

The final age is approximately244Myr. Large particle-position error does not directly specify a force-error bound, since contributions can average and cancel. Nevertheless it prevents treating the coarse cell density as a resolved formation result. The per-cell errors and every sampled direction are retained in results.json. The2 percent diagnostic is not an observational tolerance.

## Refinement underway

Every parent triangle is split at normalized great-circle edge midpoints, giving282 vertices and560 children. The summed solid angle remains4pi. Symmetry matching leaves98 distinct trajectory representatives per launch sphere:18 existing directions can be reused and80 require actual integration. The full mesh reconstruction error is1.36e-16. Adding interpolated trajectories would not constitute this refinement.

The new directions include a polar-axis case. The bar cache has even angular harmonics; its exact axis potential and axial acceleration follow from its m=0 terms, and the transverse force vanishes by symmetry. Nearby off-axis evaluations verify that limit. This is known harmonic mathematics, not a modified force law. The radial cache domain is unchanged; a trajectory outside that domain must be reported rather than hidden.

Both preparation processes remain active at this checkpoint. No improved interpolation or force result is claimed yet. After completion, check new trajectory accuracy, test interpolation against the existing finer-grid trajectories, and recompute mass-preserving volume fields with the same convergence gates.

## Scope

The source, cold-particle capture interpretation and prescribed ordinary bar remain conditional. The ordinary bar normalization is model-dependent. This audit does not fund capture, include deposited self-gravity, predict lensing or resolve the photon/time redshift mechanism. All nine goals remain active.

Reproduce the completed audit with run.py, mesh.py and verify_axis.py. prepare_nested.py takes launch radius1 or3 as a separate argument. Runtime handles and completion criteria are in progress.md.
