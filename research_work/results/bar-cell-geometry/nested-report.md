# Nested trajectory refinement: completed checkpoint

Actual new trajectories improve the cell geometry, and the remaining time-resolution discrepancy passes an additional refinement. The spatial source is still not converged: 14/30 force and 15/30 potential comparisons fail between the parent and child meshes. This is not yet a verified stellar prediction.

## What changed

The parent has 72 full-sphere directions and 140 triangles. The nested mesh has 282 vertices and 560 triangles, obtained by normalized great-circle edge subdivision. Symmetry reduces these to 98 trajectory representatives per launch sphere. For each of the 1 and 3 kpc sources, 18 trajectories were reused and 80 genuinely new trajectories integrated. Interpolating the old trajectories to more vertices would not have provided this test.

All 196 primary cached trajectories pass their Jacobi checks. The cold-particle interpretation, initial 50 km/s radial speed, ordinary-bar field, capture law and source normalization convention are unchanged. The bar normalization remains model-dependent. No physical source amplitude, self-gravity, other Galactic components or observed-star fit is introduced.

## Geometry against actual source 8 trajectories

The same reference trajectories used in the earlier audit are reused here. They are numerical references, not an astronomical holdout. The error diagnostic is not a bound on the resulting gravitational-force error.

| Launch R (kpc) | Age | Parent RMS / R (%) | Nested RMS / R (%) |
|---:|---:|---:|---:|
| 1 | 0 | 3.84 | 1.05 |
| 1 | 0.005 | 1.79 | 0.47 |
| 1 | 0.01 | 2.59 | 0.71 |
| 1 | 0.025 | 4.77 | 1.26 |
| 1 | 0.05 | 14.74 | 4.91 |
| 1 | 0.1 | 35.91 | 14.76 |
| 1 | 0.25 | 66.66 | 34.57 |
| 3 | 0 | 3.81 | 1.06 |
| 3 | 0.005 | 3.37 | 0.93 |
| 3 | 0.01 | 2.82 | 0.97 |
| 3 | 0.025 | 2.82 | 1.57 |
| 3 | 0.05 | 11.96 | 3.80 |
| 3 | 0.1 | 27.90 | 11.29 |
| 3 | 0.25 | 70.20 | 49.73 |

Age units are kpc/(km/s); 0.25 is approximately 244 Myr. Early interpolation improves substantially, while late-time errors remain 34.6 and 49.7 percent of launch radius. This does not establish chaos or reject the companion hypothesis. It shows the remaining inadequacy of the affine source map.

## Gravity and time-resolution checks

Maximum normalized source-mass discrepancy is 2.22e-16. At 256/512 age layers all 30 potential checks pass and 29/30 force checks pass. The retained failure is R 3,T 0.25 at(0.1, 0, 0), changing 8.20 percent. A 512/1024 supplement passes at all five positions; the central change falls to 1.51 percent. The matched 512 parent/child angular comparison is kept separate and still fails 15 potential and 14 force gates. Thresholds remain 2 and 5 percent, with the stated 0.01 floors; these are not observational error bars.

Sixty-digit evaluations of 65 sampled thin-cell/point combinations differ from ordinary arithmetic by at most 1.63e-12 in potential and 2.79e-12 in force per G total source mass. This does not bound every cell's error or its geometric approximation.

## Original polar-domain failure and its sensitivity test

Nine of ten selected tighter-tolerance orbit checks pass in the original cached field. The 1 kpc polar orbit reaches below the original 1e-6 kpc radial boundary during the stricter run. That failure is retained in nested-orbit-checks.json. A primary integrator not sampling that interval does not establish that the original field is defined there.

As a separate numerical test, the polar potential is continued inside r_c with a quadratic expression matching its value and derivative at the boundary:

\[\Phi_c(z)=\Phi(r_c)+\frac{\Phi'(r_c)}{2r_c}(z^2-r_c^2),\qquad |z|<r_c.\]

This is a known local matching construction, not a new gravity law or a microscopic core derivation. The exact symmetry reduction makes the polar trajectory a 1D calculation. Cutoffs 1e-4, 1e-5, 1e-6 kpc are tested with no fit to observations.

All six cutoff/source tests and all four adjacent-cutoff orbit comparisons pass. The largest departure from the cached axis trajectory is 5.97e-07 kpc and 0.00015 km/s. Matched curvature nevertheless ranges from 56289 to 184484 in(km/s)^2/kpc^2 near the cache boundary. Thus the orbit is insensitive at the declared tolerance, but the central curvature is not established as converged. These checks do not erase the original-domain failure.

## What remains

The late source map still needs spatial refinement, particularly where neighbouring trajectories separate and the map folds. Source conservation and exact tetrahedron gravity do not make an inaccurate density reliable. Local or age-dependent angular refinement is the next numerical target; extra interpolated vertices alone are insufficient.

Then the model still needs funded capture, evolving deposited gravity and ordinary-matter response, lensing from a specified metric/stress law, and joint redshift/timing tests. No physical energy-conservation closure or observational success follows from this checkpoint. All nine goals remain active.

Reproduce with prepare_nested.py 1/3, check_nested_geometry.py, check_nested_orbits.py, check_axis_core.py, nested_volumes.py 256/512 for each launch radius and 1024 for R 3,T 0.25, check_roundoff.py and export_nested.py. Numbers are separate command arguments. Runtime completion is recorded in progress.md; formulas and original failures remain in protocol.md and the raw JSON files.
