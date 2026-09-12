# Where the angular representation fails

The earlier checkpoint was progress: it isolated unresolved spatial errors after the extra time-layer check passed. This follow-up measures how broadly those errors occur and adds eight actual trajectories in previously unchecked cells. It changes no physical coefficients or source mass.

## Coverage and concentration

The existing 128 reference directions occupy only 128 of 560 triangles. Those triangles contain 23.55% and 23.56% of the normalized source mass for the 1 and 3 kpc launch spheres. These are the mass fractions of cells containing probes, not fractions of material validated: one interior ray cannot validate a whole cell. The other cells must not be treated as passing.

At the final age, 73 and 59 of the 128 sampled cells account for 90% of the measured squared interpolation error. Almost all reference weight exceeds the earlier 2%-of-launch-radius diagnostic. Thus a handful of sampled outliers cannot explain the entire failure. The ranked-cell JSON retains every nonzero contribution.

Across the full mesh, cells containing 45.0% and 68.5% of source mass have vertex-to-vertex diameters exceeding their launch radius at the final age. Large diameter is a refinement indicator, not proof of interpolation error: an affine map could be accurate even for a large cell.

## Actual new trajectories

For each sphere, select four previously unchecked symmetry classes by source mass times squared final cell diameter. Integrate the normalized angular centroid with the unchanged ordinary-bar field, launch velocity and DOP853 tolerances. Compare it with the three vertex trajectories' affine prediction. These are targeted numerical development samples, not random samples or astronomical holdouts.

| Launch radius | New trajectories | Final position error / launch radius |
|---|---:|---:|
| 1 kpc | 4 | 23.47%, 25.18%, 29.32%, 27.95% |
| 3 kpc | 4 | 60.24%, 77.64%, 61.97%, 69.56% |

All eight selected orbit invariant checks pass; maximum sampled Jacobi drift divided by 220 squared is 1.86e-8, below the unchanged 1e-5 gate. This checks numerical orbit integration, not the physical validity of the bar or capture model. The original polar-domain failure remains archived and is not erased by these off-axis checks.

## Consequence for the next calculation

Refining only the largest errors among existing reference directions would leave most cells unchecked. The next spatial refinement needs coverage across the entire sphere, with actual new trajectories and preserved source mass. Targeted centroid probes can supplement that refinement, but cannot certify all unprobed cells. Recompute the gravitational field and retain the existing 2% potential and 5% force convergence gates; improved interpolation alone is insufficient.

No trustworthy stellar-motion or lensing prediction follows yet. Energy funding, collective gravity, clock/redshift coupling and unseen-data validation remain open, as do all nine research goals.

## Reproduction and provenance

Run `python diagnose_local_refinement.py`, then `python probe_unsampled.py 1` and `python probe_unsampled.py 3` from this directory. The source caches are hash-checked. Output files retain field/mesh hashes for the new probes, all seven ages, actual and interpolated positions, and orbit-check results. Reconstructed RMS values agree with `nested-geometry.json` to 1e-12.

Affine interpolation, weighted squared-error ranking and numerical orbit integration are known mathematics. Their use here diagnoses this proposed companion-deposit model; no new physical law or literature novelty is claimed.
