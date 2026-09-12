# Adaptive angular refinement: trajectories running

The attribution audit was progress: it supplied concrete angular regions contributing most strongly to the selected force discrepancy. This refinement selects the first 66 ranked first-nested parent regions and splits all edges of their second-mesh children. Adjacent cells are split where necessary to maintain matching edges. Every unselected region remains in the sphere.

The mesh has 1,594 vertices and 3,184 triangles. All edges have two incident faces; the Euler characteristic is two; total solid angle is 4 pi; parent solid-angle discrepancy is at most 4.43e-17. Existing vertices are retained and symmetry reconstruction error is below 2.55e-16. The new mesh requires 698 representative trajectories for R=3: 480 reused and 218 new. The new paths are being integrated over the original full duration, through age 0.25, with unchanged dynamics and tolerances. No lifetime cutoff is introduced.

The geometry checker is prepared to compare with the same reference trajectories and targeted probes. The force implementation is prepared to preserve each second-mesh parent region's mass and divide it among new children by their solid angles. That first comparison isolates positional representation from changing capture quadrature. It is a numerical diagnostic; subsequent capture-weight sensitivity is still required. Explicit fourfold field averaging remains in place.

Mesh invariants and syntax checks pass. Trajectory preparation is incomplete and no adaptive force result is claimed. The earlier polar-domain limitation remains in the reused cache, and selected refinement is not a fresh observational holdout. The other failed force locations remain outstanding.

Next, after the preparation completes, run `check_adaptive_geometry.py`, then `adaptive_volumes.py 256 3 0.1` and `adaptive_volumes.py 512 3 0.1`. Compare the complete force at the declared positions with the archived second mesh, retaining all source ages and mass; do not judge success from the targeted bin alone. Additional epochs, source weighting, independent refinement coverage and all physical/observational closures remain required.

The mesh uses known local triangle subdivision, and the gravitational kernel is unchanged. All nine goals remain open. No claim of a new physical law or a verified stellar/lensing prediction follows from this preparation.
