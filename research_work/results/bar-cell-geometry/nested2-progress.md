# Second full-sphere refinement: preparation in progress

The previous turn was progress: new interior trajectories confirmed poor interpolation outside the existing reference coverage. The next calculation refines every triangle rather than assuming unchecked regions pass.

The second mesh has 1,122 vertices, 3,360 edges and 2,240 triangles. Each edge has two adjacent faces; the Euler characteristic is two; all parent vertices are preserved; solid angles sum to 4 pi within roundoff. Symmetry matching requires 480 representative trajectories per launch sphere: 98 reused paths and 382 new ones, or 764 new integrations across both sources. Representation checks do not imply gravitational convergence.

`prepare_nested2.py 1` and `prepare_nested2.py 3` are running with the unchanged source, capture and ordinary-bar dynamics. New caches use distinct nested2 filenames. Preparation JSON files are live partial outputs and are not part of this commit. No completed trajectory set or new force prediction is claimed. The earlier polar-domain limitation remains in the reused paths; this mesh does not resolve or erase it.

Once both preparations finish, `check_nested2_geometry.py` compares the actual new paths against the previously exposed source8 references and the eight targeted interior probes. Those are numerical development checks, not astronomical holdouts. Then compute the same volume gravitational fields, checking time and direction resolution separately. Preserve mass and the existing force/potential tolerance gates.

All nine goals remain open. This is progress on numerical deposition dynamics, not closure of the photon-to-companion mechanism, energy supply, self-gravity, clock law, lensing or observational fit.

Files: `mesh2.py`, `mesh2.json`, `mesh2-checks.json`, `prepare_nested2.py`, `check_nested2_geometry.py`. The mesh uses known geodesic midpoint subdivision and affine interpolation; no new physical law is introduced.

## Dependent calculations queued

`nested2_volumes.py` retains the exact tetrahedral gravity solver and unchanged target positions, source weighting and 256/512 age-layer comparison. It records the field and mesh hashes. `export_nested2.py` compares the second mesh with the first at matched 512 age layers, keeping the 2% potential and 5% vector-force thresholds.

A deliberate preflight invocation on an incomplete preparation was rejected before field calculation; this verifies the completeness guard, not numerical accuracy. Syntax checks pass.

`finish_nested2.py` waits on the two existing Windows preparation process IDs supplied at launch. It never restarts preparations. After process termination it checks that every expected record exists and passes its primary orbit check, then runs geometry, four volume calculations (at most two simultaneously), and the comparison exporter. Failed or incomplete preparations stop dependent work. Force-job logs are retained. The coordinator is currently waiting on the live preparation processes; no resulting force comparison is claimed yet.
