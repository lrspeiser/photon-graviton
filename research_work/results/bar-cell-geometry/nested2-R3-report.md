# Second refinement: completed outer-source geometry

The R=3 kpc preparation finished with all 480 representative trajectories present: 98 reused and 382 new. Cache hashes were checked during geometry evaluation. All primary orbit checks pass, with maximum sampled Jacobi error divided by 220 squared of 8.75e-8, below the unchanged 1e-5 gate. These are integration checks, not physical validation or independent strict-tolerance checks of every path.

Against the same exposed source8 reference trajectories, final interpolation RMS falls from 49.73% to 40.26% of launch radius (1.208 kpc). Earlier ages improve more: at age 0.1 it falls from 11.29% to 4.23%, and at 0.05 from 3.80% to 1.45%. Age units are kpc/(km/s); the final age is approximately 244 million years.

All final reference weight still exceeds the 2%-of-launch-radius diagnostic. The four targeted interior probes change from [60.24%, 77.64%, 61.97%, 69.56%] to [77.75%, 26.83%, 46.19%, 71.49%]. Two worsen; they are retained. More vertices do not guarantee monotonic improvement for every curved trajectory.

This establishes partial geometric improvement, not spatial convergence and not a force-error bound. The late map remains poorly represented. The R=1 preparation is still running, and the queued full-map force comparison has not yet run. No coefficients, capture law or source amplitude were changed.

Reproduce this source-only check with `python check_nested2_geometry.py 3`. The no-argument invocation still requires both complete sources and writes the combined result used by the coordinator. Source-only output is kept separately as `nested2-geometry-R3.json`, so it cannot masquerade as the combined result.

All nine goals remain open. Source funding, self-gravity, companion metric and linked redshift/observational requirements remain unresolved.
