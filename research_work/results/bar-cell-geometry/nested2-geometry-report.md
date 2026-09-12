# Second full-sphere refinement: geometry complete, gravity running

Both preparation runs completed successfully: 480 representative paths per launch sphere, including 98 reused and 382 new paths per source. All 960 primary orbit checks pass; maximum sampled Jacobi error divided by 220 squared is 4.07e-8 for R=1 and 8.75e-8 for R=3. The geometry calculation hash-checks every trajectory cache. These checks do not establish the physical validity of the ordinary-bar or capture assumptions, nor erase the retained polar-domain limitation of reused paths.

## Measured geometry improvement

The same exposed source8 reference trajectories and unchanged normalized reference weights are used at every refinement. Numbers below are weighted RMS position error divided by launch radius, in percent. They are numerical diagnostics, not measurement uncertainties or bounds on gravitational-force error.

| Source radius | Age | First nested mesh | Second nested mesh |
|---|---:|---:|---:|
| 1 kpc | 0.05 | 4.91 | 1.45 |
| 1 kpc | 0.1 | 14.76 | 6.13 |
| 1 kpc | 0.25 | 34.57 | 23.35 |
| 3 kpc | 0.05 | 3.80 | 1.45 |
| 3 kpc | 0.1 | 11.29 | 4.23 |
| 3 kpc | 0.25 | 49.73 | 40.26 |

Age units are kpc/(km/s); the final age is approximately 244 million years. At that age, 97.52% of inner reference weight and essentially all outer reference weight still exceed the earlier 2%-of-radius diagnostic. The map has improved but is not converged.

The four targeted inner probes change from [23.47%, 25.18%, 29.32%, 27.95%] to [5.11%, 11.74%, 7.78%, 11.80%]. The outer probes change from [60.24%, 77.64%, 61.97%, 69.56%] to [77.75%, 26.83%, 46.19%, 71.49%]; two worsen and remain in the result. These targeted probes are development samples, not unseen astronomical observations.

## Current execution and next decision

The coordinator has passed preparation checks and completed the combined geometry evaluation. Both 256-layer volume-force jobs have started; the 512-layer jobs follow. This uses the same exact tetrahedral gravity, source mass convention and force/potential gates. Completion of paths alone is not a force result. No new comparison is claimed until the jobs finish and the output is inspected.

Actual gravitational stability must be checked separately from interpolation improvement. Afterward, physical source funding, self-gravity, stellar response, companion metric, redshift and linked observations still require work. No new physical law or parameter fit is introduced. All nine goals remain open.

Reproduce geometry with `python check_nested2_geometry.py` once both preparations are complete. `nested2-geometry.json` contains all seven ages and all targeted-probe errors. The source-only outer result and earlier failures remain archived.
