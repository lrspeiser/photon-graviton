# Orbit-library coverage before population fitting

The 72-orbit library is numerically stable but does not yet provide adequate demonstrated coverage for a stellar-population comparison. This is a training-only diagnostic of the available saved paths, not evidence against or for the gravitational formula.

We compare 27,812 selected training stars after excluding all 72 launch stars. Of these, 27,642 pass the existing positional screen. The flagged and unassessed stars remain in a separate all-training result; none is deleted from the underlying catalog. No validation or final-test observations are opened.

For each star, we find saved orbit states within a specified three-dimensional spatial radius and ask whether any also lies within a specified three-dimensional velocity radius. We use the same Cartesian bar frame for both. These radii are diagnostic resolution choices, not measurement error bars. The inputs are conditional posterior means; distance and motion uncertainties have not been integrated into a likelihood.

The search uses 36,000 states: 500 positive-time samples from each of 72 orbits. Initial-time samples are excluded as well as the launch stars themselves. The remaining training stars influenced the original sample scales and seed-selection pool, so this is not an independent holdout.

## Coverage at several diagnostic resolutions

| Spatial radius (kpc) | Velocity radius (km/s) | Position covered | Position and velocity covered |
|---:|---:|---:|---:|
| 0.25 | 25 | 30.23% | 0.98% |
| 0.25 | 50 | 30.23% | 3.71% |
| 0.25 | 100 | 30.23% | 8.78% |
| 0.5 | 25 | 67.65% | 4.72% |
| 0.5 | 50 | 67.65% | 16.87% |
| 0.5 | 100 | 67.65% | 35.12% |
| 1.0 | 25 | 97.47% | 17.57% |
| 1.0 | 50 | 97.47% | 47.76% |
| 1.0 | 100 | 97.47% | 78.20% |

At 0.5 kpc and 50 km/s, 18,699 of 27,642 stars have nearby positions on the saved paths, but only 4,664 have a nearby position and velocity together. A fitted weighting cannot add a missing orbit state. Smooth kernels could spread support beyond these radii, but artificially broad smoothing would blur the physical prediction rather than demonstrate that the library is sufficient.

## Which regions lack coverage?

The following table uses 0.5 kpc and 50 km/s throughout. R is distance from the rotation axis; absolute z is distance above or below the disk plane. Counts are observed sample counts, not selection-corrected stellar densities.

| R (kpc) | Absolute z (kpc) | Stars | Position covered | Position and velocity covered |
|---|---|---:|---:|---:|
| 0.5–3.5 | 0–0.2 | 236 | 100.00% | 28.39% |
| 0.5–3.5 | 0.2–0.5 | 782 | 100.00% | 22.76% |
| 0.5–3.5 | 0.5–1.5 | 528 | 99.81% | 13.07% |
| 3.5–5 | 0–0.2 | 459 | 100.00% | 30.72% |
| 3.5–5 | 0.2–0.5 | 464 | 98.92% | 20.47% |
| 3.5–5 | 0.5–1.5 | 674 | 88.28% | 11.28% |
| 5–9 | 0–0.2 | 4,383 | 76.14% | 23.41% |
| 5–9 | 0.2–0.5 | 8,914 | 76.58% | 22.66% |
| 5–9 | 0.5–1.5 | 11,202 | 48.90% | 8.86% |

The upper bulge and disk regions have particularly sparse velocity coverage. This cannot yet establish that the real bulge needs a different gravity correction: incomplete orbit coverage, uncertain stellar distances, the ordinary-matter model and the hypothesized extra potential can all affect this comparison.

## Saved-time resolution and verification

Keeping every second saved positive-time sample reduces joint coverage at the middle diagnostic resolution from 4,664 to 4,519 stars (16.87% to 16.35%). The modest change suggests that simply saving twice as many nearby points is unlikely to remove the large gap. It is not a continuum bound: longer paths, additional starts and denser time sampling still require separate checks.

Twelve fixed training indices were independently checked against exhaustive distances to all 36,000 states. Search results agree exactly. Larger spatial windows never reduce coverage; thinning the saved states never increases it. Source hashes are checked before and after the calculation. These checks verify the search, not a stellar-population model.

The nearest-neighbor distances use standard Euclidean geometry. The orbit equations are the known rotating-frame Hamilton equations, and the prescribed extra potential uses known QUMOND-style mathematics with empirical coefficients. No new photon-conversion or capture formula is derived by this diagnostic.

## Consequence for the research goal

Do not tune the gravity coefficients to repair this coverage table or present these fractions as predicted stellar agreement. First expand phase-space coverage using a declared training-only rule, test duration/occupation stability, and provide equivalent orbit-population freedom under ordinary matter. Then fit nonnegative population weights with survey selection and correlated distance/motion uncertainties. Library size and regularization choices must be fixed before new validation outcomes are inspected.

For the next library, record separate coverage of signed height, bar angle and velocity components, not just nine broad radius/height bins. Retain the original results and distance flags. The two long outward paths in the current library also require an explicit outer-field/domain and duration treatment before stationary population weights can be interpreted.

The larger scientific requirements remain unchanged: a shared prediction of redshift and event timing, physically specified energy transfer and retention, capture into a gravitational source, and a common response for stellar motions and lensing. The total source-energy budget remains deferred rather than passed. Numerical orbit success alone supplies none of these missing causal links.

## Reproduction

Run `run.py`, then `report.py` from this folder with the existing catalog and trajectory caches. `results.json` records aggregate counts, input hashes and the diagnostic choices. Individual training-only distances are saved in the ignored cache `research_work/data-cache/orbit-population-coverage/training-coverage.parquet`. Source IDs retain integer precision. No source catalog, force coefficient or holdout role is changed.
