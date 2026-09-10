# Field representation matters for long-term population predictions

Changing only the numerical resolution of the candidate extra-potential correction changes the late-time orbit population more than changing the orbit solver tolerance. This rules out treating the previous small instantaneous force-refinement percentage as sufficient evidence of long-term population accuracy. It does not reject the gravity law or establish continuum error bounds.

The same targeted training start, index 12, was integrated in the existing lower-resolution correction at relative tolerances 2e-11 and 2e-13. Both original higher-resolution runs were reused. Physical parameters and ordinary-matter components were unchanged. The frozen bar was integrated for 4 kpc/(km/s), about 3.91 Gyr, as a stationary-model diagnostic, not a reconstruction of Galactic history. This single orbit remains close to the plane and cannot represent the above/below-bulge population.

## Results

The added-force difference is at most 0.160739% along the lower-resolution path and 0.113911% along the higher-resolution path. These compare two approximations to the same proposed field. The common ordinary field, outer boundary and all other approximations have not independently converged as a consequence of this check.

The lower-field solver runs differ by up to 0.010335 kpc and 0.451732 km/s over the full interval, failing the existing strict pointwise gate. Its tight-run Jacobi drift divided by 220^2 is 3.335e-10. That conserved-quantity check does not repair the trajectory gate; both results are retained.

Each interval uses exactly 2000 noninitial time samples. At the fixed illustrative 0.5-kpc / 50-km/s Gaussian-kernel resolution:

| Interval, time units | MMD between field occupations | MMD between lower-field solver occupations | MMD between higher-field solver occupations | TV between field occupations |
|---|---:|---:|---:|---:|
| (0,1] | 0.000704 | 4.90e-8 | 7.85e-8 | 0.0015 |
| (0,2] | 0.001408 | 9.57e-8 | 1.65e-7 | 0.0095 |
| (0,4] | 0.044767 | 0.000128 | 0.000832 | 0.0970 |
| (3,4] | 0.194700 | 0.000509 | 0.003278 | 0.4200 |

At 2 kpc / 100 km/s, field MMD is 0.037490 over the full four units and 0.150762 over the last unit. The corresponding higher-field solver discrepancies are 0.000355 and 0.001396. Thus the measured field-resolution sensitivity is not removed by using either of the tested solver accuracies or either fixed smoothing scale. These are empirical distances, not significance levels or percentages of wrong gravity.

## Changes in actual velocity moments

Conditional means and dispersions use fixed radial bins and separate |z| below/above 0.2 kpc. Cells need at least 50 samples in each field for a velocity comparison; all counts, including sparse and empty cells, remain in `results.json`. Samples are correlated times along one orbit, not independent observed stars.

For the full four-unit interval, the eligible 3.5–5 and 5–9 kpc near-plane cells have mean velocity component changes no larger than 0.542 km/s and dispersion changes up to 1.677 km/s. Lower-resolution sample counts are 777 and 1213 versus 698 and 1302 at higher resolution; the remaining ten lower-resolution samples occupy other cells and are not discarded from the overall distribution comparison.

In the last unit, the differences grow. In the 3.5–5 kpc cell, predicted mean rotation differs by 3.322 km/s and rotational dispersion by 7.061 km/s. Its occupancy is 1065 versus 682 samples out of 2000. In the 5–9 kpc cell, mean rotation differs by -2.882 km/s and radial dispersion by 6.894 km/s. Signs are lower minus higher resolution. Neither set is a measured stellar velocity prediction with an observational error model.

These explicit motion differences connect the abstract distribution check to quantities needed by the intended stellar likelihood. They do not say which resolution is correct or whether those differences are tolerable for any particular survey.

## Provenance and consequences

The integration, Jacobi identity, total variation, Gaussian-kernel MMD, coordinate conversion, mean and dispersion formulas are established mathematics. The gravity parameters and response form are previously fitted empirical choices, not a derived photon-capture law. No novelty is claimed for these formulas. The project contribution here is a reproducible diagnostic of the chosen implementation.

Run `python research_work/results/orbit-field-population/run.py`. Both integrations completed; result metadata retain the strict path failure. Upstream source and trajectory hashes, field hashes, solver work, force differences, all conditional cell counts and the fixed protocol are archived. No holdouts were opened or physical coefficients changed.

The next population work must include field representation as well as solver, phase and duration sensitivity across a broader paired set. Merely integrating longer or meeting a local 1% added-force criterion does not certify the target velocity distribution. Do not silently repair this by fitting gravity parameters to compensate for numerical differences. The full bulge/disk likelihood and unified redshift/companion theory remain unfinished.
