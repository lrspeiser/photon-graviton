# NK-1: a kernel outside the positive Gaussian mixtures, shell footprints
20 September 2026. [Protocol](protocol-nk1.md) efb6c2f, declared before any calculation; [amendment 1](protocol-nk1-amendment-1.md) a59f3d4 (after the first complete run, preserved as `nk1-results-as-declared.json`: both galaxy solves finished on their free face before the certificate is read; gate K2 stays as declared and failed). Archive `nk1-results.json`, runtime 33 s; suite job `nk1_checks.py`. Inputs: CL-2 stage 2's galaxy cache and archive.

**The question, in plain language.** Stage 2 exhausted the fields that are a positive blur of the matter. The one linear family left outside is a footprint written at a distance from each mass element, a shell of radius d_0 and thickness w. NK-1 asks, galaxies first, whether adding 18 such shells (d_0 from 1 to 50 kpc, w from a tenth to a half of d_0) to the 25 Gaussian widths lets a linear field describe the galaxies.

**Outcome.** b: no nonnegative combination of Gaussian and shell footprints up to 50 kpc describes the galaxies; certified minimum 31.28 km/s against 21.9.

| status | finding |
|---|---|
| **reproduction** | **passes.** With every shell amplitude forced to zero the certified minimum equals stage 2's to 4.9e-13. |
| **numerical verification** | **FAILED on K2_quadrature, C1_solvers, C2_kkt; every other gate passes.** The shell ring average and its derivative reproduce the library's Gaussian ring kernel at zero offset to 3.3e-15 (the 5-kpc offset differs by 97%, the control); the spherical shell average to 3.0e-15; source-grid doubling changes the columns by 7.3e-04. **The primary solve (reconstructed source) is certified:** its two convex solvers agree to 9.1e-14 with an optimality residual of 4.7e-16 and positive directional curvature at the optimum and at random feasible points. **The tabulated-force sensitivity is not certified:** its solvers agree only to 2.7e-05 with an optimality residual of 1.9e-02, so gates C1 and C2, which cover both solves, are recorded as failed; its 30.39 km/s is a bound from above on that sensitivity's minimum, not a certified value, and the reading of the decision rests on the primary alone. Gate K2 required the shell columns to change by less than 1e-8 between 128 and 256 quadrature nodes; they change by 4.9e-06 on the thinnest shells, recorded as failed and not re-thresholded (the shells' whole contribution to the certified minimum is 0.05 km/s against a margin of 9.4 km/s). |
| **scientific outcome** | **b: no nonnegative combination of Gaussian and shell footprints up to 50 kpc describes the galaxies; certified minimum 31.28 km/s against 21.9** |

## The galaxies
| source | certified minimum, training | validation | test | active shells (d_0, w) in kpc | active Gaussian widths | outer slope model / observed |
|---|---|---|---|---|---|---|
| reconstructed (primary) | **31.28 km/s** | 34.66 | 29.07 | (1, 0.1), (1, 0.25), (2, 1), (5, 0.5), (5, 2.5), (10, 5), (20, 2), (20, 5), (50, 5), (50, 12.5) | 6 | 0.406 / 0.292 |
| tabulated (sensitivity) | 30.39 km/s | 33.81 | 28.50 | (1, 0.1), (1, 0.25), (1, 0.5), (2, 0.2), (2, 0.5), (2, 1), (5, 0.5), (5, 1.25), (5, 2.5), (10, 1), (10, 2.5), (10, 5), (20, 2), (20, 5), (20, 10), (50, 5), (50, 12.5), (50, 25) | 25 | 0.399 / 0.292 |

Stage 2's certified minimum on the Gaussian family alone was 31.33 km/s; the shells improve it by 0.05 km/s. Part 2 (clusters and lenses) was not run, as the protocol declares for decision (b).

## What is not shown
Shells beyond 50 kpc and thinner than a tenth of their radius; the galaxy sources are stage 2's razor-thin reconstruction; exposed data.

## Next decision
The linear route is closed: no nonnegative combination of Gaussian and shell footprints up to 50 kpc describes the galaxies with these sources, a certified negative that extends stage 2 beyond the positive Gaussian mixtures. A response not proportional to its source remains the only route left open by the galaxies.
