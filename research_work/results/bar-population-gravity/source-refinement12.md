# Direction-sampling test of the gravitational response

This holds the source law and ordinary bar fixed and compares8x16 with12x24 angular sampling. The source is still a negligible test population, separately normalized on1 and3kpc launch spheres. No energy supply, self-gravity, observed star fit or new physical interaction is introduced.

All 208 combined orbit/age status checks pass: **True**. Of30 position/time/source comparisons, 2 fail the2 percent potential gate and 19 fail the5 percent vector-force gate. These development thresholds are not measurement error bars. A passing pair is useful evidence but not proof of convergence.

## Final-epoch comparison

At approximately244Myr, values below are per G times total injected rest mass. Potential units are1/kpc; acceleration units are1/kpc^2. Changes use the refined magnitude with the predeclared0.01 floor.

| Source R | Position | Refined potential / GM | Refined acceleration / GM | Potential change (%) | Force change (%) |
|---:|---|---:|---|---:|---:|
| 1 | (0.1, 0, 0) | -1.70891 | (-0.19761, -0.015515, 1.5679e-18) | 0.067 | 85.733 |
| 1 | (1, 0, 0) | -1.09278 | (-0.85232, -0.0022328, 8.1734e-19) | 0.326 | 26.397 |
| 1 | (1, 0, 1) | -0.696644 | (-0.30264, 0.00018131, -0.36577) | 0.170 | 0.751 |
| 1 | (0, 1, 1) | -0.690256 | (0.0003847, -0.30407, -0.35028) | 0.077 | 0.525 |
| 1 | (3, 0, 0) | -0.336564 | (-0.11437, -1.3496e-06, -1.6421e-21) | 0.021 | 0.055 |
| 3 | (0.1, 0, 0) | -0.60723 | (-0.046579, -0.013931, 6.6104e-19) | 0.190 | 110.017 |
| 3 | (1, 0, 0) | -0.52253 | (-0.11821, -0.014841, -2.5826e-19) | 0.133 | 11.014 |
| 3 | (1, 0, 1) | -0.454543 | (-0.06702, 0.010561, -0.086553) | 0.965 | 20.406 |
| 3 | (0, 1, 1) | -0.454596 | (-0.0090843, -0.052124, -0.088751) | 1.639 | 24.919 |
| 3 | (3, 0, 0) | -0.351118 | (-0.065071, -0.00078347, -2.3794e-19) | 1.011 | 40.999 |

## Interpretation and next step

Potential averages contributions over the moving population. Local force is more sensitive to where that material sits relative to a test point, so a stable potential or enclosed central fraction cannot substitute for a force check. Resolve failed locations before interpreting directional force differences as a bulge prediction. Keep actual source anisotropy distinct from sampling error.

The equations are the known Newtonian Green functions applied conditionally to the project cold-rest-mass source. Their derivation, units, age quadrature and analytic checks remain in [the original report](report.md). The refined force calculations do not include a fresh L40/L64 comparison, full ordinary Galactic components or collective evolution.

All raw coarse and refined outputs are retained. The next required step is further force-source refinement or a demonstrably convergent integration method, followed by source normalization and evolving capture/self-gravity. A metric/stress law is needed for lensing. Joint redshift/timing predictions and all nine overall goals remain incomplete.
