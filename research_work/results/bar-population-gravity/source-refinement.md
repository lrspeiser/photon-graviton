# Direction-sampling test of the gravitational response

This holds the source law and ordinary bar fixed and compares6x12 with8x16 angular sampling. The source is still a negligible test population, separately normalized on1 and3kpc launch spheres. No energy supply, self-gravity, observed star fit or new physical interaction is introduced.

All 100 combined orbit/age status checks pass: **True**. Of30 position/time/source comparisons, 0 fail the2 percent potential gate and 21 fail the5 percent vector-force gate. These development thresholds are not measurement error bars. A passing pair is useful evidence but not proof of convergence.

## Final-epoch comparison

At approximately244Myr, values below are per G times total injected rest mass. Potential units are1/kpc; acceleration units are1/kpc². Changes use the refined magnitude with the predeclared0.01 floor.

| Source R | Position | Refined potential / GM | Refined acceleration / GM | Potential change (%) | Force change (%) |
|---:|---|---:|---|---:|---:|
| 1 | (0.1, 0, 0) | -1.70776 | (-0.34083, 0.075951, 8.4898e-18) | 0.643 | 29.714 |
| 1 | (1, 0, 0) | -1.08922 | (-1.0705, -0.057314, 1.7461e-18) | 1.614 | 9.698 |
| 1 | (1, 0, 1) | -0.69783 | (-0.30563, 0.00027518, -0.36771) | 0.285 | 1.162 |
| 1 | (0, 1, 1) | -0.689722 | (0.0004449, -0.30399, -0.34785) | 0.121 | 0.552 |
| 1 | (3, 0, 0) | -0.336633 | (-0.11443, 1.0873e-05, -5.2596e-22) | 0.074 | 0.213 |
| 3 | (0.1, 0, 0) | -0.606078 | (-0.0033716, 0.017598, 7.3009e-20) | 0.611 | 159.490 |
| 3 | (1, 0, 0) | -0.521833 | (-0.12548, -0.025762, 1.3669e-19) | 0.457 | 16.906 |
| 3 | (1, 0, 1) | -0.458931 | (-0.081118, -0.0036867, -0.076461) | 0.276 | 16.405 |
| 3 | (0, 1, 1) | -0.462046 | (-0.011611, -0.064794, -0.066479) | 0.488 | 6.897 |
| 3 | (3, 0, 0) | -0.347568 | (-0.091051, 0.0052884, 5.3745e-19) | 0.727 | 14.609 |

## Interpretation and next step

Potential averages contributions over the moving population. Local force is more sensitive to where that material sits relative to a test point, so a stable potential or enclosed central fraction cannot substitute for a force check. Resolve failed locations before interpreting directional force differences as a bulge prediction. Keep actual source anisotropy distinct from sampling error.

The equations are the known Newtonian Green functions applied conditionally to the project cold-rest-mass source. Their derivation, units, age quadrature and analytic checks remain in [the original report](report.md). The refined force calculations do not include a fresh L40/L64 comparison, full ordinary Galactic components or collective evolution.

All raw coarse and refined outputs are retained. The next required step is further force-source refinement or a demonstrably convergent integration method, followed by source normalization and evolving capture/self-gravity. A metric/stress law is needed for lensing. Joint redshift/timing predictions and all nine overall goals remain incomplete.
