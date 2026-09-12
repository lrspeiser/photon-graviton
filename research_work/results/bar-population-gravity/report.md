# Gravity of a continuously supplied moving test population

This calculation turns the prescribed captured trajectories into potential and three-dimensional force coefficients. It is a linear response in a fixed ordinary bar, not a self-gravitating formed galaxy. The injection rate remains unspecified and unfunded. The source-direction grid is not established as converged for these new force quantities.

## Assumptions and provenance

The [population-to-gravity derivation](../bar-angular-refinement/population-to-gravity.md) gives the equations. Newtonian Green functions, Gauss quadrature and rotating-frame mechanics are known mathematics. The cold captured-particle interpretation, thin source weighting, launch speed and source history are project assumptions. No formula is claimed unique.

The retained ordinary bar component and its normalization are model-dependent. Other Galactic components, deposited self-gravity and ordinary-bar backreaction are absent. A rotating prescribed bar can exchange orbital energy with the particles; this calculation does not fund that exchange with an evolving ordinary-matter energy reservoir. The source is stationary in bar coordinates. Each of the 1 and 3 kpc launch spheres has its own unspecified constant injection rate; they are not combined into a fitted radial source.

## Results

The table is at T=0.25 kpc/(km/s), approximately 244 Myr. Multiply potential coefficients by G times the total injected rest mass of that source to obtain potential. Multiply acceleration coefficients by the same factor to obtain acceleration. Equivalently, for a constant rate, total injected mass is injection rate times T. These are predictions per unit source mass, not measured stellar accelerations.

| Source radius (kpc) | Position (kpc) | Potential / GM (1/kpc) | ax / GM (1/kpc²) | ay / GM (1/kpc²) | az / GM (1/kpc²) |
|---|---|---:|---:|---:|---:|
| 1 | (0.1, 0, 0) | -1.71874 | -0.425679 | 0.135675 | -2.45538e-18 |
| 1 | (1, 0, 0) | -1.07164 | -0.978179 | -0.00943273 | -9.4423e-19 |
| 1 | (1, 0, 1) | -0.695842 | -0.302095 | 0.000508955 | -0.36343 |
| 1 | (0, 1, 1) | -0.690556 | -0.00052629 | -0.303899 | -0.350204 |
| 1 | (3, 0, 0) | -0.336383 | -0.114188 | 2.30219e-05 | -5.04164e-21 |
| 3 | (0.1, 0, 0) | -0.609779 | -0.0319418 | 0.0182308 | 1.27949e-19 |
| 3 | (1, 0, 0) | -0.51945 | -0.111857 | -0.0425969 | 4.23037e-19 |
| 3 | (1, 0, 1) | -0.457665 | -0.0644196 | -0.00730247 | -0.0830094 |
| 3 | (0, 1, 1) | -0.4643 | -0.0114336 | -0.0591444 | -0.0695901 |
| 3 | (3, 0, 0) | -0.350095 | -0.0833628 | -0.00559343 | 7.73273e-22 |

## Verification and limits

The ring and eccentric-orbit analytic checks pass (maximum ring error 2.66e-15; inverse-radius error 0). All source orbit and age-quadrature gates pass: **True**. 0 of 108 age-quadrature comparisons remain failed. Original and refined comparisons remain in the raw data.

Age quadrature uses the integrator step intervals with 4/8 nodes and retained 16-node refinement when needed. It measures numerical time-integration error on each chosen trajectory. It does not test source-direction quadrature, field-order errors at new positions, or the continuum limit near trajectories. The earlier L40/L64 orbit agreement is useful evidence but not a whole-population force verification.

Reflection across the plane forces the vertical acceleration in the plane to zero. Above the plane, a vertical component is expected. This symmetry statement does not validate an observed Milky Way vertical force, and the listed positions are diagnostic points, not observed star bins. Bar rotation can produce an angular offset in the source distribution; individual nonzero transverse components must be checked against finer direction sampling before interpretation.

## Next calculation

Apply the same force integrals to refined direction grids before using their values in a dynamical fit. A smooth-looking central-residence estimate alone cannot establish convergence of the more locally sensitive three-dimensional force. Then evolve capture and deposited gravity together, with a funded source and ordinary-matter response. A metric or stress-source rule is also needed to predict lensing from the same population. The redshift/time mechanism and all nine overall goals remain incomplete.

Reproduce with verify_kernel.py, run.py, export.py, then document.py using the existing hashed bar cache and retained angular launches. No observational data, halo density or withheld sample is used in this checkpoint.
