# Shared wave propagation and exploratory curvature comparison

This package accompanies Appendix F of temporal_redshift_paper.docx. It retains both curvature signs and introduces shared photon/gravitational-wave propagation as a candidate requirement. No rotating-geometry or CMB likelihood fit has been performed.

To reproduce, extract the archive with both sibling directories intact and run:

```bash
OPENBLAS_NUM_THREADS=1 python shared_curved_update/curvature_test.py
```

Requires Python, NumPy, SciPy and pandas. The script imports the earlier brightness calculation from `shared_interaction_test`, reuses its full Pantheon+ covariance and Cepheid-only calibration, and fits one curvature radius for each sign. The existing kappa remains 0.000077315 per million light-years. Importing the earlier calculation regenerates its output files identically.

These are exploratory fits to the same previously inspected 960 supernovae, not held-out validation. The positive-curvature scan stops before the first conjugate point and chooses the flat boundary. Its chi-squared remains 898.4746. Negative curvature gives radius 4547.07 Mpc and chi-squared 873.0413, improving by 25.4333 with one added parameter. It does not establish a static background solution and does not fully remove the brightness offset. The original fixed reference statistic is 838.4834.

The negative branch scans a=1/(kappa*R_c) from 0 to 10; the positive branch scans from 0 to 0.999*pi/max(log(1+zHD)). Bounded minimization is compared with the flat boundary and a 201-point coarse scan. These are numerical search ranges, not physical priors. Boundary optima do not justify ordinary unconstrained likelihood-ratio significance claims.

Curvature is an isotropic beam-area change in this test. Gödel-style rotation is an additional property requiring directional ray propagation and polarization predictions. Published CMB bounds for expanding Bianchi models are not directly transferred to this fictional nonexpanding model. Mean CMB temperature does not identify curvature. The paper records a map-level temperature/polarization test plan, including radiation-source assumptions, axis-selection effects, and the limitations of previously used bandpowers.

Data provenance, checksums and the original fixed-law protocol are in the included `shared_interaction_test` directory. `curvature_results.json` records the new fits and their limitations. No coefficient, absolute magnitude, opacity or luminosity-evolution adjustment was introduced in this curvature calculation.
