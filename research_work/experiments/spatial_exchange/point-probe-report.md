# SE-PI: consistent shrinking interpolation for point probes

All 72 Hamiltonian/cone/partition controls and nine continuum consistency fixtures pass.

| Position | Mass | n24 error | n48 error | n96 error |
|---|---:|---:|---:|---:|
| [0.23, 0.37, 0.19] | 0.0 | 0.00048102716 | 0.00012511503 | 3.1530332e-05 |
| [0.23, 0.37, 0.19] | 0.1 | 0.0004809911 | 0.00012510581 | 3.1528005e-05 |
| [0.23, 0.37, 0.19] | 1.0 | 0.00041971381 | 0.00010946569 | 2.7580321e-05 |
| [1.13, -0.67, 0.41] | 0.0 | 0.00029887715 | 7.7500612e-05 | 1.9546685e-05 |
| [1.13, -0.67, 0.41] | 0.1 | 0.00029887994 | 7.7501414e-05 | 1.9546893e-05 |
| [1.13, -0.67, 0.41] | 1.0 | 0.00031111883 | 8.0776352e-05 | 2.0380097e-05 |
| [-1.71, 0.83, -0.29] | 0.0 | 0.00048531212 | 0.00012525677 | 3.1573109e-05 |
| [-1.71, 0.83, -0.29] | 0.1 | 0.00048529862 | 0.00012525336 | 3.1572252e-05 |
| [-1.71, 0.83, -0.29] | 1.0 | 0.00046473147 | 0.00012008994 | 3.0277544e-05 |

Reducing grid spacing by4 reduces RHS error by at least15.2179, consistent
with second-order convergence on these analytic fields. The numerical kernel
shrinks with the grid and has no independently fitted photon radius. Positive
weights retain the Hamiltonian averaged-cone bound. This does not imply exact
rotation invariance at finite grid spacing or converged bending in live fields.

The existing fixed-radius probe campaign is preserved. Changing probe averaging
there changed the effective response; refining the new interpolation instead
approaches the declared continuum point Hamiltonian while leaving the physical
matter-source radius unchanged. Source-radius sensitivity remains unresolved.

Seven SE-PB evolving runs are declared separately, comparing emission-only and
direct-companion cases at matched n32,48,64, plus a time control. Their results
must establish both absolute path accuracy and the much smaller paired emitter
difference. No stronger lensing or observational solution is claimed here.

Spline interpolation and Hamiltonian differentiation are established numerical
methods. The particle response is the candidate law. Raw controls, source
hashes and commit pins are in point-probe-controls-v1. All twelve goals remain active.
