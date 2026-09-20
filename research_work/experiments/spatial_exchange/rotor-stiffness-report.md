# SE-Rot3: frozen-source stiffness result

All eight declared cases pass their spectral residual, symmetry and bound checks.

| Grid | Coupling | Largest computed eigenvalue | Conservative upper bound | Half RK4 timestep bound |
|---:|---:|---:|---:|---:|
| 12 | 0.0 | 27.04 | 27.04 | 0.27196415 |
| 12 | 2.0 | 27.04 | 27.107216 | 0.27162675 |
| 12 | 10.0 | 27.04 | 28.720402 | 0.26388806 |
| 12 | 100.0 | 27.04 | 195.08023 | 0.10125311 |
| 24 | 0.0 | 108.04 | 108.04 | 0.13605757 |
| 24 | 2.0 | 108.04 | 108.3032 | 0.13589215 |
| 24 | 10.0 | 108.04 | 114.61994 | 0.13209455 |
| 24 | 100.0 | 108.04 | 766.03394 | 0.051096481 |

The added curl stiffness is a nonnegative quadratic form for fixed source Q.
The uniform oscillator term gives a strictly positive lower bound in this
frozen subsystem. The largest computed eigenvalue remains the discrete
checkerboard value here: centered curl vanishes on that grid mode. This does
not mean that the coupling leaves the rest of the spectrum unchanged.

The conservative timestep bound applies only to the declared frozen vector
system. It cannot be adopted as a complete nonlinear stability or causal-speed
guarantee. Dynamic rotor coordinates, source motion, scalar reaction and the
finite-kernel nonlocal response still require tests. In particular, a vacuum
field-front speed does not establish propagation inside the coupled source.

Code was committed before execution at2737a48. Raw eigenvectors, all values and
source hashes are preserved in rotor-stiffness-v1. Quadratic-form bounds,
eigensolvers and RK4 stability are established mathematics, not novel gravity.
All twelve original requirements remain active.
