# SE-H: finite-time homogeneous perturbations

All twelve trajectories and 48 Hamiltonian/generator controls pass the declared
numerical gates. The independent raw-archive audit passes 136 checks.

| chi | mixing k0 | Initially excited channel | Maximum tangent amplification |
|---:|---:|---|---:|
| 0 | 0 | X | 20.049876 |
| 0 | 0 | Y | 20.04902 |
| 0 | 0.5 | X | 20.049876 |
| 0 | 0.5 | Y | 23.189784 |
| 50 | 0 | X | 20.049876 |
| 50 | 0 | Y | 20.027947 |
| 50 | 0.5 | X | 20.049876 |
| 50 | 0.5 | Y | 14.878425 |
| 200 | 0 | X | 20.049876 |
| 200 | 0 | Y | 19.974575 |
| 200 | 0.5 | X | 20.049876 |
| 200 | 0.5 | Y | 21.719401 |

The amplification is a Euclidean singular value in the declared canonical
coordinates over 20 model time units. It is coordinate dependent. The uncoupled
free X coordinate has shear map [[1,T],[0,1]] and singular value
(sqrt(T^2+4)+T)/2 = 20.0498756 at T=20. Thus amplification around 20 is
not by itself evidence of exponential instability or of enhanced gravity.

The interaction changes perturbation evolution, but this horizon and these
initial states do not establish an asymptotic growth rate or nonlinear stability.
The calculation evolves the homogeneous scalar sector with the vector field
zero. Spatial perturbations couple to the vector sector and were not tested.
There are no particles, emitter or boundary here. This does not validate the
stability of the emitted 3D configurations or a permanent swirl.

Maximum absolute energy drift: 1.513591294942558e-14.
Maximum endpoint symplectic residual: 1.2280361934667184e-10.
Maximum shadow-trajectory relative discrepancy: 7.206564088578388e-05.

The tangent map is checked against independently evolved +/- perturbations at
two amplitudes. The archive audit reconstructs energies and tangent diagnostics;
it does not independently re-integrate the shadow trajectories.

Protocol and code were committed before execution. Raw trajectories, tangent
maps, controls, source hashes and all twelve results are preserved in
homogeneous-stability-v1. Run audit_homogeneous.py to regenerate this report.
Hamiltonian, symplectic, tangent and numerical integration methods are established
mathematics. The field constitutive laws are candidate assumptions; no physical
novelty or observational success is claimed. All twelve goal items remain active.
