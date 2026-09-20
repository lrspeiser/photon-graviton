# SE-H: nonlinear homogeneous perturbation diagnostic

Declared 20 September 2026 before execution. Positive field energy and frozen
propagation speeds do not prove nonlinear stability. Test the exact homogeneous
phi,X,Y sector of SE-1/SE-2 with A=0, no particles and no spatial gradients.
This sector has no field current and consistently retains A=0. Spatially varying
perturbations generally couple to A and are explicitly outside this test.

For q=(phi,X,Y), p their canonical momenta, define U=.08 phi, a=exp(4U),
k=k0 tanh(U/.001), m2=.75^2 exp(2 chi U), r=Y-kX. The Hamiltonian is

    H = a |p|^2/2 + .2^2 phi^2/2 + m2 r^2/2.

Derive qdot=a p and pdot_phi=-2g a |p|^2-.2^2 phi
-g m2(chi r^2-k_U Xr), pdot_X=m2 k r, pdot_Y=-m2 r.
Evolve the nonlinear background and its 6x6 tangent matrix simultaneously,
Mdot=J_rhs M, M(0)=I. Evaluate J by complex-step differentiation at 1e-25.

Twelve fixtures: chi=0,50,200; k0=0,.5; initial X=.01 or Y=.01 with all
other coordinates and momenta zero. Evolve to T20 with DOP853, rtol1e-10,
atol1e-12, 201 recorded times. Report energy drift, maximum tangent singular
value over time, endpoint singular values and log(sigma_max)/T. This last
quantity is a finite-time coordinate-dependent amplification rate, not an
asymptotic Lyapunov exponent. A free coordinate already has secular growth.
No amplification threshold is a physical rejection gate.

Numerical gates: absolute energy drift <1e-11 +1e-8 |H0|; Hamiltonian tangent
symplectic residual max|M^T J M-J| <1e-6 at the endpoint. Independently check
each tangent map using one seeded unit direction (seed20260925), two nonlinear
shadow solutions displaced by +/-1e-7 and +/-5e-8, central difference, requiring
relative Euclidean discrepancy <1e-4. Preserve any failures and do not adjust
their thresholds. Zero-energy equilibria are allowed and reported explicitly.

Audit analytic RHS with direct Hamiltonian directional differences on 48 seeded
states, error <1e-7 absolute; check tangent generator satisfies J_rhs^T J+
J J_rhs=0 to 1e-10. Attribute Hamiltonian/tangent methods to established math.
This tests the candidate constitutive law, not agreement with earlier gravity.
It cannot establish stability of emitted 3D fields, long-lived swirl or any
observational prediction. Preserve code/protocol hashes and all twelve results.
