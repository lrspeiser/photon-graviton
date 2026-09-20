# SR-2: resolve spatial-response ray accuracy

Declared 20 September 2026, before this campaign's calculations. SR-1 measured an 11.17% bend difference between dx=.5 and .375. Preserve that evidence. The continuum Hamiltonian, s=1, lambda=0, g=.08, eta=.08, kappa=.5, omega=.2, source radius .9, initial state and T=2.5 remain unchanged. No observations or extra optical factors are introduced.

## Equivalent compact deposition

Compute the same normalized spherical W=max(1-|x-q|^2/radius^2,0)^3 and its exact q derivative only on its nonzero bounding box. Zeros outside that box contribute nothing. Preserve self-source terms and deposition/interpolation adjointness. Reject sources whose support crosses the box edge rather than silently wrapping/truncating them. Test 36 seeded fixtures (seed 20260921), including s=0/1, eta=0/.08 and n=12/24, against the unmodified SR-1 implementation: maximum scaled RHS, energy and cone differences <1e-11. No lambda feedback is supported in the optimized implementation. This is computational optimization, not altered physics.

## Higher-order positive discrete Hamiltonian

Let D4_j f=(-f(x+2h)+8f(x+h)-8f(x-h)+f(x-2h))/(12h), a skew-adjoint fourth-order centered derivative. Let T_j=(Dminus_j Dplus_j)^2, a self-adjoint fourth derivative. Use sigma=1/64 and

G4 = sum_j,a (D4_j F_a)^2/2 + sigma h^6 sum_j,a (T_j F_a)^2/2,

J4_j = sum_a Pi_a D4_j F_a.

Replace only SR-1's discrete gradient energy/current by G4,J4. Keep its coefficients a,d,beta and all their derivatives. Thus the spatial canonical terms are sum_j[D4_j(d D4_j F)-sigma h^6 T_j(d T_j F)-D4_j(beta_j Pi)]. Fdot=a Pi-beta.D4F. Positivity of the main kinetic/gradient block follows from |beta|<sqrt(a d); the added stabilizer is nonnegative and removes centered-derivative checkerboard zero stiffness. It vanishes at order h^6 for a smooth fixed field and does not alter the continuum equations. For constant coefficients its squared frequency is a d [sum_j sin(theta_j)^2(4-cos(theta_j))^2/(9h^2)+sigma h^6 sum_j(16 sin(theta_j/2)^4/h^4)^2], plus the nonnegative mass contribution. Taylor consistency and this Fourier symbol must be tested. This is a numerical regulator, not an invented gravitational law.

For both second/fourth discretizations, test 24 random full-state directional Hamiltonian derivatives per order (n=12) against central differences, error <2e-6. Test nonzero Fourier stiffness on all nonzero n=16 3D modes, and fourth-order derivative convergence on fixed smooth modes n=32/64/128 (ratio >12). Do not infer nonlinear physical stability solely from these checks.

## Declared coupled runs

All runs start with empty fields and six rotating matter sources of momentum/mass .2 and one photon of momentum 1e-4. Source support radius stays .9. Match old source placement; use zero rotation except the explicit rotation case. In order:

1. Second order n=28,L=14,dt=.02: reproduce SR-1 enlarged-box fast case.
2. Fourth order n=28,L=14,dt=.02.
3. Fourth order n=40,L=15,dt=.01.
4. Fourth order n=56,L=14,dt=.01.
5. Fourth order n=70,L=14,dt=.01.
6. Second order n=56,L=14,dt=.01.
7. Second order n=70,L=14,dt=.01.
8. Fourth order n=56,L=14,dt=.005.
9. Fourth order n=56,L=14,dt=.01, rotated pi/3 about (1,2,3).

Save initial/final raw states and metrics every .05 model time where compatible (otherwise every integer number of steps closest to .05, always include endpoints). Track maximum sampled characteristic speed at every RHS evaluation and source extent at every completed step. Record edge amplitude at metric samples. Gates: scaled H+Q drift <1e-5, averaged photon-cone residual <1e-10, edge amplitude <1e-5, positive measured-speed geometric clearance to damping entrance. The higher stencil has wider numerical coupling: no exact discrete causal-cone or absorbing-boundary claim follows from that clearance. Fixed-resolution original/enlarged SR-1 agreement does not automatically validate the new stencil; report this limitation.

Accuracy claims require more than these consistency gates. Require the finest two fourth-order bend magnitudes to agree within 1%, a difference smaller than the preceding refinement difference, and second/fourth order at the finest shared spacing to agree within 1%. Require time-refined bend difference <0.1% and rotation-back-transformed probe velocity direction difference <1% of the unrotated bend. Report failures rather than weaken thresholds. No Richardson extrapolation is treated as measured ground truth. Rotated bends must be computed relative to the rotated initial direction; absolute atan2(vy,vx) is not the bend.

Independently reconstruct final energies with a separate stencil implementation, verify hashes and archived decisions, and publish signed/unoriented bend definitions. Compare source energy loss, angular residual and wall/CPU performance. Archive all first runs. Further grids, box tests, or changed methods need a new declared amendment if these nine runs do not establish the target precision.

## Attribution and scope

Centered finite differences, summation by parts on a periodic grid, Hamiltonian adjoints, Fourier stability analysis and positive high-order regularization are established numerical mathematics. The specific combination is a numerical implementation choice, not a claim of a new physical theory. Metric optics and the effective field assumptions retain SR-1's attribution. Physical source-size justification, local source microphysics, long-time open boundaries, affordable outer support, galaxy/cluster observations and all other twelve-goal requirements remain incomplete.
