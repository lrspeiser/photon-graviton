# PF-0: local phase-feedback source-budget screen

Declared20 September2026 before implementation/calculation. This tests a proposed modification from next-candidate-requirements.md without changing the running CC-2 equations. It is a homogeneous single-amplitude diagnostic, not a 3D swirl, galaxy or lensing solution.

The candidate lapse contains U=-lambda |A|^2/2 in a scalar-decoupled local limit. With spatial-response parameter s=1, field kinetic coefficient is exp(4U). Rescale B=sqrt(lambda) A, P=sqrt(lambda) Pi/Omega, time tau=Omega*time and energy h=lambda H/Omega^2; let t=lambda rho/Omega^2. The reduced Hamiltonian is

h=.5 exp(-2B^2) P^2 + .5 B^2 + t exp(-B^2/2).
Bdot=exp(-2B^2)P.
Pdot=2B exp(-2B^2)P^2-B+t B exp(-B^2/2)-gamma Bdot.
Qdot=gamma Bdot^2; h+Q is conserved.

All energy terms are nonnegative. Matter rest-energy is the final exponential term. The seed B=1e-8,P=0 has a specified initial energy; this is not formation from exact zero. A zero seed stays zero. Gamma is an explicit phenomenological energy outlet, not a demonstration of physical outgoing radiation. Scalar response, source motion, orientation and spatial gradients are omitted in this diagnostic.

Analytic stationary branches: B=0, and B_*^2=2 log(t) for t>1. Curvature at zero is1-t; curvature at B_* is2 log(t). At the condensed branch alpha=1/t, matter-energy fraction=1/t, field-potential fraction=log(t)/t, and residual releasable energy fraction=1-(1+log(t))/t. A weak |U| band[1e-8,1e-5] occupies t in[exp(1e-8),exp(1e-5)]. This is a local sensitivity statement, not a measured galactic density interval.

Numerical checks: finite differences of h in B,P against conservative equations at300 seeded states with B uniform[-2,2], P uniform[-1,1], t logarithmically uniform[.1,10], tolerance2e-5 scaled. Check damped h+Q derivative at gamma=.05, tolerance2e-7. Root and curvature checks for t={.5,.99,1,1+1e-6,1.01,1.1,2,10}, using numerical minimization/root methods independent of the explicit B_* formula. Preserve the critical zero-curvature case rather than call it strictly stable.

Execute those8 t values at gamma=0 and.05:16 trajectories, seed as above. T=min(10000,max(200,40/max(r_growth,1e-3))), where r_growth=(sqrt(gamma^2+4(t-1))-gamma)/2 for t>1 and0 otherwise, computed stably near threshold. Use solve_ivp Radau, rtol1e-9,atol1e-12,500 output times. Scaled energy-ledger drift<=1e-6. Record amplitude, matter/field energy exchange and whether |B| reaches half the analytic nonzero branch; horizon-limited nonformation is not an exclusion. The gamma=0 trajectories need not settle at a minimum; do not confuse oscillation with asymptotic stability. Check exact-zero seed separately.

Credit spontaneous scalarization/vectorization precedents (Damour-Esposito-Farese gr-qc/9602056; Ramazanoglu1706.01056). This is an attributed candidate adaptation; no historical novelty claim, no imported neutron-star solution, no dark matter or expansion. Preserve protocol/source hashes and first results. Do not claim this reduced potential validates the coupled3D extension.
