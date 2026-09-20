# SE-PI: shrinking numerical interpolation for the point-probe limit

Declared 20 September 2026 before calculation. SE-B uses a fixed physical
averaging radius; changing that radius changes the model response. Its26.59%
bend sensitivity is not something a finer time step can remove. Introduce a
separate numerical approximation to the same continuum point Hamiltonian,
keeping the matter source's physical spherical radius unchanged.

Use tensor-product cardinal cubic B-spline weights on the four grid nodes
around each probe coordinate. They are nonnegative, sum to1 and have support
2h per axis; their width shrinks with grid spacing h. Do not deconvolve them
or introduce a free photon-size parameter. Evaluate H_h=sum W_h h_particle
and differentiate these same weights for force. The averaged cone bound holds
because the weights are positive; pointwise accuracy must be established by
grid convergence. Cartesian interpolation can still introduce finite-grid
directional errors and is not claimed rotationally invariant at finite h.

For fractional coordinate t in[0,1], weights at nodes i-1,i,i+1,i+2 are
[(1-t)^3, 3t^3-6t^2+4, -3t^3+3t^2+3t+1, t^3]/6.
Differentiate analytically with respect to q, dividing by h. H uses the
existing alpha,z,beta at those nodes, not separately interpolated gradients.

Controls: 72 seeded fixtures seed20261001, n24/32,L16,m=0,.1,1, smooth analytic
fields as SE-P and random q in[-2,2]^3,p normal. Hamiltonian directional error
at step1e-5 <1e-7, cone excess/light equality<1e-12. Partition/derivative sums
within1e-12. Verify uniform-field force<1e-12. All derivative and cone checks
must use this new interpolation's own H, not the fixed-radius implementation.

Continuum consistency: analytic fields phi=-.04 exp(-r^2/8),
A=.1 exp(-r^2/8)(-y,x,0), q=(.23,.37,.19),(1.13,-.67,.41),(-1.71,.83,-.29),
p=(.7,.2,-.1), m=0,.1,1. Evaluate velocity/force at n24,48,96,L16. Compare to
the continuum candidate point Hamiltonian, with its position derivative by
complex-step1e-25. Require combined Euclidean RHS error at n96 <n24 error/8
(second-order trend with margin), preserving all coarse/fine errors. This
checks numerical consistency with our declared law, not old gravity.

No new trajectory claim follows until evolving bundle runs verify space/time
convergence and compare with fixed-radius results. The physical source radius
and its earlier sensitivity remain separate requirements. Splines and
Hamiltonian interpolation are established numerical methods, not invented
gravity formulas. All twelve original research requirements remain active.
