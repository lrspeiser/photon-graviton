# SE-P: common matter/light probe equations

Declared 20 September 2026 before implementation. Derive probe trajectories
from the same regularized particle Hamiltonian used in SE-1/SE-2, taking the
probe-energy limit so the finite source evolves independently. This is not a
new force fit or an inserted light-bending multiplier.

For a probe with rest parameter m>=0 and canonical momentum p, use

    H(q,p,t) = sum_grid W(q) [alpha sqrt(m^2+z |p|^2) + beta dot p],
    qdot = sum W [alpha z p/sqrt(m^2+z |p|^2) + beta],
    pdot = -sum (partial W/partial q) [alpha sqrt(m^2+z |p|^2)+beta dot p].

Here alpha=exp(U), z=exp(2U), U=.08 phi and beta is the already declared bounded
vector shift. W is the normalized spherical cubic source kernel, with its exact
normalization derivative. Light uses m=0 and nonzero p. Its speed relative to
the averaged shift is exactly the averaged front speed Cbar=sum W exp(2U).
Massive probes are bounded by Cbar. X/Y influence probes through evolved phi,A;
no additional attachment force is inserted. Electromagnetic polarization and
literal graviton attachment remain outside this model.

First validate 72 seeded local fixtures (seed20260926), cycling m=0,.1,1,
probe radii .9,1.2 and grids n24,32 in L16, using smooth nonuniform phi,A and
random q,p. Directional Hamiltonian derivatives in q,p must agree with the
Hamiltonian vector field within 1e-7 absolute, with difference step1e-5.
Require cone excess<1e-12 and massless equality within1e-12.
Massless p scaling by .5 and 2 must leave qdot unchanged and scale pdot,H
linearly within1e-12. Check four uniform-field trajectories against their
constant Hamiltonian velocities at T4 using RK4 dt.02, position error<1e-10,
momentum change<1e-12. This is an analytic control of our candidate equations.

Next evolution stage (must receive its own run declaration before calculation):
integrate probes and source field at the same RK stages, trace a transverse
bundle and massive test bodies, compare finite bundle derivatives and steps,
and test probe-radius dependence. Use the evolving field, not a frozen lens
snapshot. A fixed-time bundle map is not automatically an observed lens-plane
map or time delay; any detector convention must be declared. No observed galaxy
or cluster prediction is claimed from this implementation alone.

Hamiltonian ray tracing, test-body limits and numerical derivatives are
established methods. The shared particle/field response is a candidate choice,
not a consequence of naming the field a graviton or an old-theory pass target.
