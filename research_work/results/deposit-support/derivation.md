# A lasting galactic reservoir needs mechanical support

## Why we need this

The earlier point-source calculation gave a promising density shape if captured energy stayed where it was deposited: energy density proportional to 1/r^2, which can give a flat circular-speed contribution under ordinary gravity. But a density map is not a mechanism that keeps the energy there. A stored excitation moves with its receiver; free massive deposits accelerate; radiation has pressure and can escape.

We therefore test the motion and stress needed to preserve that shape. Without this step, multiplying local capture power by billions of years could put extra gravity in places where the energy no longer resides. This pass compares optional material/orbital and radiation descriptions. It does not adopt new matter, dark matter, expansion, a cosmic chronology or a special void-time law.

## 1. A declared diagnostic potential

On a bounded radial interval, choose the Newtonian potential Phi(r)=v_c^2 ln(r/r_ref), so the radial acceleration is -v_c^2/r and circular speed is constant. Here v_c is a specified diagnostic velocity, not a mass fitted from an assumed dark-matter halo. Any completed photon-supplied model must derive its potential and source energy; this calculation asks what support would be required **if** the desired flat-speed potential were present.

The singular center and infinite outer extent of the logarithmic form are not accepted physical predictions. Free-fall tests stop at one-tenth the starting radius. Numerical examples use v_c=200 km/s and an initial radius of 10 kpc as illustrations, not new observations.

## 2. Deposits at rest do not remain at rest

A free massive receiver released at r0 with zero velocity obeys

    d2r/dt2 = -v_c^2/r
    v_r^2 = 2 v_c^2 ln(r0/r)
    t(r0 -> r) = (r0/v_c) sqrt(pi/2) erf[sqrt(ln(r0/r))].

The code integrates the orbit and checks both this fall time and conservation of kinetic plus potential energy. In the example, going from 10 kpc to 1 kpc takes about 59 million years. This is much shorter than the illustrative multibillion-year accumulation periods used in previous budgets. It is a test-particle result in a fixed potential, not a many-body collapse simulation.

Consequently, permanently retaining the energy internally is insufficient to keep its position fixed. Existing bound matter can supply its existing orbital motion and stresses, but the captured energy then follows that matter's changing distribution. It cannot be counted again as an independently stationary reservoir.

## 3. What random orbital support would require

For a spherical collisionless receiver distribution with no mean radial flow, the radial moment of the collisionless Boltzmann equation is

    d(rho sigma_r^2)/dr + 2 beta rho sigma_r^2/r = -rho dPhi/dr
    beta = 1 - (sigma_theta^2+sigma_phi^2)/(2 sigma_r^2).

The [spherical Jeans analysis of Diakogiannis et al.](https://arxiv.org/html/1406.2542) provides this general dynamical framework. We use the equation, not their fitted mass distributions or any dark-matter premise. The solutions and diagnostics below are our conditional derivations.

For rho proportional to r^-2, constant beta<1 and radius-independent dispersion, the equation gives

    sigma_r^2 = v_c^2/[2(1-beta)].

With isotropic motions, sigma_r=v_c/sqrt(2), about 141 km/s at the illustrative circular speed. Its mean orbital kinetic energy divided by receiver rest energy is 3 v_c^2/(4 c^2), about 3.34e-7. Thus orbital support need not consume a large fraction of rest energy, but its velocities and angular momenta must have an origin. This small fraction does not solve the much larger stored-rest-energy supply problem.

For an outer radius R with zero radial stress there, the moment solution instead is

    sigma_r^2(r) = v_c^2 [1-(r/R)^(2(1-beta))]/[2(1-beta)], beta != 1
    sigma_r^2(r) = v_c^2 ln(R/r), beta=1.

The second expression matters: radial orbits are not universally forbidden. Boundary conditions and varying dispersion change the result. Numerical integration checks both branches for five anisotropies. A nonnegative dispersion satisfying one moment equation is not proof of a positive distribution function or dynamical stability.

## 4. An explicit positive equilibrium and the missing formation process

In the fixed logarithmic potential an isotropic distribution

    f(x,v) = A exp[-(v^2/2+Phi(r))/sigma^2]

is positive and satisfies the source-free collisionless equation, because it depends only on the conserved orbital energy. Integrating over velocity yields rho proportional to exp(-Phi/sigma^2). Setting sigma^2=v_c^2/2 therefore gives rho proportional to r^-2.

The code independently checks the velocity integrals and the collisionless equation at 100 random phase-space points. This is an explicit formal equilibrium, beyond a moment-only check. It still has an unphysical scale-free center and outer extent if extrapolated; a finite receiver population, distribution cutoff and self-consistent potential must replace it. No formation or stability claim is made.

Continuous deposition needs the equation with its actual source:

    df/dt + v.grad_x f - grad_x Phi.grad_v f = S_capture - S_loss.

Injecting receivers at rest supplies a sharply different velocity distribution from the equilibrium above. In a fixed potential, simply increasing its normalization would require a source with the same phase-space shape, or a demonstrated relaxation/transport process. When the reservoir contributes significantly to gravity, Phi changes too, so this fixed-potential solution cannot be scaled up indefinitely.

Purely radial incoming companions from a central source carry zero angular momentum about it. Spherical forces preserve that angular momentum. They cannot by themselves give each newly created receiver the transverse orbital motion of the isotropic solution. Existing orbiting absorbers, nonspherical torques or interactions could provide it, but their momentum and energy must be included. A separate circular-orbit integration confirms that orbital support works when the necessary angular momentum is supplied; it does not generate that momentum from nothing.

## 5. Ordinary isotropic radiation is a different store

For isotropic massless radiation the pressure is P=u/3. In a static spherical metric with lapse N(r), local stress conservation gives

    dP/dr = -(u+P) dlnN/dr
    dlnu/dlnr = -4 dlnN/dlnr.

This is consistent with the [Tolman–Ehrenfest equilibrium relation](https://journals.aps.org/pr/abstract/10.1103/PhysRev.36.1791): temperature times lapse is constant and blackbody energy density scales as temperature to the fourth power. The stress relation itself does not require a blackbody spectrum, but it does assume static isotropic radiation with the stated equation of state.

For an ordinary weak galactic field, dlnN/dlnr is approximately v_c^2/c^2. At 200 km/s the required density slope has magnitude only 1.78e-6. Across a factor of ten in radius, u changes by a factor about 0.9999959, whereas an r^-2 profile changes by 0.01. Gravity this weak therefore does not support the desired steep profile as a static isotropic radiation fluid under these assumptions.

This is not an escape-time calculation or a confined radiation solution. Boundary pressure, trapping, scattering and anisotropy require separate laws. Nor is it a prohibition on relativistic self-gravitating radiation structures: those differ from this weak-field galactic comparison, as noted in the earlier [candidate specification](../microphysics/action-specification.md). No new spatially varying time mechanism is adopted here; N denotes the ordinary static gravitational lapse used to test this comparison.

## Consequences for the unified model

The earlier rho_deposit=T q_capture/c^2 expression is conditional on a location-preserving support mechanism. It cannot be used as a dynamical solution by itself. A long internal lifetime and a long residence at a given galactic radius are different properties.

A completed candidate needs one receiver/field identity, a velocity-and-position injection law, supporting orbital motion or field stresses, finite boundaries and escape, and a time-dependent potential with source depletion. Its full stress tensor must predict both motion and lensing. For internal energy stored in ordinary matter, the complete receiver moves and gravitates; its original and added energy cannot be evolved as unrelated species without justification.

This pass establishes conditional support requirements and an ideal positive distribution, not a photon-funded, stable galaxy. It preserves the full 20-task and 32-observation scope. Run `python -X utf8 research_work/results/deposit-support/check_deposit_support.py`; saved evidence is in `deposit-support-results.json`. Mathematical checks are not astronomical validation.
