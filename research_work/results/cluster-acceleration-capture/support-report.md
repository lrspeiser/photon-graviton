# Can the deposited profile remain in place?

## Result

The calculated central depletion cannot be maintained by a steady, bound, isotropic collisionless particle population in the fixed Plummer well. A circular-orbit ensemble can support the same spherical density as a test-particle construction, but the capture/settling mechanism and self-gravity remain unspecified. This narrows the possible reservoir models without rejecting all companion or field-storage ideas.

## Why isotropic random motion fails

Assume massive nonrelativistic collisionless particles, steady spherical equilibrium, no continuing phase-space source, bound velocities, and density rho(r) proportional to the calculated q(r). In the finite-depth Plummer potential, relative potential Psi=-Phi is largest at the center and decreases outward. The known isotropic distribution-function relation is

    rho(Psi)=4 pi sqrt(2) integral_0^Psi f(E) sqrt(Psi-E) dE.

For f(E)>=0 this is nondecreasing in Psi. Our rho instead starts at zero and increases as r^2 away from the center, while Psi decreases. Therefore no nonnegative isotropic f(E) produces this profile under these assumptions. This is a phase-space consistency problem, not a numerical fitting failure. Related central density/anisotropy restrictions are derived by [An and Evans](https://arxiv.org/abs/astro-ph/0511686); the isotropic argument above is direct and does not require applying their theorem outside its stated conditions.

The same problem is visible in the established isotropic Jeans moment solution with vanishing pressure at infinity:

    sigma_r^2(r)=[integral_r^infinity rho(s) g(s) ds]/rho(r).

The numerator tends to a positive finite constant while rho(r)~r^2, so sigma_r^2 diverges as r^-2. Escape speed remains finite. At r=0.001a the required radial variance exceeds escape-speed-squared by factors approximately 28,595, 28,894, 42,531 and 4.84e10 for A=0.1,1,10,100. These are impossible bound-particle variances, not predictions that actual particles attain those speeds. Density normalization cancels, so reducing deposited mass does not remove this particular inconsistency.

Continuous injection, transient evolution, stresses or interactions can invalidate the equilibrium assumptions, but their transport/evolution equations then need calculation. A Jeans moment solution alone cannot establish a valid particle distribution.

## A limiting tangential construction exists

Place particles on circular orbits at each radius with speed v_c^2=r g(r), distributing orbital planes and phases uniformly, with equal senses of circulation. Weight each radial shell by 4 pi r^2 rho(r) dr. Each orbit remains at its assigned radius and the ensemble is spherical without net rotation. This is a positive but singular phase-space measure, not proof that a smooth or naturally formed distribution exists.

Its radial stress is zero and each tangential stress component is p_t=rho v_c^2/2. The known radial Jeans equation is satisfied:

    dp_r/dr+2(p_r-p_t)/r=-rho g,

because -2p_t/r=-rho g. The finite-core density and stress remain regular near the center. In the fixed Plummer potential the squared radial epicyclic frequency is

    omega_r^2=GM(r^2+4a^2)/(r^2+a^2)^(5/2)>0,

so individual circular orbits are stable to small radial perturbations in this imposed field. This is not a collective stability result for a self-gravitating deposit population.

Numerical test-particle integrations at 0.3a, a, 3a and 10a for ten orbital periods retain radii within 7.54e-11 relative and orbital energy within 4.16e-11 relative. The code checks bound circular speeds and positive epicyclic frequency. Units set G=M=a=1; no observed galaxy mass or time is fitted.

## What the capture model must now provide

Incoming companions cannot simply be declared stationary pressureless particles while retaining this profile. A viable version must either:

- Populate sufficient tangential orbital motion, with the associated momentum, energy and angular-momentum accounting;
- Supply stresses or a binding mechanism for a field/medium reservoir; or
- Evolve the deposited distribution dynamically rather than equating its density with the cumulative local injection.

Circular-orbit support has kinetic energy and gravitational binding; these must be included consistently when interpreting the deposited energy as mass. No process producing the required orbits has been derived. If deposits gravitate significantly, both the capture profile and orbital field change and must be solved together. The prior projected shear signature consequently remains conditional on retaining the injection shape.

Formula status: these are established distribution-function, Jeans and orbital-dynamics relations applied to the candidate deposition law. No new first-principles support physics is claimed. No observational holdout was opened and all six objectives remain open.

Run `python research_work/results/cluster-acceleration-capture/support.py` after run.py.
