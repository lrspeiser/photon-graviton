# Photon energy exchange with a bound cloud

## Question and limits

Can the same interaction that stretches light also increase the population in the proposed bound store? This calculation joins light and the preceding massive complex-scalar binding proxy. It does not adopt a massive particle or a new force as a law of the fictional universe. It tests one density-dependent interaction, not every possible companion-production mechanism. It is an original conditional derivation, not an observational inference.

## Occupation is an exact constraint

Take the Schrödinger-Poisson Hamiltonian H_SP from the [self-binding derivation](../companion-self-binding/derivation.md), including constant seed rest energy N m c0^2. Add photon rays:

```
H = N m c0^2 + H_SP[psi] + sum_i c0 |P_i| / n(X_i; |psi|^2).
N = integral |psi|^2 d^3x.
```

A density-dependent n is invariant under the global phase change psi -> exp(i theta) psi. The kinetic and gravitational terms share that symmetry. The corresponding canonical field equations therefore conserve N: dN/dt = {N,H} = 0, assuming no occupation flux through the boundary. This is a symmetry proof; fixing N in a simulation cannot independently verify it.

Light can change the cloud's kinetic, binding and interaction energy without changing N. An initially empty classical field remains empty for a regular density coupling of this type. Thus this model cannot derive the seed population from photons. Quantum pair creation in a relativistic completion is a separate possibility, not calculated here. For a charged field, net charge and total particle-plus-antiparticle occupation are different quantities; conserving charge need not prohibit neutral pair production. Our nonrelativistic single-species model has no antiparticle sector.

## Gaussian collective dynamics

Use a spherical Gaussian density with fixed mass M=N m and radius R:

```
rho(r) = M exp(-r^2/R^2)/(pi^(3/2) R^3).
V(R) = A/R^2 - B/R,
A = 3 hbar^2 M/(4m^2), B = G M^2/sqrt(2pi).
```

A homologous radial flow v(r)=Rdot r/R has kinetic energy 3M Rdot^2/4. The radius therefore has inertia I=3M/2 and canonical momentum Pi=I Rdot. This one-radius approximation retains breathing motion but omits other modes, fragmentation, radiation escape, capture and changes of occupation. It is not full Schrödinger-Poisson evolution or a lifetime calculation.

Choose radius unit R0=2A/B, energy unit E*=B/R0 and time unit R0/c0. In the dimensionless equations below the potential is 1/(2R^2)-1/R and inertia kappa=I c0^2/E*. The seed rest energy in these units is M c0^2/E*=2 kappa/3. The protocol chooses kappa=10000; it does not fit a galaxy.

Let n=1+eta R^(-3) exp(-X^2/R^2). This is a positive density coupling, with the density normalization absorbed into eta at fixed N. It is not a derived weak-gravity/void coupling. Oppositely directed equal photon pairs cross through the cloud center. Reflection symmetry cancels net momentum transfer; one representative photon per pair is evolved with weight two. This avoids hiding a net recoil in a fixed-center support, although the radial Gaussian ansatz still excludes nonspherical deformations.

```
H/E* = 2 kappa/3 + Pi^2/(2 kappa) + 1/(2R^2) - 1/R + sum_i 2 P_i/n_i.
Rdot = Pi/kappa.
Pidot = 1/R^3 - 1/R^2 + sum_i 2 P_i n_R/n_i^2.
Xdot_i = 1/n_i.
Pdot_i = P_i n_X/n_i^2.
```

On the symmetric pair subspace the canonical one-form is Pi dR + sum_i 2 P_i dX_i. This factor of two in both photon energy and symplectic form gives Xdot=1/n, not 2/n.

The initial cloud is at R=1, Pi=0. Every photon is present upstream initially, with P_i chosen to give its specified initial energy; no external source or pump injects unrecorded energy. Density tails are included at every point. Negative X is the representative ray's initial side; its partner follows -X with opposite physical momentum.

## Work and energy destination

Along each ray E_i=P_i/n_i:

```
dE_i/dt = -P_i n_R Rdot/n_i^2.
dH_cloud/dt = sum_i 2 P_i n_R Rdot/n_i^2.
```

Photon loss and cloud gain cancel exactly. The separate work accumulator integrates the second expression and is compared with independently evaluated cloud energy. The constant seed rest energy is reported but subtracted from the numerical conservation residual to avoid masking small transfers.

The unforced initial cloud already minimizes this one-radius potential. Any net positive energy left in it after light exits is excitation above that minimum, not additional constituent rest energy or a more tightly bound ground state. It can still add to gravitational mass through total energy in a compatible gravity theory. No joint lensing or metric calculation is performed here. A stationary configuration with larger N would require a production/capture channel and its own binding-energy account.

The cloud can give energy back to later photons as it breathes. A positive final transfer in the tested finite trains does not establish a one-way interaction, permanent storage, a universal redshift rate or homogeneous event stretching. The original propagation relation E_dot/E=-partial_t n/n applies here with partial_t n=n_R Rdot; independent weak-signal timing and matter-clock tests for this new background remain to be done.

## What a production mechanism must add

The next derivation must specify whether it produces neutral companion excitations, charged particle/antiparticle pairs, or a real spin-2 field with different conserved quantities. Its interaction must predict energy and momentum transfer, allowable frequencies and reverse processes. If it changes N, it cannot simply reuse this fixed-N Hamiltonian and insert an independent positive capture rate. This requirement prevents counting an unexplained preexisting cloud as photon-created gravity.

Reuse the existing [vacuum interaction constraints](../interaction-rate/derivation.md) and [matter-assisted energy/momentum calculation](../matter-assisted/derivation.md) when constructing that channel. Those studies do not yet supply the production/capture interaction for this massive bound cloud; they prevent treating energy-only bookkeeping as a complete reaction.
