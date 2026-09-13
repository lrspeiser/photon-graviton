# Orbital support: an existence result with an angular-momentum requirement

13 September 2026. No new stellar-speed fit. The original and shared-settled profiles, exact-third inventory and ordinary-matter baselines remain fixed.

## What changed

The prior pressure-law equilibria distorted the successful empirical profile. A tangential orbital population can support that profile in a stipulated spherical potential without changing its density. This is a possible mechanical realization, not an explanation of why capture creates it. It replaces the simple thermal-condensate interpretation with an explicitly different bound-state hypothesis.

The settled profile cannot be supported by a nonnegative isotropic energy-only distribution under the tested spherical assumptions: a necessary positivity condition fails in both matter baselines. The original unmodified profile passes that necessary condition, which alone does not establish a physical isotropic distribution.

Circular-orbit shells with random orbital planes and phases provide a nonnegative, singular phase-space construction for either profile. In 0.1-30 kpc the individual circular orbits have negative binding energy and positive radial epicyclic frequency squared. This is not a test of collective stability, a smooth distribution function or the full flattened Galaxy.

## Isotropic positivity condition

For an isotropic energy-only distribution f(E), with relative potential Psi and binding energy E,

\[
\rho(\Psi)=4\pi\sqrt{2}\int_0^\Psi f(E)\sqrt{\Psi-E}\,dE,
\qquad
\frac{d\rho}{d\Psi}=2\pi\sqrt{2}\int_0^\Psi\frac{f(E)}{\sqrt{\Psi-E}}\,dE.
\]

These are known distribution-function identities; see the [Eddington derivation](https://galaxiesbook.org/chapters/I-04.-Equilibria-of-Collisionless-Stellar-Systems_6-Spherical-distribution-functions.html). If f is nonnegative, the derivative cannot be negative. Since Psi falls outward in this potential, a resolved outward increase in density rejects the class. We use this necessary test before attempting an ill-conditioned second-derivative numerical inversion; failing it is sufficient for rejection, passing it is not proof of positivity.

The settled baseline-I profile has an outward-rising interval near 0.297-0.938 kpc, with maximum logarithmic density slope about 0.089. Baseline II has a shallower interval near 9.94-10.25 kpc, with refined maximum slope about 0.024. Both persist when the source-phase grid and evaluation grid are refined. Fine subdivisions near baseline II's transition depend on interpolation, so the result is reported as a resolved broad interval, not separate physical shells. These features belong to our constructed profiles, not observed density measurements.

## A tangential alternative

In a spherical potential, choose the known circular speed and specific angular momentum

\[
v_c^2(r)=r g(r),\qquad j_c(r)=r v_c(r).
\]

Populate each positive-mass shell with circular orbits distributed uniformly over orbital planes and phases, with opposite senses of rotation. There is no required net galaxy-wide spin. Locally, using tangential speed vt, the singular velocity distribution can be written

\[
f_r(v_r,v_t)=\frac{\rho(r)}{2\pi v_c(r)}\delta(v_r)\delta(v_t-v_c(r)),
\]

normalized with velocity measure 2 pi vt dvt dvr. All orbit weights are nonnegative. Uniform occupation keeps each shell's density stationary in the assumed potential. This is known orbital mechanics applied to the proposed companion profile, not new physics or a uniquely predicted distribution.

The tangential Jeans balance is rho times vt squared divided by r = rho g. Individual radial perturbations are checked through

\[
\kappa_r^2=\frac{dg}{dr}+\frac{3g}{r}.
\]

The tested profiles have positive kappa-r squared throughout the reporting region. Because the density was fixed, keeping the earlier gravitational fit is a construction identity, not a new predictive success. A continuum of exactly circular orbits is an extreme tangential limit; broadening the orbits or applying a nonspherical perturbation is substantive further work.

## Settling must transfer angular momentum

If the initially captured population were also on circular orbits, moving the selected material from r to sr requires comparing

\[
j_i=\sqrt{r^3g_i(r)},\qquad j_f=\sqrt{(sr)^3g_f(sr)}.
\]

Using the initial/final potentials and the previously selected mobile fractions gives:

| Diagnostic for moved mass initially at 0.1-30 kpc | Baseline I | Baseline II |
|---|---:|---:|
| Mean fractional reduction in specific angular momentum | 17.13% | 17.89% |
| Range across contributing shells | 16.13-33.72% | 16.41-34.74% |
| Fraction of moved mass requiring reduction | 100% | 100% |
| Mean single-orbit specific energy decrease | 10,844 (km/s)^2 | 10,119 (km/s)^2 |

The energy numbers compare individual orbit energies in different potentials. They must not be added as a global cooling budget: companion self-energy would be double-counted. The earlier properly counted settling ledger remains the appropriate global budget.

The angular-momentum reductions are conditional on circular initial states. They quantify a transport requirement rather than demonstrate capture. Isotropic incoming flux can have zero mean angular momentum while still carrying a distribution of individual angular momenta. Zero net spin neither guarantees nor forbids the required circular population.

## What this suggests next

An interaction that transfers angular momentum to an extended companion population while another channel carries away energy is a possible next branch. It must conserve the combined angular momentum and energy. A separate orbital population could supply support while a smaller dissipative population settles. Neither mechanism should be accepted just because it can be described verbally.

The next calculation should compare required torque with what radiation and companion exchanges can carry, followed by finite-width orbital populations and their response to the real flattened potential. This is more specific than adding an arbitrary pressure correction. The photon-conversion, propagation/timing, absolute supply, lensing and nuclear/diffuse radiation tasks remain unsolved and active.

## Verification

[Protocol](orbital-support-protocol.md), [code](orbital-support.py) and `orbital-support-results.json` preserve all classifications, intervals, orbital quantities and transport diagnostics. Doubling the radial grid, doubling the phase-source grid and increasing angular quadrature preserve the positivity classifications; sampled orbital speeds change by less than 0.0006 km/s. Point-mass and uniform-density-core checks reproduce epicyclic-to-orbital frequency-squared ratios 1 and 4 within 1e-5. Tangential Jeans balance passes to numerical precision.

These checks support a conditional existence result, not physical formation or collective stability. Ordinary matter is spherically averaged for the orbital construction; randomly oriented circular planes are generally not preserved in the true flattened/barred potential. The broad research goal is not complete.
