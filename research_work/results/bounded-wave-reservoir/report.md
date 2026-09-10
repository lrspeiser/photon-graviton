# A bounded wave-conversion region with ordinary endpoints

**In this explicitly localized candidate, the frequency reduction survives exit into an ordinary propagation region, while the acquired reservoir energy travels out with the wave.** The executed example gives a wavelength/event-duration factor of 1.161834 and transfers 13.9292% of the initial wave energy to the traveling reservoir.

This advances the prior matched-endpoint ray examples by including the local wave and reservoir equations. It does not derive a boundary from gravity wells, solve matter dynamics or supply a complete relativistic momentum/stress theory. The localization is an optional test assumption, not an adopted physical screening mechanism.

## Hypothesis and known equations

Retain the preceding directional Hamiltonian, but make its coupling spatially bounded. For the right-moving sector:

`h = v(s,x)*e + P*[1-v(s,x)*s_x]`, `e=(D+B)^2/4`,

`v=1/[1+(b+a*s)W(x)]`.

Choose `W=sin^4(pi*x/4)` for 0<x<4 and zero outside. The profile and coefficients are hypothetical. The trigonometric window, Hamiltonian method and characteristic equations are established mathematics; no novelty or first-principles origin is claimed.

The electromagnetic equations preserve a right-moving chiral amplitude F=(D+B)/sqrt(2), with e=F squared/2. Starting at s=0 everywhere gives the exact branch s=t, even with this explicit spatial profile. Known Hamiltonian variation then gives

`F_t + partial_x(v*F)=0`,

`P_t + partial_x(v*P) = -v_s*e`,

`partial_t(v*e+P) + partial_x[v*(v*e+P)] = 0`.

Here `v_s=-a*W/[1+(b+a*t)W]^2`. For a>=0, the source transfers energy into the reservoir. Outside the profile, v=1, v_s=0 and both components propagate without further exchange. The explicit x dependence breaks translation invariance; energy conservation does not imply that the fixed environment's momentum has been accounted for.

## Characteristic calculation and checks

Let xi label initial wave positions and J=partial x/partial xi. Along a characteristic:

`xdot=v`, `Jdot=v_x*J`, `Fdot=-v_x*F`,

`Pdot=-v_x*P-v_s*F^2/2`.

Two independent exact identities are `J*F=F_initial` and `J*(v*F^2/2+P)=e_initial`. The latter verifies the local energy carried by each initial wave element, not just the integrated global total.

The initial compact carrier wave occupies -2<=x<=-1 in an ordinary region. Its initial reservoir energy is zero under this classical ansatz. The detector is at x=8, and evolution continues to t=24. All values are dimensionless; the interval is not assigned a physical length or cosmic age. Numerical characteristic sets of 129 and 257 rays, followed by tighter integration at 257 rays, are retained.

For a=0.1 and b=0, the tight calculation gives:

| Quantity | Result |
|---|---:|
| Remaining electromagnetic energy fraction | 0.8607079764 |
| Acquired reservoir energy fraction | 0.1392920236 |
| Central endpoint frequency ratio | 0.8607079764 |
| Central arrival-interval stretch | 1.1618342428 |
| Maximum characteristic energy-identity error | 1.19e-9 |
| Maximum amplitude-identity error | 3.11e-15 |
| Arrival-derivative versus frequency/Jacobian error | 2.88e-11 |

At t=24, all wave/reservoir characteristics are beyond x=21.35, so none remains within the 0<x<4 conversion region. The field clock variable still exists there; this statement concerns the acquired P energy associated with the passing pulse, not the disappearance of every postulated field.

Tiny negative reservoir values at the 1.8e-13 level are retained as integration error. They are not treated as a physical negative-energy source. Refinement checks and the analytic solution support the positive-transfer result.

## Exact frequency and timing result

The ray equation can be written as a linear equation for travel time versus position:

`dt/dx = 1+b*W(x)+a*W(x)*t`.

Differentiate with respect to emission time. With ordinary endpoints, the exact arrival sensitivity is

`S = dt_arrive/dt_emit = exp[a*integral W(x) dx]`.

For this window the integral is 1.5, so S=exp(0.15)=1.1618342427. The same characteristic map transports the carrier phase, giving received frequency divided by emitted frequency equal to 1/S. Because propagation speed is one at both endpoints, coordinate wavelength there stretches by S as well. This exponential accumulation is known mathematics, not a new redshift formula.

The static control a=0, b=0.4 delays the wave but gives zero reservoir gain and no lasting endpoint frequency shift. The dynamic case loses energy during transit and does not regain it merely because W becomes zero on exit. Once outside, v_s=0, so the photon and reservoir energies separately stop changing.

These are wave/characteristic controls, not a supernova or galaxy fit. Endpoint atomic standards are assumed ordinary because the optical coupling is stipulated zero there; their actual matter coupling and boundary dynamics have not been derived.

## The remaining momentum and relativity problem

For the canonical fields used here, the spatial-translation generator density is `pi=D*B-P*s_x` in the right-moving case. On s=t, s_x=0, so the reservoir contribution to this canonical density vanishes even though reservoir energy has a nonzero flux.

At the ordinary exit, `J_energy=e+P`, while the displayed canonical momentum density is e. Their ratio is S=1.161834 in this example. In c=1 units, an ordinary complete symmetric relativistic stress tensor would identify energy flux with momentum density. Therefore this canonical subsystem cannot simply be declared an isolated Lorentz-invariant collection of photons and ordinary massless companions. Additional physical fields, an explicit preferred-frame medium or a different completion would have to explain the missing stress/momentum relation. A standard improvement term cannot be assumed to supply an unaccounted net momentum for a localized pulse.

The static spatial profile can exchange momentum without explicit time-dependent work, but its source and dynamics are not modeled. This is a genuine remaining requirement, not resolved by the small numerical energy error.

## Relation to the nonexpanding requirement

This test does not introduce a changing spatial scale factor or use expansion to obtain the endpoint frequency shift. It relies instead on a time-dependent, spatially localized optical interaction. It therefore lies outside the homogeneous universal-metric cancellation established in the previous turn.

That distinction is not a completed physical explanation. If the profile is tied to a gravitational environment, its matter forces, backreaction and clock/ruler behavior inside the region must be calculated. Earlier prescribed lapse profiles already exposed possible force-sign difficulties. The new traveling reservoir does not automatically cure those effects.

No rule determining W from observed ordinary matter, gravity, void structure or companion density has been established. The source of a, ideal clock preparation, arbitrary three-dimensional directions, quantum stability, capture into wells and lensing all remain open. This example is a conditional way to preserve endpoint redshift, not permission to tune separate profiles along each observed sightline.

## Reproduction

Run `run.py`. The archived results include static and dynamic controls, numerical refinement, endpoint phase/timing quantities, local characteristic energy identities, the exact window integral and input hashes. No observational dataset or holdout is read, and no prior empirical coefficient is changed.
