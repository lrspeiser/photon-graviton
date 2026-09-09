# Locally normal clocks, lasting signal stretch and the accompanying force

This optional candidate answers a specific part of the user's premise: clocks can behave normally in local experiments while a changing region produces a lasting received redshift between matched endpoint standards. It also predicts a force whose sign depends on whether the void clock is faster or slower relative to clocks in wells. The test supplies no self-consistent source for the prescribed time field, no galaxy fit and no completed companion energy reservoir.

## A common matter and light clock

Established lapse/metric and Hamiltonian structures, applied here as a candidate matter completion; not new formulas. The standard geometric framework is described by [Gourgoulhon, 3+1 Formalism](https://arxiv.org/abs/gr-qc/0703035). We use its clock/particle kinematics, without adopting Einstein's field equations as mandatory for the fictional universe.

ds^2 = -c^2 q(t,x)^2 dt^2 + dx^2 + dy^2 + dz^2, q>0.

For observers fixed in these coordinates, d tau = q dt. A pointlike stationary clock with ordinary local Hamiltonian H_local is assigned H_reference = q H_local. All internal energy terms receive the same factor, rather than changing only the electromagnetic coupling while leaving other matter parameters fixed. In a spatially uniform patch the Schrodinger equation i hbar d psi/dt = q H_local psi becomes i hbar d psi/d tau = H_local psi. A constant H_local also commutes with itself at different times, so a changing common multiplier does not by itself excite this ideal clock.

Consequently local dimensionless atomic frequency ratios and ideal local ticking stay normal in this candidate. Extending this statement to finite laboratories requires q to vary negligibly across them; tidal effects, spatial gradients and moving observers need their usual treatment. A traveler does not see their own watch slow. Spatially separated signal comparisons can still differ. This construction is a stipulated universal matter coupling, not a derivation from first principles or confirmation that real atoms obey a new field.

Established relativistic particle dispersion with the same lapse factor:

H_m = q sqrt(c^2 p^2 + m^2 c^4), H_gamma = q c |p|.

Hamilton's equations give photon dx/dt=cq in the positive-x direction and dx/d tau=c. Fixed physical spatial lengths are used; there is no spatial expansion in this ansatz. The earlier electromagnetic-medium atomic-spectrum failure held different matter parameters fixed, so it does not automatically apply to this different universal prescription. Neither prescription is thereby established as the correct physical model.

## What this predicts for gravity

Derived from the established particle Hamiltonian within this candidate; no originality claim:

a_local,at_rest = -c^2 grad ln(q).

Here a is the instantaneous acceleration of a freely released massive particle relative to the fixed local observers, evaluated at rest. It is not the proper acceleration experienced by an ideal freely falling particle. The formula follows by differentiating v_local = (1/q) partial H_m/partial p with respect to local time and using dp/d tau = -(1/q) partial H_m/partial x at p=0. The symbolic check returns exactly the expression above.

For an additional q that approaches one in wells, a lower q in voids pulls matter toward voids. A higher q in voids pulls matter toward wells. These are additional-force signs in this restricted coupling; an ordinary background gravitational field, independently derived modified-gravity dynamics or deposited sources can contribute other forces. They cannot simply be assumed to cancel an inconvenient term. The sign result is not a theorem forbidding every slower-void mechanism.

## Stretch between matching clocks

Established Hamiltonian energy identity, applied to the positive-x photon ray:

d ln(E_reference)/dt = partial_t ln(q),

ln(S_time) = - integral_ray partial_t ln(q) dt,

alpha_time = -partial_t q/(c q^2).

The last expression uses dx=cq dt and is the same established fractional-change rate as the earlier n=1/q representation, not a newly discovered redshift formula. At endpoints with q=1, the reference energy ratio equals the ratio measured against the endpoint clock standards. A decreasing q at the places/times traversed causes redshift. Merely having q smaller than one does not determine the sign of the accumulated effect.

The redshift/motion/endpoint accounting remains necessary for observations. In a general case with unequal endpoint clock factors, measured S includes q_observer/q_emitter multiplying the reference stretch. In a homogeneous universal-clock change this can cancel the reference shift, as previously noted. The localized matched-endpoint case tested here is distinct.

## Prescribed examples and numerical evidence

Example profiles proposed for this diagnostic only; originality unverified:

q=1+sigma epsilon(t) I(x), with I=sin(pi x)^2 for 0<x<1 and I=0 outside.

All quantities below are dimensionless. Rays go from x=-1 to x=2. Source/detector q is always one. The same ray law propagates carrier energy and nearby arrival-time spacing; an independent finite difference of arrival times checks the timing result. The solver is repeated with tighter tolerances and a fourfold smaller maximum step.

| Void history | epsilon(t), sigma | Received z, launch at t=0 | Received z, launch at t=2 |
|---|---|---:|---:|
| Static faster | 0.2, +1 | 0 | 0 |
| Static slower | 0.2, -1 | 0 | 0 |
| Faster, relaxing toward normal | 0.2 exp(-0.2t), +1 | 0.012188 | 0.008693 |
| Slower, slowing further | 0.2 - 0.1 exp(-0.2t), -1 | 0.009057 | 0.006316 |
| Slower, returning toward normal | 0.2 exp(-0.2t), -1 | -0.018494 | -0.011490 |

In every run the instantaneous event stretch agrees with the wavelength stretch; finite-difference arrival timing passes the declared 2e-6 absolute tolerance. This checks nearby wavefronts and neighboring launch times. The different results at launch times zero and two show why an entire finite event need not have one constant stretch factor. It is not a successful supernova light-curve prediction, a line-width calculation or a three-dimensional image-fidelity test.

## A bound on the inferred gearing

Derived here for the separable monotone examples using established calculus; originality unverified:

If partial_t q<=0, 0<=I<=1 and q_void(t)=1+sigma epsilon(t)>0, then

1 <= S_time <= q_void(t_emit)/q_void(t_arrive).

Proof: -partial_t ln(q) = -sigma epsilon_dot I/(1+sigma epsilon I). For a positive denominator, I/(1+sigma epsilon I) increases with I. Its maximum is at I=1; integrating that maximum gives the stated ratio. This bound spans the whole journey in a common time-dependent amplitude and a fixed spatial profile; it is not automatically valid for independent regions with different histories.

Thus observed stretch can impose a minimum required change in the hypothetical void clock contrast under this model. It does not measure the traveler's subjective clock rate or a unique constant slowing throughout the trip. When the journey is nearly all fully affected void, the bound can approach equality. Its time evolution still requires a physical explanation.

## Energy and unfinished physical work

Established autonomous Hamiltonian conservation, conditional on a proposed field coupling:

H_total = H_field[phi,Pi] + integral q(phi,x) h_local d^3x.

If this Hamiltonian has no explicit external time dependence and the field/matter equations all follow from it with controlled boundaries, total energy is conserved. The field source includes -partial_phi q times the entire coupled matter/radiation energy density, not just starlight. This is a framework requirement, not a chosen or solved field Hamiltonian. It would be inconsistent to require universal clocks but omit the associated matter source and force without a derived suppression mechanism.

In the numerical examples q(t,x) is prescribed externally. The photon energy loss is computable, but its receiver and the work needed to maintain q are not solved. Naming that receiver a companion field does not establish conservation of the full system. A self-consistent modified-gravity/field action must generate the environmental profile and evolution, include all matter sources, and reproduce the desired net attraction and deposited reservoir. The earlier pure-flat-slice GR source constraint remains a conditional comparison, not a prohibition on the user's modified gravity.

Research decision: preserve local normality and constant local c, keep both sign branches explicit, and use this force/clock connection when constructing the physical model. Do not fit a time-contrast amplitude solely to redshift and ignore the resulting matter force. Keep the 26 pending data candidates unscored while the independent-data and physical-law work continue.
