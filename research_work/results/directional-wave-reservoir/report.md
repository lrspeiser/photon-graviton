# A one-dimensional wave model with traveling energy reservoirs

**A local Hamiltonian candidate can shift electromagnetic mode frequencies while transferring energy into reservoirs that travel in the corresponding wave directions.** The executed one-dimensional example conserves total energy and has an explicit local energy current. It does not establish a three-dimensional, relativistic, quantum or observationally viable theory. In particular, ordinary local light-speed and clock behavior remain unresolved.

## What is postulated, and what mathematics is known?

Temporal frequency conversion is known physics; it is not a new effect claimed here. For example, [Bohn, Luk and Hendry (2021)](https://arxiv.org/abs/2107.05970) discuss frequency shifts from temporal refractive-index changes and the importance of boundaries. Such material results do not establish a cosmological mechanism or this proposed reservoir.

The following Hamiltonian coupling is an exploratory project postulate using standard canonical-field mathematics. Its physical realization and originality in the literature have not been established. It is not derived from graviton microphysics or void observations.

Use one-dimensional transverse electromagnetic variables D and B in normalized reference units c0=1. Define the two nonnegative directional energy terms

`e_sigma = (D + sigma*B)^2/4`, with `sigma=+1,-1`.

Their sum is the usual normalized electromagnetic quadratic energy. Introduce separate canonical clock/reservoir pairs `(s_sigma,P_sigma)` for the two directions, and set

`v_sigma = 1/(1+a*s_sigma)`.

The proposed Hamiltonian density is

`h = sum_sigma [v_sigma*e_sigma + P_sigma*(1-sigma*v_sigma*partial_x s_sigma)]`.

The extra spatial-gradient term is essential: it transports acquired reservoir energy rather than leaving it at the transfer location. The directional reservoirs are additional degrees of freedom. They are not identified with ordinary atomic clocks, measured particles or proven gravitational waves.

## Equations derived from the postulate

Using a vector potential with B=partial_x A and canonical momentum -D, the known Hamilton equations give

`D_t = -partial_x(h_B)`, `B_t = -partial_x(h_D)`,

`s_sigma,t = 1-sigma*v_sigma*s_sigma,x`,

`P_sigma,t = -sigma*v_sigma*P_sigma,x - v_sigma'(s_sigma)*e_sigma`.

Here `v_sigma'=-a*v_sigma^2`. Therefore the reservoir source is nonnegative for a>=0 and nonnegative electromagnetic energy. These equations are conditional deductions from the coupling, not independently measured laws.

Known local Hamiltonian energy accounting supplies

`partial_t h + partial_x J = 0`,

`J = sum_sigma sigma*v_sigma*[v_sigma*e_sigma + P_sigma*s_sigma,t]`.

An independent numerical substitution with spatially varying clock fields and nonzero reservoir fields gives a maximum local identity residual of 2.08e-14 in the chosen dimensionless test. That checks more than the spatially homogeneous special case; it is not a proof of nonlinear stability.

The Hamiltonian is not asserted globally bounded below for arbitrary P and clock gradients. On the declared branch P>=0 and `1-sigma*v*s_x>0`, its displayed terms are nonnegative. For smooth clock characteristics starting from s0, their coordinate Jacobian is `1+sigma*s0'*(v(t)-v0)`; if `1-sigma*v0*s0'>0` initially, this Jacobian remains positive as v decreases, and the clock-rate factor remains positive. The reservoir momentum increases along its own characteristic. These conditional forward-time properties do not prove a globally stable interacting quantum theory or justify restricting all physical states to this branch.

## Executed counterpropagating-wave example

Initialize s_plus=s_minus=0 and P_plus=P_minus=0. Then both clock fields obey s=t and v(t)=1/(1+a*t). With directional amplitudes `F_sigma=(D+sigma*B)/sqrt(2)`, the equations reduce to

`F_sigma,t + sigma*v*F_sigma,x = 0`,

`P_sigma,t + sigma*v*P_sigma,x = a*v^2*F_sigma^2/2`.

For each initial wave profile, let X(t)=ln(1+a*t)/a. The exact conditional solution is

`F_sigma(x,t) = F_sigma,0[x-sigma*X(t)]`,

`P_sigma(x,t) = [1-v(t)]*F_sigma,0[x-sigma*X(t)]^2/2`.

Thus each gained-energy profile travels with its own wave. Using only one right-moving reservoir would not handle a left-moving wave; the two-direction construction explicitly addresses that one-dimensional issue. It does not solve arbitrary three-dimensional propagation, polarization or interactions between rays.

The numerical test uses two Gaussian-envelope waves with different carrier wavenumbers, a periodic length-20 domain and final t=4. The nonzero-rate case has a=0.1. These are dimensionless diagnostic choices with no astronomical calibration. The waves stay away from the numerical boundary over this interval.

At 2,048 cells:

| Quantity | Result |
|---|---:|
| Remaining propagation energy / initial energy | 0.714285705 |
| Acquired reservoir energy / initial energy | 0.285714280 |
| Maximum fractional total-energy error | 1.53e-8 |
| Maximum wave-amplitude error against exact translation | 4.03e-7 |
| Maximum reservoir-profile error | 2.86e-7 |
| Carrier-frequency error against time-dependent prediction | 5.71e-9 |

The exact remaining fraction is 1/1.4 and the exact acquired fraction is 1-1/1.4. Tiny negative reservoir samples, approximately 1e-16, are retained as spectral roundoff rather than interpreted as negative physical energy. Runs at 512, 1,024 and 2,048 cells and an a=0 control are archived. The mode-frequency estimates use a short interval near the endpoint; its midpoint prediction is used, rather than treating an interval estimate as an exact endpoint derivative.

Unlike the earlier packet model with a constant E0 offset, this is a different classical field ansatz initialized at P=0. No offset energy is added in this test. That does not prove that a physical synchronized clock field can be prepared or operated at zero cost; the existence and preparation of this ideal clock sector remain assumptions.

## Frequency and event-arrival timing from the same characteristics

In the homogeneous branch, an existing spatial Fourier mode has conserved wavenumber k and frequency magnitude `omega=|k|/(1+a*t)`. Its coordinate spatial wavelength does not stretch: its frequency falls because its reference-coordinate propagation speed falls. This distinction matters for the user's desired locally constant measured c.

For a boundary emitter a fixed coordinate distance d away, characteristic travel gives

`1+a*t_arrive = exp(a*d)*(1+a*t_emit)`.

Consequently `dt_arrive/dt_emit=exp(a*d)`. A boundary phase proportional to fixed reference-clock emission time gives `omega_arrive/omega_emit=exp(-a*d)`. Whole source-time features and the carrier phase therefore have reciprocal transformations in this branch; the timing factor was not independently appended to an energy-only loss rule.

Nine ray-integral controls over three distances and three emission times reproduce these relations to 6.22e-15. They are analytic/numerical controls, not supernova fits. The exponential accumulation law is known; no novelty is claimed for it. Boundary matching to actual matter emitters and detectors remains unspecified.

## What this accomplishes—and what it does not

This candidate joins a one-dimensional electromagnetic wave phase, local energy transfer and directional reservoir transport within one postulated Hamiltonian. It advances beyond assigning a photon energy to a point packet while leaving the wave equation unspecified.

It still does not satisfy the full desired universe. The homogeneous speed is 1/(1+a*t) in fixed reference units. Local atomic clocks, rulers, emission standards and detectors have not been coupled. Rescaling all matter clocks to match the optical factor cannot be assumed harmless: the prior clock analysis showed how that can cancel the observed shift. No solution to that issue is claimed here.

The construction also lacks a three-dimensional covariant action, a complete physical stress-energy tensor, quantum stability and carrier identification, a creation mechanism for the directional clock states, and a law for capture into the required gravitational source. A preferred propagation decomposition is supplied in one dimension; isotropy and ordinary relativity are not thereby established. Captured energy is not automatically curvature or lensing mass.

The prescribed parameter a has not been derived from low-gravity environments, stellar data or microscopic interactions. No existing redshift, galaxy or lens parameter is refitted. The photon supply budget remains deferred. This is a candidate worth subjecting to the clock and stress-energy requirements, not evidence that the full theory matches observations.

## Reproduction

Run `run.py`. `results.json` records all resolutions, frequency and energy errors, the inhomogeneous local identity control, boundary timing controls and the script hash. No catalog or holdout data is read.
