# A working transfer law with separated travel and capture

The user's proposal is adopted here as a **conditional working transfer rule**, not dismissed because its microscopic origin is unfinished. We can state its conservation properties and predictions now. It is separate from the prescribed time-field candidate until an interaction derives the same rate and energy receiver.

## Spatial form of the energy account

**Provenance: proposed phenomenological transfers written with established continuity equations; no originality claim.** In fixed nonexpanding spatial coordinates and one declared common energy/time convention, define photon, traveling-companion and deposited energy densities u_gamma, u_c and u_d:

partial_t u_gamma + div(F_gamma) = j - Q,

partial_t u_c + div(F_c) = Q - C,

partial_t u_d = C.

Here Q is local energy transferred per volume per time, C is capture into deposits, and j is injected photon power density. The illustrative closures are F_gamma = c n_hat u_gamma, F_c = v_c u_c, Q = c alpha u_gamma and C = Gamma_cap u_c. Rates are nonnegative, deposits stay in place, and secondary loss, reverse transfer and scattering are absent by assumption. These prescriptions are not derived graviton physics.

**Provenance: derived by adding the established continuity equations.** For a fixed volume V:

d/dt integral_V (u_gamma + u_c + u_d) dV = integral_V j dV - integral_boundary (F_gamma + F_c) dot n dA.

Including cumulative exported energy and the fuel supplying j closes this scalar ledger. Leaving the volume is not energy destruction. Capture may occur far from conversion and need not occur along the surviving photon's path. The closures above use a straight companion path only for a calculable example; steering and capture momentum require their own dynamics.

In an evolving q field, a local energy density and a common-reference energy density cannot silently be interchanged. The appropriate field-work terms and receiver must be derived before identifying this ledger with the arrival-test Hamiltonian. Nor can arbitrary gravitational-field energy be assumed to be a positive local fluid u_d merely by naming it 'warping space'.

## Concrete path benchmark

**Provenance: standard transport and exponential survival solutions under the declared closures.** Photons cross a conversion interval with integrated fractional loss A = integral alpha dx. Companions made there travel independently at speed v_c through a later capture interval with K = integral Gamma_cap dx/v_c. Assume every companion crosses that interval, no intervening loss, and enough time for transit. For unit emitted energy:

photon fraction at detector = exp(-A),

deposited fraction = [1-exp(-A)] [1-exp(-K)],

escaped companion fraction = [1-exp(-A)] exp(-K).

The three sum to one. The missing photon energy does not automatically all become locally deposited energy; some can escape without losing its energy. At A = ln(2) and K = 1, photons retain 0.5, deposits receive approximately 0.31606 and escaping companions retain approximately 0.18394. These are illustrative dimensionless choices, not fitted astronomical rates. The photon remains at half energy after exiting the conversion interval.

The script checks 18 cases: three integrated conversion depths, three capture depths and two conversion profiles with the same integrated depth. Numerical integration of energy generated at each conversion position agrees with the analytic fractions. Conversion-free and capture-free controls are retained. This is a spatial extension of the older one-region ledger, not a new galaxy supply estimate.

## What energy transfer alone does not determine

For this particular stationary, fixed-speed closure, each photon emitted at time t_e reaches a detector a fixed distance L away at t_arr = t_e + L/c. **Provenance: established constant-speed travel relation.** Successive arrival intervals equal emission intervals, even though each photon's energy decreases. In the half-energy example the wavelength stretch is two, while the whole-event stretch is one.

That is a concrete demonstration of underdetermination: scalar energy conservation and an energy-versus-frequency rule do not fix the propagation of an entire signal. It is not a theorem that every companion mechanism fails to stretch events. The separate changing-q calculation does stretch neighboring arrival times, but its receiver remains incomplete. A combined theory must obtain both results from one consistent interaction rather than choose energy loss from one closure and timing from another without deriving their compatibility.

The number-retaining, phase-unspecified transfer closure is a packet-energy model. It does not yet define a coherent wave equation, spectral linewidth, image fidelity or microscopic momentum exchange. Those gaps are why this is a workable hypothesis rather than a finished first-principles solution. They do not prevent writing or testing the transfer law.

Run `python research_work/results/transport/spatial-transfer.py`. All examples and assumptions are saved in spatial-transfer-results.json. No galaxy parameter, catalog distance, redshift outcome or environmental input was fitted or changed. The next mechanism task is a common energy-and-phase propagation completion, retaining the locally normal clocks and lasting-arrival shift already required by the user.
