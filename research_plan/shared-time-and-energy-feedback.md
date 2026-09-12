# Selected research branches: shared temporal transport and radiation-fed temporal field

Current closure qualification: the [finite-response transport audit](../research_work/results/finite-response-transport/report.md) tests an explicit moving field-energy budget. Its homogeneous checks pass only with a stated availability rule, and that rule still fails a spatial empty-companion boundary. Neither the initial cell demonstration below nor this extension is a completed energy/clock mechanism. Shared messenger transport remains a postulate; observable redshift must retain endpoint clock factors.

User decision, 2026-09-12: proceed with option 1 (light and gravitational waves share temporal transport), and add a branch where photon energy enters whatever physically changes time and can subsequently be captured by deep gravity wells. Both remain hypotheses in a static-spatial-geometry universe. This records research directions, not validated laws. Earlier failed candidates remain archived.

## A. Shared evolving temporal transport

Postulate a dimensionless positive field n(x,t) in fixed spatial coordinates and the same geometric-optics Hamiltonian for both messenger types:

    omega_gamma = c0 |k| / n(x,t)
    omega_GW    = c0 |k| / n(x,t).

The Hamiltonian form is known wave/medium mathematics; the equal electromagnetic and tensor-wave coupling in this nonexpanding setting is a proposed physical assumption. No new tensor action or gravitational polarization law is derived here.

Known ray identities give group speed c0/n and, along a ray,

    d ln(omega)/dt = -partial_t ln(n)
    d ln(omega)/ds = -partial_t n / c0.

Thus a static spatially varying n changes propagation times but does not accumulate frequency loss in a stationary medium. Temporal evolution is essential in this particular branch. Shared coupling does not mean photons literally become detected gravitational waves.

For spatially homogeneous n, define t_o by integral[t_e,t_o] c0/n(t) dt = D. Differentiation and phase transport give

    A = dt_o/dt_e = n(t_o)/n(t_e) = omega_e/omega_o = 1+z.

These are conditional mathematical consequences, not independently invented redshift formulas. Both waveform frequencies and closely spaced event intervals transform together. The same ray has the same propagation delay for both messengers. A source emission lag remains an arrival lag transformed by their shared arrival map; simultaneous source events have zero relative delay in the ideal specialization. This does not predict the observed 1.74 seconds without a source-emission model.

For n(t)=n_e+alpha*c0*(t-t_e) over a finite applicable epoch, the ray equation recovers 1+z=exp(alpha D) and the prior enormous absolute delay. It is shared by both messengers rather than being an EM-only relative delay. Alpha remains the old empirical calibration; the n evolution is prescribed, not explained. The finite-epoch qualification avoids asserting positive n for an arbitrarily extended linear history.

Ordinary endpoint clock rates C_e,C_o change measured stretch to A_measured=A*C_o/C_e. Equal stable clock standards must be physically specified; merely renaming a coordinate time does not suffice. With fixed physical rods/clocks, c0/n is a varying measured wave speed. Declaring all local clocks to scale as 1/n instead can cancel the homogeneous measured redshift. Prior clock/ruler-completion results therefore still apply. Do not claim simultaneously unchanged local c, unchanged rods/clocks, nonexpansion and this observable redshift without a new matter-coupling completion.

## B. Radiation-fed temporal field and capture feedback

Interpret 'energy converted into time' as energy converted into a physical receiving field whose state changes temporal transport. Time itself is not assigned energy units. Let u_gamma, u_GW, u_T and u_d denote EM radiation, gravitational radiation, temporal-field/companion energy, and deposited energy densities, all measured in one explicitly declared reference frame. u_T is not yet identified with ordinary gravitons.

An exploratory local closure is

    partial_t n = kappa u_T
    h = (partial_t n)/n = kappa u_T/n
    partial_t u_gamma + div(F_gamma) = j_gamma - h u_gamma
    partial_t u_GW    + div(F_GW)    = j_GW    - h u_GW
    partial_t u_T     + div(F_T)     = h(u_gamma+u_GW) - Gamma_cap(B) u_T
    partial_t u_d                    = Gamma_cap(B) u_T.

kappa has units inverse time per energy density. n is dimensionless; h and Gamma have inverse-time units. The local response/capture closures are new project postulates, with no established novelty claim. The continuity equations and their sum are standard energy bookkeeping. A possible capture rate is Gamma_cap=Gamma0*B/(B+Bstar), with B>=0 an escape/binding-depth proxy relative to a stated surrounding boundary, not the magnitude of acceleration: acceleration can vanish at the center of a deep well. This shape and Bstar are uncalibrated choices, not a fitted deposition map.

The sum obeys

    partial_t(u_gamma+u_GW+u_T+u_d) + div(F_gamma+F_GW+F_T) = j_gamma+j_GW.

Capture moves energy into the well; it does not create a second copy of energy already held by the temporal field. To test the user's traveling-companion condition, retain the temporal-field flux F_T. In a representative open cell, model it as u_T/t_res and record that exported energy separately. This is transport out of the cell, not energy decay. Do not silently accumulate all temporal energy permanently in voids.

This closure yields alpha_eff=partial_t n/c0=kappa*u_T/c0. It offers an actual reason for a changing redshift coefficient, rather than setting alpha to the observed value on every ray. It also predicts radiation-history dependence and feedback: stronger radiation can replenish u_T, export removes it locally, and capture transfers it into deposits. In an initially empty field u_T=0 there is no conversion at all: a seed or an additional physically specified initiation mechanism is required. Seeding is not a derivation of first cause.

## Scope of the immediate calculation

Use nondimensional open-cell examples to test bookkeeping, feedback and common signal stretching; no physical parameter fit or observational pass is claimed. Set c0=kappa=t_res=1, initial n=1 and radiation energies u_gamma=1, u_GW=0.1. Compare empty seed, seed 0.01, tenfold radiation, and capture Gamma=1. Integrate exported field energy alongside all other reservoirs. The cell approximates residence/energy exchange, not a solved spatial field or galaxy halo. Trace test packets through the computed homogeneous n(t), with a fixed test path D=1 and emissions at t_e=0.1 and 0.101. Record equal messenger arrivals, frequency ratio, and event-stretch derivative. The radiation reservoirs are background drivers; no extra probe energy is counted as an independent physical source. Numerical note: initial forward differencing at t_e=0 failed the derivative tolerance in the rapidly changing bright case; the test was changed to centered differencing at t_e=0.1, allowing both neighboring emissions within the integrated domain. The physical closure and all four cases were unchanged.

Required checks: total energy including exports, nonnegative reservoirs, independent integrated optical path and finite-difference event stretching, and no seed/no conversion. These cannot verify gravitational field equations or local-clock consistency.

## Requirements before either branch can be adopted

1. Supply a local causal field/matter action or equivalently complete evolution equations, including field energy and momentum. An energy ledger does not prove that n can evolve without an additional kinetic/interaction-energy cost; such terms may alter the closure above.
2. Derive observable clock/rod behavior and local measured wave speed. Equal messenger coupling alone leaves this unsolved.
3. Derive tensor-wave amplitude/polarization/energy transport. Equal frequency stretching is not enough to reproduce a merger waveform or luminosity-distance inference.
4. Specify spatial transport, capture and the response of the well to deposited energy. Neither the capture shape nor u_d establishes the amount/location of gravitational acceleration or lensing.
5. Test radiation-history/environment dependence against independently estimated distances, without fitting a separate temporal field to each galaxy. Preserve survey selection requirements from stage 2.
6. Use the same parameters for spectra, event timing, brightness and relative delays; freeze only after these predictions and physical assumptions are complete, then evaluate genuinely withheld data.

The source history/initial seed, coupling kappa, capture scale and clock coupling are open physical parameters, not newly locked empirical constants. This is a concrete formulation of the user's additional hypothesis, not evidence that energy creates time in the real universe.
