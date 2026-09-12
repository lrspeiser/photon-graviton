# Radiation-fed spatial time field: the tested realization gives blueshift

We coupled the selected radiation/receiving-energy response to spatial transport and the operational clock rule. The field now evolves from an initial receiving-energy seed and passing radiation, rather than from a prescribed sinusoidal time profile. It produces deposits when capture is enabled. However, the tested probes are measurably blueshifted after endpoint clocks are included. This particular realization does not supply the desired redshift explanation.

## Explicit model

In reference units c0=kappa=1, define n_t=T, v=1/n and h=T/n. The postulated transport/transfer equations are

    G_t + (vG)_x = -hG
    H_t + (vH)_x = -hH
    T_t + (vT)_x = h(G+H) - Gamma*T
    deposit_t = Gamma*T.

G and H are reference-energy densities of EM and GW radiation; T is a receiving field, not yet ordinary gravitons. All three travel with coordinate speed v in this specialization. The response/capture assignments are project postulates, while the continuity equations, ray calculation and clock transformations use known mathematics.

The domain is [0,10], with constant radiation inflow G=0.2,H=0.02, no incoming receiving field, and an initial T seed centered at x=2 with width 0.4. Capture uses an illustrative Gaussian rate centered at x=8; it is not a dynamically generated gravity well. The seed is not derived from photons. Radiation amplifies the pre-existing receiving field; a zero-seed control stays inactive. Background radiation drives the fields and negligible probes measure transport without adding counted energy.

## Refined result

The 640-cell second-order reconstruction gives:

| Case | Probe path | Measured z |
|---|---|---:|
| No seed | Both paths | 0 |
| Traveling receiving field | 0 to 10 | -0.034954 |
| Traveling receiving field | 2 to 9 | -0.025020 |
| Receiving field with capture | 0 to 10 | -0.014567 |
| Receiving field with capture | 2 to 9 | -0.007909 |

Negative z means a blueshift. These numbers use demonstration units and are not an astronomical fit. Identical EM/GW transport is imposed; simultaneous probes have no relative messenger delay.

The receiver's clock matters. Measured stretch is S=J*n_e/n_o, where J is the coordinate arrival-map derivative and d tau=dt/n is the clock rule. The traveling field leaves a changed n behind it. For these probes the endpoint factor overcomes the coordinate stretching, producing event compression and higher locally measured frequency. Energy loss defined in reference coordinates must not be mislabeled as measured photon redshift.

For the probe beginning at x=0, receiving energy is initially concentrated ahead of it and moves in the same direction. The refined coordinate stretch is nearly one, while the receiver has already acquired a different clock rate. This is a particularly clear example of an absolute journey delay not guaranteeing event stretching or redshift.

Capture reduces the receiving field but n_t=T only increases n: removing T stops further clock change without reversing the accumulated state. Our earlier prescribed node-to-node example assumed unchanged endpoint clocks. The present source evolution does not automatically produce such endpoints. That difference is physical within the postulates, not a reason to suppress the endpoint factor.

## Numerical failures were retained and resolved to a stated tolerance

First-order upwind transport smeared the traveling field. Its interior-start result changed sign on refinement and still failed the initial 0.003 absolute-z convergence gate at 320 versus 640 cells. All first-order outputs remain archived. Minmod MUSCL reconstruction was then used with the same physical equations and 80,160,320,640 cells. Its 320-to-640 absolute-z differences are at most 0.000659 and pass the original gate. This supports the negative signs in these examples, not arbitrary precision or correctness of the physical theory.

All refined calculations preserve reference energy, including boundary inflow/outflow, within 2.47e-13 absolute units. Reservoirs are nonnegative to numerical roundoff. Neighboring-ray proper-clock intervals agree with the phase/arrival relation within 1.36e-5 relative across the runs. Shared arrival equality follows from the same solver/law and is not an independent gravity-wave simulation.

The reference-energy balance omits a full physical cost for the accumulated temporal memory and its gravitational/momentum effects. It must not be presented as a complete conserved proper-energy or stress-energy theory. No local action is required as the next prerequisite, but an operational energy and force account remains necessary.

## Next physical question

The next revision should specify whether capture also changes the accumulated temporal state, and what energy/force accompanies that response. It should also specify the seed's origin/distribution and whether receiving energy truly has the same propagation law as both messengers. Those choices must come from explicit postulates and then predict all endpoints, rather than manually resetting receivers or choosing only paths that redshift.

This calculation advances the spatial coupling of roadmap goals 2–4 and reveals a concrete failure of the current initialization/response/transport combination. It does not disprove all radiation-fed time models, but this version is not a successful redshift candidate. No fitted alpha, observational data or holdout changed. Reproduce with `python research_work/results/spatial-temporal-feedback/run.py`; append `--first-order` for the archived lower-order method. All nine goals remain active and incomplete.
