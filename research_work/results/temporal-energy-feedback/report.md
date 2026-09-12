# First calculation of the selected shared-time feedback branches

The user selected common light/gravitational-wave temporal transport and added radiation energy feeding a temporal field that wells can capture. The full [branch specification](../../../research_plan/shared-time-and-energy-feedback.md) supplies equations, units, assumptions and unresolved requirements. This calculation is a nondimensional mechanism experiment, not an observation fit.

## What is now explicit

Both wave types use the same proposed evolving wave Hamiltonian omega=c0*|k|/n. Standard ray/phase mathematics then links frequency stretching and infinitesimal event stretching through A=n(arrival)/n(emission). The identical ray rule gives zero differential delay for simultaneous co-located signals along the same path, even when their common absolute delay is substantial. This equality is imposed by the shared coupling, not independently discovered by two distinct wave solvers. Actual source emission lags remain separate.

The additional branch represents 'energy into time' by a physical receiving-field energy u_T and response n_dot=kappa*u_T. Radiation loses energy at rate h=n_dot/n. Both EM and GW energy enter the receiving sector; capture and export transfer that energy onward. Deposits are not counted twice. Time itself is not assigned energy units. The response and capture laws are project postulates, while their energy sum and phase identities are known mathematics.

## Results of four open-cell cases

All values use arbitrary demonstration units. The path length, field coupling, export residence time and initial n equal one. A small seed is 0.01 energy units; baseline EM/GW driver energies are 1 and 0.1. A test packet is emitted at t=0.1. These choices are not calibrated astrophysical constants.

| Case | Frequency/event derivative stretch A | EM minus GW arrival time | Deposited energy at t=12 | Field energy exported by t=12 |
|---|---:|---:|---:|---:|
| No initial receiving field | 1.000000 | 0 | 0 | 0 |
| Seeded open void cell | 1.010627 | 0 | 0 | 0.128058 |
| Tenfold background radiation | 9.574273 | 0 | 0 | 10.010577 |
| Seeded capture cell | 1.006018 | 0 | 0.010966 | 0.010966 |

The field energy exported from a cell is retained in the total ledger as energy elsewhere, not lost energy. The model has no permanent void deposit. Capturing some receiving-field energy reduces the subsequent temporal-change rate and places that energy in the deposit reservoir. No galaxy potential or deposited-energy-to-lensing conversion was solved.

The dramatic response to tenfold background radiation is a warning about this particular feedback closure, not an astronomical prediction. Here alpha_eff=kappa*u_T/c0 varies as the field is supplied and exported. It is not automatically the approximately constant coefficient in our existing redshift fit. Even at equal path length, different background histories can yield markedly different redshifts. Such a dependence must be calibrated or excluded with data, rather than hidden in object-specific initial conditions.

An exactly empty field remains empty in these equations. Radiation alone does not initiate conversion, because the transfer rate contains u_T. This makes the seed/initiator requirement explicit. A different initiation law is possible but is not derived in this calculation.

## Event intervals, memory and clocks

For small neighboring emissions the derivative stretch agrees with the phase stretch. A finite source interval need not have a single constant multiplier if n changes substantially during the event: in the bright case the 0.001-time-unit interval stretches by 9.580208, versus the local derivative 9.574273. This is a potential signal-shape distortion and cannot be dismissed as a numerical error.

n integrates past field energy in this proposed response. Export/capture makes n_dot smaller without resetting n to its original value. Thus the example has temporal memory; this is an assumption about a physical state, not an explanation of why that state costs no additional energy. A complete action must account for any kinetic, interaction or memory-state energy and momentum. The present positive-energy ledger is necessary but not sufficient for such a completion.

The largest unresolved issue is matter-clock coupling. With fixed rods and clock standards the wave speed is c0/n. A universal coordinate lapse can instead cancel the measured homogeneous redshift. Neither result can be repaired simply by renaming time. The user's desired ordinary local c, unchanged local clocks and nonexpanding spatial geometry need a separate physical completion, still absent here.

## Checks and next step

All four integrations completed with nonnegative energy sectors and total energy, including exports, conserved to at most 1.25e-14 absolute in these units. Independent optical-path integration supplies arrival times; centered finite differences verify the event-stretch derivative. The earlier forward-difference check failed in the bright case and was replaced by centered differencing, with that numerical change recorded in the branch specification. None of these checks verifies an action, spin-2 coupling or observed merger waveform.

Reproduce with `python research_work/results/temporal-energy-feedback/run.py`. Keep the prescribed shared-n branch and the radiation-fed branch separately: the first reproduces the old exponential calculator by construction, while the second predicts a changing rate that has not been fitted. Next derive the matter/field coupling and feedback stability, then jointly test brightness dependence, spectra and timing. Existing independent-distance/selection work remains required. No full-model freeze or fresh holdout was performed; the four-stage goal remains active.
