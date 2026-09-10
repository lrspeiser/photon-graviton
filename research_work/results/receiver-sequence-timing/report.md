# Sequential receiver phase can retain a stretched Fourier event sequence

The stipulated quantum mode map can stretch the centers of a sequence of independent photon profiles while the receiver evolves freely and is not reset between encounters. The sequence does not require a coherent superposition between different emission events. This extends the previous separate-trial result, but remains a conditional spectral construction, not a spacetime-local propagation or supernova timing solution.

In plain language, the receiving system's evolving quantum phase can retain a common timing reference. It is incorrect simply to insert a positive elapsed time into the old receiver-reference parameter: that parameter has the opposite sign to ordinary free Schrödinger evolution in the adopted convention. Both local and common-epoch calculations were executed and agree.

## Assumptions and known mathematics

Retain the hypothetical band map E_out = E_in/S for S = 2 or 1.5 and the receiving spectrum P_R = E_R - 100 in units c = hbar = 1. Neither the frequency factor nor this spectrum is derived from a local field action. The initial receiver is a coherent Gaussian with mean energy 100 and width 4, truncated exactly to 60 <= E_R <= 140. Its initial energy and coherence are counted resources.

The receiver evolves by exp(-i E_R t), the known Schrödinger free-evolution rule for its stipulated Hamiltonian. After each photon it releases a fixed forward companion energy L = 2.5 or 1.75, as in the earlier release alternative. The receiver retains the difference delta(E)-L. Four encounters remain safely within the positive-energy support, so no lower-boundary release failure occurs. Fixed L is NOT exact eventwise transfer of each photon's lost energy and is not adopted as the user's final transport rule.

Tracing earlier outgoing photons leaves a mixture of translated receiving states. Its probabilities evolve by convolution with the fresh photon energy distribution; no state reset is imposed. Free evolution of an energy translation differs from a translated freely evolved state only by a global phase. This known identity permits the receiving overlap matrix to be evaluated while carrying the mixture forward. Support-extreme and peak-probability translation checks are included, along with direct partial traces of every receiving mixture component in two-frequency controls.

## Correct timing-coordinate relation

The prior Fourier convention used input spectral phase exp(i E t_b) and receiver phase exp(i E_R t_R), giving

    t_profile,out = S t_b - (S-1) t_R.

This is known Fourier/overlap algebra applied to the hypothetical mode map. In a local encounter description, an incoming packet is centered at zero. A receiver freely evolved from its preparation epoch to encounter time b has reference parameter t_R = -b, not +b. Consequently

    local output profile center = (S-1)b,
    common-coordinate output profile center = b + (S-1)b = S b.

Applying the common-epoch channel and then freely evolving the output to the encounter gives the same reduced state. This coordinate identity is explicitly verified. Holding the receiver phase fixed in the local encounter frame instead produces output centers b; it is a different control requiring phase intervention, not the free-evolution history.

## Executed sequence and incoherent-envelope results

| Chosen S | Encounter times | Output Fourier-profile centers | Fixed-local-phase control centers |
|---|---|---|---|
| 2 | 0, 30, 60, 90 | 0, 60, 120, 180 | 0, 30, 60, 90 |
| 1.5 | 0, 30, 60, 90 | 0, 45, 90, 135 | 0, 30, 60, 90 |

All three profile-center interval ratios equal the chosen S in each sequence. Eight direct two-frequency partial-trace controls have maximum discrepancy 3.89e-16. Energy and forward momentum bookkeeping checks pass. Receiver variance grows from 16 to 16.0196 or 16.0044444 after four uses, despite mean energy remaining approximately 100 for the centered bands. Its state is not restored.

An incoherent mixture of the four events with weights 0.1, 0.2, 0.4 and 0.3 was also propagated. Input mean profile time is 57. Output means are 114 and 85.5. Input standard deviations are 28.526393 and 28.740216; output values are 57.052923 and 43.110369. The conditional Gaussian identity

    Var(t_out) = S^2 Var(t_in) + (S-1)^2/(4 sigma_R^2)

is verified to the stated numerical tolerance, with hbar = 1. This follows from known mixture/Fourier mathematics, not a new universal time law. Removing the common evolving phase instead leaves the envelope mean at 57; broadening individual packets alone does not produce the same sequence scaling.

## What remains unresolved

The source encounter schedule is supplied. The channel acts instantaneously on energy modes and has no physical input/output positions, local interaction Hamiltonian, causal boundary conditions or detector response. Fourier profiles are periodic finite-band functions; the selected windows avoid recurrence artifacts but do not prove arbitrary-time signal causality. An output mode profile is not by itself a demonstrated travel-time distribution.

The receiver's spectrum, preparation, motion and stress need physical explanation. The fixed-release rule still exchanges energy with that receiver rather than passing each photon's exact loss directly to its own companion. There is no derivation of the distance/environment-dependent S, protected atomic clock behavior, graviton identity, capture, deposited gravity or full source-energy budget. No real supernova light curve or observational holdout was fitted.

The next substantial step is to realize or rule out this phase behavior within a finite-history local propagation model. This result narrows a mathematical concern about receiver resetting; it does not complete priority #2 or justify treating the unified theory as validated.

Run `python research_work/results/receiver-sequence-timing/run.py`. `results.json` retains upstream and script hashes, sequence controls, mixture moments and resource accounting. The full research goal remains active.
