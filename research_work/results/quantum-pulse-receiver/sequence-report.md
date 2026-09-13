# Reusing the coherent receiver without resets

The stipulated energy-conserving spectral converter can act on two distinct photon modes through the same receiver, including intervening free evolution, without destroying its predicted relative Fourier-profile stretching. In its continuum Gaussian limit, finite receiver coherence gives a common timing offset rather than independent jitter for the two outputs. This closes a specific repeated-use gap in the earlier pulse calculation. It does not close the spatial locality, causal arrival, preparation, or escaping-companion gaps.

## Retained postulate and conditional derivation

Retain the existing energy/momentum-preserving map on designated input bands, in c=1 units:

    |E,r> -> |E/S, r+(1-1/S)E>.

The receiver momentum increases by the same amount as its energy; its stipulated spectrum has E_R-P_R=U>0. Reverse swaps complete the finite map to a unitary. S is still chosen, not derived from distance or a void interaction. This construction and the derivation below use known quantum translation, partial-trace and Fourier mathematics. No fundamental novelty is claimed. The role of temporal coherence under energy-conserving operations is discussed by [Lostaglio et al.](https://arxiv.org/abs/1410.4572).

Let delta=1-1/S. For N labeled photon modes, successive forward-branch translations add, giving the joint output

    integral [product_i a_i(E_i)] |E_1/S,...,E_N/S>
             tensor |R shifted by delta sum_i E_i> dE_1...dE_N.

Receiver gain is delta times the sum of incoming photon energies, with no reset or missing accumulated energy. This equation holds while every occupied state stays on the specified forward branch; it is not a proof of indefinite operation for a finite receiver.

For an initially coherent Gaussian receiver with energy probability variance sigma_R^2 and phase reference t_R, tracing the receiver multiplies each joint photon density-matrix element by

    G(E,E') = exp[-delta^2 (sum_i E_i - sum_i E'_i)^2/(8 sigma_R^2)]
              exp[-i delta (sum_i E_i - sum_i E'_i) t_R/hbar].

This is the usual Gaussian overlap identity. It depends on the difference in TOTAL transferred energy, not a product of independent single-photon overlaps. The receivers have not been reinitialized between photons.

In output-frequency variables this is equivalent to a mixture of ideal dilated profiles translated by one shared Gaussian time offset xi, with

    Var(xi) = [(S-1) hbar/(2 sigma_R)]^2
    t_profile,i = S t_input,i - (S-1)t_R + xi.

For initially independent profiles, the receiver contribution to Cov(t_i,t_j) is Var(xi) for i!=j. It cancels from Var(t_j-t_i). Intrinsic packet widths still contribute; the claim is not zero uncertainty in photon times. For separate independent receivers with the same parameters, the receiver contribution to a time difference would instead be 2 Var(xi). Partially correlated receivers give 2 Var(xi)(1-r) for correlation coefficient r. These are conditional joint-profile predictions, not yet astronomical predictions: the model has not specified which physical receiver different photons encounter.

## Free evolution does not mean resetting a temporal reference

For each ideal swap U_i, [U_i,H_total]=0. Including all photons and the receiver in free evolution F(Delta)=exp(-i H_total Delta/hbar) gives

    U_2 F(Delta) U_1 = F(Delta) U_2 U_1.

We verified this identity on an explicit joint state, rather than freezing the receiver phase between uses. It says a common reference-time description is consistent with the sequential calculation. It does not say operation times are physical arrival times, supply the interaction's switching mechanism, or make arbitrary operations spacetime-local. Modes are labeled and already represented in the joint Hilbert space; their actual emission/local encounter dynamics are absent.

## Executed calculation

sequence.py uses nine input energies between 9.6 and 10.4, S=2, receiver mean energy 55 and U=50, all in arbitrary units. Two photon modes each have mean energy 10 before conversion and 5 after. Receiver gain is 10 in total. Its state is translated after the first interaction and carried into the second; no trace-and-reset replacement is made. Three initial coherence widths and waiting intervals 0.3, 2 and 11 were evaluated.

- Joint-state partial traces agree with the summed-transfer formula to 3.5e-17.
- Moving the included free evolution outside the sequence changes amplitudes by at most 1.8e-14.
- Total mean energy residual is at most 1.5e-14; momentum is conserved by the same basis-state increments and fixed U.
- The first photon's marginal state is unchanged by the second use, to 8.4e-17, while joint correlations remain possible.

For input phase centers 0 and 4 and receiver reference 1, the continuum formula predicts output centers -1 and 7, separated by 8. These are analytic Fourier-profile coordinates, not measured or simulated causal arrival times. In particular, the negative first center must not be advertised as a physically emitted photon arriving before emission. The finite calculation verifies the spectral state and free evolution; it does not independently simulate a detector clock or fit these profiles to data.

## What this changes for the six objectives

The receiver cannot yet be rejected merely because consecutive uses require resetting its coherent state: this ideal branch does not. It also supplies a precise distinction between shared and independent receiver noise. The already exposed timing observations require approximately equal wavelength and event stretching; the ideal profile relation is compatible with that requirement, but no new observational fit or validation was performed here.

The missing step is now a local encounter model that realizes this operation and sends acquired energy onward as companions. Without it, we cannot identify the Fourier centers with supernova arrival stages, predict S from a galaxy's distance, or supply the cluster capture source. Energy remains in the receiver in this calculation. A prescribed frequency scaling plus a special receiver spectrum is not a first-principles derivation. All six objectives remain open, and final observational holdouts remain unopened.
