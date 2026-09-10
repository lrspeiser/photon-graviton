# Does the small-transfer candidate stretch coherent pulses?

2026-09-10. Follow-up to the fixed-energy Markov model. This is an explicit mathematical quantum-channel completion and a phase-coordinate test, not a local photon/graviton interaction, real detector model or astrophysical fit. [Calculation](coherence.py), [results](coherence-results.json).

**Result:** the mean photon energy still falls by q=exp(-alpha D), but translated pulse centers retain their original separation. This remains true after retaining off-diagonal frequency coherence. Making epsilon smaller cannot remove this obstruction while keeping the same stationary channel.

## Conditional completion and energy bookkeeping

**Known loss-channel/binomial-isometry mathematics, reinterpreted as the proposed frequency-energy ladder; no novelty claim:**

\[
|m\rangle_\gamma|0\rangle_c\longmapsto
\sum_{k=0}^{m}\sqrt{\binom{m}{k}}
q^{(m-k)/2}(1-q)^{k/2}
|m-k\rangle_\gamma|k\rangle_c.
\]

The previous postulate treated m as photon energy in units epsilon. Here |k>_c records transferred energy k*epsilon. It is not an assertion that m counts the incident photons of a monochromatic beam, and it is not an identification of the receiving modes as physical gravitons. Distinct input states have orthogonal total-energy sectors, so the map is an isometry; tracing the companion sector gives a trace-preserving quantum channel and reproduces the earlier binomial populations.

Every branch conserves energy m*epsilon=(m-k)*epsilon+k*epsilon. If all momentum is assigned forward with P=E/c, it also balances that one momentum component. This is kinematic bookkeeping only: it establishes neither an allowed angular-momentum transition nor local field stress, a nonzero decay amplitude or an astrophysical interaction rate. A particular stationary bath realization and Markov approximation would still require physical justification.

The zero-energy state is retained in the channel. For the pulse-center diagnostic it is excluded as a non-detection, and the nonzero-frequency subspace is normalized. Its probability is negligible in the tested packets but not removed from the energy ledger.

## Why pulse separation remains unchanged

Let U(delta) attach phase exp(i*m*delta) to input frequency level m; this shifts the Fourier pulse coordinate by delta under the stated sign convention. **Conditional deduction from the isometry:**

\[
\mathcal E_q[U(\delta)\rho U^\dagger(\delta)]
=U(\delta)\mathcal E_q[\rho]U^\dagger(\delta).
\]

In each term surviving the companion trace, both photon frequency indices lose the same k. Their difference, and hence the coefficient of delta in their relative phase, is unchanged. This is time-translation covariance: a shifted input yields the same shifted output, not an output shifted by delta/q. Conditioning on the nonzero-energy subspace does not change that covariance.

Thus this stationary receiving state can reduce energy and partially reduce coherence without universally stretching event intervals. It may change individual pulse shapes; a shape change is not a multiplication of the separation of arbitrary copies. The earlier population-only energy calculation could not establish this phase result by itself.

## Numerical check

A coherent packet spanning 97 energy levels has mean level 60 and a Gaussian input amplitude with probability width three levels. For q=0.99,0.8,1/1.885801, test source-coordinate pulse centers 0.2,0.5,1.2 in units hbar/epsilon. The output mean energy ratios are the specified q, while output centers remain 0.2,0.5,1.2. For example, the last case would need centers approximately 0.3772,0.9429,2.2630 for universal dilation; it does not produce them.

Trace preservation, positivity, mean-energy scaling, branch energy/forward-momentum identities and full matrix covariance pass. Covariance errors are below 1e-15. These values verify the finite channel and Fourier-coordinate calculation, not observed light curves or a causal source-to-detector flight law. The reported purity refers to the complete reduced state and documents coherence loss; it is not a claimed observational signal.

## What this eliminates and what remains

The simple stationary small-transfer completion can reproduce mean redshift and has a linewidth prediction, but it is not a complete candidate for the combined redshift/event-timing goal. It cannot be repaired by adjusting only the transfer size, alpha or a static spatial coupling profile: changing total attenuation changes q, not the covariance relation.

A receiving state with a changing phase reference, a time-dependent interaction or different clock/source dynamics would be a different completion. Such modifications must predict energy exchange, clocks, pulse timing and messenger delays together. The earlier coherent-receiver work explored a nonstationary reference and found resource/propagation gaps; those results remain relevant and unresolved. We do not attach a free event-stretch multiplier to this channel.

No observations were fitted, no new source was opened and no holdout changed. The redshift calculator is intact; this narrows its possible causes. Independent physical environment inference and a complete surviving model before withheld evaluation remain required. The full goal remains active.
