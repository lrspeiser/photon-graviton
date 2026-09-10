# Energy-conserving conversion needs a specified receiving state

**A balanced photon/companion energy ledger does not determine whether the outgoing light retains its temporal information.** We now have an explicit, energy-conserving quantum example that separates these issues. A stationary receiving state loses the tested interference pattern; a coherent, nonstationary receiving state can preserve it to high accuracy while gaining exactly the energy lost by the photons.

This is a synthetic mechanism audit, not an astronomical fit, a demonstrated graviton interaction, or a solution to supernova timing. It advances the main physical question identified in the [theory-priority review](../theory-priority-review/report.md): what receiving/driver dynamics could supply the desired optical transformation?

## What was already known, and what this adds

The earlier [frequency-transfer calculation](../frequency-transfer/conversion-kernels.md) showed that whole-photon removal, stochastic partial transfer and continuous energy drift can share the same energy equation while producing different spectra. The [electromagnetic specification](../../../research_plan/electromagnetic-transfer-specification.md) already distinguishes a fixed path delay from event stretching. We have not rerun those as if they were new findings.

The additional calculation here explicitly tracks the joint photon–receiver quantum state. It tests whether the receiver records enough information about the energy transfer to destroy coherence between photon frequencies. The source/receiver resonance calculations remain separate: they did not specify this coherent state transformation.

## A known symmetry constraint

**Established quantum mechanics, not a unique formula or a new theorem:** let a joint unitary interaction U conserve the sum of free photon and reservoir energies, and let the initially uncorrelated reservoir state be stationary:

\[
[U,H_\gamma+H_R]=0,\qquad [\rho_R,H_R]=0,
\qquad
\mathcal E(\rho)=\operatorname{Tr}_R U(\rho\otimes\rho_R)U^\dagger.
\]

Then

\[
\mathcal E(e^{-iH_\gamma t/\hbar}\rho e^{iH_\gamma t/\hbar})
=e^{-iH_\gamma t/\hbar}\mathcal E(\rho)e^{iH_\gamma t/\hbar}.
\]

This means the channel respects time translations. The standard result follows by inserting the reservoir's unchanged free evolution, commuting the joint free evolution through U, and tracing the reservoir. Thermal states are one example; the stationarity condition is the relevant assumption here, not a requirement that every receiver be thermal. Energy-conserving operations and their time-translation constraints are discussed in [Lostaglio et al., Physical Review X (2015)](https://arxiv.org/html/1410.4572).

For a coherent superposition with energies E1 and E2, its relative phase evolves at (E2-E1)/hbar. A perfect coherent fractional remapping E→E/S changes this phase-evolution rate to (E2-E1)/(S hbar). For S≠1, this is not the same time-translation behavior. Consequently, the specified stationary channel cannot implement that ideal coherent remapping for arbitrary incident superpositions. A photon energy eigenstate can still lose energy. The constraint concerns the full transformation, including coherence, not just its energy populations.

This conditional restriction assumes standard quantum dynamics, a fixed photon Hamiltonian at the compared endpoints, an initially independent stationary receiver and an energy-preserving joint interaction. A modified endpoint clock relation, correlations, an evolving field or a more general propagation geometry changes the assumptions and must be treated explicitly. It is not a prohibition against every fictional time/interaction law.

## An exact finite model we executed

**New project diagnostic built from known unitary/energy-eigenspace mathematics; no claim of fundamental novelty:** use photon energy levels 1, 2 and 4 in arbitrary units, and a receiver with a finite ladder from 0 to 80 in steps of 0.02. Within each allowed fixed-total-energy block, the code permutes basis states. This is exactly unitary and preserves each total-energy eigenvalue, a stronger check than conservation of mean energy alone. The inverse permutation exists; the operation is not claimed to be an irreversible fundamental one-way law.

On the populated input subspace the operation gives

\[
|2\rangle_\gamma|e\rangle_R\longrightarrow|1\rangle_\gamma|e+1\rangle_R,
\qquad
|4\rangle_\gamma|e\rangle_R\longrightarrow|2\rangle_\gamma|e+2\rangle_R.
\]

The initially equal photon superposition has mean energy 3. Its final mean energy is 1.5, and the receiver gains 1.5. If the receiver starts with a sharp energy, its two possible final energies distinguish the branches. Tracing over that receiver leaves the photons with the desired lower energy populations but no coherence between the two output frequencies.

We then use a receiver whose **amplitudes**, not merely probabilities, form a coherent Gaussian across energy levels. The shifted receiver states overlap. For an untruncated Gaussian with energy standard deviation sigma, their overlap is

\[
\mathcal V=\exp[-(\Delta E_{\rm loss})^2/(8\sigma^2)].
\]

**Known Gaussian overlap identity:** V is the visibility of this two-frequency beat diagnostic, and the two branches' lost energies differ by one unit here. It is not a supernova-duration formula or a measured cosmological parameter.

| Receiver state | Photon energy after / before | Receiver energy gained | Beat visibility |
|---|---:|---:|---:|
| Sharp initial energy | 1/2 | 1.5 units | 0 |
| Coherent Gaussian, sigma=0.25 | 1/2 | 1.5 units | 0.1353 |
| Coherent Gaussian, sigma=0.5 | 1/2 | 1.5 units | 0.6065 |
| Coherent Gaussian, sigma=1 | 1/2 | 1.5 units | 0.8825 |
| Coherent Gaussian, sigma=2 | 1/2 | 1.5 units | 0.9692 |
| Coherent Gaussian, sigma=5 | 1/2 | 1.5 units | 0.9950 |

Replacing each coherent Gaussian with an incoherent mixture having the **same energy probabilities** gives zero beat visibility in all six tested widths. A large energy spread or a hot reservoir is therefore not an interchangeable substitute for the coherent state in this construction. The highest-visibility case slightly differs from the untruncated formula because the finite ladder/support is explicit.

Numerical checks verify the basis permutation is one-to-one and preserves total energy; mean-energy residuals are at floating-point precision. The stationary-receiver channel obeys time-translation covariance. Holding a coherent receiver state fixed while shifting only the photon input breaks that reduced-channel symmetry; evolving **both** systems restores the closed system's free-evolution symmetry. Energy conservation has not been abandoned. The special initial receiver state supplies a temporal reference.

## What this does and does not tell us about our universe hypothesis

The useful positive result is that **energy loss with largely preserved optical coherence is not ruled out by energy conservation alone**. We constructed an example. Its extra requirement is a specified nonstationary receiving state and interaction; this resource cannot be omitted from the physics account.

The beat period between the two output levels is twice the input beat period because the energy difference was halved. That is not yet a prediction that a supernova's brightness history stretches by the same factor. Supernova radiation need not retain optical phase coherence across days. Predicting its arrival-rate history requires the full spatial and temporal channel for an ensemble of emitted wave packets, including the receiver's evolution between interactions. A converter that supplies a two-frequency beat transformation has not automatically supplied that channel.

Nor does a coherent state by itself explain why it exists throughout voids, has the required broadband interaction rate, preserves local clock/c measurements, or produces captured gravity. The receiver ladder has no spatial momentum, propagation, capture, stress-energy or gravitational model. The following repeated-use check accounts for its state change in the stipulated finite interaction; preparation, dynamical field backreaction and any later radiation/storage remain unmodeled.

The operation was stipulated for a few energy levels and chosen interaction, not derived from a local Lorentz/gauge-compatible field action. It does not predict alpha, depend on distance, establish achromatic behavior over the electromagnetic spectrum, or determine the gravitational coupling. The selected finite interval is not a universal monotonic down-converter on every possible photon/receiver state.

## Repeated use without resetting the receiver

We extended the receiving ladder to 300 energy units so that the entire populated support remains away from its artificial boundaries through 100 uses. This is a declared numerical domain change, not an inferred physical capacity. The same interaction is applied twice to the full joint state, and the two outgoing photons are traced/compared jointly. No reservoir reset occurs between those applications. No intervening free-flight dynamics is modeled.

In the interior regime, after tracing out n photons the receiver is an exact binomial mixture of its initial state shifted upward by n+k energy units, for k=0,...,n. Each photon deposits 1 or 2 units with equal probability in this test. **Known probability/translation algebra, applied to the stipulated interaction:**

\[
\langle E_R\rangle_n-\langle E_R\rangle_0=1.5n,
\qquad
\operatorname{Var}(E_R)_n-\operatorname{Var}(E_R)_0=0.25n.
\]

The direct two-step unitary agrees with the translated-state construction. The mixture calculation verifies the energy ledger at n=1,2,10,100 for all three tested widths. At n=100 the receiver gains 150 units, equal to the total photon energy loss. Its energy variance increases by 25 units squared. Those results account for accumulation; they do not establish gravitational storage or an absolute source budget.

The first and second photon have the same individual beat visibility, but their **joint state is not the product of their separate states**. For receiver widths 0.5, 2 and 5, the two-photon trace distances from that product are respectively 0.2162, 0.0294 and 0.00495. These are mathematical state-comparison diagnostics, not detection significances. Shared receiver correlations remain even where individual visibility is high. Treating repeated outputs as automatically independent would lose physical information needed by a full timing/noise model. We do not claim to have simulated a 100-photon joint wavefunction or a cosmological duration from the exact reservoir marginal.

## Consequence for the next major calculation

There are now two clearly separated branches to test:

1. **A stationary inelastic companion channel.** Derive its energy/frequency redistribution and arrival-time kernel rather than interpreting its energy balance as a coherent stretch. It still needs an independently justified explanation of the observed whole-event timing behavior; ordinary fixed travel conditions do not provide it.
2. **An evolving-field-mediated channel.** Specify the field state and local interaction, then derive its backreaction, energy receiver and full arrival-time map. The previously proposed time field is one possible direction, but a changing fundamental time constant is not uniquely required by the quantum-state test. A dynamical receiving field could be another direction. Neither is adopted as a completed theory here.

The next priority is a spatial, finite-event propagation calculation with an explicit evolving receiver and no hidden resetting/driver work, followed by the existing spectral-aging and brightness gates. Only after that mechanism supplies an independently constrained deposition/response law can the strong rotation fits be interpreted as evidence of photon origin. The mathematical gravity-amplitude degeneracy identified in the preceding review remains unresolved.

## Reproduction and scope

Run `python research_work/results/coherent-receiver-audit/run.py` and then `python research_work/results/coherent-receiver-audit/repeated_use.py`. `results.json` contains the energy, coherence and covariance tests; `beat-curves.json` stores synthetic beat signals in units with hbar=1; `repeated-use-results.json` contains the no-reset energy and correlation checks. All tests in this directory are synthetic mechanism diagnostics. No stellar or supernova holdout scores were opened, no global statistical significance was computed, and no theory-completion claim is made.
