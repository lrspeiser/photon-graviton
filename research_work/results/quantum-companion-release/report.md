# Can captured photon energy leave as companions?

The finite quantum receiver can release energy into stipulated outgoing massless modes while conserving energy and forward momentum. That alone does not complete the proposed mechanism. Three explicit release rules expose different costs: exact frequency-resolved loss tags destroy the photon coherence in this construction; a constant release borrows and returns receiver energy; releasing the receiver's entire state consumes its initial energy reserve. None is adopted as a complete explanation of observed redshift or supernova timing.

In plain language, moving energy out is possible in the bookkeeping. The unresolved question is how to do that while preserving the information in the light and accounting for the state of whatever temporarily receives the energy. Balancing the average energy for one chosen spectrum does not guarantee balance for different light sources or an unchanged receiver after many uses.

## Status and reproducibility

Run `python research_work/results/quantum-companion-release/run.py` from the repository root. A fresh execution passed all embedded assertions: eight release cases, four changed-input-spectrum controls and two 100-use histories. `results.json` records the script and upstream hashes. This is synthetic finite-state analysis, not an astronomical fit. No held-out observations were opened, no gravity or redshift coefficient was tuned, and no uncertainty from data is represented by the numerical tolerances.

The previous turn's catalog work was concrete progress: it reverified the local catalogs and corrected the inventory's outdated orbit-library description. This calculation addresses a separate major physical gap, rather than adding another fit to exposed galaxy data.

## What is assumed

All energies are dimensionless test units, with c = hbar = 1. Let S be a chosen frequency-stretch factor, and let delta(E) = (1 - 1/S) E be the photon's lost energy. These are a definition and an imposed frequency map, not a derived rate of conversion over distance. S = 2 and 1.5 are test choices, not fitted redshifts.

The upstream hypothetical conversion is

    |E, R, P_R> -> |E/S, R + delta(E), P_R + delta(E)>.

The receiver has the stipulated spectrum P_R = R - U, U = 100, and R >= 50.01. Its invariant mass squared is U(2R - U) > 0, using the known relativistic mass-shell identity. Different receiver basis states have different internal masses; this is not the dispersion relation of a single fixed-mass particle. Initial receiver energy has mean 100 and Gaussian coherent widths 0.02, 0.1, 1 or 4. That prepared receiving state is a physical resource requiring explanation, not empty space supplied at no cost.

Input photon bands are 4–6 for S = 2 and 4.5–6 for S = 1.5. The receiver absorbs the different losses coherently. Conditional receiver states overlap, allowing off-diagonal photon frequency terms to survive. Those overlaps determine the finite-band Fourier pulse profile in the upstream calculation.

Known quantum mechanics supplies the trace, inner-product and conservation identities used below. The particular spectra and conversion/release maps are hypothetical constructions for this project; no uniqueness or literature novelty is established. Energy conservation and coherence impose distinct requirements, as discussed in [Lostaglio et al., Quantum coherence, time-translation symmetry and thermodynamics (2015)](https://arxiv.org/abs/1410.4572). That paper is background for the distinction, not a derivation or endorsement of this receiver.

## Three release rules and their consequences

### 1. Exact photon-loss tags

Consider a different direct conversion with a vacuum companion input and an unchanged receiver:

    |E, R, vacuum> -> |E/S, R, C_delta(E)>.

The companion has sharp energy and forward momentum delta(E). Energy and momentum balance branch by branch. Distinct input energies produce orthogonal companion-energy states. Tracing the unobserved companion gives

    rho_photon[n,m] = psi[n] psi[m]* <C_delta(m)|C_delta(n)>
                    = |psi[n]|^2 delta_nm.

This is a conditional derivation from known partial-trace mathematics, not a general theorem against companion models. It removes the frequency coherence and makes the selected finite Fourier profile uniform over its recurrence period. Computed photon purities are 0.0402993 and 0.0846284 for the two bands. A frequency distribution remains; what disappears in this construction is the localized Fourier profile.

One cannot append a receiver-only operation to the prior coherent conversion that both restores its initial state and writes perfectly distinguishable loss labels into companions. Such an operation would map overlapping conditional receiver states to orthogonal states, violating the known preservation of inner products by an isometry. At receiver width 4, the maximum overlap-matrix mismatch is about 0.9999992 for both bands.

This result assumes sharp loss tags, the specified input ancilla state, and the displayed maps. It does not examine all finite-band companion packets, correlated inputs, nonlocal resources, field theories or detector models. A supernova light curve is not a single coherent optical wavepacket extending over days, so this calculation does not itself predict or refute supernova event timing.

### 2. Fixed-energy release with an evolving receiver

The hypothetical receiver-only release uses one frozen companion energy L:

    |R, vacuum> -> |R - L, C_L>   if R - L >= 50.01;
    |R, vacuum> -> |R, vacuum>   otherwise.

Companions have E_C = P_C = L. The map has an energy/momentum-conserving swap completion between vacuum and occupied companion sectors. The lower boundary retains energy rather than creating a negative-energy receiver. It is a state map, not an interaction Hamiltonian or a demonstrated flight trajectory.

Because the release acts only on the receiver and companion sectors, tracing both sectors leaves the photon reduced state unchanged. This is a known unitary/partial-trace identity, not a newly discovered redshift effect. Direct matrix and Fourier-profile checks confirm it here.

L is fixed to the mean loss of each original centered band, not adjusted to the changed-spectrum controls. For the broad receiver (width 4):

| Quantity | S = 2 | S = 1.5 |
|---|---:|---:|
| Mean photon energy before conversion | 5 | 5.25 |
| Mean photon energy after conversion | 2.5 | 3.5 |
| Companion energy L per emission | 2.5 | 1.75 |
| Receiver energy variance initially | 16 | 16 |
| Receiver variance after one use | 16.0049 | 16.0011111 |
| Receiver variance after 100 uses | 16.49 | 16.1111111 |
| Weight of branches borrowing receiver energy | 47.15% | 44.02% |
| Mean borrowed energy per input photon | 0.0278784 | 0.0131979 |
| Mean extra retained energy per input photon | 0.0278784 | 0.0131979 |

The receiver's mean energy stays at 100 for these centered inputs to numerical precision, but its state does not reset. The one-use no-emission probabilities are approximately 3.84e-37 and 1.68e-37, not identically zero. The equality of averaged photon loss and companion output is therefore numerical at this precision, not an exact boundary-independent identity.

Away from the boundary, standard random-increment bookkeeping gives the conditional formulas

    mean(R_N) = mean(R_0) + N (mean(delta) - L),
    Var(R_N)  = Var(R_0) + N Var(delta),

for independent fresh input photon energy distributions in this release construction. The 100-use calculation propagates receiver diagonal probabilities, including the low-energy boundary, without resetting the receiver. It agrees with the interior variance formula to the asserted 1e-7 tolerance. It does not propagate the full multiphoton quantum state or establish repeated output pulse quality. Correlations between earlier photons and the receiver are not claimed absent.

Changing the input band center by -0.2 or +0.2 while leaving L frozen changes receiver energy by approximately -0.1 or +0.1 per photon for S = 2, and -0.0666667 or +0.0666667 for S = 1.5. Thus a release balanced for one source spectrum drains or charges the receiver for another. Long-term source mixing, boundaries and state preparation must be included before asserting a sustainable void mechanism.

This rule does NOT satisfy exact transfer of each photon's loss into its own companion. It is an explicit alternative diagnostic, not a silent replacement of the user's rule. Even an unchanged average receiver energy does not establish the requirement that energy simply travels through voids with no enduring receiving-sector change.

### 3. Release the whole receiving state to a low-energy anchor

Another hypothetical isometry transfers the receiver's spectrum into the companion:

    |R, vacuum> -> |R_anchor, C_(R-R_anchor)>,
    R_anchor = 50.01.

It preserves the photon reduced state by moving the receiving-state information into the companion. It also emits the receiver's initial energy above the anchor, not just the incoming photon's loss. Mean companion energies are 52.49 and 51.74, including 49.99 units drawn from the initial receiver in either case. The anchor momentum is -49.99, giving speed approximately -0.9996001 in these units. Four-momentum remains balanced because this large recoil is counted.

This resets the receiver to a low-energy sharp state, not to the initial broad coherent state. Restoring the initial state would require energy, momentum and coherence preparation. Calling this a reusable release of only redshift energy would conceal those costs.

## Verification and limits

Single-use total energy errors are at most 2.89e-14 across all 12 centered/sensitivity cases. Momentum, basis-level release differences, photon reduced-state equality, profile equality, normalization and positive invariant mass assertions pass. The largest recorded 100-use energy-balance error is 3.47e-12. Tiny boundary probabilities are retained; the finite upper cutoff has an explicit lost-probability bound in the repeated calculation.

These are checks of the constructed mathematics. None establishes local causal propagation, a graviton identity, absence of receiver gravity in voids, a distance-dependent alpha, protected atomic clocks, real source/detector arrival times, capture into galaxy wells, deposit support or lensing. The total astronomical energy-supply calculation remains deferred, not passed.

## Consequence for the full research goal

No complete redshift-and-gravity mechanism emerges from these three release rules. The useful advance is a concrete resource ledger: an outgoing companion can preserve photon coherence only under the tested alternatives with an explicitly changed receiving sector or consumed receiving resource; the exact sharp-tag control preserves energy but loses this model's photon coherence. These restricted examples are not an exhaustive impossibility proof.

The next physical milestone is a finite-history local interaction that derives the incoming-to-outgoing phase and energy maps, evolves its receiving sector without an unexplained reset, and releases companions with accounted stress and momentum. It must then predict event timings with specified source and detector clocks. A further free adjustment of L or S would not meet that milestone.

In parallel scope, the stellar task still requires a stable, sufficiently complete orbit-population representation and a survey-selection/error-aware likelihood before evaluating reserved observations. A good empirical extra-potential fit would not identify its photon origin while the production/capture response remains free. Neither the mechanism audit nor the catalog inventory constitutes a strong case for the full theory; the original goal remains active.
