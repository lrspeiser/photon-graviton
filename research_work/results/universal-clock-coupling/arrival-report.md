# Does the accumulated shift survive arrival?

**Within the existing prescribed-field candidate, yes.** Thirty fixed calculations extend the ray beyond the changing region, through a later static clock contrast, to a detector with the same clock factor as the emitter. The downstream static region does not undo the previously accumulated matched-endpoint redshift. This is a kinematic result for the candidate, not a derivation of a self-consistent time field or a companion-energy receiver.

## Energy distinctions that matter

**Provenance: established Hamiltonian and local-clock relations, applied to our previously proposed profile; no originality claim.** For a right-moving photon, H = q c p and E_local = H/q = c p. H is the energy in the common reference-time description. Locally fixed observers measure E_local and d tau = q dt, with local light speed c.

Hamilton's equations give dx/dt = c q, dp/dt = -c p partial_x q, and dH/dt = (H/q) partial_t q. A region with no explicit time dependence therefore conserves H along the photon path. If q varies spatially there, E_local = H/q can still vary while the photon crosses it. At an exit with the same q as the entrance, E_local returns to its entrance value.

For a photon arriving from an earlier energy-losing region, that entrance value is already reduced. No rule in this Hamiltonian resets it to the original emitted energy. An observer who measures a different local frequency inside the static contrast has not thereby demonstrated that the earlier accumulated reference-energy loss was erased.

**Provenance: conditional consequence of established photon relations.** At matched endpoints q_emit = q_detector = 1, E_detector/E_emit = 1/(1+z). This remains the earlier net shift after any downstream static region with matching entrance/exit clock factors. A different actual detector motion or endpoint clock factor must still be included in the observational redshift-factor account.

## Fixed numerical cases

The earlier five profiles on x=0..1 are unchanged. A separate static contrast on x=3..4 is q_extra = contrast sin^2[pi(x-3)], with contrast 0, +0.2 or -0.2, and zero outside. This is an additive profile component; the two regions do not overlap. Source x=-1 and detector x=6 have q=1. All units are dimensionless and c=1. These are illustrative prescribed profiles, not fitted cosmic environments or an Earth gravity model.

Each profile is tested at launch times 0 and 2 for all three downstream contrasts. The result is independent of the later static contrast within the numerical tolerance. At launch zero:

| Earlier changing-region case | Received redshift after downstream region |
|---|---:|
| Static faster profile | 0 within numerical precision |
| Static slower profile | 0 within numerical precision |
| Faster profile relaxing toward normal | 0.012187812 |
| Slower profile slowing further | 0.009056963 |
| Slower profile relaxing toward normal | -0.018493622 |

The third and fourth cases preserve net redshift; the fifth preserves net blueshift. This again shows that a static time contrast alone does not produce cumulative redshift in this candidate. It also does not establish that every local observer along the path measures monotonically increasing wavelength: spatial clock-factor changes contribute to local comparisons.

## Verification and energy-receiver limit

The solver independently integrates photon momentum, reference energy, reference work and the sensitivity of arrival time to launch time. At nine checkpoints per case, q c p agrees with the independently integrated H, and H-H_initial agrees with integrated reference work. H stays constant after the changing region. Final wavelength stretch matches nearby-arrival stretch and the original saved result.

An initial coarse run crossing profile boundaries without segmentation failed the 2e-9 energy-identity gate. The integration was then split at the known fixed profile boundaries; the equations and acceptance gates were unchanged. All 30 segmented cases pass, and the maximum coarse/refined local-energy difference is 2.51e-14. This refinement resolves a numerical boundary-handling issue, not a physical energy mismatch. Original numerical results are not overwritten.

Integrating the work on the photon does **not** identify a companion field. Defining a complementary account as H_initial-H simply records how much the prescribed background would need to receive. It does not prove that such a receiver exists, obeys the desired no-loss propagation law, captures in wells or supplies extra gravity. In a completed interacting theory, reverse transfer rates must follow from the same dynamics rather than be set to zero solely because permanent redshift is desired.

Run `python research_work/results/universal-clock-coupling/arrival.py`. The saved results include all checkpoints, prior-result/protocol hashes and refinement comparisons. No astronomical parameter was fitted and no fresh catalog outcome was accessed. The new result settles arrival persistence for this limited candidate; it does not satisfy the goal's independent predictive-improvement or full physical-completion requirements.
