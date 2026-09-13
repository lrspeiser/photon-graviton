# Separating angular-momentum transport from cooling

13 September 2026. Exact-third inventory, phase-selector parameters and orbital endpoint profiles remain frozen. This is a conditional transport-budget calculation, not a new data fit.

## Main result

Escaping massless rays cannot provide the previously required individual-shell braking using only the calculated settling energy: the most favorable orbital ray-emission budget is about 1,400-1,500 times larger. Internal transfer to existing outer companion orbits is energetically less demanding in the tested fixed-potential scenarios. This supports investigating two separate roles: internal torques redistribute angular momentum, while outgoing radiation removes the remaining released energy.

The result does not establish a force or capture rate. It also does not exclude coherent wave torque, intrinsic-spin emission, exchanges among oppositely moving populations or additional energy sources; those are outside the direct orbital-ray bound.

## Direct outgoing rays

For freely escaping massless particles treated as rays emitted within radius r, orbital angular momentum satisfies

\[
|L_{\rm ray}|=|\mathbf r\times\mathbf p|\leq rE/c.
\]

Thus independently braking each contributing shell by positive angular-momentum magnitude dJ requires

\[
E_{\rm ray,min}=c\int\frac{dJ}{r_{\rm source}}.
\]

These are known energy-momentum and orbital-angular-momentum relations. Their application to companion cooling is proposed. Using the original source radius is favorable to radiation: the material settles inward, where the lever arm is smaller. Emission is maximally tangential; other ray directions cost more energy. This construction assumes no additional mechanical intermediary exchanging angular momentum.

The scalar sum dJ is not the net angular momentum of the whole spherical system. Random orbital planes can cancel vectors globally. The bound applies if the separately reduced orbital magnitudes are exported by rays rather than exchanged or canceled internally. A coherent multipole-wave mechanism has different angular-momentum accounting and is not constrained by this ray calculation alone.

| Quantity | Matter baseline I | Matter baseline II |
|---|---:|---:|
| Minimum direct-ray energy | 6.31e53 J | 6.80e53 J |
| Properly counted settling energy | 4.53e50 J | 4.52e50 J |
| Required / available | 1,394 | 1,504 |
| Required / total deposited rest energy | 1.69e-5 | 1.84e-5 |

The last row matters: direct radiation is not forbidden if additional stored energy is consumed. It is incompatible with attributing all of that output solely to the small settling-energy budget. Spending rest energy introduces reservoir depletion, different output power and a revised conservation ledger. No cosmic age or release duration is chosen.

## Internal torque receivers

In a fixed potential, circular orbits obey the established relation dE/dJ=Omega. This is why outward angular-momentum transfer and energy release are distinct in accretion dynamics; [Ogilvie's disc dynamics treatment](https://www.damtp.cam.ac.uk/user/gio10/dad2.pdf) provides the standard setting. Our spherical companions are not a demonstrated accretion disc, and no magnetic or viscous transport mechanism is imported by analogy.

We considered existing companion mass in three receiver regions, assigning equal specific angular-momentum increments within each region. Receiver orbit orientations are assumed matched to the donors as a favorable construction; opposing families may balance the total vector spin. We calculated both first-order work and finite orbit changes in the frozen final potential.

| Receiver region | Mean angular-momentum increase I / II | Mean radius increase I / II | Frozen-potential energy gain / settling budget I / II |
|---|---|---|---|
| 15-30 kpc | 3.05% / 3.52% | 3.96% / 4.52% | 28.8% / 30.7% |
| 30-60 kpc | 2.94% / 3.37% | 4.71% / 5.39% | 12.0% / 12.7% |
| 60-120 kpc | 3.61% / 4.14% | 6.55% / 7.52% | 4.48% / 4.77% |

These ratios suggest sufficient mechanical-energy headroom for an internal-transfer branch. They do not prove that the real exchanges occur. Moving the recipients changes their density and the gravitational potential; therefore these energy gains cannot simply be subtracted from the old global cooling budget and called the exact final answer. The total self-energy, orbital endpoint requirements and remaining radiation must be recomputed together. No receiver region is selected as the physical solution.

## What this means for the theory

The useful analogy is an orbital relay: inner companions pass angular momentum to outer companions, allowing some to settle while others move outward. Energy released by the combined process can then emerge as light or outgoing companions. This is distinct from asking radiation alone to provide the braking.

The next physical calculation is a coupled donor/receiver endpoint solution with conserved angular momentum, recalculated self-gravity and a nonnegative release budget. It must also test whether the outward recipient movement spoils stellar speeds or lensing. A transport coefficient, interaction strength and timescale remain to be derived. Isotropic arrival does not automatically generate the required orientation matching or torque.

This branch preserves the photon-source and one-third inventory ideas while changing the proposed settling dynamics. The other tracks remain active. In particular, nothing here solves photon redshift or its timing, verifies external energy supply, or identifies observed nuclear radiation as companion-powered.

## Reproduction and checks

[Protocol](torque-protocol.md), [code](torque.py) and `torque-results.json` retain both baselines and all three receiver choices. The old global cooling comparator uses the half-weighted companion self-energy, rather than summing single-orbit energies. Scalar receiver increments balance scalar donor decrements to 1e-12 tolerance under the orientation-matched assumption. Doubling the shell resolution changes the direct-ray bound by less than 0.0009%, and changes receiver energy-to-budget ratios by less than 0.00063 absolute. The remaining uncertainty is physical closure, not this numerical refinement.

The prior goal turn produced the orbital-support existence result; this turn quantifies the transport channels it requires. The full research goal remains active and incomplete.
