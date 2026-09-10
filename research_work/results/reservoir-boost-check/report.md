# The recoil repair does not complete a relativistic reservoir

**The current wave/reservoir/support model fails a necessary center-of-energy condition if interpreted as a closed ordinary relativistic system.** This sharpens the earlier stress-tensor caveat: conserving Hamiltonian energy and canonical momentum is not enough. The discrepancy survives the moving-support repair and equals the outgoing reservoir energy in the normalized late-state balance.

This excludes that closed relativistic interpretation of the displayed subsystem. It does not exclude every nonexpanding conversion model, a fully specified additional background sector, or an explicitly non-Lorentz-invariant theory. Those alternatives require new physical equations and tests; they are not established by this result.

## Known physical requirement

For an isolated system with localized conserved energy and a symmetric relativistic stress-energy tensor, energy flux equals c squared times momentum density. The standard relationship is explained in [Feynman's field-energy and momentum discussion](https://www.feynmanlectures.caltech.edu/II_27.html). With negligible boundary flux, energy continuity implies

`X_E = integral x*u dx / E_total`,

`dX_E/dt = integral J_energy dx / E_total = c^2*P_total/E_total`.

This is known relativistic mechanics, not a project postulate or a newly discovered constraint. The test uses c=1 in the ordinary propagation region after the pulse has left the conversion profile.

## Apply it after all local conversion has stopped

Use the actual outputs of the coupled moving-support calculation. Normalize energies and momenta by the initial wave energy. Let E_gamma denote outgoing electromagnetic energy, E_R the outgoing reservoir energy, E_support the support's total energy including its rest energy, and K its momentum.

The support's velocity is K/E_support. Both outgoing wave and reservoir energy profiles translate at speed one. Their integrated energy flow is therefore

`integral J_energy dx = E_gamma + E_R + E_support*(K/E_support)`

`= E_gamma + E_R + K`.

But the model's conserved canonical momentum in the homogeneous-clock late state is

`P_canonical = E_gamma + K`.

The reservoir's canonical translation term is `-P*s_x`, which vanishes on s=t. Consequently

`integral J_energy dx - P_canonical = E_R`.

This is an exact structural discrepancy whenever conversion produces positive outgoing reservoir energy in this branch. Tighter numerical integration cannot remove it. A finite support mass does not supply the missing term because its ordinary energy-flow and momentum contributions already agree.

## Numerical result using the archived recoil runs

| Support rest energy / initial wave energy | Predicted center-energy speed | Required P/E speed | Relative discrepancy |
|---:|---:|---:|---:|
| 10 | 0.1032861064 | 0.0909090909 | 13.6147% |
| 100 | 0.0112768776 | 0.0099009901 | 13.8965% |
| 10,000 | 0.0001139175 | 0.0000999900 | 13.9289% |

Speeds are fractions of the reference c=1. The support rest energy is included in both denominators. Its large value reduces the absolute speed difference but does not make the integrated energy-flow deficit disappear. The tighter support-ratio-10 run reproduces the same result.

An independent finite-time construction of the late energy centroid verifies the derivative. Arbitrary starting centroids are used because their offsets do not affect the velocity when each outgoing profile translates rigidly; no new wave-shape simulation is claimed. A no-conversion control satisfies the identity. Every input is checked against the earlier calculation's hashes.

## Why relabeling the reservoir does not fix the model

Assigning ordinary forward massless momentum E_R to the reservoir would make its energy current and momentum agree locally. But it would also increase the outgoing total momentum by E_R while leaving the previously calculated support impulse unchanged. The present dynamics already balances incoming momentum against E_gamma+K. A physical repair must alter the coupling and its reaction forces, not simply add an interpretation after the calculation.

Likewise, relabeling canonical momentum as something other than total physical momentum is an admission that additional physical stress/momentum accounting is needed. A complete background or medium could contribute it, but those contributions must be derived. Global improvements of a local stress tensor cannot be assumed to change the total conserved charge of a localized system without boundary or additional-sector terms.

The calculation does not assert that every possible reservoir must obey this model's canonical identification. It establishes that the presently displayed variables, energy current and support are not already a complete isolated ordinary relativistic realization.

## Consequence for the research program

The one-dimensional equations remain valid as a postulated preferred-frame canonical model, and their energy-transfer and endpoint-frequency numerical controls remain valid within that model. They should not be promoted as a completed physical graviton mechanism or used as such to interpret gravity/lensing fits.

The next physical revision must either supply a complete background/matter sector with a consistent stress-energy tensor, or use a receiver interaction whose transported energy and physical momentum are jointly specified from the start. It must then be retested for energy retention, observable clocks, endpoint redshift and timing. Adding free fit coefficients cannot repair this structural identity.

No observational holdout, galaxy coefficient, total photon budget or cosmic-age assumption is changed. This is a failed consistency gate that changes which claims are defensible, not evidence that observations themselves are unreliable or that all companion hypotheses are impossible.

## Reproduction

Run `run.py`. `results.json` records each support case, both centroid velocities, the flow/momentum deficit, the unchanged-data provenance and the scope assumptions. The prior canonical-conservation calculations remain archived without modification.
