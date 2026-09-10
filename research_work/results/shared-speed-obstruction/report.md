# Why adjusting the shared-speed field alone cannot restore lossless travel

**Within the stated canonical propagation assumptions, the problem is structural rather than a poor parameter choice.** An evolving common nondispersive propagation speed and conserved companion Hamiltonian energy cannot both hold for every companion momentum in this one-sector description. This narrows the branch that should be pursued next; it does not rule out all fictional companion physics.

## The assumptions that make the conclusion valid

This argument applies to a differentiable ray Hamiltonian H_c(t,x,p), one canonical momentum p, and no additional evolving internal state. The companion energy to be retained is H_c in the same reference standard used for the photon/receiver ledger. Equal speed means dH_c/dp=c0/n(t,x) for a continuum of positive momenta. Lossless travel means that energy stays constant throughout each free trajectory, not merely that one specially chosen packet returns to its starting energy at an endpoint.

These are model assumptions, not established properties of gravitons. The conclusion changes if energy has a different operational definition, internal degrees of freedom participate, the propagation is dispersive, or the dynamics cease to be this canonical ray model. Such changes need explicit equations and a complete energy ledger.

## Conditional derivation using established Hamiltonian calculus

Integrating the equal-speed condition with respect to momentum gives

`H_c(t,x,p) = c0*p/n(t,x) + B(t,x)`.

Hamilton's equations imply that along a ray, `dH_c/dt = partial_t H_c`: the spatial and momentum derivative terms cancel. Consequently

`dH_c/dt = -c0*p*n_t/n^2 + B_t`.

To set this to zero at two distinct momenta at the same event requires

`0 = -c0*(p1-p2)*n_t/n^2`.

Thus n_t must vanish there. For a genuine evolving propagation field, one momentum-independent offset B cannot conserve every companion energy. Spatial gradients do not remove this condition because they already enter the cancellation in Hamilton's equations.

Provenance: this is a conditional deduction from known canonical mechanics and elementary calculus, not a certified new theorem or a first-principles derivation of the existence of time. It concerns this specific meaning of shared speed and energy. It does not prove that every nonexpanding model must fail.

## Executed loophole checks

The numerical controls use a homogeneous n increasing from 1 to a chosen final value. With c0=1 and no spatial force, p is conserved. We integrate the explicit time-dependent work and verify it equals each energy change. This is a prescribed-background diagnostic, not a closed model of the source performing that work.

For n increasing from 1 to 2, an offset B=2*(1-1/n) exactly preserves the p=2 packet:

| Initial momentum/energy | Final energy with this offset | Energy change |
|---:|---:|---:|
| 1 | 1.5 | +50% |
| 2 | 2.0 | 0% |
| 4 | 3.0 | -25% |

All three retain the same group speed 1/n. Protecting one chosen energy therefore creates gains or losses at others. The program also checks final n=1.01 and 1.3; the same problem persists at smaller changes. No astronomical energy distribution is fitted.

Allowing B to depend on momentum can cancel every energy change: B=p*(1-1/n) makes H_c=p. But then dH_c/dp=1, and the companion no longer follows photon speed 1/n. The numerical derivative checks that consequence explicitly. The workaround removes the desired speed behavior as well as the energy drift.

Finally, common clocks q=1/n make the homogeneous local energy H/q constant for the unmodified H=p/n. They do this for photons too. That returns the earlier measured-redshift cancellation rather than preserving photon redshift while protecting companions.

## Consequence for the research direction

We should not spend further parameter-only searches claiming that the existing one-sector evolving-index Hamiltonian can simultaneously supply all-momentum lossless companions and the same variable propagation speed. That combination conflicts with its assumptions. The prior autonomous wave-retention failure is a concrete realization of this more general restriction.

A viable next branch must state which assumption changes. The most relevant possibilities remain a distinct photon-conversion interaction with freely propagating companions, or additional internal receiving degrees of freedom whose energy, motion and detector coupling are specified together. Neither is an established solution: the existing stationary-converter timing and coherence restrictions still apply, and extra degrees of freedom cannot supply unaccounted work.

The distinction matters in plain language: we need an explanation for why light loses observable energy while its companions do not. Merely assigning both the same slowing propagation law does not provide that difference. A successful gravity fit also cannot repair this missing propagation physics by itself.

This finding closes a direction of parameter tuning, not the overall goal. The real-data gravity comparisons, bulge/orbit likelihood, redshift and transient tests remain required. No observational values were changed, no holdouts were opened, and the total photon-supply budget remains deferred.

## Reproduction and related evidence

Run `run.py`. `results.json` records all assumptions, 18 offset controls, three clock controls, the momentum-dependent compensation check and their energy-work results. These controls illustrate the analytic restriction; their success is not observational validation.

- [Autonomous wave retention](../companion-retention-test/report.md): energy transfers from the traveling perturbation into background motion.
- [Equal-speed field proposal](../companion-speed-consistency/report.md): limited positive-energy wave-speed construction.
- [Matter-clock closure](../matter-clock-closure/report.md): common-standard cancellation.
- [Coherent receiver audit](../coherent-receiver-audit/report.md): why an energy-conserving converter still needs a specified temporal receiving state.
