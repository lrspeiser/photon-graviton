# Local capture energy and momentum ledger

The preceding goal turn made scientific progress by calculating post-capture radial redistribution. This turn checks the previously unmodeled stopping step. All six objectives remain open.

## Assumptions and provenance

Assume a traveling companion is massless with p=E/c in a local inertial frame, and a receiver initially at rest absorbs it without outgoing products. These are optional companion/interaction assumptions. Relativistic energy-momentum conservation and the invariant-mass relation are known mathematics, not new physics claimed by this project. Photon momentum and absorption recoil are reviewed in [OpenStax, The Particle-Wave Duality Reviewed](https://openstax.org/books/college-physics-2e/pages/29-8-the-particle-wave-duality-reviewed). Our derivation applies the same kinematics to the hypothetical companion. It does not derive an allowed quantum transition, capture rate, coupling to well depth, or storage lifetime.

Let M be the inertial mass of the actual receiver that takes up recoil. It cannot be set equal to the entire cluster mass without specifying how momentum is transferred to that system. For incident total energy E and net incident momentum magnitude j E/c, where 0<=j<=1:

    M_final c^2 = sqrt((M c^2+E)^2-j^2 E^2),
    Delta_internal = (M_final-M)c^2,
    K_recoil = E-Delta_internal,
    v/c = j E/(M c^2+E).

All incoming energy is retained in the closed receiver system by assumption, split between internal/rest energy and bulk kinetic energy. Recoil energy has not disappeared and can also contribute to gravity; whether it escapes the target well, thermalizes, or is radiated later is a separate transport question. These formulas are a local ledger, not a global gravitational binding-energy calculation.

## Single-direction absorption

With e=E/(M c^2) and j=1,

    Delta_internal/E=2/(sqrt(1+2e)+1),
    K_recoil/E=1-Delta_internal/E,
    v/c=e/(1+e).

| Incident energy / receiver rest energy | Internal/rest gain as fraction of incident energy | Recoil fraction |
|---:|---:|---:|
| 0.000001 | 0.9999995000 | 0.0000005000 |
| 0.001 | 0.9995004994 | 0.0004995006 |
| 0.1 | 0.9544511501 | 0.0455488499 |
| 1 | 0.7320508076 | 0.2679491924 |
| 10 | 0.3582575695 | 0.6417424305 |

These ratios are synthetic. A massive enough receiver can store almost all the incident energy as internal/rest energy while moving very little. Momentum conservation therefore does not by itself prevent efficient storage. Conversely, an unspecified receiver mass is insufficient to predict recoil or retention.

In a weak external well, for an initially stationary receiver and no subsequent forces beyond that well, requiring post-capture speed below a prescribed escape speed gives e<b/(1-b), where b=v_escape/c, for j=1. This is an approximate Newtonian escape test applied to the local recoil velocity, not a relativistic binding calculation. The example b=0.001 gives e<0.001001001. These are dimensionless examples, not measured cluster constraints. A bound center of mass alone does not guarantee that internal excitation cannot later escape as radiation.

## Arrivals from all directions

Two equal opposite incident beams have j=0. If a common closed receiver retains both, its final bulk momentum is zero and its invariant mass increases by exactly E/c^2. Individual absorption sites can nevertheless recoil, heat, or emit energy. Statistical isotropy of a cluster-wide bath does not mean every local capture event has zero momentum. A transport model must distinguish these scales.

An isolated single massless companion cannot turn into a lone massive particle at rest while conserving its four-momentum: the initial invariant mass is zero. A receiver, additional products, or a second incoming excitation can resolve the kinematics, but requires a specified interaction. Two exactly opposite equal massless excitations have nonzero combined invariant mass and can kinematically supply a stationary massive final state if such a reaction is permitted. This is not evidence that the reaction exists or occurs at the needed rate.

## Consequence for current branches

The falling-particle calculation remains conditional. Absorption into an existing receiver adds energy to that receiver; it does not automatically create the separate at-rest particle used in redistribution.py. We must choose and derive that extra interaction before treating the orbit calculation as a complete capture model.

Receiver excitation is closer to the proposed stored well-level interpretation, but its physical degrees of freedom, spatial support and decay law are still unspecified. If ordinary matter is the receiver and retains the energy, the extra rest energy initially follows that matter. An extended halo is therefore not established merely by a kinematically allowed absorption. The opacity profile and receiver distribution must refer to the same physical system.

## Verification

capture-ledger.py and capture-ledger-results.json test 21 cases: seven energy ratios from 1e-12 to 1000 and three net momentum fractions. The closed-form result is checked against direct sums of a stationary receiver and two opposite massless four-vectors, retaining the energy and momentum equations separately. All scaled residuals are below 2e-13. Three synthetic escape-speed thresholds are bracketed with incident energies 1% below and above the limit. These tests establish algebraic conservation under the stated assumptions, not the interaction's existence, permanence, or an observational fit. No holdouts were opened.
