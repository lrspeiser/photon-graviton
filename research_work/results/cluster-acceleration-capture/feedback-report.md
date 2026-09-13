# Conditional feedback of deposited gravity on capture

## Result

Allowing deposited energy-equivalent mass to deepen the well changes both later capture and the stored profile. In the A=1 comparison, when accumulated deposit mass reaches the original ordinary mass, capture power is 2.228 times its initial value. Half the stored mass lies inside 4.100a, compared with the initial injection half-radius about 3.287a. A fixed-well injection profile therefore cannot simply be multiplied by an arbitrary mass and assumed unchanged.

| Deposited / ordinary mass | Current capture / initial capture | Half-stored-mass radius / a |
|---|---:|---:|
| 0.1 | 1.111 | 3.366 |
| 0.5 | 1.582 | 3.706 |
| 1.0 | 2.228 | 4.100 |

These are imposed dimensionless evaluation points, not measured cluster mass ratios or a fitted source history. This is a controlled feedback calculation under immobile-storage assumptions, not a supported halo solution or complete energy-conserving relativistic evolution.

## What changes in the equations

Keep the ordinary Plummer mass profile and the candidate universal coefficient chi fixed, with initial A=1. At every integration step use the known spherical field equation

    g(r,t)=G [M_b(<r)+M_d(<r,t)]/r^2
    kappa(r,t)=chi g(r,t)^2.

Recompute upstream attenuation from this opacity and the imposed external isotropic companion bath. Deposit local absorbed energy into the local shell using the conditional ordinary mass-energy identification:

    d rho_d/dt=q(r,t)/c^2.

No radial redistribution or orbital motion is applied. The physical time unit is M_b c/(u_infinity a^2); the dimensionless exposures in the output do not select a universe age. The incoming bath is assumed established and transport recomputed quasi-statically. Applying this approximation physically requires evolution slow compared with relevant propagation times; no coefficients have been established to verify that ordering for real clusters.

This holds deposited mass in place without explaining the necessary stresses. The preceding support analysis already showed why ordinary isotropic particles cannot do that for the initial profile. Accumulation also changes gravitational binding energy; stress energy, support work and binding-energy corrections are not included by the simple rho_d=E/c^2 update. The numerical absorbed-energy ledger must not be mistaken for a complete gravitational energy budget.

## Numerical construction and checks

The spherical domain extends from zero to 1000a with 128 or 256 radial shells, logarithmically spaced outside 0.001a. Opacity is piecewise constant in shells. Exact upstream chord lengths through each shell give attenuation at shell volume midpoints and 48 angular quadrature nodes. Independent full-chord absorption areas use impact quadrature split at every shell boundary.

Evolution is parameterized by accumulated mass, with physical exposure integrated alongside it. A two-stage Heun update deposits each positive mass increment according to the normalized instantaneous capture profile. This is equivalent to the conditional time evolution above and guarantees the specified absorbed mass is allocated once. It is not a way of deriving the incoming energy supply. Mass steps of 0.02 and 0.01 M_b are compared.

The initial capture coefficient is 6.8815, versus 6.9085 from the earlier infinite-domain analytic-opacity calculation, a difference about 0.39%. Doubling radial shells changes subsequent power by approximately 0.60-0.63%; halving the mass step changes power by less than 4e-6 relative and exposure by less than 2.3e-5. Independent volume/chord capture integrals agree within 7.1e-6 relative in the 256-shell, 48-angle run. Deposit allocation errors are below 1.2e-16 M_b.

A separate combined check doubles the outer boundary to 2000a and angular nodes to 96 at 256 shells. Final capture power changes by 0.0677% and half-mass radius from 4.100a to 4.105a. Because boundary and angle settings change together, this is a sensitivity check rather than an isolated angular convergence estimate or rigorous truncation bound. All runs are retained in feedback-results.json.

## Implications

The external-only collecting-area proposal has positive feedback under this candidate rule: deposits deepen the well, increasing capture and shifting subsequent deposition. We have not demonstrated runaway, stability or saturation. Receiver competition and depletion can change the external bath, while orbital/field support can move the deposits and alter the result. Those effects cannot be omitted when comparing a mature halo with data.

The next model must supply support and account for the shared bath and gravitational energy before treating this evolution as physical. The earlier frozen-profile lensing thresholds are diagnostics, not self-consistent thresholds at these accumulated masses. All equations used for field, absorption and numerical integration are established tools; the capture and immobile-storage prescriptions remain explicit hypotheses. No observational holdout was opened and all six objectives remain open.

Run `python research_work/results/cluster-acceleration-capture/feedback.py`.
