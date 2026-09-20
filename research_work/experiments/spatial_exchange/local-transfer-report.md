# SE-LR3: direct transfer supplies field angular momentum

20 September 2026. Protocol681104d; executable1d2d4c8. All48 full Hamiltonian
derivative fixtures and all five individual evolution runs pass. The campaign
fails its unchanged spatial-circulation gate. Independent reconstruction
passes137 archive checks, maximum endpoint error3.25e-19. These statements
have different scopes: a faithful archive can record a failed campaign.

## Mechanism change

The curl-only coupling generates a circulating vector pattern with very
little canonical field angular momentum. In the axisymmetric continuum limit,
J_A,z=integral Pi dot(e_z cross A-partial_phi A) vanishes because both terms
cancel. This explains why circulation itself is insufficient evidence of spin
transport; it does not prove every nonlinear curl-only solution is exactly
axisymmetric.

We extended the local kinetic square to

    K = P - lambda A - epsilon (curl A) cross Q.

It adds +lambda K to the reciprocal A acceleration, while Q_t=K and the other
canonical equations follow the same Hamiltonian. The direct term introduces
no new spatial derivatives, so the earlier principal matrix/symmetrizer is
unchanged; full nonlinear stability remains unproved. Angular momentum uses
canonical P, not K. All excitation energy is included. No energy reservoir
outside the initially compact Q/P excitation or continuous forcing is added.

## Measured short-run transfer

All cases start with empty A/Pi and the same compact rotating excitation;
numbers below are at model time1, not galactic times.

| n24 case | A energy / initial total | A angular momentum / initial total | Circulation at r1.5 |
|---|---:|---:|---:|
| Curl only, epsilon2 lambda0 |0.051515%|1.43e-7%|0.000236445584|
| Direct only, epsilon0 lambda1 |21.25619%|10.18893%|2.98e-19|
| Combined, epsilon2 lambda1 |21.28249%|10.17994%|0.000212796199|
| Combined, n32 refinement |20.89778%|10.12161%|0.000187402463|

The direct-only case transfers spin while its signed equatorial circulation
is effectively zero. The combined case supports both. Neither quantity alone
establishes stronger attraction or useful lensing. No matter or light probe
was included, and outward angular-momentum flux was not measured.

Maximum relative energy drift across runs is1.03e-8. Maximum total angular
drift is2.91e-6 relative; near-edge energy stays below1.09e-10 of initial
energy. Removing the direct coupling reproduces the prior curl-only final
state exactly. Halving the timestep changes combined circulation by2.55e-7
relative and field angular momentum by5.24e-7, within the0.1% gates.

Under n24 to n32 refinement, field angular momentum changes0.58016%, within
its5% requirement. Circulation changes13.5504%, failing its5% requirement.
Two grids do not establish an asymptotic convergence rate; even the better
spin agreement requires further resolution and rotated-source tests. Earlier
failed archives and thresholds remain unchanged.

## Remaining path to the requested solution

First resolve circulation spatially and test common rotations; distinguish
spin stored near the source from angular momentum carried outward. Then
derive how ordinary matter supplies and maintains this local excitation with
reciprocal forces and a physical energy budget. Integrate the scalar/conversion
and matter/light sectors and recheck the full characteristic system before
using trajectories to claim gravity. Stability, sustained circulation, outer
orbits, lensing and held-out observations remain outstanding.

This is a proposed model ingredient built with established momentum-shift
and Hamiltonian mathematics, not a claim to have invented those constructions
or proven historical novelty. No old gravity prediction is used as a success
criterion; no dark matter or expanding background is included. Initial/final
fields and every-step ledgers are in local-transfer-v1. The auditor separately
reconstructs endpoint energy and canonical angular momentum and uses an
independent trilinear loop calculation. Intermediate extrema are checked
from recorded traces, not unarchived intermediate full fields.
