# SE-JE1: empty-field generation works; joint bend accuracy fails

20 September 2026. Protocol7f05877; runner9144a75. All eight declared runs
finish. Individual energy/angular/cone/edge gates pass, as do timestep and
rotation comparisons. The campaign fails its5% spatial bend requirement:
n24/n32 changes20.0957%. Independent endpoint/decision reconstruction passes
99 checks with maximum discrepancy2.67e-15; it confirms the failed campaign.

## One model generates fields and moves matter/light

All fields, including Q/P, start at zero. Six moving massive particles and
one backreacting massless particle supply the complete initial Hamiltonian
energy. No internal excitation is inserted. Particle-vector coupling drives A,
and the direct coupling then excites the local transfer sector. The free
control retains empty fields and straight trajectories. With transfer disabled,
Q/P remains zero. With transfer enabled, K/Q energy becomes nonzero.

| Case | Endpoint photon bend, model radians | Final K/Q energy | Particle Hamiltonian decrease |
|---|---:|---:|---:|
|Fast baseline|-0.008719965346|0|0.0257721082|
|Fast combined|-0.008717736265|4.99128e-6|0.0257715456|
|Slow baseline|-0.008329972003|0|0.0234444558|
|Slow combined|-0.008329966874|4.33484e-10|0.0234444557|

Transfer decreases the coarse bending magnitude by about0.02556% in the
fast pair and0.00006157% in the slow pair. Neither effect is resolved against
the current spatial error, and no desired-sign requirement was used. The
source momenta are.2 and.0007 in reference-speed units; these are controlled
model fixtures, not calibrated galaxy models. In particular, slow sources
need not remain on circular orbits in this starting configuration.

K/Q energy includes the kinetic-shift interaction and is not separately
conserved radiation energy. Particle energy decrease includes field-dependent
Hamiltonian energy; it is not a measured loss of rest mass or calibrated
stellar luminosity. Bare masses stay fixed. Total Hamiltonian accounting
holds, but the physical fuel interpretation and longevity remain unresolved.

## Convergence and angular transfer

| Check | Measured | Requirement | Result |
|---|---:|---:|---|
|Fast combined timestep bend change|0.00000580%|<0.1%|Pass|
|Fast combined n24/n32 bend change|20.0957%|<5%|Fail|
|Rotated final velocity direction / bend|1.48149%|<2%|Pass|

The n32 combined bend is-0.010910217105, compared with-0.008717736265
at n24. The point interpolation width shrinks with the grid, so both source
representation and evolved fields change. More resolution and a matched
on/off comparison are needed before interpreting a small transfer effect.

In the fast combined case final A angular momentum along the source axis
is3.36e-8 and Q angular momentum is-6.60e-8, compared with initial total
angular momentum0.8399. Those tiny signed transfers differ fundamentally
from the earlier experiment's approximately10% transfer of a deliberately
inserted internal excitation. That earlier percentage cannot be carried over
to matter-driven empty-field generation. The slow case's generated K/Q energy
is over four orders of magnitude smaller than in the fast case.

These results identify the next physical issue: a mechanism intended to affect
slow stars and cluster lenses needs stronger matter-driven generation or a
different field response without adding an unaccounted energy reservoir.
Simply assuming the earlier internal excitation does not meet that need.
Any changed coupling requires a new Hamiltonian/source derivation and declared
tests. No candidate is accepted because it matches an older gravity law.

## Scope and reproducibility

This is the first joint-candidate empty-field campaign, not a demonstration
of stable orbits, sustained swirl, complete lensing, polarization, physical
source fuel or observed galaxy/cluster agreement. Light is an endpoint ray,
not a lens bundle. The finite-time cone diagnostic checks interpolated
particle velocities; full nonlinear propagation/stability remains unproved.

Full initial/final states and every-step particle states/ledgers are preserved
in joint-evolution-v1. The audit independently reconstructs field and particle
energies, K/Q energy, canonical total momentum/angular momentum and endpoint
photon velocity using full-grid cardinal basis weights. Intermediate extrema
are checked from saved traces, not unarchived intermediate field states.
No dark matter, expansion or external forcing was included. All twelve goal
requirements remain active.
