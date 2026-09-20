# SE-TB1: integrated outgoing energy balances regional loss

20 September 2026. Protocol23a5358; runner9405ac4. All three replay runs and
the timestep comparison pass. All39 independent audit checks pass; endpoint
regional energies reconstruct exactly at machine precision. The final fields
exactly reproduce the previously archived combined/time/space LR3 cases:
adding the flux diagnostic did not change the dynamics.

The RK4 stages accumulate the compatible discrete outward face power. The
maximum discrepancy in regional energy change plus accumulated outward
energy, divided by initial total energy, is1.002e-8 at n24dt.02,5.23e-10 at
n24dt.01 and6.13e-10 at n32dt.01, all below the declared1e-6 threshold.
Halving dt changes net outward energy at radius1.5 by6.95e-7 relative,
well below0.1%.

| Grid, dt | Net outgoing through r1.2 / initial total | Through r1.5 | Through r2 |
|---|---:|---:|---:|
|n24,.02|8.320579%|1.002997%|0.011076%|
|n24,.01|8.320583%|1.002998%|0.011076%|
|n32,.01|8.112187%|1.107466%|0.006200%|

These are signed net energy transfers from time0 to1 through staircase
boundaries of the spherical grid masks. They are not gross radiated energy,
irreversible loss, physical luminosities or a prediction for galactic times.
Total local-sector energy includes both A and internal Q/P excitation;
separate A-sector energy is not conserved. Energy could return later.

Spatial dependence remains evident, particularly in the small signal at
radius2. This campaign verifies flux accounting, not converged transport
magnitudes or a continuum wavefront. The later n40/48/64 transfer campaign
did not store these RK-stage powers and is not silently retrofitted with
complete flux histories. No angular-momentum surface flux is measured here.

Initial/final full fields, completed-step regional energies, cumulative
transport and all stage powers are in transport-budget-v1. The independent
audit reconstructs endpoint regional energies and sums stage powers using
a separate accumulation. It does not independently reconstruct every stage
power because intermediate full fields are not archived. The face-flux
formula has its separate DF1 identity tests.

The result establishes conservative energy transport in this finite toy
excitation. It does not establish how ordinary matter produces the excitation,
whether fuel lasts, how matter/light move in the resulting field, or whether
one shared law fits galaxies and clusters. No dark matter, expansion or
old-gravity agreement enters the acceptance tests. Established numerical and
Hamiltonian methods are credited in the protocols; no novelty claim is made.
