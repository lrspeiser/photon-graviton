# SE-TB1: time-integrated local energy transport

Declared 20 September 2026 before calculation. Replay unchanged combined LR3
equations and initial excitation at n24dt.02,n24dt.01,n32dt.01, L8,T1,
epsilon2,lambda1. Preserve the earlier archives. Add no source or forcing.

For spherical grid masks r<1.2,1.5,2 accumulate the compatible SE-DF1 outgoing
face power at all four RK4 stages using dt*(P1+2P2+2P3+P4)/6. Record each
stage's power, each completed-step region energy and accumulated transport,
initial/final full fields. Record instantaneous endpoint powers as well.
The region's discrete link-energy assignment is unchanged.

Require maximum over all times/regions of
abs(E_region(t)-E_region(0)+integrated_outward_power)/initial_total_energy
<1e-6. Require final fields reproduce the corresponding LR3 combined/time/
space archive to1e-12. Require n24 timestep refinement to change net outward
energy at r1.5 by<.1%, denominator max(abs(finer),1e-12). Report n24/n32
differences as unresolved spatial diagnostics, not a declared spatial pass:
hard masks approximate different staircase surfaces at different resolutions.

An independent audit must reconstruct endpoint region energies from full
arrays and accumulate archived stage powers independently. Intermediate full
stage fields are not archived, so the audit does not independently reconstruct
every stage's power. Their underlying formula has separate DF1 identity tests.

Report signed net outgoing energy, not gross radiated energy or irreversible
loss. Energy can return. This is total local-sector energy, not separately
conserved A energy and not an angular-momentum flux measurement. Time1 and
model units do not establish physical luminosity, galactic fuel sufficiency,
long-lived circulation, or matter/light gravity. No dark matter, expansion or
old-gravity target is introduced. Keep numerical failures if any.
