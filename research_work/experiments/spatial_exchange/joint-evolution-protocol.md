# SE-JE1: empty-field generation and coupled trajectories

Declared 20 September 2026 before evolution. Use pinned JT1 equations, no
prescribed internal excitation, external source or damping. Initially all
phi,A,Q and canonical field momenta are zero. Six particles of bare mass1
form a radius.7 ring in xy; tangential momentum magnitude is.2 for the fast
fixtures and.0007 for the slow pair. A backreacting massless particle starts
at(-1,1,0) with momentum(1e-4,0,0). Reference units have vacuum speed1;
neither spatial nor physical fuel/time calibration has been established.

Periodic L10,T2. The photon is an endpoint-deflection diagnostic, not a
complete lens bundle or detector time-delay calculation. Use source-frame
velocity angle atan2(v dot rotated_y,v dot rotated_x) and record the full
3D direction too. Massive trajectories are short responses, not stable orbits.

Eight ordered runs:
1. free: n24dt.02,epsilon2,lambda1,g=eta=0,fast;
2. baseline: n24dt.02,epsilon=lambda=0,fast;
3. combined: n24dt.02,epsilon2,lambda1,fast;
4. slow-baseline: n24dt.02,epsilon=lambda=0,slow;
5. slow-combined: n24dt.02,epsilon2,lambda1,slow;
6. time: combined n24dt.01;
7. space: combined n32dt.01;
8. rotation: combined n24dt.02, angle.573 about normalized(1,2,3), with
   both sources and photon position/momentum rotated.

Archive initial/final full states, every-step field/matter/total energy,
canonical total momentum/angular momentum, K/Q-sector energy and A/Q angular
momentum, particle states, endpoint photon velocity, cone residuals and edge
amplitude (any coordinate magnitude>=4). K/Q energy includes the interaction
through K and must not be called a separately conserved radiation energy.

Individual gates: finite states, relative total-energy drift<1e-5, relative
total-angular-vector drift<1% with denominator max(initial norm,1e-8), cone
error<1e-10,edge amplitude<1e-5. Free control must retain zero fields to1e-12
and reproduce straight particle flight to1e-10. Time bend difference<.1%;
space bend difference<5%; rotation-back photon velocity direction error
<2% of the unrotated bend. Relative bend denominator max(abs(reference),1e-8).
Record total momentum residual; no unproved exact finite-grid momentum
symmetry is assumed. Report failure without changing thresholds.

Compare on/off K/Q generation, energy loss of particles and photon bend for
both fast and slow pairs. No desired sign or minimum enhancement is a
numerical pass criterion. Do not call a small effect resolved before numerical
errors are controlled. Six massive source particles also feel their own
regularized fields; self-force/point-limit implications remain to be tested.

This moves toward source generation and joint trajectories, but does not
establish galactic longevity, realistic source profiles, global stability,
Maxwell/polarization physics or an observational fit. No dark matter,
expansion or old-gravity prediction is used as a gate.
