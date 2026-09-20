# SE-LR2: finite-excitation local circulation formation

Declared 20 September 2026, before execution. Use the unchanged SE-LR1 local
Hamiltonian and RHS. This tests actual emission and reciprocal energy transfer,
not merely its principal symbol. It remains an internal-excitation subsystem:
ordinary-matter generation of Q/P and matter/light trajectories are absent.

Periodic box L8; initially A=Pi=0. Define f=max(1-r^2/1.2^2,0)^3,
Q=.03 f e_x and P=.03 f e_y. For the nonrotating control set P=.03 f e_x.
All initial energy is counted, including Q's gradient energy. There is no
external forcing or surrounding field reservoir. c_Q=.5, Omega1, omega_A=.2.
Run RK4 to T1. No claim of a sustained or astrophysical source is permitted.

Six cases: n24 dt.02 epsilon0,+2,-2 with rotating internal excitation;
n24 dt.02 epsilon2 nonrotating; n24 dt.01 epsilon2 rotating (time control);
n32 dt.01 epsilon2 rotating (space control). Preserve every-step energy,
canonical momentum/angular momentum and sector energy, initial/final fields.

J=integral[-x cross sum_a(Pi_a grad A_a+P_a grad Q_a)+A cross Pi+Q cross P].
Use centered derivatives for continuum momentum diagnostics. Grid rotational
symmetry is imperfect, so measure rather than presume exact J conservation.
Split energy into A terms and K/Q terms; the K/Q sector includes interaction.
Measure A circulation at r1,1.5,2 in the xy plane using trilinear interpolation
and128/256 loop points. Circulation is not a lensing prediction or force.

Individual gates: finite states, maximum relative total-energy drift<1e-5,
maximum angular-vector drift / max(initial norm,1e-12)<.01 for rotating cases,
and loop quadrature absolute disagreement<1e-5. Report nonrotating angular
drift as absolute. Report energy in |x| or |y| or |z|>=3 divided by initial
energy; require<1e-5 before interpreting formation as free of boundary effects.
No strict finite-grid front is claimed.

Comparison gates: epsilon0 and nonrotating A/Pi max amplitude<1e-12;
epsilon sign reversal gives opposite A/Pi and identical Q/P to1e-12;
time-control final r1.5 circulation relative difference<.001;
spatial relative difference<.05. Relative denominator max(abs(finer),1e-10).
Report any failure without relaxing gates. An unresolved spatial result must
not be promoted as converged emission. No nonzero circulation threshold is
required to call numerical bookkeeping successful; report its actual value.

Positive local coupling and finite excitation do not establish physical fuel
availability, identified gravitons, global stability, old-gravity agreement,
or observational success. Established Hamiltonian and finite-difference/RK4
methods are credited as mathematical tools, not original gravity theories.
