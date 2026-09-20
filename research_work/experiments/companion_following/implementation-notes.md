# GF-1 implementation details declared before screening

The radius perturbation is additive: r=3+0.02 cos(2 theta).
Ring behavior gates use the worst sampled circulation, radius bounds and
radial spread over the entire run, not only a favorable final frame.
Circulation is the absolute mean normalized signed tangential velocity about
the instantaneous center of mass, using velocity relative to that center.
Closest approach is checked every integration step. Other diagnostics and all
nine state components are stored every 0.2 time units.

Chain gain is each follower's maximum absolute principal velocity angle divided
by the input amplitude; the last four gains are averaged. Correlation compares
the mean signed angle of the last four followers with the input pulse, using
nonnegative lags up to 20 time units and overlapping samples. Arrivals use the
declared absolute 0.02 angle threshold. A correlation is descriptive and does
not imply finite propagation speed in these instantaneous models.

For the conservative chain, eliminate the prescribed leader acceleration from
the follower mass-matrix equations; record the required external generalized
force on the leader. Never simply overwrite one acceleration after solving the
closed system. Stored chain external work is sampled trapezoidal quadrature,
not an independent energy-conservation assertion.

M3 stores z2dot in the z1 state slots, initialized to zero; z2 initially equals
velocity. All other memory state slots begin at velocity. Inactive batch rows
after nonfinite failure are masked in evidence and cannot pass; isolated dummy
rows prevent floating-point overflow from disrupting unrelated simulations.
No failed trajectory is resumed or presented as repaired.

A model exception stops its invocation and preserves earlier output. An
unavailable positive circular balance is an explicit failure requiring a report,
not an invitation to change the fixture. No such outcome has yet been evaluated.

The minimum-step comparison uses gain and final mean-radius ratio, with all
original behavioral and conservative gates also retained. Negative coupling is
a new follow-up, not retroactive replacement of a failing positive-coupling case.
Conservative and phenomenological selections are separate. Additional
zero-input and reversed-coupling ring controls do not change selection.
