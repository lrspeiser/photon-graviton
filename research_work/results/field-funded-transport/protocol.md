# Field-funded radiation exchange: declared test

This is a new constitutive branch, not a parameter refit. Use reference energies R (shared EM/GW radiation), T (traveling companions), M (traveling temporal energy) and D (deposits); n=1+M/K, v=1/n, clocks d tau_clock=dt/n. Define b=R/(Kn), d=1-b, Q=(T-M)/tau. Assign radiation exchange to M:

R_t+(vR)_x=-b M_t;
T_t+(vT)_x=-Q-Gamma*T;
M_t+(sigma*v*M)_x=b M_t+Q;
D_t=Gamma*T.

The sum is conservative with total flux v(R+T)+sigma*v*M. Elimination gives M_t=[-div(sigma*v*M)+Q]/d. This fixes the previous empty-T gradient term but restores a possible time-matrix singularity. The local M characteristic is sigma/(n*d) times local c0. Require d>0 and n*d>sigma for the declared numerical branch; document that this is a checked trajectory domain, not a proven universal invariant. The temporal sector is explicitly slower (sigma=0.5), while EM, GW and T share local c0. No claim of ordinary graviton dynamics.

Before calculation: nondimensional K=2,tau=1,sigma=0.5, x=0..10,t=0..24. Start R=0.22, M=D=0, T=0.1*exp[-((x-2)/0.4)^2]; include zero-seed control. Radiation inflow 0.22*[1-tanh((t-2)/0.1)]/2, no incoming T/M. Compare no capture with Gamma=1/[1+exp(-(x-8)/0.25)]. Use 80,160,320 cells and limited second-order face reconstruction. Retain all failed refinement checks; require final adjacent-grid absolute measured-z difference below 0.003, energy-sum error below 1e-8, no energy below -1e-9, and clock/event identity error below 2e-4.

Probe (start,end,emission time)=(0,10,1),(2,9,0.1),(0,10,5). Derive ray arrival and J using dt/dx=n, d log J/dx=n_t. Compare S=J*n_emit/n_observe against independently integrated proper-clock intervals using emission times te +/-0.001. EM/GW equality is imposed, not observationally established. No astronomical calibration or holdout is used. Check local empty-T derivatives against the prior counterexample and record the characteristic/positivity domain along every integration step.

Follow-up after the 160-to-320 redshift refinement failed: add 640 cells for both seeded cases without changing parameters or the 0.003 gate. At this refinement also integrate radiation-out flux and signed radiation-to-M work W=integral b*M_t dx dt as separate ledgers. Verify R_final+R_out+W=R_initial+R_in to 1e-8. W may be negative; do not presume a deposit came from radiation rather than the initial seed. These new ledgers do not alter the evolution equations. Keep earlier grid results.
