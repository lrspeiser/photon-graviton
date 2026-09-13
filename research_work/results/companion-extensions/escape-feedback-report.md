# Escaping-wave feedback on protected storage

13 September 2026. Bounded theoretical experiment; not a fit to a galaxy.

## Outcome

Waves waiting to escape can re-excite protected deposits. In these examples they delay protection and allow additional energy to leave through the more fragile state. Escape therefore cannot simply be treated as instantaneous cooling. However, with a permanently stable protected state and continued illumination, the stationary solution still fills every site eventually. This is a finite storage capacity, not unlimited accumulation or a demonstrated physical mechanism.

## Model and provenance

The familiar stimulated-emission and absorption factors are n+1 and n. Our proposed application is a homogeneous population of finite sites, each with empty, bright and protected states of energies 0, 1 and 1-epsilon. This replaces, rather than derives, the earlier unlimited bosonic pair-mode model. It neither derives nor changes the empirical one-third retention reference.

In units where the bright-state energy and relaxation coefficient are one, define

\[
J=(n+1)P_a-nP_b,
\quad \dot P_0=-AP_0+(B+\gamma_a)P_a+\gamma_bP_b,
\]
\[
\dot P_a=AP_0-(B+\gamma_a)P_a-J,
\quad \dot P_b=J-\gamma_bP_b,
\quad \dot n=\mu J-n/\tau.
\]

Here mu is the number of sites per coupled radiation mode and tau is the mean escape time. The last term is an assumed one-zone exponential escape closure; it is not spatial radiative transfer or the exact sphere residence-time distribution. A single occupation represents the coupled modes. Other radiation channels are assumed to escape without returning.

The per-site energy ledger is

\[
U=P_a+(1-\epsilon)P_b,\quad W=\epsilon n/\mu,
\]
\[
\dot E_{in}=AP_0-BP_a,\quad
\dot E_{other}=\gamma_aP_a+(1-\epsilon)\gamma_bP_b,\quad
\dot E_{escape}=\epsilon n/(\mu\tau),
\]
\[
E_{in}=U+W+E_{other}+E_{escape}.
\]

This conservation identity follows from the stated equations. Energy still in transit is included and is not counted as deposited gravity. Momentum, a spatial gravitational response and a microscopic interaction remain outside this population model.

## Numerical experiment

We use epsilon=0.5, A=1, B=1.01 cubed, gamma_a=0.001, and gamma_b either zero or 0.00001. These are illustrative rates, not measured constants. Illumination lasts 100 time units, followed by 1000 with A=B=0. Time units have no established conversion to years. Mu takes 0.001, 1 and 1000; tau takes 0.1, 10 and 1000.

For tau=1000 and gamma_b=0:

| Sites per mode | Protected fraction after illumination | Wave occupation then | Stored energy after dark interval | Wave energy still inside after dark interval |
|---:|---:|---:|---:|---:|
| 0.001 | 0.998160 | 0.000906 | 0.499413 | 0.166635 |
| 1 | 0.586288 | 0.531859 | 0.369274 | 0.091625 |
| 1000 | 0.330629 | 299.369405 | 0.299605 | 0.016793 |

The maximum energy in a fully protected site is 0.5 in these units. Similar stored energies at the end of illumination do not imply equally durable storage: energy left in the bright state is easier to lose. Different cases also absorb different net input energies, which the JSON records separately.

## What is delayed, and what is assumed permanent?

For continued illumination, finite positive tau and gamma_b=0, stationarity gives J=0, then n=0, P_a=0 and P_0=0. Thus P_b=1 is the stationary solution. This algebra identifies the candidate endpoint; finite-duration runs are not a proof of a convergence time for every parameter choice. Once filled, a site accepts no additional excitation in this model.

For gamma_b>0, stationarity instead requires

\[
n=\mu\tau\gamma_bP_b,\qquad
P_a=\frac{P_b(\gamma_b+n)}{1+n},
\]
\[
1=(1+B+\gamma_a)P_a+(1+\gamma_b)P_b
\]

for A=1. A maintained deposit then requires continued input to replace leakage. The numerical calculation independently solves this scalar stationary balance and checks all population and wave equations.

## Verification and next decision

All 18 histories pass probability normalization, endpoint nonnegativity and energy-ledger checks. Eight tighter-tolerance runs agree within 1e-5 in all recorded quantities; actual maximum difference is recorded below. Eighteen stationary roots satisfy the differential-equation residual threshold of 1e-9. The stiff integrator uses an analytic Jacobian.

The next physical task is to derive the available sites, coupled modes, escape geometry and protected-state lifetime from one interaction. Until those are specified, changing mu or tau is a sensitivity study, not an explanation of an observed halo. A dark or weakly coupled outgoing channel remains an alternative, but its reduced reverse coupling and full energy accounting must follow from its proposed interaction. The prior enormous galaxy-scale free-streaming occupation estimates cannot be replaced by these arbitrary small examples.

Executable: `escape-feedback.py`. Full output: `escape-feedback-results.json`.

Maximum refinement difference: 1.26183e-09.
