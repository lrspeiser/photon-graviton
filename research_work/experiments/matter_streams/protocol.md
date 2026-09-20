# MS-1: matter-fed moving streams and cumulative ray steering

Declared before this extension's calculations. The owner explicitly adds
emission from ordinary matter, outgoing companions following companions, and
light carried along a curving stream rather than merely crossing it.
This is a separate kinematic extension of GF-1, not a completed gravitational
theory. No dark matter, expanding geometry, fitted distance or observed fit.

## The distinction being tested
A coherent directional interaction can accumulate a bend along a curving
stream. It need not increase the local field energy on every encounter.
If encounters amplify the field itself, the required gain must be funded by
a source; that nonlinear gain is not silently included in the present test.
Standard lensing already integrates deflection along a path; its thin-lens
description is not a claim that gravity acts at only one point.
Credit: Bartelmann and Schneider, Weak Gravitational Lensing, sections 3-4,
https://arxiv.org/pdf/astro-ph/9912508 .
Direction alignment and collective turning are credited in GF-1's provenance;
the projector law is a mathematical ansatz, not a novelty claim.

## Source and prescribed transport
Dimensionless c_light=1. A compact ordinary-matter source has initial available
energy 10 and emits Q=1/30 per time unit from launch radius a=0.5. At photon
injection it has emitted for age0=30. The companion speed is cg=0.5.
For a pitch h=tan(beta), radial speed vr=cg/sqrt(1+h^2); signed azimuthal
speed is sign*h*vr. The equatorial stream tangent is
t_g=(e_r+sign*h*e_phi)/sqrt(1+h^2).
The spherical-average energy transport solution is
u_g(r,t)=Q/(4*pi*vr*r^2), a<r<a+vr*(age0+t), and zero elsewhere.
Its total energy is Q*(age0+t), and source energy is 10 minus that value.
Source depletion time is never reached in the time80 ray runs.
All pitch choices have the same emitted energy at the same time:
slower radial travel raises local density but reduces the outer reach.

This is an equatorial ray slice with an analytically propagated radial energy
front. The angular stream geometry is prescribed. Neighbor-to-neighbor
formation of this wind, three-dimensional axial regularity, emitter torque
and a closed momentum-carrying mediator are not derived. The central source
is a representative source test, not a claim to have simulated every type of
matter or a galaxy's actual source history.

## Ray equation and variations
x_dot=n=(cos(theta),sin(theta)),
theta_dot=kappa*u_g*F(cos(psi-theta))*sin(psi-theta).
This is dn/ds=kappa*u_g*F(n dot t_g)*(I-n n^T)*t_g.
It preserves photon speed and energy by construction. Ray momentum transfer
is recorded as the opposite impulse owed to the prescribed stream; that
ledger is not a completed backreaction dynamics.
The law is not yet the acceleration law for stars.

Pitches h=0,0.5,1,2,4; both swirl signs.
Three gates: F0=1; F1=((1+c)/2)^2; F2=max(c,0)^2.
Couplings kappa=-16,-4,-1,0,1,4,16.
These are 210 parameter cases, including duplicate zero-pitch sign controls,
not 210 unique theories. Negative signs are explicit exploratory controls.
For each: rays start at x=-20, impacts b=+-0.5,+-1,+-2,+-4,+-8, theta=0.
Trace to their first crossing of observer plane x=20 or time80. A nonarrival
by time80 is reported as such, not as proof of permanent capture.
The front evolves during transit; no frozen-age hidden energy is added.
Use RK4 dt0.02, record every0.2, save all trajectories. Interpolate the first
observer-plane crossing within the last integration step.

Record direction change, accumulated absolute turning, exposure integral
int u_g F ds, travel time, observer position and nonarrivals.
For each +/- impact pair, separate inward bending
A=(theta_minus-theta_plus)/2 from sideways bias
B=(theta_minus+theta_plus)/2, using principal exit angles.
A sideways shift is not automatically gravitational focusing.
No fit to observed Einstein radii, shear, stellar dynamics or cluster pressure.

## Controls and refinement
Before scoring: zero coupling, exact straight-stream relaxation
tan(theta/2)=tan(theta0/2)*exp(-kappa*u*s) for F0; zero guide work;
source-energy quadrature and total budget; the unchanged total source energy
for all pitches; and straight aligned light has zero bending.
Compare complete chirality/impact reflection on the grid.
Rank configurations by mean inward pair deflection at |b|=1,2,4 among those
with all six observer arrivals; break ties by ID. This is exposed selection,
not a holdout. Refine the top6 at dt0.01. Compare exit angle and exit position
within0.01 absolute, and exposure within1% of max(1, magnitude).
Keep failures; do not promote unconverged large bends.
Store exact source hashes, controls, scores and all signed alternatives.

The intended result is a conditional test of directional cumulative steering.
A visible enhancement under this ansatz does not establish increased gravity,
a derived graviton interaction, a self-formed swirl, or solved cluster lensing.
