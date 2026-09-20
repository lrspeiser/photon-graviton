# SE-1S: the emitted wave source is not enhanced in the completed examples

20 September 2026. All eight declared final-state diagnostics pass their
Hamiltonian derivative, energy reconstruction and null checks. Flipping the
mixing sign gives identical integrated wave energy and source. Code was pinned
at c16c07c before execution; source and state hashes are in the raw archive.
The spatial archive has eight completed runs and 121 independent audit checks;
time, space and rotation refinements remain pending.

| Case | X/Y energy including mixing | Integrated source derivative D | D/energy |
|---|---:|---:|---:|
| No emission | 0 | 0 | undefined |
| No excitation | 0 | 0 | undefined |
| Emission only | 0.00136348 | 0.00222359 | 1.63082 |
| Exchange, chi=0 | 0.00133595 | -0.00770774 | -5.76949 |
| Exchange, chi=50 | 0.00134675 | -0.00394664 | -2.93048 |
| Exchange, chi=200 | 0.00136150 | 0.00179982 | 1.32194 |
| Exchange, chi=-50 | 0.00132069 | -0.01182646 | -8.95476 |
| Sign mirror of chi=200 | 0.00136150 | 0.00179982 | 1.32194 |

D is the integral of the X/Y Hamiltonian derivative with respect to U; the
scalar momentum equation receives -g D after spatial integration. It is not
the entire matter-plus-field source, the gravitational force, or a lensing
angle. The signs above therefore do not establish net attraction or repulsion
for the complete system. All values are in the declared model units at time 4.

The largest positive exchange result is about 19.1% smaller than emission-only.
Large local eigenmode response in SE-1M has not translated into enhanced
integrated source in this generated state. The inhomogeneous instantaneous
state need not have SE-1M's homogeneous, phase-averaged mode composition. The
signed mixing derivative can oppose the positive kinetic contribution.
This diagnostic does not yet measure the mode populations, so their exact
role remains to be established.

## Consequence for the next candidate change

Test the emitter's coupling to the two channels, rather than raising chi alone.
The current internal oscillator couples only to X; the high-response branch
at U=0 is Y. A conservative generalization replaces the sampled field Xbar
by cos(theta) Xbar + sin(theta) Ybar in the existing internal-energy square.
A fixed global theta would change the emission channel while keeping a finite
funding reservoir. Both field reactions and the positional force must follow
from that same square. This is a proposed next experiment, not an executed
result or a claim that emission into Y is physically established.

Before any observational comparison, resolve spatial convergence, measure
generated mode content, evolve the resulting field and test stability and
source depletion. A positive source response alone cannot establish a lasting
swirl, stellar orbit or cluster lens. No old-gravity formula supplied a target
or rejection gate here. Hamiltonian differentiation is established math;
the particular constitutive and emitter laws are candidate assumptions.

Reproduce with source_response.py in a new archive directory; the first
source-response-v1 archive must not be overwritten. It includes shell-resolved
signed/absolute contributions, U extrema, internal-energy changes and direct
finite-difference errors for every case.
