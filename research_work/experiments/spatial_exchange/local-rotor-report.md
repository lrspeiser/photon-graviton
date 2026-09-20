# SE-LR1: local interaction with multiple propagation speeds

20 September 2026. Protocol ba97a42; executable b925c58. First results are
preserved in local-rotor-v1/results.json with the executable hash and commit.
All288 principal-symbol cases and48 discrete Hamiltonian derivative cases
pass their declared gates. These are interaction-sector consistency tests,
not a galaxy/cluster fit or a completed gravity theory.

## What changed

The finite-kernel rotor used a single internal degree of freedom shared over
a source. SE-Rot4 found an instantaneous response between separated cells.
The present candidate instead assigns Q(x),P(x) at every point, with their
own gradient energy and local curl coupling K=P-epsilon(curl A) cross Q.
There is no source-wide average in this sector. Its equations, assumptions
and fixture selection are in local-rotor-protocol.md.

The new equations have a positive Hamiltonian and a positive principal
symmetrizer diag(I,I,I,c_Q^2 I) for c_Q>0. The same symmetrizer works for
every spatial direction. This supports a symmetric-hyperbolic local
first-order reduction for this sector on smooth bounded backgrounds;
constraints equating gradient variables to actual gradients must also be
retained. It is stronger evidence for local propagation than positive
energy alone. Neither it nor the sampled tests proves global regularity,
long-time nonlinear stability, or a fixed maximum speed for arbitrary Q.
In three dimensions an integrated energy bound alone does not provide a
pointwise bound on Q.

## Results

| Check | Result |
|---|---:|
| Frozen principal cases, including common rotations |288/288 pass|
| Largest numerical/analytic speed discrepancy |1.57e-13|
| Smallest positive characteristic speed in fixtures |0.07883539|
| Largest positive characteristic speed in fixtures |3.30277564|
| Full canonical directional energy derivatives |48/48 pass|
| Largest scaled derivative error |2.55e-9|

Speeds use the uncoupled A-wave speed as unit; they are not measured physical
speeds. Coupling splits branches into faster and slower modes. The unit-speed
uncoupled reference is therefore not a universal limiting cone of this
candidate. The earlier matter/light laws cannot simply be appended and
declared consistent: their coupling and the full joint principal system
still need derivation. No dark matter source, expanding background or
agreement with an older gravity law was used as a test target.

The finite-difference Hamiltonian uses a nearest-neighbor positive gradient
energy and centered curl. Its stencil is local, unlike the earlier shared
rotor. This is not a measured finite-grid signal-front experiment; numerical
dispersion and propagation convergence remain to be tested. No evolution,
matter emission, sustained swirl, stellar orbit or observed lens is generated
by these fixtures.

## Concurrent light-bundle checkpoint

Two of seven shrinking-interpolation bundle runs are complete and pass168
independent archive checks. At n32, the control bend is -0.00997321204 model
radians and the direct-companion bend is -0.00997897351. The magnitude increase
is about0.05777%; the arrival-offset increase is about0.00000792112 model time.
These are forming-field toy results. The matched n48/n64 and smaller-timestep
runs remain pending, so neither spatial accuracy nor the paired enhancement
is established. Both successful archive checks and unfinished campaign status
are recorded in point-bundle-audit.json. These runs use the earlier model,
not the new local rotor interaction.

## Next work and attribution

Derive ordinary-matter generation and reciprocal forces for the local sector
without silently assuming an initial surrounding energy reservoir. Integrate
the scalar/conversion and light sectors and recalculate the joint principal
system. Then test empty-field formation, energy/momentum/angular momentum,
propagation, convergence and perturbations before using this candidate for
galaxy or cluster predictions. The source-energy budget is still mandatory.

Hamiltonian calculus, local vector fields, singular-value analysis and
symmetric-hyperbolic methods are established mathematics, not inventions
claimed here. The specific coupling is a candidate under investigation;
its historical novelty has not been established. Related spin/rotor prior
art is recorded in rotor-report.md. Replacing an interaction sector does
not erase earlier failures or validate the twelve outstanding goal items.
