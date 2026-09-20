# SE-JT1: reciprocal joint Hamiltonian implemented and checked

20 September 2026. Protocolb19d342; implementation82ef464. All48 full-state
Hamiltonian derivative fixtures and12 zero-feedback reductions pass. Largest
scaled derivative error is5.82e-9; largest reduction error is2.23e-16. First
results and source hashes are preserved in joint-transfer-v1/results.json.

The implementation now combines seven coordinate fields (scalar phi,
vector A, internal vector Q), their canonical momenta, and massive/massless
particles in one energy function. The drift is bounded by c_Q C and the
internal kinetic momentum is K=P-lambda A-epsilon(curl A) cross Q. It uses
the shared-drift candidate declared in ST1; it is not the earlier LR3 evolution
with a new label.

Particle forces, velocities and reciprocal scalar/vector sources come from
the same interpolated particle Hamiltonian. Field equations include derivatives
of C(phi), beta(phi,A), the direct/curl momentum shift, and all centered-current
terms. Particle mass is fixed; their Hamiltonian energy can exchange with the
fields. The massless particle has nonzero momentum and backreacts, rather than
serving as an externally prescribed passive ray.

The48 checks perturb the full field-plus-particle state and compare two
central-difference energy derivatives against its canonical RHS pairing.
Field terms carry cell volume; particle terms do not. The12 reductions turn
off g and eta and reproduce the prior local transfer equations plus a free
scalar and freely moving particles. Positive energies and interpolation
velocity bounds pass. These are implementation-consistency controls, not
an independent proof of long-time stability or observational success.

## What this enables next

With all fields initially zero, moving particles source A through the
derivative of beta dot p; the direct coupling can then drive Q/P from A.
Thus a finite internally excited Q/P reservoir need not be inserted as an
initial condition. This is a route permitted by the implemented equations,
not yet a demonstrated production run or a sufficient physical source budget.
The next declared evolution must compare empty-field formation with the
uncoupled case, track particle energy loss and every field sector together,
and calculate the actual massive and massless trajectories.

Full nonlinear propagation/stability, extraction of the implemented joint
principal system, point-interpolation convergence, rotation/source-size tests,
outgoing boundaries, stable outer orbits and lens bundles remain required.
The point-particle interpolation has a shrinking lattice stencil; it is not
an exact finite-grid isotropic or local continuum source. The photon ansatz
does not provide electromagnetic polarization or attachment/detachment physics.

No dark matter or expanding background enters the implementation. Established
Hamiltonian calculus, the dispersion ansatz, finite differences and positive
interpolation are mathematical ingredients; no historical novelty of those
structures is claimed. Real galaxy/cluster observations have not been fitted
by this new joint candidate, and the twelve-goal solution remains incomplete.
