# SE-Rot2: full Hamiltonian embedding verified at derivative level

All 48 declared full-system directional derivative fixtures pass. Maximum scaled
error is 6.50552077e-11 against the declared 2e-6 gate. All twelve zero-coupling
reductions exactly match the pinned scalar-oscillator field, positional and
internal equations and energy. Positive internal mass and the massive-particle
averaged cone bound also pass on these fixtures.

The new rotor model keeps the six-field Hamiltonian and replaces each scalar
emitter by canonical vectors Q,P. The mass square depends on sampled curl A.
Its field curl source, positional reaction and internal torque are all included;
the angular ledger now includes sum Q cross P. No previous simulation source
was modified. Scalar X/Y fields remain in the model but have no direct emitter
in this variant.

These checks establish the tested implementation's Hamiltonian derivatives,
not a successful evolution or a stable physical theory. In particular, the
derivative-dependent regularized coupling changes field response near matter.
The unchanged matter speed bound does not establish the field propagation
bound. Before production evolution, test the source-region linear response,
stiffness/time-step restrictions, and full angular accounting. A finite rotor
can transfer angular momentum but cannot supply an unlimited reservoir.

Code was pinned at a15748c before execution. Raw controls and hashes are in
rotor-full-v1/results.json. Prior-art attribution and finite source assumptions
remain those in rotor-protocol.md and rotor-report.md. No physical novelty,
galaxy orbit or cluster lensing success is claimed. All twelve goals remain active.
