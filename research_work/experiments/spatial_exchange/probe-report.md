# SE-P: common matter and light probes

All 72 declared local fixtures and four uniform-field flights pass.

| Check | Maximum error |
|---|---:|
| Hamiltonian derivative | 5.84861048e-11 |
| Cone excess | 2.22044605e-16 |
| Massless speed equality | 2.22044605e-16 |
| Massless momentum scaling | 0 |
| Uniform-flight position | 1.2485332e-14 |
| Uniform-flight momentum | 1.17721175e-16 |

The test-body Hamiltonian is the same regularized particle law as the field
campaign, with m=0 for light. This fixes probe response without adding a lens
multiplier. The cone is the spherical-kernel averaged cone; convergence toward
a point-particle interpretation still needs a probe-radius test. Momentum scaling
checks do not establish electromagnetic polarization or broadband wave behavior.

Uniform backgrounds have constant canonical momentum and constant Hamiltonian
velocity, providing exact flight controls from our own equations. Nonuniform
fixtures check the complete normalization derivative of the spherical kernel.

No bundle has yet been propagated through the evolving exchange field. No
convergence, shear, image position, magnification, time delay or stellar orbit
is inferred from these controls. The planned evolving-bundle experiment needs
its own declared detector geometry and refinement tests before execution.

Raw results and source hashes are preserved in probe-controls-v1; code was
pinned at 6e4ce69 before execution. Hamiltonian ray tracing and the test-body
limit are established methods. The constitutive law is our candidate choice;
no observational success or physical novelty is claimed. All twelve goals remain active.
