# SE-Rot4: finite-width rotor coupling is nonlocal inside matter

All six checks reproduce the factorized response and reciprocity. Nonzero
remote coupling is detected for both nonzero strengths on both grids.

| Grid | Coupling | Separated-cell stiffness response | Response / cell volume |
|---:|---:|---:|---:|
| 24 | 0.0 | -0 | -0 |
| 24 | 2.0 | -0.000157329319 | -0.00424789161 |
| 24 | 10.0 | -0.00393323297 | -0.10619729 |
| 48 | 0.0 | -0 | -0 |
| 48 | 2.0 | -2.48099719e-05 | -0.00535895393 |
| 48 | 10.0 | -0.000620249297 | -0.133973848 |

The cells are separated by4/3 model units. The original local Laplacian has
exactly zero direct matrix element between them; the finite-width sampled rotor
creates a nonzero instantaneous coupling. Dividing by cell volume shows that
the shrinking nodal response is not evidence that the finite-width kernel
becomes local: its continuum integral-kernel interpretation remains nonlocal.

This is a limitation of the claimed scope, not a hidden numerical failure.
It does not measure a signal speed or reject every effective finite-sized
source description. It does prevent citing the vacuum field cone or positive
stiffness alone as proof of strict pointwise finite-front propagation through
matter. Other shared finite-size source averages also require care; this test
isolates the new rotor contribution and is not a general locality certificate
for the earlier source regularization.

A local completion would distribute internal degrees of freedom through matter
and specify their propagation and reaction, rather than use one instantly
shared rotor over its entire support. This changes the model and must be derived
and tested separately. Full rotor simulations under the present approximation
would be effective nonlocal-source tests, not proof of a causal gravity theory.

The protocol initially had a written factorization sign error, corrected and
explained before execution. The implementation and the checked factorization
agree without changing the existing rotor model. Raw checks, source hashes and
commit pin are preserved in rotor-locality-v1. No old-gravity target was used.
All twelve original goals remain active.
