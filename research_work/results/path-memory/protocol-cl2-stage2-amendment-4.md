# CL-2 stage 2, amendment 4: finishing the joint solves on their free face
Declared 19 September 2026 after the rerun under amendment 3, before the next rerun. That run is
preserved as `cl2s2-results-amendment-3.json`.

## What happened
Under amendment 3 the three-block joint solve was certified (projected gradient 1.5e-7 against an
objective of 1426, one part in 1e10), but the transfer test's galaxies-plus-clusters joint solve,
started from the three-block solution's boundary pressures, stopped at 3.2e-2 against 998, three parts
in 1e5, above the adopted 1e-6: an interior-point polish terminates on its step size and its final
gradient depends on the start. Every number outside the joint solves is identical to the first run.

## The amendment
After the polish, each joint solution is finished by Newton steps on its identified free face with the
exact Hessian (established: an active-set Newton finish), freeing any bound variable whose gradient
points inward and binding any free variable the step drives to zero; the certificate is read on the
finished point. In development, from five starting points including the origin and randomised ones,
the finished projected gradient was below 1e-14 of the objective and the objective agreed to seven
decimals. The objective, weights, blocks and the 1e-6 gate are unchanged. The rerun writes
`cl2s2-results.json`.
