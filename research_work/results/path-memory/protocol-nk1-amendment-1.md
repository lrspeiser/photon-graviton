# NK-1 amendment 1: the free-face finish for the galaxy solves; K2 recorded as failed
Declared 20 September 2026 after the first complete archived run, before the rerun. That run is
preserved as `nk1-results-as-declared.json`.

## What failed
The primary solve (reconstructed source) is certified as declared: the two solvers agree to 9e-14 and
the optimality residual is 2.3e-10. The tabulated-force sensitivity is not: its trust-region solve
stopped on step size with all 43 columns active and the two solvers differ by 2.7e-5 with an optimality
residual of 1.9e-2, so gates C1 and C2, which cover both solves, failed. Gate K2 required the shell
columns to change by less than 1e-8 between 128 and 256 quadrature nodes; they change by 4.9e-6, on the
thinnest shells (a tenth of their radius), which 128 nodes resolve only marginally.

## The amendment
Both galaxy solves are finished on their free face by Newton steps with the exact Hessian, as stage 2's
amendment 4 does for the joint solves, before the certificate is read. K2 is not re-thresholded: it
stays as declared and failed, with its magnitude recorded; the shells' entire contribution to the
certified minimum is 0.05 km/s against a margin of 9.4 km/s, so an error of 5e-6 in those columns cannot
touch the reading. Nothing else changes.
