# CL-2 stage 2, amendment 3: the joint solves' certificate
Declared 19 September 2026 after the first complete archived run, before the rerun. The run is preserved
as `cl2s2-results-as-declared.json`.

## What happened
The protocol's S2-7 joint solves combine the galaxy block's exact speed loss, in (km/s)^2, with the
cluster and lens blocks' chi-squares per point, with stage 1's per-block 1/N weights. In the run the
galaxy term (about 1000) outweighs the other blocks (about 15 for the clusters and up to a few hundred
for the lenses) and the bounded L-BFGS solve stalled at the boundary of the loss's domain: its
projected gradient was 150 against an objective of 1426, a tenth, so the joint amplitudes are not
certified. The code's joint-certificate check (projected gradient below 1e-6 of the objective), which
the protocol did not list among its gates, failed as recorded. The galaxies-alone certificate (G-C1 to
G-C3), the cluster comparison and every other gate are unaffected, and the decision rule does not use
the joint solve.

## The amendment
Each joint solve is polished with the trust-region Newton method of S2-4 (solver B, analytic Hessian)
from the bounded L-BFGS point, in variables scaled by the square root of the Hessian's diagonal at that
point (Jacobi scaling, established). An interior-point polish leaves bound variables at tiny positive
values, so amplitudes and nuisances below 1e-9 of the largest are set to zero before the optimality
reading, and only a negative gradient counts against a variable at zero (the Karush-Kuhn-Tucker
reading). The joint-certificate check is adopted as a gate: projected gradient, so read, below 1e-6 of
the objective. In development the galaxies-plus-clusters solve reached 1.8e-11. The objective, weights and blocks are unchanged;
the unit convention that makes the galaxy term dominant is recorded in the report as a property of the
declared combination, not changed. The rerun writes `cl2s2-results.json`; the first run's every
other number is deterministic and reproduces.
