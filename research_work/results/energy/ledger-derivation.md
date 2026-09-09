# T20 initial calculation: explicit energy reservoirs

This implementation covers linear, fixed-volume ledgers with constant rates, finite fuel and a finite source interval. It does not yet derive an interaction, momentum transfer, pressure work in a changing spacetime, or an entropy law.

Let P, C, D, X and F be photon, companion, deposit, receiving-sector and stellar-fuel energies. Define nonnegative h, Gamma, lambda_c and lambda_d. While fuel remains, a source injects j until the declared source stop time:

P-dot = j - h P
C-dot = h P - (Gamma + lambda_c) C
D-dot = Gamma C - lambda_d D
X-dot = lambda_c C + lambda_d D
F-dot = -j

Adding the five equations gives zero. If lambda_c=lambda_d=0, all energy can reach permanent deposits after finite fuel is exhausted and capture is nonzero. If either loss term is positive, the missing reservoir energy moves to X; the total still conserves. The choice of what X physically is remains an open action-level question.

For one initial photon-energy unit, no subsequent source and zero losses, P=exp(-h t), C=h[exp(-h t)-exp(-Gamma t)]/(Gamma-h), D=1-P-C. At Gamma=h, C=h t exp(-h t). At h=0 conversion vanishes. At Gamma=0 deposits remain zero. These cases are checked against a matrix-exponential solution.

For a continuing source with no fuel limit, P tends to j/h, C to j/Gamma and D grows indefinitely if there is capture but no deposit loss. Finite fuel changes that conclusion: the stored energy can approach the finite initial fuel plus initial radiation. Therefore permanent deposits do not require infinite energy growth when emission eventually stops. They still require a physical trapping mechanism and an acceptable formation history.

The solver also independently integrates the finite-fuel cases with adaptive ODEs and compares them to matrix exponentials. Ten checks passed, including equal rates, zero capture, zero production, shared secondary loss, finite residence and fuel exhaustion. See ledger-checks.json for errors. These are dimensionless demonstrations, not fitted astrophysical rates.

Run: python ledger_solver.py --output ledger-checks.json

Next: supply a candidate microscopic receiving/trapping sector, derive momentum/pressure/entropy terms, and use a user-selected history or response law. T20 remains in progress until those scope limitations are incorporated into the branch comparison.
