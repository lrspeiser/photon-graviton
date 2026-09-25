# Provenance

This folder is the package "Main audit and the energy-shift direction" supplied by the user on 25 September 2026,
after the step-back audit (`STEP-BACK-AUDIT.md`, commit 67075c2). `README.md`, `checks.py`, `results.json`,
`baseline_results.json`, `reproduction.log`, `environment.json` and `MANIFEST.sha256` are as supplied, produced outside
this repository with Python 3.13.5, NumPy 2.3.5 and SciPy 1.17.0. `sha256sum -c MANIFEST.sha256` passes for all six
files.

Checked here on 25 September 2026 (Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1): `python checks.py` in a copy of the
folder reproduces every number in `results.json` to better than 1 part in 10⁹ (the eight-piece finite-occupation
table, the force-versus-energy-derivative check, the local-energy scalings and the positive-mixture exponent bound).

What it contains: an exact diagonalisation of an effective exchange model (eight two-level pieces, exchange kernel
J = −C e^{−κr}/r, every excitation number 0–8) showing attraction through pair correlations with zero one-point
amplitude, a force that vanishes when every piece is excited, and two algebraic limits (a receiver whose energy
depends only on the local intensity cannot give both exponents; a fixed-weight positive mixture of Yukawa forces cannot
fall more slowly than 1/r²). The kernel's constants come from the user's earlier "structured-medium derivation"
package, which was not supplied to this repository. The review's wording corrections to the audit and its research
protocol are taken up in README §31.
