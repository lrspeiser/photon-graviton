# Four-track checkpoint — September 20, 2026 (Pacific)

## Committed work

- `af7e7f05bfa62fb5f7c8949de4319cadc745a4f4`: repaired the quick-mode fit and froze independent quick/full multipair and fixed-Hamiltonian results.
- `21562422d1955766c09be2a1de45c5ce60b2e037`: added the exact finite-link embedding audit, shared matter-source loop, local nonlinear tests, and independent quick/full CI.
- `52b351099d16e8cbea40d4ed2ef627f2e626aec5`: froze independently regenerated integration results, corrected earlier scope claims, and recorded source/environment hashes.

The baseline workflow run `35559531103`, independent quick/full validation run
`35561129077`, and integration-freeze run `35561128863` all passed their
executable checks. This checkpoint does not claim that all physical hypotheses
passed. Numerical histories, including the original full Stage-3I JSON, remain
unchanged.

See [the audit](README.md), `results/*_full.json`, and
`results/review_provenance.json` for the actual evidence and provenance.

## What is closed and what remains open

1. **Validation repair is complete at the implemented scope.** The full four-pair
   result reproduces. Fixed-Hamiltonian state-cutoff errors are now separated
   from comparisons among rebuilt cutoff models. Quick checks are not presented
   as full checks, and skipped checks are not reported as passing.
2. **The proposed natural chain embedding is rejected.** Exact leakage is nonzero
   on L=2,3,4. The chains have been relabeled reduced models. Full-real-space
   Gaussian 3D refinement is an explicitly conditional control; it does not
   establish a finite-spin 3D dynamical phase. A controlled dressed embedding or
   finite-spin many-body volume scaling is still required.
3. **The shared free-matter source loop is calculated, not the physical cone.**
   All six slab species, mirrors and doublers contribute to both responses with
   one inherited normalization. The gauge Ward test passes, while the existing
   frame source has a nonzero uniform response on both linear and
   determinant-preserving source paths. No source-specific counterterm or
   compensating constant was introduced. Interacting photon/frame pole matching
   and asymptotic infrared convergence remain open.
4. **Initial nonlinear terms and failures are explicit.** The completed-move
   commutator forces a composite energy current. A conditional affine-density
   bootstrap fixes the cubic kinetic coefficients. Local volume conditions
   change the first-/second-class split without adding a linear physical mode.
   Linear Dirac evolution preserves constraints without repeated projection.
   The naive local central-difference nonlinear brackets do not close exactly;
   microscopic nonlinear gravitational closure is not established.

## Frequency and normalization conventions

The temporal response coefficient in `check_shared_matter_ir.py` is the
**Euclidean-frequency** derivative at zero. For an occupied-to-empty gap d the
response contributes `-w*d/(d^2+omega_E^2)` and hence `+w/d^3` to the derivative
with respect to `omega_E^2`. A retarded pole requires analytic continuation and
matched bare kinetic terms. The reported temporal/spatial coefficient ratios
are therefore not physical photon or graviton speed predictions.

The uniform source curvatures are raw filled-band energy responses in inherited
lattice units, not graviton masses or established physical instabilities. The
finite matter regulator's geometry dependence and missing frame/connection
couplings must be derived before a physical interpretation. Differencing
Pi(k)-Pi(0) measures a derivative; Pi(0) is retained, not deleted from the model.

The cubic bootstrap assumes that the frame metric and its conjugate momentum
transform as a metric and a weight-one tensor density. That assumption must be
matched to the finite junction algebra; the bootstrap is not itself that match.

## Next decisive microscopic gate

Derive the finite frame/area interaction and the complete matter-plus-regulator
stress coupling from the same completed-move algebra. Keep every generated
operator, including the composite currents. Recompute the local constraints and
stress Ward identities, then match photon and frame kinetic normalizations and
perform controlled low-momentum/volume extrapolations. Do not repair the source
response by independently choosing a photon or gravity coefficient. Preserve
failed constructions as diagnostics rather than relabeling them successful.
