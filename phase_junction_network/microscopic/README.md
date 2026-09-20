# Microscopic finite-state phase

**Status:** finite kinematics, exact first-order algebra, and independent continuum-scaling checks are complete; a stable finite quantum gravity realization remains open  
**Date:** 2026-09-20  
**Branch:** `agent/phase-junction-microscopic`

This folder replaces the infrared rotor and symmetric-tensor fields with explicit finite local Hilbert-space candidates and tests which microscopic dynamics can or cannot reach the required continuum theory.

## Results obtained

1. **Electromagnetism:** a spin-1 quantum link, with three states per oriented link, is the smallest tested representation that has exact local U(1) Gauss symmetry and a non-constant electric-flux energy. Spin-1/2 preserves Gauss symmetry but its `E^2` term is constant.
2. **Gravity kinematics:** six odd-prime qudits per lattice site realize the three vector constraints and one scalar curvature constraint as an exact commuting Weyl/CSS algebra. The periodic code leaves `2*N + 4` logical qudits: two local modes per site and four global zero modes.
3. **Finite-dimensional obstruction:** no finite matrices obey `[Q,P]=iI` on their full Hilbert space. Continuous frame variables must emerge below a cutoff or from a discrete microscopic algebra.
4. **Pure-stabilizer dynamics rejected:** exact local coordinate invariants begin at two derivatives and momentum invariants at one. Squaring them gives `omega ~ k^3`.
5. **First-order identity verified:** an independent torsion-free connection gives a local `C^2 + C*partial(h)` action whose Schur complement is the Fierz–Pauli stiffness.
6. **Exact arithmetic cross-check:** the Schur/Fierz–Pauli equality is verified over the rationals by matching all six independent quadratic momentum coefficients, not only random floating-point samples.
7. **Continuum scaling cross-check:** the lowest tensor gap fits `omega ~ L^-0.995891`, the two polarizations remain degenerate to `6.11e-16`, and cubic directional anisotropy falls as `L^-2.01111`.
8. **Naïve auxiliary Hamiltonian rejected:** a positive real auxiliary block with no bare frame stiffness always induces a nonpositive Schur complement. The current 18-component connection must be constrained/nonpropagating, integrated with a nontrivial contour, or accompanied by explicit bare frame stiffness.

## Current interpretation

The finite gravity constraint code is kinematic. A commuting-projector or manifest-local-invariant-square Hamiltonian is not the required graviton theory.

The first-order frame–connection identity remains a viable **constrained action architecture**, but it is not a conventional positive oscillator Hamiltonian for all connection components. The connection quadratic form has inertia

```text
8 negative, 0 zero, 10 positive
```

and the full unconstrained first-order block has

```text
9 negative, 3 gauge zero, 12 positive.
```

After eliminating the connection and imposing the scalar constraint, the physical transverse-traceless sector has two positive modes. The next quantum construction must implement that constrained reduction rather than treating every connection component as an independent positive-energy particle.

Three branches are now separated:

1. a constrained first-order/Palatini transfer matrix;
2. a complex Euclidean auxiliary representation followed by reflection-positivity tests;
3. a bounded real Hamiltonian with explicit bare frame stiffness, used as a benchmark or fallback.

The project does not yet derive `Z_g/Z_A`, chiral fermions, the particle-mass hierarchy, nonlinear closure, or control of the vacuum-volume term.

## Files

| File | Purpose |
|---|---|
| [`finite_hilbert_derivation.md`](finite_hilbert_derivation.md) | Finite construction, canonical no-go result, discrete constraints, and derivative-order obstruction. |
| [`finite_junction_model.py`](finite_junction_model.py) | Quantum-link test, GF(p) gravity constraint code, and finite canonical-pair diagnostics. |
| [`results.json`](results.json) | Frozen finite-link and finite-constraint results. |
| [`search_local_dynamics.py`](search_local_dynamics.py) | Polynomial search for exact local gauge invariants. |
| [`dynamics_results.json`](dynamics_results.json) | Frozen `omega ~ k^3` obstruction. |
| [`first_order_frame_connection.md`](first_order_frame_connection.md) | Local independent-connection derivation. |
| [`check_first_order_frame.py`](check_first_order_frame.py) | Floating-point Schur, gauge, and TT checks. |
| [`first_order_results.json`](first_order_results.json) | Frozen random-momentum results. |
| [`check_first_order_exact.py`](check_first_order_exact.py) | Exact rational coefficient and gauge-null proof. |
| [`first_order_exact_results.json`](first_order_exact_results.json) | Frozen exact-arithmetic output. |
| [`auxiliary_stability_no_go.md`](auxiliary_stability_no_go.md) | Proof that the naïve positive auxiliary Hamiltonian cannot generate the desired stiffness. |
| [`check_auxiliary_stability.py`](check_auxiliary_stability.py) | Numerical inertia, Schur-sign, Euclidean, and stable-completion checks. |
| [`auxiliary_stability_results.json`](auxiliary_stability_results.json) | Frozen auxiliary stability output. |
| [`check_tensor_scaling.py`](check_tensor_scaling.py) | Multi-size tensor gap, polarization, and anisotropy checks. |
| [`tensor_scaling_results.json`](tensor_scaling_results.json) | Frozen finite-size output. |
| [`verification_routes.md`](verification_routes.md) | Explains the independent ways the theory must be verified. |

## Reproduce

From the repository root:

```sh
python phase_junction_network/microscopic/finite_junction_model.py \
  --output phase_junction_network/microscopic/results.json

python phase_junction_network/microscopic/search_local_dynamics.py \
  --output phase_junction_network/microscopic/dynamics_results.json

python phase_junction_network/microscopic/check_first_order_frame.py \
  --samples 500 \
  --seed 43 \
  --output phase_junction_network/microscopic/first_order_results.json

python phase_junction_network/microscopic/check_first_order_exact.py \
  --output phase_junction_network/microscopic/first_order_exact_results.json

python phase_junction_network/microscopic/check_auxiliary_stability.py \
  --samples 500 \
  --seed 44 \
  --output phase_junction_network/microscopic/auxiliary_stability_results.json

python phase_junction_network/microscopic/check_tensor_scaling.py \
  --output phase_junction_network/microscopic/tensor_scaling_results.json
```

Every script exits nonzero if a committed claim fails.

## Immediate next implementation

Perform the complete constrained canonical reduction of the first-order frame–connection branch. Construct the reduced symplectic form and physical Hamiltonian before truncating the remaining degrees of freedom. Then build a finite transfer matrix on small periodic lattices and test:

- Hermiticity or reflection positivity;
- exactly two tensor branches;
- `1/L` gap scaling;
- no propagating connection or scalar state;
- finite-state boundary leakage;
- shared microscopic coefficients with the electromagnetic quantum links.
