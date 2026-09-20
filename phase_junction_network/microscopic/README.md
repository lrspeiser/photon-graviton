# Microscopic finite-state phase

**Status:** finite kinematics and a local first-order linear-dynamics architecture are now checked; finite quantum dynamics remains open  
**Date:** 2026-09-20  
**Branch:** `agent/phase-junction-microscopic`

This folder replaces the infrared rotor and symmetric-tensor fields with explicit finite local Hilbert-space candidates and tests which microscopic dynamics can or cannot reach the required continuum theory.

## Results obtained

1. **Electromagnetism:** a spin-1 quantum link, with three states per oriented link, is the smallest tested representation that has exact local U(1) Gauss symmetry and a non-constant electric-flux energy. Spin-1/2 preserves Gauss symmetry but its `E^2` term is constant, so it cannot supply the required electric stiffness by itself.
2. **Gravity kinematics:** six odd-prime qudits per lattice site realize the three vector constraints and one scalar curvature constraint as an exact commuting Weyl/CSS algebra. On periodic odd lattices, away from modular degeneracies, the code leaves `2*N + 4` logical qudits: two local modes per site and four global zero modes.
3. **Finite-dimensional obstruction:** no finite local matrices obey an exact continuous canonical commutator `[Q,P]=iI` over their full Hilbert space. Continuous linearized diffeomorphism symmetry must therefore emerge in a low-energy sector or be replaced microscopically by an exact discrete constraint algebra.
4. **Low-energy escape:** a `d`-level truncated oscillator has `[Q,P]=i(I-d|top><top|)`, while a collective spin has relative commutator error `n/S` at excitation number `n`. Large bundles of microscopic hinges can therefore approach continuous frame variables below a controlled boundary.
5. **Pure-stabilizer dynamics rejected:** exact local coordinate invariants begin at two derivatives and momentum invariants at one. Squaring them gives `H_h ~ k^4 h^2`, `H_p ~ k^2 p^2`, and `omega ~ k^3`, not the required linear dispersion.
6. **First-order frame–connection branch passes its algebraic gate:** introducing an independent 18-component torsion-free connection gives a local `C^2 + C*partial(h)` action. Eliminating the connection reproduces the Fierz–Pauli stiffness exactly. Across 500 random momenta, the Schur-complement, gauge-null, and two-TT-mode errors are all below `1.4e-15`.

## Current interpretation

The finite gravity constraint code is kinematic. A commuting-projector or manifest-local-invariant-square Hamiltonian is not the required graviton theory. The viable branch now keeps the frame and connection independent at the microscopic level, so an area-times-curvature or connection-times-frame-gradient interaction can form before the connection is eliminated.

This supplies a **local classical linear architecture with `omega ~ k`**, but not yet a finite quantum phase. The connection quadratic form is indefinite as expected for an auxiliary constrained variable; a finite transfer matrix or Hamiltonian must still prove a stable physical spectrum, cutoff control, and nonlinear closure.

The current work also does not derive `Z_g/Z_A`, chiral fermions, the particle-mass hierarchy, or control of the vacuum-volume term.

## Files

| File | Purpose |
|---|---|
| [`finite_hilbert_derivation.md`](finite_hilbert_derivation.md) | Finite construction, canonical no-go result, discrete constraints, and derivative-order obstruction. |
| [`finite_junction_model.py`](finite_junction_model.py) | Exact single-plaquette quantum-link test, GF(p) gravity constraint code, and finite canonical-pair diagnostics. |
| [`results.json`](results.json) | Frozen finite-link, finite-constraint, and canonical-pair results. |
| [`search_local_dynamics.py`](search_local_dynamics.py) | Polynomial nullspace search for exact local gauge invariants. |
| [`dynamics_results.json`](dynamics_results.json) | Frozen derivative-order result and `omega ~ k^3` rejection. |
| [`first_order_frame_connection.md`](first_order_frame_connection.md) | Derivation of the local independent-connection architecture and its limitations. |
| [`check_first_order_frame.py`](check_first_order_frame.py) | Numerical gamma-gamma identity, Schur-complement, gauge-null, and TT-spectrum checks. |
| [`first_order_results.json`](first_order_results.json) | Frozen 500-momentum first-order verification output. |

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
```

The scripts exit nonzero if any committed claim fails.

## Immediate next implementation

Build a finite frame-and-connection transfer matrix or noncommuting Hamiltonian on a small periodic lattice. It must retain the exact/discrete constraints or quantify emergent-symmetry leakage, show a stable two-helicity spectrum whose gap closes as `1/L`, and derive its coefficients from the same microscopic penalties and move amplitudes used by the electromagnetic quantum links.
