# Microscopic finite-state phase

**Status:** first finite-state kinematic prototype  
**Date:** 2026-09-20  
**Branch:** `agent/phase-junction-microscopic`

This folder begins the next Phase Junction phase: replacing the infrared rotor and symmetric-tensor fields with explicit finite local Hilbert spaces.

## Results obtained

1. **Electromagnetism:** a spin-1 quantum link, with three states per oriented link, is the smallest tested representation that has exact local U(1) Gauss symmetry and a non-constant electric-flux energy. Spin-1/2 preserves Gauss symmetry but its `E^2` term is constant, so it cannot supply the required electric stiffness by itself.
2. **Gravity kinematics:** six odd-prime qudits per lattice site can realize the three vector constraints and one scalar curvature constraint as an exact commuting Weyl/CSS algebra. On periodic odd lattices, away from modular degeneracies, the code leaves `2*N + 4` logical qudits: two local modes per site and four global zero modes.
3. **Finite-dimensional obstruction:** no finite local matrices can obey an exact continuous canonical commutator `[Q,P]=iI` over their full Hilbert space. The trace of a commutator is zero, while the trace of `iI` is not. Continuous linearized diffeomorphism symmetry must therefore be emergent in a low-energy sector, or replaced microscopically by an exact discrete constraint algebra.
4. **Low-energy escape:** a `d`-level truncated oscillator satisfies `[Q,P]=i(I-d|top><top|)`. It is exactly canonical inside matrix elements that avoid the top state. A collective spin has relative commutator error `n/S` at excitation number `n`, so a large bundle of microscopic hinges can approach the continuous frame variables.
5. **Dynamics obstruction identified:** an exact polynomial search finds that local coordinate observables invariant under the vector gauge symmetry first appear at two derivatives, while local momentum observables invariant under the scalar symmetry first appear at one derivative. Squaring those manifest local invariants gives `H_h ~ k^4 h^2`, `H_p ~ k^2 p^2`, and therefore `omega ~ k^3`, not the required linear dispersion.

## What this does not solve

The exact finite constraint code is kinematic. A commuting-projector or manifest-local-invariant-square Hamiltonian is not by itself the required graviton theory: the latter naturally produces `omega ~ k^3`. The next gate is a finite implementation of the globally gauge-invariant Fierz–Pauli bilinear, an area-weighted term linear in frame curvature, or another emergent-symmetry mechanism that gives exactly two helicity-2 modes with linear long-wavelength dispersion and no additional gapless scalar or vector modes.

The current work also does not derive the microscopic ratio `Z_g/Z_A`, chiral fermions, the particle-mass hierarchy, nonlinear constraint closure, or control of the vacuum-volume term.

## Files

| File | Purpose |
|---|---|
| [`finite_hilbert_derivation.md`](finite_hilbert_derivation.md) | Mathematical construction, no-go result, finite alternatives, and next dynamics target. |
| [`finite_junction_model.py`](finite_junction_model.py) | Exact single-plaquette quantum-link test, GF(p) gravity constraint code, and finite canonical-pair diagnostics. |
| [`results.json`](results.json) | Frozen output of the finite-link, finite-constraint, and canonical-pair run. |
| [`search_local_dynamics.py`](search_local_dynamics.py) | Polynomial nullspace search for the minimum derivative order of exact local gauge invariants. |
| [`dynamics_results.json`](dynamics_results.json) | Frozen derivative-order result and predicted `omega ~ k^3` obstruction. |

## Reproduce

From the repository root:

```sh
python phase_junction_network/microscopic/finite_junction_model.py \
  --output phase_junction_network/microscopic/results.json

python phase_junction_network/microscopic/search_local_dynamics.py \
  --output phase_junction_network/microscopic/dynamics_results.json
```

The scripts exit nonzero if any claimed finite-state result fails.
