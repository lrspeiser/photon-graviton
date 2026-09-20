# Microscopic finite-state phase

**Status:** finite kinematics, first-order frame–connection algebra, constrained local reduction, full real-space propagation, and a reduced finite transfer matrix are checked. A finite **local** quantum implementation and nonlinear closure remain open.  
**Date:** 2026-09-20  
**Active branch:** `main`

This folder replaces the infrared rotor and symmetric-tensor fields with explicit finite local candidates and progressively tests whether their constrained dynamics reaches the required photon and helicity-2 continuum theory.

## Results obtained

1. **Finite electromagnetic link:** a spin-1 quantum link, with three states per oriented link, is the smallest tested representation with exact local U(1) Gauss symmetry and nonconstant electric-flux energy. Spin-1/2 preserves Gauss symmetry but has constant `E^2`.
2. **Finite gravitational constraint skeleton:** six odd-prime qudits per site realize the three vector constraints and one scalar constraint as an exact commuting Weyl/CSS algebra. The periodic code leaves two local logical modes per site, with the expected global zero-mode excess.
3. **Finite canonical obstruction:** no finite matrices obey `[Q,P]=iI` on their full Hilbert space. Continuous frame variables must emerge below a cutoff or from a discrete microscopic algebra.
4. **Pure stabilizer branch rejected:** exact local coordinate invariants begin at two derivatives and momentum invariants at one. Squaring them gives `omega ~ k^3`, not a relativistic tensor branch.
5. **First-order frame–connection identity:** an independent torsion-free connection gives a local `C^2 + C*partial(h)` form whose constrained elimination is exactly the Fierz–Pauli stiffness.
6. **Independent verification:** the first-order identity passes random numerical tests, exact rational coefficient comparison, gauge-null checks, finite-size scaling, static response, and a reduced Euclidean transfer-matrix benchmark.
7. **Naïve auxiliary quantization rejected:** a positive real auxiliary block with zero bare frame stiffness induces a nonpositive Schur complement. The connection must be constrained/nonpropagating, contour-integrated, or accompanied by explicit bare frame stiffness.
8. **Constrained local reduction completed:** starting from the 48-dimensional phase space `(h,p,C,P_C)`, the 36 connection constraints are second-class and the four frame constraints are first-class:

   ```text
   48 - 36 - 2*4 = 4 physical phase-space dimensions
   ```

   The quotient produces two equal positive frequencies `|k_hat|` without constructing a TT projector. All 492 nonzero modes on `L=3,5,7` pass.
9. **Full real-space reduction completed:** the reduction was repeated on complete periodic real-space matrices, without Fourier block reduction or a TT projector. The `3^3` and `5^3` lattices leave exactly two positive modes per nonzero momentum; the reduced spectra agree with the central-difference lattice symbols to `2.21e-14` or better.
10. **Global conformal issue isolated:** the unrestricted periodic zero mode has one negative homogeneous trace-momentum direction.
11. **Fixed-volume candidate passed:** imposing total-volume and trace-momentum constraints as a second-class pair removes exactly that negative global canonical pair, leaving five positive global shear momenta and five zero-potential global shape moduli.
12. **Reduced finite quantum target passed:** after the constraints are solved, two truncated oscillators give a positive transfer matrix, a twofold first excitation, `omega ~ L^-0.993857`, negligible truncation-boundary occupation, and a positive Euclidean reflection kernel to numerical precision.

## Current interpretation

The local **linear continuum bridge is now closed**:

```text
local frame + independent connection
    -> 36 second-class connection constraints
    -> 4 first-class frame constraints
    -> 4-dimensional physical phase space
    -> 2 positive linearly dispersing modes
```

This result no longer depends on inserting a transverse-traceless projector. It has also been reproduced in full real space.

The remaining microscopic problem is more specific: implement the second-class connection reduction and the four first-class frame constraints in a finite local Hilbert space, then derive the positive reduced transfer matrix from those finite local variables.

The connection is not an ordinary collection of positive-energy oscillators. The constrained first-order/Palatini route remains primary; complex Euclidean auxiliary and positive-bare-stiffness models remain independent cross-checks.

The global periodic conformal mode is controlled at linear order by a fixed-volume second-class pair. A microscopic reason for that rule and its nonlinear closure are still required.

## Files

| File | Purpose |
|---|---|
| [`finite_hilbert_derivation.md`](finite_hilbert_derivation.md) | Finite construction, canonical no-go result, discrete constraints, and derivative-order obstruction. |
| [`finite_junction_model.py`](finite_junction_model.py) / [`results.json`](results.json) | Quantum-link, GF(p) gravity-constraint, and finite canonical-pair checks. |
| [`search_local_dynamics.py`](search_local_dynamics.py) / [`dynamics_results.json`](dynamics_results.json) | Exact-local invariant search and the `omega ~ k^3` rejection. |
| [`first_order_frame_connection.md`](first_order_frame_connection.md) | Independent frame–connection derivation. |
| [`check_first_order_frame.py`](check_first_order_frame.py), [`check_first_order_exact.py`](check_first_order_exact.py) | Floating-point and exact rational first-order checks. |
| [`auxiliary_stability_no_go.md`](auxiliary_stability_no_go.md) | Schur-sign obstruction for naïve positive auxiliary quantization. |
| [`check_tensor_scaling.py`](check_tensor_scaling.py) | Tensor gap, polarization, and cutoff-anisotropy scaling. |
| [`reduced_tt_transfer.md`](reduced_tt_transfer.md) | Reduced finite transfer-matrix and reflection-positivity target. |
| [`constrained_local_reduction.md`](constrained_local_reduction.md) | Full Dirac/symplectic reduction of `(h,p,C,P_C)` without a TT projector. |
| [`check_constrained_local_reduction.py`](check_constrained_local_reduction.py) | Complete nonzero-mode reduction on `L=3,5,7`. |
| [`real_space_reduction.md`](real_space_reduction.md) | Full periodic real-space reduction and global-sector diagnosis. |
| [`check_real_space_reduction.py`](check_real_space_reduction.py) | Dense real-space checks on `3^3` and `5^3` lattices. |
| [`global_volume_constraint.md`](global_volume_constraint.md) | Fixed-volume removal of the homogeneous conformal pair. |
| [`check_global_volume_constraint.py`](check_global_volume_constraint.py) | Exact global constraint and inertia check. |
| [`verification_routes.md`](verification_routes.md) | Independent verification methods and acceptance criteria. |

Every calculation has a neighboring frozen `*_results.json` file.

## Reproduce the latest checks

From the repository root:

```sh
python phase_junction_network/microscopic/check_constrained_local_reduction.py \
  --sizes 3,5,7 \
  --output phase_junction_network/microscopic/constrained_local_reduction_results.json

python phase_junction_network/microscopic/check_real_space_reduction.py \
  --sizes 3,5 \
  --output phase_junction_network/microscopic/real_space_reduction_results.json

python phase_junction_network/microscopic/check_global_volume_constraint.py \
  --output phase_junction_network/microscopic/global_volume_constraint_results.json

python phase_junction_network/microscopic/check_reduced_tt_transfer.py \
  --output phase_junction_network/microscopic/reduced_tt_results.json
```

The earlier scripts remain independently reproducible. Every script exits nonzero if a committed claim fails.

## Immediate next implementation

1. Encode `P_C=0` and `C-Gamma[h]=0` as finite local constraints using the odd-prime Weyl or protected collective-spin variables.
2. Perform the finite constrained reduction before any continuum approximation and verify that the reduced transfer matrix matches the two-mode benchmark.
3. Quantify constraint leakage and boundary-state occupation as local dimension increases.
4. Extend the fixed-volume pair to the first nonlinear order and verify closure with gravitational self-energy.
5. Derive `U_A`, `K_A`, `U_g`, and `K_g` from one microscopic move set, rather than choosing `Z_g/Z_A` independently.

A branch fails if it needs a nonlocal TT projector, retains a propagating connection/scalar mode, has a negative physical energy, or obtains the right spectrum only by independently tuning the electromagnetic and gravitational sectors.