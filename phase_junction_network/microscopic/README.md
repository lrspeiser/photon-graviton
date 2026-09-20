# Microscopic finite-state phase

**Status:** finite electromagnetic kinematics, a finite chiral charged-endpoint prototype, first-order frame–connection algebra, constrained local gravity reduction, full real-space propagation, and a reduced finite transfer matrix are checked. A deconfined finite QED phase, a complete finite local gravity Hamiltonian, mirror/continuum completion, shared coefficients, and nonlinear closure remain open.  
**Date:** 2026-09-20  
**Active branch:** `main`

This folder replaces infrared rotor, symmetric-tensor, and continuum-matter placeholders with explicit finite local candidates, then tests whether their constrained dynamics reaches the required photon, helicity-2, and chiral-matter continuum targets.

## Results obtained

1. **Finite electromagnetic link:** a spin-1 quantum link, with three states per oriented lane, is the smallest tested representation with exact local U(1) Gauss symmetry and nonconstant electric-flux energy. Spin-1/2 preserves Gauss symmetry but has constant `E^2`.
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
13. **Finite chiral charged endpoint passed at prototype scope:** finite spin-1 flux bundles and exact odd-strand binding produce fermionic endpoint operators. An exhaustive odd-charge search through six species and `|q|<=11` yields the primitive anomaly-free chiral set `(-11,-5,-1,-1,9,9)`. Composite hopping is exactly Gauss covariant; a finite domain-wall slab has one Weyl cone per wall, seven gapped physical doublers, opposite wall chirality, particle/antiparticle pairing, and one universal frame derivative. A common two-parameter localization rule produces exponentially protected finite-width gaps spanning `3.245e7`. The result remains in the one-cone phase throughout twenty nearby Wilson-mass samples. The remote mirror wall, observed spectrum, and interacting continuum remain open.

## Current interpretation

The local **linear gravity continuum bridge is closed under its stated assumptions**:

```text
local frame + independent connection
    -> 36 second-class connection constraints
    -> 4 first-class frame constraints
    -> 4-dimensional physical phase space
    -> 2 positive linearly dispersing modes
```

This result no longer depends on inserting a transverse-traceless projector and has been reproduced in full real space. The remaining gravity problem is more specific: implement the second-class connection reduction and four first-class frame constraints in a finite local Hilbert space, then derive the positive reduced transfer matrix from those finite local variables.

The first matter bridge is also now explicit:

```text
spin-1 finite flux lanes
    + odd microscopic fermionic strands
    + exact equal-occupation binding band
    -> finite charged fermionic endpoint
    -> anomaly-free chiral wall spectrum
    -> exponentially protected opposite-wall overlap gap
```

The matter result is a regulator/prototype, not a Standard Model derivation. It supplies a finite charged defect for issues #6 and #7 and a concrete mass-protection mechanism for issue #3. Issue #8 still owns mirror completion, interacting Ward identities, radiative stability, and continuum Lorentz recovery.

The connection is not an ordinary collection of positive-energy oscillators. The constrained first-order/Palatini route remains primary; complex Euclidean auxiliary and positive-bare-stiffness models remain independent cross-checks.

The global periodic conformal mode is controlled at linear order by a fixed-volume second-class pair. A microscopic reason for that rule and its nonlinear closure are still required.

## Files

| File | Purpose |
|---|---|
| [`finite_hilbert_derivation.md`](finite_hilbert_derivation.md) | Finite construction, canonical no-go result, discrete constraints, and derivative-order obstruction. |
| [`finite_junction_model.py`](finite_junction_model.py) / [`results.json`](results.json) | Quantum-link, GF(p) gravity-constraint, and finite canonical-pair checks. |
| [`chiral_matter_defect.md`](chiral_matter_defect.md) | Finite charged endpoint, anomaly-free charge search, domain-wall chirality, protected gaps, and claim boundary. |
| [`check_chiral_matter_defect.py`](check_chiral_matter_defect.py) / [`chiral_matter_results.json`](chiral_matter_results.json) | Frozen finite-Hilbert, Gauss, anomaly, dispersion, hierarchy, frame-coupling, and stability checks. |
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

Every executable calculation has a neighboring frozen `*_results.json` file.

## Reproduce the latest checks

From the repository root:

```sh
python phase_junction_network/microscopic/check_chiral_matter_defect.py \
  --slab-width 12 \
  --output phase_junction_network/microscopic/chiral_matter_results.json

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

1. Integrate the finite defects into the actual spin-1 many-link Hamiltonian and test the finite deconfined Coulomb/QED phase under issue #6.
2. Build the finite matter-assisted photon–`varphi` transition, including reverse channel and recoil, under issue #7.
3. Derive `r0`, `eta`, `U_A`, `K_A`, `U_g`, and `K_g` from one microscopic move set under issue #3 rather than choosing each sector independently.
4. Encode `P_C=0` and `C-Gamma[h]=0` as finite local constraints and complete the finite dressed gravity Hamiltonian under issue #2.
5. Test mirror-wall completion, interacting anomaly accounting, Ward identities, common-cone recovery, and radiative stability under issue #8.
6. Extend the fixed-volume pair to first nonlinear order and verify closure with gravitational self-energy under issue #5.

A matter branch fails if its charge is only a label, its hopping violates the finite Gauss law, it obtains a small gap from an arbitrary tiny onsite energy, or it hides a mirror/doubler sector. A gravity branch fails if it needs a nonlocal TT projector, retains a propagating connection/scalar mode, has negative physical energy, or obtains the desired spectrum only by independently tuning electromagnetic and gravitational sectors.
