# Microscopic finite-state phase

**Status:** finite electromagnetic kinematics and a finite static Coulomb ground state, a finite chiral charged-endpoint prototype, a finite reversible photon–companion bridge, first-order frame–connection algebra, constrained local gravity reduction, the finite connection lock, and a finite dressed linear tensor Hamiltonian are checked. Relativistic photon dynamics, interacting QED, observed matter, shared coefficients, nonlinear closure, common-cone recovery, and interacting regulator universality remain open.  
**Date:** 2026-09-20  
**Active branch:** `main`

This folder replaces infrared rotor, symmetric-tensor, continuum-matter, and phenomenological conversion placeholders with explicit finite local candidates, then tests whether their constrained dynamics reaches the intended photon, helicity-2, charged-matter, and companion targets.

## Results obtained

1. **Finite electromagnetic link:** a spin-1 quantum link, with three states per oriented link, is the smallest tested representation with exact local U(1) Gauss symmetry and nonconstant electric-flux energy. Spin-1/2 preserves Gauss symmetry but has constant `E^2`.
2. **Finite static electromagnetic Coulomb ground state:** the spin-1 links now form a declared local 3+1D Hamiltonian with a positive weighted Rokhsar–Kivelson point. The exact periodic `2^3` zero-charge/zero-winding component contains `146,327` gauge-allowed states and `1,236,144` undirected plaquette transitions. Spin-1 `L=4,6,8` runs show Coulomb winding scaling, an equal-time rank-two transverse tensor with longitudinal fraction below `5.25e-33`, perimeter-favored Wilson-shift overlaps, and a fixed-charge response whose effective string slope falls by 56% from `L=6` to `L=8`. Spin 2 supplies a larger-representation control. This passes the **static Stage 6A** gate; a stable `z=1` photon phase and interacting QED remain open.
3. **Finite gravitational constraint skeleton:** six odd-prime qudits per site realize three vector constraints and one scalar constraint as an exact commuting Weyl/CSS algebra. The periodic code leaves two local logical modes per site plus the declared global sector.
4. **Finite canonical obstruction:** no finite matrices obey `[Q,P]=iI` on their full Hilbert space. Continuous frame variables must emerge below a cutoff or from a discrete microscopic algebra.
5. **Pure stabilizer branch rejected:** exact local coordinate invariants begin at two derivatives and momentum invariants at one. Squaring them gives `omega ~ k^3`, not a relativistic tensor branch.
6. **First-order frame–connection identity:** an independent torsion-free connection gives a local `C^2+C*partial(h)` form whose constrained elimination is exactly the Fierz–Pauli stiffness.
7. **Independent verification:** the first-order identity passes random numerical tests, exact rational coefficient comparison, gauge-null checks, finite-size scaling, static response, and a reduced Euclidean transfer benchmark.
8. **Naïve auxiliary quantization rejected:** a positive real auxiliary block with zero bare frame stiffness induces a nonpositive Schur complement. The accepted branch treats the connection as constrained/nonpropagating.
9. **Constrained local reduction completed:** the 48-dimensional phase space `(h,p,C,P_C)` has 36 second-class connection constraints and four first-class frame constraints:

   ```text
   48 - 36 - 2*4 = 4 physical phase-space dimensions
   ```

   The quotient produces two equal positive frequencies `|k_hat|` without constructing a TT projector. All 492 nonzero modes on `L=3,5,7` pass.
10. **Full real-space reduction completed:** complete periodic real-space matrices, without Fourier reduction during the constraint calculation, leave exactly two positive modes per nonzero momentum on `3^3` and `5^3` lattices.
11. **Global conformal issue isolated:** the unrestricted periodic zero mode has one negative homogeneous trace-momentum direction.
12. **Fixed-volume candidate passed:** total-volume and trace-momentum constraints remove exactly that global canonical pair at linear order.
13. **Reduced finite quantum target passed:** two truncated physical oscillators give a positive transfer matrix, a twofold first excitation, `omega ~ L^-0.993857`, negligible truncation-boundary occupation, and a positive Euclidean reflection kernel to numerical precision.
14. **Finite chiral charged endpoint passed at prototype scope:** finite spin-1 flux bundles and exact odd-strand binding produce fermionic endpoint operators. The first primitive anomaly-free chiral set in the declared search is `(-11,-5,-1,-1,9,9)`. Hopping is exactly Gauss covariant; a finite domain-wall slab has one Weyl cone per wall, seven gapped physical doublers, opposite wall chirality, and one universal frame derivative. A shared two-parameter localization rule produces gaps spanning `3.245e7`. The mirror wall and interacting continuum remain open.
15. **Finite reversible photon–companion bridge passed at architecture scope:** a source, two photon-frequency bins, neutral `chi`, matter recoil, a computed bound `chi` mode, capture recoil, a receiver, and one constrained frame form a reversible finite Hamiltonian. Energy and component ledgers close to roundoff and motion/lensing share the same frame. The finite speed ratio `v_chi/v_gamma=0.3571428571` keeps common-cone recovery open.
16. **Finite dressed linear gravity Hamiltonian passed:** the exact scalar/vector stabilizer quotient produces two logical Weyl pairs at `p=5,7,11` without a TT projector. Two positive clock/Villain copies in the exact dressed frame algebra produce signed-momentum gaps

   ```text
   Delta = c_g |k_hat| + O(|k_hat|^3)
   ```

   with exact two-polarization degeneracy, maximum fit residual `1.67e-3`, minimum lock/tensor gap ratio `2.49`, connection-band leakage `5.16e-16`, compact-cut occupation falling to `1.89e-4` at `p=11`, and full-symbol polarization splitting below `2.84e-16`.

## Current interpretation

The electromagnetic result now closes a previously missing **static finite phase** gate:

```text
finite spin-1 links
    -> exact U(1) Gauss law
    -> positive local RK Hamiltonian
    -> exact 2^3 gauge-sector closure
    -> Coulomb winding response
    -> rank-two equal-time transverse tensor
    -> perimeter-favored Wilson shifts
    -> decreasing fixed-charge string slope
```

This is not yet a photon spectrum. The equal-time rank-two tensor cannot be relabeled as two linearly dispersing quantum branches, and the winding-sector free-energy response is not a `1/L` photon gap. Issue #6 remains open for a finite coupling interval with `z=1` dynamics, charged-matter propagation, Ward identities, vacuum polarization, and QED scaling.

The local **linear gravity continuum and finite-regulator bridges are closed under their stated assumptions**:

```text
local frame + independent connection
    -> 36 second-class connection constraints
    -> exact dressed frame Weyl algebra
    -> 4 first-class finite frame constraints
    -> exact finite stabilizer quotient
    -> 2 logical Weyl fields
    -> positive compact Hamiltonian / constrained transfer matrix
    -> 2 linearly dispersing tensor branches
```

No TT projector is inserted. The logical fields are derived by modular kernel/row-space reduction. The parent finite moves use

\[
\overline Z_i=Z_{h_i},\qquad
\overline X_i=X_{h_i}\prod_aX_{C_a}^{A_{ai}},
\]

so they commute exactly with the connection lock. The accepted object is a constrained transfer-matrix route; the retained frame-only invariant-square branch still gives `omega ~ k^3` and remains rejected.

The self-dual clock ratio and overall gravity scale are regulator inputs. Issue #3 must derive shared coefficients, issue #5 must add nonlinear self-coupling and constraint closure, and issue #8 must establish common-cone recovery and regulator universality.

The matter and companion constructions are also finite architecture prototypes rather than complete particle physics or astrophysics. The matter mirror, observed spectrum, bound-state clocks, many-body deposits, continuum rates, and radiative stability remain explicit.

## Files

| File | Purpose |
|---|---|
| [`finite_hilbert_derivation.md`](finite_hilbert_derivation.md) | Finite construction, canonical no-go result, discrete constraints, and derivative-order obstruction |
| [`finite_junction_model.py`](finite_junction_model.py) / [`results.json`](results.json) | Quantum-link and GF(p) kinematic checks |
| [`finite_em_coulomb_phase.md`](finite_em_coulomb_phase.md) | Issue-#6 Stage 6A Hamiltonian, derivation, frozen diagnostics, claim boundary, and next rejection test |
| [`check_finite_em_coulomb.py`](check_finite_em_coulomb.py) | Exact cube, RK flux/transverse, Wilson-shift, and fixed-charge calculations |
| [`finite_em_coulomb_results.json`](finite_em_coulomb_results.json) | Master static Coulomb-stage result and unresolved gates |
| [`finite_em_exact_cube_results.json`](finite_em_exact_cube_results.json) | Exact periodic `2^3` gauge-sector closure |
| [`finite_em_flux_results.json`](finite_em_flux_results.json) | Spin-1/spin-2 winding, transverse, and Wilson-shift finite-size results |
| [`finite_em_charge_results.json`](finite_em_charge_results.json) | Combined `L=6,8` fixed-charge potential and string-slope scaling |
| [`chiral_matter_defect.md`](chiral_matter_defect.md) | Finite charged endpoint, anomaly search, domain-wall chirality, protected gaps, and claim boundary |
| [`check_chiral_matter_defect.py`](check_chiral_matter_defect.py) / [`chiral_matter_results.json`](chiral_matter_results.json) | Executable matter prototype checks |
| [`check_companion_bridge.py`](check_companion_bridge.py) / [`companion_bridge_results.json`](companion_bridge_results.json) | Finite reversible conversion, binding, recoil, ledger, and frame checks |
| [`search_local_dynamics.py`](search_local_dynamics.py) / [`dynamics_results.json`](dynamics_results.json) | Exact-local invariant search and the `omega ~ k^3` rejection |
| [`first_order_frame_connection.md`](first_order_frame_connection.md) | Independent frame–connection derivation |
| [`check_first_order_frame.py`](check_first_order_frame.py), [`check_first_order_exact.py`](check_first_order_exact.py) | Floating-point and exact rational first-order checks |
| [`auxiliary_stability_no_go.md`](auxiliary_stability_no_go.md) | Schur-sign obstruction for naïve positive auxiliary quantization |
| [`check_tensor_scaling.py`](check_tensor_scaling.py) | Tensor gap, polarization, and cutoff-anisotropy scaling |
| [`reduced_tt_transfer.md`](reduced_tt_transfer.md) | Reduced finite transfer-matrix and reflection-positivity benchmark |
| [`constrained_local_reduction.md`](constrained_local_reduction.md) | Full Dirac/symplectic reduction of `(h,p,C,P_C)` without a TT projector |
| [`check_constrained_local_reduction.py`](check_constrained_local_reduction.py) | Complete nonzero-mode reduction on `L=3,5,7` |
| [`real_space_reduction.md`](real_space_reduction.md) | Full periodic real-space reduction and global-sector diagnosis |
| [`global_volume_constraint.md`](global_volume_constraint.md) | Fixed-volume removal of the homogeneous conformal pair |
| [`finite_second_class_lock.md`](finite_second_class_lock.md) | Finite relative-connection oscillator and exact dressed frame algebra |
| [`check_finite_second_class_lock.py`](check_finite_second_class_lock.py) | Lock gap, factorization, and dressed-Weyl checks |
| [`finite_dressed_frame_hamiltonian.md`](finite_dressed_frame_hamiltonian.md) | Issue-#2 finite constrained Hamiltonian/transfer construction and acceptance accounting |
| [`check_finite_dressed_frame_hamiltonian.py`](check_finite_dressed_frame_hamiltonian.py) / [`finite_dressed_frame_results.json`](finite_dressed_frame_results.json) | Exact logical reduction, sparse finite spectra, lock leakage, boundary control, and full-symbol checks |
| [`verification_routes.md`](verification_routes.md) | Independent verification methods and acceptance criteria |

Every executable calculation has a neighboring frozen `*_results.json` file.

## Reproduce the latest checks

From the repository root:

```sh
python phase_junction_network/microscopic/check_finite_em_coulomb.py \
  --quick \
  --output /tmp/finite_em_quick_results.json

python phase_junction_network/microscopic/check_finite_em_coulomb.py \
  --exact-only \
  --output /tmp/finite_em_exact_cube_results.json

python phase_junction_network/microscopic/check_chiral_matter_defect.py \
  --slab-width 12 \
  --output phase_junction_network/microscopic/chiral_matter_results.json

python phase_junction_network/microscopic/check_companion_bridge.py \
  --output phase_junction_network/microscopic/companion_bridge_results.json

python phase_junction_network/microscopic/check_finite_dressed_frame_hamiltonian.py \
  --output phase_junction_network/microscopic/finite_dressed_frame_results.json

python phase_junction_network/microscopic/check_constrained_local_reduction.py \
  --sizes 3,5,7 \
  --output phase_junction_network/microscopic/constrained_local_reduction_results.json

python phase_junction_network/microscopic/check_real_space_reduction.py \
  --sizes 3,5 \
  --output phase_junction_network/microscopic/real_space_reduction_results.json

python phase_junction_network/microscopic/check_global_volume_constraint.py \
  --output phase_junction_network/microscopic/global_volume_constraint_results.json
```

Every script exits nonzero if a committed claim fails.

## Immediate next implementation

1. Detune the finite spin-1 Hamiltonian around the RK point and demonstrate or reject a finite interval with two dynamical transverse gaps scaling as `1/L`, no scalar branch, and consistent spin-2 behavior.
2. Embed the actual issue-#4 endpoint hopping in that dynamical phase and test finite Ward identities, charge renormalization, and vacuum polarization.
3. Derive `r0`, `eta`, `U_A`, `K_A`, `U_g`, `K_g`, the gravity clock ratio, and companion coefficients from one microscopic move set under issue #3.
4. Extend the frame/connection constraints to nonlinear order, include gravitational self-energy, and derive or reject the fixed-volume rule under issue #5.
5. Test mirror completion, interacting anomaly accounting, common-cone recovery, and regulator universality under issue #8.
6. Derive many-body companion capacity, lifetime, release, and self-gravity before any galaxy-scale interpretation.
7. Maintain the parameter and observable ledger under issue #9 and freeze predictions before data fitting.

A later branch fails if it hides a mirror, violates finite Gauss symmetry, mistakes an equal-time Coulomb tensor for a photon spectrum, reintroduces a TT projector, allows a low-energy connection/scalar/vector state, loses positivity, or obtains a common photon–matter–companion–gravity result only by independently tuning every sector.
