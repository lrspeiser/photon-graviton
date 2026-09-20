# Phase Junction Network

**Status:** exploratory candidate foundation. Four finite architecture gates now have explicit implementations: issue #2 closes the linear dressed-frame regulator, issue #4 supplies a finite chiral charged-endpoint prototype, issue #7 supplies a reversible photon–companion bridge, and issue #6 now has a finite **static** spin-1 Coulomb ground-state anchor. Relativistic photon dynamics and interacting QED, observed matter, shared microscopic coefficients, nonlinear gravity, common-cone recovery, interacting-continuum universality, many-body deposits, and empirical predictions remain open.  
**Created:** 2026-09-20  
**Architecture audit:** [`architecture_audit.md`](architecture_audit.md)  
**Issue #2 closure:** [`microscopic/finite_dressed_frame_hamiltonian.md`](microscopic/finite_dressed_frame_hamiltonian.md)  
**Issue #6 Stage 6A:** [`microscopic/finite_em_coulomb_phase.md`](microscopic/finite_em_coulomb_phase.md)  
**Issue #7 closure:** [`companion_integration.md`](companion_integration.md)

## Working picture

The physical comparison is placed on links between neighboring matter–geometry junctions rather than in a smooth scalar gradient. Electromagnetism is compact scalar phase holonomy; gravity compares full local frames and clocks; matter is represented by finite endpoint defects; and the selected companion is a neutral relative-phase quantum.

| Object | Role |
|---|---|
| Photon | Transverse quantum of the compact U(1) link connection |
| Companion `chi=varphi` | Neutral relative-phase quantum `theta_m-theta_g` |
| Deposit | Bound state of the same `chi` sector |
| Graviton | Transverse dressed-frame excitation, distinct from the companion reservoir |
| Frame connection | Constrained auxiliary comparison variable with a gapped relative sector |
| Matter | Finite charged or neutral defects; the observed spectrum remains open |

## Established internally

### Finite static electromagnetic Coulomb ground state

A spin-1, three-state link is the smallest tested finite representation with exact local Gauss symmetry and nonconstant electric energy. Issue #6 now goes beyond that kinematic statement with the declared finite Hamiltonian

\[
H_A(u,t,v)=\frac{u}{2}\sum_\ell E_\ell^2
-t\sum_p(W_p+W_p^\dagger)+v\sum_pD_p,
\]

where

\[
D_p=\sqrt{W_p^\dagger W_p}+\sqrt{W_pW_p^\dagger}.
\]

At `u=0, v=t`, every plaquette term is a positive weighted graph Laplacian and commutes exactly with Gauss law. The committed checks find:

- an exact periodic `2^3` spin-1 component with `146,327` gauge-allowed states and `1,236,144` undirected plaquette transitions;
- zero sampled Gauss residual and zero equal-amplitude RK residual;
- spin-1 winding scaling over `L=4,6,8`, with only 6.05% variation in `<W^2>/L`;
- an equal-time rank-two transverse correlation tensor with longitudinal fraction below `5.25e-33`;
- perimeter-favored finite Wilson-shift overlaps at every tested size;
- a fixed opposite-charge response consistent with the periodic lattice Green function;
- an effective linear-string slope that falls 56% from `L=6` to `L=8`;
- persistence of the static scaling in a spin-2 control.

This passes the **finite static Coulomb Stage 6A** gate. It does **not** yet demonstrate a stable relativistic photon phase: the equal-time transverse tensor is not a dynamical spectrum, the winding-sector free energy is not a photon gap, and the solvable RK point need not have `z=1`. Issue #6 therefore remains open for a detuned coupling interval with two linearly dispersing transverse branches, a `1/L` photon gap, no scalar, gauge-covariant dynamical matter, Ward identities, vacuum polarization, and QED scaling.

### Linear gravity and its finite regulator

The continuum, exact-arithmetic, finite-state, full-real-space, and finite-size checks establish under their linear assumptions:

- three vector and one scalar frame constraints;
- a local independent connection whose constrained elimination gives the Fierz–Pauli stiffness;
- 36 second-class connection constraints and four first-class frame constraints;
- four physical phase-space dimensions, or two positive linearly dispersing modes, per nonzero momentum;
- no inserted transverse-traceless projector;
- a finite odd-prime Weyl connection lock and exact dressed-frame algebra.

Issue #2 adds the missing finite dynamical bridge. The exact finite scalar/vector stabilizer quotient yields two logical Weyl pairs. On that quotient, two identical positive compact clock/Villain Hamiltonians give the two tensor polarizations:

\[
H_{\rm tensor}=\sum_{a=1}^{2}\frac{\lambda_g}{2}\sum_z
\left[(2-X_{a,z}-X_{a,z}^{\dagger})+
(2-Z_{a,z+1}Z_{a,z}^{\dagger}-Z_{a,z}Z_{a,z+1}^{\dagger})\right].
\]

The parent moves use only the exact dressed operators

\[
\overline Z_i=Z_{h_i},\qquad
\overline X_i=X_{h_i}\prod_aX_{C_a}^{A_{ai}},
\]

and therefore commute with the connection lock. The finite checks find:

- exactly two logical Weyl pairs at `p=5,7,11`;
- a twofold first tensor gap;
- signed-momentum gaps `Delta=c_g|k_hat|+O(|k_hat|^3)`;
- maximum linear-plus-cubic fit residual `1.67e-3`;
- connection-lock leakage `5.16e-16`;
- minimum lock-gap/tensor-gap ratio `2.49`;
- compact-cut bond occupation falling to `1.89e-4` at `p=11`;
- full three-dimensional polarization split below `2.84e-16` and anisotropy recovering approximately as `L^-2`.

This closes the **linear finite constrained regulator** gate. The overall gravitational scale and self-dual clock ratio remain regulator inputs, not derived predictions. Nonlinear closure and universality remain issues #5 and #8.

The periodic homogeneous conformal pair is removed at linear order by fixed-volume and trace-momentum constraints. That is not a vacuum-energy or cosmological-constant solution.

### Finite charged matter prototype

Issue #4 supplies a finite regulator/prototype using spin-1 flux bundles and bound odd-strand endpoint operators. Within the declared search class, the first primitive anomaly-free chiral spectrum is

\[
q=(-11,-5,-1,-1,9,9),\qquad \sum q=\sum q^3=0.
\]

Composite hopping changes exactly `|q|` finite flux lanes and commutes with the endpoint Gauss generators. A finite open Wilson/domain-wall slab supplies one Weyl cone on each wall and gaps the other seven physical Brillouin corners. A shared two-parameter localization rule gives four finite-width gaps spanning approximately `3.245e7`.

This is not the observed matter sector. The remote mirror wall, Standard Model gauge structure, generations, bound-state clocks, observed masses, interacting Ward identities, and radiative stability remain open. See [`microscopic/chiral_matter_defect.md`](microscopic/chiral_matter_defect.md).

### Finite photon–companion bridge

Issue #7 is closed at the finite architecture threshold by [`companion_integration.md`](companion_integration.md), [`microscopic/check_companion_bridge.py`](microscopic/check_companion_bridge.py), and its frozen output.

The calculation contains a finite source, two frequency bins of one photon field, a neutral `chi`, converter recoil, a computed bound `chi` mode, capture recoil, a receiver, and an exactly reduced constrained frame. Every transition has its Hermitian reverse. It closes energy and component ledgers to roundoff and uses no independent motion or lensing multipliers.

The finite benchmark also finds

```text
maximum chi speed / maximum photon speed = 0.3571428571
```

so a common relativistic cone is not established. The large finite conversion probability is an architecture stress test, not an astrophysical rate.

## What remains open

The project does not yet contain one finite local theory that simultaneously supplies:

- a stable `z=1` two-polarization photon phase and controlled interacting QED limit;
- completion of the matter prototype into the observed charged-matter spectrum without an unwanted mirror;
- one microscopic move set deriving electromagnetic, companion, matter, and gravitational coefficients;
- nonlinear gravitational self-coupling and constraint closure;
- a unitary or reflection-positive interacting continuum with regulator universality;
- a common causal cone for photons, companions, tensor modes, and matter;
- a microscopic volume/vacuum and background mechanism;
- many-body deposit capacity, lifetime, release, and self-gravity;
- a frozen cross-sector prediction that survives data.

The finite electromagnetic Stage 6A result is a static ground-state phase anchor, not completed QED. The finite gravity result is a regulator of the established **linear constrained theory**, not a unique ultraviolet completion.

## Work program

| Issue | Pillar | Status |
|---:|---|---|
| #2 | Finite dressed frame–connection Hamiltonian | **Linear finite regulator closed**; coefficient derivation, nonlinear completion, and universality delegated to #3, #5, and #8 |
| #3 | Shared microscopic coefficients and `Z_g/Z_A` | Open |
| #4 | Chiral matter defects and protected mass hierarchy | Finite prototype completed; observed spectrum and continuum completion open |
| #5 | Nonlinear gravity, self-coupling, and volume/vacuum term | Open |
| #6 | Finite electromagnetic Coulomb/QED phase | **Static RK Coulomb Stage 6A passed**; detuned `z=1` photon dynamics, dynamical matter, Ward identities, and QED continuum open |
| #7 | Photon–companion integration | **Finite architecture closed**; physical scaling and many-body completion open |
| #8 | Continuum quantum consistency and Lorentz/common-cone recovery | Open |
| #9 | Frozen predictions and empirical test ladder | Open |

```text
finite gravity regulator (#2 closed)
finite matter prototype (#4 scope passed)
finite photon-companion architecture (#7 closed)
finite static EM Coulomb anchor (#6 Stage 6A)
    -> detuned z=1 photon phase and dynamical defects (#6)
    -> shared microscopic coefficients (#3)
    -> nonlinear gravity and many-body binding (#5)
    -> interacting continuum, mirror completion, and common cone (#8)
    -> frozen predictions and data (#9)
```

## Relationship to the wider repository

The root project remains governed by `CURRENT-STATUS.md`, `research_plan/active-goal.md`, `research_plan/universe-contract.md`, and the solution-goal ledger. It requires fixed published distances, an operationally nonexpanding universe, complete energy accounting, no inserted dark halo, and joint predictions across redshift, timing, brightness, motion, lensing, gravitational waves, and background observations.

The finite constructions above supply architecture pieces. They do not yet derive a phenomenological redshift coefficient, continuum capture rate, galaxy-scale deposit, universal weak-field correction, or observational fit.

## Key files

| File | Purpose |
|---|---|
| [`architecture_audit.md`](architecture_audit.md) | Current pillar audit, dependencies, claim boundaries, and stop rules |
| [`microscopic/finite_em_coulomb_phase.md`](microscopic/finite_em_coulomb_phase.md) | Issue-#6 Stage 6A finite Hamiltonian, static Coulomb diagnostics, and rejection boundary |
| [`microscopic/check_finite_em_coulomb.py`](microscopic/check_finite_em_coulomb.py) | Exact cube, finite-size flux/transverse, Wilson-shift, and fixed-charge checks |
| [`microscopic/finite_em_coulomb_results.json`](microscopic/finite_em_coulomb_results.json) | Machine-readable Stage 6A status and remaining gates |
| [`companion_integration.md`](companion_integration.md) | Issue-#7 finite integration contract and dependency mapping |
| [`microscopic/finite_dressed_frame_hamiltonian.md`](microscopic/finite_dressed_frame_hamiltonian.md) | Issue-#2 finite Hamiltonian/transfer construction and acceptance accounting |
| [`microscopic/finite_dressed_frame_results.json`](microscopic/finite_dressed_frame_results.json) | Frozen finite logical, spectrum, lock, boundary, and symbol results |
| [`microscopic/chiral_matter_defect.md`](microscopic/chiral_matter_defect.md) | Finite chiral charged endpoint prototype and claim boundary |
| [`microscopic/`](microscopic/) | Finite constructions, executable checks, and frozen outputs |
| [`validation_protocol.md`](validation_protocol.md) | Acceptance, rejection, freezing, and data-testing rules |
| [`manifest.json`](manifest.json) | Machine-readable status summary |

## Stop rules

Do not call the equal-time rank-two electromagnetic tensor a completed photon spectrum or the winding free energy a photon gap. Do not treat one solvable RK point as a stable `z=1` phase without detuned spectral scaling. Do not treat the finite matter prototype as observed particle physics. Do not treat the finite companion event probability as an astrophysical rate, the one-particle bound mode as a halo, or closed-system coherence as a cosmological lifetime. Do not retune the finite gravity regulator independently to force a photon match. Do not call its linear success nonlinear quantum gravity or a unique ultraviolet completion. Derive the common coefficients, common cone, nonlinear closure, and regulator universality—or reject the combined branch.
