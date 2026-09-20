# Phase Junction Network

**Status:** exploratory candidate foundation. The linear gravity mode/constraint structure is internally well tested. Issue #7 now has a finite reversible photon–companion integration architecture. Finite QED, microscopic matter, nonlinear gravity, common-cone recovery, many-body deposits, and empirical predictions remain open.  
**Created:** 2026-09-20  
**Architecture audit:** [`architecture_audit.md`](architecture_audit.md)  
**Issue #7 closure:** [`companion_integration.md`](companion_integration.md)

## Working picture

The physical comparison is placed on links between neighboring matter–geometry junctions rather than in a smooth scalar gradient. Electromagnetism is compact scalar phase holonomy; gravity compares full local frames and clocks; matter is intended to arise as finite defects.

The selected companion ontology is:

| Object | Role |
|---|---|
| Photon | Transverse quantum of the compact U(1) link connection |
| Companion `chi=varphi` | Neutral relative-phase quantum `theta_m-theta_g` |
| Deposit | Bound state of the same `chi` sector |
| Graviton | Transverse frame excitation, distinct from the companion reservoir |
| Frame connection | Constrained auxiliary comparison variable |
| Matter | Finite charged or neutral defects; the complete chiral spectrum is still open |

## Established internally

### Electromagnetic kinematics

A spin-1, three-state link is the smallest tested finite representation with exact local Gauss symmetry and nonconstant electric energy. This does **not** yet demonstrate a finite 3+1-dimensional deconfined Coulomb phase or interacting QED.

### Linear gravity

The committed continuum, finite-state, exact-arithmetic, finite-size, and full-real-space checks establish under their linear assumptions:

- three vector and one scalar frame constraints;
- a local independent connection whose constrained elimination gives the Fierz–Pauli stiffness;
- 36 second-class connection constraints and four first-class frame constraints;
- four physical phase-space dimensions, or two positive linearly dispersing configuration modes, per nonzero momentum;
- no inserted transverse-traceless projector;
- a finite odd-prime Weyl connection lock with exact dressed-frame algebra;
- a reduced positive transfer-matrix target;
- linear removal of the periodic homogeneous conformal pair by fixed-volume and trace-momentum constraints.

The fixed-volume result controls one global linear mode. It is not a vacuum-energy or cosmological-constant solution.

### Finite charged matter prototype

Issue #4 has a committed finite regulator/prototype using spin-1 link bundles and bound odd-strand endpoint operators. Within the declared search class, the first primitive anomaly-free chiral spectrum is

\[
q=(-11,-5,-1,-1,9,9),
\qquad \sum q=\sum q^3=0.
\]

The same signed strand count fixes charge and fermionic exchange sign. Composite hopping changes exactly `|q|` finite flux lanes and commutes with both endpoint Gauss generators. A finite open Wilson/domain-wall slab supplies one Weyl cone on each wall, gaps the other seven physical Brillouin corners, and preserves one bare frame operator for every species. At slab width 12, a two-parameter charge-localization rule gives four distinct gaps spanning approximately `3.245e7`.

This is not the observed matter sector. The remote mirror wall, Standard Model gauge structure, generations, bound-state clocks, interacting Ward identities, observed masses, and radiative stability remain open. See [`microscopic/chiral_matter_defect.md`](microscopic/chiral_matter_defect.md).

### Finite photon–companion bridge

Issue #7 is closed at the finite architecture threshold by:

- [`companion_integration.md`](companion_integration.md), the ontology, Hamiltonian, ledgers, active-goal map, and prior-branch disposition;
- [`microscopic/check_companion_bridge.py`](microscopic/check_companion_bridge.py), the executable calculation;
- [`microscopic/companion_bridge_results.json`](microscopic/companion_bridge_results.json), the frozen output.

The calculation contains a finite source, high- and low-frequency bins of one photon field, a neutral `chi`, converter recoil, a computed bound `chi` mode, capture recoil, a receiver, and an exactly reduced constrained frame. Every transition has its Hermitian reverse.

Principal results:

```text
matter-assisted target conversion                 0.9743285451
receiver + bound chi + recoil probability         0.9742646744
energy-expectation drift                           9.77e-15
component-ledger error                             8.88e-16
forward/reverse probability difference             0
field-momentum/reciprocal-force residual           4.45e-5
frame equation residual                            2.78e-17
independent motion multipliers                      0
independent lensing multipliers                     0
```

The finite dispersion also records an unresolved result:

```text
maximum chi speed / maximum photon speed = 0.3571428571
```

A common relativistic propagation cone is therefore not established and is owned by issue #8. The large finite conversion probability is a resonant architecture stress test, not an astrophysical rate.

## What remains open

The project does not yet contain one finite local theory that simultaneously supplies:

- a demonstrated deconfined photon/Coulomb phase;
- the complete finite dressed-frame Hamiltonian;
- completion of the finite chiral prototype into the observed charged-matter spectrum and protected masses;
- shared microscopic electromagnetic, companion, matter, and gravitational coefficients;
- nonlinear gravitational self-coupling and constraint closure;
- a unitary/reflection-positive interacting relativistic continuum;
- a common causal cone for photons, companions, tensor modes, and matter;
- a microscopic vacuum/background mechanism;
- many-body capture capacity, lifetime, release, and self-gravity;
- a frozen cross-sector empirical prediction.

## Work program

| Issue | Pillar | Status |
|---:|---|---|
| #2 | Finite dressed frame–connection Hamiltonian | Open |
| #3 | Shared microscopic coefficients and `Z_g/Z_A` | Open |
| #4 | Chiral matter defects and protected mass hierarchy | Finite prototype completed; integration/continuum work delegated to #3, #6, and #8 |
| #5 | Nonlinear gravity, self-coupling, and volume/vacuum term | Open |
| #6 | Finite electromagnetic Coulomb/QED phase | Open |
| #7 | Photon–companion integration | **Finite architecture closed; physical scaling delegated** |
| #8 | Continuum quantum consistency and Lorentz/common-cone recovery | Open |
| #9 | Frozen predictions and empirical test ladder | Open |

```text
finite integration contract (#7 closed)
    -> finite EM + finite dressed gravity (#6, #2)
    -> integrate the finite matter prototype into both sectors (#4, #6, #2)
    -> shared coefficients + nonlinear closure (#3, #5)
    -> interacting continuum and common cone (#8)
    -> frozen predictions and data (#9)
```

## Relationship to the wider repository

The root project remains governed by `CURRENT-STATUS.md`, `research_plan/active-goal.md`, `research_plan/universe-contract.md`, and the twelve-item `research_plan/solution-goal-ledger.md`. It requires fixed published distances, an operationally nonexpanding universe, complete energy accounting, no inserted dark halo, and joint predictions across redshift, timing, brightness, motion, lensing, gravitational waves, and background observations.

The issue-#7 bridge supplies a finite microscopic integration architecture for those roles. It does not derive the phenomenological redshift coefficient, a continuum capture rate, a galaxy-scale deposit, or any observational fit. All twelve active goals and the R01–R32 observational program remain open except for the narrower finite architecture statements in the integration contract.

## Key files

| File | Purpose |
|---|---|
| [`architecture_audit.md`](architecture_audit.md) | Pillar audit, dependencies, and stop rules before issue-#7 closure |
| [`companion_ontology.md`](companion_ontology.md) | Ontology decision selecting `varphi`; retained as the decision record |
| [`companion_integration.md`](companion_integration.md) | Issue-#7 closure contract and complete dependency mapping |
| [`electromagnetic_derivation.md`](electromagnetic_derivation.md) | Link/loop electromagnetic target |
| [`frame_gravity_derivation.md`](frame_gravity_derivation.md) | Frame-valued gravity target and constraints |
| [`microscopic/chiral_matter_defect.md`](microscopic/chiral_matter_defect.md) | Finite chiral charged endpoint prototype and claim boundary |
| [`microscopic/`](microscopic/) | Finite constructions, executable checks, and frozen outputs |
| [`validation_protocol.md`](validation_protocol.md) | Acceptance, rejection, freezing, and data-testing rules |
| [`manifest.json`](manifest.json) | Machine-readable status summary |

## Stop rules

Do not treat the finite matter prototype as a derivation of observed particle physics. Do not treat the finite event probability as an astrophysical rate, the one-particle bound mode as a galaxy halo, closed-system coherence as a cosmological lifetime, or the synthetic frame solve as observational success. Do not impose equal photon and companion speeds after the finite benchmark found otherwise. Derive the common cone or reject the branch.
