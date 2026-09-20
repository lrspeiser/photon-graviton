# Phase Junction Network

**Status:** exploratory candidate foundation with a well-tested linear gravity branch, finite electromagnetic and gravity kinematics, a finite connection-lock construction, and a first finite chiral charged-matter prototype. It is **not** yet a complete unified quantum theory, a demonstrated finite QED phase, an integrated photon–companion mechanism, or an empirical fit.  
**Created:** 2026-09-20  
**Big-picture review:** [`architecture_audit.md`](architecture_audit.md)

## Core idea

A local matter–geometry phase difference is not itself light. For a smooth scalar phase,

\[
A_\mu=\frac{\hbar}{q}\partial_\mu\phi
\quad\Rightarrow\quad
F_{\mu\nu}=0,
\]

so the original scalar-gradient proposal is pure gauge. The corrected model moves the physical comparison to **links between neighboring junctions**:

- each link carries a compact phase comparison and conjugate matter–geometry imbalance;
- electromagnetism is the non-closing scalar phase accumulated around loops;
- photons are intended to be quantized transverse normal modes of that constrained link network;
- gravity compares full local frames and clocks, requiring a richer connection and four constraints;
- matter appears as finite charged endpoint defects rather than as a continuum field inserted after the fact;
- the provisional companion is a neutral relative-phase quantum, distinct from photons, gravitons, and charged matter.

The working statement is:

> Reality is modeled as a network of microscopic matter–geometry junctions. Electromagnetism is scalar phase holonomy. Gravity is frame-and-clock holonomy. Matter appears as network defects and endpoints. Local constraints enforce charge and stress-energy conservation.

## What is established internally

### Electromagnetic structure

The infrared target is

\[
H_A=\frac{U_A}{2}\sum_\ell E_\ell^2-K_A\sum_p\cos B_p,
\]

with

\[
c_\gamma=\frac{\ell\sqrt{U_AK_A}}{\hbar},
\qquad
\alpha_{\rm bare}=\frac{1}{4\pi}\sqrt{\frac{U_A}{K_A}}.
\]

Committed checks establish exact finite Gauss symmetry and identify a spin-1, three-state link as the smallest tested representation with nonconstant electric energy. They do **not** yet establish that the complete finite 3+1-dimensional spin-1 model lies in a deconfined Coulomb phase or reproduces interacting QED.

### Linear gravity structure

The gravity work remains the most mature pillar. It establishes, under the stated linear assumptions:

- three vector constraints and one scalar constraint;
- a local independent frame connection whose constrained elimination gives the Fierz–Pauli stiffness;
- 36 second-class connection constraints and four first-class frame constraints;
- exactly four physical phase-space dimensions, or two configuration modes, per nonzero momentum;
- two equal positive frequencies proportional to lattice momentum;
- the same result from full real-space reduction without a transverse-traceless projector;
- a synthetic long-distance `1/r` Green function;
- finite-size tensor gap closing and vanishing lattice anisotropy;
- a reduced positive transfer-matrix target;
- a finite odd-prime Weyl lock implementing the connection second-class pair with an exact dressed frame algebra.

The unrestricted periodic zero mode contains one negative homogeneous trace-momentum direction. A fixed-total-volume plus trace-momentum second-class pair removes that direction at linear order. This controls one global mode; it is not a completed vacuum-energy or cosmological mechanism.

### Finite charged matter prototype

Issue #4 now has a committed finite regulator/prototype rather than only an infrared hopping ansatz. The construction uses finite spin-1 link bundles and bound odd-strand endpoint operators. Within the declared search class, the first primitive anomaly-free chiral charge spectrum is

\[
q=(-11,-5,-1,-1,9,9),
\qquad
\sum q=\sum q^3=0.
\]

The same signed strand count fixes both the defect charge and its fermionic exchange sign. Composite hopping changes exactly `|q|` finite electric-flux lanes and commutes with both endpoint Gauss generators. A finite open Wilson/domain-wall slab supplies one Weyl cone on each wall, gaps the other seven physical Brillouin corners, preserves an identical bare frame operator for every species, and yields exponentially protected residual gaps. At slab width 12, one two-parameter charge-localization rule gives four distinct gaps spanning a factor of approximately `3.245e7`.

This is a finite chiral endpoint prototype, not the observed matter sector. The remote mirror wall remains explicit; the Standard Model gauge group, color, weak isospin, generations, observed masses, bound-state clocks, interacting Ward identities, and radiative stability remain open. See [`microscopic/chiral_matter_defect.md`](microscopic/chiral_matter_defect.md).

### Verification infrastructure

The module preserves exact arithmetic checks, numerical checks, finite-size scaling, negative branches, frozen JSON outputs, issue-level acceptance criteria, and a GitHub Actions regression workflow. See [`microscopic/README.md`](microscopic/README.md) for the calculation inventory.

## What is not established

The project does not yet contain one finite local model that simultaneously supplies:

- a demonstrated deconfined photon/Coulomb phase;
- a complete finite dressed-frame gravity Hamiltonian;
- nonlinear gravitational constraint closure and universal self-coupling;
- an observed matter spectrum, mirror-sector completion, or derived particle-scale parameters;
- a microscopic calculation of `Z_g/Z_A`, `G`, `alpha`, and the matter localization parameters from fewer independent inputs;
- an interacting unitary or reflection-positive relativistic continuum limit;
- a microscopic vacuum/background mechanism;
- a completed source/propagation/capture Hamiltonian for the provisional companion;
- a frozen cross-sector empirical prediction distinguishing the model from QED plus general relativity with selected coefficients.

## Big-picture work program

The architecture audit found that further linear-gravity polishing is no longer the sole critical path. The project is organized around these top-level issues:

| Issue | Pillar | Current role |
|---:|---|---|
| [#2](https://github.com/lrspeiser/photon-graviton/issues/2) | Finite dressed frame–connection Hamiltonian and two tensor branches | Open |
| [#3](https://github.com/lrspeiser/photon-graviton/issues/3) | Shared microscopic coefficients and `Z_g/Z_A` | Open; now also owns derivation of matter localization parameters |
| [#4](https://github.com/lrspeiser/photon-graviton/issues/4) | Chiral matter defects and protected mass hierarchy | Finite prototype completed; downstream interaction/continuum work handed to #3, #6, and #8 |
| [#5](https://github.com/lrspeiser/photon-graviton/issues/5) | Nonlinear constraints, self-coupling, and volume/vacuum term | Open |
| [#6](https://github.com/lrspeiser/photon-graviton/issues/6) | Finite electromagnetic Coulomb and QED phase | Open; must include the new dynamical defects |
| [#7](https://github.com/lrspeiser/photon-graviton/issues/7) | Photon–companion integration | Provisional identity selected; common Hamiltonian open |
| [#8](https://github.com/lrspeiser/photon-graviton/issues/8) | Continuum consistency, Lorentz recovery, and regulator universality | Open; owns mirror completion and radiative stability |
| [#9](https://github.com/lrspeiser/photon-graviton/issues/9) | Frozen predictions and staged empirical tests | Open |

The immediate dependency order is:

```text
provisional ontology/integration contract
    -> finite EM phase + finite dressed gravity
    -> new finite matter prototype integrated into both sectors
    -> shared coefficients + nonlinear closure
    -> interacting continuum consistency
    -> frozen predictions and data
```

Exploratory tasks may run in parallel, but later claims cannot be closed while upstream quantities remain arbitrary.

## Relationship to the wider photon–companion repository

The root project is governed by `CURRENT-STATUS.md`, `research_plan/active-goal.md`, and `research_plan/universe-contract.md`. It requires fixed published distances, a nonexpanding operational universe, complete energy accounting, no independently inserted dark halo, and joint predictions across redshift, timing, brightness, galaxy motion, lensing, gravitational waves, and background observations.

The provisional ontology is now explicit:

- photon: transverse quantum of the compact U(1) link connection;
- graviton: transverse frame excitation;
- companion: neutral relative-phase quantum `varphi`;
- deposit: bound collective configuration of the same `varphi` sector;
- matter: finite charged or neutral defects, distinct from the companion.

This decision removes a naming ambiguity, but issue #7 still must derive one reciprocal Hamiltonian covering production, propagation, capture/storage, reverse transitions, recoil, stress-energy sourcing, clocks, motion, and lensing. Until then, the Phase Junction and active phenomenological programs are not a completed unified theory.

## Key files

| File | Purpose |
|---|---|
| [`architecture_audit.md`](architecture_audit.md) | Full pillar-by-pillar audit, dependency graph, stop rules, and work order |
| [`companion_ontology.md`](companion_ontology.md) | Provisional companion identity, field dictionary, energy ledger, and falsification order |
| [`electromagnetic_derivation.md`](electromagnetic_derivation.md) | Corrected link/loop electromagnetic construction and QED target |
| [`frame_gravity_derivation.md`](frame_gravity_derivation.md) | Frame-valued gravity target, constraints, static limit, and parameter relations |
| [`microscopic/chiral_matter_defect.md`](microscopic/chiral_matter_defect.md) | Finite chiral charged endpoint, anomaly, hierarchy, and claim boundary |
| [`microscopic/`](microscopic/) | Finite-state constructions, constrained reductions, exact checks, and frozen outputs |
| [`validation_protocol.md`](validation_protocol.md) | Acceptance, rejection, freezing, and data-testing rules |
| [`manifest.json`](manifest.json) | Machine-readable status and result summary |

## Current decision

Continue only narrow work that closes an architecture gate:

1. build the finite matter-assisted photon–`varphi` transition in issue #7;
2. demonstrate or rule out the finite spin-1 Coulomb phase with the new dynamical defects in issue #6;
3. construct the finite dressed frame Hamiltonian in issue #2;
4. derive the electromagnetic, gravitational, and matter localization coefficients from a common move set in issue #3;
5. test mirror completion, interacting anomalies, Ward identities, and Lorentz recovery in issue #8;
6. continue the parameter ledger required by issue #9.

Do not interpret the finite matter prototype as a derivation of observed particle physics, and do not interpret another successful linear tensor residual as progress on redshift, cosmology, or empirical distinctiveness.
