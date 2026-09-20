# Phase Junction Network

**Status:** exploratory candidate foundation with a well-tested linear gravity branch, finite electromagnetic and gravity kinematics, and a finite connection-lock construction. It is **not** yet a complete unified quantum theory, an integrated photon–companion mechanism, or an empirical fit.  
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
- matter is intended to arise as finite charged or massive defects rather than being inserted as a separate continuum field.

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

The gravity work is the most mature pillar. It now establishes, under the stated linear assumptions:

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

The unrestricted periodic zero mode contains one negative homogeneous trace-momentum direction. A fixed-total-volume plus trace-momentum second-class pair removes that direction at linear order. This is a control of one global mode, not a completed vacuum-energy or cosmological mechanism.

### Verification infrastructure

The module preserves exact arithmetic checks, numerical checks, finite-size scaling, negative branches, frozen JSON outputs, issue-level acceptance criteria, and a GitHub Actions regression workflow. See [`microscopic/README.md`](microscopic/README.md) for the current calculation inventory.

## What is not established

The project does not yet contain one finite local model that simultaneously supplies:

- a demonstrated deconfined photon/Coulomb phase;
- a complete finite dressed-frame gravity Hamiltonian;
- nonlinear gravitational constraint closure and universal self-coupling;
- finite chiral charged matter and protected particle masses;
- a microscopic calculation of `Z_g/Z_A`, `G`, and `alpha` from fewer independent inputs;
- an interacting unitary or reflection-positive relativistic continuum limit;
- a microscopic vacuum/background mechanism;
- an unambiguous Phase Junction identity for the active photon–companion excitation;
- a frozen cross-sector empirical prediction distinguishing the model from QED plus general relativity with selected coefficients.

## Big-picture work program

The architecture audit found that further linear-gravity polishing is no longer the sole critical path. The project is now organized around these top-level issues:

| Issue | Pillar |
|---:|---|
| [#2](https://github.com/lrspeiser/photon-graviton/issues/2) | Finite dressed frame–connection Hamiltonian and two tensor branches |
| [#3](https://github.com/lrspeiser/photon-graviton/issues/3) | Shared microscopic electromagnetic/gravitational coefficients and `Z_g/Z_A` |
| [#4](https://github.com/lrspeiser/photon-graviton/issues/4) | Chiral matter defects and protected mass hierarchy |
| [#5](https://github.com/lrspeiser/photon-graviton/issues/5) | Nonlinear constraints, self-coupling, and volume/vacuum term |
| [#6](https://github.com/lrspeiser/photon-graviton/issues/6) | Finite electromagnetic Coulomb and QED phase |
| [#7](https://github.com/lrspeiser/photon-graviton/issues/7) | Integration with the active photon–companion/redshift/deposition program |
| [#8](https://github.com/lrspeiser/photon-graviton/issues/8) | Continuum quantum consistency, Lorentz recovery, and regulator universality |
| [#9](https://github.com/lrspeiser/photon-graviton/issues/9) | Frozen predictions and staged empirical tests |

The immediate dependency order is:

```text
ontology/integration
    -> finite EM + finite gravity + finite matter
    -> shared coefficients + nonlinear closure
    -> interacting continuum consistency
    -> frozen predictions and data
```

Exploratory tasks may run in parallel, but later claims cannot be closed while upstream quantities remain arbitrary.

## Relationship to the wider photon–companion repository

The root project is governed by `CURRENT-STATUS.md`, `research_plan/active-goal.md`, and `research_plan/universe-contract.md`. It requires fixed published distances, a nonexpanding operational universe, complete energy accounting, no independently inserted dark halo, and joint predictions across redshift, timing, brightness, galaxy motion, lensing, gravitational waves, and background observations.

Phase Junction is presently a candidate microscopic foundation, not yet the active program's completed mechanism. Issue #7 must decide whether the companion is:

- the neutral relative-phase mode;
- a frame/tensor excitation;
- a topological or matter defect;
- a bound collective junction state;
- or a distinct field outside this framework.

Until that choice and the corresponding emitter–propagation–capture–gravity Hamiltonian are derived, the two programs must not be treated as already unified.

## Key files

| File | Purpose |
|---|---|
| [`architecture_audit.md`](architecture_audit.md) | Full pillar-by-pillar audit, dependency graph, stop rules, and work order |
| [`electromagnetic_derivation.md`](electromagnetic_derivation.md) | Corrected link/loop electromagnetic construction and QED target |
| [`frame_gravity_derivation.md`](frame_gravity_derivation.md) | Frame-valued gravity target, constraints, static limit, and parameter relations |
| [`microscopic/`](microscopic/) | Finite-state constructions, constrained reductions, exact checks, and frozen outputs |
| [`validation_protocol.md`](validation_protocol.md) | Acceptance, rejection, freezing, and data-testing rules |
| [`manifest.json`](manifest.json) | Machine-readable status and result summary |

## Current decision

Continue only narrow work that closes an architecture gate:

1. construct the finite dressed frame Hamiltonian in issue #2;
2. demonstrate or rule out the finite spin-1 Coulomb phase in issue #6;
3. write the companion/field dictionary and shared energy ledger in issue #7;
4. prototype a finite charged matter defect in issue #4;
5. begin the parameter ledger required by issue #9.

Do not interpret another successful linear tensor residual as progress on matter, redshift, cosmology, or empirical distinctiveness.