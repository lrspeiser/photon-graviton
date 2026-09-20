# Phase Junction Network: Big-Picture Architecture Audit

**Audit date:** 2026-09-20  
**Issue #4 update:** finite chiral charged-endpoint prototype committed on 2026-09-20  
**Scope:** determine whether the project addresses every element required for a new fundamental theory, identify overdeveloped and neglected areas, and set the work order before further narrow tuning.

## Executive conclusion

The project has made unusually strong progress on the **linear gravitational constraint and mode problem**. Mutually reinforcing continuum, finite-state, real-space, exact-arithmetic, finite-size, transfer-matrix, and regression checks show how a local frame/connection system can reduce to two positive linearly dispersing tensor modes. That work should be preserved, but it is no longer the main uncertainty.

The first microscopic matter placeholder has now also been replaced. Issue #4 supplies a finite charged endpoint, exact Gauss-covariant hopping, a finite anomaly-free chiral charge spectrum, one universal bare frame operator, and an exponentially protected finite-width gap hierarchy. This closes the first-defect/regulator scope of the matter pillar. It does **not** derive the Standard Model, remove the regulator mirror, establish interacting QED, or derive the matter parameters from the common microscopic move set.

The full theory remains architecturally incomplete because the following pillars are still decisive:

1. a finite electromagnetic deconfined Coulomb/QED phase containing the new dynamical defects;
2. a complete finite dressed-frame gravity Hamiltonian;
3. a unified move set fixing electromagnetic, gravitational, companion, and matter coefficients;
4. nonlinear gravitational closure and a microscopic volume/vacuum mechanism;
5. an interacting continuum with mirror completion, Ward identities, unitarity, Lorentz recovery, and radiative stability;
6. one completed photon–companion production/propagation/capture Hamiltonian;
7. a frozen cross-sector prediction and staged data program.

The correct response is not to return to broad linear-gravity polishing or to identify the new matter toy spectrum with observed particles. The next work must integrate the now-explicit finite sectors and close the remaining architecture gates.

## Relationship to the wider repository

The repository's active program is governed by `CURRENT-STATUS.md`, `research_plan/active-goal.md`, and `research_plan/universe-contract.md`. It requires fixed published distances, a nonexpanding operational universe, complete energy accounting, no inserted dark-matter halo, and joint predictions for redshift, timing, brightness, motion, lensing, gravitational waves, and background observations.

`phase_junction_network/` is a candidate microscopic foundation. Its provisional ontology is now explicit:

- photon: transverse quantum of the compact U(1) link connection;
- graviton: transverse frame excitation;
- companion: neutral relative-phase quantum `varphi`;
- deposit: bound collective configuration of the same `varphi` sector;
- matter: finite charged or neutral endpoint defects;
- frame connection: constrained auxiliary comparison variable, not an additional low-energy particle.

The matter endpoint is not the companion. Stored companion energy must source the same frame as every other energy contribution, and the same geometry must predict massive-body motion and light bending. Until issue #7 derives the reciprocal source/propagation/capture dynamics, progress in the Phase Junction module does not close the active program's observational or energy-supply requirements.

## Master gate matrix

Status meanings:

- **Established internally:** the stated mathematical result has multiple committed checks under explicit assumptions.
- **Prototype pass:** an explicit finite construction passes its declared algebraic or free-regulator gates, while interaction or continuum completion remains open.
- **Partial:** a target structure or kinematic construction exists, but the physical phase or full mechanism is not demonstrated.
- **Open:** no construction yet closes the requirement.

| Pillar | Minimum requirement | Current status | Evidence already present | Decisive missing work |
|---|---|---|---|---|
| Scientific contract and ontology | One declared microscopic state space and unambiguous meanings for matter, geometry, photon, graviton, and companion | **Partial** | Photon, graviton, `varphi` companion, deposit, connection, and finite matter roles are separated in `companion_ontology.md` | Derive one common Hamiltonian and prevent all energy double counting; issue #7 |
| Finite electromagnetism | A finite 3+1D model with a deconfined Coulomb phase, two transverse photons, dynamical charge, and controlled continuum limit | **Partial** | Exact finite Gauss symmetry; spin 1 is the smallest tested link with nonconstant electric energy; finite gauge-covariant defects now exist | Demonstrate deconfinement, `1/r`, two transverse branches, Ward identities, and scaling with dynamical defects; issue #6 |
| Linear gravity | Local constraints leaving exactly two positive linearly dispersing tensor modes without a TT projector | **Established internally** | Continuum spectrum, exact coefficient identity, constrained local reduction, full real-space reduction, finite Weyl connection lock, finite-size scaling | Finish the finite dressed-frame Hamiltonian at finite local dimension; issue #2 |
| Nonlinear gravity | Closed nonlinear constraints, universal self-coupling, strong-field stability, no extra scalar/ghost | **Open** | Linear first-/second-class structure and fixed-volume linear candidate | Nonlinear frame/connection action, self-energy source, constraint closure, strong-field spectrum; issue #5 |
| Matter | Finite charged fermionic defects with chirality, anomaly control, declared doubling treatment, and protected gaps | **Prototype pass** | Odd-strand endpoint, exact flux-bundle hopping, anomaly-free `(-11,-5,-1,-1,9,9)`, one Weyl cone per wall, protected overlap gaps, universal frame derivative | Derive observed content or a declared reduced target; remove/symmetrically gap the mirror; add interactions, clocks, and bound states; issues #3, #6, #8 |
| Unification and parameter closure | One move set deriving `U_A`, `K_A`, `U_g`, `K_g`, companion couplings, matter localization, and fewer ratios than outputs | **Open** | Infrared formulas expose the required impedance ratio; matter now gives a two-parameter charge-to-gap relation | Derive all coefficients from common amplitudes and produce a held-out relation; issue #3 |
| Common causal geometry | One metric/cone for photons, tensor modes, and every matter species beyond leading order | **Partial** | Shared tetrad/frame target; tensor cutoff effects measured; all finite matter species use one bare frame derivative | Demonstrate the finite EM phase on the same frame and interacting common-cone stability; issues #6 and #8 |
| Quantum continuum consistency | Unitary or reflection-positive continuum limit, positive residues, Ward identities, regulator universality, anomaly control | **Open** | Reduced free tensor transfer matrix; finite free matter anomaly cancellation and topology checks | Full interacting continuum, mirror completion, renormalization, Lorentz recovery, regulator comparison, no-go audits; issue #8 |
| Vacuum, thermodynamics, and background | Stable vacuum, controlled volume term, thermal/statistical state, admissible nonexpanding background | **Partial at one linear global mode only** | Fixed-volume pair removes the periodic homogeneous conformal instability at linear order | Microscopic reason, local vacuum energy, nonlinear closure, background solutions; issue #5 |
| Photon–companion integration | One Hamiltonian from emission through propagation, capture/storage, gravity, clocks, and reverse transitions | **Partial ontology only** | Provisional `varphi` identity and energy-ledger rules; finite matter participant now exists | Finite reversible matter-assisted photon–`varphi` transition, propagation, capture, stress-energy, clocks, and lensing; issue #7 |
| Empirical distinctiveness | Fewer calibration inputs than independent outputs and at least one frozen cross-sector prediction | **Open** | Validation protocol postpones fitting; finite matter supplies a candidate internal charge/gap relation | Parameter ledger, weak-field dictionary, held-out prediction, staged test; issue #9 |
| Reproducibility and negative controls | Frozen outputs, independent checks, CI, retained failures, explicit rejection gates | **Established internally** | Exact/numerical checks, result manifests, failed branches, GitHub Actions | Keep status files synchronized; add tests as new pillars become executable |

## Issue #4 resolution: what changed

The matter pillar previously contained only an infrared expression of the form

\[
\psi_x^\dagger e^{iQa_{xy}}\psi_y,
\]

which assumed the field, its charge, its statistics, and its mass. The new finite construction replaces that placeholder with:

1. eleven parallel spin-1 quantum-link lanes per geometric link;
2. a charge-`q` endpoint built from `|q|` microscopic fermionic strands;
3. an exact equal-occupation binding band with sixteen composite Fock states per site/layer;
4. a composite hopping operator that moves exactly `|q|` finite flux lanes and has zero endpoint Gauss residual on every nonzero bundle transition;
5. odd charges, so the same endpoint multiplicity gives fermionic exchange `(-1)^(q^2)=-1`;
6. an exhaustive finite anomaly search selecting `(-11,-5,-1,-1,9,9)` as the first primitive odd chiral spectrum through six species and `|q|<=11`;
7. a finite open Wilson/domain-wall slab with one light corner per wall, opposite wall chirality, seven gapped physical corners, and particle/antiparticle spectral pairing;
8. one identical frame derivative for every species;
9. a common relation `r_q=r_0+eta(|q|-1)` with finite-wall gap `Delta_q ~ r_q^L_s`, producing four gaps from two shared parameters rather than one tiny onsite energy per species;
10. finite-neighborhood tests showing the one-cone phase survives twenty nearby Wilson-mass samples.

The frozen `L_s=12` hierarchy is approximately `3.245e7`. This is a protected regulator hierarchy, not an observed mass fit. The two shared localization parameters remain inputs until issue #3 derives them from the common move set.

The mirror wall is explicit. That satisfies the requirement to declare and test the doubling treatment, but it does not complete a strictly 3+1-dimensional interacting chiral gauge theory. Mirror decoupling or symmetric gapping remains a hard gate in issue #8.

## Main strategic finding

The project is **not** at a stage where only small tuning remains. Two foundations are now relatively concrete—linear gravity and a first finite matter endpoint—but the complete theory is still structurally unbalanced.

More decimal precision on the tensor spectrum or more variants of the free matter slab would not address the largest risks:

- the finite electromagnetic model may not possess a deconfined Coulomb phase once dynamical defects are included;
- the finite dressed gravity Hamiltonian is not finished;
- the matter mirror may fail to decouple in a healthy interacting continuum;
- the electromagnetic, gravitational, companion, and matter coefficients remain independently adjustable;
- nonlinear interactions may reintroduce forbidden modes or nonuniversal propagation;
- the provisional companion has no completed production/capture mechanism;
- no distinct prediction currently separates the model from QED plus general relativity with selected coefficients.

These are architecture gaps, not residual-tolerance gaps.

## Work order

### Priority 0: complete the integration dynamics

Issue #7 now has both a provisional companion identity and a finite charged matter participant. The next calculation is the smallest reversible matter-assisted photon–`varphi` transition. It must conserve charge and total energy, include recoil, contain the reverse channel, and keep photon loss, companion energy, binding energy, and gravitational sourcing in one ledger.

### Priority 0: establish the finite electromagnetic phase

Issue #6 must advance from exact Gauss kinematics to a demonstrated finite Coulomb/QED phase using the actual spin-1 links and the new dynamical endpoints. Required outputs include a phase diagram, two transverse branches, expanding-range `1/r`, finite-representation scaling, charge propagation, and Ward identities.

### Priority 0: finish the finite dressed gravity Hamiltonian

Issue #2 should remain narrowly focused on expressing frame dynamics in the exact dressed finite operators and recovering the two tensor branches at finite local dimension. Additional continuum Fierz–Pauli checks count only if required by that finite construction.

### Priority 1: derive shared coefficients

Issue #3 must place the matter relation inside the same microscopic move set as the electromagnetic and gravitational stiffnesses. It now owns at least

\[
U_A, K_A, U_g, K_g, r_0, \eta,
\]

plus matter hopping and companion couplings. The finite charge/gap relation becomes predictive only when those parameters are derived or their calibration roles are frozen.

### Priority 1: nonlinear and interacting continuum closure

Issue #5 owns nonlinear constraints, self-energy, strong-field behavior, and the volume/vacuum term. Issue #8 owns unitarity/reflection positivity, positive residues, common-cone recovery, regulator universality, mirror completion, Ward identities, interacting anomaly accounting, radiative stability, and the precise assumptions behind emergent-spin-2 no-go results.

### Priority 2: freeze predictions before flexible astrophysics

Issue #9 should maintain the parameter ledger now, but observational unblinding must wait until the common action and calibration roles are frozen. Local and weak-field tests precede galaxy, cluster, and cosmological flexibility.

## Dependency graph

```text
                 +--------------------------------+
                 | ontology + integration #7      |
                 | provisional identity selected  |
                 +---------------+----------------+
                                 |
          +----------------------+----------------------+
          |                      |                      |
+---------v---------+  +---------v---------+  +---------v---------+
| finite EM phase #6|  | finite gravity #2 |  | matter prototype  |
| open              |  | open              |  | issue #4 passed   |
+---------+---------+  +---------+---------+  +---------+---------+
          |                      |                      |
          +----------------------+----------------------+
                                 |
                 +---------------v---------------+
                 | shared coefficients #3        |
                 | nonlinear closure/vacuum #5   |
                 +---------------+---------------+
                                 |
                 +---------------v---------------+
                 | interacting continuum #8      |
                 | includes mirror completion     |
                 +---------------+---------------+
                                 |
                 +---------------v---------------+
                 | frozen predictions/data #9    |
                 +-------------------------------+
```

The dependencies are not entirely serial; exploratory work can proceed in parallel. No later layer should be declared complete while its upstream inputs remain arbitrary.

## Specific integration questions that remain

### How is the companion produced?

A bare photon–scalar vertex does not by itself establish free single-photon conversion. The finite construction must identify the additional participant—ordinary matter recoil, a background field, an inhomogeneous environment, a multiparticle process, or another explicit degree of freedom—and include the reverse transition.

### How does stored energy gravitate?

The same stored `varphi` excitation must enter total stress-energy and therefore the shared frame equations. A separately fitted galaxy-force field and lensing field would fail the integration requirement.

### What makes a material clock?

The finite endpoints supply charged one-particle states, but not atoms, chemistry, nuclei, or bound-state clocks. Clock response cannot be asserted from the one-particle frame derivative alone. It requires interactions and bound states built from the same defects.

### What removes the mirror?

The domain-wall regulator has an opposite-chirality remote wall. A complete theory must remove, symmetrically gap, or demonstrably decouple it while preserving gauge symmetry, anomaly cancellation, locality at the intended scale, and a healthy continuum.

### What fixes the background?

The linear fixed-volume pair controls one periodic zero-mode instability. It does not define the universe's large-scale statistical state, thermal background, boundary conditions, age, size, or origin. Those remain explicit hypotheses under the universe contract.

## Stop rules

Until the remaining architecture gates close:

1. Do not interpret another successful linear tensor eigenvalue test as progress on matter, redshift, or cosmology.
2. Do not identify the issue #4 charge set with Standard Model particles without deriving a mapping and the additional gauge structure.
3. Do not treat the two matter localization parameters as predicted constants until issue #3 derives them or issue #9 records them as calibration inputs.
4. Do not hide or silently discard the remote mirror wall.
5. Do not infer a deconfined QED phase merely because exact Gauss-covariant defects exist.
6. Do not derive `alpha` from an arbitrary mixing angle; derive the stiffness ratio from microscopic moves.
7. Do not call the neutral `varphi` mode a photon, graviton, or charged matter defect.
8. Do not use the fixed-volume condition as a solved cosmological-constant mechanism.
9. Do not fit galaxy or cosmological data before local weak-field predictions and parameter roles are frozen.
10. Do not count recovery of QED and general relativity after freely selecting all coefficients as a distinct prediction.

## What narrow work is justified now

The following tasks close architecture gates rather than merely polishing them:

- build the finite reversible matter–photon–`varphi` transition in issue #7;
- demonstrate or rule out a spin-1 Coulomb phase with dynamical defects in issue #6;
- build the finite dressed frame/Fierz–Pauli Hamiltonian in issue #2;
- derive shared gauge, gravity, companion, and matter coefficients in issue #3;
- close nonlinear constraints and the volume/vacuum mechanism in issue #5;
- solve mirror completion and interacting continuum consistency in issue #8;
- maintain the machine-readable parameter ledger required by issue #9.

## Decision

The linear gravitational foundation is strong enough to stop re-litigating its basic mode count. The finite matter endpoint is explicit enough to close issue #4 at its declared first-defect/regulator scope and to become an input to the electromagnetic, integration, shared-parameter, and continuum programs.

The next claim threshold is:

> A single finite local model contains a deconfined two-polarization photon sector, a finite healthy two-helicity frame sector, the protected charged matter endpoint without an unwanted mirror, and the declared `varphi` companion; it derives shared low-energy coefficients, preserves constraints nonlinearly, and produces at least one frozen cross-sector prediction.

Nothing currently committed meets that full threshold. The updated issue structure makes the remaining dependencies explicit and prevents either the strongest subproblem or the new matter prototype from being mistaken for the complete theory.
