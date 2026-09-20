# Phase Junction Network: Big-Picture Architecture Audit

**Audit date:** 2026-09-20  
**Scope:** determine whether the project is addressing every element required for a new fundamental theory, identify overdeveloped and neglected areas, and set the work order before further narrow tuning.

## Executive conclusion

The project has made unusually strong progress on one difficult slice: the **linear gravitational constraint and mode problem**. The repository now contains mutually reinforcing continuum, finite-state, real-space, exact-arithmetic, finite-size, transfer-matrix, and regression checks showing how a local frame/connection system can reduce to two positive linearly dispersing tensor modes. That work is valuable and should be preserved.

It is no longer the main uncertainty.

The theory remains incomplete at the same architectural level because five other pillars are much less developed:

1. a finite electromagnetic Coulomb/QED phase;
2. a microscopic matter and charge spectrum;
3. a unified move set that fixes electromagnetic and gravitational couplings;
4. nonlinear and quantum continuum consistency;
5. integration with the repository's active photon–companion/redshift/deposition program and its observational goals.

The correct response is not to abandon the current gravity work, but to stop treating additional linear-gravity precision as the critical path. Further gravity work should be accepted only when it closes the finite dressed Hamiltonian, nonlinear constraint algebra, universal self-coupling, or common microscopic parameter problem.

## Relationship to the wider repository

The repository's authoritative active program is governed by `CURRENT-STATUS.md`, `research_plan/active-goal.md`, and `research_plan/universe-contract.md`. It requires fixed published distances, a nonexpanding operational universe, complete energy accounting, no inserted dark-matter halo, and joint predictions for redshift, timing, brightness, motion, lensing, gravitational waves, and background observations.

`phase_junction_network/` is currently a separate candidate foundation. It supplies a possible microscopic origin for electromagnetism and gravity, but it has not yet been connected to the active program's companion excitation, emitter interaction, propagation law, capture/storage mechanism, clock response, or galaxy/cluster source model.

Until that bridge exists, progress in the Phase Junction module does not close the active program's observational or energy-supply requirements, and success in the phenomenological companion program does not establish the Phase Junction ontology.

## Master gate matrix

Status meanings:

- **Established internally:** the stated mathematical result has multiple committed checks under explicit assumptions.
- **Partial:** a target structure or kinematic construction exists, but the physical phase or full mechanism is not demonstrated.
- **Open:** no construction yet closes the requirement.

| Pillar | Minimum requirement | Current status | Evidence already present | Decisive missing work |
|---|---|---|---|---|
| Scientific contract and ontology | One declared microscopic state space and an unambiguous meaning for matter, geometry, photon, graviton, and companion | **Partial** | Link phase/imbalance, frame/connection, and neutral relative-phase candidates are distinguished in the derivations | Choose the companion identity; prevent energy double counting; specify which variables are fundamental versus effective |
| Finite electromagnetism | A finite 3+1D model with a deconfined Coulomb phase, two transverse photons, dynamical charge, and controlled continuum limit | **Partial** | Exact finite Gauss symmetry; spin 1 is the smallest tested link with nonconstant electric energy; infrared Maxwell target | Demonstrate the finite Coulomb phase, charged defects, `1/r` force, Ward identities, and QED continuum; tracked by issue #6 |
| Linear gravity | Local constraints leaving exactly two positive linearly dispersing tensor modes without a TT projector | **Established internally** | Continuum spectrum, exact coefficient identity, constrained local reduction, full real-space reduction, finite Weyl connection lock, finite-size scaling | Finish the finite dressed frame Hamiltonian at finite local dimension; tracked by issue #2 |
| Nonlinear gravity | Closed nonlinear constraints, universal self-coupling, strong-field stability, no extra scalar/ghost | **Open** | Linear first-class/second-class structure and a fixed-volume linear candidate | Nonlinear frame/connection action, self-energy source, constraint closure, strong-field spectrum; issue #5 |
| Matter | Finite charged fermionic defects with chirality, anomaly control, generations or declared reduced target, and protected mass gaps | **Open** | Gauge-covariant hopping is written only as an infrared target | Explicit finite defect model, statistics, chirality/doubling solution, mass protection; issue #4 |
| Unification and parameter closure | One microscopic move set deriving `U_A`, `K_A`, `U_g`, `K_g`, defect gaps, and fewer ratios than outputs | **Open** | Infrared formulas expose the required impedance ratio and common-speed products | Shared perturbative/direct derivation and at least one cross-sector relation; issue #3 |
| Common causal geometry | One metric/cone for photons, tensor modes, and every matter species beyond leading order | **Partial** | Shared tetrad/frame gives the desired infrared target; lattice cutoff effects are measured in the tensor branch | Demonstrate the finite electromagnetic branch on the same frame and radiative stability of common propagation; issues #6 and #8 |
| Quantum continuum consistency | Unitary or reflection-positive continuum limit, positive residues, Ward identities, regulator universality, anomaly control | **Open** | Reduced free tensor transfer matrix and reflection kernel pass; exact finite constraints exist | Full interacting continuum analysis, renormalization, Lorentz recovery, Weinberg–Witten assumption audit; issue #8 |
| Vacuum, thermodynamics, and background | Stable vacuum, controlled volume term, thermal/statistical state, and admissible nonexpanding background | **Partial at one linear global mode only** | Fixed-volume pair removes the periodic homogeneous conformal instability at linear order | Microscopic reason for the constraint, local vacuum energy, nonlinear closure, background/cosmological solutions; issue #5 |
| Photon–companion integration | One Hamiltonian from emission through propagation, capture/storage, gravity, clocks, and reverse transitions | **Open** | The wider repository has extensive phenomenological and energy-ledger work; Phase Junction has candidate excitations | Variable dictionary and end-to-end common action; issue #7 |
| Empirical distinctiveness | Fewer calibration inputs than independent outputs and at least one frozen cross-sector prediction | **Open** | Validation protocol correctly postpones data fitting; wider repository has mature data controls | Parameter ledger, weak-field dictionary, held-out predictions, staged data test; issue #9 |
| Reproducibility and negative controls | Frozen outputs, independent checks, CI, retained failures, explicit rejection gates | **Established internally** | Exact and numerical checks, result manifests, failed branches retained, GitHub Actions workflow | Keep documentation synchronized; add tests as new pillars become executable |

## Main strategic finding

The project is **not** at the stage where only small tuning remains. The linear gravity branch is close to a finite microscopic benchmark, but the complete theory is structurally unbalanced.

More decimal precision on the existing tensor spectrum would not address the largest risks:

- the finite electromagnetic model may not possess the required Coulomb phase;
- no microscopic fermion/matter spectrum exists;
- the electromagnetic and gravitational couplings remain independent;
- nonlinear interactions may reintroduce forbidden modes or nonuniversal propagation;
- the active companion/redshift mechanism has no Phase Junction identity;
- no distinct prediction currently separates the model from QED plus general relativity with chosen coefficients.

These are architecture gaps, not thoroughness gaps.

## Work order

### Priority 0: define one theory rather than parallel branches

Complete issue #7 first at the level of a **theory integration contract**, even before the full equations are solved. It must identify:

- what the companion is;
- which Phase Junction variables are fundamental;
- what carries transferred photon energy;
- what is captured or stored;
- what sources the shared frame;
- how matter clocks are represented;
- which existing phenomenological fields are effective limits and which are incompatible alternatives.

This prevents continued work from building two theories that cannot later be joined.

### Priority 0: balance the gauge and gravity sectors

Run issue #6 in parallel with the remaining finite gravity work. The electromagnetic sector must advance from exact Gauss kinematics to a demonstrated finite Coulomb/QED phase. Until then, calling the framework an electromagnetic–gravitational unification is premature.

For gravity, issue #2 should now focus narrowly on the finite dressed-frame Hamiltonian and its spectrum. Do not spend additional effort rechecking the already established continuum Fierz–Pauli kernel unless a new finite construction depends on it.

### Priority 0: construct matter before deriving precision constants

Issue #4 is required before the theory can claim QED, clocks, universal free fall, or particle masses. Matter should be built early enough that anomaly, species-doubling, and radiative-stability problems can shape the microscopic architecture rather than appear after it is frozen.

### Priority 1: derive shared coefficients and nonlinear closure

Issues #3 and #5 should converge on one common microscopic move set. The key deliverables are:

- nonlinear constraint closure;
- universal coupling to total stress-energy, including field self-energy;
- a controlled volume/vacuum mechanism;
- derived electromagnetic and gravitational stiffnesses;
- a fixed or calculated `Z_g/Z_A` rather than an equality assumption.

### Priority 1: establish the interacting continuum theory

Issue #8 tests whether the finite network actually flows to a healthy relativistic quantum theory. This is where common-cone stability, radiative corrections, regulator universality, positive residues, anomalies, and the emergent-spin-2 no-go assumptions must be addressed.

### Priority 2: freeze predictions before flexible astrophysics

Issue #9 should begin with a parameter ledger now, but observational unblinding should wait until the action and calibration roles are frozen. Local and weak-field tests precede galaxy, cluster, and cosmological flexibility.

## Dependency graph

```text
                    +-----------------------------+
                    |  ontology + integration #7 |
                    +--------------+--------------+
                                   |
             +---------------------+---------------------+
             |                     |                     |
   +---------v---------+ +---------v---------+ +---------v---------+
   | finite EM phase #6| | finite gravity #2 | | finite matter #4  |
   +---------+---------+ +---------+---------+ +---------+---------+
             |                     |                     |
             +---------------------+---------------------+
                                   |
                    +--------------v--------------+
                    | shared coefficients #3      |
                    | nonlinear closure/vacuum #5 |
                    +--------------+--------------+
                                   |
                    +--------------v--------------+
                    | quantum continuum #8        |
                    +--------------+--------------+
                                   |
                    +--------------v--------------+
                    | frozen predictions/data #9  |
                    +-----------------------------+
```

The dependencies are not entirely serial; exploratory work can proceed in parallel. But no later layer should be declared complete while its upstream inputs remain arbitrary.

## Specific integration questions that must be answered

### What is the companion?

The theory currently contains several candidates that must not be conflated:

- a photon, which is a transverse quantum of the compact link connection;
- a neutral local relative-phase mode, which is scalar and generally gapped;
- a gravitational tensor quantum, which is a frame excitation;
- a connection-relative excitation, which the finite lock intentionally gaps and removes from the physical low-energy band;
- charged or neutral defects and possible bound collective states.

The active companion cannot be called a graviton in one calculation, a scalar phase oscillation in another, and stored link imbalance in a third without an explicit conversion theory.

### How does redshift arise?

The finite electromagnetic Hamiltonian by itself propagates and quantizes light; it does not automatically transfer photon energy into another sector. A redshift mechanism requires an explicit interaction and must predict:

- frequency dependence;
- coherence and linewidth effects;
- momentum and recoil;
- reverse transitions and noise;
- duration/brightness consequences;
- source fuel and receiver energy accounting.

### How does stored energy gravitate?

The same stored excitation must enter total stress-energy and therefore the shared frame equations. A separately fitted galaxy-force field and lensing field would not satisfy the integration requirement.

### What fixes the background?

The linear fixed-volume pair controls one periodic zero-mode instability. It does not define the universe's large-scale statistical state, thermal background, boundary conditions, age, size, or origin. Those remain explicit hypotheses under the universe contract.

## Stop rules

Until the architecture gaps above are closed:

1. Do not interpret another successful linear tensor eigenvalue test as progress on matter, redshift, or cosmology.
2. Do not derive `alpha` from an arbitrary mixing-angle ansatz; derive the stiffness ratio from microscopic moves.
3. Do not call the neutral sine-Gordon mode a photon or a graviton.
4. Do not use the fixed-volume condition as a solved cosmological constant mechanism.
5. Do not fit galaxy or cosmological data with independent Phase Junction response functions before local weak-field predictions are frozen.
6. Do not insert a conventional halo, capture law, lifetime, or mass hierarchy under new names.
7. Do not count recovery of QED and general relativity after freely selecting all of their coefficients as a distinct prediction.

## What narrow work is still justified now

The following focused tasks are worth continuing because they close architecture gates rather than merely polishing them:

- build the finite dressed frame/Fierz–Pauli Hamiltonian in issue #2;
- demonstrate or rule out a spin-1 Coulomb phase in issue #6;
- write the companion/field dictionary and common energy ledger in issue #7;
- prototype the smallest charged matter defect in issue #4;
- create the machine-readable parameter ledger required by issue #9.

## Decision

The linear gravitational foundation is strong enough to stop re-litigating its basic mode count. The project should now operate as a coordinated multi-pillar program.

The next claim threshold is not “another gravity check passes.” It is:

> A single finite local model contains a deconfined photon sector, a healthy two-helicity frame sector, at least one protected charged matter defect, and a declared companion excitation; it derives shared low-energy coefficients, preserves its constraints nonlinearly, and produces at least one frozen cross-sector prediction.

Nothing currently committed meets that full threshold. The new issue structure makes every missing part explicit and prevents the most developed subproblem from being mistaken for the complete theory.