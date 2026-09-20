# Phase Junction Prior-Art and Novelty Boundary

**Date:** 2026-09-20  
**Purpose:** prevent known ingredients from being mistaken for original results, give explicit credit, and define the narrower claims that the Phase Junction program must actually establish.

## Important limitation

This is a technical literature boundary, not a legal novelty, patentability, or exhaustive-prior-art opinion. The search is necessarily incomplete. The project must continue to update this file whenever a closer construction is found.

## Established ingredients that are not Phase Junction novelty

| Ingredient used in this repository | Prior art that must be credited | Consequence for our claims |
|---|---|---|
| Finite-dimensional gauge fields on links | Chandrasekharan and Wiese, [*Quantum Link Models: A Discrete Approach to Gauge Theories*](https://arxiv.org/abs/hep-lat/9609042) | A finite spin-1 U(1) link and exact Gauss law are not by themselves new. |
| Compact plaquette/ring exchange and RK-type positive Hamiltonians | Quantum-link, quantum-dimer, and Rokhsar–Kivelson constructions | A positive graph-Laplacian point and equal-amplitude ground state are techniques, not the core novelty. |
| Emergent helicity-2 modes from finite spin/qubit models | Gu and Wen, [*Emergence of helicity ±2 modes (gravitons) from qbit models*](https://arxiv.org/abs/0907.1203) | The existence of a finite model with two tensor modes is not by itself new; our construction must differ in mechanism, closure, or predictions. |
| First-order frame/connection gravity on a discrete complex | Gionti, [*Discrete Gravity as a Local Theory of the Poincaré Group in the First Order Formalism*](https://arxiv.org/abs/gr-qc/0501082), plus Regge/Palatini/Cartan traditions | Independent frame and connection variables, torsion constraints, and area-times-curvature structure are established ideas. |
| Gauge–gravity unification through enlarged connection or Plebanski/BF structures | Alexander, Marcianò, and Tacchi, [*Towards a Loop Quantum Gravity and Yang-Mills Unification*](https://arxiv.org/abs/1105.3480); Neiman, [*Unified Lagrangians for gravity and gauge theory*](https://arxiv.org/abs/2410.18377) | Merely placing gauge and gravitational connections in one formalism is not sufficient novelty. |
| Simultaneously emergent photons and gravitons | Chkareuli, Jejelava, and Kepuladze, [*Lorentzian Goldstone modes shared among photons and gravitons*](https://arxiv.org/abs/1709.02736) and [*Emergent photons and gravitons*](https://arxiv.org/abs/1811.09578) | “Both photon and graviton emerge” is not a unique claim. The microscopic mechanism and derived relations must be different and explicit. |
| Gauge-field and loop-gravity Hilbert-space techniques | Delcamp, Dittrich, and Riello, [*Fusion basis for lattice gauge theory and loop quantum gravity*](https://arxiv.org/abs/1607.08881) | Gauge-invariant quotient bases, ribbon/flux operators, and coarse-graining language have prior art. |
| Relative-phase collective modes | Leggett-mode literature, e.g. Sharapov, Gusynin, and Beck, [*Effective action approach to the Leggett's mode in two-band superconductors*](https://arxiv.org/abs/cond-mat/0205131) | A neutral relative-phase scalar is a known type of collective excitation. Calling it a companion does not make the mode itself new. |
| Domain-wall chiral fermions | Kaplan/domain-wall-fermion literature, including Shamir, [*The Euclidean Spectrum of Kaplan's Lattice Chiral Fermions*](https://arxiv.org/abs/hep-lat/9212010) | A fifth-direction wall, mirror wall, and exponentially small overlap gap are regulator techniques, not a new matter mechanism by themselves. |
| Perturbatively generated superexchange/ring exchange | Schrieffer–Wolff and degenerate-perturbation methods | Deriving a fourth-order plaquette term from virtual local moves is standard methodology. |

## Claims we must not make

The repository should not claim originality for any of the following in isolation:

- a compact phase field;
- a Josephson or sine-Gordon relative-phase mode;
- a finite U(1) quantum link;
- exact Gauss-law stabilizers;
- a Rokhsar–Kivelson point;
- domain-wall chiral fermions;
- a tetrad/frame connection;
- two linearized helicity-2 modes;
- an emergent photon and graviton in one broad framework;
- a Plebanski, BF, Palatini, Regge, spin-foam, or enlarged-gauge-group rewriting;
- a Schrieffer–Wolff derivation of ring exchange;
- recovery of Maxwell theory or linearized Einstein gravity after their coefficients are freely selected.

## Candidate original contribution

The potentially original scientific content is the **specific combined mechanism**, if it can be completed:

1. one finite matter–geometry junction Hilbert space contains both compact scalar phase comparison and frame comparison;
2. one local constraint structure produces the electromagnetic Gauss sector and the frame scalar/vector constraints without inserting continuum gauge fields;
3. one elementary junction move and one spectrum of virtual penalties generate both electromagnetic and gravitational stiffnesses;
4. the same derivation fixes, rather than separately selects, the photon/gravity impedance ratio and causal normalization;
5. a neutral relative-phase excitation of that same junction supplies the declared companion, with reversible production, propagation, binding, release, and stress-energy;
6. finite matter defects use the same link and frame variables for charge, clocks, recoil, free fall, and source energy;
7. one sourced frame predicts both massive-body motion and light response;
8. the completed theory produces at least one cross-sector prediction not already present in QED plus general relativity with independently chosen coefficients;
9. the same microscopic theory has a controlled nonexpanding-background and many-body limit consistent with the repository's universe contract.

No single current file establishes all nine points.

## Nearest-neighbor distinction

The Phase Junction route is not intended to be:

- a quantum-link reformulation of ordinary QED alone;
- a Gu–Wen-type emergent-graviton spin model alone;
- a Plebanski or spin-foam unification obtained by enlarging a continuum gauge group;
- a spontaneous-Lorentz-breaking photon/graviton Goldstone model;
- a Kaluza–Klein or string compactification;
- a conventional scalar–tensor theory with an added dark field;
- a relabeling of a Leggett mode as dark matter;
- a domain-wall-fermion model with gravity appended afterward.

The distinguishing test is whether the **same finite microscopic moves** calculate the relative electromagnetic, gravitational, matter, and companion responses with fewer inputs than outputs.

## Originality protocol for every new stage

Every proposed “new” mechanism must include:

1. the nearest known constructions and explicit citations;
2. a sentence identifying which part is borrowed method and which result is project-specific;
3. an independent derivation or executable calculation rather than copied prose or equations;
4. at least one negative control showing that the claimed relation disappears when the shared microscopic assumption is removed;
5. a parameter count before and after the mechanism;
6. a claim boundary separating architecture, continuum physics, and empirical success;
7. preservation of failed or superseded branches;
8. a search note recording the terms and literature classes checked.

## Current originality bottleneck

The project has credible finite implementations of the separate electromagnetic, linear-gravity, matter, and companion interfaces. Their separate existence is not the novel result.

The originality bottleneck is issue #3:

> Can one elementary matter–geometry junction move and one virtual-state spectrum derive `U_A`, `K_A`, `U_g`, `K_g`, matter localization, and companion couplings—leaving fewer microscopic ratios than independent observables?

The first committed calculation addressing that question is [`microscopic/shared_junction_move_stage3a.md`](microscopic/shared_junction_move_stage3a.md). It is deliberately labeled a prototype, not proof of formal novelty or completion of the theory.

## Search terms used for this audit

The limited audit included combinations of:

- quantum link model finite-dimensional gauge theory;
- emergent helicity-2 qubit/spin model;
- emergent photon and graviton same model;
- first-order discrete Palatini/tetrad/connection gravity;
- Plebanski/BF gravity Yang–Mills unification;
- Berry-connection emergent electromagnetism and spacetime;
- relative-phase/Leggett mode;
- domain-wall chiral fermions;
- lattice gravity and gauge-theory shared Hilbert space;
- finite shared microscopic photon–graviton Hamiltonian.

The audit did not find an exact match to the complete nine-part candidate contribution above. That is evidence for continued investigation, not proof that no such work exists.