# Phase Junction Network Validation Protocol

**Status:** linear gravitational structure, finite electromagnetic kinematics, and a finite chiral charged-endpoint prototype are internally checked. The deconfined electromagnetic phase, finite dressed gravity Hamiltonian, nonlinear theory, shared coefficients, interacting continuum, companion integration, clocks, and empirical gates remain open.  
**Date:** 2026-09-20  
**Architecture review:** [`architecture_audit.md`](architecture_audit.md)

## 1. Research rule

No single successful spectrum, fit, or regulator implementation verifies the theory. The required sequence is:

\[
\text{ontology and finite variables}
\rightarrow
\text{physical phases and constraints}
\rightarrow
\text{matter and shared coefficients}
\rightarrow
\text{nonlinear/quantum continuum closure}
\rightarrow
\text{frozen predictions}
\rightarrow
\text{real-data tests}.
\]

Real measurements are initially exclusion walls, not flexible targets used to choose the action. A claim advances only when its upstream dependencies are declared and its free parameters are counted.

## 2. Current gate status

| Gate | Status | Evidence or issue |
|---|---|---|
| Exact finite electromagnetic Gauss symmetry | **Pass kinematically** | Spin-S quantum-link checks |
| Smallest tested link with nonconstant electric energy | **Pass at spin 1** | Spin-1/2 has constant `E^2`; spin 1 is nontrivial |
| Finite gauge-covariant charged endpoint | **Pass at prototype level** | Odd-strand composite and exact flux-bundle hopping; issue #4 evidence |
| Finite deconfined Coulomb/QED phase | **Open** | Issue #6 |
| Four finite gravitational constraints and local two-mode count | **Pass kinematically** | Odd-prime Weyl/CSS skeleton |
| Local frame/connection reduction to two positive linear modes | **Pass at the linear level** | Constrained local and full real-space reductions; no TT projector |
| Finite connection second-class lock | **Pass for the linear map** | Relative Weyl oscillator and exact dressed frame algebra |
| Finite dressed gravity Hamiltonian at finite local dimension | **Open** | Issue #2 |
| Nonlinear gravitational constraint closure and self-coupling | **Open** | Issue #5 |
| Linear periodic homogeneous conformal direction | **Controlled by a candidate** | Fixed-volume and trace-momentum second-class pair |
| Microscopic vacuum/volume mechanism | **Open** | Issue #5 |
| Finite chiral matter and protected gaps | **Pass at finite free-regulator prototype scope** | Issue #4; mirror/interacting completion handed to #3, #6, and #8 |
| Observed matter spectrum, clocks, and bound states | **Open** | Issues #3, #6, #8, and #9 |
| Shared microscopic `U_A,K_A,U_g,K_g` and matter localization parameters | **Open** | Issue #3 |
| Companion identity and bridge to active redshift/deposition program | **Partial** | Provisional `varphi` identity selected; common Hamiltonian open in issue #7 |
| Interacting continuum, unitarity, Lorentz recovery, anomalies | **Open** | Issue #8 |
| Frozen cross-sector predictions and empirical ladder | **Open** | Issue #9 |
| Reproducible regression suite | **Pass for current checks** | GitHub Actions and frozen outputs |

## 3. Ontology and integration gate

Before claiming a unified theory, retain one state/field dictionary containing:

- compact electromagnetic link phase and electric imbalance;
- frame/tetrad variables and connection variables;
- neutral relative-phase mode `varphi`;
- finite charged and neutral matter defects;
- the active companion excitation and any bound deposit;
- global or boundary variables such as volume constraints.

The provisional identity is fixed for present derivations: the companion is the neutral `varphi` quantum, a deposit is a bound collective `varphi` configuration, photons are transverse U(1) link quanta, gravitons are transverse frame excitations, and matter defects are distinct objects.

Every physical energy contribution must appear exactly once. Excitations may transform between sectors only through explicit reciprocal interactions. Issue #7 passes only when one common action or Hamiltonian covers emission, propagation, capture/storage, reverse transitions, stress-energy sourcing, clocks, motion, and lensing.

## 4. Electromagnetic gates

A finite electromagnetic candidate must demonstrate:

1. exact or controlled local Gauss symmetry;
2. a stable deconfined phase in 3+1 dimensions;
3. two and only two massless transverse branches;
4. positive residues and energy;
5. static charge energy approaching `1/r` over an expanding physical range;
6. gauge-covariant dynamical charged defects;
7. finite-size and finite-representation scaling;
8. Ward identities and a declared QED continuum normalization;
9. photon mass protection;
10. cutoff dispersion, birefringence, and Lorentz-violating operators quantified.

The finite-link work satisfies gate 1 and identifies the first tested nonconstant electric representation. The issue #4 construction supplies a finite endpoint and exact gauge-covariant composite hopping for gate 6. It does not establish gates 2–5 or 7–10. Do not infer a Coulomb phase from one plaquette, an infinite-rotor target, or the existence of a covariant defect.

## 5. Gravitational gates

A gravity candidate must demonstrate:

1. three momentum constraints and one scalar constraint;
2. correct first- and second-class classification;
3. exactly two positive physical tensor modes;
4. no unconstrained scalar/vector pole;
5. linear gap closing and regulator-anisotropy control;
6. long-range attraction with the correct source sign;
7. one shared frame coupling universally to total stress-energy;
8. gravitational self-energy in the nonlinear source;
9. nonlinear constraint preservation without manual projection;
10. stable strong-field and background solutions;
11. no low-cutoff strong coupling or regenerated ghost;
12. a microscopic treatment of the volume/vacuum term.

The linear mode-count and positive finite-momentum spectrum gates are internally passed. The finite matter prototype uses the same bare frame derivative for every charge species, which is a useful gate-7 input, but universal interacting stress-energy coupling is not yet demonstrated. The nonlinear, self-coupling, finite dressed-Hamiltonian, and vacuum gates remain open.

## 6. Matter gates

A matter construction must provide:

1. an explicit finite local Hilbert-space defect or endpoint;
2. exact gauge-covariant hopping;
3. universal coupling through the same frame;
4. fermionic exchange statistics or a declared alternative matter target;
5. a controlled treatment of chirality and species doubling;
6. anomaly cancellation;
7. particles and antiparticles;
8. charge assignments derived or constrained by the microscopic algebra;
9. protected small gaps rather than arbitrary tiny on-site energies;
10. material clocks and bound-state energies calculable from the same defects.

The finite issue #4 prototype passes gates 1–9 within its declared free-regulator scope:

- a charge-`q` endpoint is an exact bound composite of `|q|` finite fermionic strands;
- the same signed strand number terminates `|q|` spin-1 flux lanes and fixes the exchange sign;
- the selected primitive chiral spectrum `(-11,-5,-1,-1,9,9)` satisfies `sum(q)=sum(q^3)=0` and has no vectorlike pair;
- composite hopping commutes exactly with both endpoint Gauss generators;
- all species share one frame derivative;
- one Wilson/domain-wall cone appears on each wall while the other seven physical corners remain gapped;
- the remote mirror wall is retained explicitly rather than hidden;
- positive/negative energy pairing supplies particle and antiparticle branches;
- same-wall bilinear masses are forbidden, and the remaining finite-width gap is exponentially protected by wall separation;
- one common charge-localization rule gives four distinct gaps from two shared parameters.

This does not pass gate 10. It also does not derive Standard Model quantum numbers, observed masses, mirror-wall decoupling, or an interacting chiral gauge theory. The finite prototype closes issue #4 at its first-defect/regulator scope; issues #3, #6, and #8 own coefficient derivation, dynamical QED, mirror completion, Ward identities, and radiative stability.

Inserting a continuum Dirac field remains acceptable only as an interim infrared comparison, not as microscopic closure.

## 7. Shared-parameter gate

One finite microscopic model must derive or constrain

\[
U_A,\ K_A,\ U_g,\ K_g,\ \ell,
\]

plus matter hopping, the matter localization parameters, defect gaps, interaction strengths, and any volume parameters.

Define

\[
Z_A=\sqrt{U_A/K_A},
\qquad
Z_g=\sqrt{U_g/K_g}.
\]

A common limiting speed constrains products, not the impedance ratio. Equality of `Z_A` and `Z_g` is optional until derived.

The issue #4 prototype uses

\[
r_q=r_0+\eta(|q|-1),
\qquad
\Delta_q\propto r_q^{L_s},
\]

which is a protected common relation with fewer parameters than distinct gaps. It is not yet parameter closure because `r_0` and `eta` are not derived from the same move amplitudes that set the gauge and frame stiffnesses.

Parameter closure requires:

- fewer independent microscopic ratios than low-energy outputs;
- all perturbative denominators and combinatorial factors recorded;
- omitted-order estimates;
- no use of measured `G`, `alpha`, or particle masses to select microscopic amplitudes before the relation is derived;
- at least one held-out cross-sector relation.

## 8. Continuum quantum gate

A finite network does not automatically define a healthy quantum field theory. Require:

1. a declared scaling limit or fixed point;
2. physical propagators with positive residues;
3. unitarity or reflection positivity;
4. gauge/constraint Ward identities;
5. regulator universality across at least two implementations;
6. Lorentz and rotational recovery for photons, gravitons, and matter;
7. radiative stability of zero or small masses and universal couplings;
8. complete anomaly accounting;
9. leading irrelevant-operator inventory;
10. explicit audit of the assumptions behind emergent-spin-2 no-go results.

The reduced free tensor transfer matrix and finite free matter slab are useful benchmarks. Neither establishes the interacting continuum theory. In particular, the remote matter mirror must be removed, symmetrically gapped, or shown to decouple without violating the intended gauge symmetry.

## 9. Photon–companion bridge gates

The bridge to the active repository program must derive from the same microscopic theory:

- emitter conversion or excitation production;
- photon frequency and linewidth response;
- companion speed, dispersion, polarization, and coherence;
- capture/binding/storage and escape;
- reverse conversion and fluctuations;
- source fuel, recoil, momentum, angular momentum, and boundary flow;
- gravitational response of stored excitation;
- common matter and light trajectories;
- clock response;
- redshift, duration, brightness, and angular-size predictions;
- galaxy and cluster motion/lensing from the same sourced frame;
- background/microwave observables under the universe contract.

The new charged defects provide the finite matter participant required for a matter-assisted photon–`varphi` transition. A variable-name mapping or an interaction written without recoil and reverse channel does not satisfy this gate.

## 10. Rejection criteria before data fitting

Reject or explicitly revise a candidate if any of the following occurs:

- a physical negative-energy or negative-residue mode;
- `omega^2<0` near the vacuum;
- extra unsuppressed photon or tensor polarizations;
- a finite electromagnetic model that remains confined or gapped;
- tensor dispersion `omega ~ k^2` or `k^3` in the intended relativistic phase;
- uncontrolled finite-state boundary occupation;
- loss of charge or stress-energy conservation;
- species-dependent gravity without a complete additional-force theory;
- separately tuned photon, graviton, and matter cones;
- a volume term deleted by hand;
- a companion lifetime, capture law, halo, or particle hierarchy inserted through independently adjustable species energies;
- a matter regulator that hides an ungapped mirror or doubler;
- one field fitted to motion and another independently fitted to lensing;
- no stable continuum phase;
- more adjustable parameters than independent outputs.

A protected hierarchy from one declared local rule is a valid microscopic candidate. It cannot be used empirically until the rule's coefficients are derived or explicitly counted as calibration inputs.

## 11. Synthetic test ladder

Run these before observational unblinding:

### A. Vacuum and phase diagram

Map gaps, residues, order parameters, flux sectors, wall/mirror sectors, and phase boundaries across microscopic couplings and local dimensions.

### B. Static sources

Use finite conserved electric charge and energy sources. Measure sign, range, anisotropy, finite-size corrections, saturation, and species dependence.

### C. Wave packets and correlators

Measure dispersion, speed, polarization, leakage, birefringence, tensor/vector/scalar content, chiral propagation, and reflection positivity.

### D. Matter response

Propagate multiple finite defect species through identical electromagnetic and frame backgrounds. Test charge response, universal free fall, mirror leakage, clocks, and binding energy.

### E. Radiation and reactions

Drive electric dipole and gravitational quadrupole sources. Verify reciprocal backreaction, emitted power, forbidden multipoles, and energy/angular-momentum ledgers.

### F. Nonlinear evolution

Evolve physical and constraint-violating data separately. Physical data must preserve constraints without repeated manual projection.

### G. End-to-end companion experiment

Start from finite matter and photons, generate the declared companion excitation, propagate it, capture/store it, source the frame, and detect both matter and light response with one shared parameter set.

## 12. Freeze before empirical testing

Commit a protocol containing:

- exact microscopic Hamiltonian or action;
- vacuum/phase definition;
- complete parameter ledger;
- derived versus calibrated quantities;
- calibration observables;
- held-out predictions;
- source/data provenance;
- numerical and statistical tolerances;
- prohibited post-unblinding modifications;
- code commit and input/output hashes.

The first empirical package should test local and weak-field consequences before flexible galaxy or cosmological models:

1. photon, tensor, and matter propagation cones;
2. free-fall universality;
3. Solar-System motion, light bending, delay, and redshift;
4. gravitational-wave polarization and dispersion;
5. precision electromagnetic cutoff operators and coupling variation;
6. binary radiation and strong-field consistency.

Only a branch passing those gates should proceed to galaxies, clusters, redshift phenomenology, or cosmological/background fits.

## 13. Distinctive-prediction requirement

The theory is not empirically distinct if it reproduces QED and general relativity by independently selecting all of their coefficients. It must provide at least one relation of the form

\[
\mathcal R(\alpha,G,m_s,c,\text{cutoff coefficients},\text{companion response})=0
\]

that is derived before testing and was not used for calibration.

Preferred predictions link more than one sector, for example:

- photon and tensor cutoff dispersion;
- coupling variation and a neutral junction resonance;
- a defect mass ratio and `Z_g/Z_A`;
- a redshift-transfer coefficient and a gravitational/storage response;
- a fixed relationship between motion and lensing corrections.

The issue #4 charge-to-gap relation is an internal candidate relation. It is not yet a prediction because its shared parameters have not been derived or assigned calibration roles.

## 14. Current decision

The linear gravity mode problem and the first finite matter endpoint are sufficiently developed to stop treating either as an undefined placeholder. Continue only calculations that close an architecture dependency:

1. build the finite matter-assisted photon–`varphi` transition with reverse channel and recoil in issue #7;
2. demonstrate or reject the finite Coulomb/QED phase with the new dynamical defects in issue #6;
3. finish the finite dressed gravity Hamiltonian in issue #2;
4. derive shared gauge, gravity, and matter-localization coefficients in issue #3;
5. close nonlinear constraints and self-coupling in issue #5;
6. establish mirror completion, interacting consistency, Lorentz recovery, and radiative stability in issue #8;
7. freeze and test predictions under issue #9.

No one pillar may be used as evidence that another has been solved.
