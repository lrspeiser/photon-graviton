# Phase Junction Network Validation Protocol

**Status:** linear gravitational structure and several finite kinematic constructions are internally checked. The full electromagnetic, matter, nonlinear, continuum, integration, and empirical gates remain open.  
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
| Finite deconfined Coulomb/QED phase | **Open** | Issue #6 |
| Four finite gravitational constraints and local two-mode count | **Pass kinematically** | Odd-prime Weyl/CSS skeleton |
| Local frame/connection reduction to two positive linear modes | **Pass at the linear level** | Constrained local and full real-space reductions; no TT projector |
| Finite connection second-class lock | **Pass for the linear map** | Relative Weyl oscillator and exact dressed frame algebra |
| Finite dressed gravity Hamiltonian at finite local dimension | **Open** | Issue #2 |
| Nonlinear gravitational constraint closure and self-coupling | **Open** | Issue #5 |
| Linear periodic homogeneous conformal direction | **Controlled by a candidate** | Fixed-volume and trace-momentum second-class pair |
| Microscopic vacuum/volume mechanism | **Open** | Issue #5 |
| Finite charged chiral matter and protected masses | **Open** | Issue #4 |
| Shared microscopic `U_A,K_A,U_g,K_g` and `Z_g/Z_A` | **Open** | Issue #3 |
| Companion identity and bridge to active redshift/deposition program | **Open** | Issue #7 |
| Interacting continuum, unitarity, Lorentz recovery, anomalies | **Open** | Issue #8 |
| Frozen cross-sector predictions and empirical ladder | **Open** | Issue #9 |
| Reproducible regression suite | **Pass for current checks** | GitHub Actions and frozen outputs |

## 3. Ontology and integration gate

Before claiming a unified theory, declare one state/field dictionary containing:

- compact electromagnetic link phase and electric imbalance;
- frame/tetrad variables and connection variables;
- neutral relative-phase mode, if retained;
- charged and neutral matter defects;
- the active companion excitation used in redshift/capture/deposition work;
- global or boundary variables such as volume constraints.

Every physical energy contribution must appear exactly once. The same excitation may transform between sectors only through an explicit interaction with reciprocal source terms.

Issue #7 passes only when one common action or Hamiltonian covers emission, propagation, capture/storage, reverse transitions, stress-energy sourcing, clocks, motion, and lensing.

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

The existing finite-link calculation satisfies only the first gate and identifies the first tested nonconstant electric representation. Do not infer a Coulomb phase from one plaquette or from the infinite-rotor continuum target.

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

The linear mode-count and positive finite-momentum spectrum gates are internally passed. The nonlinear, self-coupling, finite dressed-Hamiltonian, and vacuum gates remain open.

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

Inserting a continuum Dirac field for convenience is acceptable as an interim infrared comparison, not as closure of the microscopic matter problem.

## 7. Shared-parameter gate

One finite microscopic model must derive or constrain

\[
U_A,\ K_A,\ U_g,\ K_g,\ \ell,
\]

plus matter hopping, defect gaps, interaction strengths, and any volume parameters.

Define

\[
Z_A=\sqrt{U_A/K_A},
\qquad
Z_g=\sqrt{U_g/K_g}.
\]

A common limiting speed constrains products, not the impedance ratio. Equality of `Z_A` and `Z_g` is an optional hypothesis until derived.

Parameter closure requires:

- fewer independent microscopic ratios than low-energy outputs;
- all perturbative denominators and combinatorial factors recorded;
- omitted-order estimates;
- no use of measured `G` or `alpha` to select the microscopic amplitudes before the relation is derived;
- at least one held-out cross-sector relation.

## 8. Continuum quantum gate

A finite network does not automatically define a healthy quantum field theory. Require:

1. a declared scaling limit or fixed point;
2. physical propagators with positive residues;
3. unitarity or reflection positivity;
4. gauge/constraint Ward identities;
5. regulator universality across at least two implementations;
6. Lorentz and rotational recovery for photons, gravitons, and matter;
7. radiative stability of zero masses and universal couplings;
8. anomaly accounting;
9. leading irrelevant-operator inventory;
10. explicit audit of the assumptions behind emergent-spin-2 no-go results.

The reduced free tensor transfer matrix is a useful benchmark. It does not establish the interacting continuum theory.

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

A variable-name mapping does not satisfy this gate.

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
- a companion lifetime, capture law, halo, or mass hierarchy inserted without dynamics;
- one field fitted to motion and another independently fitted to lensing;
- no stable continuum phase;
- more adjustable parameters than independent outputs.

## 11. Synthetic test ladder

Run these before observational unblinding:

### A. Vacuum and phase diagram

Map gaps, residues, order parameters, flux sectors, and phase boundaries across microscopic couplings and local dimensions.

### B. Static sources

Use finite conserved electric charge and energy sources. Measure sign, range, anisotropy, finite-size corrections, and species dependence.

### C. Wave packets and correlators

Measure dispersion, speed, polarization, leakage, birefringence, tensor/vector/scalar content, and reflection positivity.

### D. Matter response

Propagate multiple defect species through identical electromagnetic and frame backgrounds. Test charge response, universal free fall, clocks, and binding energy.

### E. Radiation and reactions

Drive electric dipole and gravitational quadrupole sources. Verify reciprocal backreaction, emitted power, forbidden multipoles, and energy/angular-momentum ledgers.

### F. Nonlinear evolution

Evolve physical and constraint-violating data separately. Physical data must preserve constraints without repeated manual projection.

### G. End-to-end companion experiment

Start from ordinary matter and photons, generate the declared companion excitation, propagate it, capture/store it, source the frame, and detect both matter and light response with one shared parameter set.

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

1. photon and gravitational propagation cones;
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

## 14. Current decision

The linear gravity mode problem is sufficiently developed to stop treating it as the whole project. Continue only gravity calculations that close issue #2 or #5.

The immediate program is:

1. define the Phase Junction/companion ontology and common ledger in issue #7;
2. demonstrate or reject the finite Coulomb/QED phase in issue #6;
3. finish the finite dressed gravity Hamiltonian in issue #2;
4. build the first finite charged matter defect in issue #4;
5. derive shared coefficients in issue #3;
6. test nonlinear and continuum consistency in issues #5 and #8;
7. freeze and test predictions under issue #9.

No one pillar may be used as evidence that another has been solved.