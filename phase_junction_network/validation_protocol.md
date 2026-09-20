# Phase Junction Network Validation Protocol

**Status:** finite kinematic stage partially completed; no observational fit has been performed for this branch.  
**Date:** 2026-09-20

## 1. Research rule

Internal consistency must be established before fitting observations. Real measurements are initially used as exclusion walls, not as flexible targets for choosing the action.

The sequence is:

\[
\text{finite junction rules}
\rightarrow
\text{effective constraints and spectrum}
\rightarrow
\text{parameter closure}
\rightarrow
\text{frozen predictions}
\rightarrow
\text{real-data tests}.
\]

## 2. Current stage status

| Gate | Status | Evidence |
|---|---|---|
| Electromagnetic finite local link with exact Gauss symmetry | **Pass** | Spin-S quantum-link plaquette; exact commutators in `microscopic/finite_junction_model.py`. |
| Nontrivial electric stiffness in the smallest tested link | **Pass at spin 1** | Spin-1/2 has constant \(E^2\); spin-1 is the first three-state nontrivial link. |
| Four finite gravitational constraints with correct local count | **Pass kinematically** | Exact odd-prime Weyl/CSS algebra leaves two local logical modes per site plus four periodic global modes. |
| Exact continuous canonical gravity algebra on a full finite local Hilbert space | **Impossible as stated** | Trace obstruction \(\operatorname{tr}[Q,P]=0\ne\operatorname{tr}(iI)\). Must emerge below a cutoff. |
| Manifest-local exact-invariant-square gravity dynamics | **Rejected for the relativistic target** | Polynomial search forces two coordinate derivatives and one momentum derivative, giving `omega ~ k^3`. |
| Local finite Hamiltonian with two linearly dispersing helicity-2 modes | **Open** | Requires a global Fierz–Pauli bilinear, frame-curvature term, auxiliary variables, or emergent symmetry. |
| Nonlinear constraint closure and universal self-coupling | **Open** | Not tested by the linear scripts. |
| Microscopic calculation of \(Z_g/Z_A\) | **Open** | Symmetry fixes operator structure, not all coefficients. |
| Chiral matter defects and protected mass hierarchy | **Open** | No finite defect construction yet. |
| Vacuum-volume term | **Open** | No balancing or unimodular mechanism yet. |

## 3. Zero-data gates

### Electromagnetic sector

The candidate must demonstrate:

1. an exact or controlled local Gauss constraint;
2. two and only two gapless transverse modes;
3. positive quadratic energy;
4. a long-distance Coulomb kernel;
5. gauge-covariant matter coupling;
6. absence of a photon mass unless the local rule is explicitly broken;
7. a stated continuum limit and cutoff corrections.

The finite-link pass establishes item 1 and the minimal nontrivial local representation. It does not by itself establish a three-dimensional Coulomb phase for spin 1.

### Gravitational sector

The candidate must demonstrate:

1. three momentum constraints and one curvature/energy constraint;
2. four first-class or exact discrete commuting constraints at the microscopic level;
3. exactly two positive transverse-traceless modes in the continuum phase;
4. no propagating negative-energy scalar;
5. a massless long-range static kernel approaching \(-1/r\);
6. one shared frame coupling to total conserved stress-energy;
7. a common low-energy characteristic cone for light and gravity;
8. controlled nonlinear constraint closure.

The current scripts establish the linear continuum target and an exact finite discrete constraint skeleton. They do not establish the nonlinear finite dynamics.

## 4. Finite-dimensional rule discovered in this phase

No finite local matrices obey an exact continuous canonical relation \([Q,P]=iI\) on their full Hilbert space. A microscopic gravity model must therefore choose one of two declared architectures:

1. **Exact discrete route:** odd-prime Weyl/qudit constraints are exact; continuous frame symmetry is emergent at long wavelength.
2. **Protected truncation route:** truncated oscillators or collective spins approximate the canonical algebra below a boundary penalty, with all cutoff leakage explicitly bounded.

A submission that silently treats finite matrices as an exact continuum canonical pair fails this protocol.

## 5. Immediate microscopic dynamics task

Use the finite constraint skeleton and enumerate local gauge-compatible operators. The first polynomial enumeration is complete: exact local coordinate invariants begin at derivative order 2 and exact local momentum invariants at order 1. A Hamiltonian made from their positive squares is therefore excluded for the relativistic target because it gives `omega ~ k^3`.

For the gravity qudits, let A denote the scalar Z-constraint matrix and B the vector X-constraint matrix. Candidate operators must satisfy:

\[
Bz=0\pmod p
\]

for Z-type coordinate terms and

\[
Ax=0\pmod p
\]

for X-type momentum terms.

For every candidate local move set:

1. enumerate support and symmetries;
2. verify exact commutation with the finite constraints;
3. determine whether the dynamical terms commute mutually or can sustain a gapless phase;
4. calculate the large-p or large-S quadratic kernel;
5. decompose all modes into scalar, vector, and tensor sectors;
6. reject all candidates without exactly two positive modes satisfying \(\omega\propto k\);
7. derive the effective coefficients by perturbation theory or direct spectrum matching.

For every perturbative candidate, calculate

\[
H_{\rm eff}=PVP+PV\frac{Q}{E_0-H_0}VP+\cdots
\]

through the first nonzero order. Record all intermediate-state energies, matrix elements, signs, and combinatorial factors. Do not promote the toy coefficient \(40t^4/\Delta^3\) to a universal result.

## 6. Required output of the microscopic stage

A candidate advances only if it produces, from one finite model:

\[
U_A,\ K_A,\ U_g,\ K_g,\ \ell,
\]

plus the defect hopping and gap parameters, with fewer independent microscopic ratios than low-energy observables.

The critical dimensionless outputs are

\[
Z_A=\sqrt{U_A/K_A},
\qquad
Z_g=\sqrt{U_g/K_g},
\qquad
\frac{Z_g}{Z_A},
\]

and the defect gaps

\[
\mu_s=\frac{m_sc\ell}{\hbar}.
\]

A model that simply assigns independent values to these quantities has not unified them. A common speed condition constrains the products \(U_AK_A\) and \(U_gK_g\), but does not by itself determine the impedance ratio.

## 7. Rejection criteria before using data

Reject or explicitly revise a microscopic candidate if any of the following occurs:

- a physical negative-energy mode;
- \(\omega^2<0\) near the proposed vacuum;
- more than two unsuppressed photon or tensor polarizations;
- a longitudinal photon in the exact massless phase;
- an unconstrained gravitational scalar with order-one coupling;
- a gravity construction that has only a gapped commuting-projector phase;
- tensor dispersion \(\omega\propto k^2\) or \(k^3\) when the intended target is relativistic;
- uncontrolled occupation of the finite canonical-pair boundary state;
- loss of local charge or stress-energy conservation;
- species-dependent frame coupling without a consistent additional force sector;
- unrelated photon and gravity propagation cones requiring numerical tuning;
- absence of a stable or well-defined continuum phase;
- a vacuum-volume term that is merely deleted rather than controlled by a mechanism;
- a microscopic local Lorentz-covariant conserved stress tensor that makes the proposed emergent spin-2 construction conflict with the assumptions of the Weinberg–Witten result.

## 8. Synthetic tests after finite-model closure

### Vacuum spectrum

Diagonalize the quadratic kernel over representative momenta and verify the complete scalar/vector/tensor spectrum, residues, and signs.

### Finite-size scaling

Run multiple odd lattice sizes and local dimensions. A putative gapless mode must show a gap closing with the declared power of L, while unwanted sectors remain gapped or constrained.

### Static sources

Place conserved charge and energy sources on large finite networks. Measure the Green functions, finite-size corrections, anisotropy, and crossover from microscopic to continuum behavior.

### Wave packets

Propagate photon and tensor packets to measure speed, dispersion, birefringence, mode leakage, and energy conservation.

### Cutoff leakage

For truncated oscillator or collective-spin realizations, measure occupation of the boundary states and confirm the predicted \(n/S\) or top-state error bound throughout evolution.

### Matter response

Propagate at least two distinct defect species through the same static frame field. Check universality of free fall and coupling to electromagnetic binding energy.

### Radiation

Drive electric dipole and gravitational quadrupole sources. Calculate emitted power and verify that forbidden lower multipoles are absent for the appropriate conservation laws.

### Nonlinear closure

Evolve constraint-violating perturbations and physical perturbations separately. Physical initial data must preserve the full nonlinear constraints without manual projection at every step.

## 9. Freeze before observational testing

Before loading detailed observations, commit a protocol containing:

- exact microscopic Hamiltonian or action;
- vacuum and phase definition;
- complete independent parameter list;
- calibration observables, if any;
- held-out predictions;
- numerical tolerances and failure thresholds;
- code commit and input hashes;
- prohibited post-unblinding changes.

A fit is informative only when the number of independent observables exceeds the number of adjustable parameters and when the held-out quantities were declared in advance.

## 10. First real-data test package

After the equations and parameters are frozen, use independent local/weak-field tests before galaxy or cosmological fits:

1. equality of electromagnetic and gravitational propagation cones;
2. composition dependence of free fall;
3. Solar-System light propagation and orbital weak-field response;
4. gravitational-wave polarization and dispersion;
5. laboratory or atomic-clock variation of dimensionless couplings;
6. precision QED bounds on higher-dimension electromagnetic operators.

Only a branch that passes those tests should be extended to binary pulsars, strong-field waveforms, galaxy dynamics, clusters, lensing, or cosmology.

## 11. Distinctive prediction requirement

The branch is not empirically distinct if it reproduces QED and general relativity by freely choosing all low-energy coefficients. It must eventually provide at least one overconstrained relation such as

\[
\mathcal R\left(\alpha,G,m_s,c,\text{correction coefficients}\right)=0
\]

that follows from the finite junction model and was not used for calibration.

The preferred target is a correlated set of effects controlled by the same microscopic ratio, for example a fixed relationship among:

- the photon/gravity cutoff dispersion;
- a variation of \(\alpha\) through a neutral junction mode;
- a tensor-sector correction;
- one or more protected matter-defect gaps.

## 12. Current decision

The finite kinematic stage has passed for the electromagnetic link and for a discrete version of the four gravity constraints. The next decisive work is a local finite Hamiltonian whose **dynamics**, not just its constraint count, produces a linearly dispersing two-helicity tensor phase and computes a shared electromagnetic/gravitational coupling ratio.
