# Phase Junction Network Validation Protocol

**Status:** prospective protocol; no observational fit has been performed for this branch.  
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

## 2. Current zero-data gates

### Electromagnetic sector

The candidate must demonstrate:

1. an exact or controlled local Gauss constraint;
2. two and only two gapless transverse modes;
3. positive quadratic energy;
4. a long-distance Coulomb kernel;
5. gauge-covariant matter coupling;
6. absence of a photon mass unless the local rule is explicitly broken;
7. a stated continuum limit and cutoff corrections.

### Gravitational sector

The candidate must demonstrate:

1. three momentum constraints and one curvature/energy constraint;
2. four first-class constraints at the linear level;
3. exactly two positive transverse-traceless modes;
4. no propagating negative-energy scalar;
5. a massless long-range static kernel approaching \(-1/r\);
6. one shared frame coupling to total conserved stress-energy;
7. a common low-energy characteristic cone for light and gravity;
8. controlled nonlinear constraint closure.

The included script checks the linear mode count, the tensor stiffness, the unwanted scalar branch when the fourth constraint is removed, and the synthetic lattice \(1/r\) kernel. It does not check nonlinear closure.

## 3. Immediate microscopic construction task

Define a finite Hilbert space on sites and/or links with:

- matter-like and geometry-like internal channels;
- compact scalar phase comparisons;
- frame/tetrad comparison variables;
- local penalties that enforce the electromagnetic and gravitational constraints;
- allowed local swaps that generate plaquette and frame-holonomy terms perturbatively.

For every candidate move set, calculate:

\[
H_{\rm eff}=PVP+PV\frac{Q}{E_0-H_0}VP+\cdots
\]

through the first nonzero order that generates both sectors. Record all intermediate-state energies, matrix elements, signs, and combinatorial factors. Do not promote the toy coefficient \(40t^4/\Delta^3\) to a universal result.

## 4. Required output of the microscopic stage

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

A model that simply assigns independent values to these quantities has not unified them.

## 5. Rejection criteria before using data

Reject or explicitly revise a microscopic candidate if any of the following occurs:

- a physical negative-energy mode;
- \(\omega^2<0\) near the proposed vacuum;
- more than two unsuppressed photon or tensor polarizations;
- a longitudinal photon in the exact massless phase;
- an unconstrained gravitational scalar with order-one coupling;
- loss of local charge or stress-energy conservation;
- species-dependent frame coupling without a consistent additional force sector;
- unrelated photon and gravity propagation cones requiring numerical tuning;
- absence of a stable or well-defined continuum phase;
- a vacuum-volume term that is merely deleted rather than controlled by a mechanism;
- a microscopic local Lorentz-covariant conserved stress tensor that makes the proposed emergent spin-2 construction conflict with the assumptions of the Weinberg–Witten result.

## 6. Synthetic tests after finite-model closure

Run these without observational fitting:

### Vacuum spectrum

Diagonalize the quadratic kernel over representative momenta and verify the complete scalar/vector/tensor spectrum, residues, and signs.

### Static sources

Place conserved charge and energy sources on large finite networks. Measure the Green functions, finite-size corrections, anisotropy, and crossover from microscopic to continuum behavior.

### Wave packets

Propagate photon and tensor packets to measure speed, dispersion, birefringence, mode leakage, and energy conservation.

### Matter response

Propagate at least two distinct defect species through the same static frame field. Check universality of free fall and coupling to electromagnetic binding energy.

### Radiation

Drive electric dipole and gravitational quadrupole sources. Calculate emitted power and verify that forbidden lower multipoles are absent for the appropriate conservation laws.

### Nonlinear closure

Evolve constraint-violating perturbations and physical perturbations separately. Physical initial data must preserve the full nonlinear constraints without manual projection at every step.

## 7. Freeze before observational testing

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

## 8. First real-data test package

After the equations and parameters are frozen, use independent local/weak-field tests before galaxy or cosmological fits:

1. equality of electromagnetic and gravitational propagation cones;
2. composition dependence of free fall;
3. Solar-System light propagation and orbital weak-field response;
4. gravitational-wave polarization and dispersion;
5. laboratory or atomic-clock variation of dimensionless couplings;
6. precision QED bounds on higher-dimension electromagnetic operators.

Only a branch that passes those tests should be extended to binary pulsars, strong-field waveforms, galaxy dynamics, clusters, lensing, or cosmology.

## 9. Distinctive prediction requirement

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

## 10. Current decision

The link-gauge and frame-gravity constructions are suitable as target infrared theories. The next decisive work is not another phenomenological data fit. It is the finite microscopic derivation that either generates their constraint structures and coupling ratio together or demonstrates that the Phase Junction interpretation cannot close as proposed.
