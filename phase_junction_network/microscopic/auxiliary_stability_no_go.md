# Auxiliary-Connection Stability Result

**Status:** one microscopic implementation rejected; three viable verification branches separated  
**Date:** 2026-09-20

## 1. Result

The first-order identity

\[
V_{\rm FP}(h)=-4Q(\Gamma[h])
\]

and its local frame–connection representation remain algebraically correct. However, the connection cannot be promoted naively to an ordinary real positive-energy auxiliary oscillator while retaining a zero bare \(h\)-\(h\) block.

The obstruction is general.

Let

\[
H(h,C)=\frac12 C^T A C+C^T B h,
\]

where \(A\) is positive definite. The stationary connection is

\[
C_*=-A^{-1}Bh.
\]

Substitution gives

\[
H_{\rm eff}(h)=-\frac12h^TB^TA^{-1}Bh.
\]

Since

\[
B^TA^{-1}B\succeq0,
\]

the induced frame stiffness is negative semidefinite:

\[
\boxed{H_{\rm eff}\preceq0.}
\]

Therefore:

> A bounded real auxiliary field with no bare frame stiffness cannot generate the positive Fierz–Pauli tensor stiffness solely through a real linear coupling.

This is independent of the detailed Phase Junction model.

## 2. Application to the current connection form

For the implemented torsion-free connection quadratic form \(J\),

\[
Q(C)=C^TJC,
\]

the numerical inertia is

\[
\operatorname{inertia}(J)=(8\text{ negative},0\text{ zero},10\text{ positive}).
\]

The complete first-order frame–connection block has, at every tested nonzero momentum,

\[
\operatorname{inertia}(\mathcal V_1)=(9\text{ negative},3\text{ zero},12\text{ positive}).
\]

The three zero directions are the spatial gauge directions. Eliminating the connection gives the unconstrained Fierz–Pauli kernel with inertia

\[
(1\text{ negative},3\text{ zero},2\text{ positive}).
\]

The fourth scalar constraint removes the one negative scalar frame direction, leaving two positive transverse-traceless modes. It does not turn the 18 connection variables into ordinary positive-energy oscillators. They must remain constrained or auxiliary.

## 3. Three distinct continuation branches

### A. Constrained first-order phase-space branch

Treat the connection as a nonpropagating constrained variable, as in a Palatini or canonical first-order action. Quantize only after imposing or solving the second-class and first-class constraints.

Required checks:

- complete Dirac constraint classification;
- reduced physical symplectic form;
- positive physical Hamiltonian;
- finite implementation of the reduced transfer matrix;
- nonlinear closure;
- no hidden connection excitation.

This branch best preserves the claim that the frame stiffness is generated through first-order frame–connection dynamics.

### B. Euclidean imaginary-auxiliary branch

For \(A>0\),

\[
\int dC\,\exp\left[-\frac12C^TAC-iC^TBh\right]
\propto\exp\left[-\frac12h^TB^TA^{-1}Bh\right].
\]

The imaginary mixed coupling produces a positive Euclidean effective action. This is an exact Hubbard–Stratonovich identity.

The unresolved issues are substantial:

- the pre-integration measure is complex;
- a sign problem may prevent straightforward Monte Carlo;
- reflection positivity is not automatic;
- a local positive transfer matrix has not been derived;
- the finite gauge constraints still have to be imposed.

This route is useful as an independent path-integral verification even if it does not become the preferred microscopic ontology.

### C. Bounded Hamiltonian with bare frame stiffness

Use

\[
H=\frac12C^TAC+C^TBh+\frac12h^TPh,
\qquad A\succ0.
\]

Eliminating \(C\) gives

\[
H_{\rm eff}=\frac12h^T\left(P-B^TA^{-1}B\right)h.
\]

A desired positive physical stiffness \(M\) is obtained by choosing

\[
P=M+B^TA^{-1}B.
\]

The full block can then be positive. This produces a stable finite Hamiltonian, but the connection only localizes or renormalizes an already-present frame stiffness. It no longer explains the entire stiffness from a zero bare \(h\)-\(h\) block.

## 4. Verification results

The accompanying script performs 500 tests for each construction.

- Positive-auxiliary no-go failures: **0**
- Actual connection inertia: **8 negative, 10 positive**
- Full first-order block inertia: **9 negative, 3 zero, 12 positive**
- Unconstrained Fierz–Pauli inertia: **1 negative, 3 zero, 2 positive**
- TT relative eigenvalue error: \(1.06\times10^{-15}\)
- Stable-completion failures: **0**
- Stable-completion maximum Schur error: \(1.33\times10^{-15}\)
- Imaginary-auxiliary positive-kernel failures: **0**

## 5. Decision

Do not build the next finite model by assigning an ordinary positive kinetic and potential energy to all 18 connection components while leaving the bare frame block zero.

Advance the branches in this order:

1. **Constrained first-order reduction**, because it best matches the intended frame–connection mechanism.
2. **Euclidean auxiliary calculation**, as an independent propagator and correlation-function check.
3. **Positive bare-stiffness Hamiltonian**, as a stable benchmark and fallback microscopic model.

The first-order identity is still useful. What has been rejected is one particular quantization of its auxiliary variable, not the identity or the broader frame-holonomy mechanism.
