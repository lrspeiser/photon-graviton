# First-Order Frame–Connection Dynamics

**Status:** local classical linear architecture verified; finite quantum realization and nonlinear closure remain open  
**Date:** 2026-09-20

## 1. Why this branch is needed

The finite constraint skeleton leaves the desired two local gravitational degrees of freedom, but the first exact-local dynamics search found a structural problem. A coordinate observable invariant under

\[
\delta h_{ij}=\partial_i\xi_j+\partial_j\xi_i
\]

requires at least two derivatives, while a momentum observable invariant under the scalar momentum shift requires at least one. Squaring those individually invariant local observables produces

\[
H_h\sim k^4h^2,
\qquad
H_p\sim k^2p^2,
\qquad
\omega\sim k^3.
\]

The intended relativistic tensor phase therefore cannot come from a simple sum of positive squares of individually gauge-invariant stabilizers.

The natural escape is to use the original frame idea more literally: make the **frame connection an independent microscopic variable**. Curvature then contains one derivative of the connection, and the action can be linear in curvature rather than quadratic in a two-derivative invariant.

## 2. Variables

For one spatial Fourier mode, let

\[
h_{ij}=h_{ji}
\]

be the symmetric frame deformation. Introduce an independent torsion-free connection variable

\[
C^a{}_{bc}=C^a{}_{cb}.
\]

In three spatial dimensions this has

\[
3\times\frac{3(3+1)}2=18
\]

components.

The connection induced by a frame deformation is the linearized Levi-Civita expression

\[
\Gamma^a{}_{bc}[h]
=\frac12\left(
\partial_bh^a{}_c+
\partial_ch^a{}_b-
\partial^ah_{bc}
\right).
\]

In Fourier space, replace each derivative by its lattice momentum component.

## 3. Connection quadratic form

Define

\[
Q(C)=\delta^{ij}
\left(
C^a{}_{bi}C^b{}_{aj}
-C^a{}_{ij}C^b{}_{ab}
\right).
\]

Let \(B(C,D)\) be its symmetric polarization,

\[
B(C,D)=\frac12\left[Q(C+D)-Q(C)-Q(D)\right].
\]

On the 18-dimensional torsion-free connection space, the implemented quadratic matrix has full rank. Its eigenvalues have both signs; the connection is an auxiliary constrained variable, not an additional set of positive-energy propagating oscillators.

## 4. Exact Fierz–Pauli identity

The spatial Fierz–Pauli stiffness used in the preceding frame derivation is

\[
\begin{aligned}
V_{\rm FP}(h)={}&
(\partial_kh_{ij})^2
-2(\partial_ih_{ij})^2\\
&+2(\partial_ih_{ij})(\partial_jh)
-(\partial_kh)^2.
\end{aligned}
\]

Direct contraction gives the identity

\[
\boxed{
V_{\rm FP}(h)=-4Q\!\left(\Gamma[h]\right).
}
\]

The accompanying code verifies this identity over 500 random nonzero momenta and random symmetric tensors, with maximum relative error below \(9\times10^{-16}\).

## 5. Local first-order potential

Define

\[
\boxed{
V_1(h,C)=4\left[
Q(C)-2B\!\left(C,\Gamma[h]\right)
\right].
}
\]

This expression contains:

- an on-site algebraic \(C^2\) term;
- a mixed \(C\,\partial h\) term with only one spatial derivative;
- no explicit \((\partial h)^2\) term.

Varying with respect to the independent connection gives

\[
C=\Gamma[h]
\]

because the connection quadratic form is nondegenerate in the tested torsion-free representation. Substitution yields

\[
V_1\bigl(h,C=\Gamma[h]\bigr)=V_{\rm FP}(h).
\]

Equivalently,

\[
V_1(h,C)=4Q\!\left(C-\Gamma[h]\right)+V_{\rm FP}(h).
\]

Thus the connection enforces metric compatibility/torsion freedom dynamically, while eliminating it reproduces the required two-derivative frame stiffness.

## 6. Matrix and Schur-complement check

Let \(J\) be the 18-by-18 matrix of \(Q\), and let \(A(k)\) map the six frame components into the 18 connection components of \(\Gamma[h]\). The first-order quadratic block is

\[
\mathcal V_1=
\begin{pmatrix}h\\C\end{pmatrix}^{\!T}
\begin{pmatrix}
0 & -4A^TJ\\
-4JA & 4J
\end{pmatrix}
\begin{pmatrix}h\\C\end{pmatrix}.
\]

Eliminating \(C\) gives the Schur complement

\[
-4A^TJA=M_{\rm FP}(k).
\]

The numerical verification gives:

- connection quadratic rank: `18 / 18`;
- maximum relative Schur-complement error: `7.97e-16`;
- maximum relative gauge-null error: `1.31e-15`;
- maximum relative transverse-tensor eigenvalue error: `9.99e-16`;
- failures over 500 random momenta: `0`.

On the transverse-traceless subspace, the two eigenvalues are

\[
\lambda_+(k)=\lambda_\times(k)=k^2,
\]

so coupling this potential to the ordinary transverse kinetic term gives

\[
\omega^2=U_gK_gk^2.
\]

## 7. How it evades the cubic-dispersion obstruction

The previous obstruction applied to Hamiltonians built from **squares of individually local gauge invariants**. The first-order frame–connection density is different:

- the connection transforms together with the frame;
- individual mixed densities need not be invariant in isolation;
- the full periodic sum is invariant through the same integration-by-parts structure as the Fierz–Pauli action;
- the connection can be eliminated locally, leaving a globally gauge-invariant two-derivative bilinear.

The architecture is therefore local without requiring curvature squared.

## 8. Finite microscopic interpretation

A finite Phase Junction realization can assign:

- collective-spin or truncated-oscillator frame amplitudes to local frame components;
- finite connection variables to oriented frame-comparison channels;
- strong local compatibility terms implementing \(C\simeq\Gamma[h]\);
- the existing four constraints as exact discrete Weyl constraints or as protected emergent continuous constraints.

The new ingredient is that the frame and connection cannot be collapsed into one local curvature stabilizer. They must remain independent at the microscopic level long enough for an area-times-curvature or connection-times-frame-gradient interaction to form.

## 9. Remaining gates

This calculation does **not** yet establish a complete finite quantum phase.

1. The connection quadratic form is indefinite. In the constrained first-order action this is expected, but a finite Hamiltonian or transfer-matrix construction must demonstrate a stable physical spectrum and unitary evolution.
2. Exact continuous canonical and diffeomorphism algebras still cannot hold on the full finite Hilbert space. The symmetry must be discrete microscopically or emergent in a protected low-energy sector.
3. The nonlinear frame–connection constraints and their algebra have not been implemented.
4. The finite move amplitudes have not yet calculated \(U_g\), \(K_g\), or \(Z_g/Z_A\).
5. A shared finite connection for all matter species has not yet been used to construct chiral defects.
6. The allowed volume term remains uncontrolled.

## 10. Decision

The first-order branch passes the next algebraic gate:

\[
\boxed{
\text{local frame} + \text{independent connection}
\longrightarrow
\text{one-derivative microscopic coupling}
\longrightarrow
V_{\rm FP}
\longrightarrow
\omega\propto k.
}
\]

The project should advance this branch rather than the pure local-invariant-square/stabilizer branch. The next concrete implementation is a finite frame-and-connection transfer matrix or Hamiltonian on a small periodic lattice, followed by a full physical-spectrum and cutoff-leakage calculation.

## References used for comparison

- T. Regge, “General Relativity Without Coordinates,” *Il Nuovo Cimento* 19, 558–571 (1961).
- Z.-C. Gu and X.-G. Wen, “Emergence of helicity ±2 modes (gravitons) from qbit models,” arXiv:0907.1203.
- Standard Hilbert–Palatini / Einstein–Cartan first-order gravity, whose connection equation imposes torsion freedom in the spinless sector.
