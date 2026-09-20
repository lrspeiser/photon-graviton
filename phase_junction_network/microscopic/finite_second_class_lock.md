# Finite Weyl Lock for the Connection Second-Class Pair

**Status:** finite local second-class locking architecture passed  
**Date:** 2026-09-20

## 1. Problem

The continuum first-order reduction uses the second-class constraints

\[
P_C=0,
\qquad
C-Ah=0,
\]

where \(A\) is the discrete linearized connection map.

At finite Hilbert-space dimension, conjugate coordinate and momentum constraints cannot both be imposed as exact sharp eigenvalue equations. The finite model instead needs a gapped minimum-uncertainty lock whose low-energy band retains the frame degrees of freedom.

## 2. Finite Weyl variables

For an odd-prime qudit, let

\[
ZX=\omega XZ,
\qquad
\omega=e^{2\pi i/p}.
\]

Use frame qudits \(h_i\) and connection qudits \(C_a\). For a linear map \(A_{ai}\in\mathbb F_p\), define the relative-coordinate Weyl operator

\[
R_a
=Z_{C_a}\prod_i Z_{h_i}^{-A_{ai}}
\]

and the connection-momentum operator

\[
M_a=X_{C_a}.
\]

The local constraint-lock Hamiltonian is

\[
\boxed{
H_{\rm sc}
=\frac{p}{4\pi}
\sum_a
\left[
2-M_a-M_a^\dagger
+2-R_a-R_a^\dagger
\right].
}
\]

The first bracket penalizes nonzero connection momentum. The second penalizes mismatch between the connection coordinate and \(Ah\).

## 3. Exact factorization

The finite linear change of variables

\[
r_a=C_a-A_{ai}h_i
\]

is a bijection over \(\mathbb F_p\). In these variables,

\[
R_a=Z_{r_a},
\qquad
M_a=X_{r_a},
\]

so

\[
H_{\rm sc}
=I_{\rm frame}\otimes
\sum_a H_{\rm oscillator}(r_a).
\]

Therefore:

- each relative connection qudit becomes a gapped finite oscillator;
- the frame Hilbert space is left exactly intact;
- the ground-band dimension is
  \[
  p^{N_h},
  \]
  where \(N_h\) is the number of frame qudits;
- the connection excitations are separated by a nonzero gap.

This is the finite analogue of eliminating the second-class connection pair.

## 4. Dressed physical frame algebra

The ordinary frame coordinate operator is

\[
\overline Z_i=Z_{h_i}.
\]

A frame shift must move the connection with the frame so that \(C-Ah\) remains unchanged. Define

\[
\boxed{
\overline X_i
=X_{h_i}
\prod_a X_{C_a}^{A_{ai}}.
}
\]

These dressed operators satisfy

\[
[H_{\rm sc},\overline Z_i]=0,
\qquad
[H_{\rm sc},\overline X_i]=0,
\]

and retain the exact frame Weyl algebra,

\[
\overline Z_i\overline X_j
=\omega^{\delta_{ij}}
\overline X_j\overline Z_i.
\]

Thus all physical finite frame dynamics can be built from \(\overline Z_i\) and \(\overline X_i\) while remaining exactly inside the connection ground band.

This is the finite counterpart of extending the continuum vector generator so that it transforms both \(h\) and \(C\).

## 5. Constraint quality and gap scaling

For one relative qudit, define

\[
K_P=2-X-X^\dagger,
\qquad
K_Q=2-Z-Z^\dagger.
\]

The chosen scaling

\[
H_p=\frac{p}{4\pi}(K_P+K_Q)
\]

keeps the oscillator gap finite as \(p\to\infty\).

For primes

\[
p=3,5,7,11,17,23,31,47,67,97,
\]

the numerical fits over \(p\ge11\) give

\[
\langle K_P+K_Q\rangle
\propto p^{-0.971906},
\]

\[
\sqrt{\langle K_Q\rangle}
\propto p^{-0.485953},
\]

and

\[
1-\Delta_p
\propto p^{-0.983141},
\]

where \(\Delta_p\) is the scaled gap.

Therefore the coordinate and momentum constraint violations vanish in the large-\(p\) limit while the scaled excitation gap approaches one.

At \(p=97\):

- coordinate chord cost: \(0.0321260\);
- momentum chord cost: \(0.0321260\);
- total cost: \(0.0642520\);
- scaled gap: \(0.983892\).

## 6. Coupled finite examples

### Two frame and two connection qudits at \(p=3\)

Using

\[
A=
\begin{pmatrix}
1&1\\
1&2
\end{pmatrix},
\]

the full Hilbert dimension is 81.

Results:

- ground-band dimension: 9, exactly \(3^2\);
- gap above the band: \(0.826993\);
- maximum \([H,\overline Z]\) residual: 0;
- maximum \([H,\overline X]\) residual: \(1.33\times10^{-15}\);
- maximum dressed Weyl residual: \(6.50\times10^{-16}\).

### One frame and two connection qudits at \(p=5\)

Using

\[
A=
\begin{pmatrix}
1\\2
\end{pmatrix},
\]

the full Hilbert dimension is 125.

Results:

- ground-band dimension: 5, exactly \(5^1\);
- gap above the band: \(0.821373\);
- maximum \([H,\overline Z]\) residual: 0;
- maximum \([H,\overline X]\) residual: \(8.88\times10^{-16}\);
- maximum dressed Weyl residual: \(1.64\times10^{-16}\).

## 7. What this solves

The finite model now has a concrete local implementation of the continuum second-class pair:

\[
\boxed{
P_C=0,\quad C-Ah=0
\quad\longrightarrow\quad
\text{gapped relative Weyl oscillator}.
}
\]

It also provides an exact physical frame algebra that commutes with the connection lock. This means finite frame dynamics need not leak into connection excitations if it is written using the dressed operators.

The connection ground band is not merely approximately degenerate; its frame degeneracy is exact because the Hamiltonian factorizes in the relative variables.

## 8. What remains open

This does not yet complete the finite graviton Hamiltonian.

The next work is to:

1. express the finite frame kinetic and Fierz–Pauli interactions using \(\overline Z_i\) and \(\overline X_i\);
2. place those operators on the full lattice constraint skeleton;
3. verify two linearly dispersing tensor branches at finite \(p\) and finite \(L\);
4. measure corrections as \(p\) and \(L\) increase;
5. extend the relative lock to the nonlinear frame-dependent connection map;
6. calculate the resulting gravitational stiffness from the same microscopic amplitudes used in the electromagnetic link sector.

At finite \(p\), the connection constraints are minimum-uncertainty locks rather than simultaneous exact coordinate and momentum eigenconditions. Their violations vanish with increasing \(p\), while the dressed physical algebra remains exact at every tested \(p\).