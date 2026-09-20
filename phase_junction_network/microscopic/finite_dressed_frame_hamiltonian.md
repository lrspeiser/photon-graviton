# Finite Dressed-Frame Hamiltonian and Transfer Matrix

**Status:** linear finite-regulator gate passed  
**Date:** 2026-09-20  
**Issue:** [#2](https://github.com/lrspeiser/photon-graviton/issues/2)

## 1. Result

The linear frame/connection theory now has a declared finite positive Hamiltonian—or, equivalently, a positive constrained transfer matrix—built from the exact dressed frame Weyl algebra.

The construction has four parts:

1. the connection-relative variables are locked by the already established finite Weyl oscillator;
2. all frame moves use
   \[
   \overline Z_i=Z_{h_i},
   \qquad
   \overline X_i=X_{h_i}\prod_aX_{C_a}^{A_{ai}},
   \]
   so they commute exactly with that lock;
3. the exact finite scalar and vector stabilizers are reduced over \(\mathbb F_p\), producing two—not assumed two—logical Weyl pairs for every nonzero propagation block;
4. those two physical fields carry the minimal positive self-dual clock/Villain Hamiltonian.

No transverse-traceless projector is inserted into the Hamiltonian or into the finite logical reduction.

For one physical polarization on an odd-prime lattice,

\[
\boxed{
H_{p,L}
=\frac{\lambda_g}{2}\sum_z
\left[
2-X_z-X_z^\dagger
+2-Z_{z+1}Z_z^\dagger-Z_zZ_{z+1}^\dagger
\right].
}
\]

The complete propagating Hamiltonian is two identical copies,

\[
H_{\rm tensor}=H_{p,L}^{(1)}+H_{p,L}^{(2)}.
\]

Every term is positive semidefinite. The corresponding transfer matrix

\[
T=e^{-\Delta\tau H_{\rm tensor}}
\]

is positive. In the unreduced redundant frame description, the same object is represented as a locally constrained transfer step,

\[
T_{\rm phys}=P_{\rm fc}\,T_0\,P_{\rm fc},
\]

where \(P_{\rm fc}\) is the product of finite group averages over the local scalar and vector first-class stabilizers. The physical logarithm need not be ultralocal after quotienting; locality belongs to the parent transfer step and its local constraints.

This distinction is important. The previously proved derivative-order obstruction still rejects an **unconstrained commuting-projector Hamiltonian made from squares of exact frame-only invariants**: that branch gives \(\omega\sim k^3\). The accepted branch is instead the constrained transfer-matrix/Dirac-reduction route already selected by the first-order analysis.

## 2. Exact finite logical reduction

For a principal-axis nonzero momentum block, use the six frame components

\[
(xx,yy,zz,xy,xz,yz).
\]

The finite scalar \(Z\)-stabilizer and three vector \(X\)-stabilizers have ranks one and three. The script does not prescribe plus and cross polarizations. It computes

\[
\ker G/\operatorname{row}C
\quad\text{and}\quad
\ker C/\operatorname{row}G
\]

over \(\mathbb F_p\), then symplectically dualizes the two quotient bases.

For every tested prime \(p=5,7,11\), the result is exactly

\[
6-1-3=2
\]

logical qudits, with canonical pairing matrix

\[
\begin{pmatrix}1&0\\0&1\end{pmatrix}
\pmod p.
\]

One selected basis is equivalent on the constraint surface to

- a cross-like field \(h_{xy}\);
- a plus-like field represented by \(h_{yy}\), whose dual shift is \(p_{yy}-p_{xx}\).

That interpretation is reported only after the modular quotient has selected the basis.

The full real-space axial constraint matrices also pass:

| \(L\) | \(p\) | frame qudits | scalar rank | vector rank | logical qudits | expected |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 5 | 18 | 2 | 6 | 10 | \(2L+4=10\) |
| 5 | 7 | 30 | 4 | 12 | 14 | \(2L+4=14\) |

All scalar/vector Weyl commutators vanish exactly modulo \(p\). The four excess periodic qudits are the already isolated global sector. The propagating test removes \(k=0\) and retains exactly \(2L\) logical qudits along the tested direction.

## 3. Finite spectrum

The calculation diagonalizes one physical clock chain in the global-shift-invariant sector of dimension

\[
p^{L-1}.
\]

Translation is then diagonalized inside each degenerate energy eigenspace. The first signed-momentum pair has lattice indices \(+1\) and \(-1\), and the second identical polarization makes the tensor first gap twofold degenerate.

The exact lattice momentum is

\[
\widehat k=2\sin\frac{\pi}{L}.
\]

With the declared test normalization \(\lambda_g=0.25\), the fitted form is

\[
\Delta_p(\widehat k)
=c_g(p)\widehat k+d_3(p)\widehat k^3.
\]

| \(p\) | tested \(L\) | power in \(\Delta\propto\widehat k^s\) | \(c_g(p)\) | \(d_3(p)\) | max relative fit residual |
|---:|---|---:|---:|---:|---:|
| 5 | 4–8 | 1.070492 | 0.219555 | 0.006830 | \(1.66\times10^{-3}\) |
| 7 | 4–7 | 1.037242 | 0.165438 | 0.002478 | \(1.78\times10^{-4}\) |
| 11 | 3–6 | 1.026897 | 0.117794 | 0.000916 | \(6.48\times10^{-4}\) |

Thus every tested finite local dimension has a linearly closing signed-momentum branch with the expected cubic lattice correction. The overall speed depends on the declared finite regulator normalization and is not yet a prediction.

Representative finite gaps are:

| \(p\) | \(L\) | \(\widehat k\) | tensor gap | lock gap / tensor gap |
|---:|---:|---:|---:|---:|
| 5 | 8 | 0.765367 | 0.170819 | 4.81 |
| 7 | 7 | 0.867767 | 0.145206 | 5.57 |
| 11 | 6 | 1.000000 | 0.118787 | 7.27 |

Across all tested rows, the connection-relative excitation remains at least 2.49 times above the tensor gap.

## 4. Finite-dimension boundary control

The compact clock regulator has a branch cut at the maximally separated bond value. The script measures both the mean fraction of bonds on that cut and the probability that any bond is on it.

For the first tensor state at the largest tested \(L\) for each prime:

| \(p\) | \(L\) | mean cut-bond fraction | probability of any cut bond |
|---:|---:|---:|---:|
| 5 | 8 | 0.045526 | 0.295958 |
| 7 | 7 | 0.007380 | 0.049316 |
| 11 | 6 | 0.000189 | 0.001129 |

The low-state compact-boundary error falls by more than two orders of magnitude between \(p=5\) and \(p=11\). The finite-\(p\) model is therefore not obtaining its low tensor branch from wraparound states.

## 5. Exact connection-lock compatibility

A dense coupled test at \(p=3\) uses two frame and two connection qudits with

\[
A=\begin{pmatrix}1&1\\1&2\end{pmatrix}.
\]

The full Hilbert dimension is 81. The connection lock has a nine-dimensional frame ground band, exactly \(3^2\). A positive frame kinetic-plus-bond operator made only from \(\overline Z\) and dressed \(\overline X\) gives

- maximum lock/frame commutator: \(1.67\times10^{-16}\);
- ground-band leakage operator norm: \(5.47\times10^{-16}\);
- lock gap: 0.826993.

The zero leakage is structural, not a fitted suppression: every dressed frame monomial commutes with the relative connection lock.

## 6. Complete three-dimensional linear symbol

The compact exact diagonalizations above test finite \(p\) and finite \(L\) nonperturbatively in a propagation sector. The complete three-dimensional parent kernel is checked independently by reducing the local Fierz–Pauli lattice symbol with all four constraints.

That reduction constructs the null space of the scalar/vector constraint matrix; it does not insert a formula for a TT projector. For \(L=12,16,24,32,48,64\):

- maximum relative frequency error from \(|\widehat{\mathbf k}|\): \(2.15\times10^{-16}\);
- maximum relative polarization split: \(2.83\times10^{-16}\);
- equal-continuum-norm cubic anisotropy scales as \(L^{-2.01663}\).

This verifies the complete linear spatial symbol while the finite clock spectra verify that a compact local-dimension regulator retains the massless branch.

## 7. Acceptance accounting

Issue #2’s linear finite gate is satisfied as follows:

| Requirement | Result |
|---|---|
| Declared finite Hamiltonian or positive transfer matrix | Two-copy positive clock/Villain Hamiltonian and constrained transfer step |
| Finite first-class treatment | Exact odd-prime stabilizer quotient |
| Exactly two physical tensor branches | Two modularly derived logical Weyl pairs; twofold first tensor gap |
| No scalar, vector, or connection branch in the physical spectrum | Removed by exact quotient; connection-relative band remains gapped |
| Linear finite-size gap | \(\Delta=c_g\widehat k+O(\widehat k^3)\) at \(p=5,7,11\) |
| No TT insertion | GF(\(p\)) constraint quotient and full constrained null-space reduction |
| Connection leakage | Below \(5.5\times10^{-16}\) in the dense test; algebraically zero |
| Finite-dimension control | Compact-cut occupation measured and rapidly decreasing with \(p\) |
| Global sector separated | \(k=0\) and fixed-volume choice remain explicit |
| Common frame ready for EM matching | Dynamics is written entirely in the exact dressed frame algebra |

## 8. What this does not solve

This closes the **linear finite-regulator** task, not the full theory.

1. The self-dual clock ratio and \(\lambda_g\) are regulator choices. Issue #3 must derive gravitational and electromagnetic stiffnesses from shared microscopic moves.
2. The compact completion has not been shown to lie in the same interacting universality class on every regulator. That belongs to issue #8.
3. Nonlinear constraints, gravitational self-energy, strong fields, and the volume/vacuum mechanism remain issue #5.
4. The common photon/gravity cone cannot be tested until issue #6 supplies a demonstrated finite electromagnetic Coulomb phase.
5. The calculation establishes the finite version of the already derived linear frame theory. It is not a unique ultraviolet completion and makes no new cosmological or galaxy-fit claim.

## 9. Reproduction

From the repository root:

```sh
python phase_junction_network/microscopic/check_finite_dressed_frame_hamiltonian.py \
  --output phase_junction_network/microscopic/finite_dressed_frame_results.json
```

The script exits nonzero if the logical count, modular commutators, signed-momentum identification, linear dispersion, lock separation, connection leakage, boundary suppression, polarization degeneracy, or three-dimensional anisotropy tests fail.
