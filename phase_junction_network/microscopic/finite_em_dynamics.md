# Stage 6B: Detuned Finite Electromagnetic Dynamics

**Status:** **PASS at the pure-gauge dynamical stage.** A finite detuned region has two transverse branches with approximately linear long-wavelength gaps; the RK point remains the required approximately quadratic control. Dynamical charged matter, Ward identities, vacuum polarization, and interacting QED remain open.  
**Date:** 2026-09-20  
**Issue:** #6  
**Executable:** `check_finite_em_dynamics.py`  
**Frozen output:** `finite_em_dynamics_results.json`

## 1. Question

Stage 6A established a finite spin-1 static Coulomb ground-state anchor at the RK point

\[
u/t=0,\qquad v/t=1.
\]

That did not establish photons: RK dynamics may have \(z=2\), while a relativistic Coulomb phase requires two dynamical transverse branches with

\[
\Delta_\gamma(\mathbf k)
=c_\gamma|\widehat{\mathbf k}|+O(|\widehat{\mathbf k}|^3).
\]

Stage 6B asks whether detuning the same finite link Hamiltonian produces a non-isolated region with that behavior.

## 2. Parent finite Hamiltonian

The Stage-6A model is

\[
H_A(u,t,v)
=\frac{u}{2}\sum_\ell E_\ell^2
-t\sum_p(W_p+W_p^\dagger)
+v\sum_pD_p,
\]

where

\[
D_p=\sqrt{W_p^\dagger W_p}+\sqrt{W_pW_p^\dagger}.
\]

All links remain finite integer-spin quantum links. No continuum rotor replaces them.

## 3. Exact transverse propagation block

For propagation along one principal axis, Gauss reduction leaves two independent transverse electric-flux fields. One polarization is the exact finite chain

\[
\boxed{
H_\perp/t
=\frac{u/t}{2}\sum_z E_z^2
-\sum_z(W_z+W_z^\dagger)
+\frac vt\sum_zD_z,
}
\]

with

\[
W_z=U_{z+1}U_z^\dagger.
\]

The zero-total-flux sector is used. The second polarization is an identical copy. A longitudinal link variable is absent because it is removed by the Gauss constraint before diagonalization.

This block uses exactly the same finite raising/lowering amplitudes and RK flippability term as the three-dimensional parent Hamiltonian.

## 4. Declared detuned region

The frozen full run tests

\[
u/t\in\{0,0.1,0.2\},
\qquad
v/t\in\{0.2,0.4,0.6\}.
\]

Thus 18 spin/coupling combinations are tested:

- spin 1 on \(L=6,8,10\);
- spin 2 on \(L=4,5,6\).

The RK controls use the same size sets at \((u/t,v/t)=(0,1)\).

For every detuned block, the lowest signed-momentum gap is fit both to

\[
\Delta\propto\widehat k^s
\]

and to

\[
\Delta=c_\gamma\widehat k+d_3\widehat k^3.
\]

The acceptance range for the finite-size exponent is

\[
0.94\le s\le1.16,
\]

with a maximum relative linear-plus-cubic residual below \(2\times10^{-3}\).

## 5. Detuned spectrum result

All 18 detuned combinations pass.

Across spin 1 and spin 2:

\[
\boxed{0.978964\le s\le1.152567}
\]

and

\[
\boxed{
\max\frac{|\Delta_{\rm fit}-\Delta|}{\Delta}
=4.31\times10^{-4}.
}
\]

Representative fitted speeds in the declared lattice normalization are:

| Spin | \(u/t\) | \(v/t\) | power \(s\) | \(c_\gamma\) | max fit residual |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.0 | 0.2 | 1.01181 | 1.43244 | \(2.01\times10^{-4}\) |
| 1 | 0.1 | 0.4 | 1.06409 | 1.21406 | \(2.03\times10^{-4}\) |
| 1 | 0.2 | 0.6 | 1.14597 | 0.97611 | \(3.27\times10^{-4}\) |
| 2 | 0.0 | 0.2 | 0.97896 | 0.87510 | \(2.91\times10^{-4}\) |
| 2 | 0.1 | 0.4 | 1.02071 | 0.78226 | \(3.00\times10^{-4}\) |
| 2 | 0.2 | 0.6 | 1.08097 | 0.68297 | \(4.06\times10^{-4}\) |

The speed is representation- and coupling-dependent in bare units. Stage 6B establishes the linear branch and its finite region; it does not derive the physical normalization of \(c_\gamma\).

The \(+k\) and \(-k\) gaps are degenerate to the reported numerical precision. Two identical transverse copies therefore give a twofold photon gap.

## 6. RK negative control

The same fit at the RK point gives

| Spin | fitted power |
|---:|---:|
| 1 | 1.96037 |
| 2 | 1.82408 |

Thus the solvable static RK anchor is not relabeled as the relativistic photon phase:

\[
\boxed{
\text{RK control: }z\approx2,
\qquad
\text{detuned region: }z\approx1.
}
\]

This is the principal Stage-6B result.

## 7. Spectral residue and scalar control

At the representative point

\[
u/t=0.1,\qquad v/t=0.4,
\]

the transverse electric operator

\[
E(k)=\sum_z e^{-ikz}E_z
\]

places almost all its spectral weight in the first signed-momentum branch:

| Spin | \(L\) | photon gap | electric residue in first branch | scalar residue in first branch | scalar gap / photon gap |
|---:|---:|---:|---:|---:|---:|
| 1 | 10 | 0.764792 | 0.992529 | \(6.93\times10^{-33}\) | 2.84960 |
| 2 | 8 | 0.600370 | 0.998534 | \(5.46\times10^{-31}\) | 2.91992 |

The scalar control operator is the Fourier transform of \(E_z^2\). It has negligible residue in the photon branch, and its first material spectral weight occurs at more than 2.8 times the photon gap.

## 8. Exact periodic-cube anchor

The full spin-1 zero-charge, zero-winding component of the periodic \(2^3\) cube has

\[
146{,}327
\]

gauge-reduced states. The complete three-dimensional sparse Hamiltonian is diagonalized at four points:

\[
(0,1),\quad(0,0.6),\quad(0.1,0.4),\quad(0.2,0.2).
\]

The detuned ground state remains continuously connected to the RK equal-amplitude state:

| \(u/t\) | \(v/t\) | RK-state overlap |
|---:|---:|---:|
| 0.0 | 1.0 | 1.000000 |
| 0.0 | 0.6 | 0.882720 |
| 0.1 | 0.4 | 0.782560 |
| 0.2 | 0.2 | 0.672287 |

For momentum \((\pi,0,0)\):

- the longitudinal electric operator has exactly zero norm in the gauge-reduced basis;
- the low multiplet carries transverse electric residue;
- the scalar electric-energy operator has at most \(2.83\times10^{-21}\) of its weight in that multiplet.

This is a finite full-cube anchor. It does not replace an infinite-volume three-dimensional spectral calculation.

## 9. Three-dimensional tensor structure

The cubic long-wavelength completion uses the exact lattice momentum

\[
\widehat k_i=2\sin(k_i/2)
\]

and the Gauss-reduced kernel

\[
K_{ij}(\mathbf k)
=\widehat k^2\delta_{ij}-\widehat k_i\widehat k_j.
\]

For every tested nonzero momentum, its spectrum is

\[
\boxed{\{0,\widehat k^2,\widehat k^2\}},
\]

with one constrained longitudinal direction and two equal transverse directions.

For equal continuum-norm modes \((3,0,0)\) and \((2,2,1)\), cubic anisotropy scales as

\[
L^{-2.01111},
\]

consistent with an \(O(a^2k^2)\) regulator correction.

## 10. Acceptance accounting

| Stage-6B gate | Result |
|---|---|
| Finite detuned region rather than one point | Pass: 18 spin/coupling combinations |
| Gap closes approximately linearly | Pass: powers 0.979–1.153 |
| Linear plus cubic cutoff fit | Pass: max residual \(4.31\times10^{-4}\) |
| Two transverse copies | Pass |
| Signed momentum degeneracy | Pass |
| Positive photon spectral residue | Pass: minimum 0.992529 |
| Scalar rejection | Pass: residue below \(5.47\times10^{-31}\), gap ratio above 2.8496 |
| RK point rejected as \(z=1\) | Pass: powers 1.824–1.960 |
| Full-cube continuity anchor | Pass |
| Longitudinal rejection | Pass |
| Spin-2 representation control | Pass |
| Cubic anisotropy recovery | Pass: approximately \(L^{-2.01}\) |

## 11. What Stage 6B establishes

The current finite electromagnetic Hamiltonian has a detuned pure-gauge region whose exact principal-axis transverse sectors support two positive photon-like branches with

\[
\Delta_\gamma
=c_\gamma|\widehat k|+O(|\widehat k|^3).
\]

The RK point is a static Coulomb anchor with approximately quadratic dynamics; detuning toward smaller \(v/t\) produces the linearly dispersing branch.

This completes the pure-gauge dynamical step requested after Stage 6A.

## 12. Claim boundary

Stage 6B does **not** establish:

- the complete infinite-volume three-dimensional many-body spectrum away from RK;
- dynamical issue-#4 matter inside the detuned phase;
- finite interacting Ward identities;
- vacuum polarization or charge renormalization;
- a regulator-independent QED fixed point;
- precision Lorentz-violation bounds;
- a common photon, companion, matter, and tensor cone;
- the physical normalization of \(c_\gamma\), \(\alpha\), or its relation to \(G\).

Issue #6 therefore remains open for its interacting-QED continuation.

## 13. Next gate

Embed the minimal \(|q|=1\) issue-#4 endpoint in the detuned region and test:

1. exact finite current continuity;
2. a finite Ward identity;
3. photon spectral dressing and charge renormalization;
4. vacuum-polarization response;
5. the issue-#7 reversible photon–\(\chi\) vertex using actual photon eigenmodes;
6. whether the matter, photon, companion, and tensor velocities can share one continuum cone without independent retuning.

## 14. Reproduction

```sh
python phase_junction_network/microscopic/check_finite_em_dynamics.py \
  --output phase_junction_network/microscopic/finite_em_dynamics_results.json
```

A smaller regression run is available with `--quick`.
