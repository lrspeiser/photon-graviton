# Finite electromagnetic Coulomb stage

**Issue:** #6 — finite electromagnetic Coulomb and QED phase  
**Status:** **Stage 6A passed: finite static Coulomb ground state.** The relativistic `z=1` photon dynamics, dynamical charged matter, and interacting QED limit remain open.  
**Date:** 2026-09-20

## 1. Result in one sentence

The project now has an explicit **finite spin-1, 3+1-dimensional, exactly gauge-invariant Hamiltonian** with a positive Rokhsar–Kivelson (RK) point whose ground-state support passes winding, transverse-correlation, Wilson-shift, static-charge, finite-size, and larger-representation Coulomb diagnostics without replacing the links by continuum rotors.

This is the first actual finite electromagnetic phase result in the project. It is stronger than the previous four-link algebra check, but it is not yet the complete solution of issue #6 because an RK Coulomb ground state does not by itself establish two linearly dispersing quantum photon branches or QED radiative consistency.

## 2. Finite Hamiltonian

Put one integer-spin representation on every oriented spatial link:

\[
E_\ell=S^z_\ell,
\qquad
U_\ell=\frac{S^+_\ell}{\max |S^+_\ell|}.
\]

The local Hilbert-space dimension is finite,

\[
d_\ell=2S+1,
\]

and the exact link algebra is

\[
[E_\ell,U_\ell]=U_\ell.
\]

At a vertex \(x\),

\[
G_x=\sum_{\ell\,\mathrm{out\ of}\,x}E_\ell
-\sum_{\ell\,\mathrm{into\,}x}E_\ell-Q_x.
\]

For an oriented plaquette \(p\), define

\[
W_p=U_1U_2U_3^\dagger U_4^\dagger.
\]

At finite representation, \(W_p\) is a partial, generally weighted shift. Define its diagonal weighted-flippability operator

\[
D_p=\sqrt{W_p^\dagger W_p}+\sqrt{W_pW_p^\dagger}.
\]

The declared electromagnetic Hamiltonian is

\[
\boxed{
H_A(u,t,v)=
\frac{u}{2}\sum_\ell E_\ell^2
-t\sum_p(W_p+W_p^\dagger)
+v\sum_pD_p.
}
\]

The selected RK point is

\[
u=0,\qquad v=t.
\]

At that point,

\[
h_p=t(D_p-W_p-W_p^\dagger).
\]

If \(W_p\) connects two electric-basis configurations \(|C\rangle\) and \(|C'\rangle\) with positive amplitude \(a_{CC'}\), the corresponding contribution is

\[
t a_{CC'}
(|C\rangle-|C'\rangle)
(\langle C|-\langle C'|).
\]

Therefore every local term is positive semidefinite, every term commutes with every Gauss generator, and the equal-amplitude superposition within each connected flux sector is an exact zero-energy ground state.

For spin 1, the nonzero matrix elements of the normalized raising operator are exactly one, so the local term is an ordinary graph Laplacian. For spin 2, it is a weighted graph Laplacian. The fixed spin representation is the finite rishon-number choice; no separate determinant term is required for this U(1) model.

## 3. Exact finite three-dimensional check

The full zero-charge, zero-winding connected component was enumerated on a periodic \(2^3\) cube:

| Quantity | Result |
|---|---:|
| Sites | 8 |
| Oriented links | 24 |
| Plaquettes | 24 |
| Link representation | spin 1, dimension 3 |
| Gauge-reduced connected states | 146,327 |
| Undirected plaquette-transition edges | 1,236,144 |
| Minimum / maximum graph degree | 4 / 48 |
| Maximum sampled Gauss residual | 0 |
| Equal-amplitude RK residual | 0 |

This is an exact finite Hamiltonian construction, not a harmonic or rotor replacement. The small cube is an algebra and closure test; it is not used as the finite-size proof of the Coulomb phase.

Frozen output: `finite_em_exact_cube_results.json`.

## 4. Static Coulomb diagnostics

The finite-size runs sample the exact equal-amplitude RK ground-state support with plaquette moves. Noncontractible loop proposals are used only to sample distinct winding sectors; they do not modify the Hamiltonian.

### 4.1 Winding-sector response

For a three-dimensional Coulomb phase, a uniform winding flux \(W\) has Gaussian free energy

\[
F(W)-F(0)=\frac{\kappa W^2}{2L}+\cdots,
\]

so

\[
\frac{\langle W^2\rangle}{L}\to\frac{1}{\kappa},
\qquad
L\,[F(1)-F(0)]\to\frac{\kappa}{2}.
\]

For spin 1:

| \(L\) | \(\langle W^2\rangle/L\) | \(L\Delta F_{W=1}\) |
|---:|---:|---:|
| 4 | 0.82276 | 0.58358 |
| 6 | 0.80761 | 0.62864 |
| 8 | 0.85777 | 0.52474 |

The fractional range is 6.05% for \(\langle W^2\rangle/L\) and 17.95% for \(L\Delta F_{W=1}\). The latter is an **equal-time sector free energy**, not a dynamical photon gap.

For the spin-2 control, \(\langle W^2\rangle/L\) changes by 4.10% between \(L=6\) and \(L=8\). Its different absolute value is a finite-representation change in the effective stiffness, while the scaling law persists.

### 4.2 Transverse correlation tensor

With link-centered Fourier phases, exact Gauss law gives

\[
\widehat q_i E_i(\mathbf q)=0,
\qquad
\widehat q_i=2\sin(q_i/2).
\]

The measured equal-time tensor has the Coulomb form

\[
C_{ij}(\mathbf q)\propto
\delta_{ij}-\frac{\widehat q_i\widehat q_j}{\widehat q^2}.
\]

Across the spin-1 \(L=4,6,8\) runs:

- maximum longitudinal power fraction: \(5.25\times10^{-33}\);
- maximum mean split of the two axis-mode transverse eigenvalues: 0.219.

The longitudinal null is exact to numerical precision. The residual transverse split is sampling noise and finite-size anisotropy. This establishes two equal-time transverse components and excludes a static longitudinal component; it does **not** yet establish two `z=1` quantum branches.

### 4.3 Wilson-shift overlap

For a closed rectangle \(C\), the finite Wilson shift is the product of link-raising and link-lowering operators around the loop. In the equal-amplitude state, its overlap is the fraction of finite flux configurations on which the shift is allowed.

For every tested spin-1 size, the perimeter fit is preferred to the area fit. At \(L=8\), for example:

| Rectangle | Shift overlap |
|---|---:|
| \(1\times1\) | 0.3096 |
| \(1\times2\) | 0.1567 |
| \(1\times3\) | 0.1058 |
| \(1\times4\) | 0.0529 |
| \(2\times2\) | 0.0933 |
| \(2\times3\) | 0.0538 |

This is a finite-operator deconfinement diagnostic. The winding distribution is the complementary global flux or 't Hooft-sector response.

Frozen output: `finite_em_flux_results.json`.

## 5. Static charged sector

A fixed \(+1\) charge is placed at the origin and a labeled \(-1\) endpoint is allowed to move. Moving the endpoint changes exactly one adjacent finite link, so

\[
G_x=Q_x
\]

is maintained exactly, including recoil of the electric string. Plaquette moves sample the flux configurations at fixed endpoints, and the endpoint histogram measures the relative charged-sector partition function.

The fitted potential is

\[
V_L(\mathbf r)=c_L-A_LG_L(\mathbf r),
\]

where \(G_L\) is the periodic cubic-lattice Green function.

| \(L\) | Histogram samples | \(A_L\) | weighted \(R^2\) | effective linear-string slope |
|---:|---:|---:|---:|---:|
| 6 | 1,038,000 | 0.8440 | 0.4721 | 0.007553 |
| 8 | 3,278,400 | 0.6045 | 0.5264 | 0.003309 |

The sample counts include endpoint measurements between full flux sweeps.

The Coulomb amplitude is positive at both sizes, the lattice-Green model improves substantially over a constant no-force model, and the effective linear-string slope falls by 56% from \(L=6\) to \(L=8\) rather than approaching a nonzero string tension. On the larger lattice the Coulomb fit also beats the linear-string fit directly.

The finite-size evidence supports a static deconfined charge response. It is not yet a dynamical charged-particle sector: no matter hopping Hamiltonian, antiparticle content, fermion statistics, or vacuum polarization is included here.

Frozen output: `finite_em_charge_results.json`.

## 6. What Stage 6A establishes

The following gates now pass:

1. A declared 3+1D finite Hamiltonian with no continuum electromagnetic field inserted by hand.
2. Exact local U(1) Gauss symmetry.
3. A positive finite spin-1 RK parent and an exactly enumerated periodic three-dimensional gauge sector.
4. Coulomb winding-sector scaling over \(L=4,6,8\).
5. A rank-two equal-time transverse correlation tensor with no longitudinal component.
6. Perimeter-favored finite Wilson-shift overlaps.
7. A static opposite-charge response consistent with the periodic lattice Green function and a string slope that decreases with size.
8. Persistence of the static phase diagnostics in a larger finite representation.

This is sufficient to replace the prior status **“finite electromagnetic kinematics only”** with **“finite static Coulomb ground state demonstrated.”**

## 7. Why issue #6 remains open

The RK point is a solvable ground-state anchor, not yet the desired relativistic quantum electrodynamics phase. The following acceptance criteria are still unpassed:

- a stable interval of microscopic couplings, not only the RK surface;
- two and only two **dynamical** photon branches with
  \(\omega=c_\gamma|\widehat{\mathbf k}|+O(k^3a^2)\);
- a lowest photon gap closing as \(1/L\);
- exclusion of an additional gapless scalar in the dynamical spectrum;
- a gauge-covariant charged-matter hopping Hamiltonian;
- charge renormalization, vacuum polarization, and finite Ward identities;
- regulator-universal QED scaling and quantitative bounds on Lorentz-violating operators.

Accordingly, `issue_6_closed` remains `false` in the master result. Calling the equal-time rank-two tensor a completed photon spectrum would violate the architecture audit's claim thresholds.

## 8. Next falsifiable calculation

The next calculation is now narrow:

1. Scan \((u/t,v/t)\) on periodic three-dimensional clusters around the RK point.
2. Measure the lowest gauge-invariant spectral gaps in the zero-flux and one-flux sectors.
3. Require a finite coupling interval with
   \[
   \Delta_\gamma(L)\propto L^{-1}
   \]
   and two transverse spectral residues.
4. Measure a longitudinal/scalar rejection channel under the same fit.
5. Compare spin 1 and spin 2 before adding matter.
6. Only after that pass, add the first gauge-covariant charged defect and test the finite Ward identity.

The rejection rule is explicit: if the only gapless point is RK-like with nonrelativistic dynamics, or if the detuned spin-1 model confines or remains gapped, this branch is not the electromagnetic continuum of the theory and must be revised.

## 9. Reproduction

From the repository root:

```sh
# Integrated deterministic CI regression.
python phase_junction_network/microscopic/check_finite_em_coulomb.py \
  --quick \
  --output /tmp/finite_em_quick.json

# Exact 2^3 gauge-reduced cube. Run separately to release its large tuple arena.
python phase_junction_network/microscopic/check_finite_em_coulomb.py \
  --exact-only \
  --output /tmp/finite_em_exact_cube.json

# Frozen full flux/transverse finite-size run.
python phase_junction_network/microscopic/check_finite_em_coulomb.py \
  --flux-only \
  --output /tmp/finite_em_flux.json

# Frozen charged-sector runs. Each starts in a fresh process.
python phase_junction_network/microscopic/check_finite_em_coulomb.py \
  --charge-only 6 \
  --output /tmp/finite_em_charge_L6.json
python phase_junction_network/microscopic/check_finite_em_coulomb.py \
  --charge-only 8 \
  --output /tmp/finite_em_charge_L8.json
```

Every mode exits nonzero when its committed acceptance threshold fails.

## 10. Related finite-state literature

The architecture follows finite quantum-link and three-dimensional U(1)-liquid precedents while retaining separate project-specific acceptance gates:

- S. Chandrasekharan, *Confinement, Chiral Symmetry Breaking and Continuum Limits in Quantum Link Models*, arXiv:hep-lat/9809084.
- M. Hermele, M. P. A. Fisher, and L. Balents, *Pyrochlore Photons: The U(1) Spin Liquid in a S=1/2 Three-Dimensional Frustrated Magnet*, arXiv:cond-mat/0305401.
- O. Sikora et al., *A quantum liquid with deconfined fractional excitations in three dimensions*, arXiv:1105.1322.
- D. Banerjee, E. Huffman, and L. Rammelmüller, *Exploring Bosonic and Fermionic Link Models on (3+1)-D Tubes*, arXiv:2201.07171.

Those results motivate the phase target; they do not substitute for the finite spin-1 calculations and claim boundaries recorded here.
