# Finite Hilbert-Space Construction for the Phase Junction Network

**Status:** kinematic construction completed; dynamical graviton phase remains open  
**Date:** 2026-09-20

## 1. Objective

The preceding Phase Junction work identified infrared target theories:

- a compact U(1) link field for electromagnetism;
- a symmetric frame deformation with three vector constraints and one scalar curvature/energy constraint for gravity.

The next question is whether those targets can arise from a finite local Hilbert space. This note gives a concrete partial answer.

> Electromagnetic Gauss symmetry has an exact, minimal three-state quantum-link realization. The four linear gravitational constraints have an exact finite **discrete** realization with the correct local degree count. An exact finite **continuous canonical** realization is impossible on the full local Hilbert space, so continuous frame variables must emerge below a cutoff.

This is progress on the microscopic architecture, not yet a complete finite quantum-gravity Hamiltonian.

## 2. Exact finite electromagnetic link

Assign a spin-S Hilbert space to every oriented spatial link:

\[
E_\ell=S^z_\ell,
\qquad
U_\ell=\frac{S^+_\ell}{u_{\max}}.
\]

The link operator is not unitary at finite S, but it obeys the exact algebra

\[
[E_\ell,U_\ell]=U_\ell.
\]

At a vertex x define

\[
G_x=\sum_{\ell\,\mathrm{out\ of}\,x}E_\ell
-\sum_{\ell\,\mathrm{into}\,x}E_\ell-Q_x.
\]

Then

\[
e^{i\sum_x\lambda_xG_x}
U_{xy}
e^{-i\sum_x\lambda_xG_x}
=e^{i(\lambda_x-\lambda_y)}U_{xy}.
\]

Thus a finite link carries an exact continuous compact U(1) gauge representation. This is the quantum-link mechanism rather than a finite clock approximation.

The minimal Hamiltonian is

\[
H_A=\frac{U_A}{2}\sum_\ell E_\ell^2
-\frac{K_A}{2}\sum_p(U_p+U_p^\dagger),
\]

where

\[
U_p=U_1U_2U_3^\dagger U_4^\dagger.
\]

Every plaquette term commutes with every Gauss generator.

### 2.1 Why spin-1/2 is too small

For spin-1/2,

\[
E_\ell^2=\frac14I.
\]

The entire electric term is a constant. It can enforce neither an electric stiffness nor a tunable photon impedance. The plaquette still has gauge-invariant dynamics, but it does not reproduce the required rotor Hamiltonian by itself.

For spin-1,

\[
E_\ell\in\{-1,0,+1\},
\]

so \(E_\ell^2\) is nontrivial. A neutral single plaquette has three flux-loop states,

\[
|m,m,-m,-m\rangle,
\qquad m=-1,0,+1,
\]

and the ring operator steps between adjacent m values. Therefore

\[
\boxed{S_A=1,\quad d_A=3}
\]

is the smallest useful quantum-link realization tested here.

### 2.2 Exact single-plaquette check

With \(U_A=1\) and \(K_A=0.4\), the included calculation gives:

| link spin | link dimension | neutral states | distinct electric energies | max \([H,G_x]\) entry |
|---:|---:|---:|---:|---:|
| 1/2 | 2 | 2 | 1 | 0 |
| 1 | 3 | 3 | 2 | 0 |
| 3/2 | 4 | 4 | 2 | 0 |
| 2 | 5 | 5 | 3 | 0 |

The result does not prove a three-dimensional Coulomb phase at S=1. It establishes the smallest local representation that can possibly reproduce both terms of the proposed electromagnetic Hamiltonian.

## 3. Obstruction to an exact finite canonical gravity variable

The infrared gravitational construction uses canonical variables

\[
[h_{ij},\pi^{kl}]
=i\,\delta_{(i}^{k}\delta_{j)}^{l}.
\]

No finite-dimensional matrices Q and P can satisfy

\[
[Q,P]=iI
\]

on their full Hilbert space. Taking the trace gives

\[
\operatorname{tr}[Q,P]=0,
\]

but

\[
\operatorname{tr}(iI)=id\ne0.
\]

This is a structural obstruction, not a failed parameter choice. Consequently, the continuum constraint algebra in `frame_gravity_derivation.md` cannot be copied literally into finite local matrices while retaining an exact canonical pair everywhere in the spectrum.

There are two viable responses.

1. Use an exact finite **discrete** Weyl algebra and ask for the continuous frame symmetry to emerge at long distance.
2. Use a finite truncation or collective spin whose canonical algebra is accurate only below a boundary cutoff.

Both are implemented below.

## 4. Exact discrete gravity constraint skeleton

Choose an odd prime p. Put six p-level qudits at each site, corresponding to

\[
(xx,yy,zz,xy,xz,yz).
\]

Each qudit has Weyl operators

\[
ZX=\omega XZ,
\qquad
\omega=e^{2\pi i/p}.
\]

Use canonical off-diagonal momentum variables

\[
p_{ij}=2\pi_{ij}\quad(i<j),
\]

so every listed coordinate has one canonical momentum.

Let \(D_i\) be commuting periodic central differences over \(\mathbb F_p\). The three X-type vector constraints are

\[
\begin{aligned}
G_x&=D_xp_{xx}+\frac12D_yp_{xy}+\frac12D_zp_{xz},\\
G_y&=D_yp_{yy}+\frac12D_xp_{xy}+\frac12D_zp_{yz},\\
G_z&=D_zp_{zz}+\frac12D_xp_{xz}+\frac12D_yp_{yz}.
\end{aligned}
\]

The Z-type scalar constraint is

\[
\begin{aligned}
C={}&-(D_y^2+D_z^2)h_{xx}
-(D_x^2+D_z^2)h_{yy}
-(D_x^2+D_y^2)h_{zz}\\
&+2D_xD_yh_{xy}
+2D_xD_zh_{xz}
+2D_yD_zh_{yz}.
\end{aligned}
\]

This is the discrete form of

\[
C=(\partial_i\partial_j-\delta_{ij}\partial^2)h_{ij}.
\]

### 4.1 Exact commutation

The Weyl commutation phase is determined by the modular symplectic product. For example, the commutator with \(G_x\) is proportional to

\[
-(D_y^2+D_z^2)D_x
+(2D_xD_y)\frac12D_y
+(2D_xD_z)\frac12D_z=0.
\]

The same identity holds for y and z. Therefore

\[
\boxed{[C(\mathbf x),G_i(\mathbf y)]=0}
\]

exactly as a finite Weyl-operator statement, not merely in a large-p expansion.

The code verifies that every entry of the scalar/vector symplectic commutator matrix is zero modulo p.

### 4.2 Degree count

For \(N=L^3\) periodic sites, with odd L and a prime p chosen away from modular degeneracies, the measured ranks are

\[
\operatorname{rank}C=N-1,
\qquad
\operatorname{rank}G=3N-3.
\]

The four missing ranks are the summed periodic redundancies: one scalar and one for each vector component.

Starting from \(6N\) qudits, the number of logical qudits is

\[
\begin{aligned}
k&=6N-(N-1)-(3N-3)\\
&=2N+4.
\end{aligned}
\]

Thus there are

\[
\boxed{2\text{ local logical modes per site}+4\text{ global modes}.}
\]

The included exact checks give:

| L | p | sites N | scalar rank | vector rank | logical qudits | nonzero commutators |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 5 | 27 | 26 | 78 | 58 = 2N+4 | 0 |
| 5 | 7 | 125 | 124 | 372 | 254 = 2N+4 | 0 |

This solves the finite-dimensional **kinematic** constraint problem.

### 4.3 What the discrete code does not prove

A Hamiltonian made only from commuting constraint projectors is generically gapped. The degree count alone does not establish:

- helicity labels;
- a gapless phase;
- linear dispersion;
- the Fierz-Pauli stiffness;
- nonlinear graviton self-coupling.

The code is therefore a finite constraint skeleton, not yet a finite graviton theory.

## 5. Emergent continuous canonical variables

### 5.1 Truncated oscillator

In a d-level oscillator truncation,

\[
a|n\rangle=\sqrt n|n-1\rangle,
\]

with the raising operation terminated at \(|d-1\rangle\). Define

\[
Q=\frac{a+a^\dagger}{\sqrt2},
\qquad
P=\frac{a-a^\dagger}{i\sqrt2}.
\]

Then

\[
\boxed{[Q,P]=i\left(I-d|d-1\rangle\langle d-1|\right).}
\]

The obstruction is concentrated at the top state. Matrix elements entirely inside the subspace

\[
\{|0\rangle,\ldots,|d-2\rangle\}
\]

have the exact canonical commutator. A large boundary penalty can therefore protect a low-energy canonical sector, although Q and P can still connect the penultimate state to the boundary and the full dynamics must be checked.

The numerical diagnostic finds low-subspace operator error below \(6\times10^{-16}\) for d=3,4,5, while the full-space defect has norm d.

### 5.2 Collective spin made from microscopic hinges

For a collective spin S, define

\[
Q=\frac{S_x}{\sqrt S},
\qquad
P=\frac{S_y}{\sqrt S}.
\]

Then

\[
[Q,P]=i\frac{S_z}{S}.
\]

Near the polarized state, with n spin-wave excitations,

\[
\frac{S_z}{S}=1-\frac nS.
\]

The relative canonical error is therefore

\[
\boxed{\epsilon_n=\frac nS.}
\]

This is directly compatible with the original idea that many two-state matter/geometry hinges form one collective frame variable. It also quantifies the cost: a one-excitation sector needs \(S\ge100\) for a one-percent canonical error, corresponding to at least 201 states in the symmetric collective-spin representation.

## 6. The remaining dynamics obstruction

The next model must add local operators that preserve or energetically enforce the four constraints and produce

\[
\omega^2=c_g^2k^2+O(k^4\ell^2)
\]

for exactly two modes.

### 6.1 Exact derivative-order search

The accompanying `search_local_dynamics.py` solves the polynomial invariance equations without assuming a particular curvature basis.

For a coordinate observable

\[
O_h=L^{ij}(\mathbf k)h_{ij},
\]

vector-gauge invariance requires

\[
L^{ij}(\mathbf k)(k_i\xi_j+k_j\xi_i)=0
\]

for every \(\mathbf k\) and \(\boldsymbol\xi\). A homogeneous polynomial ansatz gives:

| derivative order of \(L\) | unknown coefficients | rank | nullity |
|---:|---:|---:|---:|
| 0 | 6 | 6 | 0 |
| 1 | 18 | 18 | 0 |
| 2 | 36 | 30 | 6 |

Thus a nonzero manifestly local coordinate invariant first appears at **two derivatives**. The six degree-two solutions are the expected linearized-curvature/Einstein-tensor family.

For a momentum observable

\[
O_p=M^{ij}(\mathbf k)p_{ij},
\]

invariance under the scalar momentum shift

\[
\delta p_{ij}=(k_ik_j-\delta_{ij}k^2)\eta
\]

first appears at one derivative:

| derivative order of \(M\) | unknown coefficients | rank | nullity |
|---:|---:|---:|---:|
| 0 | 6 | 6 | 0 |
| 1 | 18 | 10 | 8 |

The ranks were reproduced over two independent large prime fields, eliminating accidental floating-point nullspaces.

### 6.2 Consequence for a manifest-local square Hamiltonian

A positive Hamiltonian made from squares of the lowest-order exact local invariants has

\[
H_p\sim k^2p^2,
\qquad
H_h\sim k^4h^2.
\]

Hamilton's equations then give

\[
\dot h\sim k^2p,
\qquad
\dot p\sim-k^4h,
\]

so

\[
\boxed{\omega^2\sim k^6,\qquad\omega\sim k^3.}
\]

This recovers, from the constraint algebra alone, the cubic-dispersion tendency found in the manifestly local L-type qubit construction. It is stronger than merely observing that \(R^2\) has four derivatives: exact scalar-momentum gauge invariance also forces at least one derivative into the local kinetic observable.

### 6.3 What can still produce linear dispersion

The result is not a no-go for a finite emergent graviton. It rules out the simplest architecture in which the Hamiltonian is a sum of positive squares of individually local exact gauge invariants.

A linear branch requires at least one of the following:

1. a globally gauge-invariant Fierz–Pauli bilinear \(hMh\) whose local densities change by a discrete divergence;
2. an area/tetrad factor multiplying curvature so that the action is linear rather than quadratic in curvature;
3. a protected truncated-oscillator Hamiltonian where the continuous symmetry is emergent rather than exact on the full finite Hilbert space;
4. a strongly interacting phase in which one conjugate sector is dynamically frozen and the low-energy mode reorganizes;
5. auxiliary link/frame variables that localize the global bilinear without converting it into curvature squared.

The two-derivative Fierz–Pauli energy is invariant as a full periodic sum, but it is not a sum of individually gauge-invariant commuting projectors. A successful finite model must therefore be noncommuting, use additional frame structure, or realize the symmetry only in its low-energy sector.

This is consistent with the existing qubit-model literature: finite local models with emergent helicity-2 modes have been constructed, but the linearly dispersing candidate is strongly interacting and was not established by a controlled solution. The Phase Junction project should not treat the finite constraint count as having solved that gate.

## 7. Why the impedance ratio is still free

The finite symmetries determine allowed operators but not their coefficients. At the present stage the low-energy Hamiltonian still contains

\[
U_A,\ K_A,\ U_g,\ K_g.
\]

A common causal cone imposes

\[
U_AK_A=U_gK_g,
\]

but it does not determine

\[
\frac{Z_g}{Z_A}
=\sqrt{\frac{U_g/K_g}{U_A/K_A}}.
\]

The constraint architecture alone therefore cannot derive \(G\) and \(\alpha\). Equal impedance remains an extra hypothesis unless one local move set or an exact exchange symmetry calculates both ratios.

The next microscopic Hamiltonian must expose its elementary amplitudes and energy penalties and derive all four stiffnesses by perturbation theory or direct spectrum matching. Parameter counting should reject any construction that merely renames four independent couplings.

## 8. Next Hamiltonian search

The next concrete search should use the finite gravity code as its constrained Hilbert space and enumerate short-range operators in the two symplectic kernels:

- Z-type coordinate operators z satisfying \(Gz=0\);
- X-type momentum operators x satisfying \(Cx=0\).

For every local candidate pair:

1. verify exact commutation with all constraints;
2. determine whether the terms commute with each other or can support a gapless phase;
3. calculate the small-momentum quadratic kernel in a large-p or large-S expansion;
4. reject k-squared, k-cubed, unstable, scalar-contaminated, or anisotropic branches;
5. retain only candidates with two degenerate positive modes and \(\omega\propto k\);
6. derive their coefficients from the same local amplitudes used by the electromagnetic ring term.

The primary acceptance condition is no longer merely “two modes remain.” It is

\[
\boxed{
\text{finite local Hilbert space}
+\text{four constraints}
+\text{two helicities}
+\omega\propto k
+\text{shared microscopic coefficients}.
}
\]

## 9. Current decision

The finite-state phase has produced three hard results:

1. \(d_A=3\) is the smallest useful electromagnetic link dimension in the tested quantum-link family.
2. A finite odd-prime qudit lattice can impose the four gravitational constraints exactly and leaves two local modes per site, plus four periodic global modes.
3. Exact continuous canonical frame variables cannot exist on the full finite local Hilbert space; they must be emergent below a cutoff.

The unresolved centerpiece is now sharply defined: construct and verify local noncommuting dynamics on that finite skeleton that realizes a linear helicity-2 phase and simultaneously fixes the electromagnetic/gravitational stiffness ratio.

## Primary literature used for comparison

- S. Chandrasekharan and U.-J. Wiese, “Quantum Link Models: A Discrete Approach to Gauge Theories,” arXiv:hep-lat/9609042.
- U.-J. Wiese, “From Quantum Link Models to D-Theory: A Resource Efficient Framework for the Quantum Simulation and Computation of Gauge Theories,” arXiv:2107.09335.
- Z.-C. Gu and X.-G. Wen, “Emergence of helicity ±2 modes (gravitons) from qbit models,” arXiv:0907.1203.
- A. T. Schmitz, “Gauge Structures: From Stabilizer Codes to Continuum Models,” arXiv:1809.10151.
