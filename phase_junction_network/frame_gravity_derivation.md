# Phase-Junction Frame Gravity

**Status:** zero-data candidate solution for the gravitational sector of the Phase Junction Network  
**Date:** 2026-09-20  
**Scope:** derive a healthy long-range tensor field, a Newtonian static limit, universal coupling, and a common electromagnetic/gravitational causal cone without fitting observations.

## 1. Result in one paragraph

A scalar matter–gravity phase cannot serve as gravity any more than its gradient could serve as electromagnetism. The smallest viable gravitational junction variable is a **local frame comparison**. At the linear level its physical coordinate is a symmetric rank-2 field \(h_{ij}\), with conjugate junction imbalance \(\pi^{ij}\). Three local momentum-balance constraints and one local curvature/energy constraint are required. They form four first-class constraints and reduce the six symmetric components to two propagating transverse-traceless modes. The unique local two-derivative quadratic stiffness compatible with the relabeling symmetry is the Fierz–Pauli stiffness. Its physical spectrum is positive and linear, its static Green function is \(1/r\), and coupling every excitation through one common frame yields universal gravity. A frame-holonomy or tetrad/connection implementation supplies the nonlinear completion. The remaining unsolved problem is microscopic parameter closure: the framework relates \(G\), \(\alpha\), the network spacing, and defect mass gaps, but does not yet calculate their dimensionless values.

## 2. Junction variables

At each spatial cell use a symmetric geometric phase/strain variable

\[
h_{ij}=h_{ji},
\]

and a conjugate geometry-channel imbalance

\[
\pi^{ij}=\pi^{ji},
\qquad
[h_{ij}(\mathbf n),\pi^{kl}(\mathbf m)]
=i\,\delta_{(i}^{k}\delta_{j)}^{l}\delta_{\mathbf n\mathbf m}.
\]

A microscopic interpretation is obtained from a link frame \(e_i^{\ a}\):

\[
e_i^{\ a}=\ell\left(\delta_i^{\ a}+\frac12 h_i^{\ a}+\cdots\right),
\qquad
h_{ij}=h_{ji}=\frac12\left(h_i^{\ a}\delta_{aj}+h_j^{\ a}\delta_{ai}\right).
\]

The antisymmetric part of \(h_i^{\ a}\) changes the local frame orientation and is not a physical strain.

## 3. The four indispensable constraints

Let \(D_i\) denote a lattice derivative, or \(\partial_i\) in the continuum.

### 3.1 Vector or momentum constraint

\[
\mathcal G_i=-2D_j\pi_i^{\ j}-P_i=0.
\]

In vacuum, \(P_i=0\). This generates local spatial relabelings:

\[
\delta_\xi h_{ij}=D_i\xi_j+D_j\xi_i.
\]

Physically, geometry imbalance cannot begin or end at a cell unless matter carries the corresponding momentum.

### 3.2 Scalar curvature/energy constraint

\[
\mathcal C=D_iD_jh_{ij}-D^2h-\gamma\rho=0,
\qquad h=h_i^{\ i}.
\]

In vacuum, \(\rho=0\). Its canonical transformation shifts the momentum by

\[
\delta_\eta\pi^{ij}
=\left(D^iD^j-\delta^{ij}D^2\right)\eta.
\]

This is the clock-like or temporal junction redundancy. It is not an additional radiative scalar. It removes the scalar mode that would otherwise have negative quadratic energy.

### 3.3 Constraint algebra

At the linear level,

\[
\{\mathcal G_i,\mathcal G_j\}
=\{\mathcal G_i,\mathcal C\}
=\{\mathcal C,\mathcal C\}=0.
\]

The constraints are preserved by the Hamiltonian. In particular,

\[
\dot{\mathcal C}=-\frac{U_g}{2}D_i\mathcal G_i.
\]

Thus the scalar and vector constraints are parts of one consistent local conservation structure.

## 4. Unique quadratic junction Hamiltonian

The minimal local Hamiltonian is

\[
\begin{aligned}
H_g^{(2)}={}&\frac{U_g}{2}\sum_{\mathbf n}
\left(\pi_{ij}\pi_{ij}-\frac12\pi^2\right)\\
&+\frac{K_g}{2}\sum_{\mathbf n}
\Big[
(D_kh_{ij})^2
-2(D_ih_{ij})^2
+2(D_ih_{ij})(D_jh)
-(D_kh)^2
\Big]\\
&+\sum_{\mathbf n}\left(N\mathcal C+N^i\mathcal G_i\right),
\end{aligned}
\]

where \(\pi=\pi_i^{\ i}\). The multipliers \(N\) and \(N^i\) are the clock and shift junction variables. They enforce the four constraints and do not add propagating modes.

The apparent negative trace direction in the unconstrained kinetic form is gauge, not physical. On the constrained transverse-traceless subspace the Hamiltonian is positive.

## 5. Degree-of-freedom proof

A symmetric spatial tensor has six coordinates and six conjugate momenta, for 12 phase-space variables.

There are four first-class constraints:

\[
3\times\mathcal G_i+1\times\mathcal C.
\]

Each first-class constraint removes one constraint direction and one gauge direction. Therefore,

\[
12-2(4)=4
\]

physical phase-space dimensions remain, corresponding to

\[
\boxed{2\text{ propagating configuration modes}.}
\]

For a wavevector along \(z\):

- \(\mathcal G_i=0\) removes \(\pi_{zz},\pi_{xz},\pi_{yz}\).
- Spatial relabeling removes \(h_{zz},h_{xz},h_{yz}\).
- \(\mathcal C=0\) imposes \(h_{xx}+h_{yy}=0\).
- The scalar gauge direction removes \(\pi_{xx}+\pi_{yy}\).

The two remaining pairs are

\[
h_+=\frac{h_{xx}-h_{yy}}{\sqrt2},
\qquad
h_\times=\sqrt2h_{xy},
\]

and their conjugate momenta.

## 6. Why the fourth constraint is mandatory

With only the three vector constraints, the transverse symmetric sector contains three modes. The stiffness eigenvalues are

\[
\{-k^2,+k^2,+k^2\}.
\]

The negative eigenvalue is the transverse trace or scalar breathing mode. The scalar curvature constraint removes exactly this mode. It is therefore the mechanism that makes the tensor theory healthy.

## 7. Wave spectrum

On the physical transverse-traceless subspace,

\[
H_{TT}=\frac{U_g}{2}\pi_{TT}^{ij}\pi_{ij}^{TT}
+\frac{K_g}{2}\widehat{k}^{\,2}h_{TT}^{ij}h_{ij}^{TT},
\]

with cubic-lattice momentum

\[
\widehat{k}^{\,2}
=4\sum_{a=1}^{3}\sin^2\left(\frac{k_a\ell}{2}\right).
\]

The dispersion is

\[
\omega_g^2(\mathbf k)
=\frac{U_gK_g}{\hbar^2}\widehat{k}^{\,2}.
\]

At long wavelength,

\[
\omega_g^2=c_g^2k^2+O(k^4\ell^2),
\qquad
\boxed{c_g=\frac{\ell\sqrt{U_gK_g}}{\hbar}}.
\]

Both tensor polarizations have the same positive eigenvalue.

## 8. Continuum normalization and Newtonian coupling

Eliminating \(\pi^{ij}\) and taking the continuum limit gives the transverse action

\[
S_{TT}
=\frac{M_*^2}{8}\int d^4x
\left[(\partial_t h_{ij}^{TT})^2-(\nabla h_{ij}^{TT})^2\right]
\]

in units \(\hbar=c=1\), with

\[
\boxed{M_*^2=\frac{4}{U_g\ell^3}=\frac{4K_g}{\ell}}
\]

when \(c_g=1\).

Coupling the same field to total stress-energy gives

\[
M_*^2G_{\mu\nu}^{(1)}=T_{\mu\nu}.
\]

The static weak-field equation is

\[
\nabla^2\Phi=4\pi G\rho_m,
\]

where

\[
\boxed{G=\frac{1}{8\pi M_*^2}}
\]

in natural units. Restoring \(\hbar\) and \(c\),

\[
\boxed{
G=\frac{c^4\ell}{32\pi K_g}
=\frac{U_gc^2\ell^3}{32\pi\hbar^2}.
}
\]

Therefore a point source gives

\[
\Phi(r)=-\frac{GM}{r},
\qquad
V(r)=-\frac{Gm_1m_2}{r}.
\]

## 9. Universal coupling

All matter fields must use the same link displacement/frame variable \(e_\mu^{\ a}\). This is the **single-soldering rule**:

\[
S_m=S_m[e,\omega,A,\Psi_s]
\]

for every species \(s\). A separate geometry for each species would be a different theory and would generically violate universality.

At the linear level,

\[
S_{\rm int}=\frac12\int d^4x\,h_{\mu\nu}T^{\mu\nu}.
\]

Under

\[
\delta h_{\mu\nu}=\partial_\mu\xi_\nu+\partial_\nu\xi_\mu,
\]

the variation is proportional to

\[
\partial_\mu T^{\mu\nu}.
\]

Thus the field couples to conserved total energy-momentum. When interacting sectors exchange energy, separate species-dependent gravitational coefficients are not compatible with one conserved total source. The geometry field's own energy must also enter the source at nonlinear order.

## 10. Frame-holonomy microscopic implementation

The linear symmetric-tensor theory has a natural junction-network realization.

Each oriented link carries:

- a displacement or tetrad-like vector \(e_\ell^{\ a}\);
- a local frame comparison \(\Lambda_\ell\in \mathrm{Spin}(1,3)\);
- the electromagnetic phase comparison \(U_\ell=e^{ia_\ell}\).

Around a plaquette or hinge,

\[
U_f=\prod_{\ell\in\partial f}U_\ell=e^{iF_f},
\]

\[
\Lambda_f=\prod_{\ell\in\partial f}\Lambda_\ell
=e^{R_f^{ab}J_{ab}}.
\]

The oriented area bivector is

\[
B_f^{ab}=e_{\ell_1}^{[a}e_{\ell_2}^{b]}.
\]

The minimal gravitational loop term is linear in frame curvature and weighted by area:

\[
S_g^{\rm disc}
=\beta_g\sum_f
\epsilon_{abcd}B_f^{ab}R_f^{cd}
+S_{\rm constraints}.
\]

In the continuum this becomes

\[
\boxed{
S_g=\frac{M_*^2}{4}
\int\epsilon_{abcd}\,e^a\wedge e^b\wedge R^{cd}(\omega).
}
\]

Variation with respect to the frame connection sets torsion to zero in the spinless vacuum. Eliminating that connection gives the metric two-derivative action whose quadratic expansion is the Hamiltonian above.

The important Phase Junction interpretation is:

- electromagnetic field = scalar phase holonomy;
- gravity = frame/clock holonomy;
- matter lives at endpoints and defects of the same network;
- the time-like junction variables enforce the scalar and vector constraints rather than creating extra polarizations.

## 11. Common speed of light and gravity

The electromagnetic action is built using the same tetrad:

\[
S_A=-\frac{1}{4g_A^2}
\int d^4x\sqrt{-g}\,
F_{\mu\nu}F^{\mu\nu},
\qquad
g_{\mu\nu}=e_\mu^{\ a}e_{\nu a}.
\]

Because both the photon and tensor mode use the same effective metric, their long-wavelength characteristic cone is the same:

\[
\boxed{c_g=c_\gamma=c}
\]

at leading two-derivative order. This is stronger than tuning two unrelated speeds to match.

In lattice parameters the same condition is

\[
U_AK_A=U_gK_g=\left(\frac{\hbar c}{\ell}\right)^2.
\]

## 12. Relation among \(G\), \(\alpha\), and junction impedances

For the electromagnetic junction sector,

\[
c=\frac{\ell\sqrt{U_AK_A}}{\hbar},
\qquad
\alpha=\frac{1}{4\pi}\sqrt{\frac{U_A}{K_A}}.
\]

Define dimensionless junction impedances

\[
Z_A=\sqrt{\frac{U_A}{K_A}},
\qquad
Z_g=\sqrt{\frac{U_g}{K_g}}.
\]

Then

\[
\boxed{Z_A=4\pi\alpha}
\]

and, using \(\ell_P^2=G\hbar/c^3\),

\[
\boxed{Z_g=32\pi\frac{\ell_P^2}{\ell^2}}.
\]

A genuine microscopic unification must calculate \(Z_A\) and \(Z_g\). The simplest optional closure is a universal junction impedance,

\[
Z_A=Z_g.
\]

That gives

\[
\boxed{
\ell=\sqrt{\frac{8}{\alpha}}\,\ell_P
\simeq33.11024\,\ell_P.
}
\]

This is a candidate closure, not yet a derived fact. A more realistic microscopic calculation may produce a fixed ratio \(Z_g/Z_A\) rather than equality.

## 13. Particle masses

A relativistic junction defect can have lattice dispersion

\[
E_s^2(\mathbf k)
=\Delta_s^2
+\left(\frac{\hbar c}{\ell}\right)^2
\sum_i\sin^2(k_i\ell),
\]

so at long wavelength

\[
E_s^2=m_s^2c^4+\hbar^2c^2k^2,
\qquad
\boxed{m_sc^2=\Delta_s}.
\]

Writing

\[
\Delta_s=\mu_s\frac{\hbar c}{\ell}
\]

gives

\[
\boxed{m_s=\mu_s\frac{\hbar}{c\ell}}.
\]

The dimensionless defect gaps \(\mu_s\) must be obtained from junction topology, symmetry breaking, or collective dynamics. They are not fixed by the tensor spectrum. Under the equal-impedance illustrative closure, the network energy scale is approximately \(3.69\times10^{17}\,\mathrm{GeV}\); the electron would require \(\mu_e\approx1.39\times10^{-21}\). The particle-mass hierarchy is therefore a major unsolved part of the theory, not a minor parameter fit.

## 14. Prototype verification

The accompanying script performs the following no-data checks:

1. Draw 500 random nonzero wavevectors.
2. Construct the three vector gauge directions and scalar curvature constraint.
3. Verify that the scalar constraint annihilates spatial gauge directions.
4. Verify that its momentum-space gauge direction satisfies the vector constraints.
5. Project the stiffness onto the transverse-traceless subspace.
6. Confirm two equal positive eigenvalues \(k^2,k^2\).
7. Remove the scalar constraint and confirm the spectrum \(-k^2,+k^2,+k^2\).
8. Solve a point-source lattice Poisson problem on a \(96^3\) periodic grid and fit the radial potential to \(A/r+B\).

Results:

- spectrum failures: **0 / 500**;
- maximum relative TT-eigenvalue error: **\(1.13\times10^{-15}\)**;
- maximum constraint/gauge residual: **\(1.09\times10^{-14}\)**;
- static \(A/r+B\) normalized RMS residual over lattice radii 4–20: **0.00241**.

These are algebraic and synthetic checks only. They do not establish that nature uses the model.

## 15. What is solved and what is not

### Solved at the linear structural level

- exactly two tensor polarizations;
- no physical negative-energy scalar;
- linear long-wavelength dispersion;
- a \(1/r\) static field and \(1/r^2\) force;
- one universal source under the single-soldering rule;
- a structural common photon/gravity cone;
- explicit formulas mapping \(G\) and \(\alpha\) to junction parameters.

### Not yet solved

- a finite-state microscopic Hamiltonian that demonstrably flows to this exact frame-holonomy action;
- calculation of the impedance ratio \(Z_g/Z_A\);
- derivation of fermion generations, chirality, charge assignments, and defect gaps;
- the enormous particle-mass hierarchy;
- nonlinear strong-field solutions and quantum consistency;
- a mechanism controlling the allowed vacuum-volume term;
- an explicit accounting of how the microscopic network evades the assumptions behind the Weinberg–Witten limitation;
- whether the same microscopic model also produces the previously explored long-memory or galaxy-scale behavior without adding an independent phenomenological field.

## 16. Decision

The viable gravitational branch is therefore:

\[
\boxed{
\text{scalar phase junction}
\;\longrightarrow\;
\text{frame-valued junction network}
\;\longrightarrow\;
4\text{ local constraints}
\;\longrightarrow\;
2\text{ tensor modes}.
}
\]

The low-energy gravitational action is not freely selectable once locality, one shared frame, two derivatives, and the absence of extra modes are imposed. Those requirements force the Phase Junction network into the same low-energy tensor universality class as linearized Einstein gravity. The open scientific content is the microscopic junction origin, the coupling relations, the vacuum-energy/volume term, and any calculable corrections beyond that infrared limit.

Two additional gates must be made explicit. First, the same local symmetries permit a volume term proportional to \(\epsilon_{abcd}e^a\wedge e^b\wedge e^c\wedge e^d\); a balanced flat junction vacuum must explain why its effective coefficient is absent or extremely small rather than simply deleting it. Second, an emergent massless spin-2 mode must evade the assumptions of the Weinberg–Witten limitation. A fundamentally discrete junction network can do so only if exact Lorentz covariance and a local gauge-invariant stress tensor are emergent rather than microscopic; this must be demonstrated in the finite-state model, not presumed.

## References used for comparison

- S. Weinberg, *Photons and Gravitons in S-Matrix Theory: Derivation of Charge Conservation and Equality of Gravitational and Inertial Mass*, Physical Review 135, B1049 (1964).
- S. Deser, *Self-Interaction and Gauge Invariance*, General Relativity and Gravitation 1, 9–18 (1970); arXiv:gr-qc/0411023.
- T. Regge, *General Relativity Without Coordinates*, Il Nuovo Cimento 19, 558–571 (1961).
- Z.-C. Gu and X.-G. Wen, *Emergence of helicity ±2 modes (gravitons) from qbit models*, Nuclear Physics B 863, 90–129 (2012); arXiv:0907.1203.
- C. Xu, *Gapless Bosonic Excitation without Symmetry Breaking: An Algebraic Spin Liquid with Soft Gravitons*, Physical Review B 74, 224433 (2006); arXiv:cond-mat/0609595.
- S. Weinberg and E. Witten, *Limits on Massless Particles*, Physics Letters B 96, 59–62 (1980).
- NIST, *2022 CODATA Recommended Values of the Fundamental Physical Constants*.
