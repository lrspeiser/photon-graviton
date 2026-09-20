# Phase Junction Network: Electromagnetic Derivation

**Status:** corrected zero-data construction for the electromagnetic sector  
**Date:** 2026-09-20  
**Scope:** replace the failed identification of a single scalar phase gradient with a compact link connection whose loop curvature supports two transverse photon modes.

## 1. Starting proposal and the necessary correction

The motivating proposal assigns matter-like and geometry-like channels to each microscopic junction,

\[
|\chi\rangle=\cos\Theta\,|m\rangle+e^{i\phi}\sin\Theta\,|g\rangle.
\]

The states \(|m\rangle\) and \(|g\rangle\) cannot be interpreted as a literal electron and literal graviton, because an ordinary superposition would mix charge, spin, and particle statistics. In the candidate model they are two internal storage channels of a deeper junction.

The original electromagnetic identification,

\[
A_\mu=\frac{\hbar}{q}\partial_\mu\phi,
\]

cannot describe ordinary light. For every smooth, single-valued phase,

\[
F_{\mu\nu}=\partial_\mu A_\nu-\partial_\nu A_\mu
=\frac{\hbar}{q}(\partial_\mu\partial_\nu-\partial_\nu\partial_\mu)\phi=0.
\]

It is a pure gauge. Singular winding can produce quantized holonomy or a vortex, but that is not a general propagating electromagnetic field and its winding number is not photon occupation number.

The correction is:

> The local relative phase may be a physical junction coordinate, but electromagnetism is the curvature of the **network used to compare phases at different junctions**.

## 2. A real junction exchange

A time-dependent relative phase alone does not transfer occupation. With

\[
|\chi(t)\rangle=\cos\Theta\,|m\rangle+e^{i\Omega t}\sin\Theta\,|g\rangle,
\]

the probabilities \(\cos^2\Theta\) and \(\sin^2\Theta\) remain constant.

Introduce an off-diagonal microscopic Hamiltonian,

\[
H_h=\frac{\delta}{2}\sigma_z-t\sigma_x.
\]

For the imbalance operator

\[
\hat E=\frac12(|m\rangle\langle m|-|g\rangle\langle g|),
\]

the coherent state has

\[
\langle\hat E\rangle=\frac12\cos2\Theta,
\]

and its equation of motion contains an exchange current proportional to

\[
\hbar\frac{d}{dt}\langle\hat E\rangle
\propto-t\sin2\Theta\sin\phi.
\]

Thus the phase can drive actual transfer only through a coupling that mixes the two internal channels.

## 3. From microscopic hinges to a link rotor

One two-state junction has a finite Hilbert space and cannot by itself support an arbitrarily occupied bosonic field. Put many microscopic hinges on every oriented link \(\ell\):

\[
E_\ell=\frac12\sum_{r=1}^{N}\sigma^z_{\ell r}.
\]

A collective coherent link state is schematically

\[
|\chi_\ell\rangle=\bigotimes_{r=1}^{N}
\left(\cos\Theta_\ell|m\rangle_r+e^{ia_\ell}\sin\Theta_\ell|g\rangle_r\right).
\]

Then

\[
\langle E_\ell\rangle=\frac{N}{2}\cos2\Theta_\ell.
\]

Near the balanced state \(\Theta_\ell=\pi/4+\delta\Theta_\ell\),

\[
\langle E_\ell\rangle\simeq-N\delta\Theta_\ell.
\]

In the large-\(N\) rotor limit, the collective compact phase and imbalance obey

\[
[a_\ell,E_{\ell'}]=i\delta_{\ell\ell'},
\qquad
e^{ia_\ell}|E_\ell\rangle=|E_\ell+1\rangle.
\]

The original mixing angle is therefore reinterpreted as the local electric-flux or channel-imbalance coordinate; it is not by itself the electromagnetic coupling constant.

## 4. Link phase, gauge redundancy, and loop curvature

Let \(a_{xy}\) compare the phase convention at neighboring sites \(x\) and \(y\). Under a local redefinition,

\[
a_{xy}\rightarrow a_{xy}+\lambda_x-\lambda_y.
\]

An individual link phase depends on convention. The oriented sum around a plaquette \(p\),

\[
B_p=\sum_{\ell\in\partial p}s_{p\ell}a_\ell,
\]

is invariant because all vertex phases cancel. In the continuum,

\[
B_p\simeq \ell^2F_{ij}.
\]

This loop nonclosure is the electromagnetic curvature. It is nonzero precisely when neighboring phase comparisons cannot be integrated into one global scalar phase.

## 5. Gauss constraint and charge

Define

\[
G_x=\sum_{\ell\,{\rm out\ of}\,x}E_\ell
-\sum_{\ell\,{\rm into}\,x}E_\ell-Q_x.
\]

Physical low-energy states satisfy

\[
G_x|\Psi\rangle=0.
\]

This is the lattice Gauss law. In the continuum it becomes

\[
\nabla\cdot\mathbf E=\rho.
\]

The interpretation is:

- \(E_\ell\) is signed matter–geometry imbalance carried by a link;
- a neutral junction cannot create or destroy net imbalance;
- charged matter is an endpoint or source of imbalance flux;
- closed virtual transfers preserve every local constraint.

## 6. Minimal gauge Hamiltonian

The smallest compact Hamiltonian consistent with the local constraint is

\[
\boxed{
H_A=\frac{U_A}{2}\sum_\ell E_\ell^2-K_A\sum_p\cos B_p.
}
\]

The first term is the cost of channel imbalance. The second rewards phase comparisons that close around loops. For small loop curvature,

\[
H_A^{(2)}=\frac{U_A}{2}\sum_\ell E_\ell^2
+\frac{K_A}{2}\sum_pB_p^2+\text{constant}.
\]

The canonical equations are

\[
\hbar\dot a_\ell=U_AE_\ell,
\qquad
\hbar\dot E_\ell=-K_A(\operatorname{curl}^{\dagger}B)_\ell.
\]

## 7. Toy microscopic ring exchange

A simple square-plaquette toy model illustrates how the loop term can arise. Let an isolated link swap have amplitude \(t\), and let each pair of endpoint violations of the local constraint cost energy \(\Delta\). Four oriented swaps around a square return the state to the physical sector.

Under the simplifying assumption that adjacent-link intermediate states cost \(\Delta\) and opposite-link intermediate states cost \(2\Delta\), fourth-order perturbation theory gives

\[
J_\square=t^4\left(\frac{16}{\Delta^3}+\frac{8}{2\Delta^3}\right)
=\frac{20t^4}{\Delta^3}.
\]

The effective term is

\[
-J_\square(W_p+W_p^\dagger)=-2J_\square\cos B_p.
\]

Matching the quadratic convention above gives

\[
\boxed{K_A=\frac{40t^4}{\Delta^3}}
\]

for this toy spectrum. The numerical coefficient is not universal: it changes with the microscopic graph, matrix elements, and intermediate-state energies. The robust point is that local virtual swaps can generate a closed-loop phase stiffness.

## 8. Photon spectrum and two polarizations

On a cubic prototype lattice, the transverse modes obey

\[
\omega_\gamma^2(\mathbf k)
=\frac{4U_AK_A}{\hbar^2}
\sum_{i=1}^{3}\sin^2\left(\frac{k_i\ell}{2}\right).
\]

At long wavelength,

\[
\omega_\gamma^2=c_\gamma^2k^2+O(k^4\ell^2),
\qquad
\boxed{c_\gamma=\frac{\ell\sqrt{U_AK_A}}{\hbar}}.
\]

A spatial link field has three components at nonzero momentum. Local phase redundancy removes one longitudinal coordinate, and Gauss's law removes its conjugate longitudinal momentum. Two transverse oscillators remain. Those are the two photon polarizations.

The regular lattice also predicts higher-order dispersion and directional corrections. A viable fundamental completion must make \(\ell\) sufficiently small or show that a relational/dynamical network flows to an isotropic Lorentz-invariant continuum.

## 9. Quantization and Planck's rule

For each transverse normal mode, define a canonical coordinate \(q\) and momentum \(p=\hbar E\). The quadratic Hamiltonian is a harmonic oscillator,

\[
H_{\mathbf k\lambda}=\frac{p^2}{2M}+\frac12M\omega_{\mathbf k}^2q^2,
\qquad M=\frac{\hbar^2}{U_A}.
\]

Canonical quantization yields

\[
E_n=\hbar\omega_{\mathbf k}\left(n+\frac12\right).
\]

Therefore one added network quantum has

\[
\Delta E=\hbar\omega=hf.
\]

This does not derive quantum mechanics itself; it derives the photon energy ladder from the quantized transverse network modes rather than asserting \(E=\hbar\dot\phi\). Photon number is oscillator occupation number, not scalar winding number.

## 10. Coulomb limit and electromagnetic impedance

Minimize the electric term subject to Gauss's law. The inverse lattice Laplacian approaches \(1/(4\pi r)\) at distances much larger than \(\ell\), giving, in the prototype normalization,

\[
V(r)=\frac{U_A\ell}{4\pi}\frac{Q_1Q_2}{r}
+O\left(\frac{\ell^3}{r^3}\right).
\]

For fundamental physical charge \(e\), matching the continuum Coulomb coefficient gives

\[
\frac{e^2}{\epsilon_0}=U_A\ell.
\]

Combining this with the wave speed gives

\[
\boxed{
\alpha_{\rm bare}=\frac{1}{4\pi}\sqrt{\frac{U_A}{K_A}}.
}
\]

This relation is normalization-dependent and should be understood as the bare coupling of the stated rotor model. Quantum matter and high-energy network modes can renormalize it. The measured value of \(\alpha\) would calibrate the ratio unless the microscopic model independently fixes \(U_A/K_A\).

Define the electromagnetic junction impedance

\[
Z_A=\sqrt{\frac{U_A}{K_A}},
\]

so that in this normalization

\[
Z_A=4\pi\alpha_{\rm bare}.
\]

## 11. Coupling to matter and the QED route

A charged matter defect \(\psi_x\) couples through gauge-covariant hopping,

\[
H_\psi=-t_\psi\sum_{\langle xy\rangle}
\left(\psi_x^\dagger e^{iQa_{xy}}\psi_y+\text{h.c.}\right)+H_{\rm local}.
\]

The transformation

\[
\psi_x\rightarrow e^{iQ\lambda_x}\psi_x,
\qquad
a_{xy}\rightarrow a_{xy}+\lambda_x-\lambda_y
\]

leaves the hopping term invariant. At long wavelength,

\[
a_{x,i}\simeq\frac{e}{\hbar}\int_x^{x+\ell\hat i}A_i\,dx^i,
\]

and the hopping becomes ordinary minimal coupling. If the low-energy charged defects are Dirac fermions, the infrared action can take the form

\[
\mathcal L_{\rm IR}=\bar\psi(i\gamma^\mu D_\mu-m)\psi
-\frac14F_{\mu\nu}F^{\mu\nu}+\text{suppressed corrections}.
\]

This is a route to the QED universality class, not yet a derivation of the observed fermion spectrum, charge assignments, or precision radiative corrections.

## 12. Photon mass and exact constraint protection

Exact local gauge redundancy forbids a single-link pinning term such as \(-h\cos a_\ell\). Therefore the transverse mode is gapless:

\[
\omega_\gamma(\mathbf 0)=0,
\qquad m_\gamma=0.
\]

If the local conservation rule is only approximate, a small pinning term can give

\[
\omega^2=\frac{U_A}{\hbar^2}\left(h+K_A\ell^2k^2\right),
\]

and the corresponding vector gap satisfies

\[
m_\gamma c^2=\sqrt{U_Ah}.
\]

Thus photon mass would measure explicit or emergent violation of the constraint, not the ordinary Josephson link strength.

## 13. The scalar phase mode remains separate

A gauge-neutral local relative phase \(\varphi_x\) may coexist with the link connection. With conjugate imbalance \(n_x\), a candidate Hamiltonian is

\[
H_\varphi=\sum_x\left[\frac{U_\varphi}{2}n_x^2-J_\varphi\cos\varphi_x\right]
+\frac{K_\varphi}{2}\sum_{\langle xy\rangle}(\varphi_x-\varphi_y)^2.
\]

Its small-amplitude spectrum is

\[
\omega_\varphi^2(\mathbf k)=\frac{U_\varphi}{\hbar^2}
\left[J_\varphi+4K_\varphi\sum_i\sin^2\left(\frac{k_i\ell}{2}\right)\right].
\]

This is a massive scalar Josephson/Leggett-like collective mode. Its nonlinear equation may be sine-Gordon-like, and it may support kinks or vortices. It is not the transverse photon.

## 14. Interface with gravity

The electromagnetic connection compares a compact scalar phase. Gravity must compare a complete local frame. A link can therefore carry both

\[
U_\ell=e^{ia_\ell}
\]

and a frame transport element

\[
\Lambda_\ell\in\mathrm{Spin}(1,3),
\]

plus a tetrad-like link displacement \(e_\ell^{\ a}\). Their loop products define different curvatures:

\[
U_f=\prod_{\ell\in\partial f}U_\ell=e^{iF_f},
\qquad
\Lambda_f=\prod_{\ell\in\partial f}\Lambda_\ell=e^{R_f^{ab}J_{ab}}.
\]

The frame-gravity derivation in this folder asks whether the second structure can leave exactly two healthy tensor modes while sharing the same effective causal geometry with the photon.

## 15. Status

### Established within the stated prototype assumptions

- the original one-scalar gradient cannot be the electromagnetic field;
- a constrained compact link network has nonzero loop curvature;
- the quadratic network has two transverse gapless modes;
- canonical quantization gives a photon energy ladder;
- the static kernel approaches Coulomb form;
- minimal matter hopping yields the usual gauge-covariant coupling in the infrared.

### Still open

- deriving the link rotor from a concrete finite junction Hilbert space;
- deriving rather than calibrating \(U_A/K_A\) and \(\alpha\);
- deriving fermionic defects, their charges, and their masses;
- showing exact or sufficiently accurate Lorentz symmetry;
- computing all higher-order operators and precision-QED corrections;
- coupling consistently to the nonlinear frame-gravity sector.

## References used for comparison

- K. G. Wilson, *Confinement of Quarks*, Physical Review D 10, 2445–2459 (1974).
- X.-G. Wen, *Quantum Orders and Symmetric Spin Liquids*, Physical Review B 65, 165113 (2002).
- M. Hermele, M. P. A. Fisher, and L. Balents, *Pyrochlore Photons: The U(1) Spin Liquid in a S=1/2 Three-Dimensional Frustrated Magnet*, Physical Review B 69, 064404 (2004).
- A. J. Leggett, *Number-Phase Fluctuations in Two-Band Superconductors*, Progress of Theoretical Physics 36, 901–930 (1966).
