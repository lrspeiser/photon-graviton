# Stored deformation must supply the right stress, not just energy density

The simplest literal frozen-field interpretation does not supply the gravity assumed in the galaxy/lens pilot. A static, minimally coupled canonical scalar with nonnegative potential energy has tensions that prevent it from acting as the positive, nearly pressureless extra source used there. This is a **restricted candidate failure**, not a prohibition on stored spacetime deformation, modified gravity or oscillating fields.

The constructive consequence is to distinguish a time-independent field from a time-independent average density. A rapidly oscillating field can have approximately pressureless average stress. Such a branch requires a field mass, a supported spatial solution and a capture interaction; the empirical gravity profile does not determine them.

## Reverse the pilot metric into its required source

Use the known static weak-field metric in isotropic coordinates,

\[
ds^2=-(1+2\Phi/c^2)c^2dt^2+(1-2\Psi/c^2)d\mathbf x^2.
\]

Slow stars respond to -grad Phi, while weak light bending depends on the transverse gradient of Phi+Psi. In ordinary Einstein gravity, to first order in the metric perturbation, a spherical source with local energy density u and radial/tangential stresses P_r,P_t obeys

\[
\nabla^2\Psi=\frac{4\pi G}{c^2}u,
\quad P_r=\frac{c^2}{4\pi G r}(\Phi-\Psi)',
\quad P_t=\frac{c^2}{8\pi G}\left[(\Phi-\Psi)''+\frac{(\Phi-\Psi)'}r\right].
\]

Primes denote radius derivatives. Consequently

\[
\nabla^2\Phi=\frac{4\pi G}{c^2}(u+P_r+2P_t).
\]

These are **known linearized Einstein equations**, not a new gravity law. They include stresses of the same order as energy density; stresses smaller by the weak-field parameter only affect higher orders. All sources must also obey their equations of motion and stress conservation. Prescribing arbitrary functions in these equations does not guarantee a static solution.

The pilot assumes equal extra metric potentials, Phi_c=Psi_c, and Phi_c'=g_c. It therefore requires, at this order,

\[
u_c=\frac{c^2}{4\pi G}\left(g_c'+\frac{2g_c}{r}\right)>0,
\qquad P_{r,c}=P_{t,c}=0.
\]

This is a **conditional reverse-engineering result** for the pilot metric. It does not prove that deposits are matter or that their energy budget is available. For orbit-supported matter, small supporting stresses of relative order v^2/c^2 can coexist with this leading approximation. An order-unity tension cannot be ignored in the same way.

## Test one explicit frozen-field candidate

Take a real, canonical, minimally coupled scalar with action, in c=hbar=1 units,

\[
S_\phi=\int\sqrt{-g}\left[-\frac12g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi-V(\phi)\right]d^4x.
\]

This is a **known comparison action**, not an identification of the field with gravitons. It includes no fifth-force coupling, modified Einstein term, external support or photon-capture source. For a static radial field, let K=(dphi/dell)^2/2, where ell is proper radial distance. Its known local stresses are

\[
u=K+V,\qquad P_r=K-V,\qquad P_t=-K-V=-u.
\]

A nonzero energy density therefore has tangential tension of equal magnitude. It cannot supply the pilot's negligible leading stresses. Moreover,

\[
u+P_r+2P_t=-2V,\qquad
\nabla^2\Phi_\phi=-\frac{8\pi G}{c^2}V.
\]

For V>=0 this local source of the temporal potential is nonpositive. The fitted extra field instead has

\[
\nabla^2\Phi_c=g_c\left[\frac2r-\frac{2p}{r+a}\right]>0
\]

inside its source boundary, because p=0.4624587420<1. A static canonical field with nonnegative V cannot generate that profile under these assumptions, even if equal metric potentials are relaxed. This is a statement about the local Laplacian, not a claim that every externally supported scalar configuration repels every object. Boundary sources, interactions and other stress contributions change the full problem and must be included.

The stronger equal-potential inconsistency does not depend on the sign of V: P_t=-u is incompatible with nonzero u and vanishing leading P_t. Choosing a negative potential to reverse one sign would not repair the complete required stress tensor.

## A separate equilibrium obstruction

For a static metric with lapse N>0 and spatial metric h_ij, the source-free scalar equation is

\[
\partial_i(N\sqrt h\,h^{ij}\partial_j\phi)=N\sqrt h\,V_{,\phi}.
\]

Multiply by phi and integrate. For a regular localized field with vanishing boundary term,

\[
\int N\sqrt h\,[h^{ij}\partial_i\phi\partial_j\phi+\phi V_{,\phi}]\,d^3x=0.
\]

If phi V_,phi>=0, as for a nonnegative mass term about its vacuum, both terms are nonnegative and the only localized solution is the vacuum. A static ordinary gravitational well does not remove this identity. Sources, horizons/boundaries with nonzero flux, noncanonical terms or time dependence require different arguments. This is established no-static-scalar reasoning; [Hod's static boson-star theorem](https://arxiv.org/abs/1902.05230) gives a related self-gravitating result. No new no-go theorem is claimed here.

The two restrictions have different assumptions: the local force sign uses V>=0; the integrated equilibrium identity uses phi V_,phi>=0 and the stated boundary conditions. Neither should be applied beyond its own scope.

## What time dependence can change

For a time-dependent radial canonical field, define Q=(dphi/dtau)^2/2 and K=(dphi/dell)^2/2 in natural units. Then

\[
u=Q+K+V,\quad P_r=Q+K-V,\quad P_t=Q-K-V.
\]

For the known local harmonic example phi=F cos(mt), V=m^2 phi^2/2, with constant F and negligible spatial gradients,

\[
\langle Q\rangle=\langle V\rangle=\frac{m^2F^2}{4},\quad
\langle u\rangle=\frac{m^2F^2}{2},\quad
\langle P_r\rangle=\langle P_t\rangle=0.
\]

This illustrates a possible stress structure, **not a finite galaxy solution**. With F varying over the Galaxy, gradients and gravity affect both the oscillation and envelope equation; the local oscillator cannot be pasted into an arbitrary fitted density. An averaging approximation needs periods short compared with the gravitational evolution time and a self-consistent spatial solution. Oscillating real fields and complex fields with stationary stress have a substantial existing literature, including [the boson-star review](https://doi.org/10.1007/s41114-023-00043-4). This is a known avenue, not our invention.

Calling the field a stored deformation is insufficient to determine its mass, pressure, stability or photon coupling. Minimal coupling is only a comparison choice. If we modify gravity instead, we must derive the new relation between the field, both metric potentials and detector/matter motion; the Einstein-source formulas above would no longer automatically apply.

## Reproducible checks and observational status

The companion script checks the positive required source on 1,188 interior probes across the 99 archived training configurations (33 lenses, three cutoffs). It verifies the radial Laplacian independently by differentiating r^2 g_c, and integrates the local harmonic example over a period to verify its averaged stress. These are implementation and candidate-consistency checks, not independent data points or successful galaxy predictions. The analytic sign argument covers all interior radii of the stated profile.

The earlier unfavorable isotropic lens validation remains unchanged. No additional holdout scores, mass adjustments, independent lensing multiplier or assumed photon-supply budget are introduced.

The next viable field calculation must choose a specified dynamical action, solve its supported spatial configuration, derive its motion and lensing together, and then connect conversion/capture to that configuration. A near-pressureless oscillating source is a candidate to investigate, not a completed explanation. The full redshift, event timing, broadband observations and held-out galaxy objectives remain open.
