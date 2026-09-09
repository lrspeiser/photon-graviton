# A common matter–light–gravity action: derivation and the nonexpansion obstruction

8 September 2026. This is a restricted-model derivation, not a completed successful theory of temporal redshift. It extends the earlier lapse-only calculation to universal conformal matter coupling. The purpose is to determine whether clock protection, common wave propagation and observable redshift can coexist with physical nonexpansion in this explicit family.

## 1. An action that implements common local matter physics and wave characteristics

Use units c0 = hbar = 1 in the action. Introduce a gravitational metric q, a canonical scalar phi, and a positive conformal function B(phi). Define the physical matter metric by

\[
g^{(m)}_{\mu\nu}=B(\phi)^2q_{\mu\nu}.
\]

Choose

\[
S=\int d^4x\sqrt{-q}\left[\frac{M^2}{2}R[q]
-\frac12q^{\mu\nu}\partial_\mu\phi\partial_\nu\phi-V(\phi)\right]
+S_{\rm SM}[g^{(m)},\Psi,A_\mu]. \tag{1}
\]

Here the Standard Model action includes the electromagnetic field, the charged particles that emit and detect it, and the interactions that form matter. Its parameters are constant when measured with the matter metric. There is no separate photon field reserved for traveling light. There is no screening and no imposed rotation. B and V remain functions to be specified; the obstruction below holds for any smooth positive B and any V within the stated homogeneous assumptions.

This is a familiar scalar–tensor framework, used here as a controlled candidate rather than claimed as a novel action. Conformal/disformal frame transformations and the need to track physical matter couplings are discussed by Bettoni and Liberati [1].

In a local freely falling material frame, the matter action has its usual form. It introduces no explicit varying electromagnetic or nuclear parameter that would differentially shift optical and hyperfine clocks. A material cavity with fixed proper dimensions has its usual local resonant frequency. These are properties of the stated local matter action; environmental effects, finite-size tidal corrections, gravitational binding and a quantum-gravity completion are not being claimed to vanish. Merely showing this local protection does not establish all experimental bounds on an additional scalar force.

The Maxwell action is conformally invariant in four spacetime dimensions at the classical level: its B factors cancel between the volume factor and the two inverse metrics. Its geometric-optics phase obeys

\[
(g^{(m)})^{\mu\nu}k_\mu k_\nu=0
\quad\Longleftrightarrow\quad q^{\mu\nu}k_\mu k_\nu=0. \tag{2}
\]

For Einstein gravity with a canonical scalar, the tensor modes have the q-metric null characteristics as well. Thus this action derives shared photon/tensor characteristics rather than separately setting their speeds at one epoch. It does not imply identical amplitudes or luminosity distances for electromagnetic and gravitational radiation. Scalar modes, source generation and their phenomenology remain additional questions.

## 2. The clock calculation determines what counts as expansion

Take a homogeneous nonrotating geometry

\[
ds_q^2=-c_0^2N(t)^2dt^2+a(t)^2d\Sigma_k^2.
\]

Matter proper time and material spatial scale are

\[
d\tau=B N\,dt,\qquad A(t)=B(t)a(t). \tag{3}
\]

The physical metric becomes

\[
ds_m^2=-c_0^2d\tau^2+A(\tau)^2d\Sigma_k^2.
\]

The operational definition of a nonexpanding universe is that separations of comoving sources remain constant in units of unchanged local material rulers. In this action that is A = constant. Keeping a constant while changing B does not satisfy that definition: the matter rulers and clocks couple to g_m, not to q alone.

Let Q denote the conserved comoving spatial momentum magnitude of a geometric-optics ray in this homogeneous geometry. Its coordinate-time frequency is proportional to NQ/a. A material clock converts coordinate time to proper time through BN, so the received frequency in material units is

\[
\nu_{\rm material}\propto\frac{NQ/a}{BN}=\frac{Q}{A}. \tag{4}
\]

Consequently, emission in an unchanged local atomic transition gives

\[
\boxed{1+z_{\rm observed}=\frac{A_o}{A_e}}. \tag{5}
\]

Both the lapse and any arbitrary naming of the time coordinate have canceled. Tensor phases transported in the geometric-optics limit obey the same local frequency scaling. Ordinary peculiar velocities and gravitational redshift from inhomogeneities must be added separately; they are not a universal distance-only temporal shift in the homogeneous background.

Therefore

\[
\boxed{\dot A=0\ \Longrightarrow\ z_{\rm cosmological}=0}. \tag{6}
\]

This is the central obstruction. It is valid for the action family (1) and the homogeneous transparent propagation assumptions. It is not an impossibility theorem for all nonminimal electromagnetic interactions, nonmetric propagation, anisotropic geometries, frequency-changing collisions or other theories.

## 3. The original affine timing field makes the tradeoff explicit

Set the auxiliary spatial scale to unity, use a static spatial metric h_ij, and take N = 1/n(t). Then

\[
ds_q^2=-c_0^2dt^2/n(t)^2+h_{ij}dx^idx^j.
\]

With B = 1, the material metric is q. The change dTau = dt/n removes the homogeneous lapse. Wave frequency and atomic frequency per coordinate t both carry the same factor 1/n, and there is no measured cosmological redshift.

With B = n, the material metric instead becomes

\[
ds_m^2=-c_0^2dt^2+n(t)^2h_{ij}dx^idx^j. \tag{7}
\]

Now material proper time is t, and the observable redshift is n_o/n_e. Atomic ratios and local cavity ratios remain protected, but physical separations grow as n. A coordinate wave speed c0/n is accompanied by material rulers whose coordinate lengths scale as 1/n; the local measured light speed remains c0. This is expansion relative to material rulers, even though h_ij is static.

Normalize n_o = 1 and impose the earlier affine history n(t) = 1 + gamma(t-t_o), gamma = c0 kappa. The present-normalized radial geometric distance is

\[
R=\int_{t_e}^{t_o}\frac{c_0\,dt}{n(t)}
=\frac{1}{\kappa}\ln\frac{1}{n_e}.
\]

Thus

\[
1+z=\frac1{n_e}=e^{\kappa R}. \tag{8}
\]

The desired exponential survives, but it has become a linearly expanding material geometry. The coefficient remains an imposed history parameter, not a prediction of a chosen V(phi). The choice B=n and the affine history in this section are kinematic examples; no claim is made that the chosen static q and matter content solve (1). The conclusion that they violate operational nonexpansion already follows before solving the field equations.

## 4. Varying the action does not remove the obstruction

Define alpha = d ln B/d phi and let T_m^(q) be the trace of the matter stress tensor obtained by varying with respect to q. The equations from (1) are

\[
M^2G_{\mu\nu}[q]=T^{(\phi)}_{\mu\nu}+T^{(m,q)}_{\mu\nu},
\qquad \Box_q\phi=V_{,\phi}-\alpha T_m^{(q)}. \tag{9}
\]

In q proper time t_q, with H_q = (1/a) da/dt_q and c0=1,

\[
3M^2(H_q^2+k/a^2)=\rho_q+\tfrac12\dot\phi^2+V,
\]
\[
-2M^2(\dot H_q-k/a^2)=\rho_q+p_q+\dot\phi^2. \tag{10}
\]

Dots in (10) denote q proper-time derivatives. The material expansion rate is

\[
H_m=\frac{H_q+\alpha\dot\phi}{B}. \tag{11}
\]

A physically static solution could in principle have H_q = −alpha dot(phi), depending on B, V and matter. But imposing H_m = 0 still makes the measured frequency (4) constant along the homogeneous expansion history. Changing the background potential cannot make (5) yield redshift while retaining A constant. Full stability cannot rescue a branch that already lacks the target observable.

For the separate Einstein-gravity subclass B = constant with a static physical metric, (10) also gives rho+p = 2M² k/a² in natural units. Negative spatial curvature requires negative total enthalpy. Ordinary matter, radiation and a canonical scalar cannot supply it, regardless of scalar potential. This sign argument must not be indiscriminately applied to the variable-B physical-frame equations as if they were unmodified Einstein equations. The more general redshift obstruction (6), however, already covers the entire conformal family considered here.

## 5. Luminosity and angular distances from the same construction

Let chi be the dimensionless radial coordinate of a source and S_k(chi) equal sin(chi), chi or sinh(chi), for normalized positive, zero or negative curvature. In the material metric, the physical transverse size corresponding to a small observed angle is A_e S_k(chi) times that angle. Hence

\[
D_A=A_e S_k(\chi).
\]

With conserved photons and transparent propagation, each photon loses an energy factor 1+z and the arrival interval grows by the same factor. The receiving sphere has area 4 pi A_o² S_k(chi)², yielding

\[
F=\frac{L}{4\pi A_o^2 S_k(\chi)^2(1+z)^2},
\qquad D_L=(1+z)A_o S_k(\chi).
\]

Combining with (5),

\[
\boxed{D_L=(1+z)^2D_A}. \tag{12}
\]

This is derived here from the material metric and conservation assumptions rather than imported as a universal law for every candidate. It agrees with the usual metric distance reciprocity relation; its assumptions and modifications are discussed in [2]. It applies to this classical transparent electromagnetic construction, not automatically to gravitational-wave amplitudes or the earlier fixed-material nonmetric transport branch.

In the earlier fixed-material transport branch, the corresponding derivation instead used a fixed source-size geometry and gave

\[
D_A^{\rm static}=D_M,\qquad D_L^{\rm static}=(1+z)D_M.
\]

For the material coasting geometry above, using the same present-normalized path and curvature,

\[
D_A^{\rm metric}=\frac{D_M}{1+z},\qquad D_L^{\rm metric}=(1+z)D_M. \tag{13}
\]

Both use D_M = Rc sinh[ln(1+z)/(kappa Rc)] for the negative-curvature comparison. Their brightness laws are identical but their angular-distance laws differ by 1+z. In the flat case the common luminosity law is (1+z)ln(1+z)/kappa, also found in existing coasting cosmology work [3]. This equivalence does not imply that either model explains all other observations.

Using the previous exploratory fit Rc = 4.547 Gpc and fixed kappa, without a new fit:

| z | Common luminosity distance, Gpc | Fixed-material transport angular distance, Gpc | Metric coasting angular distance, Gpc |
|---|---:|---:|---:|
| 0.1 | 0.416 | 0.378 | 0.344 |
| 0.5 | 2.462 | 1.642 | 1.094 |
| 1.0 | 5.839 | 2.919 | 1.460 |
| 2.0 | 15.164 | 5.055 | 1.685 |

These are ideal homogeneous-frame predictions. The earlier Pantheon+ calculation's zHEL/zHD conventions are retained in its fit; they are not independently rederived here. A new action must also justify source calibration. Equal luminosity-distance curves would have equal supernova fit statistics only under the same adopted luminosity and calibration assumptions. The table is not a new angular-distance measurement or likelihood.

## 6. What has actually been derived

The action yields protected ordinary local matter physics and shared photon/tensor characteristics, but physical nonexpansion forces zero homogeneous redshift. A conformal choice that restores the exponential gives operational expansion instead. The clock and redshift calculation establishes this for a broader class than the earlier lapse-only example; it does not provide the requested successful nonexpanding theory.

The paired distance derivation also shows why the supernova brightness fit cannot settle this question alone. Different physical interpretations can share the same brightness law while predicting angular distances different by a factor of two at z=1. Testing angular sizes would require an independently justified physical ruler or source model; standard CMB/BAO rulers cannot simply be imported into the fictional nonexpanding model.

To keep the stipulated nonexpansion, a next candidate must depart from the universal conformal metric construction—for example through a specified nonminimal interaction that changes radiation frequency relative to matter. It must then demonstrate atomic protection, common tensor behavior, energy conservation and stable background dynamics rather than assign those properties. No such successful interaction has been derived in this note. This result narrows the mathematical search and prevents a relabeling of expansion from being mistaken for a new mechanism.

## Reproduction and sources

`check_derivation.py` checks the redshift cancellation, its recovery with expanding material scales, and both distance identities at four redshifts. These are numerical checks of analytic formulas, not new data validation. Its only empirical inputs are the previously fitted curvature radius and fixed redshift coefficient; both are explicitly recorded in the code and output.

1. Bettoni and Liberati, Disformal invariance of second order tensor-scalar theories: framing the Horndeski action (2013), https://arxiv.org/abs/1306.6724 .
2. More et al., Modifications to the Etherington Distance Duality Relation and Observational Limits (2016), https://arxiv.org/abs/1612.08784 .
3. Sultana, The Rh = ct universe and quintessence (2016), https://arxiv.org/abs/1601.05474 .
