# The lens source needs an explicit support law

The finite companion source used in the lens pilot cannot be a spherical, collisionless, isotropic particle population in equilibrium. Its density ends abruptly while still positive, which conflicts with the available particle orbits near the edge. This is a conditional failure of that particle interpretation, **not a rejection of stored spacetime deformation or all possible companion reservoirs**.

A nonnegative construction using circular orbits can reproduce the very same source profile. Thus the force law alone does not tell us how the reservoir is supported. We have to derive a capture/transport law that creates an admissible population, or a field equation supplying the required stresses. Declaring that deposits remain where they land does not supply either.

This extends the earlier [support audit](../deposit-support/derivation.md) to the actual empirical force and finite boundaries used for the lens predictions. The [positive stellar-tracer result](../lens-orbit-admissibility/report.md) remains valid: the luminous stars and the hypothetical source are different populations with different density profiles.

## Assumptions and provenance

The particle comparison assumes Newtonian spherical gravity, a static collisionless source, no confining wall, and an effective source density that also acts as ordinary inertial mass density. These are additional assumptions, not consequences of photon conversion. They do not automatically apply to a spacetime field or relativistic radiation. No total photon-supply calculation is made here.

The pilot stipulates

\[
g_b=\frac{GM}{(r+a)^2},\qquad
g_c=Aa_*\left(\frac{g_b}{a_*}\right)^p,\qquad p=0.4624587420
\]

inside a radius r_t, and sets g_c(r)=g_c(r_t)(r_t/r)^2 outside. The baryonic Hernquist field is known; the extra-force power law is an **empirical fit**, not a derived photon-graviton interaction. Its source radii 5, 20 and 100 R_e are sensitivity choices, not measured capture boundaries.

Applying the **known spherical Poisson relation** gives the following conditional source:

\[
\rho_c(r)=\frac{g_c(r)}{4\pi G}
\left[\frac2r-\frac{2p}{r+a}\right]>0\quad(r<r_t),
\qquad \rho_c(r)=0\quad(r>r_t).
\]

The force and enclosed mass are continuous at r_t, so there is no mass sheet there. The volume density nevertheless has a finite jump. Removing infinite outer mass by keeping enclosed mass fixed did not by itself solve dynamical support.

## Why isotropic particles cannot realize this boundary

Define the available potential above the edge as

\[
\psi(r)=\Phi(r_t)-\Phi(r)=\int_r^{r_t}(g_b+g_c)\,ds.
\]

An isotropic nonnegative equilibrium f(E), with no particles outside r_t, must have its energy support above the boundary threshold. In terms of binding energy relative to that threshold, the **known velocity integral** is

\[
\rho(\psi)=4\pi\sqrt2\int_0^\psi f(E)\sqrt{\psi-E}\,dE.
\]

Differentiating yields the necessary inequality

\[
\frac{d\rho}{d\psi}\ge\frac{\rho}{2\psi},\qquad
\mathcal D(r)\equiv
-\frac{2\psi}{g_b+g_c}\frac{d\ln\rho}{dr}\ge1.
\]

Each nonnegative term obeys this bound since 1/(psi-E)>=1/psi. Singular distributions must still give finite density; a nonnegative measure at the threshold also cannot produce the stipulated finite edge jump. Equivalently, rho/sqrt(psi) must be nondecreasing with psi. At the pilot's edge, psi tends to zero but the density tends to a positive finite value with finite inward logarithmic slope, so D tends to zero. This proves a failure near the edge for every positive amplitude in this family, independently of a numerical grid.

This is an application of established distribution-function mathematics, with no novelty claim. Related finite-radius consistency problems and viable tangential alternatives are studied by [Baes (2022)](https://academic.oup.com/mnras/article/512/2/2266/6548135). The standard isotropic inversion is also documented in [galpy](https://docs.galpy.org/en/v1.11.0/reference/dfeddington.html).

In plain language: at a free outer edge, an isotropic population includes outward-moving particles. To prevent any from going beyond that edge, its allowed speeds and population must thin out appropriately. Our prescribed density does not thin out. A positive velocity dispersion from the averaged Jeans equation would not fix this problem.

## Size of the inconsistency in the archived models

The script reads the 33 **training** lens masses already inferred from aperture dispersions, updated I-band sizes and conditional static geometry. It evaluates the three cutoffs for the archived isotropic stellar case with 1.5-arcsecond seeing, for 99 source configurations. Stellar isotropy does not imply companion isotropy; the latter is the optional interpretation being tested here.

| Source boundary | Median radius where D falls below 1, as fraction of boundary | Median equivalent source mass in the outer failing region |
|---|---:|---:|
| 5 R_e | 0.758 | 27.9% |
| 20 R_e | 0.767 | 25.4% |
| 100 R_e | 0.772 | 24.4% |

Crossings are numerically located after inspecting 2,000 radial probes per configuration; the global impossibility follows from the analytic edge argument. These percentages describe where this necessary condition fails in the specified profiles. They are **not a required percentage of mass to remove**, a full distribution-function inversion, or observed missing mass. Passing this inequality at smaller radius is not sufficient to establish admissibility there.

No parameters are fitted, and no validation or test predictions are newly evaluated. Existing lens residuals and their unfavorable isotropic validation result are unchanged.

## A positive alternative that preserves the force

For each spherical shell, populate circular orbits with speed

\[
v_{\rm circ}^2(r)=r[g_b(r)+g_c(r)].
\]

Distribute orbital planes and phases uniformly, with equal senses of circulation, and assign shell weight dM=4pi r^2 rho_c dr. Every weight is nonnegative. Each orbit remains on its shell in the fixed spherical potential, producing the desired density with no radial flux even at the cutoff. This is a **known circular-orbit construction**, not a new deposition mechanism; it is singular in velocity space rather than a smooth isotropic distribution.

Its stresses satisfy the known radial momentum equation:

\[
P_r=0,\quad P_\theta=P_\phi=\frac{\rho_c v_{\rm circ}^2}{2},\qquad
\frac{dP_r}{dr}+\frac{2P_r-P_\theta-P_\phi}{r}=-\rho_c g.
\]

For the representative 20 R_e cutoff, the required circular speed at the source edge ranges from 204 to 368 km/s across the 33 archived models, with median 274 km/s. These are conditional model speeds, not detected companion velocities. They illustrate how substantial directed orbital motion must replace the picture of stationary particles simply staying at their landing sites.

The radial single-orbit restoring coefficient is kappa^2=g'+3g/r>0 on either side of the boundary in this force family. That addresses small radial disturbances of an individual circular orbit in the fixed potential. It does **not** prove collective stability when the source and gravity respond together.

Six direct orbit integrations, spanning all three boundaries and two radii per boundary, preserve radius to better than 8e-12 over ten periods. Independent integration of source density recovers the prescribed enclosed mass within 5e-16 relative error; potential integration agrees within 2e-15, and finite-difference density slopes within 7e-11. These checks verify the implementation and construction, not formation, long-term collective stability or observational agreement.

## Consequence for the next physical model

Three distinct branches remain possible, but none is silently adopted:

1. An isotropic particle reservoir needs a density that tapers consistently with a nonnegative energy distribution. Solving that distribution and its gravitational field together will generally change the lens predictions near the boundary; an arbitrary force cutoff is insufficient.
2. An orbit-supported reservoir can keep a sharper boundary, but capture and transport must supply its angular-momentum distribution. Random arrival directions alone do not demonstrate circularization. Momentum, energy and any dissipation must balance.
3. A deposited spacetime field needs its own dynamical stress and metric equations. The particle density and pressure formulas above cannot simply be assigned to that field without deriving the correspondence.

The physically decisive next task is a shared source/support law followed by renewed motion and lensing predictions. A new per-galaxy cutoff or orbital choice selected to improve each residual would hide this missing connection. No strong case for the full photon-companion theory is claimed from this consistency audit.

Reproduce with `python research_work/results/deposit-boundary-admissibility/run.py`. Input hashes, all 99 configurations and numerical checks are archived beside this report. Dependencies: Python, NumPy and SciPy.
