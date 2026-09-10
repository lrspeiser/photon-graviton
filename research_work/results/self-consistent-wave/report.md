# Nine supported wave-source profiles, with motion and lensing from one density

We have constructed nine numerical stationary wave configurations in a fixed spherical ordinary-matter galaxy, including the wave's own gravity. They satisfy the coupled field equations and independent equilibrium checks. Unlike the empirical source previously imposed, these profiles are obtained from a support equation; they have regular central density and rapidly decaying outer tails.

This is a constructive candidate calculation, **not observed galaxy agreement, a photon-capture solution or a stability proof**. It does not turn the earlier unfavorable lens validation into a success. The deposited mass and wave parameter are still specified inputs with no derived connection to redshift energy transfer.

## Equations and status

The free, nonrelativistic stationary Schrodinger-Poisson equations are **known physics**, not a new photon-graviton law:

\[
-\frac{\hbar^2}{2m}\nabla^2\psi+m(\Phi_b+\Phi_c)\psi=\mu\psi,
\qquad \nabla^2\Phi_c=4\pi Gm|\psi|^2.
\]

The fixed ordinary-matter potential is the known Hernquist form Phi_b=-GM_b/(r+a). We use this spherical ordinary galaxy as a controlled starting point, not a complete model of a rotating Galactic bar. The weak-field, slowly varying envelope approximation must hold in a physical application. The field is not an ordinary massless graviton: calling it a companion remains a proposed identity requiring a conversion/capture interaction.

The stationary-equilibrium approach and wave/fluid equations are established methods, as in [Construction and Evolution of Equilibrium Configurations of the Schrodinger-Poisson System](https://doi.org/10.3390/universe8080432). No expansion or cosmological dark-matter formation assumption is used in this calculation.

Define x=r/a, f=M_c/M_b, the dimensionless self-potential v=Phi_c/(GM_b/a), and

\[
\eta=\frac{\hbar^2}{m^2GM_ba},\qquad
\rho_c=\frac{M_b}{4\pi a^3}\frac{u(x)^2}{x^2},\qquad
\mathcal M(x)=M_c(<r)/M_b.
\]

Then the boundary-value problem is

\[
u''=\frac2\eta\left[-\frac1{1+x}+v-\epsilon\right]u,
\qquad v'=\frac{\mathcal M}{x^2},\qquad \mathcal M'=u^2.
\]

The unknown eigenvalue epsilon=mu/(mGM_b/a) is solved along with the profile. Regularity requires u proportional to x at the center, so u/x stays finite. At infinity u tends to zero, the enclosed mass tends to f, and v tends to -f/x. These are **a nondimensionalization and numerical application of known equations**, not a novel fundamental derivation.

## Numerical method and independent checks

A collocation boundary-value solver follows a nodeless branch by increasing f in small steps. Overall sign of u is physically irrelevant; a genuine sign change away from numerical tail noise is rejected. Nodelessness by itself is not a proof that this is the global energy minimum or collectively stable.

The initial domain is 0.001<=x<=60 with tolerance 2e-7. Each reported case is recomputed on 0.0005<=x<=90 with tolerance 2e-9. At the inner boundary the regular expansions u=epsilon_radius*u' and M=epsilon_radius*u^2/3 are used. At the numerical outer boundary u=0, M=f and v=-f/x. This finite numerical boundary is an approximation to the decaying tail, not a physical wall adopted in the theory. The very small outer mass and wider-domain agreement support this approximation for these cases.

The solver had two unsuccessful early approaches: placing the inner boundary at 1e-6 caused mesh exhaustion near the singular coordinate origin at tight tolerance, and large normalization jumps could converge to a nodal state. Those failures were not accepted as solutions. The regular inner expansion, explicit inner-radius check and smaller continuation steps resolve the numerical issues in the reported family. The physical kinetic-energy integral uses (u'-u/x)^2 on the finite domain; substituting u'^2 without its boundary term gave an erroneous virial residual and was corrected.

The independent virial relation in these units is

\[
2T+W_c-\int u^2\frac{x}{(1+x)^2}\,dx=0,
\quad T=\frac\eta2\int(u'-u/x)^2dx,
\quad W_c=\frac12\int u^2v\,dx.
\]

The external-potential term is not the total ordinary-matter potential energy: it is the appropriate radial force moment. Treating the Hernquist galaxy as fixed excludes its dynamical response and any combined-system stability claim.

Across all nine refined cases:

- Maximum mass-normalization error is 2.8e-9 relative to specified source mass.
- Maximum normalized virial residual is 4.0e-9.
- Maximum change in enclosed mass between numerical setups is 1.45e-7 of total source mass, over 0.05<=r/a<=30.
- Maximum relative change in half-mass radius is 1.59e-7.

These checks establish numerical consistency for this boundary-value problem. They do not test perturbation growth, non-spherical modes, formation or observables.

## Resulting sizes

| Wave parameter eta | Source/ordinary mass f | Radius containing half the source mass, in units a |
|---:|---:|---:|
| 0.3 | 0.1 | 1.3844 |
| 0.3 | 1 | 0.8146 |
| 0.3 | 3 | 0.3720 |
| 1 | 0.1 | 2.8582 |
| 1 | 1 | 1.9642 |
| 1 | 3 | 1.0891 |
| 3 | 0.1 | 6.1145 |
| 3 | 1 | 4.4940 |
| 3 | 3 | 2.7644 |

For the Hernquist ordinary light profile R_e=1.8153a. At fixed eta, adding source mass strengthens self-gravity and makes the source more compact in this family. Increasing eta spreads the wave out. These are derived parameter sensitivities, not measured source radii. The parameter choices sample the equation; they were not selected to fit any galaxy.

In particular, eta cannot be held fixed across arbitrary galaxies if the particle mass m is universal: eta varies as 1/(M_b a). The source fraction f is also not independently predictable until the capture and accumulation mechanism is specified. Fitting both freely per galaxy would not establish photon origin.

## One source for stellar motion and light bending

The conditional circular speed is

\[
\frac{v_{\rm circ}^2}{GM_b/a}=\frac{x}{(1+x)^2}+\frac{\mathcal M(x)}x.
\]

For the leading weak-field, nonrelativistic gravitational source, the extra metric potentials are equal to the order retained. Gradient stresses and relativistic corrections must be small for this approximation; they are not declared zero in the full scalar theory. Under that approximation the **known weak-lensing integral** gives

\[
\frac{\widehat\alpha_c(b)}{4GM_b/(ac^2)}
=\int_0^{\pi/2}\frac{\mathcal M(\beta/\cos\theta)\cos\theta}{\beta}\,d\theta,
\quad \beta=b/a.
\]

The code evaluates four impact parameters per case. An independent integral of density inside the projected cylinder gives the same deflection, with maximum relative difference 5.2e-11. No separate lensing multiplier is fitted. The bending angle here is an asymptotic deflection, not an observed Einstein angle; source/lens geometry, luminous stellar orbits, aperture effects and data uncertainties are still needed for that comparison.

For a physical M_b,a mapping, one must check GM_b/(ac^2) times the dimensionless potential depths is small, and that the envelope energies are small relative to mc^2. Dimensionless equilibrium alone does not establish those conditions at arbitrary physical scales.

## What this changes in the research plan

We now have a working forward solver for a specific supported source, rather than only an empirical force or an inverse inconsistency check. This makes it possible to test an explicit wave branch against motions and lensing. The next observational analysis must use a shared physical mass, proper ordinary-matter uncertainties and the derived profile, while keeping source-normalization assumptions visible. It must not silently inherit the old aperture-derived masses as correct for the new potential.

Before promoting any success as support for the full theory, we still need photon conversion and event timing, capture and energy/momentum balance, the origin of f, lifetime and stability, and broad observational/held-out checks. The total photon-supply question remains deferred, not passed. No validation or test scores were opened in this pass.

Reproduce with `python research_work/results/self-consistent-wave/run.py`. Requires NumPy and SciPy. `results.json` contains equilibrium and lensing checks; `profiles.json` contains 150 sampled radii for each case. No observational catalog or fitted redshift/gravity coefficients are inputs to this numerical family.
