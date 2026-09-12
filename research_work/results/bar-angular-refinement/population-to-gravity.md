# From moving captured energy to a gravitational prediction

## Why central visits are not the final target

A trajectory that just grazes a diagnostic sphere can switch from “never entered” to “entered” with an arbitrarily small change of direction. Its time inside tends to zero at grazing, but its ever-entry flag jumps to one. For a continuous source, its ever-entry history contribution can also jump by a finite amount if grazing occurs well before the observation time. Therefore central encounter fractions can converge more slowly with direction sampling than actual residence fractions. These are mathematical properties of these diagnostics, not new physics or proof that a particular quadrature is adequate.

Gravity depends on where the deposited material is **now**. Counting all material that has ever passed a location as if it stayed there silently adds a retention rule. We have not adopted that rule for these moving-particle tests.

## Conditional forward model

**Project assumptions:** use the declared cold captured-particle interpretation, a stationary source in rotating bar coordinates, a fixed ordinary potential, and an initially negligible captured population. Let dot M_R be the rest-mass injection rate on launch sphere R. The current two spheres have separate, unspecified rates; do not choose their mixture by fitting a desired halo. Let w_R(Omega) be the normalized source probability density on launch directions. X_R(a,Omega) is the trajectory at age a in the same rotating coordinates.

**Known transport mathematics, conditionally applied:** a zero-initial-population source switched on at time zero gives

\[
\rho_c(\boldsymbol{x},T)=\sum_R\dot M_R
 \int d\Omega\,w_R(\Omega)\int_0^T da\,
 \delta^3\!\left[\boldsymbol{x}-\boldsymbol{X}_R(a,\Omega)\right].
\]

The normalization is total injected mass, sum_R dot M_R T, when all positions including escaped material are included. No central sink or permanent in-place deposit is implicit. With time-dependent injection, dot M_R becomes dot M_R(T-a); with a changing potential, the trajectory depends on birth time as well as age. Rest-mass formation must still be connected to the physical energy/momentum budget.

Integrating this expression over a sphere gives the residence statistic already measured:

\[
M_c(<r,T)=\sum_R\dot M_R T f_{{\rm inside},R}(r,T).
\]

For example, a 0.5 percent residence fraction means 0.5 percent of that source's **total injected rest mass** is inside at that time. It is not an absolute mass until the injection rate is specified and funded.

**Known Newtonian potential and force, conditional on ordinary gravitation of captured rest mass:**

\[
\Phi_c(\boldsymbol{x},T)=-G\sum_R\dot M_R
\int d\Omega\,w_R\int_0^T\frac{da}{|\boldsymbol{x}-\boldsymbol{X}_R|},
\]
\[
\boldsymbol{a}_c(\boldsymbol{x},T)=-G\sum_R\dot M_R
\int d\Omega\,w_R\int_0^T da\,
\frac{\boldsymbol{x}-\boldsymbol{X}_R}{|\boldsymbol{x}-\boldsymbol{X}_R|^3}.
\]

These are the standard superposition integrals, not novel force laws. They implement the user's “deeper well where energy collects” interpretation: the force follows from the spatial gradient of the resulting potential. Incoming wave direction does not itself specify the gravitational force direction. If the retained state instead has non-particle stress, these rest-mass formulas require an explicit replacement.

The integrals require appropriate spatial resolution near a source trajectory. A finite trajectory quadrature is not a set of literal point masses. Introducing smoothing must be declared and convergence checked rather than interpreted as a physical core by default.

**Known exact angular-average identity:** even for nonspherical density, the monopole potential is

\[
\overline\Phi_c(r,T)=-G\sum_R\dot M_R
\int d\Omega\,w_R\int_0^T\frac{da}{\max(r,|\boldsymbol{X}_R|)},
\qquad
\langle a_{c,r}\rangle=-\frac{GM_c(<r,T)}{r^2}.
\]

This is the shell-theorem kernel for the spherical average. It is useful for a first mass-to-force diagnostic, but cannot predict bar-direction or above/below-disk differences for individual stars. Those require the full spatial integral and the ordinary matter components.

## What the next solver must retain

1. Trajectories or age-integrated spatial kernels, rather than only first-entry times.
2. Source rates tied to transported energy, rather than an arbitrary normalization at each radius.
3. Full three-dimensional forces for stellar and bulge comparisons; enclosed mass alone is insufficient in a bar.
4. Backreaction once deposited mass is no longer negligible. The fixed-field expression is a linear response, not a completed formation model.
5. A specified metric/stress response before converting the same deposits into a lensing prediction. A Newtonian force integral alone does not establish light bending.

This derivation is sufficient to define the next calculational target without claiming a microscopic explanation for why time or gravity exists. It adds no expansion or dark-halo term, and does not validate the still-unresolved photon/time redshift mechanism.
