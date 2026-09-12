# Radial capture changes the deposited profile

The initially bound weak-opacity source does not stay near its capture positions. A small continuously injected population becomes substantially more centrally concentrated as particles move through the well. Its kinetic energy increases through infall while its orbital energy remains conserved. Thus permanent deposition at the capture point is not the dynamical consequence of the tested cold-merger rule.

## Scope and normalization

This is the linear test-population response in a frozen spherical potential. Use the existing p6,C1,kappa0=.1,R100 snapshot with c/v0=1000. Earlier work found its sampled instantaneous cold source to be initially bound. Ordinary matter and the old deposited profile supply a fixed background potential; we do not simultaneously claim that the old profile is the result of these new trajectories.

The added source has amplitude eta, taken infinitesimal for this diagnostic. Its shape and injection velocities come from the previous capture calculation. All normalized mass fractions and mean specific energies are independent of eta. Reported source rates and accumulated masses are coefficients per unit eta; they are not an actual massive addition for which self-gravity has been omitted. Eta also absorbs the unspecified conversion between the earlier normalized exposure rate and dynamical time. No physical luminosity or age is assigned.

Finite accumulation requires changing self-gravity, source attenuation and the background response. Central concentration can make those effects important locally before the added mass becomes large globally. The test-population approximation is not a justification for neglecting such effects indefinitely or arbitrarily close to the center.

## Equations and provenance

Newtonian radial trajectories and their energy integral are known mathematics. Applying them to the proposed cold-capture source is a conditional project calculation, not a novel gravity formula. Define a signed coordinate x along each orbit's line through the center:

`dx/dt=v`, `dv/dt=-sign(x) g(|x|)`,

`g(r)=[M_b(r)+M_D(r)]/r^2`, `e=v^2/2+Phi(|x|)`.

Each particle starts at r0 with inward velocity `v0=-c beta(r0)`. Its rest-mass weight is the captured energy weight multiplied by `sqrt(1-beta^2)`, under the earlier merger postulate. All launch nodes are bound in the frozen field. The signed-coordinate treatment allows genuine center crossing onto the opposite side of the orbit; there is no sticky center, reflecting wall, circularizing torque or cooling.

For constant injection, a source launched at all times produces the age integral:

`delta M(<r,T)=eta integral_0^T da integral d r0 S_rest(r0) 1[|x(a;r0)|<r]`.

Dividing by `eta T integral S_rest dr0` gives the reported enclosed fractions. This is a transport response, not a final equilibrium density. Permanent in-place capture would instead retain the source's initial radial fractions at every age.

Time units are a/v0, radii are in a, and specific energies are in v0^2=G M_b/a. These are not kpc or billions of years until physical scales are specified.

## What the calculation establishes

The innermost fraction increases greatly, so deposited gravity cannot be calculated from capture positions alone under this particle interpretation. At fixed added total mass, an enclosed fraction directly determines the added spherical acceleration at that radius. This does not mean that the reported acceleration change has been fed back into the trajectories; that is the next nonlinear calculation.

Kinetic energy grows at the expense of potential energy. Therefore a low injected kinetic energy is not, on its own, proof that orbital support lacks an energy source. Conversely, energy conservation alone does not supply tangential motion: all orbits here still have zero angular momentum. A spherical ensemble with radial velocities can be statistically spherical while being strongly anisotropic in velocity. No collective stability or stationary limit is established from three finite exposure times.

The same trajectories also determine projected enclosed mass. A spherical shell at radius r contributes fraction1 inside projected radius b>=r, or `1-sqrt(1-b^2/r^2)` for b<r. This is known geometry; it is integrated over the same source and age weights. If the deposited component follows the usual weak-field spherical light-deflection law, its contribution is `delta alpha(b)=4 G delta M_projected(<b)/(c^2 b)`. That lensing law is an explicit conditional choice, not a new formula or a validation of the temporal-field model. The projected fractions therefore provide a linked lensing input without identifying them as observed deflections. Force and projected-mass enhancements refer only to the added component, not the entire galaxy.

## Verification and reproducibility

The [protocol](protocol.md) predates calculation. A kick-drift-kick leapfrog integrates all launch radii simultaneously; a trapezoidal age integral represents continuous injection. Coarse512-shell fields use two launch nodes per shell and dt=.004; a dt=.002 control isolates time refinement. The fine1024-shell calculation uses four launch nodes per shell and dt=.002. Launch quadrature is in r^3 so shell-volume weights are represented directly.

All declared mass-fraction, kinetic-energy, orbital-energy and mass-normalization gates pass. The largest enclosed-fraction refinement difference is0.0002065 absolute, or0.02065 percentage points; the largest mean kinetic-energy difference is0.01468 percent. The fine run's maximum individual orbital-energy error is7.48e-7 in v0^2 units, below the1e-4 gate. No mass is discarded at center crossings or at a numerical outer wall; the external deposited potential continues as a point-mass potential outside the source domain.

Run `python research_work/results/radial-capture-response/run.py` with OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 for practical runtime. It reads retained source profiles and no observations. The code and results.json preserve each resolution and every comparison. `export.py` writes the comparison table and checkpoint below.

## Next requirement

An adopted particle branch must evolve its source and retained population together, including changing gravity and escape. The initial-motion calculation and this response kernel provide explicit starting equations, but not a completed galaxy explanation. A different stored-field branch must specify its stress and energy law rather than borrowing orbital support implicitly. All nine goals remain active, including the unresolved photon/time/redshift mechanism and the joint observational and withheld tests.


## Final checkpoint

| Population | Inside r0.1 | Inside r1 | Inside r3 | Inside r10 | Mean kinetic energy per added rest mass |
|---|---:|---:|---:|---:|---:|
| Source at injection | 0.02905% | 16.11324% | 69.77784% | 97.50857% | 0.0488102 |
| Moving, T=1 | 1.18917% | 25.17914% | 71.44785% | 97.50950% | 0.0993956 |
| Moving, T=10 | 2.99945% | 33.31296% | 82.85856% | 97.67979% | 0.1770752 |
| Moving, T=30 | 2.78729% | 31.42511% | 80.06104% | 98.19142% | 0.1650653 |

The source row gives the spatial baseline for permanent in-place capture and the kinetic energy at injection. It is not a stationary population with that kinetic energy.

At T30 the innermost enclosed fraction is 95.96 times the in-place value. The mean mechanical energy remains approximately -0.502400 in all moving-population samples; kinetic growth is balanced by a more negative potential energy. The finite-time response is not a stationary solution.

[comparison.csv](comparison.csv) retains all values. The accumulated mass coefficient at T30 is approximately10.2752 per unit source amplitude; actual added mass is eta times this coefficient. Setting eta=1 and treating this as negligible added mass would be inconsistent with the approximation.


Projected enclosed fractions (the input to a conditional spherical weak-lensing calculation):

| Population | Inside b0.1 | Inside b1 | Inside b3 | Inside b10 |
|---|---:|---:|---:|---:|
| Source at injection | 0.42676% | 29.32198% | 79.04967% | 98.45247% |
| Moving, T=1 | 2.27108% | 36.82130% | 80.01768% | 98.45680% |
| Moving, T=10 | 4.74024% | 45.68769% | 87.27735% | 98.53978% |
| Moving, T=30 | 4.41425% | 43.72108% | 86.67587% | 98.81204% |

At fixed added mass, the innermost projected fraction grows 10.34-fold by T30, compared with 95.96-fold for true spherical enclosed mass. These are changes in the added component, not multipliers for the total galaxy force or lensing.
