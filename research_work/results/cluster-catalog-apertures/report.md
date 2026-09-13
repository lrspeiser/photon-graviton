# Full catalog aperture predictions

## Result

All 814 usable THEMIS catalog emitters now contribute to physical spherical and observer-projected aperture calculations. The measured finite M87 approximation replaces its point emitter; the other 813 retain point emission. This extends the common-response calculation beyond M87 without assuming the complete cluster is spherical.

| Radius (kpc) | Mean spherical power (million Lsun) | Mean projected aperture power (million Lsun) | Conditional angular-mean response ratio |
|---|---:|---:|---:|
| 100 | 1.544 | 5.649 | 14.63 |
| 300 | 11.299 | 35.191 | 12.46 |
| 500 | 30.517 | 80.007 | 10.49 |
| 900 | 118.824 | 189.287 | 6.37 |
| 1000 | 217.923 | 217.923 | 4.00 |

At 500 kpc the corresponding M87-only ratio was 5.57. Including other emitters changes the shape relationship, not merely the total power. This is a conditional model prediction, not agreement with measured stellar motions or lensing. The catalog remains incomplete, and no accumulated history is inferred.

## What is averaged

For a generally nonspherical Newtonian source, Gauss's law gives the spherical average inward radial acceleration G M3(<r)/r^2. Define v_monopole^2=r times this acceleration; this is a summary of the field, not a claim that individual stars follow circular orbits.

In the thin-lens equal-potential weak-field branch, the azimuthal mean radial deflection around a circle of impact radius b is 4G M2(<b)/(c^2 b), even when the projected density is not circular. This follows by the two-dimensional divergence theorem applied to the lens potential. Tangential deflection and angular structure are not specified by this average. At b=r, the known relations give

    mean_radial_deflection / (v_monopole^2/c^2) = 4 P2(<r)/P3(<r)

if stored energy shares the calculated power shape with one common accumulation normalization and ordinary mass-energy response. This is the same conditional branch as [the M87 response report](../cluster-m87-finite-emitter/response-report.md), generalized to angular averages. These are established field/geometry identities, not newly derived companion physics. Deflection is physical asymptotic bending, not observed image separation; distance geometry is still needed.

Common exposure and any common constant response multiplier cancel, but spatially varying source histories, motion, stress or unequal potentials can change the prediction. Treating a pressureless reservoir as stationary still requires a support/evolution model. No such model is supplied by this calculation.

## Numerical method and checks

The unchanged three-dimensional directional kernel from cluster-projected-deposition is integrated over a uniformly sampled 1 Mpc sphere centered on the catalog position of M87. Four independent scrambled Sobol sequences use seeds 4701–4704. Each has 131,072 points and a nested 32,768-point comparison. Multiply the mean volumetric power by 4 pi/3; spherical and projected aperture indicators select the appropriate points. Coordinates retain the prior Earth-facing observer basis. No dark-matter profile or lensing normalization is fitted.

M87's point contribution is subtracted at every point and the verified finite-source radial convolution is added. Point singularities for other emitters are integrable in these finite volumes and no sample lands exactly on a source. This does not resolve their pointwise projected peaks or establish convergence for arbitrarily small apertures. The smallest aperture here is 100 kpc.

Every fine-run total differs from the independent emitter/ray total 217,926,911.95 Lsun by less than 0.01187%. The mean differs by -0.00162%. A second check compares spherical apertures against the separately computed ray-based radial catalog profile with the same finite-M87 correction: differences of the sampled means are -0.1454% at 100 kpc and less than 0.0085% at the larger radii.

The 100 kpc aperture is least resolved: its ratio ranges from 14.45 to 14.73 across the four fine runs, while the 500 kpc ratio ranges from 10.4839 to 10.4904. The coarse/fine mean spherical power changes by about 2.1% at 100 kpc, so fine-run scatter alone should not be presented as a rigorous accuracy bound. Replacing the sampled spherical denominator with the independent ray value gives 14.613 at 100 kpc and 10.4863 at 500 kpc. All run values, coarse means and comparisons are retained.

Scramble ranges and refinement are numerical diagnostics, not observational error bars or confidence intervals. Luminosity and distance uncertainties, completeness, source evolution, capture competition, deprojection and field-response uncertainty are not covered. The calculations inherit the original model-based SED luminosities and stipulated published distances. No final holdouts are opened.

## Consequence and next step

Surrounding emitters can substantially alter the relationship between projected and three-dimensional deposits. They cannot be represented solely by multiplying the M87-only profile. We now have catalog-wide aperture targets for the equal-potential branch, but actual motion/lensing fitting still requires measured observables, ordinary matter, source/receiver histories and a physically supported reservoir. Full angular lensing maps, collisions and the other research objectives remain unresolved. All six objectives remain open.

Run `python research_work/results/cluster-catalog-apertures/run.py`, then `python research_work/results/cluster-catalog-apertures/verify.py`.
