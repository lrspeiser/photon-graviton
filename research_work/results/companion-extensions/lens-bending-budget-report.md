# Where the six-lens bending discrepancy originates

13 September 2026. Decomposition at fixed stellar nuisance fits and conditional geometry.

## Main result

In three of the six systems, the fitted stellar component alone already exceeds the bending needed at the catalogue lens angle. This is true under both population proxies and all four tested profiles: original reference, amplitude-adjusted reference, local filling and full recycling. The affected systems are J0037-0942, J1204+0358 and J1402+6321.

For the adjusted reference's Chabrier case, stars alone supply 114.39%-117.65% of the required bending in those systems. Even removing every companion contribution would leave too much bending. This is a conditional incompatibility of the fitted stellar model, geometry and target angle; it is not an independent measurement that the real stars are too massive. The stellar masses were inferred from inner motions under the same model assumptions.

Consequently, changing only nonnegative companion density while holding these stellar fits and geometry fixed cannot resolve those three cases. A joint refit of the stellar component, orbital structure or geometry could change that conclusion and has not been performed here.

## Formula and provenance

At the catalogue angle theta_E, let b_E=D_l theta_E and define the geometry-scaled stellar and companion deflections

\[
\theta_\star=\frac{D_{ls}}{D_s}\widehat\alpha_\star(b_E),\qquad
\theta_c=\frac{D_{ls}}{D_s}\widehat\alpha_c(b_E).
\]

The known spherical weak-field mass integral supplies each deflection, consistently with the previous lens calculation. The required lens equation at this impact parameter is theta_star+theta_c=theta_E. Define S=theta_star/theta_E and Q=theta_c/theta_E. If only the companion amplitude is allowed to change, its inverse-required multiplier is

\[
\boxed{m_c=\frac{1-S}{Q}.}
\]

This is algebraic inversion of the existing effective lens equation, not a new gravity law or a fitted predictive parameter. For positive Q, S>1 requires m_c<0. Such a negative value is a diagnostic flag: it is not adopted as negative energy or repulsive companion gravity.

At the same fixed geometry and stellar light-profile shape, the largest stellar mass compatible with the target after removing companions is

\[
M_{\star,\max}=\frac{M_{\star,\rm fit}}{S}.
\]

If the companion contribution is instead held fixed, the required stellar-mass multiplier is (1-Q)/S. These relations quantify conditional changes only. They do not show that changing those masses would still fit the stellar motion data; that would require another explicit nuisance fit and consistency check.

## Per-system diagnostic

Values below use the refined amplitude-adjusted reference, Chabrier proxy. Fractions refer to bending at the fixed catalogue angle; they are not fractions of the model's predicted Einstein radius and should not be substituted into the previously reported angle RMS.

| System | Stellar bending / required | Companion bending / required | Required companion multiplier | Minimum stellar mass reduction even with companions removed |
|---|---:|---:|---:|---:|
| J0037-0942 | 1.14922 | 0.08387 | -1.77918 | 12.98% |
| J1112+0826 | 0.87935 | 0.09212 | 1.30969 | 0% |
| J1204+0358 | 1.14388 | 0.07057 | -2.03867 | 12.58% |
| J1402+6321 | 1.17651 | 0.07465 | -2.36463 | 15.00% |
| J1621+3931 | 0.98200 | 0.08874 | 0.20280 | 0% |
| J1630+4520 | 0.85039 | 0.10990 | 1.36128 | 0% |

A zero in the last column means stellar bending alone does not exceed the target; it does not mean the full profile fits. J1621+3931 would require much less companion bending at its current stellar fit, whereas J1112+0826 and J1630+4520 would require more. This mixed pattern explains why one universal increase or decrease of companion strength is inadequate at fixed nuisance fits.

Both population cases give three stellar-only excesses in every branch. Across systems, the adjusted-reference stellar fraction ranges from 0.85039 to 1.17651 (Chabrier) and 0.84849 to 1.17384 (Salpeter). Local filling reduces the maximum to 1.14530/1.14351, and full recycling gives 1.15121/1.15001, but neither removes the three excesses. The complete 48-case decomposition is recorded, including inverse-required multipliers and conditional stellar-mass ceilings.

## What can change next

This diagnostic separates several physical and modeling possibilities:

- Stellar mass normalization and radial population gradients: the current constant mass-to-light profile may allocate too much projected mass inside the lens aperture. Changing it must still match stellar motions and available population information.
- Orbital structure and geometry: constant anisotropy and spherical modeling determine which stellar mass is inferred from inner velocities. A more general model must be constrained by the data rather than chosen solely to fix the angle.
- Distance/optical mapping: the adopted nonexpanding geometry remains conditional. Any change must also recompute physical light-profile radii and dynamical fits, and satisfy the propagation tests; an arbitrary per-lens distance factor is not a theory.
- Joint redistribution: a changed companion profile could alter the stellar mass inferred from the inner motions. That requires refitting the coupled motion/lensing problem. The present fixed-stellar argument does not exclude that possibility, but simply adding positive projected companion mass cannot lower the bending floor while holding stars fixed.

These are alternatives to assess, not newly selected solutions. The result also applies to any additional nonnegative spherical gravitating component at these fixed stellar fits, so it is not evidence against dark matter specifically or for a special companion effect.

## Scope and checks

The catalogue angles are SIE model summaries; no raw-image likelihood or alternative image reconstruction is fitted. The same conditional angular distances and source/receiver ratios are used as in the previous lens comparison. Luminosities remain population-mass-derived proxies. No observational confidence interval or independent stellar mass measurement is supplied by this inversion.

Run `python research_work/results/companion-extensions/lens-bending-budget.py`. It reconstructs the stellar deprojection and companion profiles using the original/refined numerical grids as appropriate, but performs no new fit. Separately integrating stellar and companion bending reproduces every existing predicted lens root to relative error below 3.93e-10. This verifies the decomposition against the implemented lens equation, not the physical validity of its inputs. Raw/derived input and code hashes are recorded in `lens-bending-budget-results.json`.

The main conclusion is conditional and specific: half of these lens targets require changes beyond companion-only amplitude/redistribution at the existing stellar fits. The broader theory remains incomplete, and the next joint model must address that stellar/geometry dependency explicitly.

Related: [matched capacity lens transfer](capacity-lensing-report.md), [earlier target-profile compatibility](../isotropic-galaxy-transfer/lens-profile-compatibility.py).
