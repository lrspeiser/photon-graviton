# J1621+3931: stored mass is weakly identified by the inner fit

13 September 2026. Fixed-scale conditional profile scan; all existing motion bins fitted and catalogue lens angle consumed as calibration. No new outer observations.

## Result in plain language

The inner measurements do not uniquely demand the previously selected enormous reservoir. A capture scale ten times the half-light radius gives virtually the same motion fit as a scale one hundred times that radius, while storing about 486 times less effective mass. Even removing the companion component entirely worsens raw chi-squared by only 0.332, despite using fewer fitting freedoms.

This does not prove that a smaller reservoir forms, nor that companions are absent. It demonstrates that the earlier extreme mass was not a measured requirement of these inner data. The successful outer predictions remain predictions of that particular frozen extreme branch, not of every acceptable companion fit.

## Executed profile scan

| ac/Re | All-motion chi-squared | Companion bending fraction | Stored mass inside grid (solar masses) | Energy equivalent (J) |
|---|---:|---:|---:|---:|
| 0.1 | 3.093851234 | 0.000000 | 0 | 0 |
| 0.3 | 2.948645172 | 0.095112 | 4.5879e+10 | 8.1992e+57 |
| 1 | 3.093851234 | 0.000000 | 0 | 0 |
| 3 | 3.093851234 | 0.000000 | 0 | 0 |
| 10 | 2.762839496 | 0.086619 | 1.7348e+14 | 3.1004e+61 |
| 30 | 2.761992001 | 0.084749 | 3.3482e+15 | 5.9837e+62 |
| 100 | 2.761989333 | 0.084744 | 8.4243e+16 | 1.5056e+64 |

The stellar-only score is 3.093851234. At ac/Re=10 the difference from the extreme free-fit score is only 0.0008502; at 30 it is about 0.00000267. The 100-to-10 mass ratio is approximately 486. Both still represent substantial conditional mass/energy inventories. Zero-amplitude solutions at 0.1, 1 and 3 have no identifiable companion scale. The compact 0.3 branch supplies a separate modest improvement with about 4.59e10 solar masses. The nonmonotonic sequence reflects the changing attenuated density shape and local optimization; it is not a proof of all profile minima between sampled scales.

## Model, fitting and provenance

The retained diagnostic density is rho_d=D J(r;ac,k0)/[1+(r/ac)^2]^2. Attenuation, mass integration, Jeans moments and lens projection are known mathematics; the capture shape and interpretation are project hypotheses. k0 remains fixed. Each chosen ac changes the optical depth as well as the spatial scale. D is fitted through the fraction of lens bending assigned to companions, so this experiment cannot validate the one-third exponent or predict source supply.

At each of seven prespecified scales, refit five parameters: stellar gradient h, central/outer anisotropy, orbital transition radius and deposited bending fraction. Stellar mass is calibrated from the remaining lens bending. The light profile, full covariance, parameter bounds and necessary slope condition are inherited from free-companion.py. Fourteen starts are used initially and fifteen subsequently, including the prior scale solution. All trial success flags and objective values are recorded. Direct profiles replace the earlier capture-scale interpolation table. These are bounded numerical minima found by local searches, not a global optimum proof or confidence interval.

The inventory is the integrated positive density inside the inherited grid. The remaining mass is bounded by 4*pi*D*ac^4/R using J<=1. Conditional energy equivalents use E=M c^2 under the reference effective-density closure, not a unique gravitational-field energy density. No universe age or size is assumed.

## Numerical checks and limits

Zero-companion reproduction, full-covariance decomposition, lens closure and the final 8193-point necessary orbital condition pass. The largest-scale score reproduces the prior free result. All seven scales were refitted after doubling incoming angular quadrature from 96 to 192; the largest score difference is 6.24e-11, and maximum motion-prediction difference is 3.49e-06 km/s. The repeated search supports the reported plateau at its stated precision.

Adaptive bending integration emits a roundoff warning at the requested 1e-9 tolerance for an extended profile. Angular refinement does not independently establish that adaptive tolerance; the shared radial grid, stellar projection and bending method remain inherited approximations. Earlier split-projection refinement tested the extreme branch separately. These checks support a distinction much larger than numerical noise in mass, not exact parameter inference or physical stability. Necessary orbital conditions do not establish a positive distribution function.

## What changes next

The immediate problem is identification of the profile, not automatically a shortfall of 8.42e16 solar masses. Three distinguishable routes remain:

1. Profile the bending fraction and fit stellar populations/orbits with independent constraints. This can determine whether any deposited component is needed without selecting its mass through a weakly constrained radius.
2. Test a physically derived outer stopping or support scale, determined from the incident field or environment. A cutoff chosen only to reduce inventory would be another fit, not a solution; attenuation and lensing must be recomputed together.
3. Compare several near-equivalent inner-fit branches with independent outer shear and tracer measurements. The tiny inner-score differences mean that outer evidence can discriminate far more strongly. Source geometry and environmental matter must be included.

The shared exact-third reference is unchanged. This scan does not solve propagation, timing, capture support or source supply, and it does not establish an advantage over freely fitted dark halos. It prevents treating the extreme optimizer endpoint as a uniquely required physical reservoir.

## Reproduction

Run `python research_work/results/companion-extensions/reservoir-scale-scan.py` and repeat with `--angles 192`. The two JSON files record all seven fits, measured/predicted motions, inventories, necessary-condition checks and optimizer outcomes. The script records its direct inputs and the inherited model-input hashes. The pre-execution scope is in reservoir-scale-scan-protocol.md.
