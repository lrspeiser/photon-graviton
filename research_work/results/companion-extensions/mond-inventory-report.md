# Real SPARC data: MOND-guided redistribution improves rotation predictions

13 September 2026. All 149 saved SPARC galaxies, 3150 accepted measured radii, original 89 training / 29 validation / 31 test split. Observations were previously exposed. Original source-data hashes were verified against the third-radiation-retention manifest. No galaxies, stellar mass-to-light ratios, distance inputs or observed speeds were changed.

## Result

A shared mixture trained on proportional speed errors selects **f=0.923557**. It shifts the existing companion density toward a positive, finite-inventory spherical equivalent of the simple MOND acceleration target. This substantially improves the original rotation prescription, but is not a derived capture law or independent discovery of MOND behavior.

| Model | Training RMSE km/s | Validation RMSE km/s | Test RMSE km/s |
|---|---:|---:|---:|
| Original exact-third companions | 29.025 | 32.495 | 23.591 |
| Known simple MOND, training-fitted a0 | 19.890 | 26.876 | 16.398 |
| Fully redistributed, inventory-capped target | 20.474 | 27.384 | 16.604 |
| Training-selected shared mixture | 20.468 | 27.136 | **16.241** |

The selected test RMSE improves about **31.2%** over the original companions. Its slight km/s advantage over simple MOND is not a superiority claim: MOND has better logarithmic error in all three splits and better training/validation RMSE. The shared mixture uses three inherited companion constants, the inherited MOND a0 and one fitted mixture parameter; the simple MOND control has only its fitted a0. No uncertainty or significance estimate supports the small test RMSE difference.

| Model | Training log RMS | Validation log RMS | Test log RMS |
|---|---:|---:|---:|
| Original companions | 0.138489 | 0.115394 | 0.091058 |
| Simple MOND | 0.109172 | 0.095427 | 0.078337 |
| Selected mixture | 0.119156 | 0.104234 | 0.083155 |

## Formula and provenance

Known simple MOND mathematics supplies the extra acceleration target:

    gc = [sqrt(gb^2 + 4 a0 gb) - gb] / 2
    a0 = 8.563335193921255e-11 m/s^2

This a0 was previously fitted on the training sample and was not retuned. The original disk-and-gas acceleration gb remains unchanged. The scalar MOND formula is the existing rotation prescription, not a solved three-dimensional modified-gravity field.

Known Newtonian spherical inversion gives:

    Mtarget(<r) = r^2 gc(r)/G

The hypothetical fixed-inventory redistribution rule is:

    Mend(<ri) = min(Mtotal, max over j<=i of Mtarget(<rj))
    Mnew(<r) = (1-f) Mreference(<r) + f Mend(<r)
    vnew(r) = sqrt(vb(r)^2 + G Mnew(<r)/r)

Mtotal is calculated by integrating the original exact-third capture density, including its outer tail. It is not fitted to each galaxy's observed speed or target halo. The cumulative maximum ensures a nondecreasing enclosed-mass endpoint; the cap prevents borrowing more mass-energy than the reference inventory. Neither operation is a microscopic capture or migration derivation. The MOND target is explicitly borrowed physics, while interpreting its finite-inventory redistribution as companions is the project hypothesis. The one-third factor survives in Mtotal and the original component; it no longer sets most of the fitted radial shape.

All mass remaining beyond the final measured radius can be given a positive exterior completion, for example in a shell from 2 to 3 times that radius. That construction proves a positive finite-inventory spherical completion exists, but its stopping radius, width and support are arbitrary here. No particular exterior shell was observed or tested with lensing. It uses the measurement boundary only for a diagnostic completion, not a proposed universal physical galaxy boundary. Redistribution includes both outward and inward movement, and f is a profile mixture fraction, not an inferred fraction moving solely inward.

## Positivity and inventory

No adjacent sampled target-mass interval decreases in any of the 149 galaxies. Thus no running-maximum repair was necessary at these radii. A piecewise-positive interpolation can fill the shells; this does not establish the density at unmeasured radii or verify a full nonspherical reconstruction.

**21 galaxies, covering 151 measured radii, exceed the original inventory in the uncapped MOND target.** They remain in every score, with the cap applied. NGC3741 needs approximately 39.33 times its current total model inventory to reach the uncapped target at its most demanding sampled radius; UGC05721 needs 8.27, DDO154 6.43, NGC2915 4.67 and UGC04483 3.92. These are shortfalls of this specified model inventory, not an absolute cosmic photon-supply impossibility. The target is capped instead of adding energy.

The endpoint inventory fraction outside the last measured star ranges from 0 to 99.04%. An exterior reservoir is therefore compatible with the rotation diagnostic, but its lensing, energy source and confinement have not been established. Orbital/binding work during redistribution still needs an energy ledger even though the deposited rest-energy inventory is held fixed.

## Inner/outer behavior

Each galaxy contributes equally within its radial bin. Negative means speeds too low; positive means too high.

| Region in stellar disk scale lengths | Original mean error km/s | Selected mean error km/s |
|---|---:|---:|
| Inner, below 1 | -6.726 | +1.912 |
| Middle, 1 to below 3 | -7.486 | -4.991 |
| Outer, at least 3 | +10.412 | -4.831 |

The modification reduces the original inner deficit and outer excess, but leaves a middle/outer deficit. It is a material improvement over the original shape, not a complete fit.

## Verification and implications

Zero mixture reproduces archived companion velocities within 1e-10 km/s. Raw MOND reproduces the saved control within 1e-9 km/s. Total-inventory quadrature at 96 and 192 incoming angles agrees within 3.47e-7 fraction; selected speeds change by at most 0.000025 km/s. Positive shell increments and inventory caps pass. These verify the numerical calculation, not observational uncertainties, physical energy supply or dynamical stability.

Keep this as a promising MOND-guided redistribution diagnostic alongside the original reference. It demonstrates that a fixed reference inventory can support substantially better sampled rotation predictions for this ensemble, with specific supply-limited cases. It does not demonstrate that incoming waves naturally form the required profile, nor that the mechanism predicts lensing, vertical gravity, clusters or redshift. A next meaningful test is to derive an acceleration-dependent capture/retention/transport rule that produces a similar profile without using the target as its definition, and check the finite exterior completion against independent outer gravity or lensing.

Reproduce with `python research_work/results/companion-extensions/mond-inventory.py`. Predictions, the full training scan, positivity/inventory diagnostics and source hashes are in mond-inventory-results.json. The v1.5 PDF predates this report.
