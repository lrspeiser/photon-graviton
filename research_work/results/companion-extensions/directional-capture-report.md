# Directional capture: disk geometry changes the errors, but does not resolve both signs

13 September 2026. Four geometries, each with two normalization choices. One shared mixture per branch, fit only to the original 89 training galaxies and frozen for 29 validation and 31 test galaxies. All samples exposed; no new observational validation.

## Actual three-dimensional calculation

The incoming far-field bath is isotropic. The capture region is oblate with m squared=R squared+z squared/q squared and opacity kappa=k0 W(n_z)/(1+m squared/a squared)^2. Top and bottom are symmetric, but paths through the disk plane differ from vertical paths. Unlike the earlier spherical calculation, deposited density depends on both radius and height, and its equatorial gravitational force includes nonspherical multipoles. This is an idealized capture geometry tied to the observed stellar disk scale, not a measured disk/bulge absorber or mapped nearby/distant stellar radiation sky.

The four endpoint choices are q=0.5 with W=1 (flat isotropic capture); q=0.5 with W=3 n_z squared (vertical-sensitive capture); q=0.5 with W=1.5(1-n_z squared) (side-sensitive capture); and q=0.25 with W=1 (thinner isotropic capture). Each W has angular mean one. W represents an unproved directional interaction response, not a claim that more photons are observed arriving from those directions. The same W enters both local absorption and path attenuation.

    rho_q(R,z) = 2 C0 eta(X) [1+m squared/a squared]^-2 <W exp(-tau)>
    tau = integral_along_incoming_path kappa dl
    rho_mix = (1-f_mix) rho_reference + f_mix rho_endpoint

The exact-third eta and all original reference constants remain unchanged. Ray integration, exponential attenuation, angular quadrature, linear source superposition and the Legendre solution of Poisson gravity are known mathematics. The ellipsoidal opacity, W choices and storage interpretation are proposed hypotheses. Standard nonspherical potential methods are described in [Tremaine, potential theory lectures](https://www.ias.edu/sns/tremaine/lectures/ast513/potential).

Two normalization choices separate effects. Raw branches retain the calculated capture-derived density amplitude. Fixed-mass branches rescale the endpoint density to the spherical inventory over the same numerical domain: this isolates redistribution, but is not derived from the transport law. That deterministic per-galaxy rescaling uses no measured rotation residual, yet it inherits the previously fitted reference inventory. The global mixture is an additional fitted parameter, not a derived capture fraction.

## Frozen-mixture rotation scores

Lower RMSE is better. Values are equal-galaxy km/s errors. Fitting minimizes mean squared log10 speed, so a km/s gain need not improve the actual objective. Eight searched families and one additional parameter must be acknowledged.

| Branch | Fitted mixture | Training RMSE | Validation RMSE | Test RMSE |
|---|---:|---:|---:|---:|
| Spherical reference | 0 | 29.025 | 32.495 | 23.591 |
| flat_isotropic_raw | 0.00000 | 29.025 | 32.495 | 23.591 |
| flat_isotropic_fixed_mass | 0.01926 | 29.145 | 32.581 | 23.728 |
| flat_vertical_raw | 0.00028 | 29.024 | 32.493 | 23.589 |
| flat_vertical_fixed_mass | 0.01189 | 29.079 | 32.537 | 23.643 |
| flat_sides_raw | 0.01067 | 28.969 | 32.437 | 23.531 |
| flat_sides_fixed_mass | 0.02724 | 29.157 | 32.594 | 23.730 |
| thin_isotropic_raw | 0.03073 | 28.783 | 32.230 | 23.334 |
| thin_isotropic_fixed_mass | 0.03791 | 29.454 | 32.760 | 24.069 |
| Existing fitted simple MOND | Not this mixture | 19.890 | 26.876 | 16.398 |

## Inner versus outer signed error

Negative means too slow; positive means too fast. Bins are below one, one to below three, and at least three stellar disk scale lengths. These are not necessarily the outer radii of the companion reservoir. Each galaxy is weighted equally within a bin.

| Branch | Inner mean km/s | Middle mean km/s | Outer mean km/s |
|---|---:|---:|---:|
| Spherical reference | -6.726 | -7.486 | +10.412 |
| flat_isotropic_raw | -6.726 | -7.486 | +10.412 |
| flat_isotropic_fixed_mass | -6.651 | -7.180 | +10.815 |
| flat_vertical_raw | -6.726 | -7.488 | +10.407 |
| flat_vertical_fixed_mass | -6.699 | -7.370 | +10.593 |
| flat_sides_raw | -6.738 | -7.545 | +10.233 |
| flat_sides_fixed_mass | -6.650 | -7.164 | +10.861 |
| thin_isotropic_raw | -6.758 | -7.686 | +9.642 |
| thin_isotropic_fixed_mass | -6.377 | -6.127 | +11.919 |

The fixed-mass flattened branches can raise inner speeds, but also raise outer speeds that were already too high. The raw thinner branch lowers the outer excess while making the inner/middle shortage slightly worse. Thus direction-dependent capture can alter the ratio, but these simple variants do not repair both errors together. The trained mixtures remain small; the fit does not support replacing the entire spherical component by these flattened endpoints.

The summary JSON also reports signed paired logarithmic-loss changes, per-split improvement counts and descriptive galaxy-bootstrap intervals. These do not include full distance, inclination, stellar-mass or selection uncertainty, and are not blind significance or model-selection evidence. Numerical score improvements are not proof of a physical capture mechanism.

## Inventory and energy

| Endpoint | Raw total / reference total, across galaxies |
|---|---:|
| flat_isotropic | 0.5024 to 0.5455 |
| flat_vertical | 0.4993 to 0.5108 |
| flat_sides | 0.4993 to 0.5139 |
| thin_isotropic | 0.2523 to 0.3071 |

The raw endpoints store less, rather than creating or deleting energy. Relative to the spherical capture assumption, the uncaptured fraction remains in traveling companions in the bookkeeping interpretation; a common self-consistent source history and collision law are not derived. Mixture inventories are (1-f_mix)+f_mix times the raw ratio and are recorded separately. Fixed-mass branches preserve the numerical inventory by construction but require redistribution/support that has not been calculated. No cooling spectrum, stability, vertical stellar-motion or lensing success follows from this comparison.

## Matched inventory-only control

As a post-fit diagnostic, reduce the spherical density by exactly the same per-galaxy total-inventory ratio as each raw mixture, with no further fitting. This separates lower retention from directional shape. It is not an additional trained competitor or an untouched test. The summary JSON retains all split speed/log errors.

| Raw branch | Directional validation / test RMSE | Same-inventory spherical validation / test RMSE |
|---|---:|---:|
| flat_isotropic_raw | 32.495 / 23.591 | 32.495 / 23.591 |
| flat_vertical_raw | 32.493 / 23.589 | 32.493 / 23.589 |
| flat_sides_raw | 32.437 / 23.531 | 32.419 / 23.504 |
| thin_isotropic_raw | 32.230 / 23.334 | 32.190 / 23.243 |

For the thin raw branch, the directional test RMSE is 23.334 km/s versus 23.243 for its same-inventory spherical control. The directional branch has slightly better logarithmic scores than that control, so the verdict depends on the stated metric; it does not establish a general advantage. Against the original reference, its validation logarithmic score slightly worsens, and descriptive paired log-loss intervals span zero in all three splits. The speed gain therefore cannot be cleanly attributed to a successful directional mechanism.

## Verification and limitations

Coarse resolution was {'nr': 384, 'nt': 24, 'nv': 16, 'nphi': 32, 'lmax': 12}; refined resolution was {'nr': 768, 'nt': 48, 'nv': 32, 'nphi': 64, 'lmax': 24}. Mixtures were not refitted on comparison data or at higher resolution. The largest frozen-mixture speed change on refinement is 0.029255 km/s. The maximum reference force-reconstruction speed discrepancy is 0.027225 km/s; the mixture baseline itself uses the exact archived spherical predictions.

Independent direct numerical ray integrals agree with the analytic expression within 1.2e-14 relative error over 45 cases. Smooth oblate density fields checked against a separate homoeoid integral agree to 0.000504 relative error at refined resolution. These validate the numerical machinery, not the proposed physics.

The grid extends from 1e-6 to 3000 reference capture radii. Positive exterior density bounds give omitted mass at most 0.0718% of the grid inventory. This finite domain is a numerical approximation, not an assumed age or universe size. Fixed-mass normalization is exact only within this common finite domain.

The model is axisymmetric and reflection symmetric. It does not model the Milky Way bar, independent bulge opacity, individual stellar source positions or environmental asymmetry. Direction-sensitive cross sections and a long-lived deposited state are assumptions. Rotation was computed in the equatorial plane; predicting motion above/below the disk and light bending requires additional calculations. No outward-force term was clipped to simulate a positive circular speed.

## Consequence

Keep the original shared law as the reference. Assess the small raw thinner-capture improvement against the inventory-only control, rather than attributing any gain automatically to direction. None of these fitted mixtures repairs both radial signs. Further work should specify a physically motivated disk/bulge capture relation or use vertical-motion constraints to distinguish geometries, rather than enlarge this sweep solely to find a lower residual.

The user also proposed making capture respond to stronger gravity deeper in the galaxy. This is distinct from the present prescribed opacity profile, although that profile already has greater local opacity toward the center. A gravitational-depth rule must distinguish potential depth from local acceleration, calculate the incoming energy surviving interception along each path, and include feedback from deposits. A deep central well can have small local acceleration by symmetry. Greater central deposition still contributes to outer gravity. Earlier binding-feedback branches should be reviewed before constructing a directional version, rather than labeling all depth feedback untested or repeating an old prescription. No depth-feedback experiment was performed in this run.

## Reproduction

Run directional-capture.py, directional-capture.py --refine, directional-capture-check.py and directional-capture-report.py in that order. Protocol, source hashes, density-derived force arrays, inventories, scans, frozen predictions, refinement and summary are saved alongside this report. The paper supplement indexes this result; PDF v1.5 predates it.
