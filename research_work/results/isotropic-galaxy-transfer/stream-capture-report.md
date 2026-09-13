# Do companion rivers deposit the required halo shape?

Fixed capture parameters and the exact one-third law were used throughout. Straight incoming rays were attenuated along their paths inside each finite 5 Re diagnostic region. This is a conditional external-input shape test, not an absolute energy supply test or a new fit to observed lens images or stellar velocities.

## Result

Under the standard comparison geometry, source mixing reduces cumulative-profile differences from 14.46 to 2.44 percentage points for J1112, from 18.83 to 5.52 for J1621, and from 18.91 to 5.83 for J1630 (Chabrier proxy). It does not improve J0037, J1204 or J1402: their differences remain 17.97, 49.94 and 65.11 points. Their required profiles are more centrally concentrated than these straight-stream deposits. Under our retained geometry none of the six source mixtures improves on distant illumination in this basis. These are conditional shape comparisons, not observed lens residuals.

The successful shape adjustments require the nearest tested shell: J1112 assigns about 93.1% of deposited power to D/R=1.1, while J1621 and J1630 put all power there. Thus they hit the source-distance boundary and have no independently measured source justification. Sources closer than 1.1 R, internal sources, different capture regions or subsequent migration were not tested; no general exclusion follows.

## Executed transport

Known transport: dF_c/ds=-kappa F_c after geometric dilution is factored out, with tau=integral kappa ds and q_ret=eta kappa F_c. Proposed capture profile: kappa=k0/[1+(r/a)^2]^2, retaining the previously calibrated k0 and a. The Hill expression eta=X^(1/3)/(1+X^(1/3)) is known mathematics with our proposed retention interpretation. Unretained absorbed energy belongs to another reservoir; no energy creation is assumed.

For a point source, F_c(x)=L[1-exp(-alpha*d_entry)] exp(-tau)/(4*pi*d(x)^2). The entry distance is source to boundary; d(x) is source to deposition point. These rays include finite-distance divergence. The external companion component is supplied at the boundary. Conversion inside the region and internal stars are excluded from this particular test. Exterior opacity is unspecified; the boundary field is conditional. R=5 Re is a diagnostic boundary, not a measured capture cutoff.

For each of 24 geometry/population combinations, sources were placed at 1.1, 2, 5, 10 and 100 boundary radii, plus a distant isotropic field. Spherical averages also represent uniformly populated source shells. A single beam can make an asymmetric deposit; a shell averages that asymmetry away. The one-third factor multiplies all rates and cancels only when comparing normalized radial shapes. An isotropic distant field has zero dipole.

## Halo shape comparison

Both curves are normalized to their mass or deposited energy inside R. The error below is the largest difference between cumulative fractions, in percentage points; it is not a velocity error, lens error, confidence level or statistical significance. Constant illumination history with no migration is required to identify deposited-power shape with stored-energy shape.

| Galaxy | Geometry | Distant field error, pp | Best distance-mixture error, pp | Target half-mass radius / R | Distant half-deposit radius / R |
|---|---|---:|---:|---:|---:|
| J0037-0942 | companion_regular | 62.25 | 62.25 | 0.112 | 0.592 |
| J0037-0942 | standard_flat_FLRW | 17.97 | 17.97 | 0.480 | 0.578 |
| J1112+0826 | companion_regular | 5.98 | 5.98 | 0.570 | 0.554 |
| J1112+0826 | standard_flat_FLRW | 14.46 | 2.44 | 0.646 | 0.540 |
| J1204+0358 | companion_regular | 49.43 | 49.43 | 0.167 | 0.526 |
| J1204+0358 | standard_flat_FLRW | 49.94 | 49.94 | 0.156 | 0.512 |
| J1402+6321 | companion_regular | 56.57 | 56.57 | 0.148 | 0.577 |
| J1402+6321 | standard_flat_FLRW | 65.11 | 65.11 | 0.077 | 0.563 |
| J1621+3931 | companion_regular | 16.47 | 16.47 | 0.491 | 0.578 |
| J1621+3931 | standard_flat_FLRW | 18.83 | 5.52 | 0.700 | 0.564 |
| J1630+4520 | companion_regular | 7.63 | 7.63 | 0.562 | 0.565 |
| J1630+4520 | standard_flat_FLRW | 18.91 | 5.83 | 0.690 | 0.551 |

The mixture is deliberately generous: six nonnegative deposited-power weights are fitted separately to each target. It is an inverse diagnostic, not a measured source population or shared prediction. No source luminosities, durations, or energy budgets follow from these weights. Both population cases, every distance, angular dipoles, normalized curves and weights are retained in the JSON.

These NFW profiles are our restricted target fits, not unique measured halos. Some have parameter boundaries and poor stellar fits. Standard geometry remains a comparison only; original capture scales are retained when evaluating that target, so it is not a recalibration under FLRW. Matching this enclosed spherical profile would still not prove projected lens agreement, support, or a viable energy history.

## Interpretation and next decision

Changing source distance changes which paths illuminate the receiver, but leaves its local capture coefficient fixed. A finite smooth capture law with nonsingular external sources gives finite central deposited density, whereas NFW has a central cusp. This does not alone exclude agreement over a finite measured range, but source mixing cannot arbitrarily prescribe the inner profile. The table quantifies the discrepancy over this declared diagnostic region.

No bending, scattering, binding, or subsequent inward migration was added. If distance mixtures leave a deficit of central deposition, the next physical variable is an explicit transport/capture change that delivers energy inward. It must account for momentum and retained-state support; simply renaming rays rivers does not supply it. A target-fit mixture is not a reason to change the one-third exponent.

Actual nearby-source positions for these six receivers are not available in this calculation. The cached M87 catalog was not reassigned to them. Thus this test is a controlled distance-family test against existing halo targets, not an observed source-to-halo prediction. No new galaxy, cluster, redshift, or lensing success is claimed.

## Verification

Independent optical-depth quadrature agrees within 7.33e-15. Doubling radial/angular resolution changes cumulative fractions by at most 0.00126; the D/R=100000 limit agrees with distant illumination within 3.77e-07. All cumulative profiles are monotonic, and mixture weights are nonnegative and sum to one. See stream-capture-protocol.md, executable stream-capture.py and input hashes in the results. All six project goals remain open.
