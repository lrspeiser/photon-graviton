# Where the catalog-fed cluster deposition occurs

## Result

The same catalog and capture assumptions now determine a radial cumulative deposit-power profile without fitting a desired gravity distribution. External illumination deposits mostly in outer layers, while the internal sources add substantial inner deposition.

| SED model | Internal deposits outside R/2 | External deposits outside R/2 | Combined deposits outside R/2 |
|---|---:|---:|---:|
| themis | 64.752% | 95.168% | 85.994% |
| dl14 | 64.757% | 95.179% | 86.007% |

The outer half of the radius occupies 87.5% of a sphere's volume. Thus the external component is preferentially outer-loaded relative to uniform deposition, but the combined 86.0% fraction is slightly below that volume reference. A large percentage near the edges is not automatically a greater local density there. This cumulative comparison does not establish monotonic local density or the location of a projected lensing peak.

## Frozen physical assumptions and provenance

Inputs are the preceding cluster-catalog-pilot/results.json, hashed in the output. Both 814-source SED variants retain the same NGC4486 center, stipulated R=1 Mpc boundary, archived empirical alpha=0.0002488993265191759 per Mpc, and hypothetical kappa=10 per Mpc. There is no parameter fitting in this calculation. Published distances are stipulated; SED luminosities, point emitters, steady illumination, no intervening capture, permanent stationary retention, and neglected backreaction retain all the preceding report's limitations. No actual cluster age, accumulated energy or target mass is introduced.

Known Euclidean ray geometry and the known coupled linear transfer equations give the deposited fraction D(s) along each ray, as derived in the preceding internal/external source reports. For a radial aperture r0, a ray's distance from the center obeys

    r(s)^2 = b^2 + (s-s_closest)^2.

The portion inside the aperture is the intersection of [0,L] with

    [s_closest-sqrt(r0^2-b^2), s_closest+sqrt(r0^2-b^2)]

when r0>=b; otherwise it contributes zero. Its deposited fraction is D(s_hi)-D(s_lo). Integrating over the emitter's ray directions and multiplying by its luminosity produces the enclosed deposited power. This uses established geometry and transfer mathematics conditional on our postulates, not a new fundamental equation. The resulting radial profile is calculated from source geometry, not copied from a dark-matter profile.

## Numerical checks

All 814 sources are integrated at 256 and 512 angular nodes for each SED version, at 21 radial apertures. The maximum change in cumulative power between resolutions is 5.67e-5 of the total. The final aperture recovers the preceding independent total-power calculation to within 1e-14 relative. Cumulative power begins at zero, remains nonnegative and increases monotonically. All source records are retained in the upstream hashed input. The reported refinement is an estimate of integration sensitivity, not an observational uncertainty or rigorous error bound.

## Consequences for the six objectives

This closes a spatial-accounting step but is not a full 3D deposit field. It gives radial cumulative power for the fixed illumination, not the supported energy density after dynamical evolution. Internal and external source contributions are both included, with no independent density fitting.

A spherical monopole approximation could use cumulative stored energy to predict angle-averaged radial force under an explicitly chosen gravity law, but the actual source illumination is anisotropic. A lensing map requires the spatial distribution and metric response; radial averaging loses the directional information needed to locate peaks. Deposits must also move consistently as receivers or clusters move. Changing the total normalization cannot by itself repair an incorrect shape.

Next retain directional deposition information, propagate luminosity/distance and history uncertainty, and evolve a supported common source before testing motion and lensing. No final holdouts were opened; all six physical objectives remain open.

Run `python research_work/results/cluster-catalog-deposit-profile/run.py`.
