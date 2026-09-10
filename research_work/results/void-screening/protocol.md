# Void screening, local path lengths and clock direction: diagnostic protocol

Written before running this diagnostic, 9 September 2026. This is a conditional mathematical test of the latest proposed clock factor, not a fit to new observations or a completed gravitational theory. Energy-supply normalization and halo formation are deferred at the user's request.

## Questions

1. How small would unscreened accumulated photon loss be over Solar System paths, compared with million-light-year paths?
2. Does multiplying an ordinary clock lapse by a screened slower-void factor preserve its ordinary direction and local gradients?
3. Does redshift require clocks to slow, stop, or reverse? Distinguish environmental amplitude from temporal evolution and from the rate of photon loss.

## Fixed assumptions and choices

- Use the previously exposed constant conversion coefficient from `../conversion-first/results.json`, alpha=0.0002488993286382367/Mpc. Do not refit or claim a calibrated microscopic void rate. Applying it at local distances is a hypothetical extrapolation. For the screened examples treat it as alpha_void solely to hold the photon-loss scale fixed across comparisons.
- Standard constants: c=299792458 m/s; au=149597870700 m; Julian year=31557600 s; parsec=648000/pi au; nominal IAU GM_sun=1.3271244e20 m^3/s^2, GM_Earth=3.986004e14 m^3/s^2 and equatorial Earth radius=6378100 m. The mass parameters and radius are nominal conversions, not precision orbit data.
- Compute unscreened z=expm1(alpha L), energy-loss fraction=-expm1(-alpha L), and c*alpha. Paths: 1 m, 20200 km, 1 au, 100 au, 1 light-year, 1 million light-years, 100 million light-years. No visibility threshold or device sensitivity is inferred from length alone.
- Proposed environmental factor: S=1/[1+(W/W_star)^n], F=1+sigma*epsilon*S. Sigma=-1 is the recommendation being tested; sigma=+1 is an explicitly different faster-void control, not an adopted model. For either control use alpha=alpha_void*S; for the faster control this requires changing the old link to alpha=kappa*abs(F-1), not retaining a negative conversion coefficient.
- In the weak field, use q_base=sqrt(1-2W/c^2) as a static lapse surrogate, with q_total=q_base*F. It has the correct leading potential dependence; for combined Earth/Sun/background potentials it is not an exact multi-body GR solution. Compare stationary hypothetical clocks, omitting real orbit velocities, tides, rotation and signal-transfer corrections. Use differences and ratios; a common rescaling of coordinate time has no observable significance.
- Define W_Earth=GM_Earth/R_Earth+GM_sun/au+B. Define W_GPSheight=GM_Earth/(R_Earth+20200km)+GM_sun/au+B. For a 50-au illustrative location use GM_sun/(50au)+GM_Earth/(49au)+B. Background B takes 0 and 3e10 m^2/s^2. The second is an illustrative external well depth, NOT a measured Milky Way potential or a dark-matter map. The first is a sensitivity control, not a claim the Milky Way has no potential.
- Scan epsilon in {1e-12,1e-9,1e-6,0.01,0.1}; W_star in {1e6,1e8,1e10} m^2/s^2; n in {1,2,4}; both signs and both backgrounds: 180 fixed cases. All scan parameters are illustrations, not fitted values. Report all cases; do not select a successful theory from them.
- For each case report local suppression, extra Earth-to-GPS-height and Earth-to-50-au clock-ratio changes, and the derivative contribution relative to the ordinary baseline. Compare the derivative also at W=W_star, to expose transition forces that local screening can conceal.
- Under the earlier universal lapse coupling with locally Euclidean spatial metric, the extra acceleration of released slow matter relative to static observers is a_extra=-c^2 grad ln F. For an outward-decreasing W, sigma=-1 contributes outward acceleration; sigma=+1 contributes inward acceleration. This is a conditional coupling diagnostic, not a constraint on all modified-gravity theories. The extra-to-baseline acceleration ratio equals the ratio of their log-lapse derivatives with respect to W, since the common gradient and minus sign cancel. Absolute forces for real galaxies require a real field and spatial metric.
- Include a uniform unscreened F control. It cancels from local clock ratios and has no spatial-gradient force in this idealized coupling even if F differs from 1. Conversion can remain nonzero by separate postulate. This prevents confusing a common coordinate factor with an observable local clock anomaly.

## Checks and limits

Verify path conversions using standard unit identities; check the analytic derivatives symbolically with SymPy and all endpoint ratios/derivatives with independent 60-digit mpmath calculations. Check the mode signs and positive lapse over the declared grid, epsilon=0 recovery and exact uniform-factor cancellation. Use stable log1p/expm1 arithmetic for tiny differences. Relative errors in nonzero extra clock ratios and derivatives must be below 1e-7; path formulas below 1e-12. No astronomical likelihood, confidence interval, validated screening threshold or claim of compatibility with all Solar System tests follows from these checks.

Relevant existing evidence: [universal-clock force and timing diagnostic](../universal-clock-coupling/report.md), [spatial packet-transfer timing result](../transport/spatial-transfer.md). A stationary loss law does not by itself stretch event arrival intervals; a static clock contrast does not by itself produce a lasting matched-endpoint shift.

External references: [IAU nominal constants](https://arxiv.org/abs/1510.07674), [NIST GPS relativity](https://www.nist.gov/atomic-clocks/a-powerful-tool-for-science/putting-einstein-test), [Galileo gravitational-redshift test](https://arxiv.org/abs/1812.03711), [environmental screening precedent](https://arxiv.org/abs/1001.4525). The Galileo fractional-deviation constraint applies to that experiment's varying relativistic redshift signal; it is not a universal bound on F-1 or on the one-way photon-loss fraction.
