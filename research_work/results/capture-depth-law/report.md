# Capture-depth dependence: radial forces and lensing from one rule

## Result

A common, steeper capture law can avoid the earlier extreme central deprivation and produce a nearly flat deposit-only circular-speed contribution over a finite radial interval. However, the favorable example declines farther out and has low total interception efficiency. It is a candidate for further testing, not a validated galaxy model or a demonstrated radiation supply.

The calculation compares all declared cases at matching central optical depth. It does not fit a separate density at each radius, use observed velocities as inputs, or certify a microscopic reason for the preferred-looking exponent.

## Postulate and comparison design

**Project postulate, not claimed unique:** in scale-radius units use

\[
\kappa(r)=A\left[\frac{2}{1+\sqrt{1+r^2}}\right]^p,\quad p=1,2,3.
\]

The bracket is a normalized saturating function of a fixed baryonic well's binding depth. A sets inverse-length opacity; p changes how strongly capture favors deeper regions. Choosing p=2 or3 is not a microscopic derivation. In particular, a squared depth dependence is not justified merely because some quantum probabilities involve squared amplitudes.

For each shape, set A by the specified diameter optical depth tau_d=0.1,1 or10:

\[
\tau_d=2\int_0^{R}\kappa(r)dr.
\]

This is a normalization for controlled comparison, not an astronomical calibration. Repeat at sphere radii10 and30 with unchanged scale radius. Isotropic incident companions, straight rays, absorption and in-place retention are the same explicit assumptions as the [preceding capture calculation](../isotropic-capture-profile/report.md). The resulting density follows from q=kappa integral I dOmega and rho_D proportional to q, not from an imposed halo profile.

**Conditional consequence of known spherical gravity:** when optically thin and far outside the core, kappa proportional to r^(-p) implies rho_D proportional to r^(-p). If p<3 and that outer contribution dominates enclosed mass, v_c is proportional to r^(1-p/2). Hence p=2 approaches a flat contribution asymptotically. That limit does not guarantee flatness over the finite core-to-boundary interval actually tested; attenuation also changes the profile.

## Finite-range predictions

The ratio v(9)/v(3) measures how much the deposit-only circular speed changes over a factor of three in radius. One means equal speeds. It is independent of the incident-intensity normalization.

| Depth exponent p | Thin: tau_d=0.1 | Moderate: tau_d=1 | Strong: tau_d=10 |
|---|---:|---:|---:|
| 1 | 2.058 | 2.212 | 5.647 |
| 2 | 1.451 | 1.563 | 3.246 |
| 3 | 1.076 | 1.144 | 1.920 |

These use outer radius10. Strong absorption continues to shield the interior and steepens the rising speed profile. At small opacity, p=3 gives only a7.6% speed increase across3..9, with fitted logarithmic speed slope0.052. It is the closest-to-flat member of this declared grid over that interval, not a universal or data-selected optimum.

The p=2 result still rises45% over3..9 because the finite profile has not reached its asymptotic regime. Looking only at the formal inverse-square limit would have overstated its success.

## Extension beyond the favorable interval

With the boundary moved to30, the p=3 thin case retains v(9)/v(3)=1.0762, so that middle-range behavior is not primarily caused by the original radius10 cutoff. But it then gives

\[
v(27)/v(9)=0.8275.
\]

That is a17.2% decline farther out. Its moderate-opacity counterpart declines14.8%; the strongly absorbing case is much flatter over9..27 but rises sharply over3..9. No grid case is shown to give a universally flat curve across all these radii.

These are additional-force contributions, not total velocities after adding a real disk, gas, bulge and their uncertainties. A declining contribution is not automatically excluded by real galaxies, but it cannot be advertised as an indefinitely flat halo prediction.

## Lensing and supply remain coupled constraints

The projected-to-spherical enclosed mass ratio at radius5 gives an intensity-independent diagnostic of lensing relative to local circular force. For p=3 and tau_d=0.1 it changes from1.313 at boundary10 to1.410 at boundary30. For p=1 at the same depth it changes from2.047 to3.563. Thus outer material can substantially affect projected lensing even when its effect on inner spherical forces is limited. The ratio is not itself an observed lensing mass measurement.

Interception is also different. The p=3 thin case captures0.948% of incident power at boundary10, and0.205% at boundary30. The incident power denominator grows with boundary area, so these fractions should not be mistaken for a loss of total absorbed power under fixed incident intensity. They do show why profile shape alone cannot establish energy sufficiency. We must calculate the actual external companion intensity, illumination duration and total capture power before claiming the required gravity is funded.

All intensities and gravitational normalizations remain unfitted. Retention/support, companion focusing, changing capture opacity under deposit gravity and nonspherical geometry are unresolved.

## Evidence and next use

[Protocol](protocol.md), [solver](run.py), [36 cases](cases.csv) and [results](results.json) are archived. Reproduce with:

```text
python research_work/results/capture-depth-law/run.py
```

Fine grids use641 radial nodes and96 nodes for angle/ray integration; coarse grids use321 and48. Volume absorption agrees with independent projected-chord absorption to relative1.58e-5 or better on fine grids. The central intensity agrees with the independent exp(-tau_d/2) identity. All declared speed-ratio/slope refinement and projected-mass quadrature gates pass. The largest coarse/fine speed-ratio difference is0.000725 and slope difference0.000317. These are numerical checks, not observational confidence intervals.

The p=2 andp=3 alternatives merit comparison with real ordinary-matter geometry and motion data, with one shared depth law and explicit scale/opacity parameters. Retain their full predicted curves and lensing consequences, rather than choosing a favorable radial segment. This experiment changes neither the astronomical redshift calibration nor any holdout exposure. It also does not resolve persistent redshift or its joint energy mechanism; all nine goals remain incomplete.
