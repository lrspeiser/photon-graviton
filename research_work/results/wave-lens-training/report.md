# First observed comparison of supported wave sources

With one shared field mass and one shared deposited-mass fraction, the supported wave model modestly reduces training RMS errors relative to the same photometric ordinary-matter baseline, but **does not give a satisfactory joint match**. It changes the typical prediction from too little gravity to too much: stellar aperture dispersions are median 11.1% high and lens Einstein angles median 28.7% high.

This is a conditional 32-system training comparison, not validation of the wave model or evidence that photons produced the source. No held-out outcomes were newly opened. The earlier negative validation of the empirical profile remains unchanged.

## What was specified before this run

The archived protocol fixes a shared particle mass m=10^-24 eV/c^2 and total source-to-stellar mass ratio f=1. The field mass is an illustrative numerical choice that keeps dimensionless wave parameters near the previously checked family; it was not estimated from this run's velocity or lens residuals. The mass fraction is also illustrative, not derived from the photon supply. Neither is optimized here.

For each training lens, we use the archived Salpeter stellar-population estimate rescaled under the energy-loss-plus-event-stretch propagation branch, and the updated I-band effective radius. The same stellar mass and radius enter both the ordinary-matter benchmark and the wave model. Neither is adjusted to fit stellar velocity dispersion or lens angle. Of the 33 earlier eligible training lenses, J0109+1500 lacks the required population mass, leaving 32 matched systems. The eligibility restriction follows missing input data, not residuals.

The population masses remain **conditional estimates**. Their rescaling changes luminosity normalization while retaining published stellar population, age, dust, metallicity, stellar evolution and IMF assumptions; it is not a full population-posterior refit in our fictional universe. This is a more independent check than obtaining mass from the very velocity being tested, but it is not assumption-free mass measurement.

## Prediction pipeline

The source is solved separately in each ordinary galaxy with the same physical field mass and f. Its dimensionless parameter varies correctly as

\[
\eta=\frac{\hbar^2}{m^2GM_ba},\qquad a=R_e/1.8153.
\]

Across these systems eta ranges approximately 0.134 to 3.81. Holding eta fixed across galaxies would have silently changed the particle mass. The free stationary Schrodinger-Poisson equations, regular boundary conditions and branch-following solver are given in the [source-equilibrium report](../self-consistent-wave/report.md). These are **known equations applied to a proposed companion identity**, not a new fundamental interaction.

The total spherical acceleration is

\[
g(r)=\frac{GM_b}{(r+a)^2}+\frac{GM_c(<r)}{r^2}.
\]

The second term is omitted in the ordinary benchmark. The wave density is obtained from the equations rather than prescribed by the old empirical power law. A Hernquist luminous tracer in each potential supplies the known isotropic Jeans prediction

\[
j(r)\sigma_r^2(r)=\int_r^\infty j(s)g(s)\,ds.
\]

The predicted line-of-sight variance is projected and averaged through a 3-arcsecond-diameter aperture with Gaussian seeing FWHM 1.5 arcseconds. This seeing remains a shared sensitivity assumption, not an independently measured value for every exposure. Uniform stellar mass-to-light ratio, spherical geometry, isotropic stellar orbits and omission of separate gas and central black-hole components remain limitations. No velocity is predicted for a particular star solely from its location.

At leading weak-field order, the nonrelativistic wave source contributes the usual metric gravity to stars and photons. Integrating its enclosed mass and the ordinary density along the ray gives the deflection. Combining it with the archived static Euclidean redshift-derived distances gives a predicted Einstein angle. The target is the published intermediate-axis SIE image-model summary; it is not a raw image likelihood. The same density sets both motion and lensing, with no independent light-bending multiplier.

The pilot's distances are still conditional on our redshift rule, extrapolated beyond its nearby calibration. This run does not validate that geometry, event stretching or the propagation mechanism.

## Observed versus predicted summaries

All values below refer to the same 32 training systems and fixed photometric stellar masses.

| Model | Dispersion RMS error | Median predicted/measured dispersion | Einstein-angle RMS error | Median predicted/SIE angle |
|---|---:|---:|---:|---:|
| Ordinary matter | 62.77 km/s | 0.790 | 0.4628 arcsec | 0.663 |
| Ordinary matter plus supported wave | 52.50 km/s | 1.111 | 0.4400 arcsec | 1.287 |

Mean residuals change from -54.44 to +27.42 km/s for dispersion and from -0.4085 to +0.3777 arcseconds for angle. Thus a somewhat smaller RMS is not a claim of a good fit: substantial mismatches remain, and the shared wave source generally overshoots the photometric baseline's shortfall.

The full table in `predictions.json` contains measured and predicted quantities for each object. Quoted spectroscopic uncertainties are preserved, but these RMS summaries do not combine all stellar-mass, distance, orbital, seeing and imaging uncertainties. No reduced chi-square, statistical rejection probability or likelihood-ratio significance is claimed.

These RMS values should **not** be compared directly with the earlier velocity-normalized lens pilot as if only the gravity law changed. That earlier calculation chose stellar mass to reproduce aperture dispersion, used 33 systems and did not test the dispersion independently. The current comparison fixes photometric masses and uses 32 systems; its paired benchmark is the appropriate comparison.

## Checks and scope

All 32 source equilibria pass normalization and virial checks, with maximum relative errors below 1.8e-9 and 2.9e-9 respectively. The largest computed central potential depth divided by c^2 is 1.86e-5, supporting the weak-field/nonrelativistic use of these configurations. This does not establish their collective stability or how they are populated.

Doubling radial and aperture angular quadrature resolution for the first system and the two extremes in eta changes predicted dispersion by at most 0.000026 km/s. The ray integral is separately checked at each predicted Einstein angle against an independent projected-density calculation. Resolution of the aperture integration does not itself test the ray integral.

An additional cross-check uses the exact square-root mass scaling of the fixed-shape ordinary Jeans problem: these ordinary-matter dispersions agree with values reconstructed from the archived dispersion-normalized pilot to within 4.4e-11 km/s. This checks normalization and the aperture calculation, not the adopted population masses.

Numerical precision is therefore much smaller than the quoted discrepancies. The discrepancies belong to the chosen physical and observational assumptions, rather than being removable through finer aperture integration.

## Next scientific decision

This supplies an actual data comparison for a source with an explicit support equation. It does not yet supply a shared photon-conversion and capture theory. A carefully constrained training study could vary a common field mass and a common source-normalization rule, while propagating stellar-population and orbital uncertainty. Any candidate selected that way must then be evaluated on reserved data under a frozen protocol. The present shared f=1 must not be replaced by independent per-galaxy values merely to reproduce each target.

Even a successful future galaxy fit would not determine why companion energy forms this state or how much accumulates. Redshift, supernova event timing, broadband light behavior, lossless companion travel, capture, conservation, stability and photon supply remain parts of the full goal. The deferred total photon budget is not passed by treating f as an input.

Run `python research_work/results/wave-lens-training/run.py`, then `run.py --verify-resolution` in this directory and `verify.py`. Inputs and protocol are hashed in `results.json`. The source solver was made importable without triggering its nine-example run; its equations and example parameters are unchanged. Dependencies: NumPy and SciPy.
