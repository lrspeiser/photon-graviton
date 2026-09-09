# Temporal-field hypothesis checks

These calculations explore fictional non-expanding physics against existing observations. They do not establish a mechanism or replace a complete cosmological model.

## Local reanalysis

checks.py uses the existing 164 Cosmicflows-4 SBF galaxy-group representatives: 104 train, 35 validation, 25 test. Each variant fits one slope on the training subset using the original Gaussian objective, propagated distance errors and 300 km/s velocity floor. The test labels have previously been examined: new comparisons are exploratory, not a fresh blind validation. Original input SHA256 hashes are in results.json. Calibration, selection, peculiar-motion and source-physics assumptions remain inherited from the previous experiment.

Put the included source_inputs directory alongside checks.py and run Python with numpy/scipy installed. No network needed.

Full-signal model: S=exp(k*r/c), D_L=r*S, z=S-1. Energy-only model: same energy shift but no envelope stretching, D_L=r*sqrt(S). General inversion with a=1 or 1/2 is z=exp(W(a*k*D_L/c)/a)-1.

Positive-history alternative: n(t)=exp(gamma*(t-t0)), fixed Euclidean matter positions, unchanged atomic reference, n(t0)=1 and v=c0/n. Integrating the ray gives r=c0*(S-1)/gamma, hence z=k*r/c and D_L=r*(1+z). Its inversion is z=2*x/(1+sqrt(1+4*x)), x=k*D_L/c. This avoids n=0 at finite time but does not establish screened boundaries, energy backreaction, causal consistency or atomic coupling. It is a different hypothesis from the previously proposed bounded sigmoid. The full sigmoid is not fitted: its potential and transition history are unspecified, and nearby data cannot justify extra parameters.

## External summary-measurement checks

DES (1504 supernovae): b=1.003 +/-0.005 statistical +/-0.010 systematic for duration stretching (1+z)^b. Source https://arxiv.org/abs/2406.05050 and https://academic.oup.com/mnras/article/533/3/3365/7738388 . We combine quoted errors in quadrature for a rough consistency calculation, not a new raw light-curve analysis. Method two constructs a reference with assumed dilation and is explicitly partly a consistency check; their method one independently minimizes scatter. Do not claim a new 90-sigma or 200-sigma rejection from extrapolating this Gaussian summary. Source evolution and passband modeling remain assumptions.

Clock comparison: alpha_dot/alpha=(1.0 +/-1.1)e-18/year, https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.126.011102 . Define a conditional coupling q via alpha_dot/alpha=q*gamma, gamma=7.731496595524618e-11/year. An approximate conservative 95% absolute bound is |q| < (|mean|+1.96*sigma)/gamma =4.082e-8. This is NOT a general bound on time variation, nor on a common rescaling of clocks. In the original epsilon=n*epsilon0, mu=n*mu0 illustrative model, c=c0/n and epsilon*c is constant: one must not assume alpha varies as 1/n. Detailed spectral ratios and atomic coupling still need calculation.

Present CMB temperature: T0=2.72548 K (Fixsen https://arxiv.org/abs/0911.1955). Assuming the local usual blackbody energy density, u=a*T0^4=4.17466e-14 J/m3. If all this radiation loses energy at the PRESENT fractional rate gamma, Q=gamma*u=1.02278e-31 W/m3. Extending that fractional rate constantly for 1 Gyr in a fixed volume with no replenishment transfers 7.44017% of initial radiation energy to the field. This illustrative energy accounting does not solve field backreaction, thermalization, photon production or gravitational stress; the fractional rate need not stay constant in other histories.

Distant CMB: H2O absorption at z=6.34 implies 16.4--30.2 K (published 1-sigma range), https://arxiv.org/abs/2202.00693 . T0*(1+z)=20.005 K lies inside that range; an everywhere identical 2.725 K background does not. This tests a constant-temperature hypothesis, not every static geometry. In a changed atomic/radiative theory the excitation inference must be recalculated. No new raw spectral fit was performed.

## Untested or underidentified hypotheses

The 164-group features contain distances, their errors, sky position, group identifiers and split/tile membership. They do not contain independent line-of-sight density maps. No void-screening or gravity-cliff test is claimed. A constant active path fraction f is exactly degenerate with k in k*f*r. Spatially variable f requires independently reconstructed path environments before fitting.

A global distance-scale change is exactly degenerate with k in the original k*D_L prediction. The earlier 8.744% rate shift is equivalent in magnitude to 0.1820 mag of distance zero point. It is not evidence for environmental variation by itself.

Rotation and lensing cannot be newly predicted from these photon hypotheses until a spatial field equation, matter coupling, and photon trajectories are specified. The reservoir energy is not directly measured by a redshift fit. No claim of discovery or elimination of all dark-matter alternatives is made.
