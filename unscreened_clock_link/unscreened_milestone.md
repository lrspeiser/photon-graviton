# Unscreened time theory: clock response and a measurable radio link

The next milestone is to specify how the proposed temporal field affects the clock that generates a signal, the signal during its journey, and the clock that measures it. Removing density screening does not settle those relationships. The equations below define an unscreened branch and derive a two-way observable before attempting a fit to spacecraft observations.

## 1. Working assumptions

The primary calculation assumes fixed material-coordinate endpoints separated by R, unchanged ordinary material rulers, a nondispersive optical propagation speed c/n(t), and a coherent transponder that immediately retransmits the received phase with a fixed frequency ratio q. Ordinary gravitational and kinematic effects are omitted in this benchmark and must be included in the actual tracking model. Gravitational waves are not fitted here; retaining their shared propagation with light remains a separate required completion.

Use n(t)=exp(γt), with γ=7.7314965955×10⁻¹¹ per year and n=1 at the reference epoch. This is the previously retained rolling-history branch. Unlike a distance-only fitted curve, it specifies how conditions evolve between emission and reception. Material atomic reference frequencies are parameterized as f_atom(t)=f0 n(t)^(−p). The exponent p is a placeholder for a physical clock coupling, not a newly measured constant. The main unscreened candidate takes p=0: material clock frequencies do not change with n.

Ruler and material-dynamics responses must eventually be derived too. If atomic lengths change, numerical changes of a dimensional light speed cannot be interpreted without comparing consistent material units. The fixed-ruler, fixed-endpoint assumptions make this calculation concrete; they do not establish a complete matter action or planetary dynamics.

## 2. The two-way observable

Let t0 be transmission, t1 the coherent retransmission, and t2 reception. Conservation of phase in the homogeneous optical model gives the returned coordinate frequency

f_return = q f_atom(t0) n(t0)/n(t2).

Comparing it with the Earth reference at reception gives

\[
\boxed{1+y\equiv\frac{f_{\rm return}}{qf_{\rm atom}(t_2)}
=\left[\frac{n(t_2)}{n(t_0)}\right]^{p-1}.}
\]

Here y is a fractional frequency residual; y<0 is a redshift. The coherent turnaround ratio q cancels. An ideal coherent transponder therefore does not require comparing an onboard free-running clock with an Earth clock. Actual DSN observables also require uplink frequency ramps, finite count intervals, station clock conversions and known transponder behavior.[1,2]

For fixed endpoints, the integrated round-trip path satisfies ∫c/n dt=2R. With rolling n=exp(γt) and n(t2)=1,

\[
S\equiv n(t_2)/n(t_0)=1+2\gamma R/c,
\qquad y=S^{p-1}-1.
\]

For the retained affine history n=1+γt, the corresponding S is exp(2γR/c). Both histories predict the same leading solar-system term,

\[
\boxed{y\simeq-2(1-p)\gamma R/c.}
\]

This is a derived two-way prediction under stated assumptions, rather than automatically doubling a one-way fitted shift. At solar-system distances the difference between these histories is higher order and extremely small, so this experiment primarily tests the effective coupling, not which cosmological curvature of the distance law is correct.

## 3. Three retained clock-response cases

| Clock response | Measured consequence at fixed positive γ | Implication |
|---|---|---|
| p=0: unchanged atomic frequencies | Redshift and event stretching survive | Primary unscreened branch to test |
| p=1: f_atom proportional to 1/n | Frequency shift cancels exactly | Also cancels this mechanism's astronomical redshift under identical endpoint standards |
| p=2: stronger atomic response | Net measured blueshift | Retain as a separate coupling branch; changing γ or other physics changes its interpretation |

A universal p=1 cannot hide the effect only in the laboratory while keeping it in distant galaxies under the same assumptions. For p≠1, the leading observed redshift slope determines γ_eff=(1−p)γ. If γ is then refitted to preserve the astronomical slope, the leading solar-system signal depends on that same γ_eff. Merely adjusting a universal p therefore does not independently suppress the local signal. Higher-order history, material dynamics and different source physics are separate questions.

No choice of p has been derived from atomic physics here. A constitutive electromagnetic model may change atomic frequencies and lengths together, so assigning p=0 requires a physical completion. These cases expose the requirement rather than resolve it by naming a new kind of time.

## 4. Range and Doppler must agree

The station measures elapsed time in its own clock units. Define dτ=n^(−p)dt and inferred radar range ρ=c[τ(t2)−τ(t0)]/2. In the fixed-endpoint rolling model,

\[
\rho=\frac{c}{2\gamma}n(t_2)^{-p}\frac{S^p-1}{p}\quad(p\ne0),
\qquad \rho=\frac{c}{2\gamma}\ln S\quad(p=0).
\]

At p=1 this equals R exactly, matching frequency cancellation. At p=0 and to leading order,

\[
\boxed{\frac{d\rho}{d\tau_2}\simeq\gamma R,\qquad y=-\frac{2}{c}\frac{d\rho}{d\tau_2}.}
\]

The last identity is exact in the ideal fixed-endpoint model. It ties frequency and ranging predictions to one parameter rather than fitting unrelated anomalies.

| Fixed one-way distance | Two-way y, p=0 | Equivalent range rate | Apparent range change per year |
|---|---:|---:|---:|
| 1 AU | −2.4451×10⁻¹⁵ | 0.0003665 mm/s | 11.57 m |
| 10 AU | −2.4451×10⁻¹⁴ | 0.0036651 mm/s | 115.66 m |
| 30 AU | −7.3353×10⁻¹⁴ | 0.0109953 mm/s | 346.98 m |
| 150 AU | −3.6676×10⁻¹³ | 0.0549764 mm/s | 1,734.92 m |

These are calculated signals, not observed unexplained drifts. The annual figures assume fixed true separation, fixed rulers and unchanged material clocks. They do not predict that a moving spacecraft simply gains this amount of physical orbital distance every year. Baseline orbit parameters, clock conventions and changing geometry can absorb or modify an apparent drift.

At the reference epoch the inferred range offset itself is much smaller: about −0.0183 m at 10 AU in the rolling branch. Its derivative can nevertheless accumulate over a long observing interval because n keeps evolving. Reporting only an accumulated annual number without the normalization and evolving clock model would be misleading.

## 5. First calculation using Cassini geometry

A JPL Horizons query retrieved 13 reconstructed geocentric Cassini states at 30-day spacing from 2005-01-01 through 2005-12-27. The source identifies a final post-mission reconstructed trajectory. The geometric Earth–Cassini ranges span 8.0944–10.0858 AU. These are genuine JPL ephemeris products, which incorporate standard dynamics and data fitting; they are not raw Doppler measurements or observations independent of a dynamical model.[3]

Evaluating the fixed-endpoint formula separately at those ranges gives y from −1.9792×10⁻¹⁴ to −2.4661×10⁻¹⁴. Its equivalent velocity scale is 0.00297–0.00370 mm/s. After subtracting the mean of these 13 forecast values, their RMS is only 1.754×10⁻¹⁵. This subtraction illustrates that much of the raw forecast is nearly constant. It is not a full projection against orbital parameters or a fit to actual data.

Published Cassini experiments achieved approximately 3×10⁻¹⁵ fractional frequency stability at 1,000 seconds under favorable conditions.[4] Comparing a forecast with that noise scale motivates investigation, but does not establish detectability: a nearly constant offset may be absorbed by velocity or calibration parameters, and the observations have correlated errors. The monthly sampling here does not represent the actual tracking schedule or independent 1,000-second data points.

A moving-endpoint light-time solver must evaluate the spacecraft at turnaround and Earth at transmission and reception, not insert an instantaneous range into the static formula. If a propagation-induced range term depends on time times R(t), its derivative also contains a time times range-rate contribution. The CSV deliberately labels its entries stationary-endpoint forecasts so that this geometry exercise is not mistaken for a full model of moving Cassini.

## 6. Numerical verification completed

The reference calculation uses 90-digit decimal arithmetic because ordinary floating-point subtraction can lose a fractional signal around 10⁻¹⁴. An independent one-second pulse propagation calculation at 10 AU agrees with the differential frequency-stretch factor to approximately 3×10⁻³² relatively; the remaining difference is consistent with finite pulse spacing. The p=1 cancellation of both frequency and radar range was checked for both histories. A symmetric finite-difference derivative checks the range–Doppler relation with a residual around 4×10⁻⁴² m/s. These checks validate the implemented ideal equations, not their applicability to nature.

## 7. Empirical test protocol now defined

The next observational calculation should fit a single effective temporal slope jointly to coherent Doppler and range, with standard nuisance parameters. Cassini ODF archives contain Doppler, range and frequency-ramp information, and identify associated trajectory, Earth-orientation, ionosphere and troposphere products. The inspected SCC1 archive covers July 2005 solar-corona sessions; its strong plasma environment makes it a challenging diagnostic sample rather than an automatically clean first sample.[1]

The required sequence is:

1. Implement moving-endpoint light-time and station clock mappings, preserving both rolling and affine histories as branches. Include finite transponder delay and the actual turnaround/ramp conventions.
2. Reproduce the standard model's observables before adding the temporal term. Use consistent station positions, reference frames, relativistic delays and media corrections.
3. Inject a known temporal signal into a controlled residual simulation using real observing times. Verify that joint orbit and clock fitting can recover it; report how much is absorbed by nuisance parameters.
4. Fit the no-screening prediction with its rate fixed by the astronomical branch, then separately allow one fitted effective rate. Keep those tests distinct.
5. Reserve time arcs and, where possible, another mission for prediction. Avoid giving each arc an independent time coefficient that would erase the universality test.
6. Report agreement, discrepancy, or lack of identifiable sensitivity, retaining all physical branches and showing the effect of nuisance assumptions.

The criterion for this milestone is an operational prediction from clock through propagation to readout, supported by independent pulse and range checks. That part is completed for ideal stationary endpoints. The retrieved Cassini geometry provides a first forecast. No ODF observations have been decoded, no observed residuals have been fitted, and no local detection or exclusion is claimed.

## Sources and files

[1] NASA PDS. Cassini RSS Raw Data Set SCC1 V1.0. https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=CO-SS-RSS-1-SCC1-V1.0 . Archive overview, observables and ancillary-data requirements.

[2] NASA. Basics of Space Flight, Chapter 10: Telecommunications. https://science.nasa.gov/learn/basics-of-space-flight/chapter10-1/ . Coherent spacecraft transponder operation.

[3] JPL. Horizons API documentation. https://ssd-api.jpl.nasa.gov/doc/horizons.html . Exact request, returned states and source description are bundled in data/horizons_request.json and data/horizons_response.json. Target −82, center 500@399, geometric vectors, KM-S, 30-day sampling in TDB.

[4] Armstrong et al. Reducing antenna mechanical noise in precision spacecraft tracking (2008). https://doi.org/10.1029/2007RS003766 . Published Cassini Doppler stability benchmark; not a uniform noise model for this forecast.

Reproduce with `python calculate.py`, using Python's standard library only. Input geometry is bundled, so no network is required. `fetch_geometry.py` records the acquisition query and requires network only if refreshed. The outputs are static_link_results.json, cassini_geometry_forecast.csv and geometry_summary.json. The scientific scope is an unscreened branch with explicit provisional clock response, not a completed universal time theory.
