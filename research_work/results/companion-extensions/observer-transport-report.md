# Observer geometry: angular spread and geometric delay

13 September 2026. Paraxial moments for prescribed independent kicks.

## Result

The previous accumulated-direction diagnostic is not an observed image width. Source-to-observer geometry weights deflections by their location. For the earlier 30.660139 Mpc path, the fixed-angle example gives an apparent directional RMS of about 0.520 arcseconds and a mean extra travel time of 5025 seconds. A source-size-tracking kernel gives much smaller values, but that kernel is prescribed rather than derived. Neither branch is an observational fit or a complete radiation-transport solution.

## Derivation and provenance

Place a source at s=0 and observer at s=D. Let a small transverse kick dtheta occur at s. For a ray to end at the observer, its initial slope changes by -(D-s)dtheta/D. Its final angular change is therefore (s/D)dtheta. Integrating the two straight segments' squared slopes gives the geometric excess path length dtheta^2 s(D-s)/(2D).

For independent zero-mean isotropic kicks with total angular variance rate q(s), these standard small-angle geometric relations give

\[
\left\langle\theta_{obs}^2\right\rangle=\int_{s_0}^D q(s)(s/D)^2\,ds,
\]
\[
\left\langle\Delta t\right\rangle=\frac1{2c}\int_{s_0}^D q(s)\frac{s(D-s)}D\,ds.
\]

Here q includes both transverse components; there is no additional factor of two. The mathematical geometry is known, while our choice of kicks and energy-loss process is hypothetical. The endpoint construction assumes the source supplies the needed nearby launch directions and no additional angle-dependent weighting of rates; it is not a full radiative-transfer boundary-value solution.

For a fixed positive energy step delta and prescribed mean fractional loss alpha,

\[
E(s)=E_0e^{-\alpha(s-s_0)},\qquad
\nu(s)=\alpha E(s)/\delta,\qquad q(s)=\nu(s)\theta_{step}(s)^2.
\]

This fixes the mean event-rate law to the existing redshift prescription; it does not derive alpha. Reverse transfers and their extra angular traffic are omitted, favorably to the model. The exponential energy law is known fractional-loss mathematics, not a novel formula. All sources are treated as static in this diagnostic.

## Two deliberately distinct kernels

1. **Fixed angle:** theta_step=d R_sun/(1 pc), with d=0.157241775 imported from the earlier dilute r=0.9 angular construction. This freezes the angle locally but does not maintain the same empty-output fraction along the ray.
2. **Apparent-source tracking:** theta_step(s)=d R_sun/s. This maintains the earlier geometrical scaling for an isolated source, but assumes the interaction can respond to apparent source size. We have not supplied a local physical cause for that dependence.

Both use E_0=2 eV, delta=1e-8 eV and the earlier alpha=0.0002488993286 per Mpc. Source radius is 6.957e8 m. The three start distances are illustrative domain choices (1 AU, 1 pc, 1000 pc), not proposed physical screening thresholds. The calculation does not model transfer before the chosen start.

For the source-size-tracking kernel, cancellation of the s^2 factors yields the independent check

\[
\left\langle\theta_{obs}^2\right\rangle
=(dR_{sun}/D)^2\frac{E_0}{\delta}
[1-e^{-\alpha(D-s_0)}].
\]

## Numerical comparison

At D=30.660139 Mpc and s_0=1 pc:

| Prescribed kernel | Apparent directional RMS | Mean geometric delay |
|---|---:|---:|
| Fixed angle | 0.520086 arcseconds | 5025.43 seconds |
| Tracks source size | 2.94087e-8 arcseconds | 5.22728e-10 seconds |

The unweighted fixed-angle direction proxy was about 0.902 arcseconds. Its reduction here is consistent with source-observer weighting, roughly a factor sqrt(3) for a nearly constant variance rate over the path. The adaptive kernel's much smaller image moment cannot be credited to an established interaction.

Starting the tracking kernel at 1 AU instead gives a mean delay of 9.16735e-10 seconds; starting at 1000 pc gives 3.00311e-10 seconds. The geometric delay is more sensitive than the apparent angular moment to near-source events. Where such events are rare, the mean need not describe a typical photon's delay; no Gaussian delay distribution is asserted.

## Consequences for the broader theory

A stationary scattering kernel applies the same delay distribution to successive otherwise identical pulses. It can shift and broaden an event, but does not by itself derive the common redshift-proportional dilation of its temporal structure. The timing problem therefore remains separate from these angular moments.

These are optical examples. They cannot be directly compared with gamma-ray/gravitational-wave arrival differences without deriving the wavelength dependence and the response of both messengers. Likewise, angular RMS alone does not establish an acceptable telescope image, surface brightness or flux selection.

Next, a local interaction must determine the angular kernel, with its total rate and recoil, and the coupled radiation field must establish whether the required output modes stay relatively empty. The fixed and source-size-tracking branches cannot be combined by silently taking the favorable feature of each. The traveling-companion and protected-storage mechanisms remain unfinished; the exact-third reference is unchanged.

## Verification

Eighteen cases cover three path lengths, three start distances and two kernels. Adaptive quadrature agrees with an independent 8193-point logarithmic-grid trapezoidal integration within 1e-4 relative tolerance for all three moments (observer angle, delay and raw direction). The source-size-tracking angular integral agrees with its analytic expression within 1e-10. Three explicit two-segment ray constructions verify endpoint closure and the single-kick delay factor. These tests validate the stated moment calculation, not a physical or observational success.

Executable: `observer-transport.py`; output: `observer-transport-results.json`. Predecessor: [angular outlet](angular-outlet-report.md).
