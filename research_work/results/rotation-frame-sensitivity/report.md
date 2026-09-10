# Reference-frame sensitivity of the two Milky Way rotation curves

**The different adopted solar frames can explain part of the offset between the published red-giant and Cepheid curves in a simple approximation, but not all of it.** The RMS difference between those two observed-summary curves falls from 12.95 to 8.68 km/s in the central-ray, unchanged-drift approximation below. All twelve Cepheid comparisons remain higher than the interpolated red-giant curve.

These numbers are **curve-to-curve differences**, not revised errors for our gravity model. The frozen Cepheid model scores remain 20.46 km/s for the primary completion and 14.82 km/s for its previously calibrated mass variant. Neither published catalog nor any earlier score is changed.

## Published reference frames

| Source | Adopted solar Galactic radius | Adopted solar azimuthal speed |
|---|---:|---:|
| [Eilers et al., section III](https://arxiv.org/html/1810.09466v2) | 8.122 kpc | 245.8 km/s |
| [Feng et al., section 2.2](https://doi.org/10.1093/mnras/stag011) | 8.275 kpc | 250.2 km/s |

The differences are 0.153 kpc and 4.4 km/s. Stellar motions are measured from the Solar System; converting them to a Galactic frame requires an adopted position and velocity for the Sun. A different choice can therefore shift an inferred rotation curve without changing the underlying gravity.

## What was calculated

**Known coordinate transformation specialized to a conditional geometry; not a new physical law:** on the Sun–Galactic-center line, with the star's heliocentric position and velocity held fixed and both positions on the same side of the Galactic center,

```
R_Eilers_frame = R_Feng_frame - 0.153 kpc
vphi_Eilers_frame = vphi_Feng_frame - 4.4 km/s.
```

The radial coordinate identity is checked by reconstructing the same heliocentric position in the two frames. Away from that central line, at fixed cylindrical direction axes, the azimuthal velocity contribution of this tangential solar-speed change is 4.4 cos(phi) km/s. Across |phi|<=30 degrees that contribution lies between 3.81 and 4.4 km/s. **This is not a bound on the complete change in the circular-speed estimator.** A shifted Galactic center also changes the direction axes, radial moments, bin assignments and density gradients.

The published quantity is circular speed Vc, not an individual star's vphi. Applying the same velocity offset to Vc additionally assumes an unchanged additive asymmetric-drift correction. We use it only as an illustrative central-ray sensitivity. It is not an exact re-reduction of either catalog, and it does not capture distance calibration or non-axisymmetric stellar motion.

**Known interpolation and summary statistics:** the archived Eilers bins are interpolated with a shape-preserving cubic interpolator at the twelve Cepheid radii. All queries stay inside the archived radius range. We compare the published curves directly, then compare after the stated radius/velocity transformation. No offset is fitted to the differences: both shifts come from the published frame constants.

| Curve-to-curve diagnostic | RMS difference (km/s) | Mean Cepheid minus Eilers (km/s) | Positive differences |
|---|---:|---:|---:|
| Published radii/speeds, no frame adjustment | 12.95 | 12.03 | 12/12 |
| Central-ray, constant-drift approximation | 8.68 | 7.24 | 12/12 |

Replacing cubic interpolation with piecewise linear interpolation changes individual unadjusted differences by at most 0.357 km/s, and adjusted differences by at most 0.483 km/s. This tests interpolation choice, not measurement covariance. No significance or complete uncertainty bound is inferred.

## Consequence for the theory test

The frame mismatch should not be ignored, but this approximation does not remove the population-to-population difference. It also does not establish that the remaining difference must arise from gravity: distances, selection, streaming motions, asymmetric drift and other modeling choices remain possible contributors.

Before a rigorous combined inference, the individual stars must be transformed in one declared frame, with uncertainties propagated, and each population's velocity distribution modeled under the same potential. That is different from shifting an already-compressed circular-speed table. The current check keeps this distinction explicit and does not quietly replace the frozen new-catalog result with a more favorable post-exposure score.

This is a post-exposure diagnostic. The Cepheid table remains exposed for future development, and the original protocol and scores remain authoritative. There is no new photon-conversion, capture, lensing or field-response formula in this pass.

`run.py` reproduces `results.json` and all twelve rows in `curve-differences.json`, with input hashes. No external data acquisition or parameter optimization is performed by that script.
