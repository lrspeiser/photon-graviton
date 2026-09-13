# Added-noise test of the object-excluded spectral clock

The tested independent Gaussian noise increases phase errors but does not
generate the roughly 25-percent median aging slowdown required by the earlier
group-specific no-stretch fit. This is a limited processed-library experiment,
not a constraint on every distant-observation systematic.

## Fixed experiment

The previous 269 target spectra and 2,210 comparison spectra are unchanged.
All templates belonging to the target object remain excluded. The method still
uses the median phase of the best matches from five distinct comparison objects.
The vectorized implementation reproduces both previous zero-added-noise outputs
exactly before noisy cases are scored.

Eight declared seeds generate unit Gaussian noise independently per processed
wavelength sample. Noise standard deviation is the target's centered spectral
RMS in the broad window divided by a selected contrast SNR of 3, 1 or 0.3.
This definition concerns flattened spectral features: it is not the photon SNR
or reported continuum SNR of an instrument. Existing template noise is already
present and has not been subtracted or assigned a covariance.

For each seed, the same random draw is scaled across the three noise levels and
reused in the overlapping wavelength windows. These paired cases are not
independent experiments. The noiseless comparison library and the eligible
object/phase sample remain fixed; no rejection threshold is tuned after adding
noise. Gaussian noise and normalized-correlation matching are established
methods, used here as declared development choices rather than a new physics law.

## Results

| Contrast SNR | Rest wavelength window, Angstroms | Mean phase RMS over seeds, days | Mean of seed median object slopes | Range of seed median slopes |
|---|---|---:|---:|---|
| 3 | 4000–7000 | 1.962 | 1.021 | 1.011–1.038 |
| 3 | 4000–5500 | 2.165 | 1.027 | 1.013–1.050 |
| 1 | 4000–7000 | 2.224 | 1.018 | 1.001–1.033 |
| 1 | 4000–5500 | 3.187 | 1.026 | 0.998–1.045 |
| 0.3 | 4000–7000 | 4.566 | 1.027 | 0.980–1.065 |
| 0.3 | 4000–5500 | 5.880 | 1.074 | 1.013–1.163 |

For each seed and window, the aging slope is estimated separately in 17 objects
with at least three epochs spanning 15 labeled days. The table first takes the
median across those objects, then averages that median over eight seeds. Phase
RMS instead weights spectra. Every object's slope and every recovered phase are
retained, including failures to recover the correct phase accurately. Five
objects remain in the phase-error calculation but fail the declared slope
span/count requirement.

These seed ranges are not confidence intervals or uncertainty coverage. A
median near one does not imply every supernova is accurately calibrated.
At the strongest tested noise, some individual phase estimates are poor even
though the median slope remains near or above one.

## What this resolves and what it does not

Neither this noise prescription nor the earlier wavelength restriction supplies
the required 25-percent population clock shift in this estimator and library.
We should therefore not invoke those tested effects as an explanation of the
observed distant aging-rate pattern. The result does not prove the absence of
source evolution or all calibration bias.

The experiment does not reproduce observed flux uncertainties, sky subtraction,
spectral resolution, host-galaxy contamination, source selection or correlation
quality cuts. It does not shift and resample actual observer-frame spectra.
The assigned phases still come from the existing library and share its
calibration. Their agreement under added noise is not an independent clock
measurement or an implementation of the full SNID algorithm.

The next empirical step requires actual distant spectra and observing metadata,
with errors propagated through the same object-excluded analysis. The separate
160-case artificial light-curve calibration is unchanged and remains a different
test. No reserved scientific observations or prediction scores were opened.

Run `run.py` to reproduce the 48 seed/noise/window cases. `protocol.json` records
the choices; `results.json` preserves predictions, per-object summaries, zero-noise
checks and dependency hashes. The loader checks the cached archive hash before
reading its spectra.
