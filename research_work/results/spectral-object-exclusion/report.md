# Object-excluded spectral age recovery and wavelength sensitivity

A simple independent matching implementation recovers the existing labels of
nearby processed spectra with about 1.8 days RMS error over 4000–7000 Angstroms,
worsening to 2.2 days over 4000–5500 Angstroms. All spectra from the tested
supernova are excluded from its comparison library. This establishes a usable
development baseline, not an independently calibrated supernova clock.

## Fixed method and sample

The inspected expanded SNID library supplies 2,210 eligible spectral columns
from 274 Type Ia objects after the declared phase and wavelength-support cuts.
The targets are 269 spectra of all 22 nearby objects in the exposed published
aging table. These are not the exact 145 spectra used in the historical paper.

Training phases range from -15 to 50 days; target phases from -10 to 30 days.
The same target spectra and comparison spectra are used in both wavelength
windows. Coverage is defined conservatively using the first and last nonzero
samples of the processed spectrum; it is not an instrument-derived validity
mask. Rejected support cases remain recorded in `results.json`.

For each target, all columns with its object identity are removed. The matcher
centers and normalizes each spectrum, then computes its normalized dot product
with comparison spectra at fixed rest wavelengths. For each comparison object
it retains the highest-scoring epoch; the median phase of the five best distinct
objects is the estimate. Thus a prolific reference object cannot occupy all five
positions. Every selected identity and score is recorded and checked for exclusion.

Normalized correlation and nearest-template phase estimation are established
methods. This particular configuration is a development choice with no novelty
claim. It is **not a SNID replication**: SNID's full processing, Fourier filtering,
quality statistic and uncertainty machinery are not implemented here. The input
library is itself processed and carries previously assigned phase labels.

## Results

| Rest wavelength range | Target spectra | Mean phase error | Median absolute error | RMS phase error |
|---|---:|---:|---:|---:|
| 4000–7000 Angstroms | 269 | +0.755 days | 1.10 days | 1.799 days |
| 4000–5500 Angstroms | 269 | +0.935 days | 1.30 days | 2.161 days |

These aggregate phase errors weight individual spectra, so objects with more
epochs contribute more. Individual-object summaries are retained rather than
treating the 269 spectra as independent experimental repetitions.

For 17 objects having at least three eligible epochs spanning at least 15 days,
an unweighted line fits recovered phase against labeled phase. The median slope
across objects is 1.0139 in the broad window and 1.0115 in the narrow window.
The median paired change in slope is -0.00390; individual changes differ. The
remaining five objects are not scored for slope because of the declared span or
count requirement, not because of unfavorable residuals.

An additive phase bias is different from an aging-rate bias: fitting an intercept
allows a constant age offset without forcing it to change the slope. These
slopes concern the library's phase coordinate, not independent observer dates.
They are not measurements of the propagation exponent b.

## Meaning for the timing question

This limited wavelength restriction does not produce a median 25-percent aging
slowdown in this nearby library. It therefore does not supply the group-dependent
clock change needed by the earlier no-stretch fit. This is not a general bound
on distant-sample bias: noise, spectral resolution, source-population differences,
host contamination and selection were not varied.

The recovery also does not validate the input phase labels. All targets and
comparison objects inherit the library's processing and calibration. Object
exclusion prevents direct self-matching but does not remove shared reductions,
calibration errors or source-population dependence.

The next useful calibration must preserve object exclusion while introducing
known noise/resolution and observer-frame wavelength transformations. Actual
distant spectra with dates are still required for an empirical aging measurement.
The photon interaction must eventually predict redshift, event stretch and
brightness together; accurate label recovery alone does not derive that physics.

`protocol.json` fixes the choices used here. Run `run.py` to reproduce the
calculation; archive and input hashes are checked/recorded. No original SNID
binary was executed and no reserved scientific observations were scored.
