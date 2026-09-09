# Observer-time evidence for testing the conversion model

## Why this is needed

Changing photon energy and stretching an event's arrival-time pattern are different predictions. The current stationary conversion kernels establish no common event-duration stretch. We therefore need the actual arrival-time photometry, not only widths produced by another analysis, to test this part of the fictional universe independently of an assumed expansion history.

This pass acquires and audits inputs. It does not estimate a time-dilation exponent, select a final supernova sample, or claim observational success. The previous goal turn made progress by deriving inverse companion transfers; this turn advances the separate timing requirement that remains mandatory.

## What the earlier timing evidence does and does not say

The recovered DES width products were built against time-scaled reference curves. Refitting those products alone is a consistency calculation, as the existing `time_revision/README.md` states. This does not make the entire published timing result circular: the authors also varied the temporal scaling exponent directly in stacked photometry, separately from the width fit. Their analysis matches emitted wavelength ranges and uses SALT estimates for peak normalization and timing. Their population and selection assumptions remain relevant. See [White et al., sections 4.1–4.3](https://arxiv.org/html/2406.05050v2).

Keep the observations and methodological checks. Do not import their expanding-universe interpretation into this project, and do not dismiss the measurements merely because the publication uses that interpretation.

## Inputs now available locally

The pinned [DES public release](https://github.com/des-science/DES-SN5YR/tree/c9a4fcafc4cbd19bd750dee47fc76194a45c181f) supplies calibrated observer-time photometry. This is not raw detector imagery. Instrumental corrections and error models remain part of its provenance.

| Product | What it supplies | Use in this pass |
| --- | --- | --- |
| DES HEAD FITS | Event identifiers, observation pointers, type flags, redshift and metadata | Validate indexing and inventory fields |
| DES PHOT FITS | MJD, band, calibrated flux, error and measurement flags | Validate measurement rows without time rescaling |
| Classification CSV | Published classifier probabilities | Inventory schema only; no probability-based selection |
| Release README files | Column definitions and calibration/update notes | Interpret units, pointers and data limitations |
| Authors' Methods.py and license | Exact reference-curve and input dependencies at a pinned commit | Inspect source; do not execute it |

All eight files were downloaded and verified against their Git blob hashes at pinned commits; SHA-256 hashes are in `input-audit.json`. About 70 MB of source files are kept in the ignored `research_work/generated/timing-inputs` cache, rather than duplicated in Git. The retrieval script reconstructs and checks that cache.

The FITS inventory contains **19,706 transient entries and 1,779,030 pointed-to photometry rows** in g, r, i and z. It includes more than the timing study's selected supernovae. All pointers are in bounds, nonoverlapping, and agree with NOBS; event identifiers are unique. Pointer-selected times, fluxes and errors are finite, with positive errors. The table also contains one separator row per transient, excluded by the pointers.

The release labels **353 events as spectroscopically confirmed SNe Ia**, all with positive finite heliocentric redshifts, spanning 0.0176 to 0.85. This is an inventory, not a selected or completeness-corrected timing sample. Using that subset would avoid relying on the photometric classifier for type identity, but would introduce its own spectroscopic-selection limitations and a smaller redshift range.

The acquisition applied no time dilation correction, flux normalization, width fit, redshift cut or quality optimization. No distance, expansion rate, cosmological age, dark-matter estimate or luminosity-distance simulation was used to generate these counts.

## Dependencies that remain before a valid timing inference

The downloaded `Methods.py` reads `FITOPT000.FITRES`, `FITOPT000.LCPLOT` and `hubble_diagram_wPIa.txt`. It uses zHEL, fitted peak dates, model peak fluxes and a probability-selected event list. Its `get_reference_curve` supports varying the time exponent and also defaults to scaling by 1+z. The acquired FITS files are not those exact fit-product files. Consequently, this acquisition is not a claim to reproduce the original paper's complete pipeline.

The pinned data release includes updates later than the 2024 paper. Its README documents calibration/metadata changes; its observed band labels must be read from the FITS table rather than assumed from an update note. Exact historical-paper reproduction would require matching historical inputs. A new observer-time reanalysis can instead use this explicitly versioned release and state the differences.

Before fitting, freeze a separate inference protocol with:

1. **Event selection and redshift provenance.** Decide and record whether the primary sample is spectroscopic or classifier-selected. A field named heliocentric redshift can still contain a photometric estimate when a spectroscopic value is unavailable; audit source flags. Preserve the user's fixed-distance rule without needing distance for this timing comparison.
2. **Measured-wavelength matching.** Compare emission bands using observed redshifts and actual filter throughput. Matching wavelength with 1+z does not by itself impose an event-duration factor. Central-wavelength shortcuts must be tested for bias.
3. **Observer-time normalization.** Fit peak times and flux scales without pre-dividing time by 1+z, or explicitly quantify the dependency introduced by reusing template-derived estimates. Photometric errors and selection cuts need documented handling.
4. **A variable timing law.** Fit or scan b in duration proportional to (1+z)^b, including b=0 and b=1. Rebuild the relevant reference times for every tested b, rather than holding a b=1-derived reference fixed and calling the result independent.
5. **Injection and null tests.** Inject both unstretched and stretched synthetic events into the actual cadence/error pattern, with matching selection and normalization. Recover both before interpreting the real-data result. These are pipeline checks, not synthetic confirmation of the theory.
6. **Correlations and source evolution.** Shared references and multiple bands are not independent events. Resample at event level and propagate reference construction. Predeclare sensitivity to intrinsic duration evolution, wavelength dependence and selection.

An exact degeneracy must remain visible: if mean intrinsic duration evolves as (1+z)^e while propagation contributes (1+z)^b, widths alone can measure b+e. Additional source information or a justified restricted evolution model is needed to separate them. Freely assigning a duration trend to match each observation would remove predictive content.

The DES events and aggregate results have already been examined in this project. Newly downloading their photometry does not make them an unexposed validation set. This work will be exploratory; genuinely independent validation remains a later requirement.

## Reproduce the acquisition

From the repository root:

```sh
python -m pip install -r research_work/results/timing-foundation/requirements.txt
python research_work/results/timing-foundation/acquire_and_audit.py
```

The script requires network access to verify the pinned Git trees and downloads missing cache files. It reads FITS data through Astropy and does not execute downloaded Python or unpickle external objects. Existing cached files must match the pinned blob size/hash. It writes only the timing cache and this pass's generated input audit; it does not overwrite the recovered historical products.

This is an acquisition command, separate from the offline physics diagnostic suite. The last numerical suite remains the 31-job companion-bath verification; no new timing fit has passed an observational test.
