# Spectral clock template acquisition and identity audit

The official SNID archives are now cached locally and inventoried by content
hash. We can inspect the underlying template spectra instead of relying solely
on the published aging-rate table. This acquisition does not yet recover the
distant observed spectra or independently calibrate the template phase labels.

## Verified local contents

| Inspected set | Template files | Spectral columns | Objects shared with our aging-rate table |
|---|---:|---:|---:|
| templates-2.0 | 349 | 3,760 | 22 |
| snid-5.0/templates | 111 | 1,515 | 8 |

Counts include the set's non-Ia and non-supernova templates. They are not counts
of independent Type Ia supernovae or independent clock measurements. Every
parsed template has the declared wavelength/sample dimensions, finite values
and increasing wavelengths. Phase vectors have the declared number of entries.

All 22 overlaps in templates-2.0 are the nearby members of the existing 35-object
aging sample. There are no matched distant-sample identities in either inspected
set. Identity matching uses the old two-digit and newer four-digit year names;
it is not a full alias or underlying-exposure audit.

## Source and release distinction

The archives come from the [author's SNID page](https://people.lam.fr/blondin.stephane/software/snid/index.html).
The templates-2.0 README describes an expanded set incorporating later CfA
observations. The current software archive bundles several template directories,
including alternate wavelength grids; those variants are not additional
independent objects. No claim is made that either downloaded directory is the
exact frozen library used in the 2008 timing paper.

`manifest.json` preserves archive URLs, byte sizes, SHA256 hashes, and per-file
hashes. The downloaded sizes are 6,452,995 and 28,163,032 bytes. These are hashes
of the retrieved contents, not publisher-signed release checksums. The archives
remain under `research_work/generated/snid-template-audit/`, outside tracked
scientific results. No downloaded executable or installer was run.

## What the phase and spectrum fields do not establish

The [SNID manual](https://people.lam.fr/blondin.stephane/software/snid/howto.html)
describes template age selection, object exclusion and continuum processing.
Our parser records the phase and age-flag fields as data; it does not infer a
new physical clock from them. Template phases cannot simply be treated as
independent raw observer dates. Their time origin, frame corrections and source
metadata must be audited before testing a propagation law.

The template spectral arrays are processed matching inputs, not a complete set
of calibrated photon fluxes with observer dates, noise covariance, throughput
and extinction corrections. They cannot by themselves establish the joint
redshift, event-duration and absolute-brightness prediction required by goal 1.

## Necessary next test

When evaluating one nearby object, exclude every template from that object,
not merely its selected epoch. The 22 overlaps make this essential. Changing
wavelength coverage, degrading resolution/noise and shifting wavelengths can
then measure how the estimator's recovered phase depends on observing conditions.
Such transformations must preserve the known input phase and record that they
are artificial tests, not additional supernova observations.

To evaluate actual spectral aging, obtain the distant observed spectra and
their independent observation dates, plus provenance for nearby reference
dates. Fit wavelength shift without using an assumed event-stretch factor to
restrict allowed ages. Use object-level separation between template training,
method development and validation, and account for shared template errors.

This audit does not assert that the original authors made a self-matching error.
It identifies a condition our own reproduction must satisfy. No new fit or
reserved scientific prediction has been scored. The separate frozen synthetic
light-curve calibration continues unchanged.

Run `acquire.py` and `inspect_templates.py` to reproduce the inventory and array
checks. The acquisition script reads archive members in memory and does not
extract archive paths into the workspace. Future downloads must be compared
against the saved hashes before claiming identical inputs.
