# Data foundation for the fictional-universe theory

The roadmap now maps all 32 requirements to concrete observables or derivations, units, responsible tasks, local inputs and missing work. Twenty-one selected input products were checked for presence and hashed. This is a selected input audit, not a claim to possess every astronomical observation or a star-by-star catalog of the entire universe.

The previous goal turn made concrete progress by adding candidate-action and microphysical evidence. This step addresses the still-open T02 evidence protocol rather than selecting a physical law.

## What the galaxy distances actually mean

The local SPARC table defines distance method 1 as Hubble-flow distance using H0=73 km/s/Mpc and a Virgo-centric infall correction. Its method counts are:

| Distance method | Full 175-galaxy table | Reused 149-galaxy sample | 26 galaxies in energy comparison |
|---|---:|---:|---:|
| Hubble flow | 97 | 81 | 4 |
| Tip of the red giant branch | 45 | 37 | 11 |
| Cepheid magnitude–period relation | 3 | 3 | 3 |
| Ursa Major group distance | 28 | 26 | 8 |
| Supernova light curve | 2 | 2 | 0 |

These counts come directly from the recovered catalog and frozen object lists. Distance estimates from the other methods also have calibration and stellar-physics assumptions; “not Hubble flow” does not mean assumption-free. Group distances share correlated uncertainty.

The user has now fixed published galaxy distances as fictional facts. The Hubble-flow origin of some values does not invalidate their use or authorize expansion as an explanation. No distance values changed. See the governing [universe contract](../../../research_plan/universe-contract.md).

Historical sensitivity note, not an active permission to re-infer fixed distances: if a distance were changed, the same observations must be propagated through all dependent quantities. For fixed angular size, radius scales with distance; under a fixed luminosity-flux relation, luminosity scales with distance squared. Gas and stellar mass templates, inclination-corrected dynamics, source separations and external illumination also depend on the adopted calibrations. The inferred extra mass depends on both the total dynamics and the baryonic subtraction. A revised distance cannot be used only where it improves the photon supply.

The actual catalog rows use whitespace-separated fields whose layout differs from the byte offsets printed in its header. The audit parser therefore reads the 19 data columns and verifies the 175 names, the 149 unique selected names and the 26 unique matched names. A naive fixed-offset parse failed during development and was corrected before the reported results were produced.

## What is available, and at what level

| Local product | Evidence available | Interpretation that must be preserved |
|---|---|---|
| SPARC catalog and 175 rotation-model files | Galaxy properties, processed rotation curves, baryonic templates and published errors | These are processed observations/models. Physical radii, deprojected speeds and baryonic mass contributions carry assumptions. They are not original telescope images or spectral cubes. |
| DustPedia THEMIS and DL14 tables plus schema | Published SED-model outputs and two alternative reconstructions | Bolometric luminosities depend on distance and SED modeling. They are not directly measured photon-energy histories. |
| Archived 26-galaxy energy table | Exact benchmark quantities and source selections | Required mass and energy shortfalls are outputs of the earlier assumptions, not independent observations. |
| Cluster source-paper text and comparison table | Published mass/gas summaries and the previous derived comparison | Hydrostatic masses, gas fractions and R500 are inferred quantities; a new gravity theory needs profile-level observations and a consistent reconstruction. |
| Pantheon+ table and full covariance | 1,701 processed supernova entries and a 1,701 × 1,701 covariance | Entries are not necessarily unique supernovae. Corrected magnitudes, distance moduli and redshift-frame choices contain calibration/model processing. |
| DES event-averaged widths and spectral-aging summary | Previously processed timing information | Source templates, rest-band matching, selection and shared errors still need a forward model. |
| FIRAS spectrum table | Frequency, reconstructed monopole spectrum, residuals and published uncertainties | The header states that the supplied monopole is a 2.725 K blackbody plus residuals. It is not an untouched detector timestream. Units differ between the spectrum and residual/error columns. |
| Planck TT and EE binned tables | Angular-spectrum estimates with errors and a BestFit column | BestFit is a model prediction, not an additional observation. The binned tables do not supply a complete foreground/beam/bin-window/covariance likelihood. |
| Cosmicflows catalog and schema | Published redshift/distance products | Distance indicators, velocity corrections and group calibrations must be audited for the new model. |
| Atom/cavity paper text | A documented experimental comparison | Paper text is not the raw frequency time series or a ready-made secular-drift likelihood. |

All 21 selected paths exist. The Pantheon+ data hash matches its archived manifest. The compressed covariance was decompressed in memory; its uncompressed hash matches the archived original, and its 2,893,401 matrix entries agree with the table's 1,701 rows. This verifies the recovered products' identity and dimensions, not the suitability of every reduction assumption for a new theory.

## Prospective validation rules

The accompanying `validation-protocol.json` freezes process rules for future work:

1. All recovered and previously used products remain exploratory, including old train/validation/test partitions. Re-splitting exposed objects does not make a blind test.
2. Before a new candidate-selection run, freeze the action, finite parameter set, source/capture assumptions, data selection, transformations, likelihood, nuisance treatment, comparator and computational stopping rule.
3. Report predictions with their physically derived normalization separately from fits with free amplitudes. A fitted response factor is not a derivation of an energy source.
4. Use relevant covariance and shared calibration parameters. Disclose approximations and report galaxy-level alongside point-level performance.
5. Preserve all invalid/stability failures in the candidate's record. If an assumption changes, create a new candidate version and recompute dependent quantities.
6. Select genuinely unused validation data only after checking overlap and prior exposure. Until then, report exploratory evidence honestly.

No universal observational pass threshold has been invented. Candidate-specific likelihoods, sample selection and stopping budgets must be registered before the relevant fits. Existing numerical test tolerances check arithmetic and integration, not agreement with the universe.

## Missing evidence and next work

The registry identifies acquisition or derivation gaps for every requirement. High-priority gaps include resolved luminosity profiles and underlying photometry, galaxy distance/calibration provenance, lensing measurements, cluster profiles, local gravity and gravitational-wave data, background-radiation covariances, and a candidate-consistent cosmic history. A file being present is not sufficient if it contains a published fit rather than the observation needed to test a new law.

T02 remains incomplete because the detailed matter-frame convention, candidate-specific likelihoods and genuinely unused validation sample are not finalized. T03's candidate comparison work and the existing energy/transport diagnostics remain preliminary; none has been relabeled as observationally confirmed.

## Deliverables

- `input-audit.json`: hashes, presence, distance-method counts, per-object matched metadata and covariance checks.
- `observable-registry.json`: all 32 requirements with units, inputs, dependencies, ownership and remaining gaps.
- `validation-protocol.json`: prospective process rules and explicitly unfinished registration steps.
- `registry-checks.json`: coverage and reference-integrity checks.
- `audit_observables.py`, `enrich_observable_registry.py`: reproducible audit and registry builder using the current local project paths.

No archived input or original scientific result was modified.
