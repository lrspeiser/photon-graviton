# Recovering the survey-membership information for the vertical test

**The current training sample mixes substantially different targeting populations.** Targeting metadata were recovered for all 27,884 existing training stars. There are 16,221 EXTRATARG=0 stars and 11,663 with other flags; no unmatched or conflicting EXTRATARG classifications remain. Original catalogs, training definitions and holdout roles are unchanged.

## Documented classification

SDSS identifies EXTRATARG=0 as its convenient Main Red Star Sample indicator. This flag does not supply a selection probability and can overlap special targeting categories. Selection fractions require photometric parent counts by field, color and magnitude cohort. Converting these to distance-dependent completeness additionally requires dust and stellar-population modeling. [SDSS selection documentation](https://www.sdss4.org/dr17/irspec/targets/selection-biases/), [targeting flags](https://www.sdss4.org/dr17/irspec/targettingbits/).

We matched the saved training source IDs to APOGEE_ID, FIELD and TELESCOPE, then recovered targeting flags and J/H/K and targeting-limit metadata from the local original allStarLite file. Duplicate candidates would be retained and conflicting flags marked ambiguous, rather than silently selecting a row. Cached candidate and membership tables are hash recorded. No validation or final-test kinematics were analyzed.

| Existing region | Stars | Main Red Star Sample | Other EXTRATARG |
|---|---:|---:|---:|
| all | 27884 | 16221 | 11663 |
| bulge_plane | 299 | 261 | 38 |
| bulge_offplane | 577 | 501 | 76 |
| disk_plane | 4425 | 3352 | 1073 |
| disk_offplane | 11214 | 5544 | 5670 |

Region definitions are unchanged from the earlier vertical decomposition: bulge R<3.5 kpc, disk R>=5 kpc, plane |z|<0.2 kpc, off-plane 0.5<=|z|<=1.5 kpc within the existing overall cuts. The four displayed regions do not partition the entire training set. This table reports targeting composition, not velocity dispersions or physical density ratios.

## What is now possible and what remains

The recovered flags allow a documented main-sample branch while retaining the original mixed sample for sensitivity comparison. They do not establish a representative sample: field/cohort completion, the photometric parent population, extinction, Gaia matching/quality, StarHorse availability, chemical selection and association uncertainties still affect inclusion. The original sample is not overwritten or silently filtered.

Next, retrieve the field-specific targeting parent products and cohort/design definitions for the main-sample branch, and quantify the extra selection introduced by this project. Raw main-sample counts also cannot yet be treated as intrinsic density. This is necessary groundwork for comparing the competing vertical-gravity predictions with observations; no new gravitational fit or companion mechanism is validated. All six scientific goals remain open.
