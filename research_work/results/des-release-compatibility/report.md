# Archived and pinned DES light-curve compatibility

The 2024 SNDATA_ROOT archive and our pinned later DES release contain identical light-curve measurements for all matched events. This resolves the measurement-version ambiguity in using the recovered historical calibration workflow. It does not establish that the complete calibration error budget or source-emission model is correct.

## Comparison performed

The full SNDATA_ROOT archive SHA-256 was reverified before extracting its DES HEAD, PHOT and README files. Current input hashes were checked against the pinned acquisition manifest (commit c9a4fcafc4cbd19bd750dee47fc76194a45c181f). Events were paired by SNID, with unique identifiers required on both sides. Within each event, point counts, observation times and canonical band labels had to agree in their existing order before measurement values were compared. No nearest-neighbor pairing or interpolated flux comparison was used.

| Check | Result |
| --- | --- |
| Archived events | 19,706 |
| Pinned-release events | 19,706 |
| Matched event IDs | 19,706 |
| Events present in only one version | 0 |
| Events with unequal point counts | 0 |
| Events with mismatched time/band alignment | 0 |
| Aligned measurement rows | 1,779,030 |
| Changed FLUXCAL entries | 0 |
| Changed FLUXCALERR entries | 0 |
| Changed entries in the other 16 shared measurement fields | 0 |

These are exact comparisons of stored values, not approximate statistical agreement. Separator rows outside event pointers are not counted as measurements. The compared fields include photometry flags, image identifiers, observing conditions, zero-point fields and pixel coordinates. No fields exist only in one PHOT schema.

## What did change

Eight host-galaxy magnitude/error columns differ in the HEAD tables. HOSTGAL_MAG_g differs for 16,697 events, HOSTGAL_MAG_i for 16,703, and each of the other six magnitude/error fields for 16,704. The archived header has IAUC; the newer header instead includes NAME_IAUC and adds NAME_TRANSIENT, LENSDMU and LENSDMU_ERR. This is a schema inventory, not a claim that these added quantities should be used in our physical model.

No other shared header field differed in this comparison after excluding identifier/pointer fields handled separately. Host magnitudes are scientifically relevant to population analyses, so the corrected pinned metadata should be retained with its provenance. They must not silently be replaced by archive values. Host magnitudes are not measurements of supernova intrinsic luminosity.

## Consequence for the joint-light work

The later release has not altered the measured light curves relative to those packaged with the recovered historical calibration resources. We can retain the later corrected metadata while developing the original-calibration forward comparison on the identical light curves. This does not authorize applying an unrelated newer recalibration, omitting cross-filter covariance, or treating a trained supernova template as an independently known clock.

The known calibration-offset convention, measured zero-point covariance and real passbands are now supported by direct measurement-version matching. Additional passband/source-model systematics and source-evolution constraints remain necessary. There is no redshift/timing/brightness fit in this audit and no new blind observational test; these DES data were already exposed.

The separate 160-case artificial timing calibration is unchanged and must finish with its declared failures retained. All six scientific demonstrations remain incomplete.

```powershell
python research_work/results/des-release-compatibility/run.py
```

The existing original-archive acquisition and pinned timing-input acquisition supply the two verified source products. The script retains extracted files only in ignored generated storage and commits comparison counts and hashes, not a duplicate light-curve release.
