# New calibration, previously represented galaxy groups

Decision: the audited targets in [Jensen et al. 2025, TRGB-SBF Project III](https://arxiv.org/abs/2502.15935) do not provide a fresh group-independent prediction sample. Keep this release as a possible calibration-sensitivity comparison. No new redshift prediction was scored and no adopted distance was changed.

The paper recalibrates 61 earlier distant SBF measurements and describes an expanded HST sample as future work. Its Table 1 supplies calibration galaxies, not a new large distant evaluation catalog. The calibration is tied to NGC 4258. The metadata and selected methods paragraphs are exposed; the NGC 4874 paragraph also exposed its old/revised distances. None of this release may now be represented as completely unread.

## Completed identifier audit

All 14 Table 1 galaxy names plus NGC 4258 and NGC 4874 were selected before alias lookup. SIMBAD primary identifier responses resolve all 16 to PGC/LEDA identifiers. Exact CF4 identifier joins give:

| Targets | Previously excluded CF4 group |
|---|---:|
| IC 2006; NGC 1344, 1374, 1375, 1380, 1399, 1404 | 13418 |
| NGC 4458, 4472, 4489, 4552, 4649, 4697 | 41220 |
| NGC 4636 | 42734 |
| NGC 4258 | 39600 |
| NGC 4874 | 44715 |

Thus 16/16 targets overlap previously excluded groups, while only 1/16 matches a PGC object directly in the structured historical registry. Checking only exact object names would substantially understate shared environments and calibration exposure. This uses existing CF4 grouping as a conservative exclusion convention, not as a proof of theory-independent gravitational membership.

The JSON retains per-object identifier-section evidence, source URLs and the hashes of the PDF, decision protocol and pre-audit registry. verify_sbf_identity_audit.py replays all aliases and group joins against the registry in commit 539fe59, preventing later additions from manufacturing historical overlap. Only the first 15 characters of each raw CF4 row enter the group join; no target redshift or distance affects selection.

## Limits and next source decision

This audit covers the calibration table and the two named anchor/Coma objects, not an individually crossmatched list of all 61 revised distant measurements. It therefore does not assert that every revised galaxy was previously used. Nonetheless the publication's recalibration of an older sample and already exposed aggregate results do not satisfy the requested genuinely unexposed final test. A separately released expanded sample would require its own metadata and identity audit.

No new unexposed validation sample is certified. Broaden the source search to published geometric or standard-candle distance releases, including nearby objects outside the current distance range, while keeping distant-sample acquisition open. Reject target-redshift-derived distances as independent validation and preserve uncertainty about the physical environment and motions. Do not substitute calibration changes or old targets for a fresh predictive test.
