# Fresh-test eligibility: five excluded, 24 pending

This is the current stage-2 eligibility screen. It supersedes the initial three-exclusion/26-pending count for new work, while preserving the earlier identity, frame and radio-quality snapshots for reproduction. It does not change their historical results or certify new identities.

## Two ambiguities no longer block candidate screening

| Target | Evidence | Current fresh-test decision |
|---|---|---|
| LV J1017+2922 | 1.726 arcsec from CF4 PGC4231240, which is a previously exposed object/group; explicit alias remains unconfirmed | Precautionary exclusion |
| AGC740112 | 1.918 arcsec from CF4 PGC5808772, also previously exposed; explicit alias remains unconfirmed. An identity search additionally displayed a published target distance. | Precautionary exclusion plus documented pre-freeze distance exposure |

The same rule is applied to every one of the 29 originally staged candidates: an unresolved fixed-cone coincidence with a previously exposed object/group is sufficient to exclude it from a claim of fresh evidence. It is not sufficient to assert that the catalog identities are identical. The initial UGC05797, dw1046p1244 and NGC4592 exclusions remain. All 29 records, including these five exclusions, are retained with reasons.

In plain language, we do not need to prove two names refer to the same galaxy before deciding that a possible repeat is unsuitable as our clean test. This resolves their eligibility without choosing by prediction success or changing their distance measurements. They remain available for explicitly nonblind work, provided measurement provenance is established.

## Expanded historical scan

The saved single-counterpart support supplies recognized catalog alias patterns for 15 targets. A new scan of the same 799 tracked text files at pre-ELVES commit 3e71904 finds no additional matches under those aliases. This supports a narrower claim than freshness: none of those exact normalized catalog-name patterns was found in that historical text scope. It does not search all PDFs, external conversations, untracked downloads, coordinate-only identities or all possible host relationships. It does not undo any initial exclusion or count present audit text as old exposure.

## Unplanned search exposure

An identity search returned an excerpt from a [published galaxy distance table](https://www.sao.ru/hq/dolly/bs/table_1_150.pdf). The excerpt displayed AGC740112's distance, together with neighboring table rows. No target redshift or prediction residual was shown or computed. The accompanying elves-search-exposure.json records the source, query, all displayed object names and categories of displayed fields. Those objects must not be described in future work as having wholly unseen distance information. We do not use the displayed values for any fitting or selection threshold.

No numerical pair was scored in this pass, but saying that all candidate measurements remained unseen would now be inaccurate. The distance exposure gives a second reason to keep AGC740112 outside a blind test; it does not claim that its redshift had already been scored. Broad identity searches can expose measurements in snippets, so further pending-target source work should use explicit metadata projections where available.

## Reproduction and next action

Run `python research_work/results/redshift-priority/elves_quarantine_screen.py`. It verifies complete initial target coverage and the fixed positional radius, reads saved identity evidence, scans the frozen historical checkout in one batch, and outputs all decisions with input hashes. No candidate distance or redshift values enter the decision code. The underlying CF4 grouping is used only as a conservative historical exclusion convention, not as an adopted dark-matter environment model.

Current count: **29 initial candidates, five excluded, 24 pending, zero certified fresh**. New acquisition must take its pending names from elves-quarantine-screen.json. Earlier scripts that intentionally replay 26-target snapshots should remain pinned to their original inputs; do not silently rewrite those historical audit populations.

Continue source/host checks on the remaining 24 and settle a bounded, explicitly documented eligibility scope. The final common formula and credible uncertainty remain unfrozen; neither a clean identity nor an absence of historical name hits supplies an independent galaxy-motion model. No new rate was fitted, no fresh predictive improvement is established, and the full research goal remains active.
