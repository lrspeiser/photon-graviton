# Environmental time and companion-energy manuscript

Current working draft: **v0.4, 9 September 2026**. This is an incomplete theoretical framework with a reproduced observational calibration, not a validated replacement cosmology.

- [Manuscript PDF](manuscript.pdf) and [editable source](manuscript.md).
- [All 164 observed/predicted rows](analysis/observed-predicted-redshift.csv), [metrics](analysis/metrics.json), [comparison figure](analysis/redshift-comparison.png), [residual figure](analysis/redshift-residuals.png).
- [Evidence map](evidence-map.md) and [artifact verification](verification.json).

## Reproduce

From repository root, install `papers/cumulative-time-companions/requirements.txt`, then run:

```text
python papers/cumulative-time-companions/analysis.py
python papers/cumulative-time-companions/extended_analysis.py
python papers/cumulative-time-companions/build.py
```

The analysis verifies the adopted source hash, independently refits only the inherited training partition and checks every prediction against the saved conversion-first result. It never alters the catalog. The inherited partitions are previously exposed, not fresh holdouts. The PDF builder embeds the generated figures and numbered equations.

## Revision scope

Version 0.2 clarifies the root law and its link to temporal evolution of the propagation factor; adds the full redshift comparison, residuals and limitations; corrects the schematic canonical field notation; and incorporates progress through the 45-job pair-production-balance verification. Evidence is pinned to 7b0b296. The physical conversion rate, common clock action, persistent three-dimensional field, capture, supported halo and joint gravity/lensing solution remain open.

Version 0.1 remains available in Git history. Recovered historical papers and original data are unchanged. The full research program remains incomplete.

Version 0.3 labels all 15 numbered equations by provenance and adds the predeclared sky-tile comparison, fixed-rate six-maser diagnostic and conditional homogeneous-field result. [Sky-tile predictions](analysis/sky-tile-predictions.csv), [maser comparison](analysis/maser-comparison.csv), and [extended checks](analysis/extended-checks.json) retain all rows. No fresh-validation success or originality claim is made.

Version 0.4 incorporates evidence through ee5733b: 20 provenance-labeled equations, larger-region validation, [gearing targets](analysis/observed-time-stretch-targets.csv), [observation-factor diagnostics](analysis/observation-factor-diagnostics.csv), [larger-region predictions](analysis/coarse-sky-predictions.csv), universal local-clock/force checks and fresh-sample exclusions. [Revision checks](analysis/revision-checks.json) recompute saved metrics without fitting new parameters. The 26 pending ELVES targets remain unscored; their metadata alone does not establish a fresh sample.
