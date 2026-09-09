# Environmental time and companion-energy manuscript

Current working draft: **v0.2, 9 September 2026**. This is an incomplete theoretical framework with a reproduced observational calibration, not a validated replacement cosmology.

- [Manuscript PDF](manuscript.pdf) and [editable source](manuscript.md).
- [All 164 observed/predicted rows](analysis/observed-predicted-redshift.csv), [metrics](analysis/metrics.json), [comparison figure](analysis/redshift-comparison.png), [residual figure](analysis/redshift-residuals.png).
- [Evidence map](evidence-map.md) and [artifact verification](verification.json).

## Reproduce

From repository root, install `papers/cumulative-time-companions/requirements.txt`, then run:

```text
python papers/cumulative-time-companions/analysis.py
python papers/cumulative-time-companions/build.py
```

The analysis verifies the adopted source hash, independently refits only the inherited training partition and checks every prediction against the saved conversion-first result. It never alters the catalog. The inherited partitions are previously exposed, not fresh holdouts. The PDF builder embeds the generated figures and numbered equations.

## Revision scope

Version 0.2 clarifies the root law and its link to temporal evolution of the propagation factor; adds the full redshift comparison, residuals and limitations; corrects the schematic canonical field notation; and incorporates progress through the 45-job pair-production-balance verification. Evidence is pinned to c0425ecba9f13f937c56e0806b420539dd5c2fcd. The physical conversion rate, common clock action, persistent three-dimensional field, capture, supported halo and joint gravity/lensing solution remain open.

Version 0.1 remains available in Git history. Recovered historical papers and original data are unchanged. The full research program remains incomplete.
