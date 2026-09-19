# Coupled conversion campaign (CWC-1)

Read [the report](report.md), [the pre-run plan](protocol.md),
[the attribution register](provenance.md) and
[the exploratory amendment](amendments.md).
No dark matter or expanding geometry is used in these runs.
The candidate is an effective fictional model, not a validated gravity theory.

Reproduce with Python3.13, NumPy2.2.6 and SciPy1.16.1. Matplotlib3.10.5 is
used only for the summary figure. Run from the repository root, with one BLAS
thread if running stages concurrently. Every evidence directory must be new.

```powershell
python -B research_work/experiments/coupled_conversion/run.py --stage small --output research_work/generated/cwc-small
python -B research_work/experiments/coupled_conversion/run.py --stage one_dimensional --output research_work/generated/cwc-1d
python -B research_work/experiments/coupled_conversion/run.py --stage two_dimensional --output research_work/generated/cwc-2d
python -B research_work/experiments/coupled_conversion/amend.py --stage refinement --output research_work/generated/cwc-refinement
python -B research_work/experiments/coupled_conversion/amend.py --stage explore --output research_work/generated/cwc-explore
python -B research_work/experiments/coupled_conversion/validate.py --stage mirror_audit --output research_work/generated/cwc-mirror
python -B research_work/experiments/coupled_conversion/validate.py --stage readiness --output research_work/generated/cwc-readiness
```

Refinement and mirror diagnostics compare new calculations with the committed
baseline evidence. They do not overwrite it. A scientific false result is
saved and is distinct from a failed process. Exceptions remain in manifests.

The original1D arrival-refinement failure and2D mirror-gate failure remain
archived. Follow-up diagnostics retain the original thresholds and separately
explain their results. Higher-power runs receive no hidden scalar seed.

The broad historical suite is not invoked: some historical comparison jobs
use assumptions excluded by the owner. CWC-1 has its own isolated runners and
checks loaded modules. Earlier PF5/RUT evidence is unchanged.
