# TF-1: transverse coupling of the converted field

Read [the protocol](protocol.md), [the attribution register](provenance.md) and, once the stages
have run, [the report](report.md). No dark matter, expansion, distance refit or halo model is used.
The candidates are effective fictional couplings, not validated gravity.

Run from the repository root with Python 3.11+, NumPy and SciPy (versions are recorded in each
manifest). Every evidence directory must be new; failures are preserved.

```sh
python -B research_work/experiments/transverse_coupling/run.py --stage exact  --output research_work/experiments/transverse_coupling/evidence/exact-v1
python -B research_work/experiments/transverse_coupling/run.py --stage scan   --output research_work/experiments/transverse_coupling/evidence/scan-v1
python -B research_work/experiments/transverse_coupling/run.py --stage lenses --output research_work/experiments/transverse_coupling/evidence/lenses-v1
python -B research_work/experiments/transverse_coupling/tf1_checks.py          # the suite job, about a minute
```

The exact stage reuses CWC-1's archived generated field by hash and CWC-1's own index-ray code for
the live comparison; the scan reuses CWC-1's `spatial.simulate` unchanged; the lens stage reuses
CL-2's `LensSystem` (CR-2's measurement interface) under the static geometry.
