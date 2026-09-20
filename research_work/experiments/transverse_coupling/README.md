# TF-1: transverse coupling of the converted field

Read [the protocol](protocol.md), [the attribution register](provenance.md), [amendment 1](amendment-1.md)
and [the report](report.md). No dark matter, expansion, distance refit or halo model is used.
The candidates are effective fictional couplings, not validated gravity.

Run from the repository root with Python 3.11+, NumPy and SciPy (versions are recorded in each
manifest). Every evidence directory must be new; failures are preserved.

```sh
python -B research_work/experiments/transverse_coupling/run.py --stage exact  --output research_work/experiments/transverse_coupling/evidence/exact-v1
python -B research_work/experiments/transverse_coupling/run.py --stage scan   --output research_work/experiments/transverse_coupling/evidence/scan-v1
python -B research_work/experiments/transverse_coupling/run.py --stage scan   --fit-max 0.02 --output research_work/experiments/transverse_coupling/evidence/scan-v3   # amendment 1
python -B research_work/experiments/transverse_coupling/run.py --stage lenses --output research_work/experiments/transverse_coupling/evidence/lenses-v1
python -B research_work/experiments/transverse_coupling/tf1_checks.py          # the suite job, about a minute
```

Evidence directories exact-v1 and scan-v1 are failed-process records (provenance guard tripped by an
edit during the run); exact-v2, scan-v2, scan-v3 and lenses-v1 are the cited runs.

The exact stage reuses CWC-1's archived generated field by hash and CWC-1's own index-ray code for
the live comparison; the scan reuses CWC-1's `spatial.simulate` unchanged; the lens stage reuses
CL-2's `LensSystem` (CR-2's measurement interface) under the static geometry.
