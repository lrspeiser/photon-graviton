# Photon-companion hitchhiking

Start with [the comparison report](report.md) and [actual results figure](comparison.png).
This is a fictional coupling experiment, not a derived physical graviton bound state.

## Records

- HH-1 [protocol](protocol.md), source at f9eaaed, evidence-v1: heading memory.
- HH-2 [protocol](protocol-curvature.md), source at 43d892c, evidence-curvature-v1:
  curvature memory. The changed law was declared after HH-1's results.
- [Numerical follow-up](numerical-followup.md), source at e3b2760,
  evidence-numerical-v1: smaller steps, mirrors and event ledgers.
- [Comparison data](comparison.json) and [integrity audit](integrity-audit.json).
- Three evidence SHA256 lists retain original bytes, including failed results.

Each mean campaign uses 81 parameter combinations, not 81 independent theories.
Mean trajectories and discrete event histories are different approximations.
The report keeps nonarrivals and unresolved resolution tests visible.

## Verify or reproduce

Requires Python 3.13, NumPy 2.2.6 and Matplotlib 3.10.5.
From the repository root, set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1:

    python -B research_work/experiments/photon_hitchhiking/audit.py
    python -B research_work/experiments/photon_hitchhiking/build_report.py

The audit reads evidence and pinned Git objects; it does not rerun dynamics.
The report command regenerates derived comparison JSON, Markdown and PNG only.
The audit also checks GF-1 and MS-1 evidence byte digests.

For an exact fresh dynamics reproduction, use an isolated checkout of the source
commit above, with the relevant experiment's output directory absent, then run
hitch.py, curvature.py or numerical.py as appropriate. Every dynamics driver
refuses to overwrite an existing evidence directory. Do not delete the canonical
archived results to make a rerun work. The Git-committed source is the reproduction
authority; later print-only fixes and optional curvature support do not rewrite
the source manifest of an earlier run.

The comprehensive report distinguishes completed numerical campaigns from work
still required to produce energy-conserving stellar and cluster-lensing predictions.
