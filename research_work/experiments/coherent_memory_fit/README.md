# Shared memory fits and observed comparison

Start with [the chart](observed-comparison.png), [report](report.md), [exact summary](summary.json), and [audit](audit.json).
The selected shared formula remains unsuccessful as a full galaxy/cluster solution.

## Reproduction
Python 3.13.5; NumPy 2.2.6; SciPy 1.16.1; Matplotlib 3.10.5.
Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. From the repository root:

    python -B research_work/experiments/coherent_memory_fit/audit.py
    python -B research_work/experiments/coherent_memory_fit/make_comparison.py
    python -B research_work/experiments/coherent_memory_fit/build_report.py

These verify evidence or rebuild derived artifacts, not the fitted archives.
For a fresh campaign, use an isolated checkout at the source commits:
- CMF-1: protocol f9a5928, source 9726084; run.py creates evidence-v1.
- CMF-2: protocol 7121a21, source 115f3fa; gradient_run.py uses the preserved first archive and creates evidence-gradient-v1.
Both drivers refuse to overwrite their archive directories. Keep Git history for provenance checks.

The 132 configurations comprise 84 plus 48 fits, with three starts each.
All partitions are exposed; no fresh validation, full energy conservation, microscopic attachment or historical originality is claimed.
No dark-matter model or expanding universe is active. No catalog distance is adjusted.
