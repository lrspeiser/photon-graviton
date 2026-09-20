# Observation-constrained coherent hitchhiking (IH-1)

Start with [the report and fitted equations](report.md), [figure](results.png),
[exact coefficients and scores](summary.json), and [independent audit](audit.json).

This treats the owner's mechanism as a working hypothesis and solves for an
effective deterministic response. It is not a completed field-generation or
photon-binding theory. No dark matter, expansion or distance changes are used.

## Preserved stages

| Stage | Protocol | Source commit | Archive |
|---|---|---|---|
| Positive coupling, 40 fits | [primary](protocol.md) at e521c76 | 4c23b2c | evidence-v1 |
| Opposite turning, 24 fits | [counter protocol](counter-protocol.md) at 260c8bb | 2732814 | evidence-counter-v1 |
| Large-acceleration release, 6 fits | [release protocol](recovery-protocol.md) at 5487d89 | d580cff | evidence-recovery-v1 |

Every fit has three archived starts. These are 35 response variants, not 210
different theories. Training, validation and test-labelled partitions are
historically exposed; later stages are explicitly exposed follow-ups.
Selections are written before each stage's test/Coma calculation.

## Reproduction

Python 3.13.5, NumPy 2.2.6, SciPy 1.16.1, Matplotlib 3.10.5.
Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. From the repository root:

    python -B research_work/experiments/inverse_hitchhiking/audit.py
    python -B research_work/experiments/inverse_hitchhiking/build_report.py

The first verifies archived bytes, original Git sources and the final
predictions independently. The second rebuilds derived report, summary and PNG.
Neither reruns the fits. For fresh fit reproduction use an isolated checkout
of the source commit above and run campaign.py, counter.py or recovery.py.
Those commits precede their respective result directories; the drivers refuse
to overwrite existing archives. Keep full Git history for source-hash checks.

All three fitting invocations completed normally. The first report-rendering
attempt had an unmatched-bracket syntax error; it was corrected before the
report was built. No simulation or observation archive was altered by that fix.

The final report distinguishes improved empirical matching from a statistically
adequate fit, absolute lensing calibration, conservation of a complete dynamical
field, and historical novelty. None of the latter is asserted.
