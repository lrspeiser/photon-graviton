# GF-1 companion following

Start with [the report](report.md), [results figure](results.png), and the
[618-formula catalogue](catalogue.md). Machine-readable equations are in
[catalogue.json](catalogue.json). These are candidate variants, not 618 novelty
claims or a validated theory of gravity.

The [original protocol](protocol.md) was committed and pushed before calculation.
[Implementation details](implementation-notes.md) specify measurements.
The [initial branch review and attribution](review-and-provenance.md) is
supplemented by the [completed RW-1 review](rw1-completion-review.md).
[Execution note](execution-note.md) records an unavailable circular balance;
[scoring correction](scoring-correction.md) preserves a lag-window bug and its
correction. [Diagnostic plan](diagnostic-plan.md) distinguishes exposed
follow-ups from the original test.

## Evidence

- controls-v1: 26 original mechanics checks; controls-v2: 27, including the
  lag-window regression test.
- screen-v1: all 600 phenomenological chains/rings and 18 conservative chains.
  Its first invocation stopped on an unavailable conservative circular balance.
- screen-completion-v1: 17 conservative rings and the explicit unavailable case.
- scoring-v2: corrected 20-time-unit lag scores from original saved trajectories.
  No dynamics were rerun for this rescore.
- refine-v1 and diagnostics-v1: preserved first selections and follow-ups.
- refine-v2 and diagnostics-v2: follow-ups for the corrected selections.
- audit.json and evidence-sha256.json: independent mechanical diagnostics,
  source-commit comparisons and exact evidence byte digests.

Large original NumPy files are committed as lossless 32 MiB byte chunks plus
an index. No trajectory was thinned beyond the declared 0.2 recording cadence.
archive_io.open_npz reads an ordinary local NPZ or verifies and joins its
committed chunks in memory. The intact original large files may remain locally;
only those two regenerable byte-identical copies are ignored by Git.

## Reproduction

Requires Python 3.13, NumPy 2.2.6 and Matplotlib (rendering only). Set
OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. From the repository root:

    python -B research_work/experiments/companion_following/controls.py --output NEW_CONTROLS_DIRECTORY
    python -B research_work/experiments/companion_following/audit_and_plot.py

The first command refuses to overwrite an existing evidence directory.
The second regenerates derived audit/plot/checksum files from immutable
evidence, comparing source hashes against their recorded Git commits.

For a fresh dynamics reproduction, use a new output directory with
campaign.py screen. The original batch initializer intentionally raises if any
candidate has no positive circular balance; complete_screen.py documents and
executes the recovery used in this campaign without changing its fixture.
Use the recorded source commit for an exact historical rerun. New experiments
need new evidence directories and an explicit protocol; do not overwrite these.

Model equations use dimensionless declared units and instantaneous planar
interactions. No dark matter, expanding geometry, distance fitting, star-motion
fit or cluster-lensing fit is executed by GF-1. Initial packet energy is counted;
photon formation, microscopic graviton identity and matter/light coupling
remain unestablished.
