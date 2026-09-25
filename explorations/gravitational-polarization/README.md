# Gravitational polarization explorations

**Exploratory, not adopted.** This directory is separate from the working theory, astronomical regression baseline, fitted constants and locked forecasts. No result here inherits the adopted model's galaxy or cluster scores.

## Archived GP-1 review

`gp1-review.tar.xz` preserves the complete 14-file gravitational-polarization review supplied in this conversation, including REPORT.md, checks.py, results.json, the numerical-integrity checks and the supplied independent source audit. It is a lossless repack of the uploaded ZIP, not a revision to its contents.

Extract and verify:

```sh
tar -xJf gp1-review.tar.xz
cd gravitational_polarization_review
sha256sum -c MANIFEST.sha256
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python checks.py --output reproduced.json
```

All original manifest entries were checked, and the original diagnostic results were independently rerun. They reproduce identically except elapsed execution time; all ten implementation-integrity checks pass. These are numerical checks, not evidence that a new gravitational theory is correct.

The review's main limitation remains: its cubic polarization energy is selected to obtain the cold scaling, not derived from independent microscopic physics. Motion sourcing, screening, persistent emission and a relativistic lensing metric are not established by GP-1.

Provenance is in `PROVENANCE.md`. Subsequent candidate work belongs in this directory and must retain failures and assumptions alongside successful checks.
