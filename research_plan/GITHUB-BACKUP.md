# Complete research checkpoint - 13 September 2026

The paper, scripts, canonical results and research plans are tracked on `main`.
Additional local research files excluded from normal Git tracking are preserved in
the same GitHub repository's [research backup release](https://github.com/lrspeiser/photon-graviton/releases/tag/research-backup-2026-09-13).

The adjacent `github-local-backup-2026-09-13.json` records the base commit, archive
names, SHA-256 hashes, original file paths and duplicate-content mappings. This
backup includes downloaded stellar catalogs, intermediate calculations, recovered
evidence, the Pantheon covariance file and partial computation checkpoints.
Partial results remain partial; historical direction/state files are not new
instructions or evidence that a computation finished.

## Restore

1. Clone this repository and check out `main` (or the recorded base commit for
   the precise scientific checkpoint).
2. Download every ZIP asset from the linked release. Verify its SHA-256 against
   the manifest before extraction.
3. Extract the ZIPs at the repository root, retaining their relative paths.
4. For manifest entries with `duplicate_of`, copy that restored file to the
   entry's `path`. All unique archived file contents have SHA-256 hashes in the
   manifest. For `tracked_duplicates`, copy `duplicate_of_tracked` to `path`;
   these files already exist in Git and were not archived again.

Python bytecode caches and temporary PDF rendering images are omitted: they are
reproducible working files, not unique research. Large datasets are release assets,
not Git objects on `main`; the versioned manifest links the two. The canonical
paper remains `output/pdf/theory-basis.pdf`.
