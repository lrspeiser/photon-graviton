# Pending Integration Archive

This folder preserves research artifacts that were produced against earlier repository heads but were not safely merged before `main` advanced.

They are stored on `main` so no work is lost. **Do not apply these patches blindly.** Reconcile them against the current code and retain newer implementations where they supersede archived work.

## Archive

The exact source artifacts are packaged in one compressed archive split into eight ordered binary parts:

```text
archive/pending_integration_patches.tar.xz.part-000
archive/pending_integration_patches.tar.xz.part-001
archive/pending_integration_patches.tar.xz.part-002
archive/pending_integration_patches.tar.xz.part-003
archive/pending_integration_patches.tar.xz.part-004
archive/pending_integration_patches.tar.xz.part-005
archive/pending_integration_patches.tar.xz.part-006
archive/pending_integration_patches.tar.xz.part-007
```

Reassemble and verify with:

```sh
cat archive/pending_integration_patches.tar.xz.part-* \
  > /tmp/pending_integration_patches.tar.xz

sha256sum /tmp/pending_integration_patches.tar.xz
# expected:
# e2ee617149c985d213f7aae814c0c82520881822783d046eca9e7937accd1906

tar -xJf /tmp/pending_integration_patches.tar.xz -C /tmp
cat /tmp/MANIFEST.txt
```

Individual part checksums are recorded in [`archive/SHA256SUMS`](archive/SHA256SUMS).

## Included artifacts

### `phase_junction_stage5b.patch`

Local quadratic nonlinear constraint closure, gravitational self-energy sourcing, nonlinear-background mode counting, and CI changes.

SHA-256:

```text
eee19fbbb5104fe5cc0a5a7efc20707568d9e7b2abe892c43ea4c1805f987802
```

### `phase_junction_cone_breakthrough.patch`

Stage 6C regulator-complete matter stress, Stage 6D spacetime-coherent domain-wall/PV tests, cone-renormalization checkpoint, scripts, results, documentation, and workflow changes.

SHA-256:

```text
c77929cec9f62f27bc53eebb5591c7d63aa48de1a7ba738f595088ff6f0184e6
```

**Cleanup required before application:** the archived patch was previously noted to contain a compiled `__pycache__/*.pyc` entry. Remove generated binary/cache entries from any reconciled commit.

### `alternative_gravity_branch_portfolio.patch`

A branch-portfolio document expanding the research program beyond the earlier two-branch framing, including screened companion, memory/deposition, split quasistatic/radiative, preferred-frame, modified-inertia, massive/screened, and global/nonlocal branches.

SHA-256:

```text
8368b9bb31bea3108197e042c0e1ea38350339d45dd5a403a64ef2ec890f0f49
```

This artifact predates newer screened-memory work on `main`; reconcile it with the implemented SM-1 branch instead of treating it as the current active plan.

## Safe integration procedure

For each artifact:

1. Reassemble and verify the archive checksum.
2. Create a temporary branch from the current `main`.
3. Run `git apply --check` only as an initial diagnostic.
4. Inspect every overlapping path against current implementations.
5. Port scientific results and tests manually where later code has diverged.
6. Exclude compiled files, caches, and stale generated outputs.
7. Run the complete relevant regression workflows.
8. Merge only after updating the claim boundary and current-status documents.

## Archive status

These files are preserved evidence and integration inputs. Their presence in this folder does **not** mean their diffs have been applied to the live theory implementation.
