# Pending Integration Archive

This folder preserves research artifacts that were produced against earlier repository heads but were not safely merged before `main` advanced.

They are stored on `main` so no work is lost. **Do not apply these patches blindly.** Reconcile each patch against the current code and retain newer implementations where they supersede the archived work.

## Contents

### `stage5b/`

Local quadratic nonlinear constraint closure, gravitational self-energy sourcing, nonlinear-background mode counting, and CI changes.

Original file:

```text
phase_junction_stage5b.patch
```

SHA-256:

```text
eee19fbbb5104fe5cc0a5a7efc20707568d9e7b2abe892c43ea4c1805f987802
```

The patch is split into ordered 20,000-byte fragments. Reassemble with:

```sh
cat stage5b/phase_junction_stage5b.patch.part-* > /tmp/phase_junction_stage5b.patch
sha256sum /tmp/phase_junction_stage5b.patch
```

### `cone_breakthrough/`

Stage 6C regulator-complete matter stress, Stage 6D spacetime-coherent domain-wall/PV tests, cone-renormalization checkpoint, scripts, results, documentation, and workflow changes.

Original file:

```text
phase_junction_cone_breakthrough.patch
```

SHA-256:

```text
c77929cec9f62f27bc53eebb5591c7d63aa48de1a7ba738f595088ff6f0184e6
```

Reassemble with:

```sh
cat cone_breakthrough/phase_junction_cone_breakthrough.patch.part-* > /tmp/phase_junction_cone_breakthrough.patch
sha256sum /tmp/phase_junction_cone_breakthrough.patch
```

**Cleanup required before application:** the archived patch was previously noted to contain a compiled `__pycache__/*.pyc` entry. Remove any generated binary/cache entry from the reconciled commit.

### `alternative_gravity/`

A branch-portfolio document expanding the research program beyond the earlier two-branch framing, including screened companion, memory/deposition, split quasistatic/radiative, preferred-frame, modified-inertia, massive/screened, and global/nonlocal branches.

Original file:

```text
alternative_gravity_branch_portfolio.patch
```

SHA-256:

```text
8368b9bb31bea3108197e042c0e1ea38350339d45dd5a403a64ef2ec890f0f49
```

This artifact predates newer screened-memory work on `main`; reconcile the document with the implemented SM-1 branch instead of treating it as the current active plan.

## Safe integration procedure

For each artifact:

1. Reassemble and verify the SHA-256 checksum.
2. Create a temporary branch from the current `main`.
3. Run `git apply --check` only as an initial diagnostic.
4. Inspect every overlapping path against current implementations.
5. Port scientific results and tests manually where later code has diverged.
6. Exclude compiled files, caches, and stale generated outputs.
7. Run the complete relevant regression workflows.
8. Merge only after updating the claim boundary and current-status documents.

## Archive status

These files are preserved evidence and integration inputs. Their presence in this folder does **not** mean their diffs have been applied to the live theory implementation.
