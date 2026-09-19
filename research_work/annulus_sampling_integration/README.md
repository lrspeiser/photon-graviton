# Integration record for `research_work/annulus_sampling`

The owner's corrected annulus sampler (protocol 2f1d5ed, code bde7a7c) came with a six-item handoff. This
directory is the record of item 1, kept outside their package so that nothing in it is touched and its own
file manifest stays valid.

## What was run (2026-09-18, on main at e59ef6f)

| step of item 1 | result |
|---|---|
| pull main without overwriting an active local campaign | rebased onto 2f1d5ed and bde7a7c and pushed without force; my one local edit to the historical `equilibrium.py` (a docstring) was reverted first, so that file is byte-identical to the blob the owner verified (`3f76b359`) |
| run the focused checks | `python -B -m research_work.annulus_sampling.checks`: 15 of 15 pass |
| finish the four-family audit | `four-family-audit.json` here: all four stage 6 families pass, none failed |
| run the historical suite | 73 of 73 jobs pass with the saved results preserved: the 72 historical jobs unchanged, and the new one |
| register the focused checks as a separate numerical-verification job | `run_checks.py` runs them as job 73, as a module from the repository root, outside the frozen results tree |
| do not edit archived rut6 to use the new sampler | not edited: `rut6.py` and `equilibrium.py` are unchanged and stage 6 replays against its archive with zero differences, failed H6 included |

## The four-family audit

`verification.json` in the owner's package records the two narrow families. The audit tool was run here
over all four (`--all`); it refuses to write into the frozen results tree, and `research_work/generated/`
is ignored by git, so its output is copied here unchanged.

| family | radial marginal, max abs error | deliberate extra-R control | mean-radius shift of the control | Var(R)/Mean(R) |
|---|---|---|---|---|
| narrow_cold | 5.4e-16 | rejected | 0.001733356 | 0.001733356 |
| narrow_warm | 4.7e-16 | rejected | 0.003085551 | 0.003085551 |
| wide_cold | 3.6e-16 | rejected | 0.059165297 | 0.059165297 |
| wide_warm | 3.6e-16 | rejected | 0.081275658 | 0.081275658 |

The shift of the biased sampler equals Var(R)/Mean(R) to the last digit in every family, which is the
owner's analytic statement of the defect. It also says how large the stage 6 error was: small for the
narrow annuli, and 6 to 8 percent of the mean radius for the wide ones.

The audit's own status block applies here too: numerical verification passed; physical stability and
continuum convergence **not evaluated**; the historical H6 **failed and preserved unchanged**.

## What this is not

It is not items 2 to 6. Those are stage 7 (`research_work/results/path-memory/protocol-rut7.md` and its
amendment), which uses this sampler as the sampler of record at every source-sampling call site.
