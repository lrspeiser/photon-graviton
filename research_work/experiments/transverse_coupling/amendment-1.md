# TF-1 amendment 1: the power-law fit range of the longitudinal scan
Declared 19 September 2026 after scan-v2 fell, before any rerun. Nothing in protocol.md is changed;
this note adds one refit and records two process facts.

## What failed
Stage 2 had three numerical gates: energy, reproduction, and "the fitted power laws (gain and each cost
against g_par over the scan) have residual scatter < 10% in the log", fitted over the scan points with
g_par <= 0.1. Energy and reproduction passed. The power-law gate FAILED as declared: the receiving gain
saturates between g_par = 0.1 and 0.2 and the timing discrepancy is not monotone there, so one power
law over 0.005 to 0.1 leaves log scatters of 0.144 (gain), 0.207 (timing) and 0.185 (colour); only the
fixed-ruler speed change fits at 0.053. The failed result stays in evidence/scan-v2, with its screens,
table and exponents as they fell. The point-by-point screens that carry the stage's reading do not use
the fit.

## The refit
The exponents are refitted over the weak-coupling points only, g_par in {0.005, 0.01, 0.02}, where every
cost is below one percent, with the same 10% log-scatter criterion, into a new directory scan-v3 through
`run.py --stage scan --fit-max 0.02`. The runs themselves are the same deterministic runs (scan-v2
reproduced scan-v1 bit for bit); no run, window, threshold or screen is altered. Three points give one
degree of freedom per fit, so the exponents are weak-regime descriptions, not laws. A pass of the refit
does not pass the original gate; both are reported.

## Two process facts, recorded
exact-v1 and scan-v1 completed their calculations and tripped the provenance guard, because the lens
module was edited while they ran; their results are unchanged by that edit (the exact and scan stages do
not import it) and the clean reruns exact-v2 and scan-v2 reproduce them bit for bit. The v1 directories
are preserved as failed-process records and are not cited as evidence.
