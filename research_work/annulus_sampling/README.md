# Corrected annulus sampler and local handoff

## What is implemented

The new sampler removes the extra radius factor from the phase-space weights.
It exposes its discrete node masses, radial marginal, and moments; uses radial
trapezoidal weights (including endpoints); retains the historical density
routine's equal rectangular angular-momentum weights; and draws both radial
velocity signs and a uniform angle. It validates inputs and never mutates the
source while sampling.

The adapter inherits the existing equilibrium solver and replaces **only**
`sample()`. The old `results/path-memory/equilibrium.py` and `rut6.py` are
intentionally unchanged: historical stage-6 reproduction still uses the old
sampler and still records failed H6. Follow-on studies must import the new
adapter or call the corrected standalone function. This is not a silent
correction to the old archive, a new force law, or a stability claim.

Protocol: `protocol.md`, committed before the new checks at `2f1d5ed`.
Baseline: `72aef4456aae63bbb3482706a5ecda6af77422c8`.

## Run from the repository root

```sh
python -B -m research_work.annulus_sampling.checks
python -B -m research_work.annulus_sampling.audit --cases narrow_cold narrow_warm --output research_work/generated/annulus-sampling-v2-narrow.json
python -B -m research_work.annulus_sampling.audit --all --output research_work/generated/annulus-sampling-v2-all.json
```

Use a new output filename on each run. The audit refuses to overwrite a file or
write into the frozen `research_work/results` tree. It requires NumPy and SciPy,
not the external/local observational databases. `--all` re-solves all four
archived matched-support construction families; the inherited wide-annulus
fixed-point solves can be substantially more expensive than the narrow ones.
Construction exceptions and numerical failures are recorded and make the audit
exit nonzero; omitted/failed cases cannot pass as a reduced sample.

## Use in the next physical experiment

```python
from research_work.annulus_sampling.adapter import WarmAnnulus

a = WarmAnnulus(dL=0.06, dE=0.03).solve(target_support=0.1)
check = a.verify_sampling()
if not check['numerical_verification_passed']:
    raise RuntimeError(check)
x, v = a.sample(256, seed=1)
```

For an already constructed historical object, there is no need to rebuild it:

```python
from research_work.annulus_sampling import build_distribution, verify_distribution

nodes = build_distribution(existing_annulus)
check = verify_distribution(existing_annulus, nodes)
if not check['numerical_verification_passed']:
    raise RuntimeError(check)
x, v = nodes.draw(256, seed=1)
```

Keep one `nodes` object when drawing several seeds from the same state; do not
cache it across a change to the source. Verify the equilibrium's own construction
residual independently before starting a simulation. A successful sampling
check does not imply that an equilibrium has converged.

## Verification performed here

- Fifteen focused unit/integration tests pass, including the existing
  WarmAnnulus adapter, an incomplete-gamma analytic density, nonuniform radius
  grids, deterministic moments, seeded Monte Carlo draws, input validation,
  source immutability, and output safety.
- The positive control deliberately inserts the old extra radius factor and is
  rejected. Its mean-radius bias equals Var(R)/Mean(R).
- Both full narrow-annulus constructions were re-solved and audited at the
  archived default resolution. The largest radial-marginal discrepancy was
  4.86e-16; mass and radial-moment discrepancies were at roundoff.
- Two attempts at the all-four audit exceeded this environment's execution
  budget during the inherited wide-annulus solve. No all-four output or verdict
  is claimed. The local worker should complete `--all`.
- The loaded historical equilibrium.py was byte-verified against Git blob
  `3f76b3591916c32e5d8d2a9665ca35244239dce0` before testing.
- The full 72-job historical suite was **not run here**. No existing tracked
  source, protocol, result, suite runner, or global manifest was modified by this
  package. `verification.json` is a local verification record, not a replacement
  for any historical archive or gate.

## Status policy

Keep three statuses separate:

1. Reproduction: whether an archived calculation can be reproduced.
2. Numerical verification: whether a sampler, integrator, or solver passes its
   correctness checks. These checks exit nonzero on failure.
3. Physical outcome: whether a declared model meets its scientific objective.

The audit explicitly marks physical stability and continuum convergence as
`not_evaluated`, and historical H6 as `failed_preserved_unchanged`. A reproducible
negative physical result is legitimate. A reproducible numerical bug is not a
passing numerical test.

## Remaining work for the local agent

1. Pull main without overwriting an active local campaign. Run the focused
   checks, finish the four-family audit, then run the historical suite. Register
   the new focused checks as a separate numerical-verification job in a new
   integration commit; do not edit archived rut6 to use the new sampler.
2. Declare a new experiment/output set for corrected cold/warm runs. Use the
   new adapter, or the standalone sampler, at every source-sampling call site.
   Keep H6 and the old outcomes frozen; do not regenerate them with corrected
   samples and present them as the original experiment.
3. Define invariant energy/angular-momentum cutoffs or a bound-energy taper.
   The inherited Gaussian and finite integration domain are not a demonstrated
   globally finite-mass distribution. Independently refine spatial/velocity
   quadrature and outer domain; a small fixed-point residual is not that check.
4. Compare identical corrected initial samples under live and frozen fields at
   common times, with multiple seeds and body counts. Record signed drifts and
   coherent streaming separately from residual dispersion. Do not divide two
   near-zero drift measurements and call their ratio a stability result.
5. Keep matched-support examples separate from fixed-physics temperature tests.
   The former change fitted alpha and total mass. The latter require fixed mass,
   alpha, kernel, formation and retention times, plus explicit DF normalization.
6. After construction and sampling are verified, proceed to warm-population
   modes, full-state perturbations, formation from empty fields, and reciprocal
   energy/momentum accounting. Use the locally available databases only for
   explicitly declared observational comparisons; this patch supplies no new
   galaxy fit, reservoir model, or proof of stability.

The new package has its own file manifest so it can be checked without touching
historical manifest pins. The local integration can append its new files to the
project-wide manifest without refreshing any old scientific input hash.
