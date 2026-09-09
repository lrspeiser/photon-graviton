# Current research

Read the [research checkpoint](results/RESEARCH-CHECKPOINT.md) and [roadmap](../research_plan/START-HERE.md). The theory is incomplete; passing these checks verifies the stated calculations, not observational agreement with the universe.

## Layout

- `energy/`: reusable energy-ledger and radial-transport modules.
- `results/`: canonical saved reports, numerical evidence, figures and diagnostic scripts, organized by subject.
- `run_checks.py`: portable entry point for the current diagnostics.
- `generated/`: ignored output directories from fresh runs. Committed evidence is not overwritten.
- `publication-verification.json`: verification performed while preparing the GitHub update.

Local orchestration notes and duplicate `evidence/run-*` copies are not part of the public source tree. Scientific results from those runs are consolidated under `results/`.

## Run

From the repository root, using a virtual environment with `requirements.txt` installed:

```sh
python research_work/run_checks.py
python research_work/run_checks.py --baseline
```

The first command runs 27 jobs: catalog no-loss comparison, closed and open energy ledgers, internal/external spatial capture, microphysical kinematics, frequency transfer, scalar-wave evolution, gravity/lensing response, the data/requirement audit, the conversion-first empirical pass, interaction-rate constraints, matter-assisted conversion, collective spatial response, oscillator frequency response, ordinary soft-graviton emission, reversible capture/storage, mechanical deposit support, capture momentum, moving-receiver isotropic capture, broadband receiver response, joined source/receiver spectral kinetics, dynamical matched electromagnetic waves, thermal conversion with an exploratory FIRAS comparison, background replenishment source/energy requirements, and heated thermal absorption/emission with an exploratory FIRAS comparison. The second adds the historical baseline's three scripts as one additional job, run in an isolated copy. Each job logs its output, and the runner checks that saved results and the 271 recovered original files stay unchanged.

To select a new output location:

```sh
python research_work/run_checks.py --baseline --output-dir research_work/generated/my-run
```

Existing nonempty output directories are refused. Script failures or integrity failures produce a nonzero exit code. Numerical differences in the historical baseline are recorded for review rather than hidden by relaxing tolerances; a successful process exit is not a claim of bit-identical regression. See [baseline reproduction](results/baseline/reproduction-report.md).

The `PHOTON_GRAVITON_RESULTS` environment variable is set by the runner for its child scripts. Individual scripts may require preceding outputs; use the entry point to preserve their order. No network download is required for these numerical runs after Python dependencies are installed.

The diagnostics' mathematical logic and original result records are preserved. Only current scripts were made portable; old paper builders and historical checkpoints retain their documented environment assumptions. Archived JSON containing paths describes the original calculation, not an active background process.
