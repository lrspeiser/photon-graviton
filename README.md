# Photon–Companion Research in a Hypothetical Universe

**This project explores a fictional universe that has the same observed stars, galaxies, and astronomical datasets as ours, while allowing different physical laws.** We use real measurements as a common test bed so that we can safely explore new concepts, derive their consequences, and compare their predictions without presenting speculation as established physics.

The working universe is nonexpanding. Candidate mechanisms let light redshift through changing propagation and exchange energy with hypothetical companion waves or fields. Captured energy or a modified field response might contribute to galaxy and cluster gravity. These are hypotheses to investigate, not discoveries about our universe. “Companion” does not mean an experimentally detected particle or an established graviton interaction.

The aim is to find a self-consistent emitter–propagation–receiver model that predicts color, duration, brightness, angular size, gravitational-wave propagation, galaxy rotation, and lensing with shared independently constrained parameters. No candidate in this repository has yet achieved all of those requirements.

## Start here

1. **[Current research checkpoint](research_work/results/RESEARCH-CHECKPOINT.md).** Verified results, open questions, and links to the latest energy, capture, microphysics, data, wave and gravity derivations.
2. **[Research roadmap](research_plan/START-HERE.md).** All 32 requirements and 20 tasks, with prerequisites and completion criteria. Baseline reproduction is complete with documented numerical differences; the full theory remains incomplete.
3. **[Run the current checks](research_work/README.md).** Portable diagnostics write to a fresh ignored directory and preserve saved evidence.
4. **[Recovered causal-supply analysis](companion_causal_test/report.md).** DustPedia inputs, cross-validation, cluster transfer, controls, uncertainty checks and circulation. This is the historical baseline, not the latest full research checkpoint.
5. **[Main paper, version 9](redshift_paper/temporal_redshift_paper.docx).** The historical unified paper; it has not incorporated all later findings.

The [universe contract](research_plan/universe-contract.md) fixes published galaxy distances as facts, including adopted Hubble-flow values. Dark matter, cosmic expansion and the Big Bang are excluded as active explanatory premises. We retain observations and derive the photon-supplied mechanism without importing those theories. Propagation, capture, storage and the detailed gravity law remain under investigation.

## Current findings

- **Latest: [observer-time evidence](research_work/results/timing-foundation/report.md).** Pinned DES calibrated photometry is now available through a reproducible acquisition script and verified local cache. The audit covers 1,779,030 observations without time rescaling. Timing inference still requires a frozen selection/normalization method, injection tests and shared-reference uncertainty treatment; no new time-dilation measurement is claimed.
- **Conversion remains incomplete.** The [164-group exploratory fit](research_work/results/conversion-first/report.md) gives a constant-loss coefficient about 7.63e-5 per million light-years on reused data. Scalar interaction tests constrain color dependence and spectral broadening; [ordinary soft-graviton emission](research_work/results/soft-graviton/derivation.md) has very weak energy throughput in its tested limit. No complete microscopic redshift law is established.
- **Capture is more than absorption.** [Storage calculations](research_work/results/capture-storage/derivation.md) account for release, reverse transitions and finite capacity. [Mechanical support](research_work/results/deposit-support/derivation.md) and [momentum transfer](research_work/results/capture-momentum/derivation.md) show why retaining energy does not guarantee its location or orbit. [All-direction illumination](research_work/results/isotropic-capture/derivation.md) is not automatically momentum-balanced for a moving receiver.
- **Supply remains conditional and unresolved.** In the recovered 26-galaxy comparison, converting all present luminosity for 10 billion years falls short of the ordinary deposited mass-energy requirement by a median factor about 5,400. The local-redshift-only shortfall is much larger; catalogued neighbors and repeated circulation do not create extra energy. See the [recovered supply report](companion_causal_test/report.md).
- **Shape is not normalization or stability.** [Fast local capture](research_work/results/radial-capture/capture-shape-findings.md) can yield a conditional flat outer speed contribution if deposits remain where formed. Constant absorption of uniform external illumination produces uniform or outer-heavy deposition. Neither gives a formed, supported, adequately funded galaxy.
- **Gravity must predict motion and lensing together.** Earlier empirical galaxy/cluster fits do not derive their strength from photon energy. The [scalar gravity comparison](research_work/results/gravity-response/motion-and-lensing.md) shows why strengthening stellar motion need not supply matching light bending. Clock protection, time stretching and the nonexpanding background remain mandatory checks.

The [research checkpoint](research_work/results/RESEARCH-CHECKPOINT.md) indexes every retained derivation and numerical pass. Earlier propagation, circulation and fit findings remain in their original reports; all 20 tasks and 32 observational areas remain in scope.

## Scientific boundaries

Real datasets are observational constraints; their published distance calibrations, stellar-population assumptions, hydrostatic assumptions, and other model dependencies still matter in an alternate-universe interpretation. We do not silently treat model-derived quantities as raw observations.

We distinguish:

- mathematical derivations under explicit assumptions;
- energy-accounting identities versus demonstrated interactions;
- descriptive fits versus predictions with frozen parameters;
- exploratory reuse of a dataset versus genuinely new holdouts;
- fitted gravitational response versus independently calculated companion supply;
- a fictional model's internal consistency versus evidence that it describes our universe.

A good fit is not proof. Negative results and superseded candidates are retained so that the research history remains reviewable. The model does not yet supply exact local clock protection, a complete photon–matter–gravity action, a stable nonexpanding background, or joint lensing predictions. CMB origin, rotation, and topology are exploratory branches.

## Repository map

| Directory | Work preserved |
|---|---|
| `research_plan` | Current roadmap, task queue, assumptions and original planning provenance |
| `research_work` | Current solver modules, portable verification command, canonical reports and saved results |
| `archive/original-uploads` | Original GitHub upload bundles and historical paper/report, retained for provenance |
| `redshift_paper`, `redshift_paper_sources` | Main paper, initial group data and original redshift work |
| `option3_test`, `option3_cliff` | Galaxy timing/cliff candidates and prediction audits |
| `time_revision`, `time_first_principles`, `temporal_hypothesis_checks` | Signal stretching, thermal tests and early first-principles checks |
| `screened_time_candidate`, `unscreened_clock_link` | Screened and unscreened hypotheses, Solar System checks |
| `temporal_candidate_audit` | Broad candidate audit and source manifests |
| `time_light_atoms`, `fixed_atom_time_field`, `minimal_clock_interaction` | Atomic standards and photon propagation |
| `shared_interaction_test`, `shared_curved_update` | Interaction, brightness, shared gravitational-wave propagation and curvature |
| `nonrotating_action_test`, `conformal_action_derivation`, `interaction_stability_attempt` | Action candidates and stability limitations |
| `five_candidate_tests`, `alternate_universe_extensions` | Alternative mechanisms, archived for comparison |
| `focused_history_test`, `theory_closure`, `nature_tests` | Revised history, dynamics and natural/laboratory checks |
| `companion_wave_test`, `companion_deposition_fit` | Energy supply, halo shapes and global fits |
| `companion_causal_test` | Latest source-driven tests, controls, bootstrap and circulation calculations |
| `abstract_revision`, `paper_fixed_atom_update` | Available paper revisions; five historical `unified_paper_v9` build files remain missing |

The large Pantheon covariance input is stored losslessly as `.cov.gz`; `restore_data.py` restores and checks it before running historical brightness analyses. Its manifest entry records the original uncompressed checksum.

The recovered relative directory structure is preserved because analysis scripts share inputs. `SNAPSHOT_MANIFEST.json` records the original 276-file snapshot, of which 271 files are recovered with matching checksums; it does not inventory later research additions. See [recovery details](RECOVERY_REPORT.md). Original upload bundles remain under `archive/original-uploads`; local orchestration state, duplicate working evidence and regenerated runs are ignored. `.gitattributes` preserves recovered file bytes across platforms for checksum verification.

## Reproduce the latest analysis

Use Python 3.11 or later in a virtual environment:

```sh
python -m venv .venv
# Activate the environment using your operating system's usual command.
pip install -r requirements.txt
python research_work/run_checks.py
# Include a fresh isolated reproduction of the historical three-script baseline:
python research_work/run_checks.py --baseline
```

The current command checks original source hashes, runs diagnostics, and verifies that original files and saved results remain unchanged. Each invocation creates a fresh directory under `research_work/generated/`, with logs and a verification summary. `--output-dir` selects another new or empty directory. The latest [companion-bath verification](research_work/companion-bath-verification.json) records 31 current default jobs; the optional historical baseline makes 32 jobs and was not rerun for this update. Earlier [conversion-first verification](research_work/conversion-first-verification.json) and [publication verification](research_work/publication-verification.json) retain their original suite scope and results.

The historical scripts can still run directly, but doing so overwrites their corresponding outputs. Use the isolated runner above for baseline reproduction. `python restore_data.py` restores the compressed covariance if needed by older brightness analyses; the current diagnostics read its compressed form directly. These are research tools, not a complete cosmology library.

Earlier checkpoints may depend on additional libraries, archived inputs, network services, or document-rendering utilities. Their environment records and reports are preserved; not every historical script has been rerun for this commit. Paper-building scripts include original workspace paths and need adjustment for a different machine. Running all old paper-update scripts in sequence is not a supported build workflow and can overwrite the current paper.

## Data and attribution

Primary sources include SPARC (Lelli, McGaugh & Schombert), DustPedia/CIGALE (Nersesian and collaborators), Pantheon+SH0ES, Cosmicflows-4, X-COP, DES supernova duration products, and COBE/FIRAS and Planck products. Source URLs, citations, original catalog ReadMe files and existing license notices accompany the relevant checkpoints. Cite the original datasets and measurements when using them; this project does not claim ownership or grant a new license over third-party material.

See each report for exact sample cuts, units, fixed versus fitted quantities, uncertainty limitations, and source references. Start with the current checkpoint for later findings and read the recovered causal report for the original baseline interpretation.
