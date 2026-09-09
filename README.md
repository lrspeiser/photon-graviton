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

- [Collective spatial response](research_work/results/collective-response/derivation.md) narrows the scalar candidate's scattering and weakens its rate dependence from E^6 to E^3 at large target size. Energy transfers remain about 25% per event; this is not yet an achromatic narrow-line redshift mechanism.

- [Matter-assisted conversion](research_work/results/matter-assisted/derivation.md) now has an explicit scalar/polarizability action, recoil accounting and a leading cross section. The simple contact limit yields strongly color-dependent, broad-angle scattering rather than the desired narrow-line redshift. This is a constrained comparison candidate, not an adopted theory.

- The resumed [interaction-rate study](research_work/results/interaction-rate/derivation.md) derives conditions for the same fractional shift across colors, checks why one scalar-photon vacuum operator does not supply partial redshift under standard dispersion, and shows how reverse transfers can increase spectral noise. These candidate-specific constraints do not rule out the overall conversion-first program.

- **Latest conversion-first pass:** a constant photon-energy loss law was fitted to 164 previously exposed galaxy groups with fixed published distances. It gives alpha about 7.63e-5 per million light-years; the nearby data do not distinguish it from a linear distance trend. A resonant inelastic mode toy conserves energy without special time stretching, but does not derive the physical rate, smooth redshift, capture or gravity. Read [the executed pass](research_work/results/conversion-first/report.md) and [revised goal](research_plan/active-goal.md).

- Changing the propagation history improves the earlier redshift–brightness fit, but does not establish an interaction or protect atomic/cavity frequency ratios.
- Using measured total luminosities for 26 DustPedia–SPARC overlap galaxies, converting all of their present luminosity for 10 billion years still falls short of the ordinary deposited mass-energy requirement by a median factor of about **5,400**. At the fixed local redshift rate, the shortfall is much greater.
- Catalogued neighboring sources do not solve the energy budget. A constant-bath control performs as well as or better than actual neighboring-source illumination under the primary log-residual metric.
- A depth-dependent empirical extra-acceleration law improves galaxy validation and approximately matches eleven cluster excess-mass estimates without fitting a cluster normalization. Its uncertainty is broad and **it is not yet derived from the photon energy supply**.
- Repeated light circulation can increase encounters, but cannot generate unlimited energy. Extrapolating the current positive-exponent propagation history gives a finite remaining path length even in static space; that extrapolation is a conditional prediction, not a measured horizon.

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

The current command checks original source hashes, runs diagnostics, and verifies that original files and saved results remain unchanged. Each invocation creates a fresh directory under `research_work/generated/`, with logs and a verification summary. `--output-dir` selects another new or empty directory. The latest [collective-response verification](research_work/collective-response-verification.json) records 15 current default jobs; the optional historical baseline makes 16 jobs and was not rerun for this update. Earlier [conversion-first verification](research_work/conversion-first-verification.json) and [publication verification](research_work/publication-verification.json) retain their original suite scope and results.

The historical scripts can still run directly, but doing so overwrites their corresponding outputs. Use the isolated runner above for baseline reproduction. `python restore_data.py` restores the compressed covariance if needed by older brightness analyses; the current diagnostics read its compressed form directly. These are research tools, not a complete cosmology library.

Earlier checkpoints may depend on additional libraries, archived inputs, network services, or document-rendering utilities. Their environment records and reports are preserved; not every historical script has been rerun for this commit. Paper-building scripts include original workspace paths and need adjustment for a different machine. Running all old paper-update scripts in sequence is not a supported build workflow and can overwrite the current paper.

## Data and attribution

Primary sources include SPARC (Lelli, McGaugh & Schombert), DustPedia/CIGALE (Nersesian and collaborators), Pantheon+SH0ES, Cosmicflows-4, X-COP, DES supernova duration products, and COBE/FIRAS and Planck products. Source URLs, citations, original catalog ReadMe files and existing license notices accompany the relevant checkpoints. Cite the original datasets and measurements when using them; this project does not claim ownership or grant a new license over third-party material.

See each report for exact sample cuts, units, fixed versus fitted quantities, uncertainty limitations, and source references. Start with the current checkpoint for later findings and read the recovered causal report for the original baseline interpretation.
