# Photon–Companion Research in a Hypothetical Universe

**This project explores a fictional universe that has the same observed stars, galaxies, and astronomical datasets as ours, while allowing different physical laws.** We use real measurements as a common test bed so that we can safely explore new concepts, derive their consequences, and compare their predictions without presenting speculation as established physics.

The working universe is nonexpanding. Candidate mechanisms let light redshift through changing propagation and exchange energy with hypothetical companion waves or fields. Captured energy or a modified field response might contribute to galaxy and cluster gravity. These are hypotheses to investigate, not discoveries about our universe. “Companion” does not mean an experimentally detected particle or an established graviton interaction.

The aim is to find a self-consistent emitter–propagation–receiver model that predicts color, duration, brightness, angular size, gravitational-wave propagation, galaxy rotation, and lensing with shared independently constrained parameters. No candidate in this repository has yet achieved all of those requirements.

## Start here

1. **[Latest analysis: causal supply, retention, gravity response, and circulation](companion_causal_test/report.md).** New DustPedia luminosity inputs, cross-validation, cluster transfer, controls, uncertainty checks, and the closed-universe question.
2. **[Main paper, version 9](redshift_paper/temporal_redshift_paper.docx).** The unified temporal–companion framework before the latest causal tests. Read it alongside the addendum above: the paper has not yet incorporated that addendum.
3. **[Companion deposition fits](companion_deposition_fit/deposition_fit_report.md)** and **[energy-budget tests](companion_wave_test/companion_wave_report.md)**.
4. **[Propagation-history fit](focused_history_test/report.md)** and **[clock/interaction derivation](minimal_clock_interaction/derivation.md)**.

## Current findings

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
| `unified_paper_v9`, `abstract_revision`, `paper_fixed_atom_update` | Paper construction, previous versions and revision scripts |

The original relative directory structure is preserved because analysis scripts share inputs. `SNAPSHOT_MANIFEST.json` inventories the committed research files and checksums. Regenerable page-render previews, caches, and duplicate delivery archives are omitted; scientific plots, source datasets, scripts, results, reports, and available paper versions are retained.

## Reproduce the latest analysis

Use Python 3.11 or later in a virtual environment:

```sh
python -m venv .venv
# Activate the environment using your operating system's usual command.
pip install -r requirements.txt
python companion_causal_test/run.py
python companion_causal_test/followup.py
python companion_causal_test/circulation.py
```

The latest numerical scripts use local datasets and NumPy/SciPy. The follow-up includes bootstrap refits. Results are written into `companion_causal_test/` and overwrite the corresponding output files. The scripts are research checkpoints, not a packaged cosmology library.

Earlier checkpoints may depend on additional libraries, archived inputs, network services, or document-rendering utilities. Their environment records and reports are preserved; not every historical script has been rerun for this commit. Paper-building scripts include original workspace paths and need adjustment for a different machine. Running all old paper-update scripts in sequence is not a supported build workflow and can overwrite the current paper.

## Data and attribution

Primary sources include SPARC (Lelli, McGaugh & Schombert), DustPedia/CIGALE (Nersesian and collaborators), Pantheon+SH0ES, Cosmicflows-4, X-COP, DES supernova duration products, and COBE/FIRAS and Planck products. Source URLs, citations, original catalog ReadMe files and existing license notices accompany the relevant checkpoints. Cite the original datasets and measurements when using them; this project does not claim ownership or grant a new license over third-party material.

See each report for exact sample cuts, units, fixed versus fitted quantities, uncertainty limitations, and source references. `companion_causal_test/report.md` is the authoritative interpretation of the newest results, including the tests that did not improve the hypothesis.
