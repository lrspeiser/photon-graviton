# Temporal candidate audit: reproduction

The main deliverable is `candidate_comparison_report.md`. The updated paper is included as `temporal_redshift_paper.docx` in the packaged release. Every branch is retained; no result establishes a complete temporal cosmology.

## Reproduce the calculations

Use Python with NumPy, SciPy and pandas. From this directory run:

```bash
python analyze.py
python uncertainty_checks.py
OPENBLAS_NUM_THREADS=1 python connected_field.py
```

No network access is needed: the input data are bundled. Results are written to `results/`. `prepare_inputs.py` documents the original session staging and official Planck downloads; it is not needed after extracting this package. `connected_field_grid_tied_packet.py` reproduces the earlier, grid-dependent packet-width branch under its own output filename. It should not be substituted for the fixed-width calculation.

## Evidence status

- Redshift: 164 preserved groups with training/validation/reserved counts 104/35/25. New parameters fitted on training groups only. Reserved objects were previously examined, so this is exploratory reuse, not a new blind test.
- SPARC: 149 eligible galaxies, split 89/29/31. New surface-density fits use only the 89 training galaxies. Existing power and RAR parameters are preserved. Published distances and mass-to-light ratios are not refitted.
- Supernovae: derived DES width averages and published spectral-aging measurements. Their template and shared-systematic assumptions are retained.
- FIRAS/Planck: released residuals or bandpowers, diagonal-error diagnostics rather than complete instrument likelihoods.
- GW, Galileo, atomic ratios and HFLS3: published numerical summaries, not raw-data reanalysis.
- Connected fields: synthetic dimensionless Hamiltonian experiments, not observational fits. Packet width, physical field energy, atomic clocks, gravity and the cosmic history are not established.

`input_manifest.json` identifies staged inputs and SHA-256 hashes. `results/verification.json` records split and numerical checks. `results/SPARC_uncertainty.json` and redshift result entries contain descriptive bootstrap intervals; they do not include the full research-selection history or every systematic uncertainty.

## Paper history

`paper_before_candidate_audit.docx` preserves the preceding paper. `update_paper.py` documents the appendix addition and intentionally refuses to append twice. Its original session path expects the paper in the sibling `redshift_paper/` directory; calculations above do not depend on that path or on the document script. The final paper is supplied directly in this release.

## Numerical correction recorded

For n_GW=exp(p*gamma*t), the derivative of GW arrival time at p=1 is

```text
dt_GW/dp = -[(1+x)*ln(1+x)-x]/gamma, x=gamma*R/c.
```

This was checked against a finite difference at p=1±0.0001. At the illustrative 40 Mpc distance, the relative agreement was 2.05e-8. The corresponding 10-second sensitivity is |p-1|≈4.83e-13. It is not a measured confidence bound.

Full citations and physical caveats are in the report. The archive's source tables remain the responsibility of their original authors; computed candidate predictions are distinct from those measurements.
