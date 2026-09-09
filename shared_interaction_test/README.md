# Reproduce the fixed-law interaction and brightness checks

Run from this directory with Python 3 and NumPy, SciPy and pandas installed:

```bash
OPENBLAS_NUM_THREADS=1 python brightness.py
python interaction.py
```

Both scripts use paths relative to their own location. The complete official brightness table and covariance are bundled, so no network access is needed. `data/manifest.json` pins source URLs and SHA-256 hashes; `data/commit.json` records the Git commit. `protocol.json` was written before the brightness outcomes were inspected. It is an internal analysis record, not an external preregistration.

`brightness_predictions.csv` records all 960 evaluation predictions. `brightness_results.json` includes calibrator results, full-covariance statistics, bin covariance, the constant-offset diagnostic and observer-frame sensitivity. `interaction_results.json` records the leading atomic-response and benchmark gravitational-wave calculations. The complete derivation and limitations are in `shared_interaction_report.md`.

The paper's new Appendix E summarizes this work. The report and scripts distinguish a proposed shared action from a proven microscopic theory, standardized real-data tests from raw measurements, and leading estimates from exact species-specific predictions. No coefficient, opacity, luminosity-evolution or radiation-backreaction parameter was fitted to the high-redshift brightness data.

Source papers used for reading are cited but not redistributed. Numerical covariances and released standardized photometry retain their original filenames and README.
