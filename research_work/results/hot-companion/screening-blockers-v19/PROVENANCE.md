# Provenance

This folder is the derivation package "One rate mechanism for screening and delayed release" supplied by the user on
25 September 2026 (the README, `run_checks.py`, `make_plots.py`, `results/` and `figures/` are as supplied, produced
outside this repository with Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0; `manifest.json` lists their SHA-256).

Checked here on 25 September 2026 (Python 3.11): `python run_checks.py --n-per-seed 20000` reproduces the full
master-equation check exactly (largest error in the released fraction 1.22 × 10⁻¹⁴, normalisation 1.79 × 10⁻¹⁴) and the
event simulation within sampling error (largest |z| 3.2 over 42 conditions at 100,000 realisations each; 2.17 at the
package's 1,000,000).

Carried to the Solar System by `code/blocker_release_v19.py` → `run-blocker-release-v19/` (README §29.5): the Sun's
companion is launched in the Sun's own field, so under this mechanism its release depends on its history.
