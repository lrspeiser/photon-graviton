# Observational extension of the temporal-transport paper

Run time_revision/analyze.py with Python, numpy, scipy and pandas. It reads the included derived data, exports event tables, and computes results.json. No network is required. DES pickles are read with a restricted allowlist of numpy/pandas reconstruction classes. The downloaded authors' scripts are provided as methodological provenance, not executed by this analysis; their GPL-3.0 license is included.

Sources (retrieved 2026-09-08):
- NASA FIRAS residual spectrum: https://lambda.gsfc.nasa.gov/data/cobe/firas/monopole_spec/firas_monopole_spec_v1.txt
- NASA documentation: https://lambda.gsfc.nasa.gov/product/cobe/firas_monopole_spect.html
- DES per-band widths, authors' plotting and width code: https://github.com/ryanwhite1/DES-Time-Dilation (main snapshot; hashes of downloaded inputs recorded in results.json).
- DES paper: https://arxiv.org/html/2406.05050v2
- Distant temperature article XML: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8810383/fullTextXML . Only the previously published z=6.34 temperature interval is compared; its image-based compilation was not converted into an independent dataset or fit.

EVIDENCE STATUS
DES: 1504 event width measurements are derived from reference curves already corrected for time dilation. Reproducing b=1.002877 +/-0.004825 is a consistency reproduction, not raw-photometry proof. Event bootstrap omits cross-event shared-reference covariance and systematic uncertainty. Per-band event overlap prevents combining these as independent results.
FIRAS: fit the released residuals, not the rounded reconstructed absolute spectrum. A temperature-perturbation, y-distortion, and Galaxy-template nuisance fit uses diagonal errors only. The y-template is a restricted distortion family; absence of this distortion does not establish all possible spectra fit. Coherent thermal-preserving transport predicts the same present Planck curve as many alternatives, so this does not fit kappa or show a cosmic origin. Unweighted RMS across all released residual channels is noise-dominated at the high-frequency end and must not be confused with the published peak-normalized precision statement.
Screened ray calculation: synthetic consistency test using n=1+a*y*sin(pi*u)^2 in dimensionless coordinates. It is NOT observed data. Smooth adiabatic geometric optics, unchanged endpoint matter, prescribed field, no scattering, and a sufficiently short signal are assumptions. Field backreaction and microscopic matter coupling are not solved.
CMB origin, geometry/topology, directional alignment, environmental density and rotation: not fitted. The proposed families do not yet specify independent predictions, initial conditions, collision terms, or a spatial field equation. A free aligned harmonic fit would reproduce a selected target without identifying a cause. No Planck sky-map likelihood is claimed.

PAPER REPRODUCTION
The archive contains the original paper at redshift_paper/temporal_redshift_paper.docx. Run time_revision/update_paper.py once against that original to create the revised paper (python-docx required). It adds Sections 9--11, revises abstract and conclusion, and adds references. Repeated execution would duplicate the extension, so start from the archived original for reproduction.

The earlier temporal_hypothesis_checks folder contains the 164-group exploratory alternatives and summary constraints; its README documents its separate inputs and assumptions. No old held-out labels are newly blind. The original experimental paper archive remains separately available.

After updating the paper, run time_revision/layout_fix.py to apply the final pagination correction.
