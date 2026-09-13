# Cluster observation acquisition and target audit

## New evidence

The current Virgo/M87-centered forward calculation still has no verified same-target observational lensing comparison. A search found that [NGVS XXI](https://arxiv.org/abs/1705.04329) measures stacked background clusters at 0.2<z<0.5 in the survey fields. Its name does not make it a measurement of Virgo's own shear. Those mass estimates are not eligible substitutes for observations of our 1 Mpc M87 receiver. This search is not proof that no Virgo measurement exists.

A separate candidate is Coma: [Kubo et al.](https://arxiv.org/abs/0709.0506) measured weak-lensing shear with SDSS, and [Okabe et al.](https://arxiv.org/abs/1304.2399) studied its spatial lensing substructure. Published halo masses from those analyses are model inferences, not the raw observable to fit. Their lensing measurements motivate acquiring Coma data, not transferring Virgo predictions to Coma unchanged.

We acquired [Kang et al. (2025), CDS J/ApJS/278/51](https://cdsarc.cds.unistra.fr/ftp/J/ApJS/278/51/ReadMe), a spectroscopy catalog intended for comparison with Coma lensing structure. Its wider table contains positions, corrected apparent r-band magnitudes, redshifts, redshift uncertainties, source/reference flags and author-assigned membership. The narrower weak-lensing-field counts discussed in the abstract refer to a subset, not the complete table.

## Verified download

| Check | Result |
|---|---:|
| Total table rows | 37,758 |
| Unique sequence numbers and SDSS IDs | 37,758 each |
| Extended-source flags | 26,895 |
| Point-source flags | 10,863 |
| Author-assigned Coma member flags | 1,826 |
| Redshifts attributed to this work | 2,194 |
| Rows reporting zero redshift uncertainty | 5 |

The row and flag counts agree with the byte-by-byte catalog documentation. Numeric columns parse with finite values, uncertainties are nonnegative, and positional bounds are recorded. This validates ingestion and schema, not the authors' membership decisions or the accuracy of every measurement. The five zero uncertainties must not become infinite statistical weights or exact velocities.

Source URLs, sizes and SHA-256 hashes of the decompressed files are in acquisition.json. The public table is served as table2.dat.gz; acquisition decompresses it explicitly. Raw files are cached under research_work/generated/coma-kang2025 and excluded from git; the acquisition script and manifest are committed. Cached files are reused and rehashed on rerun. This is an exploratory acquisition, not a fresh held-out test set. Existing final holdouts remain unopened.

## What this enables and what it does not

The data supply a concrete starting sample for source photometry and line-of-sight kinematics around Coma lensing structures. They do not include the individual background shape measurements, calibration weights, source-distance information, image masks or shear covariance required for a lensing likelihood. A table of inferred subhalo masses cannot replace those inputs.

Before a joint test:

1. Obtain calibrated background shapes or a documented binned shear/reduced-shear product with covariance, centers, angular bins and selection. Preserve angular observables before importing any cosmological conversion.
2. Specify Coma distances and source geometry consistently with the hypothetical model. Measured redshift may contain both propagation and motion contributions; it is not automatically distance or velocity. Author membership flags need an audit rather than automatic adoption.
3. Construct Coma's own ordinary-matter and emitter models, including flux-to-luminosity conversion, bolometric corrections, source histories and light outside the surveyed region. Apparent r-band magnitude alone is not bolometric energy supply.
4. Predict the velocity distribution and reduced shear from the same supported deposit field. Individual galaxy velocities are not circular velocities, and reduced shear is not a deflection angle or a published halo mass.
5. Freeze a comparison protocol before fitting, with training/validation separated and observational/selection errors included. The presently acquired table is already exploratory and must not later be described as untouched evidence.

The Virgo calculation remains a useful forward pilot; Coma is an additional observational candidate, not a replacement of the six objectives. No mass, age, response gain, membership cut or observational fit was optimized here. All six objectives remain open.

Run `python research_work/results/cluster-observation-readiness/acquire.py`.
