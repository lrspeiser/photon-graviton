# Measured lensing inputs acquired: SLACS

The public SLACS release is now archived locally: 131 lens candidates and 63 modeled grade-A strong lenses. Fifty-eight of those 63 have reported spectroscopic stellar velocity dispersions and uncertainties. This supplies an actual observational target for the next joint motion-and-lensing test, beyond the previous hypothetical deflection tables.

The sample is predominantly early-type galaxies: 57 early-type, five late-type and one unclassified. It therefore tests transfer to a different population; it is not a ready-made matched set of the 149 SPARC disk galaxies. Its strong-lens selection must be modeled before making population-wide claims.

## Keep observations and inferences distinct

| Archived quantity | Status and treatment |
|---|---|
| Foreground/background spectroscopic redshift | Measured redshift; not itself a distance or proof of expansion. |
| Stellar sigma and its quoted error | Spectroscopic velocity dispersion, uncorrected for aperture effects. Not a circular speed and not a lens-model dispersion. Five unavailable values remain missing. |
| Apparent I-band magnitude | HST image-model photometry; Galactic extinction is supplied separately. No silent luminosity conversion. |
| Angular effective radius and light axis ratio | De Vaucouleurs image-model summaries, not physical radii or an assumption-free luminosity profile. |
| SIE Einstein angle, mass axis ratio and position angle | Inferred from imaging using a singular isothermal ellipsoid model. Einstein angle uses the catalog's intermediate-axis normalization. Not a raw image pixel or an independent total-mass measurement. |
| LTM Einstein angle and external shear | Alternative light-traces-mass image model. Keep as a model-dependence comparison, not an independent observation of the same system. |
| Catalog rest-frame luminosity and enclosed luminosity fractions | Preserved in raw source tables, excluded from the analysis-ready field subset. Luminosity incorporates cosmological distance modulus, evolution and k-corrections; fractions depend on image/lens models. |

The catalog states a minimum quoted sigma error of five percent and withholds dispersions for insufficient spectral S/N or multiple foreground systems. We preserve its availability and quality flags rather than filling missing values with quantities calculated from lensing.

These distinctions prevent a circular test: deriving a “velocity dispersion” from the Einstein angle and then claiming that the same angle agrees with stellar motions would not provide independent evidence.

## Frozen system assignments

Before any lensing-fit score is computed, system IDs are assigned by a documented SHA256 rule:

| Role | All modeled lenses | With reported stellar dispersion |
|---|---:|---:|
| Training | 41 | 37 |
| Validation | 9 | 8 |
| Test | 13 | 13 |

The salt, exact rule and per-system assignments are archived. These are future evaluation roles, not a claim of pristine blindness: the catalog is public, its basic fields and aggregate summaries have been inspected, and the systems were selected as lenses. No model was fitted or ranked against the Einstein angles or stellar dispersions in this acquisition step. Previous stellar and SPARC roles are unchanged.

## Geometry under the current redshift proposal

**Known inversion of the project's empirical exponential redshift relation, plus an additional static-Euclidean geometry assumption:**

\[
D(z)=\frac{\ln(1+z)}{\alpha},\qquad
\frac{D_{ls}}{D_s}=1-\frac{\ln(1+z_l)}{\ln(1+z_s)}.
\]

The distance rule is not claimed as a new mathematical formula or a physical derivation of photon conversion. The per-system geometry table uses the exact fitted alpha stored in the joint-galaxy audit, 0.0002488993265191759 per Mpc. The nearby rounded/fixed reference alpha used elsewhere is not silently substituted. Foreground redshifts span 0.0629–0.5132 and source redshifts 0.1965–1.1924.

The table contains **conditional model distances**, not independently measured galaxy distances. It omits peculiar motion and endpoint redshift corrections and extends the empirical relation beyond its nearby calibration domain. Treating those distances as angular distances also assumes the stated static Euclidean geometry. A different propagation metric may change that identification and must be derived consistently.

An algebraic consequence is that alpha cancels from D_ls/D_s in this geometry. Merely retuning the overall photon-loss rate cannot adjust this distance ratio. Physical lens radii still depend on D_l, so alpha does not disappear from every mass or light-profile prediction.

No expanding-universe distance relation, inferred dark-halo mass or cosmologically corrected catalog luminosity is adopted as a fact in the fictional model.

## What is needed before a meaningful comparison

1. **Specify the companion source and both metric potentials.** Stellar motion constrains the temporal potential; light bending depends on the sum of temporal and spatial potentials. The previous equal-potential calculation was conditional. Without this step, an independently tuned lensing multiplier would hide a missing theory prediction.
2. **Construct each ordinary-matter light/mass model under declared distance and stellar-population assumptions.** Recover the relevant imaging/photometric corrections and their errors. An early-type lens cannot be assigned a SPARC disk template or a dark-halo mass from its discovery paper.
3. **Predict the measured aperture dispersion.** Include the luminosity-weighted aperture, seeing, orbital anisotropy and potential. A circular velocity cannot simply be compared with the tabulated stellar dispersion.
4. **Use an imaging likelihood or explicitly limited image-model summary likelihood.** The release does not provide a full Einstein-angle/shape/shear covariance in table5. SIE and LTM values cannot be treated as two independent measurements or their difference as a calibrated one-sigma error.
5. **Fit only the training systems, select procedure on validation, and then score reserved systems once frozen.** Account for selection and ordinary-matter uncertainty, and retain failures. Do not label a successful isothermal benchmark as a derived photon-companion result.

Why these steps matter: otherwise the apparent agreement could come from distances, masses, orbital behavior or lensing strength that were chosen using the same observations we claim to predict. The acquired data now make those requirements concrete. Acquisition is complete for these summary tables; the physical model and joint likelihood are not.

## Files, verification and source

Source tables and the publisher/CDS byte descriptions are archived under source-tables, including notes about unmodeled systems. SHA256 hashes are recorded in manifest.json. prepare.py reads the CDS schema, verifies unique IDs, checks the 63-to-131 join, preserves missing dispersions, confirms foreground redshift is below source redshift and saves the analysis-ready JSON plus conditional geometry. The literal catalog unit string `----` on the name column triggers an Astropy unit warning; this is a dimensionless identifier and does not alter the numeric angle or velocity units.

Primary publication: Bolton et al., **The Sloan Lens ACS Survey. V. The Full ACS Strong-Lens Sample**, ApJ 682, 964 (2008), [DOI: 10.1086/589327](https://doi.org/10.1086/589327). Published electronic tables and column documentation: [CDS/VizieR J/ApJ/682/964](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJ/682/964). Retrieved 2026-09-10 from the release's CDS archive. The catalog's model assumptions are retained as provenance, not endorsed as the fictional universe's physical explanation.
