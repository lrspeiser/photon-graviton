# Milky Way stellar-data readiness

**Input-provenance correction:** the [later training audit](../stellar-orbit-support/report.md) recovers StarHorse_INPUTFLAGS, the wider distance quantiles and public astrometric-fidelity scores. The preparation below had retained output flags but omitted input flags. StarHorse did not record parallax use for 19,381 of the 77,927 subsequently allocated training stars, so its use of Gaia cannot be assumed for every row. A new sidecar records the required distance/astrometry treatment and disagreements without altering this frozen parent catalog. Earlier exploratory motions and uncertainty reconstructions remain provisional; they are not validated gravity measurements.

We recovered useful observations locally, downloaded the missing official APOGEE and StarHorse catalogs, and prepared an exploratory sample of 140,407 giant stars. This establishes actual data coverage for the proposed bulge/plane comparison. It does not yet measure a companion-induced gravitational field.

**Gaia enrichment is complete:** all 140,407 requested sources were returned with astrometric correlations and quality fields. The additional stated Gaia cuts retain **128,772 candidate stars**. The highest absolute difference between Gaia's PMRA and the rounded value packaged in APOGEE is 0.0000064 mas/yr. The enriched table contains both catalogs' values and a `gaia_quality_candidate` column; it retains excluded rows for auditing. It does not silently change the original 442-cell exploratory summary.

## What is on the machine

| Product | Rows | Status |
|---|---:|---|
| Recovered BRAVA catalog | 8,585 | Original IRSA-format table; retained as a separate line-of-sight check |
| Recovered BRAVA–Gaia table | 6,189 | Existing cross-match; missing full astrometric covariance, not automatically accepted as a validated reduction |
| Recovered GIBS catalog | 5,651 | Additional bulge radial velocities and errors; provenance requires review before combining |
| Recovered Gaia sky wedges | 144,000 | Twelve files of 12,000 rows each; selected subsets, not a complete all-sky census |
| Official StarHorse APOGEE v2 | 562,424 | Downloaded and parsed; distances and uncertainty quantiles |
| Official APOGEE DR17 allStarLite, synspec_rev1 | 733,901 | Downloaded and parsed; combined-spectrum rows, not this many unique stars |

The recovered project is `C:/Users/henry/Documents/Codex/2026-08-09/sigma-theory-compiler-publish`. Its earlier velocity corrections, inferred circular speeds and theory-specific processed tables are not imported into the new sample. The inventory includes several older APOGEE/Eilers cross-matches; those are unnecessary for the new official-catalog join. The search covered the relevant Documents/Codex and Downloads locations, with an initial attachment search; it was not an exhaustive scan of every disk.

The [APOGEE access documentation](https://www.sdss4.org/dr17/irspec/spectro_data/) identifies the revised SAS catalog and notes that the CAS version differs. We use SAS synspec_rev1. [StarHorse's release page](https://data.aip.de/projects/aqueiroz2023.html), DOI 10.17876/data/2023_1, identifies the v2 product and the correction to age/mass priors. This preparation uses distance quantiles, not giant-star ages. Gaia astrometry in EDR3 and DR3 is shared, as documented by [ESA](https://www.cosmos.esa.int/web/gaia/dr3); the EDR3 astrometry packaged in APOGEE is therefore useful without substituting Gaia DR2 radial velocities for APOGEE velocities. [BRAVA's archive](https://irsa.ipac.caltech.edu/data/BRAVA/) remains the reference for the recovered independent bulge sample.

## Exact matching and initial quality cuts

We match on both APOGEE_ID and the integer Gaia source identifier. Identifiers never pass through floating point. Within repeated APOGEE pairs, the highest-SNR spectrum is selected, with field/telescope tie-breaks. No velocities are used to select the best spectrum. The join yields 554,738 pairs; 28 rows share 14 Gaia identifiers across different APOGEE IDs. Keeping one spectrum per source leaves **554,724 distinct matched stars**. Ambiguous duplicate StarHorse pairs would be excluded; none occurred. We do not fill unmatched distances with invented values.

The exploratory cuts require ordered positive distance quantiles, relative distance half-width at most 20%, SNR at least 70, zero STARFLAG and ASPCAPFLAG, blank StarHorse output warnings, giant-star surface gravity between 0 and 3.5, valid radial velocity with positive error at most 5 km/s, visit scatter at most 2 km/s, finite proper motions with positive errors at most 1 mas/yr, and finite plausible iron/alpha abundances. These conservative cuts define a selected sample; they are not a survey-completeness correction. A small measured visit scatter does not prove a star is single.

The result is **140,407 stars**, with the full sequential count audit in `cutflow.json`. A blank FITS warning string is masked by Astropy and becomes a missing value in pandas; it is explicitly restored to blank before interpreting flags.

## Can we examine the places the user proposed?

| Region | Cylindrical radius R | Absolute height above plane | Stars before further Gaia-quality cuts |
|---|---|---|---:|
| Plane beneath the bulge | 0.5–3.5 kpc | below 0.2 kpc | 758 |
| Above/below the bulge | 0.5–3.5 kpc | 0.5–1.5 kpc | 4,265 |
| Disk plane control | 5–9 kpc | below 0.2 kpc | 8,508 |
| Disk off-plane control | 5–9 kpc | 0.5–1.5 kpc | 30,988 |

After Gaia quality cuts, these counts become **717**, **4,085**, **7,703**, and **27,620**, respectively. Those cuts require RUWE below 1.4, at least nine visibility periods, a five- or six-parameter astrometric solution, no duplicated-source flag and no Gaia non-single-star flag. Non-detection of a binary is not proof that the star is single.

These are coverage counts, not measurements of extra pull. R is distance from the rotation axis, not distance from Earth. The bulge labels identify comparison regions, not confirmed membership of every star. Foreground/background contamination and distance uncertainty can move stars between regions.

The prepared table uses distance medians to calculate positions and all three velocity components. Adopted frame settings are centre distance 8.2 kpc, solar height 0.0208 kpc, solar velocity (11.1, 248, 7.25) km/s in Astropy's Galactocentric convention, and an illustrative bar angle of 27 degrees. These are declared working choices, not newly measured constants. They must be varied in the dynamical analysis.

The initial transformation uses APOGEE sky coordinates and packaged Gaia proper motions. The enriched table also supplies Gaia coordinates. The final likelihood should adopt a common position/velocity epoch and account for the observation epochs; the exploratory transformation is not an epoch-propagated orbital solution.

For inspection, `exploratory-moments.json` contains 442 cells with at least 30 stars, separated by radius, height, hemisphere, metallicity, alpha abundance and bar azimuth. Means and standard deviations are observed sample summaries at median distance, without measurement-error subtraction or survey-selection correction. They must not be used as intrinsic dispersions with falsely precise errors. The standard transformations and sample statistics are known mathematics, not formulas unique to our hypothesis.

## Ordinary-matter model audit

We inspected and pinned the [official AGAMA bar example](https://github.com/GalacticDynamics-Oxford/Agama/blob/f302756b8af2b763db58e278e30478517dc8eea3/py/example_mw_bar_potential.py). Its complete model includes an inner X-shaped bar, two long-bar components, a disk, a central mass concentration, and a separate halo. Loading its complete returned potential would include dark matter. Adding another full disk would double-count ordinary matter. The example also warns that the original model is unsuitable as a realistic description beyond about 5 kpc. Code was inspected, not executed.

The [Sormani et al. paper](https://arxiv.org/abs/2204.13114) explains its dependence on a Portail model fitted to density and kinematic observations. Its normalization is therefore not independent evidence for our proposed extra gravity. Use its bar shape as a candidate while propagating mass-model uncertainty. Construct a consistent disk/bar transition, gas distribution and central components before fitting the outer disk. A standard halo may be included as an explicitly labeled comparison, as requested, but is not an active premise of the companion hypothesis.

## What remains, why, and the consequence of skipping it

1. **Use the downloaded Gaia covariance in the uncertainty calculation.** `fetch_gaia.py` obtained all 140,407 requested sources; `gaia-download.json` records the completed file hash. Correlations describe which measurement errors move together. Ignoring them can turn distance/proper-motion errors into an apparent physical pattern. Acquisition is complete; propagating these errors into a fitted model remains necessary.
2. **Propagate distances and positions consistently.** StarHorse already uses Gaia parallax. Treating its posterior as an independent distance measurement and multiplying by the same parallax likelihood would count information twice. Marginal distance quantiles alone are not a complete joint posterior with proper motion.
3. **Match populations and survey selection.** Compare like-for-like chemistry, bar position and detection probability. Otherwise different orbit populations or obscuration can masquerade as a difference caused by deposits. Keep above/below comparisons separate and assess their overlap in coverage.
4. **Build a bar-aware orbital prediction for each candidate potential.** A potential gives accelerations; an orbital population gives velocity distributions. Circular-speed formulas cannot predict every bulge star's speed. Include bar rotation, nuisance parameters, and tracer density/velocity structure rather than fitting an arbitrary dispersion conversion factor.
5. **Specify how deposits change the scalar potential.** The prior local-well hypothesis allows deepening without an arrival-directed force. Its kernel/coupling and deposition distribution must predict one shared spatial field. Freely adjusting a correction in every observed cell would reproduce a table without testing a mechanism.
6. **Reserve complete fields or spatial/population regions for evaluation.** Neighboring stars share selection and model errors; a random star split alone overstates independence. These exploratory summaries have now been inspected and must not be advertised as a pristine blind test.
7. **Return to the joint constraints.** A successful stellar fit still needs consistent lensing, photon redshift and timing, companion transport, and conservation of energy. None is certified by the new catalog preparation.

## Reproduction and storage

Run from the repository root:

```powershell
python research_work/results/stellar-data-readiness/prepare.py --download --local-project C:/Users/henry/Documents/Codex/2026-08-09/sigma-theory-compiler-publish
python research_work/results/stellar-data-readiness/fetch_gaia.py
python research_work/results/stellar-data-readiness/verify.py
```

Raw catalogs and the prepared `matched-exploratory.parquet` are under `research_work/data-cache/stellar-catalogs/`, excluded from Git to keep the repository small. Tracked scripts, exact download URLs, SHA-256 hashes, schemas, cut counts and summary results preserve the audit trail. The Gaia query is a join on uploaded public source identifiers with no arbitrary TOP truncation. `bar-model-audit.json` pins the inspected model source. This is data preparation and model auditing; the three-potential bulge fit has not yet been performed.

The enriched deliverable is `matched-with-gaia-covariance.parquet`. `verify.py` checks unique integer identifiers, complete one-to-one Gaia enrichment and a round trip from derived cylindrical positions/velocities back to the original observables. Calculation inputs are promoted to double precision to avoid unnecessary rounding in coordinate transformations. Catalog measurement precision is not increased by that promotion.
