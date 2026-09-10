# Cepheid measurements: acquisition and readiness

The full Gaia DR3 classical-Cepheid parent join has been downloaded: **13,077 stars, 57 columns**. A separate archive count agrees exactly; source identifiers are unique integers. The query has no TOP limit, sky filter, or measurement-quality cut. This provides individual measurements for a future common-reference-frame comparison with the existing red-giant sample.

## What is available

| Measurement check, applied cumulatively | Remaining stars |
|---|---:|
| Classical-Cepheid parent | 13,077 |
| Sky coordinates and proper motions | 13,050 |
| Intensity-averaged G, BP and RP photometry | 12,956 |
| At least one positive reported pulsation period | 12,956 |
| Combined radial velocity and positive reported error | 2,632 |
| RUWE below 1.4 | 2,436 |
| At least eight combined-RV transits | 2,304 |
| Complete astrometric covariance fields | 2,304 |

These are measurement-availability flags, not a final science selection. All 13,077 rows remain in the derived local parquet. Covariance completeness does not establish a positive-semidefinite matrix. Modes comprise 7,712 fundamental, 5,110 first-overtone and 255 multimode stars; they must not share an invented period-distance calibration.

Independently of the cumulative cuts, 2,633 stars have combined radial velocities and 718 have pulsation-model mean velocities with positive errors. Both measurements are retained. The combined-RV transit count and variable-star clean-epoch count are different fields, and Cepheid pulsation can affect the interpretation of measured velocities.

There are **zero exact Gaia source-ID matches** between this parent and the prepared 140,407-star red-giant sample (including its 128,772 quality candidates). This establishes distinct catalog objects, not independent Gaia calibration, distance assumptions, or Galaxy-wide systematics.

## What this does and does not establish

This is not a reconstruction of Feng et al.'s 903-star sample. No calibrated distances, Magellanic/external-galaxy exclusions, Galactic radius/height cuts, or velocity cuts have been applied. Consequently, 2,304 must not be presented as a contradiction of their selected count. No gravity model was fitted, and no published-curve or holdout score changed.

The next analysis requires:

1. A documented period-luminosity/Wesenheit calibration for the appropriate pulsation modes, with uncertainty and extinction assumptions. Distances turn an angular movement into a physical speed; a wrong distance can imitate a gravity discrepancy.
2. A common solar position and velocity convention for both populations, propagating correlated astrometric errors. Otherwise different coordinate choices can create part of the apparent difference.
3. Declared spatial, quality and population selections, including pulsation-velocity treatment. Different populations need not have the same average orbit even in the same gravitational well.
4. An orbital-population likelihood and survey selection treatment in the three-dimensional ordinary-matter and candidate fields. A potential predicts acceleration; velocity distributions also depend on the orbits occupying it.
5. A documented validation split before selecting model changes using individual-star outcomes. Previously examined rotation curves remain exposed diagnostics.

The photon-to-companion production, permanent storage and shared gravitational-response mechanisms remain unresolved. Acquiring these data does not validate those mechanisms or settle the deferred total photon-supply budget. This report introduces no physics formula; its checks use standard catalog processing.

## Reproduction and provenance

Run `fetch.py` to retrieve or reuse the cached VOTable, then `audit.py` to produce the derived parquet, counts and overlap list. `fetch.py` resumes its recorded archive job rather than launching a duplicate. Raw data remain in the ignored `research_work/data-cache/cepheid-stars/` directory. Scripts, query and manifests are versioned.

The raw VOTable SHA-256 is `c19fb1d53e4977daca4814ccc916ec47edfc451b8d868ec9896f084cdc5d597c`. `download.json` records the query and archive job; `audit.json` records the derived-file hash and exact counts. The raw SOURCE_ID capitalization is preserved in the VOTable and normalized only in the derived table.

- [Gaia DR3 documentation](https://gea.esac.esa.int/archive/documentation/GDR3/)
- [Gaia Cepheid table metadata](https://gaia.aip.de/metadata/gaiadr3/vari_cepheid/)
- [Gaia DR3 Cepheid release paper](https://www.aanda.org/articles/aa/pdf/2023/06/aa43990-22.pdf)
- [Feng et al. rotation-curve paper](https://academic.oup.com/mnras/article/546/2/stag011/8416425)

For the already-local APOGEE, StarHorse, Gaia subset and BRAVA inventory, see [the verified catalog recheck](../stellar-data-readiness/local-recheck.md).
