# Positional association audit before gravity inference

Four extreme-speed training examples have suspect APOGEE–Gaia associations. Checking the actual 2MASS observation dates does not remove their position disagreements. This changes the next step: resolve associations before combining Gaia parallaxes with StarHorse distances, and before treating the high inferred speeds as evidence for additional gravity.

## Measurements and result

The original APOGEE positions and J/H/K photometry agree with the queried 2MASS entries. APOGEE's packaged Gaia identifiers and parallaxes agree with the Gaia entries used in our matched catalog. Thus this diagnostic does not reveal an accidental identifier swap introduced by our parquet join. It does reveal a positional problem in those inherited associations; it does not establish its exact upstream cause.

| APOGEE target | 2MASS epoch, year | Original Gaia separation after epoch correction, arcsec | Follow-up |
|---|---:|---:|---|
| 2M11112524-6128593 | 2000.06 | 2.548 | Another proper-motion source is 0.782 arcsec away after correction; association remains ambiguous. |
| 2M17452907-3227078 | 1998.63 | 2.252 | Different source 4054069937711518208 is 0.045 arcsec away after correction. Strong candidate for review, not automatically substituted. |
| 2M17353316-3239126 | 1998.63 | 2.858 | No close alternative with measured proper motion in the queried five-arcsecond neighborhood. |
| 2M16453225-4452301 | 1999.38 | 2.207 | A different position-only Gaia detection lies close to the target; no proper motion or parallax is supplied, so it cannot supply a replacement full velocity. |

The 2MASS position-error ellipse major axes are 0.13 arcsec for the first target and 0.06 arcsec for the others. These comparisons are diagnostic, not formal association probabilities or independent significance estimates. Crowding, blending, source multiplicity, orbital motion and uncertainty calibration need explicit treatment. The nearby candidates' parallaxes must not be chosen because they happen to agree better with StarHorse.

## Formula provenance and computation

**Known astrometric approximation, not new theory:** propagate the Gaia sky position using its measured angular velocity:

`angular displacement = proper motion × (2MASS observation year − 2016.0)`.

The script applies the displacement on the sphere. Gaia's reference epoch is 2016.0; the 2MASS `jdate` field supplies the actual observation date. Coordinate equinox J2000 is not assumed to be the observation epoch. Perspective acceleration, annual parallax and individual orbital motion are not fitted here. These are positional checks, not a new distance posterior.

**Known least-squares geometry:** for tangent-plane offset vector d and proper-motion vector m, the best unrestricted elapsed time is `t = (d·m)/(m·m)`, leaving residual `|d − t m|`. This asks whether any constant-motion epoch could reconcile the positions. It deliberately permits unrealistic epochs, so a small residual does not certify an association. Among 77,927 training stars, 190 have residual >0.3 arcsec, 109 >0.5, 79 >1, and 46 >1.5. These are review counts, not estimates of the full mismatch rate. This screen can miss offsets aligned with proper motion and incorrect sources close on the sky.

## Consequences for the bulge test

1. Keep these four original associations in the historical data, but label them unsuitable for inference pending association review. Their provisional orbit launches are not validated physical states.
2. Review training positional flags before any additional parallax likelihood. A high Gaia fidelity score describes the Gaia solution; it does not prove that Gaia observed the same star as APOGEE.
3. For credible replacements, retrieve full astrometric covariance and source quality, resolve duplicate mappings, and reassess which data StarHorse used. Preserve the original APOGEE identity and frozen evaluation role through any repair.
4. For unresolved or position-only matches, use an explicit missing-motion/association treatment or document exclusion and its spatial selection effect. Do not fabricate transverse velocities.
5. Only then fit comparable bulge-plane, off-plane and outer-disk populations under the shared potential/selection likelihood. These four cases do not explain all 1,070 previously flagged distance disagreements.

No original IDs, distances, frozen roles, gravity parameters or held-out scores were changed. The original matched-parent SHA256 remains `a60182c30b98fe06d7e928a3aed5f70fa59342bdc44dc65b9f8005ff0b56f959`. No claim of a successful gravity fit follows from this audit. Timing, photon–companion dynamics, capture, lensing and the deferred source-energy budget remain separate unfinished requirements.

## Reproduction and sources

Run `python research_work/results/stellar-orbit-support/association_audit.py` with the existing parent and role caches. The small downloaded query snapshots are archived in `association-inputs/`; copy them into the ignored stellar-orbit-support cache if necessary. Full results, integer-safe source IDs, candidate lists and input hashes are in `association-audit.json`. The broader training screen is an ignored parquet sidecar. No held-out velocity likelihood was evaluated.

- [ESA Gaia DR3 contents](https://www.cosmos.esa.int/web/gaia/dr3): reference epoch and catalog context.
- [2MASS point-source column documentation](https://irsa.ipac.caltech.edu/data/2MASS/docs/releases/allsky/doc/sec2_2a.html): positions, observation dates and uncertainties.
- [IRSA TAP](https://irsa.ipac.caltech.edu/TAP): four exact designation queries in `fp_psc`.
- [Gaia archive TAP](https://gea.esac.esa.int/tap-server/tap): DR3 sources within five arcseconds of those four positions, query archived verbatim. Neighbors outside this radius were not investigated.
- [APOGEE DR17 synopsis](https://www.sdss4.org/dr17/irspec/dr_synopsis/): source catalog context.

Queries retrieved 2026-09-10. This is a targeted follow-up of training anomalies selected after inspecting the training set, not a blind validation exercise.
