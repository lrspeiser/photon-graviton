# ESSENCE passband version and convention audit

The existing local SNDATA_ROOT archive contains two legacy ESSENCE response curves and a convention README. They have been recovered and hashed. They are **not adopted as the 2016 release calibration**.

## Why the distinction matters

The archive README states that its response includes an extra wavelength factor relative to a photon transmission. For photon-counting synthetic photometry, the familiar operation is an integral proportional to F_lambda times lambda times photon transmission. Thus one must divide the archived energy response by wavelength before using a routine that itself supplies the photon-counting factor. This conversion is known photometric bookkeeping, not a new physics formula.

Narayan et al. 2016 Appendix B describes the release photon throughput as the product of atmosphere, optics, filter and detector response. Table B2 lists selected numerical samples in print; the full electronic table is required for a release-matched calculation. Its full-system response must not be replaced by the filter glass transmission alone.

## Numerical comparison

After dividing the legacy curves by wavelength, fit one amplitude per band to four positive printed Table B2 samples. This tests shape compatibility while permitting arbitrary overall normalization. It is not a scientific fit or a reconstructed filter curve.

- R: residual at 7005 Angstrom is +2.06%; the other three selected points differ by roughly -0.8 to -1.0%.
- I: selected in-range samples differ by about -0.8%, +0.5% and -1.7%. At 9725 Angstrom the printed throughput is 0.0008, but the legacy array ends at 9500 Angstrom. The small tail is absent; its integrated photometric impact is not determined by this one point.
- The legacy arrays have 25 and 10 negative samples respectively. They are retained in the originals, not silently treated as physical negative transmission or clipped without documentation.

These discrepancies exceed simple rounding at some printed samples. A single normalization does not reproduce the published response. Eight sparse samples cannot recover the missing detailed curves or establish an error bound on synthetic supernova photometry.

## Access state and next step

The historical author host telescopes.rc.fas.harvard.edu failed DNS resolution. The publisher landing-page request did not provide usable table links. The published paper and the local legacy curves remain accessible. These are specific access failures, not an overall research impasse.

Next acquire Table B2 from the electronic supplement or an independently traceable copy and verify its printed anchor values and photon convention. Then integrate the same evolving source spectrum through R4m and I4m. The existing Hsiao template in the cache may support a conditional source-model pilot, but its training and temporal conventions must be audited before it is used to test time stretching.

No joint brightness prediction or new propagation fit is claimed here. All six research goals remain open.

## Reproducibility and sources

run.py recovers only the three specified archive members if the cache is missing and rejects changed cached inputs on subsequent runs. results.json contains hashes, all eight comparisons and negative-sample counts. No downloaded code is executed.

- [Narayan et al. 2016, Appendix B and Table B2](https://lss.fnal.gov/archive/2016/pub/fermilab-pub-16-402-ae.pdf)
- Local archive: research_work/generated/des-photometric-calibration/SNDATA_ROOT_2024-07-03.tar.gz; members filters/ESSENCE/README, CTIO4m_R.dat and CTIO4m_I.dat.
- [Same-object photometry](../sn2006mk-joint-photometry/report.md)
