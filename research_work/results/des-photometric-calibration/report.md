# DES brightness calibration resources and photon-count verification

The original DES five-year calibration inputs are now available locally, with full-archive checksum verification. This closes an acquisition gap for the joint light test; it is not a fit to supernova brightness or validation of companion transport.

## Acquired provenance

Source: [SNDATA_ROOT archive, DOI 10.5281/zenodo.12655677](https://zenodo.org/records/12655677), file `SNDATA_ROOT_2024-07-03.tar.gz`, 1,653,076,974 bytes. Its published MD5 b75e8fec76e53b9b9e51901a09226b46 was verified and a local SHA-256 recorded. The scripts extracted regular calibration/standard-spectrum files plus the four exact DECam filter dependencies and Hsiao spectral template named by the input. Every extracted file has a recorded SHA-256. No archive code was executed.

The archive and selected resources are in ignored `research_work/generated/des-photometric-calibration`; the acquisition scripts and manifests are committed. This preserves reproducibility without putting a 1.65 GB archive into Git.

The input declares AB magnitudes and COUNT response, and points to DES-SN3YR_DECam filter curves. Its historical calibration offsets are:

| Filter | Declared magnitude offset | Pivot wavelength (angstrom) |
| --- | ---: | ---: |
| DES-g | -0.0008 | 4809.40 |
| DES-r | -0.0116 | 6420.44 |
| DES-i | -0.0090 | 7816.03 |
| DES-z | +0.0071 | 9171.13 |

These are calibration-input offsets, not corrections we have already applied to observed flux. The sign/convention in the observation-to-model comparison must be traced before application. The archived input labels filters DES-g/r/i/z while the pinned photometry has g/r/i/z; that mapping must be made explicit.

## Independent count-rate check

**Known radiometry, not a new physical equation:** for a photon-counting band with transmission T(lambda), the rate per collecting area is

\[
C=\int T(\lambda) F_\lambda(\lambda)\frac{\lambda}{hc}\,d\lambda.
\]

For a flat-frequency AB reference, this is

\[
C_{\rm AB0}=\frac{10^{-0.4(48.6)}}{h}\int\frac{T(\lambda)}{\lambda}\,d\lambda
\]

in cgs units with a consistently chosen wavelength unit. We integrated the tabulated curves independently and compared to the original calibration program's retained log:

| Band | Our photon rate (s^-1 cm^-2) | Archived log | Relative difference |
| --- | ---: | ---: | ---: |
| g | 652180.42 | 652200 | 0.00300% |
| r | 745615.31 | 745700 | 0.01136% |
| i | 706898.36 | 706900 | 0.00023% |
| z | 554337.43 | 554400 | 0.01129% |

All satisfy the code's 0.1% numerical agreement check. Differences include numerical integration and rounded reference values; this is not a measurement of instrument calibration uncertainty. File hashes, wavelength ordering and nonnegative finite transmission were also verified.

## How this connects to the proposed propagation model

**Conditional consequence of our declared transport postulates:** for static dilution, spectral stretch S and event stretch S^b,

\[
F_{\lambda_o}(t_o)=\frac{L_{\lambda_e}(\lambda_o/S,t_o/S^b)}{4\pi D^2 S^{b+2}}.
\]

One factor 1/S accounts for photon energy, S^-b for arrival rate, and another 1/S for the emitted-to-observed wavelength-bin Jacobian. Inserting this into the band integral gives predicted counts through the actual instrument response. This is the appropriate next bridge between the model and recorded multi-band flux, rather than replacing a filter by its central wavelength. It still requires an emitted spectrum and an independently justified source-population model.

## Remaining requirements

The photometry is pinned to a later release than this 2024 archive. Compatibility with that release, exact zero-point application and calibration covariance remain to audit. This is an original-calibration resource set, not an automatic adoption of a newer recalibration. The recovered Hsiao template is a dependency, not an independently verified clock or luminosity prior; its source-model assumptions must be audited before use. No correction has been applied to the observed light curves.

Known source-duration/luminosity evolution ambiguity remains, as derived in the joint-light identifiability report. Filter acquisition alone does not resolve it. The synthetic timing calibration runs separately; none of the six scientific demonstrations is completed by this acquisition.

## Reproduce

```powershell
python research_work/results/des-photometric-calibration/acquire.py
python research_work/results/des-photometric-calibration/extract_dependencies.py
python research_work/results/des-photometric-calibration/verify.py
```
