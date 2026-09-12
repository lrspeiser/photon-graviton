# Combined redshift, timing and detector-flux forward calculator

The verified propagation bookkeeping, real DES passbands and historical calibration convention are now combined in `transport.py`. Given an independently supplied distance and an emitted spectral luminosity as a function of source time, it predicts redshift, event stretching, photons per area per time and native DES FLUXCAL. No observed redshift or per-event fitted brightness is required by this calculator.

This is a conditional effective-model implementation, not a derived photon-companion interaction or a demonstrated fit to real observations.

## Inputs and formula provenance

The source callback returns L_lambda(t_e, lambda_e) in erg/s/angstrom, with time in days and wavelength in angstrom. Source luminosity, spectrum and evolution must be specified separately rather than adjusted to each desired result. Distance is supplied in Mpc. Observer time is elapsed time relative to the arrival corresponding to source time zero; no absolute travel time is inferred.

**Existing fitted description, using known fractional-loss mathematics:**

\[
S=\exp(\alpha D),\qquad z_{\rm predicted}=S-1.
\]

The default alpha is the earlier exposed-galaxy calibration, 0.0002488993286382367 Mpc^-1. It is not a new independent determination or a universal rate proven from first principles.

**Proposed shared-transport parameterization:** event stretch is S^b, with default b=1. b=0 remains the energy-only control. The implementation does not derive this event mapping from energy loss.

**Known radiometry conditional on those assumptions:**

\[
F_{\lambda_o}(t_o)=\frac{L_{\lambda_e}(t_o/S^b,\lambda_o/S)}{4\pi D^2 S^{b+2}},
\qquad
C_j(t_o)=\int T_j(\lambda_o)F_{\lambda_o}(t_o)\frac{\lambda_o}{hc}\,d\lambda_o.
\]

The code converts Mpc to cm and keeps the wavelength unit consistent with h*c. It integrates through each full verified filter response.

**Known magnitude calibration:** if C_AB0,j is the integrated zero-AB reference photon rate and delta_j is the checked historical offset,

\[
F_{\rm CAL,native,j}=10^{11}\frac{C_j}{C_{\rm AB0,j}}10^{-0.4\delta_j}.
\]

The factor 10^11 comes from the documented FLUXCAL magnitude zero point 27.5. It is not an adjustable physical gain. The calculator also returns photon rates and AB-equivalent FLUXCAL so the unit conversions remain inspectable.

## Verification

A flat-frequency reference source at 10 pc reproduces the four archived calibration reference magnitudes to relative flux error below 3.4e-16. Doubling the distance with zero conversion reduces flux by exactly four within the numerical tolerance. The artificial spectrum is restricted to 1000-20000 angstrom; every checked passband maps inside its support.

Twenty-four pulse cases span the four bands, b=0/1 and S=1/1.2/2.2. The measured half-maximum width scales as S^b for this wavelength-independent-duration source. Its fixed-band peak count scales as S^-b because the chosen flat-frequency spectrum contributes a compensating frequency-bin factor. That band-specific result is not the bolometric attenuation law, which remains S^-(1+b). Other source spectra must be integrated and can give different band brightness ratios.

At the illustrative distance of 100 million light-years (30.660139 Mpc), the default rate gives predicted z=0.0076604805. The saved brightness for that example uses the declared artificial reference source; it is not a measured star or galaxy and must not be cited as observational agreement.

## What is still missing

There is no dust, lensing, peculiar-velocity correction, unknown peak-date fitting, source-population likelihood or physical receiving-sector calculation in this module. The measured calibration covariance is available separately and must enter a joint likelihood; it is not silently treated as independent noise here. The source-duration/luminosity evolution ambiguity remains unless restricted or constrained by independent evidence.

This creates the forward calculation needed for a joint test. It does not establish that the default alpha, event mapping or source model fits real supernova redshift, timing and brightness together. The larger timing calibration continues under its frozen protocol, with boundary failures retained. All six scientific demonstrations remain open.

```powershell
python research_work/results/joint-light-forward/verify.py
```

The earlier DES calibration acquisition supplies the hash-verified filter files. The calculator accepts observation arrays containing g/r/i/z band labels and a vectorized emitted-luminosity callback; `verify.py` gives runnable examples with declared artificial sources.
