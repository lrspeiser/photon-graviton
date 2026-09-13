# Radiation-conditioned capture gives mixed galaxy transfer

The proposed susceptibility multiplier is F=[(L_3.6/10^9 L_sun)/(R_disk/kpc)^2]^q. It multiplies both k0 and C, preserving their ratio and the assumed common external exposure. This is a postulated capture response using a band-luminosity proxy, not a derived microscopic interaction or an independent companion source.

Fitting only the 89 existing training galaxies gives q=0.4644964204, C=8.267782946e7 Msun/kpc^3, k0=0.4297539287/kpc and a/R_disk=2.252509382. Freeze all four for the 29 validation and 31 test galaxies.

| Sample | Retained RMS (km/s) | New RMS (km/s) | Retained log RMS | New log RMS |
|---|---:|---:|---:|---:|
| Training | 30.70873 | 28.82801 | 0.141004 | 0.135166 |
| Validation | 33.18419 | 35.34063 | 0.119980 | 0.115922 |
| Test | 23.99377 | 24.50764 | 0.098545 | 0.094572 |

The actual training objective and its fractional-error transfer measure improve. Absolute speed RMS worsens in both comparison groups. Both facts matter; this is a mixed candidate, not an across-the-board improvement or a rejected fit solely because a secondary metric worsens. Keep it for further physical/observable comparison without replacing the retained reference.

Three starts converge, all succeed, and no parameter reaches a bound. q=0 reproduces the old model numerically before fitting. Finer integration changes group RMS by at most 0.00110 km/s. Published stellar mass-to-light ratios, distances and ordinary-matter inputs remain fixed assumptions, and every partition was previously exposed. There is no independent significance or cause-identification claim.

Higher radiative surface-intensity proxy increases capture susceptibility in this fit. That is a different response and target from the inverse-flux effective-retention Solar System full-mass family. A single universal equation has not been found merely because both mention radiation. The galaxy model addresses additional gravity above ordinary matter, while the Solar System exercise attempted all source mass. Moreover, susceptibility and lifetime are different quantities.

Next evaluation of this candidate must use the same deposited source for motion and lensing and specify the luminosity/morphology mapping; no lens normalization or per-object photon supply may be added to rescue a mismatch. Energy origin, support, independent confinement and microscopic graviton production remain open.

Reproduce with `python research_work/results/isotropic-galaxy-transfer/run.py --radiative-flux`. Complete parameters and all radial observed/predicted speeds are in `radiative-flux-results.json` and `radiative-flux-predictions.json`. All six broader goals remain open.
