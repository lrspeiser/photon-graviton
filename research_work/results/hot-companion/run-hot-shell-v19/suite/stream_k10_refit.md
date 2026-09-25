# Regression suite: stream_k10_refit (quick tier)

Law: a = 6.1158e-11 m/s^2; g_d = 2.0273e-10 m/s^2 (lambda 3.315); u = 119.9 km/s; reach 1.59 Mpc; gradual release over 3e+04 AU; refit: a, u; distance scale alpha x0.9500 (H0-like 70.9); SPARC distances static; heard through a stream absorbing inward waves at 10/Mpc

Round 19: a stream absorbing inward waves at 10 per Mpc (100 kpc), a and u refitted.

Run 2026-09-25 03:31 UTC on commit 36e6ef5+changes, 29 s. Baseline: baseline.json (round12, full, 2026-09-24 04:37 UTC).

| | count |
|---|---:|
| pass | 32 |
| close | 10 |
| fail | 7 |
| info | 11 |
| error | 0 |
| vs baseline: regressed | 8 |
| vs baseline: improved | 3 |
| vs baseline: worse | 5 |
| vs baseline: better | 10 |
| vs baseline: changed | 23 |
| vs baseline: new | 0 |
| vs baseline: known | 5 |
| vs baseline: error | 0 |


## machinery

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| deep-regime strength a | 6.116e-11 m/s^2 | constant of the law under test | info | changed (was 6.298e-11, info) |
| release scale g_d | 2.027e-10 m/s^2 | constant of the law under test | info | (was 2.027e-10, info) |
| companion speed u | 119.9 km/s | constant of the law under test | info | changed (was 169.4, info) |
| reach u x 13 Gyr | 1594.0 kpc | constant of the law under test | info | changed (was 2252.8, info) |
| companion emission per kg of matter, l = a u / 2 | 3.666e-06 W/kg | energy budget (round 3: 6.48e-6) | info | changed (was 5.336e-06, info) |
| standing rule: not MOND, not Newton, not dark matter | 1.148 dex | no forbidden class in the guard verdict; value = spread of the pull with random speed 0-1000 km/s | pass | changed (was 0.9994, pass) |
| SPARC pull is not a function of g_N alone | 0.04377 dex | >= 0.005 dex beyond a local function of g_N | pass | changed (was 0.02749, pass) |
| extra pull at 100 g_d, relative to Newton | 0  | < 1e-12 | pass | (was 0, pass) |

## galaxies

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| SPARC typical speed miss (all galaxies) | 16.23 km/s | no worse than MOND (16.13); within 5% is close | close | REGRESSED (was 15.87, pass) |
| SPARC mean squared log residual (the fitted statistic) | 0.2129  |  | info | changed (was 0.2111, info) |
| SPARC typical speed miss (test galaxies) | 12.91 km/s | no worse than MOND (13.52); within 5% is close | pass | changed (was 12.35, pass) |
| SPARC typical speed miss (validation galaxies) | 18.92 km/s | no worse than MOND (19.80); within 5% is close | pass | changed (was 19.01, pass) |
| SPARC median residual log10(observed/predicted pull) | 0.02212 dex | 0 +- 0.02 dex (stellar mass-to-light zero point) | pass | better (was 0.02713, pass) |
| SPARC median residual where the release switches (g_N 10^-10.25 to 10^-9.5) | 0.01219 dex | 0 +- 0.02 dex | pass | better (was 0.02706, pass) |
| bulge-dominated galaxies (25): typical speed miss | 31.71 km/s | no worse than MOND (30.35) | close | REGRESSED (was 29.44, pass) |

## clusters

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| X-COP: typical mass miss, rms ln(M_hydrostatic/M_predicted) (static distances) | 0.2111  | <= 0.30 (hydrostatic masses carry 10-20% bias and error); MOND 1.06, NFW fit 0.10 | pass | changed (was 0.2215, pass) |
| X-COP: worst mean miss at any of the six radii (0.1-1 R500) | 0.2711  | <= 0.25 (the usual hydrostatic bias at R500) | close | REGRESSED (was 0.2389, pass) |
| the companion speed X-COP prefers for these a and g_d | 119.9 km/s | law uses 119.9 | info | changed (was 169.4, info) |

## lensing

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| KiDS lensing pull, all isolated lenses (mixture): median log10(observed/predicted), static distances | -0.01208 dex | 0 +- 0.025 dex (7 reliable bins) | pass | improved (was 0.06476, close) |
| KiDS lensing pull, blue lenses (k = 1.01, measured): median log10(observed/predicted), static distances | -0.03272 dex | 0 +- 0.025 dex (7 reliable bins) | pass | worse (was 0.01116, pass) |
| KiDS lensing pull, red lenses (k = 3.94, measured): median log10(observed/predicted), static distances | -0.0144 dex | 0 +- 0.025 dex (7 reliable bins) | pass | improved (was 0.07776, fail) |
| KiDS lensing pull, Sersic n < 2 (disks; k = 0.56): median log10(observed/predicted), static distances | 0.03674 dex | 0 +- 0.025 dex (7 reliable bins) | pass | improved (was 0.06627, close) |
| KiDS lensing pull, Sersic n > 2 (bulges; k = 3.95): median log10(observed/predicted), static distances | -0.06232 dex | 0 +- 0.025 dex (7 reliable bins) | close | REGRESSED (was 0.03364, pass) |
| KiDS lensing pull, GAMA spectroscopic lenses: median log10(observed/predicted), static distances | -0.06752 dex | 0 +- 0.025 dex (7 reliable bins) | close | REGRESSED (was 0.009316, pass) |
| KiDS: red lenses lens more than blue ones of the same visible mass (colour gap) | 0.1732 dex | 0.153 +- 0.04 dex (same bins) | pass | better (was 0.1263, pass) |
| KiDS: bulge-dominated lenses lens more than disk-dominated ones (Sersic gap) | 0.22 dex | 0.154 +- 0.04 dex (same bins) | pass | worse (was 0.1553, pass) |
| lensing circular speeds 50-300 kpc, spirals: rms z (static distances) | 2.729  | Mistele et al. 2024 (ApJL 969, L3), Table 1 | close | changed (was 2.859, close) |
| lensing circular speeds 50-300 kpc, ellipticals: rms z (static distances) | 7.668  | Mistele et al. 2024 (ApJL 969, L3), Table 1 | fail | REGRESSED (was 2.553, close) |
| SLACS: stellar mass from lensing minus from star speeds (mean of 6) | -0.09383 dex | 0 +- 0.022 (the lenses' own scatter) | fail | REGRESSED (was -0.02775, pass) |
| SLACS: stellar mass needed for the Einstein radii, vs Chabrier (static distances) | 0.3556 dex | Salpeter is +0.25; the IMF trend of Posacki et al. 2015 is derived with dark haloes | info | changed (was 0.4264, info) |
| Einstein Cross: companion share of the pull at the ring | 1.058e-12  | dark matter inside the ring < 15% (Trott et al. 2010; van de Ven et al. 2010) | pass | changed (was 1.073e-12, pass) |
| bulge microlensing: companion share at a star's Einstein radius | 0  | standard lensing by stars (OGLE, MOA optical depths) | pass | (was 0, pass) |

## milky_way

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| rotation speed at the Sun (8.2 kpc) | 210.8 km/s | 229-234 (Eilers 2019, Zhou 2023, Ou 2024, Jiao 2023); +-7 for the visible-matter model | close | better (was 209.2, close) |
| rotation curve 15-27 kpc: mean offset over the four Gaia analyses | -6.529 km/s | 0 +- 5 (the analyses differ by about 10) | pass | better (was -8.672, pass) |
| rotation curve 5-27 kpc: typical miss against Eilers et al. 2019 | 14.6 km/s |  | info | better (was 15.94, info) |
| vertical pull 1.1 kpc above the Sun (as surface density), counted local matter | 71.62 Msun/pc^2 | 69.8 +- 3.3 (Bovy & Rix 2013: 68 +- 4; Holmberg & Flynn 2004: 74 +- 6) | pass | better (was 72.29, pass) |
| mass inside 20 kpc (as inferred from orbits) | 1.704e+11 Msun | 1.91e+11 +- 1.8e+10 (Posti & Helmi 2019) | pass | better (was 1.661e+11, pass) |
| mass inside 50 kpc (as inferred from orbits) | 3.671e+11 Msun | 4.50e+11 +- 4.0e+10 (Correa Magnus & Vasiliev 2022) | close | better (was 3.555e+11, close) |
| mass inside 100 kpc (as inferred from orbits) | 6.765e+11 Msun | 6.90e+11 +- 7.0e+10 (Deason 2021 + Correa Magnus & Vasiliev 2022) | pass | better (was 6.533e+11, pass) |
| mass inside 200 kpc (as inferred from orbits) | 1.289e+12 Msun | 1.10e+12 +- 2.5e+11 (Correa Magnus & Vasiliev 2022) | pass | worse (was 1.243e+12, pass) |
| escape speed at the Sun (to 400 kpc) | 511.4 km/s | 445-580 (seven analyses 2014-2024) | pass | changed (was 504.4, pass) |
| inner Galaxy (bar region, 2-3 kpc): share of the rotation speed from visible matter | 0.9801  | 0.88 +- 0.07 (Wegg, Gerhard & Portail 2016, bulge microlensing) | pass | changed (was 0.9842, pass) |

## dwarfs

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| Draco: stars' speed spread | 2.718 km/s | 9.1 +- 1.2 | fail | known (was 2.788, fail) |
| Ursa Minor: stars' speed spread | 3.365 km/s | 9.5 +- 1.2 | fail | known (was 3.442, fail) |
| Sculptor: stars' speed spread | 6.544 km/s | 9.2 +- 1.4 | pass | changed (was 6.639, pass) |
| Sextans: stars' speed spread | 2.048 km/s | 7.9 +- 1.3 | fail | known (was 2.104, fail) |
| Carina: stars' speed spread | 3.38 km/s | 6.6 +- 1.2 | close | changed (was 3.472, close) |
| Fornax: stars' speed spread | 12.2 km/s | 11.7 +- 0.9 | pass | changed (was 12.25, pass) |
| Leo II: stars' speed spread | 5.191 km/s | 6.6 +- 0.7 | close | REGRESSED (was 5.234, pass) |
| Leo I: stars' speed spread | 9.228 km/s | 9.2 +- 1.4 | pass | changed (was 9.264, pass) |
| Crater II: stars' speed spread | 0.9692 km/s | 2.7 +- 0.3 | fail | known (was 0.9986, fail) |
| Antlia 2: stars' speed spread | 1.073 km/s | 5.71 +- 1.08 | fail | known (was 1.105, fail) |
| all ten dwarfs: chi-squared (MOND with the same stars: 119) | 141.6  |  | info | changed (was 137.7, info) |

## precision

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| planets Mercury-Neptune: largest extra pull relative to Newton | 0  | < 1e-12 (planetary ephemerides) | pass | (was 0, pass) |
| star S2 at the Galactic Centre (far point): extra pull relative to Newton | 0  | GR precession 1.10 +- 0.19 (GRAVITY 2020) | pass | (was 0, pass) |
| Double Pulsar: orbit change from the companion emission, in measurement errors | 0.00924  | < 0.5 of the error on the orbital decay (Kramer et al. 2021) | pass | changed (was 0.01345, pass) |
| light bending at the Sun's limb: extra relative to GR | 0  | gamma - 1 = (2.1 +- 2.3) x 10^-5 (Cassini) | pass | (was 0, pass) |
| Cassini: Galactic distortion of the Sun's field, Q2 | 4.316e-27 1/s^2 | (1.6 +- 1.8) x 10^-27 (Park et al. 2026) | pass | worse (was 4.377e-27, pass) |
| Cassini Q2 with the Milky Way model's pull and heat at the Sun | 2.970e-27 1/s^2 |  | info | worse (was 3.472e-27, info) |
| wide binaries (1 Msun, 20,000 AU): boost of the pull | 1.074  | between the two published analyses, 1.0 (Banik et al. 2024) and 1.4-1.5 (Chae 2023-24) | pass | changed (was 1.076, pass) |
