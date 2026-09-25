# Regression suite: dist_metric_eta_r12 (full tier)

Law: a = 6.6320e-11 m/s^2; g_d = 2.0273e-10 m/s^2 (lambda 3.057); u = 169.0 km/s; reach 2.25 Mpc; gradual release over 3e+04 AU; refit: a, u; distance scale alpha x0.9814 (H0-like 73.2); SPARC distances static; distance geometry metric (D_A = D/(1 + z)); path factor sqrt(1 + 0.5 z/(1 + z))

Round 21 (registered comparison): the supplied distance law, D* = (ln(1 + z)/alpha) sqrt(1 + z/(2(1 + z))), D_A = D*/(1 + z), D_L = (1 + z) D*, alpha from the supernovae at eta = 1/2 (H0-like 73.2); a and u refitted.

Run 2026-09-25 18:59 UTC on commit 1dfcc49+changes, 988 s. Baseline: baseline.json (round12, full, 2026-09-24 04:37 UTC).

| | count |
|---|---:|
| pass | 62 |
| close | 8 |
| fail | 7 |
| info | 13 |
| error | 0 |
| vs baseline: regressed | 1 |
| vs baseline: improved | 3 |
| vs baseline: worse | 16 |
| vs baseline: better | 10 |
| vs baseline: changed | 48 |
| vs baseline: new | 0 |
| vs baseline: known | 5 |
| vs baseline: error | 0 |


## machinery

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| deep-regime strength a | 6.632e-11 m/s^2 | constant of the law under test | info | changed (was 6.298e-11, info) |
| release scale g_d | 2.027e-10 m/s^2 | constant of the law under test | info | (was 2.027e-10, info) |
| companion speed u | 169.0 km/s | constant of the law under test | info | changed (was 169.4, info) |
| reach u x 13 Gyr | 2246.4 kpc | constant of the law under test | info | changed (was 2252.8, info) |
| companion emission per kg of matter, l = a u / 2 | 5.603e-06 W/kg | energy budget (round 3: 6.48e-6) | info | changed (was 5.336e-06, info) |
| standing rule: not MOND, not Newton, not dark matter | 1.001 dex | no forbidden class in the guard verdict; value = spread of the pull with random speed 0-1000 km/s | pass | changed (was 0.9994, pass) |
| SPARC pull is not a function of g_N alone | 0.02737 dex | >= 0.005 dex beyond a local function of g_N | pass | changed (was 0.02749, pass) |
| extra pull at 100 g_d, relative to Newton | 0  | < 1e-12 | pass | (was 0, pass) |

## galaxies

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| SPARC typical speed miss (all galaxies) | 15.94 km/s | no worse than MOND (16.13); within 5% is close | pass | changed (was 15.87, pass) |
| SPARC mean squared log residual (the fitted statistic) | 0.2129  |  | info | changed (was 0.2111, info) |
| SPARC typical speed miss (test galaxies) | 12.44 km/s | no worse than MOND (13.52); within 5% is close | pass | changed (was 12.35, pass) |
| SPARC typical speed miss (validation galaxies) | 19.12 km/s | no worse than MOND (19.80); within 5% is close | pass | changed (was 19.01, pass) |
| SPARC median residual log10(observed/predicted pull) | 0.03295 dex | 0 +- 0.02 dex (stellar mass-to-light zero point) | pass | worse (was 0.02713, pass) |
| SPARC median residual where the release switches (g_N 10^-10.25 to 10^-9.5) | 0.03682 dex | 0 +- 0.02 dex | pass | worse (was 0.02706, pass) |
| bulge-dominated galaxies (25): typical speed miss | 29.64 km/s | no worse than MOND (30.35) | pass | changed (was 29.44, pass) |

## clusters

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| X-COP: typical mass miss, rms ln(M_hydrostatic/M_predicted) (static distances) | 0.222  | <= 0.30 (hydrostatic masses carry 10-20% bias and error); MOND 1.06, NFW fit 0.10 | pass | changed (was 0.2215, pass) |
| X-COP: worst mean miss at any of the six radii (0.1-1 R500) | 0.2359  | <= 0.25 (the usual hydrostatic bias at R500) | pass | changed (was 0.2389, pass) |
| the companion speed X-COP prefers for these a and g_d | 169.0 km/s | law uses 169.0 | info | changed (was 169.4, info) |
| X-COP: held-out typical miss, u fitted on the other half | 0.2374  | <= 0.30 | pass | changed (was 0.2365, pass) |

## lensing

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| KiDS lensing pull, all isolated lenses (mixture): median log10(observed/predicted), static distances | 0.02684 dex | 0 +- 0.025 dex (7 reliable bins) | pass | improved (was 0.06476, close) |
| KiDS lensing pull, blue lenses (k = 0.51, measured): median log10(observed/predicted), static distances | -0.03154 dex | 0 +- 0.025 dex (7 reliable bins) | pass | worse (was 0.01116, pass) |
| KiDS lensing pull, red lenses (k = 1.98, measured): median log10(observed/predicted), static distances | 0.04238 dex | 0 +- 0.025 dex (7 reliable bins) | pass | improved (was 0.07776, fail) |
| KiDS lensing pull, Sersic n < 2 (disks; k = 0.28): median log10(observed/predicted), static distances | 0.02861 dex | 0 +- 0.025 dex (7 reliable bins) | pass | improved (was 0.06627, close) |
| KiDS lensing pull, Sersic n > 2 (bulges; k = 1.99): median log10(observed/predicted), static distances | -0.001737 dex | 0 +- 0.025 dex (7 reliable bins) | pass | better (was 0.03364, pass) |
| KiDS lensing pull, GAMA spectroscopic lenses: median log10(observed/predicted), static distances | -0.0286 dex | 0 +- 0.025 dex (7 reliable bins) | pass | worse (was 0.009316, pass) |
| KiDS: red lenses lens more than blue ones of the same visible mass (colour gap) | 0.1247 dex | 0.153 +- 0.04 dex (same bins) | pass | changed (was 0.1263, pass) |
| KiDS: bulge-dominated lenses lens more than disk-dominated ones (Sersic gap) | 0.1533 dex | 0.154 +- 0.04 dex (same bins) | pass | changed (was 0.1553, pass) |
| lensing circular speeds 50-300 kpc, spirals: rms z (static distances) | 2.49  | Mistele et al. 2024 (ApJL 969, L3), Table 1 | close | better (was 2.859, close) |
| lensing circular speeds 50-300 kpc, ellipticals: rms z (static distances) | 4.027  | Mistele et al. 2024 (ApJL 969, L3), Table 1 | fail | REGRESSED (was 2.553, close) |
| SLACS: stellar mass from lensing minus from star speeds (mean of 6) | -0.02634 dex | 0 +- 0.023 (the lenses' own scatter) | pass | changed (was -0.02775, pass) |
| SLACS: stellar mass needed for the Einstein radii, vs Chabrier (static distances) | 0.4424 dex | Salpeter is +0.25; the IMF trend of Posacki et al. 2015 is derived with dark haloes | info | changed (was 0.4264, info) |
| Einstein Cross: companion share of the pull at the ring | 1.102e-12  | dark matter inside the ring < 15% (Trott et al. 2010; van de Ven et al. 2010) | pass | changed (was 1.073e-12, pass) |
| bulge microlensing: companion share at a star's Einstein radius | 0  | standard lensing by stars (OGLE, MOA optical depths) | pass | (was 0, pass) |

## milky_way

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| rotation speed at the Sun (8.2 kpc) | 210.0 km/s | 229-234 (Eilers 2019, Zhou 2023, Ou 2024, Jiao 2023); +-7 for the visible-matter model | close | changed (was 209.2, close) |
| rotation curve 15-27 kpc: mean offset over the four Gaia analyses | -7.137 km/s | 0 +- 5 (the analyses differ by about 10) | pass | better (was -8.672, pass) |
| rotation curve 5-27 kpc: typical miss against Eilers et al. 2019 | 15.16 km/s |  | info | changed (was 15.94, info) |
| vertical pull 1.1 kpc above the Sun (as surface density), counted local matter | 72.89 Msun/pc^2 | 69.8 +- 3.3 (Bovy & Rix 2013: 68 +- 4; Holmberg & Flynn 2004: 74 +- 6) | pass | worse (was 72.29, pass) |
| mass inside 20 kpc (as inferred from orbits) | 1.688e+11 Msun | 1.91e+11 +- 1.8e+10 (Posti & Helmi 2019) | pass | better (was 1.661e+11, pass) |
| mass inside 50 kpc (as inferred from orbits) | 3.632e+11 Msun | 4.50e+11 +- 4.0e+10 (Correa Magnus & Vasiliev 2022) | close | better (was 3.555e+11, close) |
| mass inside 100 kpc (as inferred from orbits) | 6.689e+11 Msun | 6.90e+11 +- 7.0e+10 (Deason 2021 + Correa Magnus & Vasiliev 2022) | pass | better (was 6.533e+11, pass) |
| mass inside 200 kpc (as inferred from orbits) | 1.274e+12 Msun | 1.10e+12 +- 2.5e+11 (Correa Magnus & Vasiliev 2022) | pass | worse (was 1.243e+12, pass) |
| escape speed at the Sun (to 400 kpc) | 509.2 km/s | 445-580 (seven analyses 2014-2024) | pass | changed (was 504.4, pass) |
| inner Galaxy (bar region, 2-3 kpc): share of the rotation speed from visible matter | 0.9837  | 0.88 +- 0.07 (Wegg, Gerhard & Portail 2016, bulge microlensing) | pass | (was 0.9842, pass) |

## dwarfs

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| Draco: stars' speed spread | 2.817 km/s | 9.1 +- 1.2 | fail | known (was 2.788, fail) |
| Ursa Minor: stars' speed spread | 3.479 km/s | 9.5 +- 1.2 | fail | known (was 3.442, fail) |
| Sculptor: stars' speed spread | 6.709 km/s | 9.2 +- 1.4 | pass | changed (was 6.639, pass) |
| Sextans: stars' speed spread | 2.127 km/s | 7.9 +- 1.3 | fail | known (was 2.104, fail) |
| Carina: stars' speed spread | 3.511 km/s | 6.6 +- 1.2 | close | changed (was 3.472, close) |
| Fornax: stars' speed spread | 12.38 km/s | 11.7 +- 0.9 | pass | worse (was 12.25, pass) |
| Leo II: stars' speed spread | 5.293 km/s | 6.6 +- 0.7 | pass | changed (was 5.234, pass) |
| Leo I: stars' speed spread | 9.356 km/s | 9.2 +- 1.4 | pass | changed (was 9.264, pass) |
| Crater II: stars' speed spread | 1.01 km/s | 2.7 +- 0.3 | fail | known (was 0.9986, fail) |
| Antlia 2: stars' speed spread | 1.118 km/s | 5.71 +- 1.08 | fail | known (was 1.105, fail) |
| all ten dwarfs: chi-squared (MOND with the same stars: 119) | 136.0  |  | info | changed (was 137.7, info) |

## precision

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| planets Mercury-Neptune: largest extra pull relative to Newton | 0  | < 1e-12 (planetary ephemerides) | pass | (was 0, pass) |
| star S2 at the Galactic Centre (far point): extra pull relative to Newton | 0  | GR precession 1.10 +- 0.19 (GRAVITY 2020) | pass | (was 0, pass) |
| Double Pulsar: orbit change from the companion emission, in measurement errors | 0.01412  | < 0.5 of the error on the orbital decay (Kramer et al. 2021) | pass | changed (was 0.01345, pass) |
| light bending at the Sun's limb: extra relative to GR | 0  | gamma - 1 = (2.1 +- 2.3) x 10^-5 (Cassini) | pass | (was 0, pass) |
| Cassini: Galactic distortion of the Sun's field, Q2 | 4.485e-27 1/s^2 | (1.6 +- 1.8) x 10^-27 (Park et al. 2026) | pass | worse (was 4.377e-27, pass) |
| Cassini Q2 with the Milky Way model's pull and heat at the Sun | 3.560e-27 1/s^2 |  | info | worse (was 3.472e-27, info) |
| wide binaries (1 Msun, 20,000 AU): boost of the pull | 1.079  | between the two published analyses, 1.0 (Banik et al. 2024) and 1.4-1.5 (Chae 2023-24) | pass | changed (was 1.076, pass) |

## collisions

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| Bullet main cluster: outer stars needed for the galaxy speed (static distances) | 3.700e+12 Msun | 3.5-6.0 x 10^12 (Legacy Survey star count, round 6, converted) | pass | changed (was 3.711e+12, pass) |
| Bullet main cluster: lensing strength (kappa) on its galaxies | 0.6968  | at least 0.36 +- 0.06 (a stated lower bound) | pass | changed (was 0.7168, pass) |
| Bullet subcluster: lensing strength (kappa) on its galaxies | 0.1466  | at least 0.20 +- 0.05 (a stated lower bound) | pass | worse (was 0.2764, pass) |
| Bullet: leftover lensing on the main cluster's gas | 0.03251  | 0.05 +- 0.06 | pass | worse (was 0.04957, pass) |
| Bullet: leftover lensing on the subcluster's gas | 0.04205  | 0.02 +- 0.06 | pass | better (was 0.0516, pass) |
| Bullet main: lensing peak distance from its galaxies | 8.339 kpc | within a quarter of the galaxy-gas separation (200 kpc) | pass | changed (was 15.68, pass) |
| Bullet sub: lensing peak distance from its galaxies | 13.82 kpc | within a quarter of the galaxy-gas separation (184 kpc) | pass | changed (was 22.97, pass) |
| Bullet main cluster: lensing mass inside 237 kpc (250 in LCDM units) | 2.361e+14 Msun | 2.36-2.64 x 10^14 (Paraficz 2016; Bradac 2006; converted) | pass | better (was 2.969e+14, pass) |
| Bullet subcluster: lensing mass inside 237 kpc (250 in LCDM units) | 9.103e+13 Msun | 1.89-2.17 x 10^14 (Paraficz 2016; Bradac 2006; converted) | fail | worse (was 1.385e+14, fail) |
| collision stack: how far the lensing follows the gas (beta, median of 20) | 0.01782  | -0.04 +- 0.07 (Harvey et al. 2015, 72 collisions) | pass | changed (was 0.01914, pass) |
| collision stack: lensing peak distance from the galaxies (median) | 15.61 kpc | raw offset; includes the pull of the other cluster | info | changed (was 15.75, info) |
| MACS J0025 SE subcluster: lensing mass inside 300 kpc (LCDM units; 285 kpc static) | 1.390e+14 Msun | 2.31 (+0.92/-1.57) x 10^14 | pass | changed (was 2.053e+14, pass) |
| MACS J0025 NW subcluster: lensing mass inside 300 kpc (LCDM units; 285 kpc static) | 1.181e+14 Msun | 2.40 (+0.46/-1.29) x 10^14 | pass | changed (was 1.853e+14, pass) |
| MACS J0025 SE: lensing peak distance from its galaxies | 15.76 kpc | with the galaxies, not the gas (352 kpc away): lensing offset from the gas at > 4 sigma | pass | changed (was 23.61, pass) |
| MACS J0025 NW: lensing peak distance from its galaxies | 49.12 kpc | with the galaxies, not the gas (162 kpc away): lensing offset from the gas at > 4 sigma | pass | changed (was 79.2, pass) |
| MACS J0025: galaxies' speed spread (line of sight, inside 1.5 Mpc) | 669.9 km/s | 835 +- 59 (108 galaxies within 1.5 Mpc) | close | changed (was 668.7, close) |
| Abell 520 P1, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 2.667e+13 Msun | 1.99-2.67 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | changed (was 3.158e+13, pass) |
| Abell 520 P2, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 5.177e+13 Msun | 3.84-3.95 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | worse (was 6.056e+13, pass) |
| Abell 520 P3, the gas-rich, galaxy-poor centre ("dark core"): lensing mass inside 150 kpc (LCDM units) | 3.077e+13 Msun | 2.70-3.18 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | changed (was 3.717e+13, pass) |
| Abell 520 P4, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 3.527e+13 Msun | 4.01-5.31 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | better (was 4.117e+13, pass) |
| Abell 520 P5, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 2.883e+13 Msun | 2.78-3.01 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | changed (was 3.453e+13, pass) |
| Abell 520 P6, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 2.070e+13 Msun | 3.49-3.49 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | close | changed (was 2.527e+13, close) |
| Abell 520: lensing mass inside 710 kpc of the centre (LCDM units) | 4.368e+14 Msun | 4.76 +- 0.52 x 10^14 in the static law (Mahdavi et al. 2007; <D_ls/D_s> = 0.59) | pass | better (was 5.272e+14, pass) |
| Abell 520: galaxies' speed spread at P1, P2, P4, P5 (rms z) | 2.659  | Girardi et al. 2008 Table 4: 811, 749, 579, 668 km/s (6-9 galaxies each) | close | changed (was 2.697, close) |
| El Gordo: lensing mass inside 500 kpc of the centre of mass (Planck units) | 4.559e+14 Msun | 5.13e+14 +- 12% (aperture mass) | pass | changed (was 8.423e+14, pass) |
| El Gordo: lensing mass inside 1000 kpc of the centre of mass (Planck units) | 1.136e+15 Msun | 1.32e+15 +- 12% (aperture mass) | pass | worse (was 2.142e+15, pass) |
| El Gordo NW: galaxies' speed spread (line of sight, inside 1 Mpc) | 994.9 km/s | 1290 +- 134 (Menanteau et al. 2012) | close | worse (was 1013.7, close) |
| El Gordo SE: galaxies' speed spread (line of sight, inside 1 Mpc) | 881.6 km/s | 1089 +- 200 (Menanteau et al. 2012) | pass | worse (was 903.6, pass) |
| El Gordo SE: lensing peak distance from the cool gas core | 144.8 kpc | about 12 arcsec, the lensing peak nearer the merger centre (Kim et al. 2021; Ng et al. 2015) | info | changed (was 248.4, info) |
