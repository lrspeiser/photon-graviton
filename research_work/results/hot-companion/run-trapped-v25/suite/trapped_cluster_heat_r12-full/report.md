# Regression suite: trapped_cluster_heat_r12 (full tier)

Law: a = 6.2979e-11 m/s^2; g_d = 2.0273e-10 m/s^2 (lambda 3.219); u = 169.4 km/s; reach 2.25 Mpc; gradual release over 3e+04 AU; refit: u; distance scale alpha x0.9500 (H0-like 70.9); SPARC distances static; trapped companion, cluster galaxies' own glow at cluster heat (code/trapping_v25.py)

round 25, registered (round25-trapped-suite.md): as trapped_own_heat_r12, with a cluster galaxy's own glow at the cluster's heat (the law's formula taken literally inside the galaxy); u refitted on X-COP

Run 2026-09-26 18:07 UTC on commit ed61d9c+changes, 1331 s. Baseline: baseline.json (round12, full, 2026-09-24 04:37 UTC).

| | count |
|---|---:|
| pass | 57 |
| close | 9 |
| fail | 11 |
| info | 13 |
| error | 0 |
| vs baseline: regressed | 6 |
| vs baseline: improved | 1 |
| vs baseline: worse | 11 |
| vs baseline: better | 2 |
| vs baseline: changed | 16 |
| vs baseline: new | 0 |
| vs baseline: known | 6 |
| vs baseline: error | 0 |


## machinery

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| deep-regime strength a | 6.298e-11 m/s^2 | constant of the law under test | info | (was 6.298e-11, info) |
| release scale g_d | 2.027e-10 m/s^2 | constant of the law under test | info | (was 2.027e-10, info) |
| companion speed u | 169.4 km/s | constant of the law under test | info | (was 169.4, info) |
| reach u x 13 Gyr | 2252.6 kpc | constant of the law under test | info | (was 2252.8, info) |
| companion emission per kg of matter, l = a u / 2 | 5.335e-06 W/kg | energy budget (round 3: 6.48e-6) | info | (was 5.336e-06, info) |
| standing rule: not MOND, not Newton, not dark matter | 0.9995 dex | no forbidden class in the guard verdict; value = spread of the pull with random speed 0-1000 km/s | pass | (was 0.9994, pass) |
| SPARC pull is not a function of g_N alone | 0.007459 dex | >= 0.005 dex beyond a local function of g_N | pass | changed (was 0.02749, pass) |
| extra pull at 100 g_d, relative to Newton | 0  | < 1e-12 | pass | (was 0, pass) |

## galaxies

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| SPARC typical speed miss (all galaxies) | 16.03 km/s | no worse than MOND (16.13); within 5% is close | pass | changed (was 15.87, pass) |
| SPARC mean squared log residual (the fitted statistic) | 0.2116  |  | info | changed (was 0.2111, info) |
| SPARC typical speed miss (test galaxies) | 13.07 km/s | no worse than MOND (13.52); within 5% is close | pass | changed (was 12.35, pass) |
| SPARC typical speed miss (validation galaxies) | 19.21 km/s | no worse than MOND (19.80); within 5% is close | pass | changed (was 19.01, pass) |
| SPARC median residual log10(observed/predicted pull) | 0.0343 dex | 0 +- 0.02 dex (stellar mass-to-light zero point) | pass | worse (was 0.02713, pass) |
| SPARC median residual where the release switches (g_N 10^-10.25 to 10^-9.5) | 0.0334 dex | 0 +- 0.02 dex | pass | worse (was 0.02706, pass) |
| bulge-dominated galaxies (25): typical speed miss | 29.95 km/s | no worse than MOND (30.35) | pass | changed (was 29.44, pass) |

## clusters

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| X-COP: typical mass miss, rms ln(M_hydrostatic/M_predicted) (static distances) | 0.2215  | <= 0.30 (hydrostatic masses carry 10-20% bias and error); MOND 1.06, NFW fit 0.10 | pass | (was 0.2215, pass) |
| X-COP: worst mean miss at any of the six radii (0.1-1 R500) | 0.2389  | <= 0.25 (the usual hydrostatic bias at R500) | pass | (was 0.2389, pass) |
| the companion speed X-COP prefers for these a and g_d | 169.4 km/s | law uses 169.4 | info | (was 169.4, info) |
| X-COP: held-out typical miss, u fitted on the other half | 0.2365  | <= 0.30 | pass | (was 0.2365, pass) |

## lensing

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| KiDS lensing pull, all isolated lenses (mixture): median log10(observed/predicted), static distances | 0.1059 dex | 0 +- 0.025 dex (7 reliable bins) | fail | REGRESSED (was 0.06476, close) |
| KiDS lensing pull, blue lenses (k = 0.41, measured): median log10(observed/predicted), static distances | 0.02251 dex | 0 +- 0.025 dex (7 reliable bins) | pass | worse (was 0.01116, pass) |
| KiDS lensing pull, red lenses (k = 1.18, measured): median log10(observed/predicted), static distances | 0.1384 dex | 0 +- 0.025 dex (7 reliable bins) | fail | worse (was 0.07776, fail) |
| KiDS lensing pull, Sersic n < 2 (disks; k = 0.16): median log10(observed/predicted), static distances | 0.08394 dex | 0 +- 0.025 dex (7 reliable bins) | fail | REGRESSED (was 0.06627, close) |
| KiDS lensing pull, Sersic n > 2 (bulges; k = 1.21): median log10(observed/predicted), static distances | 0.09221 dex | 0 +- 0.025 dex (7 reliable bins) | fail | REGRESSED (was 0.03364, pass) |
| KiDS lensing pull, GAMA spectroscopic lenses: median log10(observed/predicted), static distances | 0.05042 dex | 0 +- 0.025 dex (7 reliable bins) | close | REGRESSED (was 0.009316, pass) |
| KiDS: red lenses lens more than blue ones of the same visible mass (colour gap) | 0.07951 dex | 0.153 +- 0.04 dex (same bins) | pass | worse (was 0.1263, pass) |
| KiDS: bulge-dominated lenses lens more than disk-dominated ones (Sersic gap) | 0.1157 dex | 0.154 +- 0.04 dex (same bins) | pass | worse (was 0.1553, pass) |
| lensing circular speeds 50-300 kpc, spirals: rms z (static distances) | 3.038  | Mistele et al. 2024 (ApJL 969, L3), Table 1 | fail | REGRESSED (was 2.859, close) |
| lensing circular speeds 50-300 kpc, ellipticals: rms z (static distances) | 1.912  | Mistele et al. 2024 (ApJL 969, L3), Table 1 | pass | improved (was 2.553, close) |
| SLACS: stellar mass from lensing minus from star speeds (mean of 6) | -0.001845 dex | 0 +- 0.025 (the lenses' own scatter) | pass | better (was -0.02775, pass) |
| SLACS: stellar mass needed for the Einstein radii, vs Chabrier (static distances) | 0.4552 dex | Salpeter is +0.25; the IMF trend of Posacki et al. 2015 is derived with dark haloes | info | changed (was 0.4264, info) |
| Einstein Cross: companion share of the pull at the ring | 1.073e-12  | dark matter inside the ring < 15% (Trott et al. 2010; van de Ven et al. 2010) | pass | (was 1.073e-12, pass) |
| bulge microlensing: companion share at a star's Einstein radius | 0  | standard lensing by stars (OGLE, MOA optical depths) | pass | (was 0, pass) |

## milky_way

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| rotation speed at the Sun (8.2 kpc) | 209.1 km/s | 229-234 (Eilers 2019, Zhou 2023, Ou 2024, Jiao 2023); +-7 for the visible-matter model | close | (was 209.2, close) |
| rotation curve 15-27 kpc: mean offset over the four Gaia analyses | -11.2 km/s | 0 +- 5 (the analyses differ by about 10) | close | REGRESSED (was -8.672, pass) |
| rotation curve 5-27 kpc: typical miss against Eilers et al. 2019 | 16.47 km/s |  | info | changed (was 15.94, info) |
| vertical pull 1.1 kpc above the Sun (as surface density), counted local matter | 72.4 Msun/pc^2 | 69.8 +- 3.3 (Bovy & Rix 2013: 68 +- 4; Holmberg & Flynn 2004: 74 +- 6) | pass | changed (was 72.29, pass) |
| mass inside 20 kpc (as inferred from orbits) | 1.603e+11 Msun | 1.91e+11 +- 1.8e+10 (Posti & Helmi 2019) | pass | worse (was 1.661e+11, pass) |
| mass inside 50 kpc (as inferred from orbits) | 3.403e+11 Msun | 4.50e+11 +- 4.0e+10 (Correa Magnus & Vasiliev 2022) | close | worse (was 3.555e+11, close) |
| mass inside 100 kpc (as inferred from orbits) | 6.225e+11 Msun | 6.90e+11 +- 7.0e+10 (Deason 2021 + Correa Magnus & Vasiliev 2022) | pass | worse (was 6.533e+11, pass) |
| mass inside 200 kpc (as inferred from orbits) | 1.181e+12 Msun | 1.10e+12 +- 2.5e+11 (Correa Magnus & Vasiliev 2022) | pass | better (was 1.243e+12, pass) |
| escape speed at the Sun (to 400 kpc) | 496.0 km/s | 445-580 (seven analyses 2014-2024) | pass | changed (was 504.4, pass) |
| inner Galaxy (bar region, 2-3 kpc): share of the rotation speed from visible matter | 0.9842  | 0.88 +- 0.07 (Wegg, Gerhard & Portail 2016, bulge microlensing) | pass | (was 0.9842, pass) |

## dwarfs

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| Draco: stars' speed spread | 2.852 km/s | 9.1 +- 1.2 | fail | known (was 2.788, fail) |
| Ursa Minor: stars' speed spread | 3.516 km/s | 9.5 +- 1.2 | fail | known (was 3.442, fail) |
| Sculptor: stars' speed spread | 6.717 km/s | 9.2 +- 1.4 | pass | changed (was 6.639, pass) |
| Sextans: stars' speed spread | 2.154 km/s | 7.9 +- 1.3 | fail | known (was 2.104, fail) |
| Carina: stars' speed spread | 3.554 km/s | 6.6 +- 1.2 | close | changed (was 3.472, close) |
| Fornax: stars' speed spread | 12.28 km/s | 11.7 +- 0.9 | pass | changed (was 12.25, pass) |
| Leo II: stars' speed spread | 5.253 km/s | 6.6 +- 0.7 | pass | changed (was 5.234, pass) |
| Leo I: stars' speed spread | 9.266 km/s | 9.2 +- 1.4 | pass | (was 9.264, pass) |
| Crater II: stars' speed spread | 1.023 km/s | 2.7 +- 0.3 | fail | known (was 0.9986, fail) |
| Antlia 2: stars' speed spread | 1.133 km/s | 5.71 +- 1.08 | fail | known (was 1.105, fail) |
| all ten dwarfs: chi-squared (MOND with the same stars: 119) | 134.4  |  | info | changed (was 137.7, info) |

## precision

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| planets Mercury-Neptune: largest extra pull relative to Newton | 0  | < 1e-12 (planetary ephemerides) | pass | (was 0, pass) |
| star S2 at the Galactic Centre (far point): extra pull relative to Newton | 0  | GR precession 1.10 +- 0.19 (GRAVITY 2020) | pass | (was 0, pass) |
| Double Pulsar: orbit change from the companion emission, in measurement errors | 0.01345  | < 0.5 of the error on the orbital decay (Kramer et al. 2021) | pass | (was 0.01345, pass) |
| light bending at the Sun's limb: extra relative to GR | 0  | gamma - 1 = (2.1 +- 2.3) x 10^-5 (Cassini) | pass | (was 0, pass) |
| Cassini: Galactic distortion of the Sun's field, Q2 | 4.377e-27 1/s^2 | (1.6 +- 1.8) x 10^-27 (Park et al. 2026) | pass | worse (was 4.377e-27, pass) |
| Cassini Q2 with the Milky Way model's pull and heat at the Sun | 3.472e-27 1/s^2 |  | info | worse (was 3.472e-27, info) |
| wide binaries (1 Msun, 20,000 AU): boost of the pull | 1.076  | between the two published analyses, 1.0 (Banik et al. 2024) and 1.4-1.5 (Chae 2023-24) | pass | (was 1.076, pass) |

## collisions

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| Bullet main cluster: outer stars needed for the galaxy speed (static distances) | 3.709e+12 Msun | 3.4-5.8 x 10^12 (Legacy Survey star count, round 6, converted) | pass | (was 3.711e+12, pass) |
| Bullet main cluster: lensing strength (kappa) on its galaxies | 0.7169  | at least 0.36 +- 0.06 (a stated lower bound) | pass | (was 0.7168, pass) |
| Bullet subcluster: lensing strength (kappa) on its galaxies | 0.2765  | at least 0.20 +- 0.05 (a stated lower bound) | pass | (was 0.2764, pass) |
| Bullet: leftover lensing on the main cluster's gas | 0.04957  | 0.05 +- 0.06 | pass | (was 0.04957, pass) |
| Bullet: leftover lensing on the subcluster's gas | 0.0516  | 0.02 +- 0.06 | pass | (was 0.0516, pass) |
| Bullet main: lensing peak distance from its galaxies | 15.68 kpc | within a quarter of the galaxy-gas separation (252 kpc) | pass | (was 15.68, pass) |
| Bullet sub: lensing peak distance from its galaxies | 22.98 kpc | within a quarter of the galaxy-gas separation (232 kpc) | pass | (was 22.97, pass) |
| Bullet main cluster: lensing mass inside 301 kpc (250 in LCDM units) | 2.969e+14 Msun | 3.09-3.46 x 10^14 (Paraficz 2016; Bradac 2006; converted) | pass | (was 2.969e+14, pass) |
| Bullet subcluster: lensing mass inside 301 kpc (250 in LCDM units) | 1.385e+14 Msun | 2.47-2.85 x 10^14 (Paraficz 2016; Bradac 2006; converted) | fail | known (was 1.385e+14, fail) |
| collision stack: how far the lensing follows the gas (beta, median of 20) | 0.01913  | -0.04 +- 0.07 (Harvey et al. 2015, 72 collisions) | pass | (was 0.01914, pass) |
| collision stack: lensing peak distance from the galaxies (median) | 15.75 kpc | raw offset; includes the pull of the other cluster | info | (was 15.75, info) |
| MACS J0025 SE subcluster: lensing mass inside 300 kpc (LCDM units; 429 kpc static) | 2.053e+14 Msun | 3.64 (+1.46/-2.48) x 10^14 | pass | (was 2.053e+14, pass) |
| MACS J0025 NW subcluster: lensing mass inside 300 kpc (LCDM units; 429 kpc static) | 1.853e+14 Msun | 3.79 (+0.73/-2.04) x 10^14 | pass | (was 1.853e+14, pass) |
| MACS J0025 SE: lensing peak distance from its galaxies | 23.61 kpc | with the galaxies, not the gas (530 kpc away): lensing offset from the gas at > 4 sigma | pass | (was 23.61, pass) |
| MACS J0025 NW: lensing peak distance from its galaxies | 79.17 kpc | with the galaxies, not the gas (244 kpc away): lensing offset from the gas at > 4 sigma | pass | (was 79.2, pass) |
| MACS J0025: galaxies' speed spread (line of sight, inside 1.5 Mpc) | 668.6 km/s | 835 +- 59 (108 galaxies within 1.5 Mpc) | close | (was 668.7, close) |
| Abell 520 P1, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 3.158e+13 Msun | 2.44-3.26 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | (was 3.158e+13, pass) |
| Abell 520 P2, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 6.056e+13 Msun | 4.70-4.83 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | (was 6.056e+13, pass) |
| Abell 520 P3, the gas-rich, galaxy-poor centre ("dark core"): lensing mass inside 150 kpc (LCDM units) | 3.715e+13 Msun | 3.29-3.89 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | (was 3.717e+13, pass) |
| Abell 520 P4, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 4.117e+13 Msun | 4.91-6.48 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | (was 4.117e+13, pass) |
| Abell 520 P5, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 3.452e+13 Msun | 3.40-3.68 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | (was 3.453e+13, pass) |
| Abell 520 P6, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 2.524e+13 Msun | 4.27-4.27 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | close | changed (was 2.527e+13, close) |
| Abell 520: lensing mass inside 710 kpc of the centre (LCDM units) | 5.270e+14 Msun | 5.84 +- 0.64 x 10^14 in the static law (Mahdavi et al. 2007; <D_ls/D_s> = 0.59) | pass | (was 5.272e+14, pass) |
| Abell 520: galaxies' speed spread at P1, P2, P4, P5 (rms z) | 2.698  | Girardi et al. 2008 Table 4: 811, 749, 579, 668 km/s (6-9 galaxies each) | close | (was 2.697, close) |
| El Gordo: lensing mass inside 500 kpc of the centre of mass (Planck units) | 8.424e+14 Msun | 9.45e+14 +- 12% (aperture mass) | pass | (was 8.423e+14, pass) |
| El Gordo: lensing mass inside 1000 kpc of the centre of mass (Planck units) | 2.142e+15 Msun | 2.43e+15 +- 12% (aperture mass) | pass | (was 2.142e+15, pass) |
| El Gordo NW: galaxies' speed spread (line of sight, inside 1 Mpc) | 1013.7 km/s | 1290 +- 134 (Menanteau et al. 2012) | close | (was 1013.7, close) |
| El Gordo SE: galaxies' speed spread (line of sight, inside 1 Mpc) | 903.6 km/s | 1089 +- 200 (Menanteau et al. 2012) | pass | (was 903.6, pass) |
| El Gordo SE: lensing peak distance from the cool gas core | 248.4 kpc | about 12 arcsec, the lensing peak nearer the merger centre (Kim et al. 2021; Ng et al. 2015) | info | (was 248.4, info) |
