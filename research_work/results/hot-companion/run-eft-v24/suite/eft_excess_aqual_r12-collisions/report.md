# Regression suite: eft_excess_aqual_r12 (full tier)

Law: a = 6.2979e-11 m/s^2; g_d = 2.0273e-10 m/s^2 (lambda 3.219); u = 169.4 km/s; reach 2.25 Mpc; gradual release over 3e+04 AU; distance scale alpha x0.9500 (H0-like 70.9); SPARC distances static; collision maps: eft_excess_aqual (code/eft_field_v24.py)

Round 24: eft_excess_r12 (the repaired Casimir-EFT equation, with the law's memory) solved exactly on the collision maps, the curl field kept (code/eft_field_v24.py mode 'eft_excess_aqual': damped Helmholtz-projection iteration, relaxation 0.3, 60 steps). Round systems: the law exactly. Registered comparison, not adopted.

Run 2026-09-26 08:48 UTC on commit 2ee5c1c, 5801 s. Baseline: baseline.json (round12, full, 2026-09-24 04:37 UTC).

| | count |
|---|---:|
| pass | 21 |
| close | 4 |
| fail | 2 |
| info | 2 |
| error | 0 |
| vs baseline: regressed | 1 |
| vs baseline: improved | 0 |
| vs baseline: worse | 7 |
| vs baseline: better | 2 |
| vs baseline: changed | 13 |
| vs baseline: new | 0 |
| vs baseline: known | 1 |
| vs baseline: error | 0 |


## collisions

| Check | Ours | Measured / required | Status | vs baseline |
|---|---:|---|---|---|
| Bullet main cluster: outer stars needed for the galaxy speed (static distances) | 3.711e+12 Msun | 3.4-5.8 x 10^12 (Legacy Survey star count, round 6, converted) | pass | (was 3.711e+12, pass) |
| Bullet main cluster: lensing strength (kappa) on its galaxies | 0.6947  | at least 0.36 +- 0.06 (a stated lower bound) | pass | changed (was 0.7168, pass) |
| Bullet subcluster: lensing strength (kappa) on its galaxies | 0.2223  | at least 0.20 +- 0.05 (a stated lower bound) | pass | changed (was 0.2764, pass) |
| Bullet: leftover lensing on the main cluster's gas | 0.06387  | 0.05 +- 0.06 | pass | worse (was 0.04957, pass) |
| Bullet: leftover lensing on the subcluster's gas | 0.0837  | 0.02 +- 0.06 | pass | worse (was 0.0516, pass) |
| Bullet main: lensing peak distance from its galaxies | 27.01 kpc | within a quarter of the galaxy-gas separation (252 kpc) | pass | worse (was 15.68, pass) |
| Bullet sub: lensing peak distance from its galaxies | 49.72 kpc | within a quarter of the galaxy-gas separation (232 kpc) | pass | worse (was 22.97, pass) |
| Bullet main cluster: lensing mass inside 301 kpc (250 in LCDM units) | 2.940e+14 Msun | 3.09-3.46 x 10^14 (Paraficz 2016; Bradac 2006; converted) | pass | worse (was 2.969e+14, pass) |
| Bullet subcluster: lensing mass inside 301 kpc (250 in LCDM units) | 1.365e+14 Msun | 2.47-2.85 x 10^14 (Paraficz 2016; Bradac 2006; converted) | fail | known (was 1.385e+14, fail) |
| collision stack: how far the lensing follows the gas (beta, median of 20) | 0.02062  | -0.04 +- 0.07 (Harvey et al. 2015, 72 collisions) | pass | changed (was 0.01914, pass) |
| collision stack: lensing peak distance from the galaxies (median) | 25 kpc | raw offset; includes the pull of the other cluster | info | changed (was 15.75, info) |
| MACS J0025 SE subcluster: lensing mass inside 300 kpc (LCDM units; 429 kpc static) | 2.032e+14 Msun | 3.64 (+1.46/-2.48) x 10^14 | pass | changed (was 2.053e+14, pass) |
| MACS J0025 NW subcluster: lensing mass inside 300 kpc (LCDM units; 429 kpc static) | 1.873e+14 Msun | 3.79 (+0.73/-2.04) x 10^14 | pass | changed (was 1.853e+14, pass) |
| MACS J0025 SE: lensing peak distance from its galaxies | 44.29 kpc | with the galaxies, not the gas (530 kpc away): lensing offset from the gas at > 4 sigma | pass | changed (was 23.61, pass) |
| MACS J0025 NW: lensing peak distance from its galaxies | 212.9 kpc | with the galaxies, not the gas (244 kpc away): lensing offset from the gas at > 4 sigma | fail | REGRESSED (was 79.2, pass) |
| MACS J0025: galaxies' speed spread (line of sight, inside 1.5 Mpc) | 668.7 km/s | 835 +- 59 (108 galaxies within 1.5 Mpc) | close | (was 668.7, close) |
| Abell 520 P1, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 3.134e+13 Msun | 2.44-3.26 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | changed (was 3.158e+13, pass) |
| Abell 520 P2, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 5.793e+13 Msun | 4.70-4.83 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | better (was 6.056e+13, pass) |
| Abell 520 P3, the gas-rich, galaxy-poor centre ("dark core"): lensing mass inside 150 kpc (LCDM units) | 3.969e+13 Msun | 3.29-3.89 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | worse (was 3.717e+13, pass) |
| Abell 520 P4, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 3.952e+13 Msun | 4.91-6.48 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | worse (was 4.117e+13, pass) |
| Abell 520 P5, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 3.377e+13 Msun | 3.40-3.68 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | pass | changed (was 3.453e+13, pass) |
| Abell 520 P6, a galaxy clump: lensing mass inside 150 kpc (LCDM units) | 2.646e+13 Msun | 4.27-4.27 x 10^13 in the static law (Jee et al. 2014; Clowe et al. 2012) | close | better (was 2.527e+13, close) |
| Abell 520: lensing mass inside 710 kpc of the centre (LCDM units) | 5.261e+14 Msun | 5.84 +- 0.64 x 10^14 in the static law (Mahdavi et al. 2007; <D_ls/D_s> = 0.59) | pass | changed (was 5.272e+14, pass) |
| Abell 520: galaxies' speed spread at P1, P2, P4, P5 (rms z) | 2.697  | Girardi et al. 2008 Table 4: 811, 749, 579, 668 km/s (6-9 galaxies each) | close | (was 2.697, close) |
| El Gordo: lensing mass inside 500 kpc of the centre of mass (Planck units) | 8.447e+14 Msun | 9.45e+14 +- 12% (aperture mass) | pass | changed (was 8.423e+14, pass) |
| El Gordo: lensing mass inside 1000 kpc of the centre of mass (Planck units) | 2.138e+15 Msun | 2.43e+15 +- 12% (aperture mass) | pass | changed (was 2.142e+15, pass) |
| El Gordo NW: galaxies' speed spread (line of sight, inside 1 Mpc) | 1013.7 km/s | 1290 +- 134 (Menanteau et al. 2012) | close | (was 1013.7, close) |
| El Gordo SE: galaxies' speed spread (line of sight, inside 1 Mpc) | 903.6 km/s | 1089 +- 200 (Menanteau et al. 2012) | pass | (was 903.6, pass) |
| El Gordo SE: lensing peak distance from the cool gas core | 286.9 kpc | about 12 arcsec, the lensing peak nearer the merger centre (Kim et al. 2021; Ng et al. 2015) | info | changed (was 248.4, info) |
