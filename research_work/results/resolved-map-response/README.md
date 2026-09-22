# JR-8: measured galaxy shapes, spatial exchange and exact response sensitivities

21 September 2026, America/Los_Angeles. Baseline `20255778863150392baef9702b6ebb5df35d1342`. Executed locally. This is a map-conditioned construction, NOT a calibrated absolute gravity fit.

## Theory first

Ordinary matter and radiation may generate a companion state whose spatial organization, currents and history affect gravity. Here actual projected source shapes enter a three-dimensional forward/return transport solve. The same simulated source supplies force and light-deflection readouts. This is one implementation, not the definition of the wider branch portfolio.

## Completed

- 288 primary stationary 3D map-driven solves, 12 selected finer-grid solves, six larger-domain controls, and separately labeled nonlinear counterfactuals.
- The same three documented MaNGA DR17 example objects and original FITS bytes from JR-7; no new outcome-selected galaxies.
- A radial projected-velocity reconstruction on actual stellar spectral bins.
- Local rate-balance calculations over 778 electron-density shells from 12 X-COP clusters.
- Final audit: 1,562 checks, zero failures, including three independently re-solved primary states. This is not every historical repository suite.

An initial primary invocation was interrupted by the execution runtime. A first audit read only 11 of 12 refinements and correctly failed its count check. Both incomplete records remain in the full archive; no threshold was changed. All scientific computation ran locally, not on GitHub Actions.

## Inputs versus assumptions

Continuum SPX_MFLUX and Halpha EMLINE_GFLUX are required positive with signal-to-noise at least five; flagged flux pixels are excluded. Named FITS channels and original stellar BINID values are used. The custom restricted IMAGE-HDU reader exactly reproduces the prior Astropy-extracted arrays. It is not a general FITS-table reader.

Positive source shapes use six smooth radial basis functions and harmonics through m=2. The observed photometric ellipse supplies R/Re and angle under a thin-disk coordinate interpretation. Apertures are 2.00, 1.078 and 2.00 Re for 7443-12703, 8082-12704 and 8081-1901. Their beam widths are 0.308, 0.407 and 0.924 Re. Small-scale clumps are not resolved by this reconstruction.

Halpha is emissivity, not total gas density. Continuum is not bolometric source power. Each source is normalized to unit source power and unit proxy-gas inventory, so these are dimensionless shape-response calculations, NOT predicted km/s, Einstein angles, physical source masses or conversion efficiencies. No gravity coefficient is fitted to observations.

Proxy shapes proportional to Halpha and to its square root are distinct sensitivity hypotheses. Gaussian vertical lifts use gas heights 0.20/0.45 Re and source height 0.20 Re. They preserve the assigned inventory, but have NOT been jointly reprojected and fitted to images; the alternatives are not established equally plausible depth reconstructions. The exterior source is unmeasured. Masked regions use the declared smooth reconstruction, not secretly measured gas.

Observed Halpha/Hbeta, [NII]/Halpha and [SII] ratios are retained without relabeling them as density/temperature. No [SII] doublet pixels for 8082-12704 meet the joint criteria in this extraction.

## Spatial construction

The fixed JR-7 rate families are retained:

    P_t = L_P P - k_plus P + k_minus C - 0.5 P + J
    C_t = L_C C + k_plus P - k_minus C - 0.05 C
    h = rho_proxy/(rho_proxy+0.15)
    k_plus = 0.5 exp(a h)
    k_minus = 0.5 exp(b h^2)

The four (a,b) pairs are (0,0), (2,0), (0,4), (2,4). Diffusion is 0.5/0.05, carrier drift is 0 or +/-0.6, radial/vertical boundaries reflect and angle is periodic. Prescribed drift is NOT measured stellar pattern speed. These are frozen source/gas snapshots, not a self-consistent moving galaxy.

Force and deflection use the same softened 3D source, G=c=1, softening 0.15 Re. Both P+C and companion-only results are retained. This closes a two-channel exchange/escape energy ledger, not recoil, gravitational binding energy or a complete relativistic stress coupling.

## Main measured result

Compare the observed-pattern configuration with a 90-degree rotation of only the gas pattern at fixed radial inventory. Positive means the original configuration produces more than the counterfactual. The table is the FINEST mixed/no-drift result, not an observed physical change:

| Source geometry | Companion content | Mean inward force, P+C | Projected deflection, P+C |
|---|---:|---:|---:|
| 7443-12703 | +0.07655% | +0.03149% | -0.00013%, unresolved sign |
| 8082-12704 | -0.02327% | -0.00235% | +0.22962% |
| 8081-1901 | -0.06367% | -0.17044% | +0.18936% |

Force is averaged at 0.6 times the aperture radius in the adopted disk plane. Deflection is averaged on a projected impact ring at the same numeric radius, at the inclination inferred from the thin-disk axis ratio. These are different observation operators. In two cases deflection increases while the local radial force decreases, WITHOUT a separate photon multiplier. The changes are small and do not resolve historical lens residuals; these three objects do not have matched observed lensing constraints here.

A sign is numerically resolved only if its magnitude exceeds three times its change between the two finest paired grids. All listed signs except the first deflection satisfy that rule. This is not an astrophysical confidence interval. Between the two finest meshes, companion content changes by at most 0.037% and mean inward force by at most 0.230%.

## Exact spatial susceptibility

Write the stationary system A y=j, y=(P,C). For a linear observable O=l^T y solve A^T psi=l. At fixed source and transport:

    dO/d ln rho_j = (psi_C-psi_P)_j
                    * [a k_plus P - 2 b h k_minus C]_j
                    * h_j(1-h_j).

This standard adjoint mathematics is applied to the exact finite-volume operator. Both positive and negative gas-response regions occur within EVERY primary mixed object. A global increase in gas-proxy inventory has a negative derivative in all three, while some local additions have a positive derivative. Total gas amount therefore loses sign information.

All twelve companion derivative checks match independent finite differences within 1.1e-7 relative error. Six force/deflection derivatives also pass the unchanged 1e-5 criterion. Fixed-inventory sensitivities include the normalization derivative rather than silently adding mass.

## A constructive same-force/different-lensing calculation

As a post-primary feasibility calculation, project the deflection sensitivity off total gas inventory and ONE radial-force mean, then solve the full nonlinear state to remove the remaining force error. The largest of the three declared perturbations gives:

| Source | Proxy gas inventory moved | Deflection increase | Relative radial-force change |
|---|---:|---:|---:|
| 7443-12703 | 2.109% | 0.210% | 5.1e-15 |
| 8082-12704 | 1.304% | 0.306% | -2.1e-12 |
| 8081-1901 | 2.550% | 0.299% | -1.9e-15 |

These are deliberately computed counterfactual gas movements, NOT discovered configurations or altered observations. Local gas multipliers can reach about 1.8. Other radii and stellar velocities are not held fixed, and no full image likelihood is fitted. This demonstrates geometric freedom, not an observational solution.

## What dominates uncertainty

Increasing the assumed gas height from .20 to .45 Re changes companion content by +2.133%, +14.377% and +6.639%, respectively, in the primary coarse-map sensitivity. For 8081-1901 it changes the phase contrast from -0.0632% to +0.0933%. This thickness alternative is not as finely certified as the primary table and is not an independently measured height. The square-root emissivity proxy also changes the response.

Doubling the radial domain, retaining inner resolution and assuming zero exterior source/proxy gas, changes companion content by -0.545%, +5.317%, -1.091%. Corresponding phase contrasts are +0.07045%, -0.01434%, -0.05329%. Real galaxies need not end at the survey footprint; this is a boundary sensitivity, not a physical exterior model.

Thus density interpretation, vertical structure and exterior continuation are larger uncertainties than these small phase effects. Low-harmonic reconstruction omits sub-beam clumping.

## Kinematics

The matched-quality apertures contain 1,084, 106 and 176 distinct stellar spectral bins, not independent measurements at every map pixel. Gas velocities are averaged onto those bins. Median inverse variance is a working weight, not a full covariance.

In 7443-12703 the relative first-harmonic phase IN THE ADOPTED ELLIPSE COORDINATES is about +4.8, -2.0, +2.3 and +27.6 degrees at 0.67, 1.07, 1.47 and 1.87 Re. The outer structure is not represented by one global gradient angle. No formal uncertainty or 3D spin interpretation is assigned, and the simple gas harmonic fit remains imperfect. These velocities are not used as both source-current inputs and held-out force validation.

## Cluster turnover

For a local cell with fixed power J and no transport,

    C/J = k_plus/[0.5 k_minus+0.05 k_plus+0.025].

For the mixed law its density derivative switches sign at

    (8h-2) * 0.5 exp(4h^2) = 0.1
    h = 0.26872794856
    n_e/n0 = 0.36748013004.

Below the threshold more density increases content per unit source; above it, return dominates that derivative. If source power itself varies with density/temperature, total production can differ. The maximum local content per unit source is 27.888% above its zero-gas-rate value; this is not a bound on every spatial theory.

Using all 778 X-COP density shells, illustrative n0=1e-3 cm^-3 puts the crossings at 644-1,029 kpc in the twelve clusters. For n0=1e-2 cm^-3 they move to 44-235 kpc where present; one cluster does not cross within the supplied shells. These are log-interpolated, conditional positions, NOT observed conversion fronts. The reference density is uncalibrated; no cluster source, temperature dependence, velocity field or lensing fit is supplied. Hydrostatic and NFW/Einasto masses are excluded.

## Reproduction and limits

The full companion archive contains original FITS/NPZ maps, cluster profiles, fields, all 288 cases, refinements, counterfactuals, hashes, incomplete execution records, 1,562-check readback and report-generation utilities. Scientific scripts and a compact result record are committed here. Full arrays are supplied in the originating conversation archive, not asserted permanently hosted by this Git tree.

With package inputs and fresh output paths:

    python code/maps.py --inputs inputs/maps --output fresh-maps
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/response.py --maps fresh-maps/maps.json --output fresh-response
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/response.py --maps fresh-maps/maps.json --output fresh-response --refine-only
    python code/cluster_response.py --input inputs/clusters/cl2-inputs-xcop-profiles.json --output fresh-clusters.json
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/supplement.py --maps fresh-maps/maps.json --output fresh-supplement.json

Saved-output readback: use a copy of the archive, preserve the old results/audit.json under another name, then run `python code/audit.py --root COPY`.

The completed link is observed shape -> assumed 3D source -> spatial exchange -> common force/lensing -> exact susceptibility. Absolute source energy, true gas density, full depth/current inference, stellar dynamics, observed lens images, Solar-system recovery, merger history, and a microscopic gravitational coupling remain open. These gaps are not filled by selecting convenient geometry after inspecting gravitational residuals.

Primary documentation: official MaNGA DAP data model https://sdss-mangadap.readthedocs.io/en/latest/datamodel.html ; official DR17 access https://www.sdss4.org/dr17/manga/manga-data/data-access/ ; prior JR-7 and the pinned X-COP extract in `research_work/results/path-memory/cl2-inputs-xcop-profiles.json`.
