# Cross-scale status of MOND-guided redistribution

13 September 2026. Treat the SPARC improvement as a candidate branch, not a global replacement. Freeze f=0.9235570945357928 and a0=8.563335193921255e-11 m/s^2. Do not reselect them for each system.

## Milky Way: executed frozen transfer

All 18 archived combinations of two ordinary-matter baselines, three disk sizes and three luminosity factors were transferred without fitting any Milky Way speed. Compare the same 38 Eilers circular-speed summaries, about 5-25 kpc. These are kinematically inferred circular speeds, not individual measured stellar orbit speeds; keep their existing baseline and systematic limitations.

| Model | Baseline I RMS km/s | Baseline II RMS km/s |
|---|---:|---:|
| Original exact-third | 6.761 | 10.410 |
| Simple MOND with frozen SPARC a0 | 9.526 | 12.044 |
| MOND-guided fixed-inventory mixture | 9.259 | 11.870 |

Both fiducial baselines worsen. Thirteen of the 18 declared sensitivity cases improve, but choosing one afterward would change the input assumptions to favor the result. Keep all cases and the original fiducials. Neither main baseline requires an inventory cap. No universal transfer success is established. The inherited reference quadrature was already refined to below 0.001 km/s changes; the algebraic mapping and nonnegative capped endpoint are checked directly here. No new physical capture solver or full barred/vertical Milky Way model is implied.

The calculation uses the archived reference total inventory and motion predictions, applies the same MOND-derived target mass and running-maximum/cap rule, then predicts the same 38 speeds. Results retain every sensitivity and prediction in mond-cross-scale-results.json. This is a real-data transfer, not a refit or a new independent dataset.

## Clusters: updated calculation specification, no new observational fit

The archived six-bin Coma comparisons are freely normalized projected-shape fits. They do not supply a matched extended ordinary gas-and-star mass model or a fixed galaxy-law companion inventory for Coma. Their point-baryon MOND control must not be promoted to a realistic cluster baseline. A new fully normalized transfer cannot be honestly reported from those shape fits alone.

Apply the candidate only after specifying M_b(<r) from extended gas and stars and the original physically defined companion capture inventory. In spherical approximation, use gb=G M_b/r^2, calculate the known MOND excess, and use exactly the same f, a0 and cap. Do not create a new cluster-specific disk luminosity/scale proxy without identifying the added assumption. Show baseline Newtonian, original companions, MOND-guided companions, galaxy-calibrated MOND and independently fitted halo comparisons with their differing fit freedoms disclosed.

For lensing, specify a complete density outside the observed region and calculate both local and mean projected surface density. Compare reduced shear when that is the observed quantity. Moving spare mass into an arbitrary outer shell can alter lensing even if interior rotation is unchanged; the SPARC shell-existence construction does not fix the answer. Do not omit the gas, fit an amplitude to restore agreement or hide cap-induced shortfalls. The next-sample specification now carries these conditions explicitly. No cluster result from this branch is claimed.

## Redshift: maintain separate conversion and redistribution modules

The accumulated conversion postulate remains dE_gamma/ds=-alpha E_gamma and 1+z_conv=exp(integral alpha ds). The exponential loss solution is known mathematics; its interpretation as a cosmological frequency shift caused by companion conversion is hypothetical. Redistribution supplies no new equation for alpha and no new event-duration prediction.

Known local conservation bookkeeping can be written:

    partial_t u_gamma + div F_gamma = j_star - Q
    partial_t u_travel + div F_travel = Q - C + R
    partial_t u_deposit + div F_deposit = C - R

Q is conversion power per volume; C is capture; R is release. Redistribution changes the deposited flux F_deposit, while the internal transfers cancel on summing the equations. These equations are bookkeeping constraints, not a derived interaction. Boundary fluxes and any field/binding or kinetic energy released during settling must be included in the complete energy ledger. No continuous creation of energy is authorized by fitting a new density profile.

With alpha and photon paths fixed, changing the deposited distribution alone leaves the conversion-only redshift predictions unchanged. Do not claim that it improves the existing redshift score or equate MOND a0 with alpha: their units and physical roles differ. The saved 164-group comparison still gives cross-validated c*z-residual RMS 448.38 km/s for constant conversion versus 448.09 for the linear control and 452.52 for the tested smooth rate. These are redshift residual units, not measured peculiar speeds.

A changed gravitational potential may affect paths, endpoint gravitational redshifts and timing in a future coupled field model. Those changes have not been calculated here. The earlier result that stationary same-frequency photon/graviton mixing can transfer energy without producing spectral redshift remains relevant: improving deposited gravity does not resolve conversion microphysics or supernova event stretching.

## Working decision

Retain original and MOND-guided branches side by side. Require frozen-parameter transfer before promoting the new branch. Its SPARC success and fiducial Milky Way degradation must both be included in the paper. Prioritize an acceleration-dependent transport/retention mechanism and a specified outer profile; a cluster fit or redshift change cannot be inferred from the SPARC improvement alone. The v1.5 PDF predates these results.
