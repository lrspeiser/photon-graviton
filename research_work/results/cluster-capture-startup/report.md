# Finite travel time in all-direction cluster capture

## Result and scope

The external companion bath can feed a large spherical receiver from every direction, but intercepted energy is not immediately deposited. A finite-flight calculation now separates deposited, still-traveling and escaped energy. This closes the startup omission identified in the earlier isotropic-capture-profile report for a simpler constant-opacity sphere. It does not fit an observed cluster or derive its photon supply.

## Assumptions and formula provenance

Project postulates for this diagnostic: fixed spherical boundary of radius R; straight rays moving at c; constant absorption coefficient kappa; an initially empty interior; a constant isotropic incoming specific intensity I0 switched on at boundary time zero; permanent stationary capture; no internal sources. These are optional effective assumptions, not claimed original physical laws. Fixed geometry neglects focusing, mechanical response and gravitational work; the ledger below is the transport-sector energy balance, not full dynamical stress-energy closure.

Known angular-flux geometry gives inward incident power

    P_in = 4 pi^2 R^2 I0 = pi R^2 c u_boundary,
    u_boundary = 4 pi I0/c.

Here u_boundary describes the unattenuated isotropic bath, not the interior density. Zero net vector flux of an isotropic bath does not mean zero inward incident power. The source maintaining the imposed bath must eventually be included in the global energy ledger.

At radius r and propagation cosine mu the backward path to the boundary is

    ell = r mu + sqrt(R^2-r^2+r^2 mu^2).

Known absorption and causal travel time then give the conditional accumulated deposit density

    u_D(r,T) = 2 pi kappa I0 integral[-1,1]
               exp(-kappa ell) max(T-ell/c,0) dmu.

This is a derivation from the stated assumptions, not a new fundamental law. At early times the interior has not yet been illuminated. For T >= 2R/c the instantaneous capture rate reaches its steady value, but the accumulated deposit retains the startup deficit.

An independent ray calculation uses chords L=2R sqrt(1-b^2/R^2) with incident area weight 2b db/R^2. Along each chord, m=min(L,cT), and the deposited, traveling and escaped energies per unit incident power are respectively

    integral[0,m] kappa exp(-kappa s) (T-s/c) ds,
    integral[0,m] exp(-kappa s) ds/c,
    max(T-L/c,0) exp(-kappa L).

Their sum is T. Averaging over chords therefore conserves every unit of incoming energy without double counting it as both traveling and captured.

## Executed results

Fractions below are relative to all energy that has entered by T. Time is expressed in R/c; no universe age or physical cluster size is imposed. Optical depth means kappa R. The last column compares deposits with the shortcut steady-capture-rate times T.

| Optical depth | T/(R/c) | Deposited | Traveling | Escaped | Deposit / steady shortcut |
|---|---:|---:|---:|---:|---:|
| 0.1 | 1 | 0.04641 | 0.87430 | 0.07929 | 0.37476 |
| 0.1 | 2 | 0.07886 | 0.61923 | 0.30191 | 0.63676 |
| 0.1 | 10 | 0.11485 | 0.12385 | 0.76131 | 0.92735 |
| 1 | 1 | 0.35621 | 0.59197 | 0.05182 | 0.50670 |
| 1 | 2 | 0.51316 | 0.35150 | 0.13534 | 0.72996 |
| 1 | 10 | 0.66503 | 0.07030 | 0.26466 | 0.94599 |
| 10 | 1 | 0.89650 | 0.09950 | 0.00400 | 0.90101 |
| 10 | 2 | 0.94575 | 0.04975 | 0.00450 | 0.95050 |
| 10 | 10 | 0.98515 | 0.00995 | 0.00490 | 0.99010 |

For example, at optical depth 10 and T=2R/c, about 94.6% of incident energy is deposited, 5.0% remains in transit, and 0.45% has escaped. Strong capture can retain most incoming energy under these assumptions. This is not evidence that the incoming intensity is large enough for any actual cluster.

## Verification and consequences

The 15 cases compare independent chord-energy integration against spatial/angular integration at two resolutions (256/512 radial nodes and 512/1024 angular nodes). Maximum absolute volume/chord fraction difference is 1.29e-07; maximum ledger error is 2.22e-16. All deposited fractions remain below the steady shortcut. These checks validate the conditional transport calculation, not capture physics or a gravitational response.

To obtain an actual cluster prediction, replace the prescribed bath with photon-funded, direction- and time-dependent supply; include internal emitters without double counting boundary inflow; derive capture and subsequent reservoir motion; and calculate lensing and matter dynamics from the same source. A larger radius increases intercepted power for the same bath, but does not determine the bath or the gravitational distribution. No observational holdouts were opened, no model parameters were tuned to cluster data, and all six research objectives remain open.

Run: `python research_work/results/cluster-capture-startup/run.py`.
