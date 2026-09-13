# Present radiation does not uniquely determine accumulated deposits

## Constructive result

Within the unchanged homogeneous conversion/capture equations, positive finite-duration emission histories can have exactly the same present photon energy and traveling companion energy but substantially different deposited energy. Even measuring both current radiation reservoirs would not uniquely identify the past stored amount.

| History | Present photon energy | Present companion energy | Deposited energy | Total past emitted energy |
|---|---:|---:|---:|---:|
| Reference | 1.272738 | 0.691277 | 0.035985 | 2.000000 |
| Alternative with older burst at alpha c T=10 | 1.272738 | 0.691277 | 0.730888 | 2.694903 |
| Alternative with older burst at alpha c T=30 | 1.272738 | 0.691277 | 8.582205 | 10.546220 |
| Alternative with older burst at alpha c T=100 | 1.272738 | 0.691277 | 9957.267845 | 9959.231860 |

All energies use the same arbitrary physical energy-density unit. These histories are synthetic mathematical constructions, not inferred stellar populations. No age in years or universe size is stipulated. The first alternative stores 20.31 times the reference energy; the additional stored energy is fully paid for by additional past emission, not energy creation.

## Construction and verification

Set beta/alpha=0.1, one of the preceding transport comparison cases. For dimensionless pulse age x=alpha c T, the known impulse responses are

    P(x)=exp(-x)
    C(x)=[exp(-x)-exp(-0.1x)]/(0.1-1)
    D(x)=1-P(x)-C(x).

Replace impulses by constant positive bursts of width 0.02 in x. The finite-burst response is the average of these kernels across each burst interval. The reference emits one unit at each age 0.1 and 1. The alternative changes these two positive burst amplitudes and adds an older positive burst. Solving the two linear equations for equality of current P and C leaves a free history direction; choose 90% of the maximum old-burst amplitude before one of the recent amplitudes reaches zero. The resulting recent amplitudes remain positive.

This is linear-algebraic nonuniqueness of the inverse problem, using established transport mathematics. It is not a new microscopic law or a fitted explanation of the observed universe. A common continuing source may be added to every history without changing their equality of present radiation or difference in stored energy; zero current emission in the displayed burst example is not essential to the construction.

An augmented matrix exponential integrates each finite burst and free propagation. Direct quadrature of the analytic pulse kernels independently agrees within 1e-12 per unit burst. The code checks nonnegative amplitudes/reservoirs, identical current radiation to absolute 1e-12, normalized pulse energy, and total emitted energy equal to the sum of present reservoirs. Complete amplitudes and results are saved in history-results.json.

## What follows and what does not

The current limited galaxy inventory or bolometric photon background alone cannot establish an upper bound on past deposits if arbitrary past source energy is allowed. More time permits old energy to have left both traveling reservoirs. Consequently the earlier steady relation P_receiver=sigma J/beta is a conditional current rate, not a general reconstruction of accumulated energy from today's J.

These alternatives are not indistinguishable under every observation. Equal bolometric photon energy does not imply equal spectra, source counts, stellar ages, remnants, element abundances, geometry, time variations or gravitational fields. The required stellar fuel is very different. No claim is made that stars can actually realize the large old burst or that its deposits remain supported as assumed.

The exact ledger supplies the next necessary constraint:

    E_deposited = E_total_emitted - E_photons - E_companions

for the closed homogeneous bookkeeping with no initial reservoir. A physically supported upper bound on total emitted energy would therefore bound the deposits. Removing a fixed universe age does not remove the need to account for fuel, remnants and past luminosity. For open regions, boundary fluxes and initial energy must also enter this ledger.

Next constrain admissible stellar histories and source energy using ordinary-matter inventories and spectral/stellar-population observations under the proposed propagation law. Then evolve transport and deposit dynamics before comparing lensing and motions. Arbitrarily adding old light until gravity matches would not constitute a prediction. All six objectives remain open; no observational holdouts were opened.

Run `python research_work/results/cluster-competing-capture/history.py`.
