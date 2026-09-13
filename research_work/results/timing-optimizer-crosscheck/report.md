# Timing optimizer cross-check

All 16 exposed boundary cases were checked with derivative-free Powell optimization, starting at each previously selected SLSQP fit. Each received a free-scatter run and an explicit sigma=0 run. Original parameter bounds and outside-population probability limit were retained. No base calibration estimate or coverage classification was changed.

Maximum improvement in negative log likelihood: 9.59161639e-09. Maximum absolute exponent change: 8.15928907e-08.

| Case | Objective improvement | Exponent change |
|---|---:|---:|
| split_gaussian-snr5-b1-seed902 | 9.59161639e-09 | 2.99316127e-12 |
| split_gaussian-snr5-b0-seed903 | 1.27897692e-13 | -1.23999525e-08 |
| split_gaussian-snr5-b0-seed904 | 1.13686838e-13 | -8.93362239e-09 |
| split_gaussian_shoulder-snr5-b0-seed904 | 2.34479103e-13 | -8.20116135e-09 |
| split_gaussian_shoulder-snr5-b1-seed904 | 1.84741111e-13 | -1.26968083e-08 |
| split_gaussian-snr5-b0-seed905 | 1.85451654e-12 | -2.80988219e-11 |
| split_gaussian-snr5-b1-seed905 | 2.84379098e-09 | -8.15928907e-08 |
| split_gaussian_shoulder-snr5-b0-seed905 | 6.60804744e-13 | -1.24766738e-08 |
| split_gaussian-snr5-b0-seed909 | 2.15058549e-09 | 9.18616572e-12 |
| split_gaussian-snr5-b1-seed909 | 0 | 0 |
| split_gaussian_shoulder-snr5-b0-seed909 | 2.08444817e-10 | 9.16963311e-11 |
| split_gaussian-snr5-b1-seed912 | 3.55271368e-14 | 2.1684512e-09 |
| split_gaussian-snr5-b0-seed916 | 3.33955086e-12 | -3.33270078e-11 |
| split_gaussian-snr5-b0-seed919 | 8.71743566e-10 | -1.64224879e-08 |
| split_gaussian-snr5-b0-seed920 | 2.02241779e-10 | 7.60486618e-10 |
| split_gaussian_shoulder-snr5-b0-seed920 | 3.25897531e-10 | 3.73437392e-12 |

This supports local stability of the selected candidates under a second optimization method. Both methods start near the same candidate; this is not an independent global search and does not resolve every alternative start. The earlier spread across SLSQP starts remains recorded. No confidence-interval calibration follows from this check.

The next experiment applies the revised estimator to all 160 original artificial samples with fixed-truth profile likelihoods. It is development recalibration on exposed samples, not a blind test. New seeds and event-grid sensitivity remain required before observational inference. All six scientific goals remain open.
