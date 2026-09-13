# Boundary-capable timing estimator: development refit

This is a selected-case development study on all 16 invalid cases from the completed frozen simulation batch. No real fluxes or untouched observations are used. Original estimates and coverage remain unchanged.

## Method and formula provenance

The population model log W ~ Normal(a + b log(1+z), sigma) is an existing statistical model, not a new physical law. We integrate a piecewise-linear event likelihood against its normal population density using known Gaussian interval integrals. At sigma=0 the integral becomes evaluation at the population mean. The proposed change is numerical/statistical, not a derivation of photon conversion.

Mean-duration bounds remain 5-100 days and b remains -2 to 3. Scatter now includes zero and retains upper bound 0.6. The normal distribution remains normalized on the original finite width support of 2-256 days; an explicit constraint retains the original maximum outside probability of 1e-6 for every event. No narrower mean/exponent domain is substituted. The acceptance check permits 1e-11 absolute numerical tolerance on that probability; all selected best fits have computed outside mass zero.

Five SLSQP starts per case are retained, including a zero-scatter start. The best feasible successful start is reported; optimizer success is not proof of global optimality.

| Case | Original b | Revised b | Revised scatter | Successful starts | Objective spread across successful starts |
|---|---:|---:|---:|---:|---:|
| split_gaussian-snr5-b1-seed902 | 1.030032 | 1.030900 | 1.68078e-08 | 5/5 | 0.00660573 |
| split_gaussian-snr5-b0-seed903 | -0.043089 | -0.043128 | 0.0288319 | 5/5 | 0.0221713 |
| split_gaussian-snr5-b0-seed904 | 0.005378 | 0.005405 | 0.0289823 | 5/5 | 0.0146321 |
| split_gaussian_shoulder-snr5-b0-seed904 | -0.012734 | -0.011031 | 0.00458509 | 5/5 | 0.000473188 |
| split_gaussian_shoulder-snr5-b1-seed904 | 0.905505 | 0.905148 | 0.0226619 | 5/5 | 0.0217325 |
| split_gaussian-snr5-b0-seed905 | 0.009280 | 0.018650 | 4.84801e-23 | 5/5 | 0.000153953 |
| split_gaussian-snr5-b1-seed905 | 0.895701 | 0.903914 | 3.54756e-08 | 5/5 | 0.00102211 |
| split_gaussian_shoulder-snr5-b0-seed905 | -0.055064 | -0.053809 | 0.0152918 | 5/5 | 0.00147991 |
| split_gaussian-snr5-b0-seed909 | -0.066992 | -0.062662 | 6.20743e-25 | 5/5 | 0.00238841 |
| split_gaussian-snr5-b1-seed909 | 0.944102 | 0.947316 | 1.13263e-23 | 5/5 | 0.00171269 |
| split_gaussian_shoulder-snr5-b0-seed909 | -0.091785 | -0.087484 | 4.12637e-09 | 5/5 | 3.23382e-05 |
| split_gaussian-snr5-b1-seed912 | 1.044052 | 1.038886 | 9.52193e-07 | 5/5 | 1.5973e-11 |
| split_gaussian-snr5-b0-seed916 | -0.122220 | -0.130584 | 6.39827e-24 | 5/5 | 4.14425e-05 |
| split_gaussian-snr5-b0-seed919 | 0.144146 | 0.147601 | 6.82573e-30 | 5/5 | 0.00208057 |
| split_gaussian-snr5-b0-seed920 | -0.273140 | -0.274231 | 8.8518e-23 | 5/5 | 1.27552e-08 |
| split_gaussian_shoulder-snr5-b0-seed920 | -0.233301 | -0.233655 | 2.07087e-09 | 5/5 | 1.3604e-07 |

## Outcome and remaining work

All 16 cases return a feasible candidate; 11 have scatter below 1e-6. Maximum exponent change from the original is 0.00936951. The maximum successful-start objective spread is 0.0221713.

The prior integration verification passes 35 comparisons to independent quadrature (maximum absolute error 4.45e-16); it establishes only the declared interpolant integral. It does not establish underlying event-grid accuracy, global fitting convergence, or confidence-interval coverage. Small-scatter interpolation corners and optimizer-start differences remain unresolved. Further optimization comparisons and grid checks must precede adopting this method. Then recalibrate all 160 simulations, including fixed-truth profiles, and test new simulations. Do not use this selected sample to claim repaired coverage or redefine the original failed gates.

No physical conclusion about redshift or companions follows from these synthetic refits. All six scientific goals remain open.
