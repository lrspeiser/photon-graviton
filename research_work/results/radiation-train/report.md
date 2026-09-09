# Finite radiation trains: smoother signals, but continuing drift

## Question and result

The previous shared-field calculation used one driving photon packet and found launch-dependent stretching. This comparison replaces that burst with a finite train of packets, keeping the incoming average power fixed while changing their spacing. It tests the mean redshift, gradual drift and remaining sampled variation together.

**More frequent, smaller packets reduce the sampled residual variation, but do not make the redshift steady in this finite window.** At the same nominal power, the residual RMS after a descriptive linear detrending falls from about 16.2% to 3.27% of the mean generated redshift as spacing decreases from 2 to 0.5. Those percentages are relative to the small redshift increment z, not to the entire duration factor 1+z. A substantial gradual drift remains. Doubling the power approximately doubles the mean probe redshift in the tested weak regime.

This is a deterministic finite-source comparison, not a validated cosmic radiation population, a stationary solution or an observational fit. Detrending is a diagnostic decomposition; it does not remove physical variation from a real signal.

## 1. A finite, energy-accounted source

The model retains the three cosine-squared conversion regions centered at 0, 3 and 6, each with half-width 0.5. It uses K=1, companion-wave speed v=0.5 in c=1 units and Gaussian packet smoothing width 0.25. All regions and all packets share one field. The Hamiltonian generalizes the previous packet model to

```
H = integral K/2 [u^2 + v^2 phi_x^2] dx + sum_i P_i/n_bar(X_i),
n_bar(X_i) = 1 + integral W_sigma(x-X_i) g(x) phi(x) dx.
```

The coupled equations are

```
phi_t = u,
u_t = v^2 phi_xx + sum_i [P_i/(K n_i^2)] W_sigma(x-X_i) g(x),
X_i_dot = 1/n_i,
P_i_dot = (P_i/n_i^2) partial_X_i n_i.
```

Each packet can change the propagation of later packets through the same evolving field. No independent receiving reservoir or field reset is assigned to a packet. Static spatial profiles exchange momentum without supplying energy in this approximation; their material realization remains unspecified. The numerical momentum-rate ledger includes the profile and grid reaction, not a claim that those unmodeled supports are a closed dynamical matter system.

Packets are initially positioned at x=-2-j Delta, with j starting at zero. Their initial energy is Q Delta, where Q is the nominal incoming power. There are 20/Delta packets, so the total initial photon energy is exactly 20Q. At Q=0.003 it is 0.06 for all three spacings; the Q=0.006 control begins with 0.12. The field initially has zero energy.

The initial photon train is the finite energy supply. There is no unaccounted continuous injection, infinite stellar fuel or arbitrary receiving energy term. Stellar production of this radiation is outside the calculation. The nominal interval Delta corresponds to propagation before the conversion regions in the common reference units; subsequent packet timing is evolved, not imposed.

The domain is -30<x<36 and the end time is 28. The initial trains fit inside it, and boundary-strip checks exclude relevant periodic-return contamination. This run is a finite illumination episode, not the infinite-time limit of a constant-power source.

## 2. Weak probes and sample statistics

Infinitesimal probes traverse x=-2 to x=9 at thirteen nonuniform launch times between 6 and 12.07. They sample the actual shared field but add negligible energy and therefore have no finite backreaction. As before, independently integrated carrier rays and arrival-map derivatives agree on the local frequency and duration stretch. Matter clock laws are not derived; the comparison uses declared common reference standards. Clock-scaling controls retained by the reused probe integrator are illustrative, not an adopted atomic response.

For each case, z is the probe frequency stretch minus one. The report gives the unweighted mean of the thirteen samples, a fitted linear drift against launch time, and the RMS remaining after subtracting that line. These are descriptive sample statistics, not ensemble estimates or bounds on unsampled launch times. Regular packet trains do not represent a random stellar population, and residual variation can include nonlinear startup drift as well as periodic ripple.

The finer-grid results are:

| Incoming power Q | Packet interval Delta | Packet count | Mean sampled z | Linear drift dz/dt | Residual RMS in z | Residual RMS / mean z |
|---:|---:|---:|---:|---:|---:|---:|
| 0.003 | 2.0 | 10 | 0.00548574 | 0.00024035 | 0.00088701 | 16.17% |
| 0.003 | 1.0 | 20 | 0.00530431 | 0.00031956 | 0.00022758 | 4.29% |
| 0.003 | 0.5 | 40 | 0.00521530 | 0.00033344 | 0.00017033 | 3.27% |
| 0.006 | 0.5 | 40 | 0.01043320 | 0.00066424 | 0.00034287 | 3.29% |

At fixed Q, the mean changes by about 5% across this spacing range; it is similar, not exactly unchanged. The residual RMS decreases much more strongly. However, for Delta=0.5 and Q=0.003, the sampled raw redshift spans about 0.004302 to 0.005996. Subtracting a trend must not be mistaken for a physical removal of that change. The increase in the fitted slope also shows why residual RMS alone is an incomplete measure of signal fidelity.

Doubling Q at Delta=0.5 approximately doubles both mean z and the absolute drift and residual RMS. This is evidence of dependence on the radiation background in these examples, not an independently calibrated universal redshift law. It is also not a test of individual photon color: packet energy here is a classical radiation-loading parameter.

## 3. Energy accounting

At time 28, the field energies at Q=0.003 are approximately 0.00031963, 0.00030003 and 0.00029173 for spacings 2, 1 and 0.5. They equal the corresponding total photon-energy losses within the numerical residual. At Q=0.006 and spacing 0.5, the field holds approximately 0.00115937, again taken from the finite initial photon supply.

The field energy in the double-power comparison is approximately four times as large, while the weak-probe redshift is approximately twice as large. Field energy is quadratic in its wave amplitude, so those are different quantities. Neither scaling can be extrapolated indefinitely while ignoring backreaction and the finite input energy. No gravitational amplification factor has been introduced or measured.

The companion field remains the recipient of all lost photon energy in this model. There is no capture channel, deposited mass distribution or derived extra lensing. The present calculation concerns production and propagation; it does not replace the separate entry, retention and gravitational-response requirements.

## 4. Numerical verification

An initial 1,056/2,112-cell pilot failed the predeclared absolute probe-stretch refinement threshold of 1e-4 at Q=0.003 and Delta=1. Its maximum difference was about 1.0592e-4. This partial failure is retained in [pilot-verification.json](pilot-verification.json); the threshold was not relaxed.

The final eight backgrounds use 2,112 and 4,224 cells. All 104 independent carrier/event comparisons pass. The total-energy drift is below 2.4e-14 of initial energy, numerical momentum-rate residual below 4.7e-14, carrier/event relative disagreement below 4.5e-15, and independent arrival-time difference below 9e-14 in reference-time units.

Refinement changes total photon loss by at most 0.028% relatively, sampled mean redshift by at most 0.031%, and an individual probe stretch by at most 2.45e-5 absolutely. The largest relative change in detrended RMS is about 3.12%, meeting the additional 10% RMS convergence gate. These statements establish numerical consistency of the stated finite model, not compatibility with astronomical measurements. Both grids and all un-subtracted probe values are retained.

## 5. What remains to be established

The experiment supports a limited smoothing effect from more frequent illumination. It does not show that a sustained population produces a stationary or acceptably slowly evolving background. The chosen launch window may include propagation and startup transients. A longer source history and later observation windows are needed to distinguish settling from continuing secular evolution. Mean redshift, drift and event distortion must be carried forward together.

The next comparison should extend the finite, energy-accounted source duration and follow later probes while controlling packet resolution and boundary returns. It should derive or test the mean-field evolution rather than infer stationarity from a reduced residual RMS. Random source histories, angular illumination and three-dimensional transport remain separate requirements; regular one-dimensional trains do not establish their behavior.

The source power, region widths and placement are still model inputs, not consequences of limited gravity. A complete theory must determine the environmental coupling and matter clocks, then join source production to capture, supported deposits and a common motion/lensing response. The 20-task, 32-area research goal remains active.

## Reproduction

```sh
python research_work/results/radiation-train/check.py
```

The [protocol](protocol.json) records the finite source supply, cases, observation window and checks. The [saved evidence](radiation-train-results.json) contains all packet/background cases, energy snapshots, raw probes, finite-event intervals, sample statistics and source hashes. NumPy and SciPy are required. The independent probe integrator is reused from the existing weak-signal timing module, and the standard diagnostic suite includes this calculation.
