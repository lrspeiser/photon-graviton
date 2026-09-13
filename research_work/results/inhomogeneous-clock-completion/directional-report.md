# Directional requirement for the evolving lapse candidate

The preceding goal turn compared stationary conversion with the published supernova timing summary. This continuation examines the already implemented evolving, nonuniform clock candidate rather than introducing another unconstrained redshift coefficient. All six objectives remain open.

## Established model and deduction

Retain the existing postulates: fixed spatial rulers, clocks at rest obeying d_tau=dt/n, and nondispersive Hamiltonian omega=c|k|/n. These are known lapse-metric kinematics used as hypothetical physical rules. Prior work already shows homogeneous cancellation and nonzero measured stretching in prescribed evolving nonuniform fields. This is not a claim that those earlier results are new.

Along a ray with unit spatial direction k_hat, Hamiltonian evolution gives d log(omega)/ds=-partial_t n/c. A local clock measures Omega=n omega. Since dt/ds=n/c, ordinary differentiation yields

    d log(Omega)/ds = k_hat dot grad(log n),
    log(S) = -integral k_hat dot grad(log n) ds.

The time derivative cancels locally between the wave and clock, but spatial gradients sampled at different times can still produce a finite path effect. The existing model identifies S with wavelength and event-duration stretch. This generalizes its one-dimensional identity; it is a conditional calculus result, not a new fundamental law.

Define the local logarithmic redshift rate alpha_local=-k_hat dot grad(log n). At the SAME spacetime event,

    alpha_local(-k_hat)=-alpha_local(k_hat),
    angular_mean(alpha_local)=0,
    angular_mean(alpha_local^2)=|grad(log n)|^2/3.

Therefore this clock rule alone cannot supply a strictly positive, direction-independent local rate alpha. Replacing the dot product with a gradient magnitude would change the physical propagation rule, not simply calibrate its coefficient. Our earlier phenomenological isotropic conversion coefficient cannot automatically be identified with this local lapse-gradient rate.

## Why averaging redshift itself can mislead

As an angular benchmark, prescribe log(S)=-x mu, with directions uniformly distributed in mu between -1 and 1. Here x is an integrated-gradient parameter. This is a specified shift distribution, not a solved three-dimensional ray field; straight paths through a finite gradient need not obey the full bent-ray dynamics.

Known angular integration gives mean(log S)=0, while mean(z)=sinh(x)/x-1>0. Half the directions blueshift. Thus a positive arithmetic mean redshift alone does not establish the desired predominantly redshifted distance trend.

| Target arithmetic mean z | x | RMS log(S) | Blueshift fraction |
|---:|---:|---:|---:|
| 0.001 | 0.07745 | 0.04471 | 50% |
| 0.01 | 0.24458 | 0.14121 | 50% |
| 0.1 | 0.76340 | 0.44075 | 50% |
| 1 | 2.17732 | 1.25708 | 50% |

For example, mean z=0.01 comes with extrema approximately -0.217 and +0.277 in this benchmark. These are synthetic predictions, not measured sky statistics. No real-data likelihood or observational exclusion is claimed for this benchmark.

## Consequence and limits

This local directional obstruction does NOT prove zero cumulative redshift in every evolving, statistically isotropic inhomogeneous universe. Opposite observed rays traverse different spacetime points; their path histories and source/observer environments can correlate with the field. Those correlations could generate a nonzero integrated mean and must be calculated, including scatter and angular dependence. Source and observer motion also changes the operational prediction.

The next adequate candidate must generate its field and actual ray paths, calculate endpoint clock readings, and predict the DISTRIBUTION of redshift and timing across directions. It cannot claim an isotropic local positive conversion rate from the present lapse rule without an additional interaction or a changed clock/propagation law. This is a substantive requirement for an all-sky explanation, not an exclusion of nonexpanding physics.

directional.py and directional-results.json verify the analytic angular means, variance, opposite-direction wavelength-ratio product, and blueshift fraction with 256-point quadrature in four cases. These checks validate this deduction, not the existence of the field. Existing radiation-powered initiation also settles rather than maintaining the needed redshift; no autonomous sustained solution is established. No observed flux or final holdout was accessed.
