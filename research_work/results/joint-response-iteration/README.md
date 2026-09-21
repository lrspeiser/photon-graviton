# JR-1: real-data iterative companion response

**21 September 2026. Twelve candidates executed, with 38 optimizer starts. Substantial improvement; not a completed or observationally validated theory.**

## Theory first

The working hypothesis is that matter and radiation generate an extended companion state whose gravitational response affects stellar motion and lensing, without appreciably disturbing internal planetary orbits. JR-1 fits an effective, source-linked spatial closure for that state. It does not derive microscopic production, its energy supply, formation history, or the nested Solar-system solution. A fitted envelope is not evidence that a physical source has produced it.

## Executed work and provenance

Scientific baseline: c22d188949874ea0ef05b7e782d341973d421385. Initial declaration: 7267c709ace23af6afa56cadd04688c85285cf12, [PROTOCOL.md](PROTOCOL.md). Residual-driven revisions are recorded in [AMENDMENT-1.md](AMENDMENT-1.md), [AMENDMENT-2.md](AMENDMENT-2.md), [AMENDMENT-3.md](AMENDMENT-3.md), and [AMENDMENT-4.md](AMENDMENT-4.md).

The fitting and numerical checks ran locally on CPython 3.13.5, NumPy 2.3.5, SciPy 1.17.0. GitHub Actions run 35663965932 exported public input files only; it was not the scientific fit. Its artifact SHA-256 is 170d3632038db72f34f70d0031f086508ecf6b7cf42a26f0903be93825fcc6ce.

R0 is a fixed initial envelope. R1-R9 contain 30 optimizer starts. R9 was selected before scoring the 31 SPARC comparison galaxies or two out-of-fit lenses in this campaign. Their scores were then exposed and preserved. All systems were historically examined, so even that first transfer is not new blind confirmation.

R10 and R11 are explicitly post-transfer development, each with four optimizer starts. Their optimization still uses the original training subset, but their formulation followed exposed transfer residuals. Do not relabel their comparison scores as independent validation. R9's selection and predictions are not overwritten. A further ordinary-matter control was fitted from three starts with matched R9 stellar/orbit nuisance freedom.

## Real observational inputs

SPARC: 149 galaxies and all 3,152 positive-radius raw rotation measurements, preserving 89 training / 29 development-validation / 31 comparison roles. These primarily gas-traced rotation curves are not 3,152 independent stellar measurements. Ordinary radial-force components are from the catalog mass models. Source stellar disk mass is separately integrated from the measured surface-brightness profile and its stated continuation; gas mass is H I times 1.33. The bulge retains the archived spherical reconstruction from the baryonic bulge contribution. Observed rotation speeds are not used to construct those source masses. Velocities and their uncertainties were not changed or assigned a fitted error floor.

SLACS: six systems with usable resolved stellar V_rms, full within-object covariance, fixed published light components, conditional population masses, and catalog Einstein-angle summaries. The 40 radial bins split into 27 fitted bins in J0037-0942, J1112+0826, J1204+0358 and J1402+6321, and 13 out-of-fit bins in J1621+3931 and J1630+4520. J0330-0020 remains excluded by the existing use flag; J1538+5817 lacks detailed light components. Missing/excluded systems are not successes.

The lens calculations retain static-Euclidean conditional distances and the energy-loss-plus-event-stretch normalization of Chabrier population masses. These are model-dependent conversions, not independent distance measurements or a new stellar-population inference. Sersic profiles are deprojected spherically, neglecting flattening, with the released 0.8-arcsecond PSF and annular covariance. The SIE Einstein angles are image-model summaries, not a full image likelihood. Their full errors and covariance are unavailable in these inputs.

## The connected spatial model

The extra inward force is

    g_chi(r) = (A/r) * (r/rc)^q/[1+(r/rc)^q] / sqrt[1+(r/rt)^2].

For q>0, the implied enclosed effective mass r^2*g_chi/G increases monotonically and the total is finite, A*rt/G. This is a proposed spherical companion envelope, including around disks whose baryonic radial force remains the catalog disk force. It is not a full source-generated three-dimensional state, AQUAL solve, or off-plane disk prediction.

R0-R10 use equal temporal/spatial potentials. The same g_b+g_chi feeds the annular Jeans/seeing calculation and the weak-field light-deflection integral. There is no independent per-lens gravity multiplier. R11 allows one universal gamma_chi, with U_t=Phi_b+chi and U_s=Phi_b+gamma_chi*chi. Its light integral uses g_b+(1+gamma_chi)*g_chi/2; stellar motion still uses g_b+g_chi. This explicitly relaxes the original equal-potential assumption, but is not a derived relativistic action.

For R10 let M_d and M_s be inferred stellar disk and spheroid masses, M_g be gas, and R_e the photometric half-light radius. Define

    M_b = u_star*(M_d+u_s*M_s)+M_g
    S = u_star*(M_d+w_s*u_s*M_s)+M_g
    f_s = u_s*M_s/(M_d+u_s*M_s)
    C = [M_b/(1e10 Msun)]/(R_e/kpc)^2
    A = Aref*(S/1e10 Msun)^p*C^dA
    rc = c0*R_e*C^dc*10^(-k_s*f_s)
    rt = t0*R_e
    q = q_disk+(q_spheroid-q_disk)*f_s.

S is a production proxy, not a derived physical energy injection rate. The source parameters are shared across galaxies, including the same spheroid efficiency for SPARC bulges and spherical SLACS systems. Stellar nuisance offsets enter baryonic force, source proxy and both observables consistently. No object is assigned an independently fitted companion amplitude.

R10 has 20 fitted parameters: ten field/source coefficients and ten stellar/orbit nuisance parameters. The nuisance parameters are two shared population scales, four individual stellar-mass offsets and four orbital anisotropies, with the declared priors/bounds. The out-of-fit lenses use no individual mass correction and the mean fitted beta as a point prediction. Full population uncertainty is not marginalized. J0037-0942 is at beta=0.35; positive Jeans moments alone do not establish a nonnegative stable distribution function.

R10 values: Aref=10091.75 (km/s)^2 approximately; c0=10^0.4069629462; t0=10^1.1596824779; p=0.526234245; dA=-0.064663394; dc=0.225231311; q_disk=1.890635427; q_spheroid=1.676008720; log10(w_s)=0.732739288; k_s=1.652351040; log10(u_star)=-0.021265419; log10(u_s)=-0.108897934. Use the exact machine-readable parameters, not rounded numbers, for replay. The mass-normalization follow-up traded a lower spheroid stellar scale against larger companion efficiency; this is a degeneracy, not an independently measured change of stellar population.

## All twelve revisions

SPARC column: mean per-galaxy velocity RMSE over the 29 development-validation galaxies. Stellar column: mean fractional RMS over four fitted lenses. Angle column: fractional RMS over those four Einstein angles. These are descriptive metrics, not a joint acceptance probability.

| Revision | Parameters | SPARC validation RMSE km/s | Fitted stellar error | Fitted angle error |
|---|---:|---:|---:|---:|
| R0 initial reservoir | 0 | 19.992 | 39.408% | 60.380% |
| R1 amplitude/core | 2 | 46.891 | 19.263% | 31.436% |
| R2 source compactness | 4 | 44.161 | 8.607% | 16.518% |
| R3 spatial shape | 7 | 44.753 | 7.086% | 16.157% |
| R4 shared stellar/orbits | 9 | 39.396 | 8.619% | 19.025% |
| R5 production/release shape | 10 | 39.419 | 8.629% | 19.038% |
| R6 radial orbits | 12 | 38.509 | 8.352% | 20.383% |
| R7 source geometry | 11 | 18.215 | 3.373% | 8.392% |
| R8 geometry/rate/orbits | 14 | 18.079 | 3.528% | 7.396% |
| R9 stellar priors/shape | 19 | 18.210 | 1.658% | 5.311% |
| R10 spheroid population | 20 | 17.684 | 1.563% | 5.257% |
| R11 universal metric ratio | 21 | 17.687 | 1.575% | 5.264% |

The important early failure was a tradeoff: R1-R6 improved lenses while worsening disk curves. Geometry-dependent source efficiency and core response removed much of that tradeoff. Later stellar/orbit nuisance freedom improved lens kinematics. The R11 optimum gamma_chi=1.005827 is nearly equal-potential and adds negligible training improvement while slightly worsening the original validation-plus-fitted-lens criterion. Keep R11 as evidence; retain the simpler R10 as the latest substantive development candidate. No uncertainty on gamma_chi was established.

The objective gives separate equal-weight mean standardized-square blocks to rotation, lens kinematics and deflection; objects are equal-weighted within a block. Kinematics use supplied errors/covariance. Lensing uses a declared 5-percent working deflection scale for optimization, NOT an observed uncertainty. Prior penalties are included in fitting. The R9 selection criterion has no formal model-complexity penalty. Neither optimizer convergence nor this balanced objective supplies a p-value or evidence ratio.

## Preserved R9 transfer and R10 follow-up

| Sample / metric | Initial R0 | Frozen R9 | Post-transfer R10 |
|---|---:|---:|---:|
| 89 training disks, RMSE km/s | 16.948 | 18.710 | 18.633 |
| 29 development-validation disks, RMSE km/s | 19.992 | 18.210 | 17.684 |
| 31 comparison disks, RMSE km/s | 14.393 | 15.258 | 15.089 |
| Four fitted lenses, stellar fractional RMS | 39.408% | 1.658% | 1.563% |
| Four fitted lenses, angle fractional RMS | 60.380% | 5.311% | 5.257% |
| Two optimizer-excluded lenses, stellar fractional RMS | 36.020% | 3.208% | 2.660% |
| Two optimizer-excluded lenses, angle fractional RMS | 60.560% | 13.570% | 11.567% |

Not every metric improves. The initial reservoir has lower comparison-disk velocity RMSE than the joint models. R10 has 22 of 31 comparison rotation curves below 20-percent fractional RMS and 12 below 10 percent. Difficult objects are retained. R10/R11 comparison numbers are post-transfer development, not independent validation.

### Actual R10 predictions

| Galaxy | Role | Observed Einstein arcsec | Predicted arcsec | Angle error | Stellar fractional RMS |
|---|---|---:|---:|---:|---:|
| J0037-0942 | fitted | 1.530 | 1.569390 | +2.575% | 1.240% |
| J1112+0826 | fitted | 1.490 | 1.348843 | -9.474% | 1.693% |
| J1204+0358 | fitted | 1.310 | 1.345144 | +2.683% | 1.368% |
| J1402+6321 | fitted | 1.350 | 1.385681 | +2.643% | 1.951% |
| J1621+3931 | optimizer-excluded; previously exposed | 1.290 | 1.136533 | -11.897% | 3.499% |
| J1630+4520 | optimizer-excluded; previously exposed | 1.780 | 1.580135 | -11.228% | 1.822% |

### Nuisance-matched ordinary-matter control

With the same R9 shared stellar normalization, individual stellar offsets and orbit parameters/priors, a zero-companion model gives fitted stellar error 3.370%, fitted angle error 9.829%, and 31-galaxy comparison RMSE 35.192 km/s. On the two excluded lenses its stellar error is 15.174% and angle error 34.949%. This is stronger than a fixed-population baseline, but not a model-evidence comparison: R9 has ten more field/source parameters. No separate R10-matched two-population-scale baryonic control was fitted.

## Why this is not yet a validity claim

R10 stellar chi-square is 75.666 over 27 fitted radial bins and 52.562 over 13 optimizer-excluded bins. The raw rotation chi-square is 211200.744 over 3152 measurements with quoted errors. The remaining differences exceed measurement errors alone; small fractional RMS is not statistical agreement. Simplified source geometry, stellar populations and distance conventions add unquantified systematic uncertainties. We did not enlarge error bars until the model passed.

The positive finite envelope is not a demonstrated energy source. Under an ordinary weak-gravity stored-energy interpretation its equivalent total mass is A*rt/G. R9 implies about 10-22 times fitted stellar mass for the four lenses, over the large outer envelope. This is not mass measured inside the Einstein radius, and its physical fuel has not been demonstrated. A different energy-to-gravity law would need its own specified derivation. Time reparameterization cannot supply missing energy.

No nested star/planet/galaxy calculation, cluster merger, source-removal evolution, gravitational-wave, full image, source-budget or transient timing success is claimed. The galaxy-level closure is not an automatically valid local rule for every star or arbitrary grouping of matter.

## Exchange time: explicit degeneracy, not a claimed discovery

The production/release model has

    dN/dt = P-Gamma1*N-2*Gamma2*N^2
    Nstar = 2*P/[Gamma1+sqrt(Gamma1^2+8*Gamma2*P)].

Scaling all rates by a common k leaves Nstar exactly unchanged, while tau=1/(Gamma1+4*Gamma2*Nstar) is divided by k. A numerical scan from k=10^-12 to 10^12 changes Nstar by at most 1.11e-16 fractionally. The instantaneous-equilibrium limit is allowed in this rate description but gives the same equilibrium gravity/lensing prediction. This does not establish that time disappears or that spatial propagation is instantaneous. A different clock or timeless-exchange proposal must specify an observable change to the source/field mapping; these static data cannot identify the common internal clock rate.

## Numerical reliability and reproducibility

R9-R11 have frozen-parameter radial/angular/deprojection resolution checks. Maximum stellar-speed changes are below 0.015%, angle changes below 0.0007%. R9 light deflection agrees with independent adaptive quadrature at observed and predicted ring impacts to 1.23e-10 relatively. Three independent stellar-source quadratures agree below 7e-16. A no-refit saved-parameter replay reproduced 120 scalar summary values for R9-R11 exactly and verified all ten input hashes. These are numerical checks, not observational validation.

Code: [run.py](run.py), [audit.py](audit.py), [population_followup.py](population_followup.py), [metric_followup.py](metric_followup.py). Complete arrays, all optimizer starts, raw input snapshots, code snapshots, logs, figures and a saved-parameter replay are in the originating conversation's research package. A compact parameter/result record is committed alongside this report. Full empirical records are not replaced by rounded tables.

From the repository root with the pinned NumPy/SciPy versions, use a fresh directory:

```sh
python research_work/results/joint-response-iteration/run.py --output-dir jr1-replay --max-stage 3 --starts 3 --max-nfev 180 --no-final
python research_work/results/joint-response-iteration/run.py --output-dir jr1-replay --max-stage 6 --starts 3 --max-nfev 180 --resume --no-final
python research_work/results/joint-response-iteration/run.py --output-dir jr1-replay --max-stage 8 --starts 4 --max-nfev 260 --resume --no-final
python research_work/results/joint-response-iteration/run.py --output-dir jr1-replay --max-stage 9 --starts 4 --max-nfev 300 --resume
python research_work/results/joint-response-iteration/audit.py --output-dir jr1-replay
python research_work/results/joint-response-iteration/population_followup.py --output-dir jr1-replay
python research_work/results/joint-response-iteration/metric_followup.py --output-dir jr1-replay
```

## Source attribution

Lelli, McGaugh and Schombert (2016), SPARC, arXiv:1606.09251; Bolton et al. (2008), SLACS V, arXiv:0805.1931; Auger et al. (2009), SLACS IX, arXiv:0911.2471. Detailed light/kinematics use the pinned TDCOSMO2025_public release d7f38db341f68be1df0d9ac1fc528c45113f94cf through the repository's prior audit. Existing slacs-component-refit, slacs-resolved-fit, lens-photometric-audit and lensing-data-readiness document the reused Abel, Jeans/PSF, covariance, population-normalization and conditional-geometry calculations. Those mathematical tools are reused, not claimed as newly invented gravitational physics.
