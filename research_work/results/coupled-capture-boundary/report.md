# Capture boundary and finite-mass audit

The arbitrary outer capture radius substantially changes the inferred deposited mass. Refitting one amplitude partly hides that change in the inner rotation curve. None of these cases establishes a physical outer boundary or a successful complete galaxy model.

## Model and provenance

The project postulate is `d rho_D/ds = C [W/(W+Wstar)]^p`, with normalized exposure `0 <= s <= 1`, zero initial deposits, `W=max(-Phi_b-Phi_D,0)`, and `Wstar=40000 (km/s)^2`. C has density units, solar masses/kpc^3. The potential is referenced to infinity. Newtonian Poisson gravity is known mathematics; this capture response, uniform optically thin illumination and permanent in-place retention are assumptions. C combines supply, opacity and exposure, not a measured energy source. Deposits feed back into W throughout the integration.

The fixed axisymmetric ordinary-matter model and 542 training stars are unchanged from the [100 kpc calibration](../coupled-capture-calibration/report.md). Twelve observed quantities are approximate Jeans-corrected training proxies, not individually measured circular velocities: tracer scales are fixed and the mixed vertical term is omitted. The objective is unweighted RMS, not a likelihood or confidence test. Ordinary matter alone gives RMS 75.24 km/s under this reduction. No reserved observations or dark-halo-derived vertical-force summaries are used. Mass-model uncertainty remains unpropagated.

The [protocol](protocol.md) specified 30 and 200 kpc as sensitivity choices before their calculation, retaining the exposed 100 kpc result as control. They are not measured galaxy edges.

## Transferring versus refitting the amplitude

Holding the amplitude previously fitted at 100 kpc fixed gives these coarse-grid RMS errors:

| Outer radius (kpc) | p=2 RMS (km/s) | p=3 RMS (km/s) |
|---|---:|---:|
| 30 | 45.88 | 54.60 |
| 200 | 26.41 | 34.84 |

Refitting the single amplitude at each boundary gives the following refined results. This is calibration on exposed data, not independent validation.

| Outer radius (kpc) | p | Fitted C (solar masses/kpc^3) | RMS (km/s) | Deposited mass (solar masses) | Deposit depth at Sun ((km/s)^2) |
|---|---:|---:|---:|---:|---:|
| 30 | 2 | 2.12739e7 | 21.23 | 7.60525e11 | 168280 |
| 100 | 2 | 1.16994e7 | 23.67 | 3.34690e13 | 2176107 |
| 200 | 2 | 9.95276e6 | 23.98 | 2.83975e14 | 9181711 |
| 30 | 3 | 3.93058e7 | 18.72 | 6.23173e11 | 144328 |
| 100 | 3 | 1.93629e7 | 22.62 | 2.97711e13 | 1953815 |
| 200 | 3 | 1.51110e7 | 23.26 | 2.65592e14 | 8609737 |

For p=3 the mass spans approximately 426-fold while RMS spans 18.72 to 23.26 km/s. Even the best of these exposed cases, p=3 at 30 kpc, underpredicts the innermost proxy by 36.51 km/s and overpredicts the outermost by 22.91 km/s. Its predicted curve rises from 205.51 to 249.13 km/s, while those proxies are 242.02 and 226.22 km/s. The [comparison table](comparison.csv) retains all 72 predictions, residuals and exposure labels.

Outer material affects the depth controlling capture, even when its direct interior force is small. In the exactly spherical limit an exterior shell adds interior potential depth without interior acceleration. Our calculation is axisymmetric, so that shell statement explains the distinction rather than asserting exact force cancellation in every component. This creates a route for remote deposits to alter inner capture through feedback. Changing the cutoff changes the physical model, not merely numerical resolution.

## Conditional outer-mass requirement

The following is a conditional derivation using known asymptotic power counting, not a unique project formula or a general impossibility claim about nonexpanding universes.

Assume an isolated finite-mass source, a nonzero uniform incoming companion intensity extending arbitrarily far, in-place retention for positive exposure, and the proposed low-depth opacity `kappa proportional to W^p`. The finite-mass assumption implies `W ~ G M_total/r` far away. With asymptotically unattenuated incoming intensity, deposition then has `rho_D proportional to r^(-p)` at leading order. Consequently:

`M_D(<R) proportional to integral r^(2-p) dr` at large radius.

For p=2 this tail grows linearly with R; for p=3 it grows logarithmically. Thus neither law under those assumptions supplies an isolated finite-mass outer profile on an infinite domain. A low-depth exponent greater than three is necessary for an integrable tail in this particular setting, but is not sufficient to establish stability, inner agreement or a physical mechanism. Feedback can invalidate the assumed finite-mass asymptotic form; that does not rescue a finite-mass solution under the same assumptions.

Other explicit possibilities are declining incident supply, an actual transport front, escape/redistribution, or a physically defined capture threshold. Finite retention age alone does not impose a spatial boundary when illumination already exists everywhere. A multi-galaxy background also changes the isolated boundary problem and must be modeled explicitly. None of these alternatives is implemented or passed here.

## Numerical verification and reproducibility

`python research_work/results/coupled-capture-boundary/run.py 30` and the corresponding `200` command reproduce the new integrations. Set OPENBLAS_NUM_THREADS and OMP_NUM_THREADS to 1 for practical runtime. `export.py` reads those results and the committed 100 kpc results to produce comparison.csv; it performs no fit and opens no additional star samples.

Calibration uses 400 radial nodes, 64 angular nodes, multipoles through 24 and 128 midpoint time steps. The selected amplitudes are checked with 256 steps on the same grid and with 800 radial nodes, 128 angular nodes, multipoles through 48 and 256 steps. All new changes in predicted speeds are below the predeclared 0.5 km/s gate; the largest is 0.20848 km/s. No amplitude optimum lies at a scanned endpoint. Full scans, optimizer outcomes, positivity/density-bound assertions and input hashes are retained in the runner and JSON results. These checks verify the implemented calculation, not its physical assumptions.

## Consequence for the active goals

Goal 5 now has direct evidence that a depth-dependent, optically thin, permanent in-place capture law plus an arbitrary cutoff cannot yet be presented as the galaxy explanation. Goal 4 still lacks a physical supply and retention budget. The next model decision must specify outer behavior and retention before further observational tuning. An explicit integrable low-depth law is a testable candidate, not a preferred result inferred from these fits. No redshift equation is changed or validated by this audit; joint timing, energy, lensing and genuinely withheld evaluation remain incomplete.
