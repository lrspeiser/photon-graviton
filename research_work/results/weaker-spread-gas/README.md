# JR-4: weaker gravity, wider companion structure, and gas

**21 September 2026, America/Los_Angeles. Executed exploratory calculation, not a validated theory or a demonstrated cause of gas richness.**

## Theory first

Matter and radiation may generate an extended companion state whose gravitational response affects stars and light. The user proposes that in the faint, gas-rich problem galaxies this response is weaker and spreads differently, potentially leaving more gas. JR-4 tests that proposal without changing the observed data, source masses, distances, inclinations, or previously fitted R10 coefficients. It separates an effective reduction of gravitational strength from redistribution of the companion state, and separates galaxy dynamics from gas-cloud self-gravity.

**Result:** individually adjusted wider companion configurations bring 27 of the 38 originally overpredicted difficult rotation curves below 20-percent RMS, compared with 22 using weaker overall gravity. The median required spreading in those 38 is 1.547 times the old radial scale. But a single low-surface-density rule helps 11 of them while pushing 21 previously satisfactory galaxies above 20 percent. That rule worsens median rotation error and joint lens agreement. The gas-collapse calculation gives a conditional route to reduced star formation, not a prediction of gas abundance.

## Provenance and scope

Baseline main: `7386427ad7228b8d57319b1be78bbea3c47f31e8`.
Protocol declared before fitting: `215fc0231056b8c9a935d42eb12bd6db04f80bc3`.
Executed implementation published at `2fd9e186deb6416fd34b056678a3f0c3a01e5833`.

Inputs: the attached original `Photon-Graviton-JR1-Complete-Research.zip`, with R10 parameters, photometric/source reconstructions, 149 SPARC galaxies and 3,152 positive-radius rotation rows, six SLACS lens galaxies and 40 resolved stellar-motion bins. The 89/29/31 disk roles are retained. All objects were examined in earlier work; no fresh blind confirmation is claimed. No clusters were modeled.

The scientific computations were executed locally using Python 3.13.5, NumPy 2.3.5 and SciPy 1.17.0. They were not GitHub Actions runs. There are 64 shared-law optimizer starts across the two declared objectives, plus per-object inverse diagnostics. Full arrays, all optimizer starts, failed alternatives and immutable input hashes are preserved in the downloadable package.

## 1. The exact change to the gravitational response

Let gb be the unchanged baryonic radial acceleration and gchi the R10 companion acceleration. The new companion force is

    gchi_new(r) = s/lambda^2 * gchi_old(r/lambda),
    0.1 <= s <= 1, 1 <= lambda <= 10.

In the existing profile this means

    A_new = s*A/lambda,
    rc_new = lambda*rc,
    rt_new = lambda*rt.

For a spherical companion, the implied enclosed mass is therefore

    Mchi_new(<r) = s*Mchi_old(<r/lambda).

When s=1, the total equivalent companion mass is exactly conserved: the same state is distributed more broadly. This is not merely moving the core radius while inadvertently adding total mass. When s<1, the gravitational normalization is reduced; this calculation does not determine whether that means a different coupling, occupation, or energy conversion. The galactic effective-strength variants are not claims that laboratory Newton's constant changed.

Companion-only variants retain gb unchanged. Overall-response variants use

    g_total_new(r) = s*gb(r) + s/lambda^2*gchi_old(r/lambda).

The luminous matter is not redistributed. The same force enters stellar Jeans projection and lens deflection under R10's inherited equal-potential assumption. No separate lensing gain or new stellar/orbital nuisance fit is introduced.

The spatially expanded static configuration is not itself a derivation of the dynamics or energy needed to expand it. Increasing its radius changes binding energy; a physical formation model must account for that work.

## 2. Independent-galaxy feasibility: does the proposed direction help?

These fits ask what each individual galaxy would require. They are not a shared predictive law. Original difficult galaxies have more than 20-percent fractional rotation RMS; 38 are overpredicted in mean signed fractional velocity, and 16 underpredicted.

| Individual allowed change | Original overpredicted 38 now below 20% | Median error of those 38 | Median requested change | All 149 below 20% |
|---|---:|---:|---|---:|
| None | 0 | 38.42% | None | 95 |
| Weaker total galaxy response | 22 | 19.44% | Strength 0.579 of baseline | 117 |
| Weaker companion only | 22 | 18.43% | Companion strength 0.448 | 117 |
| Wider companion, same total equivalent mass | 27 | 14.38% | Radial scale 1.547 times baseline | 122 |
| Both weaker companion and wider structure | 27 | 13.72% | Median strength 1; median spread 1.430 | 122 |

The combined medians do not describe one representative paired solution. Two pure-spreading fits reach the tenfold spread limit; the combined calculation has three upper-bound cases. All 16 originally underpredicted difficult galaxies prefer no reduction or spreading and none crosses below 20 percent. These changes cannot increase the enclosed gravitational influence at any radius for an increasing spherical mass profile.

Examples, with independent per-object fits:

| Galaxy | Original RMS | Weaker-total result | Pure-spreading result |
|---|---:|---|---|
| UGC07125 | 55.84% | 4.48%, strength 0.412 | 9.86%, spread 2.093 |
| IC2574 | 57.62% | 18.16%, strength 0.413 | 8.39%, spread 1.726 |
| UGC07559 | 51.81% | 9.26%, strength 0.437 | 4.03%, spread 1.749 |
| KK98-251 | 69.71% | 20.02%, strength 0.354 | 7.34%, spread 1.851 |

These show that a less centrally concentrated companion can change the shape of the curve, not only its amplitude. They do not show that the real galaxies possess the fitted new distribution. Catalog gas fraction, stellar populations, and noncircular-motion uncertainties remain relevant.

## 3. One environment-dependent law, rather than separate corrections

Use the fixed-source compactness proxy

    Sigma = Mb/(2*pi*Re^2), in Msun/pc^2,
    w = Sigma0/(Sigma+Sigma0),
    s = exp(-a*w), lambda = exp(b*w).

Mb is the original R10-normalized stellar-plus-gas mass. Re is the fixed photometric half-light radius, including when gas extends farther. Thus Sigma is not a directly measured gas surface density. Gas fraction and residuals are not explicitly used to choose each correction, although gas mass contributes to Mb.

Compare the baseline, a constant suppression, weaker companion only, pure spreading, weaker-and-spread companion, weaker total response, and weaker total response plus spreading. The common coefficients are fitted on the 89 training galaxies only. The primary objective is the mean per-galaxy fractional squared velocity residual. A separate fit using original quoted velocity errors is retained as sensitivity. The 29 validation galaxies select the primary form; the 31 comparison galaxies and all six lenses are excluded from fitting and selection. All are previously exposed data.

The selected combined-companion candidate returns

    a = 1.35e-20 (effectively zero),
    b = ln(10) (the upper bound),
    Sigma0 = 5.43217 Msun/pc^2.

It is numerically indistinguishable from pure spreading: maximum velocity difference is 4.38e-7 km/s. The selector's preference for the combined label is a numerical tie, not evidence for an extra parameter. The upper bound means the low-density limiting extent has not been identified. The bound was not expanded after seeing the outcome.

### Shared primary results

| Shared model | All-galaxy mean velocity RMSE | Median fractional RMS | Galaxies below 20% | Original overpredicted 38 below 20% |
|---|---:|---:|---:|---:|
| R10 unchanged | 17.71 km/s | 15.05% | 95 | 0 |
| Constant overall suppression | 23.87 km/s | 21.36% | 65 | 13 |
| Low-density weaker companion | 19.68 km/s | 18.72% | 77 | 7 |
| Low-density spreading | 19.13 km/s | 17.09% | 85 | 11 |
| Low-density weaker total response | 20.03 km/s | 19.96% | 75 | 6 |
| Low-density weaker total plus spreading | 19.27 km/s | 17.14% | 81 | 10 |

The selected law's mean fractional-squared error improves from 0.10315 to 0.06756 across all galaxies, because it reduces some large outliers. But the median, physical-unit RMSE and number below 20 percent worsen. It brings 11 originally difficult galaxies below 20 percent while pushing 21 previously acceptable galaxies above it. These metrics are not interchangeable and none is a calibrated acceptance probability.

For the 31 comparison galaxies, mean per-galaxy velocity RMSE changes from 15.09 to 16.67 km/s; below-20-percent count changes from 22 to 18. Their mean fractional-squared error changes only slightly, 0.04762 to 0.04597. No global-success claim is made.

Using the quoted-error objective favors much smaller effects, with pure-spreading Sigma0 about 0.464 Msun/pc^2 instead of 5.432. Its 29-galaxy validation RMSE is 17.70 versus 17.68 km/s at baseline. This sensitivity shows that the strength of the inferred environmental change is not robust to error weighting.

### Lenses remain a connected constraint

Across the six unchanged lenses, the selected shared law changes mean stellar fractional RMS from 1.93% to 2.58% and Einstein-angle fractional RMS from 7.94% to 9.13%. Stellar chi-square increases from 128.23 to 193.80 across 40 bins. The three already under-bent systems remain under-bent. A lower gravitational influence cannot provide their missing projected focusing under this fixed geometry.

The stellar population, conditional static distance convention, spherical deprojection, prior orbital parameters and equal-potential assumption are inherited rather than independently remeasured. This is not a new fit of full lens images.

## 4. Is the needed change really correlated with gas richness?

Among the 38 overpredicted outliers the median inferred gas fraction is 68.5%, versus 53.0% in the 16 underpredicted outliers, using the exact R10 normalization here. This descriptive difference is compatible with investigating gas physics but does not establish a cause.

Across all 149 galaxies, independently fitted total strength has Spearman correlation -0.162 with inferred gas fraction. Independently fitted spreading has correlation +0.094. These are weak associations. Gas fraction and the environmental proxy share source-mass inputs; the data were used to infer the corrections, and many optima sit at no-change boundaries. No causal significance or proof of varying G is claimed.

## 5. Conditional gas physics: what could leave more gas?

There are two physically different possibilities: suppress conversion of gas into stars, or improve retention of gas against escape. Weaker gravity does not automatically do both.

### Gas-cloud self-gravity

For a uniform cold cloud of fixed density, hypothesize G_cloud=s*G. Then

    t_ff = sqrt[3*pi/(32*G_cloud*rho)] proportional to s^(-1/2),
    Jeans length proportional to s^(-1/2),
    Jeans mass proportional to s^(-3/2),

with density and sound speed fixed for the Jeans comparison. At half-strength, collapse takes 1.414 times longer and the critical mass is 2.828 times larger. At quarter-strength those factors are 2 and 8. This is a possible route to less gas turning into stars; no gas is created by the change.

A direct radial shell integration for s=1, 0.5, 0.25 and 0.1 agrees with the analytic finite-radius collapse time to at most 1.60e-12 relative error. These are conditional physics controls, not simulations of the observed galaxies' molecular clouds. The mean galactic rotation field does not determine G_cloud.

### Gas retention and rotational stability can point the other way

If a fixed-shape entire potential is scaled by s, escape speed scales as sqrt(s). At half-strength it is 29.3% smaller: heating or feedback can remove gas more easily. That competes with slower consumption.

In a thin gas disk with fixed column and sound speed, Q=cs*kappa/(pi*G_cloud*Sigma_g). If only the large-scale rotational field is weakened, kappa falls as sqrt(s) while G_cloud stays fixed, so Q decreases: rotational support against local collapse is reduced. If both the rotational field and gas perturbation coupling scale by s, Q increases as 1/sqrt(s). These limits have opposite trends. They are not a claim that one Q threshold fully determines star formation.

### A direct geometric route: thicker, more diffuse gas

As a separately labeled post-fit analytic calculation, the spherical companion's vertical harmonic curvature at radius R is nu_chi^2=gchi(R)/R. Isothermal tracer gas in a confining harmonic field has thickness h=cg/nu. If baryonic vertical curvature is unchanged and nonnegative, and companion gravity weakens, then

    1 <= h_new/h_old <= sqrt[gchi_old(R)/gchi_new(R)].

At fixed gas column and local G, the corresponding free-fall time ratio is at most the square root of the thickness ratio. For the independently fitted pure-spreading configurations of the 38 overpredicted outliers, evaluated at their photometric Re, the median upper bounds are a 1.703 thickness ratio and 1.305 free-fall ratio. For the shared spreading law, the corresponding upper bounds are only 1.205 and 1.098.

These are upper bounds under external harmonic confinement, not measured thicknesses, robust star-formation predictions, or full hydrostatic solutions including gas self-gravity. Dominant baryonic confinement reduces the effect. The physical interpretation is nevertheless explicit: a more extended companion can lower vertical compression and leave gas less dense without changing the microscopic gravitational constant.

## 6. What was not measured or established

No star-formation-rate data, molecular-cloud velocity dispersions, cooling, metallicity, feedback, accretion/outflow, gas vertical profiles or source histories were fitted. Therefore JR-4 does not predict the observed gas fractions or determine which way causality runs. Low light could be a consequence of low star formation, but that was not established by fitting circular speeds.

Observational context is compatible with the broader question, not proof of this explanation: Leroy et al. (2008) find lower star-formation efficiency in atomic-gas-dominated environments, but free-fall/orbital-time changes alone do not reproduce its full decline. Bigiel et al. (2010) find very inefficient star formation in outer atomic disks. Neither study demonstrates weaker fundamental gravity. No star-formation measurements from those studies were added to the fitting data here.

A complete next construction needs the companion's response to a gas density perturbation, including spatial and temporal scale dependence, coupled to pressure, gas self-gravity, cooling and energy flow. That decides whether the successful broader background field suppresses cloud formation or merely changes rotation. Turning down a static force is not that derivation.

## Numerical reliability and preservation

Original R10 disk predictions reproduce exactly in the current environment. All original input hashes are unchanged. The companion dilation conserves total equivalent mass to 2.22e-16. Selected lens resolution doubling changes stellar velocities by at most 0.0146% and angles by at most 0.000563%; independent ray quadrature differs by at most 4.16e-10 relatively. These are smaller than observational/model residuals.

Independent readback: 4,336 checks, zero failures. It reconstructs every saved model/galaxy error measure, lens covariance chi-square, mass/dilation identities, input hashes and collapse controls. This validates arithmetic and numerical implementation, not the physical hypothesis.

One inert conditional expression in the published source differs from the original local source. Its Git blob was verified as ddbe3c1e70108036678f50e2911ac383110b04eb and the exact published code was separately rerun. All scientific prediction files, inverse results, gas controls and selection records are byte-identical. Both provenance records are included; elapsed times and source hashes naturally differ.

## Reproduce

Extract the original-JR1.zip included in the package, retaining its Photon-Graviton-JR1 folder. With NumPy 2.3.5 and SciPy 1.17.0:

    python code/run.py --jr1 Photon-Graviton-JR1 --output fresh-results
    python code/audit.py --jr1 Photon-Graviton-JR1 --results fresh-results
    python code/vertical_response.py --jr1 Photon-Graviton-JR1 --results fresh-results

Use a fresh output directory. The original data are not overwritten. The first two commands reproduce the primary campaign; the third is the separate conditional thickness calculation, not another optimization.

## Sources

- Project R10 code, parameters and inputs: original JR-1 package and its manifest, retained in the downloadable archive.
- Original data: Lelli, McGaugh & Schombert (2016), SPARC, https://arxiv.org/abs/1606.09251; pinned SLACS/TDCOSMO inputs documented in JR-1.
- Leroy et al. (2008), https://arxiv.org/abs/0810.2556, star-formation efficiency and gas environments.
- Bigiel et al. (2010), https://arxiv.org/abs/1007.3498, inefficient star formation in outer disks.

**Decision:** retain environmental redistribution as a concrete direction that explains much of the overpredicted subgroup when fitted individually. Do not promote the tested one-variable environmental law, declare lower universal G, or claim that the abundance of gas has been explained.
