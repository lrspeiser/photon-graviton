# Depletion limits capture, but does not establish a galaxy model

This calculation replaces the untouched isotropic supply in the previous feedback experiment with attenuated incoming rays. Each captured portion is removed from the energy available farther along that chord. The absorption ledger balances. Attenuation substantially changes the deposited mass, but does not universally remove outer-boundary sensitivity or supply mechanical support for deposits.

## What changed

The spherical Plummer source, depth-dependent capture family and chronological gravity feedback are retained. The original amplitude C is split into `C=B*kappa0`, where B measures incident mass-equivalent fluence and kappa0 sets absorption strength. The [derivation](derivation.md) defines these quantities and integrates both legs of every chord under isotropic illumination. This split is a proposed model specification; the absorption and Poisson mathematics are established, not claimed as new formulas.

Three original cases are tested: p4,C1; p6,C1; and p6,C100. For each, kappa0=.1,10,1000 at domains R30 and R100. Holding C fixed while changing kappa0 changes B inversely. This exposes a degeneracy hidden by the thin-supply amplitude; it does not hold the source intensity constant.

The model evolves deposits from zero for one normalized exposure. No astronomical data are accessed. All quantities are dimensionless in the same G=M_b=a=Wstar=1 convention as the [steep-tail experiment](../integrable-capture-feedback/report.md). These are not inferred Galactic masses, measured conversion rates or physical exposure times.

## Energy-accounting result

For boundary radius R and isotropic illumination, incoming mass-equivalent fluence is `B*pi*R^2`. At every opacity update, the deposited energy plus transmitted energy equals that incoming amount. Independently summing absorption over full chords also matches the mass accumulated in all shells. The solver checks both identities at relative tolerance 1e-7.

The derived finite-domain bound `M_D<=B*pi*R^2` replaces the previous thin-supply volume bound. An area-dependent bound is not an isolated finite-mass proof: it still increases as the domain grows. A highly absorbing outer layer can intercept nearly all incident energy while shielding the interior.

This is explicitly an absorption ledger. Energy in transit, momentum exchange, binding/support energy, finite-speed gravity and a calibrated photon source are absent. It must not be called a complete physical energy-conservation solution.

## Physical assessment

The low-source p6,C1 branch remains comparatively insensitive to increasing the domain from30 to100 at the tested opacities. Increasing opacity while reducing incoming fluence at fixed C reduces its mass. This preserves a conditional weak-source candidate; it does not establish an adequate amount or spatial distribution of extra Galactic gravity.

The p4,C1 and p6,C100 cases continue to show substantial domain dependence. More opacity greatly reduces their mass at fixed C, but changing the numerical boundary still changes the amount and location of capture. Strong attenuation also sharply suppresses deposition near the center. Thus the absorption budget and the stellar-force profile are separate requirements.

The smooth opacity remains positive wherever W is positive. It suppresses weak-well deposition but does not enforce strictly zero permanent collection in voids. To meet that preference literally, the retention dynamics must prevent a persistent bound population there; this model does not yet do so. The opacity is also independent of companion frequency by assumption, not by a demonstrated spectral calculation.

## Numerical status and reproduction

The [protocol](protocol.md) sets the parameter grid and numerical gates before computation. The initial128-shell/4-impact-node and256-shell/8-impact-node comparisons fail the1-percent force/depth/mass gate for10 of18 cases. All initial constant-opacity, zero-opacity and analytic Plummer gravity tests pass. Initial failed comparisons remain in results.json; they must not be reported as successful merely because their mass estimates agree better than their forces.

`run.py` reproduces the initial grid, saves results, and intentionally exits with a failed-gate assertion when any initial check fails. `refine.py 30` and `refine.py 100` continue only affected cases at512 and, if needed,1024 shells, preserving the gate and each intermediate failure. `export.py` combines the latest evidence and labels each case's numerical status in comparison.csv. Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1 for practical runtime.

The numerical comparison uses relative errors even for nearly zero forces. This demanding gate is retained; a tiny absolute force is not silently reclassified as a pass. The final checkpoint below reports any remaining unresolved cases.

An additional uniform thin-source check (`check_thin_source.py`) isolates local angular integration accuracy. Eight impact nodes have a maximum local source error of0.1914 percent, failing a supplementary0.1-percent criterion;16 and32 nodes reduce it to0.0261 and0.00341 percent. These checks preserve the eight-node failures in thin-source-check.json. This establishes a known local angular bias, not a rigorous error bound on nonlinear evolution.

To check its propagation in the weak-source candidate, `check_angular_evolution.py` reruns all three p6,C1 opacities at R100 on the same256-shell grid with16 impact nodes. Changes in mass, depth and selected squared speeds are all below0.000836 percent. The separate results are retained in angular-evolution-check.json, without replacing the primary radial-refinement records. This supports the weak-source result; high-source cases have not received this additional angular-evolution test.

## Consequence for the active goals

We now have an explicit depleted-supply version of capture with an auditable incoming/outgoing ledger. This advances goals4–5 without closing them. The next physical requirement is a retention model: specify whether the deposited state behaves as bound particles or field stress, and what prevents it from simply moving after capture. That choice must enter the force and energy equations, rather than leaving deposits permanently fixed by assumption. The incoming supply must subsequently be connected to the same photon/time mechanism used for redshift. No observational holdout or complete joint model has been passed.

There is already a useful conditional constraint on retention. If companions carry photon-like momentum E/c and absorption transfers it locally, the spherical inward supply produces an inward momentum force, not outward support. The derivation gives `F_r=-(1/r^2) integral_0^r q(u)u^2 du` and `f_rad=kappa F_r/c<=0`. A static deposited fluid would need additional pressure with `dP/dr=-rho_D g+f_rad`, or a specified particle/field stress alternative. These are conditional force-balance requirements, not a completed support model or a numerically verified momentum ledger.

## Final numerical checkpoint

| p | C | kappa0 | Mass, R30 | Mass, R100 | Mass ratio | All diagnostic gates, R30 / R100 |
|---|---:|---:|---:|---:|---:|---|
| 4 | 1 | 0.1 | 20987.85 | 252791 | 12.04463 | pass / pass |
| 4 | 1 | 10 | 212.6452 | 2416.095 | 11.36210 | pass / pass |
| 4 | 1 | 1000 | 1.361981 | 10.02173 | 7.35820 | unresolved / unresolved |
| 6 | 1 | 0.1 | 0.2315506 | 0.231797 | 1.00106 | pass / pass |
| 6 | 1 | 10 | 0.2130322 | 0.2132694 | 1.00111 | pass / pass |
| 6 | 1 | 1000 | 0.04949779 | 0.04964893 | 1.00305 | pass / pass |
| 6 | 100 | 0.1 | 2605486 | 3.047513e+07 | 11.69652 | pass / pass |
| 6 | 100 | 10 | 27457.26 | 304988.5 | 11.10775 | pass / pass |
| 6 | 100 | 1000 | 247.7481 | 2704.334 | 10.91566 | pass / pass |

16 of18 cases pass the unchanged numerical gate after refinement. Maximum final incoming/transmitted/deposited ledger relative error is 2.86e-12; maximum independent chord-versus-shell absorption discrepancy is 3.32e-12. These numerical identities do not include missing physical energy sectors.

Unresolved: R=30.0, p=4, C=1.0, kappa0=1000.0, shells=1024; maximum relative diagnostic change=0.038349, squared force contribution at r1=4.1577e-17. Do not treat the full case as numerically validated.
Unresolved: R=100.0, p=4, C=1.0, kappa0=1000.0, shells=1024; maximum relative diagnostic change=0.048918, squared force contribution at r1=3.6646e-17. Do not treat the full case as numerically validated.

Full values and status labels are in [comparison.csv](comparison.csv); initial and subsequent failed gates remain available in the JSON files. The mass ratios diagnose boundary dependence only over the tested domains; they are not an infinite-domain extrapolation.
