# Steeper capture tails: a conditional repair, with a feedback limit

A steeper low-depth capture law can remove the outer-boundary problem in weak-source examples. It does not automatically do so when deposited gravity feeds back into capture. In strong-source examples, inner force convergence conceals a domain-filling density and an unconverged mass and potential. No astronomical data were fitted or validated in this experiment.

## Equations and physical meaning

This is a spherical mechanism diagnostic, not a replacement for the axisymmetric Milky Way calculation. Use units G=M_b=a=1, where a is the ordinary-matter scale radius. Density units are M_b/a^3; potential and squared speed units are G M_b/a. The ordinary potential depth is the known Plummer expression:

`W_b(r)=1/sqrt(1+r^2)`.

Retain the project's proposed capture family:

`d rho_D/ds = C [W/(1+W)]^p`, `W=W_b+W_D`, `rho_D(r,0)=0`, `0<=s<=1`.

Wstar=1 in these units; p=4 or6 replaces the nonintegrable low-depth p2/p3 tail. These exponents are proposed alternatives, not derived microscopic laws or claimed unique formulas. C is an exposure-normalized density supply, not a measured photon-to-gravity conversion rate. Uniform thin incident supply, permanent in-place deposition and instantaneous Newtonian potential updating remain assumptions. In particular, no finite-speed gravity-field evolution or depletion of the incident supply is implemented.

Known spherical Poisson mathematics gives:

`M_D(r)=4 pi integral_0^r rho_D(u) u^2 du`,

`W_D(r)=M_D(r)/r + 4 pi integral_r^R rho_D(u) u du`,

`v_D^2(r)=M_D(r)/r`.

The origin uses the finite limiting expression for W_D. R is a computational truncation, not a physical capture cutoff. Four values R=30,100,300,1000 test whether results settle. The ordinary Plummer potential remains untruncated. In an isolated finite-mass asymptotic solution W falls as 1/r, giving rho_D proportional to r^(-p) and an integrable mass tail for p>3. This conditional result does not establish that feedback admits such a solution at every C.

## Results

The [protocol](protocol.md) was written before computation. Each of 24 parameter/domain combinations was run at two numerical resolutions, without optimization. The following compares R=300 with R=1000. Mass is quoted in ordinary-matter mass units, not solar masses or an inferred Galactic mass.

| p | C | Deposited mass at R=1000 | Mass ratio, 1000/300 | Central depth ratio, 1000/300 | Finite-range convergence gate |
|---|---:|---:|---:|---:|---|
| 4 | 0.01 | 0.0348114 | 1.009044 | 1.000069 | Pass |
| 4 | 1 | 3.40933e9 | 37.07849 | 11.12344 | Fail |
| 4 | 100 | 4.18107e11 | 37.03746 | 11.11122 | Fail |
| 6 | 0.01 | 0.00160055 | 1.000003 | 1.000002 | Pass |
| 6 | 1 | 0.232031 | 1.000004 | 1.000003 | Pass |
| 6 | 100 | 4.08377e11 | 37.03722 | 11.11114 | Fail |

The finite-range gate requires less than 1 percent change in mass, central depth and deposit-only squared circular speeds at r=1,3,10. Passing is evidence of finite-range convergence, not proof of an infinite-domain solution or long-term stability. The weak p4 case passes narrowly on mass and has a slower asymptotic tail than p6.

At p6,C=1, deposit-only squared speeds at r=1,3,10 are 0.04239,0.05639,0.02272. This gives a concrete spatial force prediction; it is not a demonstrated galaxy fit. Its exposure interval and supply do not have a physical age or radiation-budget calibration.

At p4,C=1, the density at the outer edge reaches 0.814 C; at p4,C=100 it reaches 0.998 C, and at p6,C=100 it reaches 0.975 C. Feedback makes the well sufficiently deep to nearly saturate capture far from the original ordinary source. These cases are incompatible with the intended isolated finite-mass interpretation over the tested domain range. They do not prove all time-dependent or supply-limited completions fail.

## Why apparently converged inner speeds can mislead

Known uniform-sphere mathematics gives `M_D(R)=4 pi rho R^3/3`, `W_D(0)=2 pi rho R^2`, while at fixed interior r, `v_D^2(r)=4 pi rho r^2/3`. Thus a nearly uniform deposited density can make the interior force settle while mass grows as R^3 and central depth as R^2.

The expected ratios from expanding R by 1000/300 are 37.037 for mass and 11.111 for central depth. These closely match the saturated cases above. Their squared circular speeds change by at most about 0.11 percent between the two largest domains despite the enormous mass change. Therefore checking rotation curves alone would accept a numerically misleading isolated model. Total potential and mass are required checks because potential also controls the capture rule.

Density remains bounded by C s on each finite domain since the postulated rate lies between zero and C. The failed boundary behavior is not a demonstrated finite-time density singularity. It is rapid feedback into a domain-filling profile, with no physically acceptable isolated mass established. Permanent positive capture also does not produce a stationary total mass at indefinite exposure.

## Verification

The coarse calculation uses 1000 radial nodes with adaptive ODE relative/absolute tolerances 1e-7/1e-10. The refined calculation uses 2000 nodes and 2e-9/1e-12. Every mass, depth and selected force diagnostic agrees within the predeclared 0.5 percent numerical gate; the largest relative change is 0.03534 percent. All retained solver states satisfy the nonnegative density and rho_D<=C s bounds within numerical tolerance.

An independent analytic gravity check uses a truncated Plummer density, for which `M(r)=r^3/(1+r^2)^(3/2)` and `W(r)=1/sqrt(1+r^2)-(1+R^2)^(-3/2)` inside R. The maximum tested relative error is 0.002466 percent. This checks the potential including exterior shells, rather than only the radial force.

Run `python research_work/results/integrable-capture-feedback/run.py` to reproduce all 48 evolutions and analytic checks. [results.json](results.json) retains every coarse/refined case, convergence decision and boundary ratio. These are numerical checks of the specified equations, not observational validation.

## Next physical step

Retain the steep-tail family as a conditional candidate, with both successful finite-range and failed feedback regimes recorded. Do not adopt p6 merely because one synthetic case passes. Before a new galaxy fit, the supply and retention model must specify a physically allowed exposure range and how incoming energy is depleted, transported or redistributed. A stationary deposited distribution also needs support or binding dynamics; gravity alone does not keep newly deposited material fixed at its arrival position.

This experiment advances outer capture and feedback checks in goals4–5. It leaves the unified clock/redshift mechanism, physical energy and momentum closure, galaxy/cluster lensing comparison, environmental redshift test and genuine withheld evaluation open. No nine-stage goal is declared complete.
