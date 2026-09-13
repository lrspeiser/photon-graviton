# External deposition from a candidate acceleration-squared capture law

## Result

A concrete optional gravity-to-capture rule produces finite collecting area and extended deposition using external companions alone. In a fixed Plummer ordinary-matter well, at dimensionless capture strength A=1, half the captured power is deposited outside 3.287 scale radii. The ordinary-matter half-mass radius is 1.305 scale radii. Increasing capture strengthens absorption along incoming paths and moves the deposition outward.

| A | Effective collecting radius / a | Half-deposition radius / a | Deposits outside a |
|---|---:|---:|---:|
| 0.1 | 0.484 | 3.098 | 92.545% |
| 1 | 1.483 | 3.287 | 93.433% |
| 10 | 3.903 | 4.783 | 98.327% |
| 100 | 8.745 | 9.699 | 99.99998% |

The effective collecting radius is an area-equivalent radius, not the boundary containing deposits; the two radius columns need not agree. This unbounded smooth profile has a finite but extended outer contribution. A large outer deposition fraction does not by itself identify a lensing peak or demonstrate a successful rotation curve.

## Candidate rule, not a locked-in mechanism

Use the known Plummer acceleration for ordinary mass M and scale a:

    g(r)=GM r/(r^2+a^2)^(3/2).

Postulate the capture probability per length

    kappa(r)=chi g(r)^2
            =g(r)^2/(ell_star g_star^2).

chi is one universal dimensional coefficient in a future cross-system test; ell_star and g_star are an equivalent parametrization, not two independently identifiable parameters. This is a phenomenological candidate, not a derived particle interaction, a proven novel law, or a generally covariant theory. g is defined here in the static receiver frame of the Newtonian comparison. The physical frame dependence and relation to a relativistic field require resolution before claiming fundamental physics.

The square was chosen as a testable candidate with a convergent outer collecting area, not inferred from astronomical data. In this cored symmetric well g vanishes at the center and so does local capture. It measures acceleration strength, not potential depth: a deep central potential does not guarantee strong local capture under this rule. That distinction must be tested rather than hidden in terminology.

In dimensionless radius x=r/a,

    a kappa=A x^2/(1+x^2)^3
    A=chi G^2 M^2/a^3.

No arbitrary outer cutoff, local stellar illumination or dark halo is added. The exterior rate falls as r^-4 and has finite effective area.

## Exact ray and volume calculation

For impact b measured in units of a, established straight-ray absorption yields

    tau_full(b)=A pi (1+4b^2)/[8(1+b^2)^(5/2)]
    sigma_eff/a^2=2 pi integral_0^infinity b[1-exp(-tau_full(b))]db.

For deposition at radius r, integrate the opacity upstream along each incoming direction, then average survival over solid angle:

    q(r)=c u_infinity kappa(r) * (1/2) integral_-1^1 exp[-tau_upstream(r,mu)] dmu.

The upstream optical depth is evaluated analytically by writing the opacity as A[(1+r^2)^-2-(1+r^2)^-3]. Three direct semi-infinite line integrals independently check that expression at each A. The spatial integral of q agrees with c u_infinity sigma_eff: the largest relative difference across four strengths is 1.36e-8. Doubling angular nodes from 128 to 256 changes total power by at most 7.95e-8 relative. These are numerical checks, not observational errors or proof of a supported deposit distribution.

The ordinary Plummer enclosed-mass fraction is r^3/(r^2+a^2)^(3/2), giving r_half/a=[2^(2/3)-1]^(-1/2)=1.305. The tabulated deposit half-radii come from integrating the calculated power density, not prescribing an outer shell. The comparison concerns deposition per time; equality with stored-energy shape requires a stationary common history and no redistribution.

## Testable scaling and limitations

At weak capture, sigma_eff approximately equals the volume integral of opacity:

    sigma_eff~(3 pi^2/4) a^2 A, proportional to M^2/a.

At very strong capture the exterior r^-4 tail gives leading scaling

    sigma_eff~pi Gamma(1/3) a^2 (pi A/2)^(2/3), proportional to M^(4/3).

These scalings assume the same chi and fixed external companion bath. They are conditional consequences of the candidate, not empirical galaxy relations or universal bounds. Global receiver competition changes that bath and must be solved with the same cross-sections.

The well is fixed by ordinary matter in this calculation. If deposited energy gravitates, it changes the well and potentially the capture rate; that feedback, the reservoir's support and motion, and incident trajectory bending have not been included. We cannot simultaneously neglect backreaction for this test and claim a self-consistent massive halo. The next step is to test the profile's motion/lensing shape and backreaction, followed by a common-coefficient test using measured ordinary-matter structures. All six objectives remain open; no observational holdouts were opened.

Formula status: the Plummer well, absorption integrals and mathematical solutions are established tools. The acceleration-squared capture prescription is an exploratory postulate in this project; no priority or uniqueness claim is made.

Run `python research_work/results/cluster-acceleration-capture/run.py`.
