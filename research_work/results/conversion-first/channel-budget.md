# How much redshift can come from conversion without event stretching?

The observationally preferred limit of this restricted two-channel family assigns nearly all the logarithmic redshift to an event-stretching process. Adding a large stationary conversion contribution worsens the timing prediction. This selects a required behavior; it does not supply a radiation-powered clock field or establish companions as its cause. All six research objectives remain open.

## Formula and provenance

**New modeling choice for this comparison, not a claim of unique physics:** split the accumulated logarithmic shift into nonnegative components A_t (a process that stretches phase and event time together) and A_s (stationary energy conversion with unchanged event spacing). Assume achromatic behavior, conserved photon number, and one common proportion across objects:

    A = A_t + A_s = ln(1+z)
    A_t = b A,  A_s = (1-b) A,  0 <= b <= 1.

**Known phase/time-map calculus:** for a phase-preserving signal whose arrival map is t_o=f_D(t_e), frequency scales inversely with f_D'. An affine map over the event gives

    t_o - t_o,ref = exp(A_t) (t_e - t_e,ref)
    S = dt_o/dt_e = exp(A_t)
    nu_after_time / nu_emit = exp(-A_t).

The reference arrival time includes the unknown common flight time. This is an effective event-window map, not a globally derived propagation law. It is not licensed to reset separately for each pulse or to assume a different law for each object. The existing evolving-clock candidates must generate this map with their actual endpoint clocks; fitting b does not resolve their field/source difficulties.

**Known fractional-loss solution and energy/count bookkeeping**, conditional on these assumptions:

    E_observed / E_emitted = exp(-A_t) exp(-A_s) = 1/(1+z)
    S = (1+z)^b
    spectral aging rate = 1/S = (1+z)^(-b)
    F_bol = L_bol / [4 pi D^2 (1+z)^(1+b)].

D is a stipulated geometric distance in the static Euclidean flux calculation. No expansion-derived luminosity distance is inserted. The power law in b is already used in the supernova literature; it is not our invention. The two-channel interpretation is the hypothesis being constrained. Brightness here is a prediction, not a fit to measured luminosities. Intrinsic evolution, dust, selection, peculiar velocities and lens magnification are omitted and cannot be silently assigned to b.

## Calibration and comparison with actual observations

We use the previously read [DES published timing estimate](https://arxiv.org/html/2406.05050v2), b=1.003 with statistical scale 0.005 and systematic estimate 0.010, as a summary calibration. Projecting its central estimate onto the allowed interval gives b=1, stationary fraction zero. This is a constrained fit to a published estimate, not our independent raw-photometry likelihood. Source and reduction assumptions remain those discussed in timing-verdict.md.

Without retuning b on the spectral sample, we compare with all 35 previously exposed [Blondin et al. Table 3 aging rates](https://arxiv.org/html/0804.3595v1#S4.T3). Their ages come from spectral templates, including a small time correction in the nearby reference sample. We keep the original aging-rate errors and use a diagonal residual score; shared template uncertainties are not included. These are different observations, but this retrospective comparison is **not an untouched or statistically independent validation**. No final holdout was opened.

| Stationary share of ln(1+z) | b | Spectral residual chi-square, 35 rows | Event stretch at z=1 | Same-distance flux fraction at z=1 |
|---:|---:|---:|---:|---:|
| 0% | 1.000 | 26.949 | 2.000 | 0.250 |
| 1.2% | 0.988 | 26.889 | 1.983 | 0.252 |
| 2.7% | 0.973 | 26.850 | 1.963 | 0.255 |
| 10% | 0.900 | 27.262 | 1.866 | 0.268 |
| 50% | 0.500 | 50.608 | 1.414 | 0.354 |
| 100% | 0.000 | 150.569 | 1.000 | 0.500 |

The spectral sample alone does not tightly separate 0% from 10%; the tighter constraint comes from the DES summary. We do not sum these results into a global significance.

For an explicitly descriptive sensitivity calculation, let b vary by k times the sum 0.005+0.010 around 1.003. At k=1,2,3 the maximum allowed stationary logarithmic shares are 1.2%, 2.7%, and 4.2%. **These are not confidence limits:** adding the quoted scales does not define a probability distribution or account for unknown source evolution. With a free intrinsic timing exponent b_source, the observed exponent becomes b+b_source and this channel constraint is lost unless source behavior is independently constrained.

## Energy consequence and next physical requirement

The difference between emitted and observed photon energies is E_emit[1-exp(-A)], when expressed in the stipulated endpoint standards. This is not automatically a local energy deposit into a separate field: geometric redshift compares energies measured by different observers. See the subsequent [receiver-source derivation](../inhomogeneous-clock-completion/receiver-source.md). It is incorrect to count this difference again when event stretching lowers the arrival power: the longer duration cancels the arrival-rate factor in integrated received energy.

Only in an additional common-energy bookkeeping model with simultaneous physical loss rates alpha_t=b alpha and alpha_s=(1-b) alpha, acting on the same E and assigned explicit receivers, do the channel energy integrals obey

    Q_s = integral alpha_s E ds = (1-b)(E_emit-E_obs)
    Q_t = integral alpha_t E ds = b(E_emit-E_obs).

The logarithmic share need not equal the energy share for sequential or spatially varying channel proportions. These integrals are not established physical transfer terms for our lapse-metric candidate. The earlier wording identifying Q_t with delivered work was too broad: observer-dependent energy change is insufficient to establish such work or companion production. If only stationary conversion feeds companions, the added common-energy loss model leaves them the Q_s portion. A process that produces both time stretching and companions remains possible, but its local coupling and conservation law must establish that transfer rather than infer it from endpoint redshift.

**Decision:** retain b approximately 1 as the timing requirement for a candidate joint mechanism. Do not add a large stationary loss on top of a clock mechanism already producing the full measured shift. The next physical calculation must identify the receiver of Q_t and produce the arrival map from the same interaction. The observational comparison shows why that missing step matters: a formula that transfers energy but leaves event spacing unchanged cannot be the dominant explanation under these source assumptions.

Run channel-budget.py to reproduce the scores and 210 per-object predictions. The input hash and assumptions are recorded in channel-budget.json. This work supersedes no historical measurement and claims no new blind prediction.
