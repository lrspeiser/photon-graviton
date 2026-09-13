# Observational verdict on stationary conversion alone

The preceding user-facing turn clarified that numerical passes are not evidence of a new physical theory. This turn executes a direct, retrospective comparison of an existing branch with a published observational summary. It is not a prospective holdout test or a new raw-data analysis. All six objectives remain open.

## Candidate being tested

Call this existing branch S0: static Euclidean paths, constant propagation speed, unchanged local clocks, stationary photon-energy loss into companions, conserved photon number, and no intrinsic source-duration evolution. Its redshift follows the known fractional-loss solution 1+z=exp(integral alpha ds). Its timing prediction is t_arrival=t_emission+T(D), hence duration stretch S=1, regardless of alpha. The companion interpretation is hypothetical; the loss solution and time-difference identity are known mathematics.

This is a specific branch, not a definition of all our possible nonexpanding models. The prior nearby-distance fit in report.md does not establish its timing behavior observationally or determine alpha at every distance.

## External observational comparison

[DES timing study, White et al. (2024)](https://arxiv.org/html/2406.05050v2) reports b=1.003 +/-0.005 statistical +/-0.010 systematic in duration proportional to (1+z)^b. Its analysis matches emitted-wavelength ranges across observed filters; Appendix C checks a reference construction without time de-redshifting. Appendix A discusses source evolution, partly using quantities already fitted under a time-dilation convention. Those assumptions remain relevant; this is not an assumption-free measurement of individual intrinsic clocks.

| Redshift | S0 duration ratio | Ratio from published fitted relation |
|---:|---:|---:|
| 0.1 | 1 | 1.1003 |
| 0.2 | 1 | 1.2007 |
| 0.5 | 1 | 1.5018 |
| 1 | 1 | 2.0042 |
| 1.2 | 1 | 2.2052 |

Rows evaluate one fitted relation, not five independent observed supernovae. The discrepancy is large compared with the quoted uncertainty scales, but we do not turn the systematic estimate into an independent Gaussian or claim a calibrated extreme significance. The conclusion is conditional on the study's source-comparison and reduction assumptions: **S0 does not reproduce the published timing result and should not be promoted as the complete explanation.** Matching its nearby redshift trend cannot repair this failure.

## What a revision must supply

For t_arrival=t_emission+T(t_emission,D), an actual propagation mechanism must produce

    S = 1 + partial T / partial t_emission.

A large common travel delay is insufficient; later portions must acquire a different delay. Setting S=(1+z)^b without deriving T merely imposes the observed relation. The identity does not require an expanding universe, but a nonexpanding field or interaction model must generate it and meet other timing constraints.

Alternatively, if S0 is retained, all the apparent stretching must come from intrinsic source behavior or observational selection/reduction. In the simplified multiplicative population model, b_apparent=b_propagation+b_source. With b_propagation=0, source effects must supply approximately 1.003, implying a factor about two in characteristic intrinsic duration between z=0 and 1. That is a required explanation, not evidence that such evolution exists; arbitrary source evolution would make timing alone unidentifiable.

Any propagation revision also changes brightness bookkeeping. With the same stipulated distance D, emitted luminosity L, photon-number conservation, an achromatic energy shift, and duration stretch S, the known flux ledger is

    F=L/[4 pi D^2 (1+z) S].

At z=1, S0 gives half the unshifted same-distance flux; adding the measured timing factor gives about 0.2495. These are conditional predictions, not a brightness-data fit or an adoption of cosmological luminosity distance. The added factor must be included in a joint model rather than adjusted independently.

## Recorded result and next work

timing-verdict.py writes timing-verdict.json with the source, fixed branch assumptions, predicted and reported exponents, and the table. It checks the timing identities, not the observational reduction. This comparison supersedes treating S0 timing as merely untested against the published summary, while preserving the need for our independent raw-photometry analysis.

The next substantive theory step is a specified evolving propagation/receiver model that calculates energy transfer and arrival-time mapping together. If no such candidate survives, report that failure rather than count more numerical passes as support. Cluster-capture and orbit calculations remain conditional components and do not rescue a failed propagation branch. No final observational holdouts were opened.
