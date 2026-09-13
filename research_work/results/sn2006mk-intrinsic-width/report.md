# Intrinsic source width versus propagation stretch: SN 2006mk

The exposed two-band light curve prefers an effective template width of 1.528877 with starting template phase -3.106 days in this restricted fit. Both bands use the same time map and source width; their previously chosen measured brightness anchors remain fixed. No absolute luminosity or dust model is fitted.

This quantifies a known identifiability problem on the real-data pilot. It does not discover a new degeneracy or validate an intrinsic-source explanation.

## Exact equivalent descriptions

| Propagation stretch A | Intrinsic template width w | Product A*w | Conditional Gaussian score |
| ---: | ---: | ---: | ---: |
| 1.000000 | 1.528877 | 1.528877 | 40.450391 |
| 1.475400 | 1.036246 | 1.528877 | 40.450391 |
| 1.595706 | 0.958119 | 1.528877 | 40.450391 |

Unchanged propagation requires an intrinsic duration 52.9% longer than the reference. Propagation stretch 1+z=1.4754 requires only 3.6% longer. These percentages are relative to the adopted empirical template, not measured bounds on physically allowed supernova diversity. No source-population prior is supplied by this fit.

## Formula provenance and proof

The assumed template coordinate is u=q+(t_observer-t_first_spectrum)/(A*w). Introducing w is an empirical source-time rescaling, not a physical supernova model or a companion equation. This is a known scaling degeneracy: under A -> k*A and w -> w/k, every template coordinate is unchanged. Therefore the integrated spectra used for photometry, temporal ratios, predicted fluxes and shared-anchor covariance all remain unchanged. The light-curve likelihood depends on B=A*w, not A and w separately.

The fitted q is a template coordinate at the first observation; changing w also changes its conversion to physical source days. Treating q as an independently fixed physical clock while changing w would mix definitions. The earlier spectral phase estimates are likewise template-based; this photometry calculation does not establish how a truly slower physical source would alter spectral line evolution.

Changing the photon redshift law is not required for this equivalence: the measured wavelength redshift remains fixed in every row. Absolute brightness is canceled by the two measured anchors. The result therefore concerns relative light-curve shape only, not the full requested redshift/timing/brightness solution.

## Numerical checks

Two differential-evolution searches with distinct seeds over q in [-10,5] and B in [1,2] converge to interior solutions, with score agreement better than 1e-5. The best chi-square is 35.6773 and the determinant-inclusive score is 40.4504. All 25 unanchored points retain their prior eligibility. Direct evaluation of the three decompositions agrees in score to 1e-10 and in all sampled phases to 1e-10 days. Repeated execution reproduces results.json exactly. These checks verify equivalence and optimization under the stipulated source, not source plausibility.

## Consequence for the research claim

The preceding fixed-width unchanged-time model was a poor description; allowing an unconstrained intrinsic width removes that photometric distinction exactly. We must not report the earlier score gap as a general rejection of unchanged propagation. Nor should the new width be called a successful physical alternative merely because it fits.

To distinguish the explanations, specify and independently constrain source evolution: use measured source populations, spectral diversity and physically motivated relations among luminosity, ejecta velocity, opacity and diffusion/decay times. Such constraints must be evaluated with their timing assumptions exposed and must predict spectral evolution as well as brightness. A model-dependent upper limit on plausible width cannot simply be invented.

An arbitrary time-rescaled reference light curve does not provide energy-conserving explosion dynamics or a nuclear clock. A companion-based model still needs its own causal arrival-time map, local energy/momentum accounting and joint source uncertainties. These remain independent of the galaxy-deposition requirements. All six research objectives remain open.

[Starting-phase fit](../sn2006mk-phase-profile/report.md) | [Source-template assumptions](../hsiao-source-audit/report.md) | [Original two-epoch cross-prediction](../sn2006mk-source-pilot/report.md)
