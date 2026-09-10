# Test of the screened clock recommendation

9 September 2026. Conditional calculation, not an observational validation. [Protocol](protocol.md), [executable calculation](run.py), [all numerical tables](tables.md), [machine-readable results](results.json).

## Main findings in ordinary language

The user's short-path argument is correct for accumulated photon conversion: the unchanged exploratory redshift coefficient implies very small losses across Solar System distances, even without suppressing that conversion locally. It does not establish that every relevant measurement is too insensitive to see them, and it does not bound a separate clock-rate change.

The recommended slower-void multiplier has a real sign issue. Ordinary stationary clocks run faster when moving to higher potential (outward in a simple well). The proposed slower-void factor decreases outward, so it opposes this ordinary trend. Small amplitudes can leave the total trend unchanged; larger amplitudes can reverse it. No clock has to stop or run backward. A reversal of the spatial gradient is not reversal of time.

Screening can reduce local corrections by many orders of magnitude, but it is not by itself a complete way to preserve gravity. Under the previously explored universal lapse coupling, the change in the factor across the transition produces an additional force. Slower-void factors push outward; faster-void factors pull inward. A sharp switch can protect Earth while generating a large transition effect. The complete field theory might couple differently, but that must be derived rather than assumed.

## Unscreened photon-loss calculation

**Provenance: established exponential fractional-loss mathematics applied to our proposed transfer law.**

\[
A=\alpha L,\qquad z_{\rm transfer}=e^A-1,\qquad
f_{\rm energy\ lost}=1-e^{-A}.
\]

Use alpha=0.0002488993286382367 per Mpc from the [earlier exploratory conversion fit](../conversion-first/report.md), without refitting. This is an exposed-data calibration over roughly 10-93 Mpc, not a derived void rate or a tested local prediction. Applying it unchanged at all lengths below is an extrapolation used to test the user's scale argument. The earlier fit did not prove a companion mechanism or adequate predictive uncertainties.

| Unscreened path length | Predicted fractional photon energy loss |
|---|---:|
| 20,200 km, an illustrative ground-to-satellite-height path | 1.6294e-19 |
| 1 au, approximately the Earth-Sun distance | 1.2067e-15 |
| 100 au | 1.2067e-13 |
| 1 light-year | 7.6313e-11 |
| 1 million light-years | 7.6310e-5, or about 0.00763% |
| 100 million light-years | 0.0076022, or about 0.760% |

There is no million-light-year threshold in this law: the effect exists from the start and accumulates continuously. The implied fractional-loss rate is about 2.4182e-18 per second of local travel in a uniform environment. Whether that rate is measurable depends on the actual frequency-comparison, Doppler correction, storage, source and instrument protocol; it is not determined by path length alone. In particular, an ideal year-long photon storage test is not the same as a one-metre one-way path or an ordinary cavity decay measurement.

## What was tested for clocks

**Provenance: proposed screened multiplier, using established clock-ratio and derivative mathematics.**

\[
S(W)=\frac{1}{1+(W/W_*)^n},\qquad
F(W)=1+\sigma\epsilon S(W),\qquad
q_{\rm total}=q_{\rm base}F.
\]

Sigma=-1 is the slower-void recommendation; sigma=+1 is an alternative faster-void control. The control is not a newly adopted model. The baseline q_base=sqrt(1-2W/c^2) is a static weak-field surrogate. For the summed Earth, Sun and external background potentials it is NOT an exact multi-body GR solution or a satellite timing model. Comparisons are between hypothetical stationary clocks; actual orbit velocities, tides, rotating reference frames and link corrections are not modeled.

The 180 predeclared cases use five epsilon values, three W_star values, three transition exponents, two signs, and two external-background assumptions. The background values 0 and 3e10 m^2/s^2 are sensitivity choices, not measured Milky Way potentials. The nominal Earth/Sun constants are from [IAU nominal conversions](https://arxiv.org/abs/1510.07674). All cases are reported, with no parameter selection by astronomical outcomes.

**Provenance: exact algebra for the stipulated multiplicative factors.** For clocks A and B, the additional ratio relative to the baseline is

\[
\delta_{AB}=\frac{(q_B/q_A)_{\rm total}}{(q_B/q_A)_{\rm base}}-1
=\frac{F(W_B)}{F(W_A)}-1.
\]

This observable depends on differences between local factors, not solely on F-1 at one location. A uniform factor cancels exactly from local clock ratios and gives no spatial-gradient force in this idealized coupling. Calling a uniform rescaling a large locally measured time anomaly would be incorrect. Photon conversion could still exist under its separate postulate.

Examples below all use W_star=1e8 m^2/s^2 and illustrative background B=3e10 m^2/s^2. The baseline ground-to-GPS-height ratio differs from unity by 5.28485e-10, corresponding to about 45.66 microseconds per day for stationary clocks. This is a gravity-only surrogate; the observed GPS net shift also includes orbital motion. [NIST GPS explanation](https://www.nist.gov/atomic-clocks/a-powerful-tool-for-science/putting-einstein-test).

| Slower-void amplitude epsilon | Transition n | Extra ground-to-GPS-height clock ratio | Extra force / baseline force at transition, conditional coupling |
|---|---:|---:|---:|
| 0.1 | 1 | -4.9359e-7 | -2.3651e7 |
| 0.1 | 4 | -6.7162e-14 | -9.4606e7 |
| 1e-9 | 1 | -4.9343e-15 | -0.22469 |
| 1e-9 | 4 | -6.7162e-22 | -0.89876 |

These are examples from the fixed grid, not best fits. The first produces a roughly 0.043-second-per-day additional clock-ratio effect, overwhelming the baseline in this toy setup. Sharper screening makes the local correction much smaller, but the last column shows why this does not establish harmless gravity throughout the transition. Even a one-part-in-a-billion void amplitude can be a substantial fraction of the very weak baseline gradient there.

The [Galileo redshift experiment](https://arxiv.org/abs/1812.03711) reports a fractional deviation from its modeled relativistic redshift of (+0.19 +/- 2.48)e-5 at one sigma. That is a constraint on the experiment's changing redshift signal, not a universal limit on F-1, and we have not evaluated this model against the Galileo likelihood. No row above is labeled empirically allowed or excluded by that measurement.

## Direction, force and the meaning of zero

For a simple outward-decreasing well depth W:

| Quantity | Behavior from deeper well toward void |
|---|---|
| Ordinary stationary-clock factor | Increases toward its distant reference value |
| Our proposed slower-void factor | Decreases toward 1-epsilon |
| Faster-void control | Increases toward 1+epsilon |
| Combined slower-void factor | Can increase or decrease, depending on amplitude and transition slope |

**Provenance: established Hamiltonian/lapse result within the previously tested universal matter coupling, not a universal theorem about modified gravity.**

\[
\mathbf a_{\rm extra}=-c^2\nabla\ln F.
\]

The force is relative to the static observers in that coupling, not proper acceleration experienced by a free-falling observer. With the same spatial-gradient factor for the baseline and extra term, their ratio is

\[
R(W)=\frac{d\ln F/dW}{d\ln q_{\rm base}/dW}.
\]

For the slower branch at W=W_star, cancellation of the total gradient occurs at

\[
\epsilon_{\rm critical}=
\frac{4W_*}{n(c^2-2W_*)+2W_*}.
\]

**Provenance: algebraic consequence of the stipulated functions, not an empirical bound.** At W_star=1e8 m^2/s^2, epsilon_critical is about 4.4506e-9 for n=1 and 1.1127e-9 for n=4. Above it, the total gradient reverses at the transition in this surrogate. Sixty of the 90 slower-void cases do so; this grid count is not a statistical probability or a measure of viable parameter space. See the [earlier universal-clock report](../universal-clock-coupling/report.md) for the coupling and its limitations.

Zero well depth is a chosen potential reference, not zero elapsed time. Zero additional time dilation means a relative factor of 1. With 0<epsilon<1, the slower multiplier stays positive, and all 10,980 checked lapse combinations are positive. Time never needs to stop, become negative, or reverse in these examples. A clock rate that decreases as a function of position remains a forward-running clock rate.

## Does redshift require slower void clocks?

Not under the adopted energy-transfer postulate alone. The redshift sign is set by positive alpha, which reduces photon energy. Holding alpha_void fixed while changing epsilon changes kappa=alpha_void/epsilon: the same photon redshift can coexist with different proposed clock amplitudes. The unmeasured efficiency is why the observed redshift does not determine epsilon.

For a faster-void alternative, the conversion prescription could be changed to

\[
\alpha=\kappa|F-1|=\alpha_{\rm void}S.
\]

**Provenance: a different proposed interaction law used only for the sign-control comparison.** It has exactly the same prescribed loss profile, but is not derived from the faster clocks. Keeping alpha=kappa(1-F) unchanged would instead produce a negative alpha for F>1 and would contradict the intended energy flow. The faster version aligns the clock-factor direction with the ordinary outward trend and makes the conditional extra force inward; neither feature establishes a correct force magnitude, stable deposits or a successful theory.

If redshift is to come from time geometry alone, the earlier matched-endpoint candidate depends on partial_t q, not merely on whether a region has q below or above 1. Static differences do not automatically accumulate permanent matched-endpoint stretching. Similarly, the stationary transfer law does not automatically stretch the spacing of separate packets or a whole supernova event. These requirements remain open; changing signs cannot bypass them.

## Verification and next decision

The script passed 12 symbolic derivative checks, 360 independent 60-digit clock-ratio comparisons, 720 derivative comparisons, 10,980 positivity checks, zero-amplitude recovery, and uniform-factor cancellation. Maximum relative discrepancies were 2.08e-16 for path conversion, 2.78e-11 for extra clock ratios, and 7.00e-16 for derivatives. Checks verify the mathematics of the declared model, not its physical truth. No catalog was refitted and no energy-supply calculation was performed.

The recommendation is therefore revised in strength: screening is a possible mechanism to investigate, not a demonstrated necessary or sufficient safeguard. The user's distance argument means we should calculate actual local photon-link predictions before asserting that conversion must be screened. The clock amplitude, its sign and its force coupling need independent justification. Do not lock in a substantial slowing of void clocks solely because photons redshift. Retain both sign possibilities as hypotheses until a common interaction supplies a reason to choose one.
