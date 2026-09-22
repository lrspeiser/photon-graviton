# JR-9: same-object cross-prediction, and one interaction instead of two rate laws

**22 September 2026. Executed on the six real SLACS systems with the R10 universal
coefficients frozen, plus one derived response kernel. Development on historically
examined objects; not a blind test and not a validated theory.**

Protocol declared before scoring: [PROTOCOL.md](PROTOCOL.md), committed at `6bc559b`
with no cross-prediction score computed. Every forward calculation reproduces the
published JR-1 R10 lens predictions to zero relative difference before anything
new is scored; that check is the first thing each run writes out.

## The one-paragraph summary

R10's published joint fit adjusted each lens's stellar mass and orbital anisotropy
against its stellar motions *and* its Einstein angle together. When you stop doing
that and make each observable predict the other, the agreement is materially worse
than the published numbers suggest, and the residual is not random: it separates the
six lenses almost perfectly by redshift. The companion still does most of the work —
it removes about nine tenths of the cross-observable stellar-mass inconsistency that
ordinary matter alone suffers — but a real 10 percent discrepancy survives. The
useful news is that a spherical redistribution of the companion, staying inside its
own mass budget, repairs five of the six lenses at essentially zero kinematic cost,
while a stronger companion amplitude is firmly blocked. That is a precise
instruction for the causal model: change the shape, not the strength.

## A. Kinematics predicts the lens, and the lens predicts the kinematics

`code/cross_predict.py` → [`run-v1/cross-prediction-results.json`](run-v1/cross-prediction-results.json)

The ten universal field/source coefficients and the two shared population scales
stay at their published R10 values. Per lens the only freedom is the ordinary
nuisance pair R10 already declared: a stellar log-mass offset with its published
uncertainty as prior, and a constant orbital anisotropy. All six systems are
treated identically, unlike JR-1's fitted/excluded split.

### Direction K → L: predict the Einstein angle from the resolved V_rms bins

| Lens | Observed | Predicted (median) | 68% band | Error | Kinematic sigma |
|---|---:|---:|---|---:|---:|
| J0037-0942 | 1.530 | 1.6051 | 1.5866 – 1.6248 | **+4.91%** | 0.0193 |
| J1112+0826 | 1.490 | 1.2711 | 1.2468 – 1.2957 | **-14.69%** | 0.0245 |
| J1204+0358 | 1.310 | 1.3840 | 1.3684 – 1.4007 | **+5.65%** | 0.0162 |
| J1402+6321 | 1.350 | 1.4102 | 1.3930 – 1.4276 | **+4.46%** | 0.0173 |
| J1621+3931 | 1.290 | 1.1450 | 1.1281 – 1.1630 | **-11.24%** | 0.0178 |
| J1630+4520 | 1.780 | 1.5051 | 1.4861 – 1.5264 | **-15.44%** | 0.0203 |

Angle fractional RMS **10.46%**, against **7.94%** for R10's published joint fit and
**48.19%** for the nuisance-matched ordinary-matter control. The kinematic posterior
is tight — the predicted angle carries only a 1.2 to 1.6 percent credible width — so
this is an informative test, not a vacuous one.

The comparison that matters: **the published joint fit was absorbing roughly half of
the discrepancy into the nuisance parameters.** R10's over-bent lenses read +2.6% when
both observables were fitted together and +4.5 to +5.7% under cross-prediction; its
under-bent lenses read -9.5 to -11.9% jointly and -11.2 to -15.4% under
cross-prediction. The published 5.26% and 11.57% figures overstate how well one
observable anticipates the other.

### Direction L → K: predict the V_rms profile from the Einstein angle alone

| Lens | Bins | chi2 (K → L) | chi2 (L → K) | Delta chi2 | V_rms error (K → L) | (L → K) |
|---|---:|---:|---:|---:|---:|---:|
| J0037-0942 | 6 | 6.08 | 25.72 | +19.64 | 0.79% | 2.52% |
| J1112+0826 | 7 | 25.49 | 87.01 | +61.52 | 2.17% | 4.81% |
| J1204+0358 | 7 | 7.51 | 37.44 | +29.93 | 0.64% | 2.79% |
| J1402+6321 | 7 | 15.75 | 29.55 | +13.80 | 1.51% | 2.82% |
| J1621+3931 | 6 | 9.05 | 46.37 | +37.33 | 1.87% | 5.18% |
| J1630+4520 | 7 | 5.40 | 124.57 | +119.17 | 0.94% | 6.39% |

Forcing the model to reproduce the lens costs **Delta chi2 = 281 over 40 bins** and
degrades the mean fractional V_rms error from about 1.3% to about 4.1%.

### The consistency statistic: does the lens want a heavier galaxy than the stars do?

`Delta_dm = dm(lensing) - dm(kinematics)`, in dex and in published stellar-mass sigma:

| Lens | sigma_M | dm from motions | dm from lens | Delta | In sigma | Mass factor |
|---|---:|---:|---:|---:|---:|---:|
| J0037-0942 | 0.060 | +0.0961 | +0.0566 | -0.0395 | **-0.66** | 0.91 |
| J1112+0826 | 0.090 | +0.1103 | +0.2393 | +0.1291 | **+1.43** | 1.35 |
| J1204+0358 | 0.070 | +0.1605 | +0.1169 | -0.0436 | **-0.62** | 0.90 |
| J1402+6321 | 0.070 | +0.0759 | +0.0413 | -0.0346 | **-0.49** | 0.92 |
| J1621+3931 | 0.060 | +0.0061 | +0.1041 | +0.0980 | **+1.63** | 1.25 |
| J1630+4520 | 0.070 | -0.0400 | +0.0980 | +0.1380 | **+1.97** | 1.37 |

Against the ordinary-matter control run through the identical protocol:

| | Companion (R10) | Baryons only |
|---|---:|---:|
| Mean `Delta_dm` in published sigma | **+0.54** | **+5.46** |
| Range | -0.66 to +1.97 | +3.86 to +7.85 |
| Extra stellar mass the lens demands | 0.90x – 1.37x | 1.93x – 2.96x |
| K → L angle fractional RMS | 10.46% | 48.19% |

This is the clearest quantitative statement JR-9 produces. **Ordinary matter needs
the same galaxy to be twice to three times heavier when it explains the lens than
when it explains the stars, a 4 to 8 sigma inconsistency in every system. The
companion, with universal constants frozen and no per-lens gravity parameter,
reduces that to 0.5 to 2 sigma and changes its sign from system to system.** It does
not remove it.

The surviving pattern is bimodal and matches the sign of the angle error exactly:
three lenses need slightly *less* stellar mass for the lens than for the stars, three
need 1.4 to 2.0 sigma more.

### Is the discrepancy just an unstated lensing error bar?

The catalog SIE angles ship no uncertainty, and JR-1's 5 percent was an optimization
scale. Scanning the assumed fractional angle uncertainty:

| Lens | 1% | 2% | 3% | 5% | Needed for 1 sigma |
|---|---:|---:|---:|---:|---:|
| J0037-0942 | -3.18 | -2.16 | -1.57 | -0.99 | 4.9% |
| J1112+0826 | +7.52 | +5.59 | +4.23 | +2.75 | 14.4% |
| J1204+0358 | -3.70 | -2.50 | -1.81 | -1.14 | 5.7% |
| J1402+6321 | -2.86 | -1.96 | -1.43 | -0.90 | 4.5% |
| J1621+3931 | +6.50 | +4.55 | +3.35 | +2.13 | 11.0% |
| J1630+4520 | +10.05 | +6.62 | +4.75 | +2.97 | 15.2% |

The three under-bent lenses stay at 2.1 to 3.0 sigma even with a generous 5 percent
allowance, and need 11 to 15 percent to become consistent. No plausible image-model
uncertainty covers that. The discrepancy is real.

A note on reading the tables below: sections A, B and D quote the angle error at
slightly different reference points for the same nuisance pair — the posterior median,
the posterior mode, and a direct penalty minimum respectively. They differ by a few
tenths of a percent, which is the size of that choice, not a discrepancy.

## B. The residual sorts by redshift, and standard cosmology does not explain it

`code/geometry_check.py` → [`run-v1-geometry/geometry-degeneracy.json`](run-v1-geometry/geometry-degeneracy.json)

Ranking the six by cross-predicted angle error and by lens redshift gives the same
order. Rank correlation **-0.94**, Pearson **-0.93**, with a clean split at
z ≈ 0.22: the three over-bent lenses sit at z = 0.164, 0.196, 0.205; the three
under-bent at z = 0.245, 0.248, 0.273.

Because the project uses a declared static-Euclidean distance convention rather than
a measured one, the first thing to check is whether this is a geometry artefact:

| Lens | z | Error | Needed change in D_ls/D_s | Project | FLRW(0.3, 70) | FLRW/project |
|---|---:|---:|---:|---:|---:|---:|
| J0037-0942 | 0.1955 | +4.44% | **-4.63%** | 0.6355 | 0.6532 | 1.028 |
| J1112+0826 | 0.2730 | -14.46% | **+19.21%** | 0.5057 | 0.5234 | 1.035 |
| J1204+0358 | 0.1644 | +5.63% | **-6.08%** | 0.6887 | 0.7055 | 1.024 |
| J1402+6321 | 0.2046 | +4.45% | **-4.56%** | 0.5263 | 0.5433 | 1.032 |
| J1621+3931 | 0.2449 | -11.63% | **+14.27%** | 0.5352 | 0.5531 | 1.034 |
| J1630+4520 | 0.2479 | -15.83% | **+21.11%** | 0.6208 | 0.6388 | 1.029 |

The required correction runs from -6% to +21% and correlates with lens redshift at
r = +0.93; a linear trend in z cuts its RMS from 13.5% to 4.5%. But a standard flat
FLRW geometry differs from the project convention by only 2.4 to 3.5 percent, almost
uniformly, and correlates with what is needed at only r = +0.57.

**Switching to standard cosmology would not fix this.** Averaged over the six systems
the needed change is 11.6 percent against a 3.0 percent convention difference — about
four times larger, and far steeper in redshift, since the convention offset is nearly
uniform while the requirement swings from -6% to +21%. That weakens the obvious ordinary explanation without removing the
general degeneracy between a redshift-dependent geometry error and a companion with
history.

## C. A universal history term does not rescue it

`code/history_term.py` → [`run-v2-history/history-term.json`](run-v2-history/history-term.json)

Since the hypothesis says the companion has a history, the cheapest test is one
universal coefficient, `A_chi -> A_chi (1+z_lens)^s`, with everything else frozen.
It is unusually clean to test because all 149 SPARC rotation galaxies sit at z ≈ 0,
so no rotation curve moves for any value of s.

Fitting s to the six cross-predicted angle residuals gives **s = 2.82**, a factor 1.87
boost at z = 0.25, and improves the angle RMS from 10.58% to 7.27%. It is rejected
anyway: it buys that only by driving stellar masses to -0.25 to -0.35 dex, at or
through the 5-sigma edge of their published uncertainty, and by degrading the V_rms
fit from about 1.3% to 4.2% and 6.3% in the two worst systems.

**A stronger companion is not the answer.** The resolved kinematics are measured to
1 to 2 percent, and they pin the amplitude. This is the central constraint the causal
model has to respect.

## D. But a redistribution is nearly free — and that is where the room is

`code/response_kernels.py` → [`run-v2-kernels/response-kernels.json`](run-v2-kernels/response-kernels.json)

Both the deflection integral and the Jeans pressure integral are linear in the force,
so adding a thin spherical shell of mass at radius r gives the exact linear-response
kernel of each observable, and any spherical redistribution is a superposition. For
each lens we solve for the shell mass that closes the Einstein-angle gap exactly, then
let the ordinary nuisance pair re-adjust in its own defence and measure the cost.

One methodological result first: **without a mass budget this test is vacuous.** A
shell far outside the light distribution is nearly invisible to interior kinematics
while still bending light, so the cost falls to zero as the radius grows — at the
price of a required mass that diverges (one lens reaches 4565 stellar masses at 64
half-light radii). Any such analysis has to constrain the budget or it proves nothing.

Constraining the shell to the frozen companion's own equivalent mass `A rt / G`, so
the change is a redistribution rather than new mass:

| Lens | Gap | Best shell radius | Shell mass | Of budget | Delta chi2 | Delta penalty |
|---|---:|---:|---:|---:|---:|---:|
| J0037-0942 | +4.29% | 5.82 Re | -21.2 M* | -0.94 | -0.92 | **-0.58** |
| J1112+0826 | -14.61% | 2.02 Re | +3.3 M* | +0.23 | -16.76 | **-16.92** |
| J1204+0358 | +5.68% | 4.36 Re | -5.1 M* | -0.43 | -1.15 | **-0.47** |
| J1402+6321 | +4.63% | 0.53 Re | -0.08 M* | -0.00 | -4.44 | **-3.49** |
| J1621+3931 | -11.69% | 0.10 Re | +0.19 M* | +0.01 | -3.18 | **+0.98** |
| J1630+4520 | -15.81% | 3.96 Re | +11.6 M* | +0.96 | +7.62 | **+7.94** |

**Five of six lenses can be brought to their measured Einstein angle by a spherical
redistribution inside the companion's own budget at no kinematic cost at all, and one
of them improves the kinematic chi-square by 16.8 in the process.** Only J1630+4520
resists, at Delta penalty +7.9.

So the answer to the question the portfolio has been circling — can projected lensing
be changed without disturbing the local force on stars? — is yes, on the real lenses,
quantitatively, and the required redistributions are now written down. The two
observables are not degenerate even though the Einstein radius sits inside the
kinematically measured region in all six systems (theta_E / Re runs 0.46 to 0.92).

## E. The screening statistic on the real joint vector

`code/joint_screen.py` → [`run-v1-screen/joint-screen.json`](run-v1-screen/joint-screen.json)

The starter package's diagnostic, vendored byte-identical (SHA-256
`1d43056f…88741bc02`, matching its manifest), run on the real 46-element vector: 40
seeing-convolved V_rms values in km/s plus six Einstein angles in arcsec, with each
lens's released kinematic covariance and all twelve `(dm, beta)` nuisance modes
marginalized.

| Assumed angle sigma | sqrt(t'Wt) | A (common) | A (motion) | A (lensing) | Difference |
|---:|---:|---:|---:|---:|---:|
| 1% | 10.21 | +1.32 ± 0.10 | +0.55 | +1.23 | -5.05 sigma |
| 2% | 8.08 | +1.26 ± 0.12 | +0.66 | +1.28 | -4.09 sigma |
| 3% | 6.82 | +1.18 ± 0.15 | +0.73 | +1.32 | -3.36 sigma |
| 5% | 5.71 | +1.06 ± 0.18 | +0.80 | +1.36 | -2.41 sigma |
| 10% | 5.04 | +0.93 ± 0.20 | +0.85 | +1.39 | -1.35 sigma |

The candidate shell was sized by requiring the angle to match, so its lensing
amplitude is fitted by construction and is not evidence. The informative entries are
the other two.

`sqrt(t'Wt)` between 5 and 10 says the proposed correction is **comfortably
detectable** in this data — the dataset is not inconclusive. And **A(motion) lands
between +0.55 and +0.85, not at zero**: the resolved stellar motions independently
prefer a substantial share of the same correction that repairs the lens, with the
same sign. The two observables disagree about its size by 1.4 to 5 sigma depending
on the lensing allowance, so the template's radial shape is not yet right — but this
is the opposite of the failure mode where a lensing fix leaves no kinematic trace.

## F. One interaction in place of two fitted rate laws

`code/derive_response.py` → [`run-v3-derived/derived-response.json`](run-v3-derived/derived-response.json)

JR-5 through JR-8 chose `k_plus ~ exp(2h)` and `k_minus ~ exp(4h^2)` independently.
Both now follow from one reversible event, `gamma + M <-> gamma' + chi + M'`, through
Fermi's golden rule with a factorized material operator. A Galilean boost gives
`S_v(Q, omega) = S_0(Q, omega - Q.v)` exactly, so the medium's velocity enters only
through `Q.v`. Direction is derived; no orientation multiplier appears anywhere.

**The scale correction that matters.** Q is set by the photon, not by the cloud, so
`Q.v` is the ordinary optical Doppler shift, of order 1e12 rad/s at 200 km/s, while
collisional dephasing in any astrophysical gas is of order 1e-8 rad/s — a ratio near
1e-20. The Lorentzian is therefore a delta function against the medium's own velocity
spread, and the physical rate is the Lorentzian averaged over the velocity
distribution, exactly as an absorption line samples a velocity profile:

    k = 2 pi u^2 (c / omega) f(v_res - v_bulk . nhat),    v_res = c Delta_0 / omega

Verified numerically against the full velocity integral to 7e-4 relative, and the
starter's reduced two-state example reproduced against the exact Bloch slow pole to
0.155 percent.

Three things follow with no further choices:

**1. The two exploratory laws are one kernel.** Writing `h` for the resonance offset
in profile widths, the derived shape is `exp(-h^2/2)`. Fitting `a exp(b h^2)` to it
recovers `b = -0.500000` to machine precision, which validates the algebra. Fitting
`a exp(b h)` works only locally: to 5 percent over an interval of half-width 0.34 in
`h`, with an exponent set by where you centre it rather than by any physics. So
`k_plus ~ exp(2h)` was a local approximation with no stated range, and
`k_minus ~ exp(4h^2)` **flatly contradicts the derived law** — opposite sign, eight
times the magnitude. The derived kernel says a resonance offset suppresses the return
channel; the fitted law said it amplifies it. That is a substantive disagreement, and
it is testable.

**2. Enhancement and suppression from one process.** For a Gaussian velocity profile,
`d ln k / d ln sigma_v = (dv/sigma_v)^2 - 1`. Extra dispersion raises the rate when
the resonance sits outside the profile and lowers it when it sits inside, with the
turning point at `|dv| = sigma_v`. The numerical scan brackets the sign change at
134.9 – 142.2 km/s against the analytic prediction of 140.0 km/s. Unlike the earlier
`Gamma = |Delta|` framing, this threshold is in velocity dispersion — a quantity
resolved CO cubes measure directly.

**3. The production traces velocity, not density.** The angular pattern peaks exactly
where the projected flow equals the resonance velocity. In a rotating disk the kernel
therefore selects an iso-velocity contour, so the companion production pattern must
rotate with the kinematic major axis and not with the surface-brightness or
column-density pattern. That is a morphological prediction a CO cube can check.

### The independent consequence, declared in advance

The conversion is velocity selective: light is preferentially lost along sightlines
whose gas lies within `sigma_v` of `v_res` in projected velocity. The observable is a
correlation between residual dimming and the CO or H-alpha line-of-sight velocity
field **at fixed column density**, with a Gaussian profile of the locally measured
width. An ordinary dust or opacity effect correlates with column and not with the
velocity offset at fixed column; this kernel predicts the opposite ordering. No
gravitational quantity enters, so nothing that improves a lensing or rotation fit can
tune it. The data route is the PHANGS-ALMA, PHANGS-MUSE and HST/JWST common footprint.

### How strong the selectivity comes out, and why that is useful

With a single sharp `Delta_0` the on- to off-resonance ratio reaches many orders of
magnitude, so conversion would occur only in a thin projected-velocity slice. That is
the honest consequence of one narrow resonance and it is discriminating: if resolved
data show a broad smooth dependence on projected velocity instead, a single sharp
`Delta_0` is excluded and the companion must carry a spread of energies. The
sign-change condition survives that generalisation, because it needs only that the
profile have a finite width.

## What JR-9 does not establish

* The matrix element is not derived from a specified companion identity, so there is
  no absolute rate, no absolute resonance velocity and no energy budget.
* Nothing connects section F to sections A–E. **The redistribution the lenses require
  is not yet shown to be what the derived kernel produces**, and until that
  calculation exists, the two halves of this package are a constraint and a mechanism
  that have not been joined.
* No gravitational source, stress tensor or force law follows from a conversion rate.
* The lens calculation is still six catalog SIE angle summaries, not an image
  likelihood, and the kinematics are still azimuthally rebinned radial profiles under
  a spherical non-rotating Jeans model.
* The static-Euclidean distance convention is unchanged and unvalidated. Section B
  reports the degeneracy; it does not resolve it.
* Six systems, all historically examined. The redshift split is a 1-in-20 arrangement
  by chance. It is a pointer, not a detection.

## The reserved test this sets up

The KCWI release carries 14 SLACS lenses. Eight are in this project and only six are
usable; the rest have never been scored here and would need light components and
population masses acquired first. They are a genuine holdout for the redshift trend
in section B, and the way to keep them one is to acquire those inputs and declare the
prediction **before** computing a single angle. The prediction to declare is specific:
under frozen R10, cross-predicted Einstein angles should run increasingly negative
with lens redshift, crossing zero near z ≈ 0.22.

## Reproduction

```
cd code
python cross_predict.py    --output-dir ../run-v1          # ~2 min
python response_kernels.py --output-dir ../run-v2-kernels  # ~3 min
python geometry_check.py   --output-dir ../run-v1-geometry
python history_term.py     --output-dir ../run-v2-history  # ~2 min
python joint_screen.py     --output-dir ../run-v1-screen
python derive_response.py  --output-dir ../run-v3-derived
```

Every script refuses an existing output directory, and `cross_predict.py` aborts
before scoring if the R10 reproduction drifts by more than 1e-6. Requires NumPy and
SciPy; no network access. `code/joint_template_test.py` and
`code/exchange_rate_bridge.py` are vendored unmodified from the starter package and
their SHA-256 values match its manifest.
