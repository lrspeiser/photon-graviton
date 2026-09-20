# TF-1: transverse coupling of the converted field
19 September 2026. Executed under the owner's static-universe, no-dark-matter constraints on the
branch `claude/tests-clusters-lensing-xrhttm`. [Protocol](protocol.md) eec8965, declared before any
calculation; [amendment 1](amendment-1.md); [attribution](provenance.md); [reproduction](README.md).
Evidence: exact-v2, scan-v2, scan-v3 and lenses-v1 under [evidence/](evidence/); exact-v1 and scan-v1
are failed-process records (provenance guard tripped by an edit during the run) that the clean reruns
reproduce bit for bit.

**Outcome.** A sideways coupling of the converted field passes, exactly and by construction, the three
requirements CWC-1's scalar channel failed: it changes neither light's speed nor its colour, and it
treats all matter alike. It bends light along the same paths as the index coupling, because ray
curvature under an index is the transverse gradient of its logarithm. It converts no energy, so every
unit of conversion has to come from the along-the-motion channel at a fixed exchange rate measured
here. On the six lenses it leaves CL-2's verdict unchanged, and per lens it is indistinguishable from a
well: **outcome (b)**. The magnetism-like vector coupling is rejected by its geometry. The obstruction to
one universal law is the field's shape, not the direction of its action.

| status | finding |
|---|---|
| **reproduction** | **passes.** CWC-1's archived 1024-point runs (primary, three vacua, coupling 0.05) reproduced to 2e-15; CWC-1's archived 2D rays reproduced by the live index-ray code to better than 1e-9; CL-1's archived Newtonian Einstein angles under the static geometry reproduced through RPG-1's three-dimensional route to 8.5e-9 (the projected route differs by 4.4e-6, as CL-1 recorded); CL-2's stars-only benchmarks under flat FLRW to 1.7e-7. |
| **numerical verification** | **stage 1 passes 18 of 18 gates with every control rejected; stage 2 passes energy and reproduction and FAILED its declared power-law gate (scan-v2), which passes under amendment 1's weak-regime refit (scan-v3), both quoted; stage 3 passes 3 of 3** (reproduction, monotone alternation, optimality residuals below 5e-16). |
| **scientific outcome** | **(b), by the rule fixed before the run.** The steering law passes G1 to G9 but under the joint universal family the lenses are described neither by the well (kinematics chi-square 327 where 128 was allowed, Einstein residual 13.1% where 3% was allowed) nor by the steering law (435, 15.2%). The per-lens ceilings are indistinguishable (ratio 1.015). The vector law is rejected by G7 and G8 as declared. Exposed data, static geometry, mean-field kinematics: a description, not a validated theory. |

## Plan-to-execution record

| Stage | Executed evidence | Result |
|---|---|---|
| Exact properties | 2D rays on the archived field with the live index rays; probes; Lagrangian family at four exponents; 1e6-sample velocity averages; 48 rays through the model vector field; orbits and a 2000-probe ensemble | 18 of 18 numerical gates; steering exact on light, vector law rejected |
| Longitudinal scan | 20 new 1D runs (six couplings, three colours, three vacua, one reproduction) plus two archived cases | energy and reproduction pass; power-law gate failed as declared, refit under amendment 1 passes |
| Lenses | 12 lens systems (six lenses, two geometries), 48 per-lens feasibility fits, six joint alternating solves | 3 of 3 numerical gates; outcome (b) |

## Stage 1: the exact properties (evidence/exact-v2)

| Gate | Measured | Threshold | Control |
|---|---|---|---|
| G1 work-free | steering: unit speed by construction, frequency not a variable; vector: speed drift 2.2e-14 | < 1e-9 | index rays: speed deviates by 4.06e-3; rejected |
| G2 composition | steering and vector accelerations identical for action-to-mass ratio 2 (difference 0) | < 1e-14 | CWC-1 material force ratio 2.000; rejected |
| G5 archived field | live index rays reproduce the archived angles | < 1e-9 relative | |
| G5 inward | all six steering rays inward | | zero field: exactly straight |
| G5 against the index bending | relative difference 2.4e-8 | < 2% in the weak regime (g max phi = 0.0041) | |
| G5 ray step | 4.2e-8 of the largest bend | < 1e-4 | |
| G6 Euler-Lagrange | residual at most 5.3e-9 at alpha 0, 0.5, 1, 2; formula errors at most 2.2e-15 | < 1e-6 | alpha = 0 is a potential force, speed changes |
| G6 speed | only alpha = 1 conserves speed (2.2e-15); the others change it by 0.97 to 0.99 of the acceleration | | |
| G7 vector, isotropic | mean force 2.1e-3 sigma B at N = 1e6 | < 4.2e-3 (three standard errors) | co-rotating population: v_rot B to 2.9e-4 |
| G7 steering support | quadrature and Monte Carlo agree to 5.1e-4; f(0) = 2/3 to 1.1e-16 | < 2e-3 | |
| G8 model field | divergence 1.7e-11 of the field scale | < 1e-6 | |
| G8 focusing | transverse 0.99997 in every view; vector -0.0007 | > 0.999; below 0.1 | deflections at most 0.020 radian |
| G9 circular orbit | radius drift 6.4e-11; period equal to the well's to 1.6e-11 | < 1e-8 | |
| G9 radial orbit | direction and speed change exactly 0 | < 1e-12 | |
| G9 ensemble | conservation drift per orbit 2.9e-13 (steering), 9.6e-14 (well) | < 1e-8 | |

**What the steering law is.** For a ray Hamiltonian H = a(x)|p| the direction obeys dn/ds = P_perp grad ln(1/a)
per unit path length; the speed a cancels from the path (Fermat's principle, established). The steering
law dn/dt = g P_perp grad phi with g_perp = g is therefore the index law's path at unit speed, which the
2.4e-8 agreement confirms to the integrator's tolerance. Every timing and colour cost of CWC-1 sits in
the speed and clock factors the steering law drops, so the three failed screens are passed identically,
not approximately. The price is structural: (i) no Lagrangian of the form |v|^2/2 + psi|v|^alpha
produces a speed-independent transverse force (G6: speed is conserved only at alpha = 1, where the
force scales with speed), and the derived ray argument in the protocol shows no achromatic Hamiltonian
ray law keeps unit speed while bending; (ii) with no Hamiltonian there is no reciprocal recoil term, so
the momentum the law takes from light and matter is deposited nowhere; (iii) it does no work, so it
converts nothing; (iv) it does not turn radial motion at all, and it drives an isotropic population
radial: anisotropy from -0.007 to +0.23 in 20 dynamical times where the well drives the same population
to -0.38, the median radius shrinking to 1.06 against the well's 1.22 and the largest excursion reaching
5.8 against 2.4. Nothing escaped within 20 dynamical times, but the orbits it "keeps" are the tangential
ones, converted steadily into radial ones it cannot hold.

**What the vector law is not.** With the Lorentz structure a = v x B the mean force on a non-rotating
population is zero (2.1e-3 sigma B at a million samples, within three standard errors of zero), the
deflection of light has the sign of the field rather than of the side of passage, so rays on opposite
sides deflect the same way and nothing focuses (-0.0007 across three viewing directions), and at equal
field it pushes light c/v times harder than stars: 1,499 at 200 km/s, 999 at 300 km/s. It cannot
support a pressure-held galaxy or cluster, cannot produce mass-like lensing, and cannot serve light and
stars with one field. No sourcing question needs to be asked of it.

## Stage 2: the longitudinal scan (evidence/scan-v2; refit in scan-v3)

The transverse channel is silent in one dimension, so every cost here is the along-the-motion channel's.
CWC-1's 1D fixture at 1024 points, clock exponent 2; energy error below 1e-3 in every run.

| g_par | receiving gain | timing discrepancy | colour spread | fixed-ruler speed change | screens: timing 1%, colour 1%, speed 0.1%, conversion 0.1% |
|---|---|---|---|---|---|
| 0.005 | 0.088% | 0.020% | 0.003% | 0.007% | conversion fails |
| 0.010 | 0.351% | 0.080% | 0.010% | 0.027% | **all four pass** |
| 0.020 | 1.38% | 0.315% | 0.039% | 0.107% | speed fails |
| 0.050 | 7.64% | 1.64% | 0.197% | 0.628% | timing and speed fail |
| 0.100 | 18.2% | 3.15% | 1.45% | 2.08% | all three fail |
| 0.200 | 17.9% | 1.50% | 3.29% | 6.53% | all three fail (CWC-1's primary) |
| 0.500 | 27.1% | 34.6% | archived, one colour | 32.2% | fails |

In the weak regime (0.005 to 0.02) gain, timing, colour and speed all scale as g_par to the power 1.99,
1.99, 1.97 and 1.99 with log scatters below 0.004 (scan-v3). So the costs per unit of conversion are
fixed there: **each 1% of the light's energy converted costs 0.23% in timing discrepancy, 0.03% in
colour spread and 0.08% in fixed-ruler speed.** The window in which every screen passes is narrow,
between about 0.007 and 0.02 in this fixture, and delivers 0.2% to 1.4% conversion. Over 0.005 to 0.1
the same quantities do not follow one power law (the gain saturates near 0.1 and the timing metric is
not monotone, which CWC-1's own resolution study already showed to be grid-sensitive); that is the
declared gate that failed. The observed cosmological rate is H0/c, about 2.3e-4 per Mpc; no mapping of
the fixture's units to megaparsecs is made, but the exchange rate says that a coupling tuned to that
rate would sit far below every local screen.

## Stage 3: the six lenses (evidence/lenses-v1)

Static Euclidean geometry as primary (the owner's no-expansion constraint), flat FLRW as a recorded
sensitivity; CR-2's measurement interface through CL-2's `LensSystem`; the universal family of 21
footprint widths from 0.158 kpc to 10 Mpc on each lens's deprojected stars; anisotropy free in
[-2, 0.45] per lens. The well W bends light with the factor two and moves stars with the full force; the
steering law S bends light with factor one and moves stars with f(beta) = 1 - <v_r^2/v^2>: 0.798 at
beta = -2, 0.754 at -1, 0.667 at 0, 0.584 at 0.45.

**Joint, one spectrum for six lenses.**

| geometry, IMF | coupling | kinematics chi-square (40 bins, 128 allowed) | worst Einstein residual (3% allowed) | described |
|---|---|---|---|---|
| static, Chabrier | W | 327 | 13.1% | no |
| static, Chabrier | S | 435 | 15.2% | no |
| static, Salpeter | W | 330 | 13.8% | no |
| static, Salpeter | S | 389 | 16.1% | no |
| flat FLRW, Chabrier | W | 408 | 13.1% | no (CL-2's archived joint, 407.7 and 13.1%, recovered) |
| flat FLRW, Chabrier | S | 493 | 14.8% | no |

Per lens in the static Chabrier joint fits, kinematics chi-square and Einstein residual: J0037 108 and
-4.7% (W), 116 and -4.8% (S); J1112 29 and -12.2%, 36 and -13.9%; J1204 5.6 and +3.4%, 12.0 and -1.7%;
J1402 40 and +7.3%, 44 and +5.7%; J1621 7.8 and -8.5%, 10.0 and -9.2%; J1630 137 and +13.1%, 217 and
+15.2%. The alternating solves converged monotonically in six iterations (objectives 22.0 to 17.1 for W,
27.4 to 21.7 for S). The active widths: W at 0.63, 1.0, 3.98, 25 and 40 kpc; S at 0.40, 0.63, 2.5, 3.98,
63 and 100 kpc. The steering law's smaller light factor makes it want more force, which the stars then
feel too strongly unless the anisotropy turns radial (its fitted betas run from -0.86 to +0.42 against
the well's -2.0 to +0.02); the trade does not close the gap.

**Per lens, free amplitudes (feasibility ceilings).** With 21 widths or 16 shells free per lens against
six to eight bins, every lens is fitted to a total kinematics chi-square of 13 to 15 over 40 bins with
the Einstein constraint met, under either coupling and either family; the S-to-W ratio of the ceilings
is 1.015, indistinguishable by the declared rule. The problems are underdetermined, so this is a
statement that the data cannot tell the two light-to-star ratios apart, not that either is right. Most
per-lens fits sit at the tangential anisotropy bound (six of six under the universal family, four or
five of six under the shells).

**The vector law** was not fitted: G7 gives zero mean support for a non-rotating population and G8 no
focusing; the ratio c/v is recorded as the reason.

## The dichotomy, as it now stands (derived here; originality unverified)

A transverse force that bends light closer and keeps stars in orbit with the observed ratio must be
speed-independent: a speed-dependent one gives light c/v times the push it gives stars, does not focus,
and has no mean effect on pressure-held matter. A speed-independent transverse force has no Lagrangian
in the family tested and no reciprocal recoil, converts nothing, and for light paths and circular orbits
is a well with light factor one. So any reciprocal response that does both jobs is a well, possibly with
a tensor light factor between one and two, and its remaining freedom is the well's shape and that
factor. The conversion channel can coexist with it at weak coupling, at the exchange rate above. That
returns the programme to the shape, which is CL-2 stage 2, now with two extra facts in hand: the light
factor is not pinned to two, and a transverse law's star factor is f(beta).

## What is not shown
No reciprocal formulation of the steering law; no sourcing, formation or retention of a steering profile
at any strength; a frozen periodic 2D fixture for the rays; mean-field kinematics for a non-conservative
force, whose steady state G9 shows is not the isotropic one assumed; underdetermined per-lens problems;
no galaxy or cluster block (for circular orbits and isotropic gas the steering law is the well with
factors 1 and 2/3, so CL-2's shape verdicts carry over unchanged); no Coma; exposed data; the static
geometry is a conditional premise. Passing a gate is a property of the declared candidate in the
declared fixture, not priority and not empirical gravity.

## Next decision
Return to the shape with CL-2 stage 2 (the review's six items), carrying the light factor and f(beta) as
declared options. Keep nonlinear sourcing as the candidate for the mass-to-speed slope. If the steering
law is to be pursued as physics rather than as a probe, the theoretical question is a reciprocal
formulation: a field term that receives the momentum the law takes, with the energy budget closed.

## Correction (19 September 2026; text only)
The universal family in stage 3 is described above and in the protocol as "21 widths from 0.158 kpc to
10 Mpc". It has 25 widths, five per decade, and every computation used all 25 (the archived lens rows
have 25 columns). The count was a miscount in the text; no number in the evidence changes.
