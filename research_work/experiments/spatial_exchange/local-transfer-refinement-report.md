# SE-LR4: spatial gates pass, rotation gate fails

20 September 2026. Protocol69c6786; numerical runner4b8472d. All four
declared runs are complete. The campaign fails because circulation changes
2.52915% under a common rotation, above the unchanged2% gate. All individual
run gates pass. The independent complete archive audit passes71 checks,
maximum endpoint discrepancy2.17e-19.

| Final comparison | Difference | Gate | Result |
|---|---:|---:|---|
|n48/n64 circulation|3.59587%|<5%, decreasing|Pass|
|Previous n32/n48 circulation|9.17437%|Trend reference|Decreases|
|n48/n64 field angular momentum|0.153806%|<5%, decreasing|Pass|
|Previous n32/n48 field angular momentum|0.432482%|Trend reference|Decreases|
|Common rotation, circulation|2.52915%|<2%|Fail|
|Common rotation, field angular momentum|9.78e-8%|<2%|Pass|

At n64, circulation at radius1.5 is0.000165696048, A carries10.06273% of
initial angular momentum and20.51141% of initial energy. At rotated n48,
circulation is0.000167312865 compared with unrotated0.000171654259. The field
evolution stencil and loop interpolation both have grid direction dependence;
this comparison does not by itself separate those sources of error. Both
must be investigated without retroactively changing the declared criterion.

The n40 run passes its individual finite-state, energy, angular-momentum,
boundary-energy and quadrature gates. Independent reconstruction passes20
checks with maximum endpoint discrepancy2.17e-19 at the first checkpoint;
the two-run audit now passes36 checks with the same maximum discrepancy.
The numerical equations
and parameters are unchanged from SE-LR3; this campaign refines the grid.

| n40 result at model time1 | Value |
|---|---:|
| A-sector energy / initial total energy |20.7146%|
| A-sector angular momentum / initial total angular momentum |10.0936%|
| Signed circulation at radius1.5 |0.000179942557|
| Maximum relative energy drift |4.27e-10|
| Maximum relative total angular-vector drift |1.38e-6|
| Maximum edge-energy fraction |1.56e-14|

At n48 the circulation is0.000171654259, the A energy fraction20.6133%,
and the A angular fraction10.0782%. All individual gates pass. A-sector
angular momentum outside radius1.2 is0.684015% of initial total angular
momentum; total outer angular momentum is1.754690%. These are still
exploratory masked-region diagnostics, not measured angular surface flux.

Compared with the previous n32 run, circulation decreases further. The actual
acceptance pair is preregistered as n48/n64 with a decreasing difference from
n32/n48; n40 is a trend diagnostic and cannot replace that pair.

## Angular momentum outside the initial excitation

The initial Q/P profile is zero beyond radius1.2. Initially A/Pi is zero
everywhere. At time1 the projected angular-momentum fractions are:

| Outside radius | A sector / initial total | All sectors / initial total |
|---|---:|---:|
|1.2|0.650193%|1.615571%|
|1.5|0.052617%|0.065122%|
|2.0|0.00002652%|0.00002699%|

This is an exploratory measurement of angular momentum accumulated outside
the initial support, not a measured surface flux. Both Q/P and A can spread;
outgoing Q/P can generate A locally outside the original source. Hard radial
masks introduce their own grid error, and the small outer tails must not be
treated as a demonstrated continuum propagation front. The failed rotation
gate still needs resolution. At n64 the A and total angular fractions outside
radius1.2 are0.674001% and1.706395%, respectively. The short run does not
establish sustained angular transport.

This finite prescribed excitation is not yet ordinary-matter generation with
a physical fuel budget. It supplies no stellar orbit or light trajectory and
does not demonstrate sustained gravity or a fit to observations. No dark matter,
expanding background or agreement with an older gravity law is used. All twelve
goal requirements remain active.

Files in local-transfer-refinement-v1 preserve manifest/source hashes and
completed per-run arrays and traces. The auditor reconstructs endpoint energy,
canonical angular momentum, outer-region densities and loops independently of
the evolution ledger and interpolation functions. Intermediate extrema are
checked from recorded traces rather than intermediate full field archives.
