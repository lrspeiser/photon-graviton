# SE-LR4: refinement checkpoint, campaign still running

20 September 2026. Protocol69c6786; numerical runner4b8472d. One of four
declared runs is complete: n40. The n48,n64 and commonly rotated n48 cases
remain pending. No overall spatial or rotation accuracy result is available.

The n40 run passes its individual finite-state, energy, angular-momentum,
boundary-energy and quadrature gates. Independent reconstruction passes20
checks with maximum endpoint discrepancy2.17e-19. The numerical equations
and parameters are unchanged from SE-LR3; this campaign refines the grid.

| n40 result at model time1 | Value |
|---|---:|
| A-sector energy / initial total energy |20.7146%|
| A-sector angular momentum / initial total angular momentum |10.0936%|
| Signed circulation at radius1.5 |0.000179942557|
| Maximum relative energy drift |4.27e-10|
| Maximum relative total angular-vector drift |1.38e-6|
| Maximum edge-energy fraction |1.56e-14|

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
treated as a demonstrated continuum propagation front. Finer grids and the
rotated-source test are still needed.

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
