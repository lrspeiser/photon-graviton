# Several inward-migration variations

The exact one-third retention and existing capture constants remain unchanged. Four post-capture inward transport maps were executed using both distant illumination and the nearest tested source shell, for both geometries and both population proxies. These maps conserve deposited inventory; they do not yet derive binding, support or the energy released during settling.

## What was tested

Let x=r/R, where R=5 Re is the previous diagnostic region; b is a dimensionless migration rate times duration, and s=exp(-b). No universe age is assumed.

| Variation | Where a captured element moves | Interpretation |
|---|---|---|
| Uniform | y=s*x | Every element relaxes inward |
| Partial | Fraction f moves to s*x; remainder stays | Only some captured energy is mobile |
| Supported | y=x below h; y=h+s*(x-h) above h | Settling slows toward a support radius |
| Retention-conditioned | y=exp[-b*(1-eta)]*x | Migration rate depends on the existing retention factor |

These are known exponential-relaxation and distribution-transport constructions. Their application to companion deposits, the mobile fraction and stopping prescription are proposed phenomenology, not unique first-principles physics. The supported map solves dy/db=-(y-h) above h; the uniform map solves dy/db=-y. R-based scaling is a diagnostic assumption, not a gravitational force law. b=0 reproduces no migration.

## Shared rule and transfer to an omitted galaxy

Parameters are shared across six targets within each geometry/population/source branch. A separate leave-one-out exercise fits five and predicts the sixth without changing that fit. These are already-inspected halo targets, not fresh astronomical data. Errors below are RMS cumulative-fraction differences, in percentage points, averaged equally over radius and galaxies; they are not lensing or rotation residuals.

| Geometry | Source | Variation | No migration | Shared fit | Omitted-galaxy transfer |
|---|---|---|---:|---:|---:|
| companion_regular | distant | uniform | 26.74 | 22.18 | 25.42 |
| companion_regular | distant | partial | 26.74 | 16.80 | 20.02 |
| companion_regular | distant | supported | 26.74 | 22.18 | 25.42 |
| companion_regular | distant | retention_conditioned | 26.74 | 21.76 | 24.95 |
| companion_regular | near_1.1R | uniform | 31.71 | 23.05 | 26.34 |
| companion_regular | near_1.1R | partial | 31.71 | 17.29 | 20.51 |
| companion_regular | near_1.1R | supported | 31.71 | 23.05 | 26.34 |
| companion_regular | near_1.1R | retention_conditioned | 31.71 | 22.51 | 25.70 |
| standard_flat_FLRW | distant | uniform | 23.55 | 23.24 | 25.78 |
| standard_flat_FLRW | distant | partial | 23.55 | 20.66 | 24.32 |
| standard_flat_FLRW | distant | supported | 23.55 | 23.24 | 25.78 |
| standard_flat_FLRW | distant | retention_conditioned | 23.55 | 23.25 | 25.61 |
| standard_flat_FLRW | near_1.1R | uniform | 26.23 | 23.87 | 26.74 |
| standard_flat_FLRW | near_1.1R | partial | 26.23 | 20.45 | 24.33 |
| standard_flat_FLRW | near_1.1R | supported | 26.23 | 23.87 | 26.75 |
| standard_flat_FLRW | near_1.1R | retention_conditioned | 26.23 | 23.92 | 26.65 |

All Salpeter cases, fitted parameters, boundary flags, individual targets and curves are retained in the JSON. The nearest-source branch assumes an isotropic source shell at 1.1 R; its existence has not been measured. It is not selected as the actual source distribution.

## What this can establish

Every strictly inward map increases enclosed deposit fraction at fixed radius. It can help a target requiring more central material, but cannot reduce an already excessive central fraction. Shared fitting therefore tests whether one movement rule can balance these incompatible demands. Closely fitting a target-derived halo shape still does not establish the source budget, its true shape, or correct projected lensing. Standard geometry is a comparison only.

The finite-time maps track locations and retain total positive inventory. They do not conserve the full gravitational, kinetic and reservoir energy automatically: settling releases energy that requires an explicit receiving channel, and stopping needs stress or a binding mechanism. The support-radius candidate imposes where motion slows; it does not explain that support. Self-gravity and source histories have not been evolved.

## Findings in plain language

Partial migration is the strongest of the four tested variations. For distant input under our retained geometry, a shared fit moves about 35.8% of the deposited inventory to 15.0% of its former radius, leaving the rest in place. Its cumulative-profile RMS falls from 26.74 to 16.80 percentage points in the shared fit and to 20.01 in omitted-galaxy transfer. These are substantial placement improvements, but not a solution: the improvement comes with worse central overconcentration in three other targets. This new fraction is fitted to normalized halo shapes and must not be substituted for the earlier 0.37% rotation-trained redistribution result.

Uniform settling and retention-conditioned settling provide smaller improvements under the retained geometry. The supported model chooses h=0, so its added stopping-radius parameter collapses to ordinary uniform contraction rather than establishing nonzero support. No bounds were expanded.

Under standard comparison geometry and distant input, none of the variants improves aggregate omitted-galaxy RMS: the no-migration value is 23.55 points and partial migration gives 24.32. The near-source partial case improves that branch from 26.23 to 24.32, still worse than the distant-input no-migration benchmark. Thus this is not a geometry-independent or source-independent solution.

For the retained geometry and distant input, the following maximum cumulative errors show the tradeoff (percentage points):

| Galaxy | No migration | Partial migration, omitted-galaxy prediction |
|---|---:|---:|
| J0037-0942 | 62.25 | 35.84 |
| J1112+0826 | 5.98 | 37.03 |
| J1204+0358 | 49.43 | 21.79 |
| J1402+6321 | 56.57 | 29.72 |
| J1621+3931 | 16.47 | 31.02 |
| J1630+4520 | 7.63 | 36.31 |

The next mechanism question is what independently measured property controls the mobile fraction and where movement stops. Fitting a separate fraction to each desired halo would describe the answer rather than predict it. A gravity-dependent mobility or support law must be specified from source/matter inputs, then transferred without using the held-out halo shape. These results do not yet establish such a law.

## Verification

Zero-migration identity, endpoint inventory, monotonicity and inward cumulative ordering pass. Doubling evaluation-grid resolution changes fitted RMS by at most 0.000261. Three fitting starts and the explicit no-migration endpoint are compared. Bounds were declared before execution and were not expanded after fits. All six broader project goals remain open.
