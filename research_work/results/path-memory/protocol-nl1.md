# NL-1: a response not proportional to its source, the root spectrum
Declared 19 September 2026 before any calculation, on `claude/tests-clusters-lensing-xrhttm` at d2604eb
(CL-2 stage 2 complete). Stage 2's archive and caches are inputs; nothing declared there changes.
Formula provenance follows [research_plan/formula-provenance.md](../../../research_plan/formula-provenance.md).

## The question, in plain language (the owner's communication rule)
CL-2 stage 2 certified that no field proportional to the ordinary matter, blurred by any positive
combination of widths, describes the galaxies: their extra pull grows with mass more slowly than any
such field can. The classic remedy makes the extra pull grow as the square root of the ordinary pull.
Its local form is already in the repository as PM-1's law, g = g_N + sqrt(a* g_N), and in this
programme's own numbers it passes the galaxies' acceptance by a hair on stage 2's consistent source
(21.62 km/s against 21.9) while falling ten times short on the clusters under the corrected comparison
(113 per point against the reference's 11.6). NL-1 asks whether the same root response, applied to the
*written* field at every width rather than to the local force alone, with one set of amplitudes for
every system, can keep the galaxies while bringing the clusters and lenses inside their acceptances.

Why it matters: this is the first of the two routes left open by stage 2 and TF-1, and the one the
mass-to-speed slope points at. What would stay unreliable if skipped: the claim that nonlinear sourcing
is the next ingredient would rest on the slope argument alone, and the cluster problem of root laws
would remain untested in its nonlocal form.

## The family
    g_extra(x) = sqrt(a_0 g_N(x)) + sum_j sqrt(a_j g_j(x)),   a_0, a_j >= 0
with g_N the ordinary Newtonian force (stage 2's consistent source for the galaxies) and g_j the written
force per unit amplitude at width w_j of stage 2's family (the 25 widths, 0.158 kpc to 10 Mpc). The
member a_0 alone is PM-1's law (established: Famaey and Binney's Bekenstein toy function, per the
prior-art audit); the members j > 0 are proposed here, originality unverified. Where a written force is
outward (a negative operator entry) the root term is zero there; the count and size of such entries
are recorded. The clusters use widths of 1 kpc and more on stage 1's 1500-point grid (gate G4 below)
plus the local member; the lenses use all 25 widths plus the local member through CR-2's interface under
the static geometry with the well coupling (light factor two, stars the full force, the extra force
scaling as the square root of the stellar mass inside the root). The lens bend is the three-dimensional
route, (4/c^2) int g_extra(b sec t) b sec t dt.

## Objectives and the joint combination
Galaxies: stage 1's acceptance metric, the equal-galaxy speed error, exact. Clusters: the chi-square
under the release's SZ correlations with the corrected boundary fit of stage 2, the boundary pressures
profiled by nonnegative least squares at each amplitude vector. Lenses: the Einstein constraint at 3%
and the kinematics chi-square in V with anisotropies alternated as in TF-1.
The joint objective is dimensionless and balanced, each block divided by its acceptance so that one
means "at the threshold":
    J = F/(21.9 km/s)^2 + (chi2_c/N_c)/(2 x 11.64) + chi2_K/128 + chi2_E/6
(11.64 is stage 2's release reference per point under the correlations; a block is described when its
term is at most one). The family is not convex in the amplitudes, so no global certificate is claimed:
each solve is a multi-start bounded L-BFGS (20 starts: the local law alone, the local law plus stage 2's
linear spectrum scaled, and eighteen random) followed by the trust-region Newton polish and the free-face
finish of stage 2, and the best point is reported with the spread of the twenty.

## Gates (numerical, each with a control where one exists)
R1 reproduction: the local law at PM-1's fitted a* on the tabulated force reproduces the archived
   galaxies-alone reference 20.214 km/s to 1e-3; the clusters' PM-1 references under the correlations
   reproduce stage 2's 113.06 and 98.92 per point to 1e-9; stage 2's caches load by hash.
G1 optimality: at every reported optimum the projected gradient, read on the finished point, is below
   1e-8 of the objective.
G2 multi-start: the best objective is reached within 1e-6 relative by at least three of the twenty
   starts; a failure marks that optimum "not reproduced across starts" in the reading and is recorded.
G3 the analytic gradient of the root transform agrees with finite differences at 20 random feasible
   points to 1e-6.
G4 cluster grid: the linear-law pressure columns on the 1500-point grid agree with stage 2's 6000-point
   cache to 1e-3 for widths of 1 kpc and more (control: the 0.158-kpc column, which differs by more).
G5 negative operator entries zeroed for the root law: fewer than 1e-3 of all entries.

## Decision rule (fixed now)
(a) One amplitude vector, the joint minimiser, has every term at most one: galaxies at most 21.9 km/s
    with the consistent source, clusters within a factor two of the reference, lenses within 3% and 128.
    The root spectrum is the candidate universal law; next, its energy budget and formation.
(b) The galaxies-alone best is at most 21.9 km/s but no joint minimiser brings the clusters (or the
    lenses) inside: the cluster problem of root laws persists under nonlocal sourcing. The reading
    reports the trade-off: the cluster term at the galaxy optimum and the galaxy term at the cluster
    optimum, and the best joint point.
(c) The galaxies-alone best exceeds 21.9 km/s with the consistent source: the root law's pass on the
    galaxies depended on the tabulated force; reported as such with the tabulated sensitivity.
Validation and test are scored at the training optimum; the outer mass-to-speed slope against the
integrated baryonic mass is reported; transfer as stage 2's E3. Exposed data; the three statuses stay
separate.

## Evidence
`nl1.py` writes `nl1-results.json` with input hashes; `nl1_checks.py` reruns the fast gates and anchors
the archive, registered before the stage-8 job, which stays last. Failed gates are recorded and never
re-thresholded; a corrective rerun needs a declared amendment.
