# CL-2 stage 2: the shape, with every column visible to every block
Declared 19 September 2026 before any calculation, on `claude/tests-clusters-lensing-xrhttm` at 9131422
(the text-only corrections to the stage 1 record), `main` at 9aff6b2. Stage 1's archive
`cl2-results.json`, its protocol and amendments are untouched. Formula provenance follows
[research_plan/formula-provenance.md](../../../research_plan/formula-provenance.md).

## The question, in plain language (the owner's communication rule)
Stage 1 concluded that no member of the universal family of footprints serves galaxies, lenses and
clusters at once. The review of that stage, verified line by line, found that the conclusion was the
output of a procedure with five defects: part of the family was switched off in two of the three
systems; the solver minimised a different quantity from the one the rule judged; two descriptions of the
same galaxy's matter were mixed; one sensitivity's reference carried a factor twice; and the cluster SZ
points were treated as independent when the release supplies their correlations. Stage 2 repairs all
five and asks the sharpened question: evaluated the same way in every system and judged on the
acceptance metric itself, does any member of the family describe the galaxies, and how far is the
corrected cluster comparison from the release's own reference?

Why it matters: TF-1 showed that the direction of the field's action is not the obstruction, so the
shape of the response is the whole remaining question for a universal law. A certified minimum of the
galaxy metric either produces a viable spectrum or a defensible rejection of the entire linear family,
which is what would justify new field physics. What would stay unreliable if skipped: the galaxy
exclusion would remain uncertified, the cluster factor 1.49 would rest on an approximate likelihood, and
every later stage would inherit both.

## Premises
Static Euclidean geometry PF1 (D = ln(1+z)/alpha0) is primary for the lenses, as the programme's
nonexpansion premise requires; flat FLRW is recorded as a sensitivity for continuity with stage 1's
archive. No dark matter and no halo model is fitted; the release's own NFW score for the clusters is
quoted as the release's pipeline consistency number, not as a floor. All data are exposed: every reading
is exploratory. The held-out lens remains unavailable.

## S2-2. One family, evaluated identically everywhere
The family is C = sum_j Lambda_j (rho_b * K_{w_j}), Lambda_j >= 0, on the 21 widths 0.158 kpc to 10 Mpc
(stage 1's grid without the 0.1-kpc column, per amendment 2), for every block: the 12 X-COP clusters
through the forward pressure, the SPARC training galaxies (validation and test scored, never fitted),
the six SLACS lenses, and the Milky Way baselines I and II for transfer only.
- Clusters: n_grid 6000 (stage 1 used 1500) and every width resolved by the shell kernel; no column is
  zero. Gate G4-2, on A1795 and A2319: doubling to 12000 changes the pressure columns at the pressure
  radii by < 1e-4 for widths >= 1 kpc and < 1e-3 for 0.398 to 0.631 kpc. For 0.158 and 0.251 kpc no
  convergence is required; instead gate G4-3, negligibility: at the largest amplitude any fit of this
  stage assigns them, their pressure contribution is < 1e-3 of the Newtonian pressure at every point.
- Galaxies: widths >= 0.158 kpc as amended, n_grid 3000; gate G4b rerun on three galaxies (< 1e-3).
- Lenses: the 21 widths as in TF-1, through CR-2's measurement interface.

## S2-3. One source per galaxy
The Newtonian force is computed from the same reconstructed components as the written force: the
razor-thin stellar and gas disks of RPG-1's construction (Upsilon_disk = 0.5 on the SPARC surface
brightness with an exponential tail; the gas disk as the fitted exponential of 1.33 M_HI) and the
spherical bulge. Disk force by the exact ring kernel (complete elliptic integrals, established) with the
uniform-sheet subtraction: g(R) = int 2 pi s [Sigma(s) - Sigma(R)] F(R, s) ds - Sigma(R) int_{R_max}^inf
2 pi s F(R, s) ds, using the identity that an infinite uniform sheet exerts no in-plane force (derived
here from established mathematics; the singular ring sum without subtraction was found not to converge
and is the control). Gates: G-N1 Freeman's exponential-disk closed form reproduced to 1e-3 at eight
radii from 0.3 to 30 kpc (control: the plain ring sum, off by more than 10%); G-N2 the bulge force
equals G M(<r)/r^2 to 1e-12; G-N3 a compact ring set reproduces the point-mass force to 1e-6 at ten
times its radius. Recorded, not gated: the ratio of the reconstructed to the tabulated SPARC Newtonian
force, point by point. The primary galaxy comparison uses the reconstructed force; the tabulated force
is a sensitivity. The mass in the slope diagnostic is the components' integrated baryonic mass.

## S2-4. The acceptance metric minimised exactly
F(Lambda) = mean over galaxies of the mean over points of (sqrt(u_i) - v_i)^2, u_i = R_i (g_N,i +
A_i . Lambda), Lambda >= 0: the square of stage 1's equal-galaxy RMSE. On u > 0 it is convex, since its
Hessian is sum_i w_i v_i/(2 u_i^(3/2)) R_i^2 a_i a_i^T (established; the review's identity, checked).
Solver A: L-BFGS-B with bounds, analytic gradient, the objective +inf where any u_i <= 0, started from
the nonnegative least-squares solution of stage 1's acceleration problem. Solver B: projected Newton
with the analytic Hessian, an independent code path, started from the origin.
Gates: G-C1 the two solvers agree to 1e-6 relative in F; G-C2 optimality: for active amplitudes the
gradient's magnitude < 1e-8 of its magnitude at the start, for inactive amplitudes the gradient >= -1e-8
of that scale; G-C3 convexity: the smallest directional curvature over 200 random directions at the
optimum and at 20 random feasible points >= 0 (numerical Hessian-vector products); G-C4 reproduction:
the nonnegative least-squares start with the tabulated force reproduces stage 1's archived
galaxies-alone score, 30.99 km/s, to 1e-6.
The certified minimum RMSE* = sqrt(F*) is judged by stage 1's rule, <= 21.9 km/s on training;
validation and test are scored at the training optimum. The slope diagnostic is reported with the
integrated mass. The joint solves of S2-7 use the same convex objective for the galaxy block and
weighted least squares for the others, combined with stage 1's per-block 1/N weights (convex).

## S2-5. Lenses
TF-1 already ran the two-step test on the static geometry: per lens, free amplitudes describe every lens
under either coupling; one common spectrum describes none. Those results are cited, not rerun. Here the
lens block enters the three-block joint solve through TF-1's rows with the well coupling (light factor
two), the factor-one coupling recorded as a sensitivity; and the per-lens ratio-hull intervals of the
review (the ratio of extra second moment to extra bending per width, whose convex hull must contain the
ratio the data need) are recorded as a diagnostic without a gate.

## S2-6. Clusters corrected
The boundary-pressure fit applies the thermal fraction once (stage 1's `fit_boundary` applied it twice
to reference models); stage 1's non-thermal sensitivity references are recomputed. The SZ points' 
covariance is taken from the release: the per-cluster pressure correlation matrix (the COVMAT block of
the Y-PROF-COVMAT files, read from the pinned tarball) applied to the extract's SZ errors, Cov =
corr o (err err^T); the X-ray points keep independent errors. Gates: G-S1 the block's radii agree with
the extract's SZ radii to 2% (1.6% on A2319 and A644 was found in preparation and is recorded);
G-S2 every correlation matrix is symmetric and positive definite with condition number < 1e6;
G-S3 with the identity in place of the correlation matrix, the chi-square of stage 1's archived
spectrum reproduces stage 1's clusters-alone chi-square to 1e-9. The references (baryons only, the
release's NFW, PM-1's law) are recomputed under the covariance and the corrected boundary fit; the
spectrum is refitted; stage 1's sensitivities (mu, no stars, non-thermal fractions) are repeated.

## S2-7. Joint and transfer
Stage 1's E1 (each block alone; all three jointly) and E3 (calibrated on galaxies, on clusters, on both;
scored on the others and on the Milky Way) are repeated with every column visible and the galaxy metric
exact.

## Decision rule (fixed now)
(a) The certified minimum on training with the reconstructed source is <= 21.9 km/s: the family
    describes the galaxies; universality is then decided by the three-block joint solve under stage 1's
    thresholds (clusters within a factor two of the release reference under the covariance comparison;
    lenses by TF-1's acceptance: Einstein residual <= 3%, kinematics chi-square <= 128).
(b) The certified minimum exceeds 21.9 km/s: no member of the universal linear family describes the
    galaxies with these sources, evaluated consistently and judged on the metric itself. That is the
    certified negative the review asked for; the next ingredient is a response not proportional to its
    source or a kernel outside the positive Gaussian mixtures, each in its own protocol.
The cluster reading, whatever the galaxy outcome: the spectrum's score against the corrected references
under the covariance comparison, superseding stage 1's factor 1.49.
The three statuses stay separate: reproduction, numerical verification, scientific reading.

## Evidence and checkpoints
This protocol is committed alone. The stage runs into `cl2s2-results.json` with input and source hashes
and the runtime; `cl2s2_checks.py` reruns the fast gates and anchors the archive, registered before the
stage-8 job, which stays last. Failed gates are recorded and never re-thresholded; a corrective rerun
needs a declared amendment. The report keeps the three statuses apart.
