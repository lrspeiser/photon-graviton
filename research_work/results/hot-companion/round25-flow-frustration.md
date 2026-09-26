# Round 25: flow frustration (test registered before it was run)

Registered 26 September 2026, before any value of the measure below was computed. The results go in README
section 36; this file is not edited after the run except to add a pointer to them.

## The proposal (project owner, 26 September 2026)

The law does well where the companion's energy has one clean place to go, and it falls short where that energy is
*directionally frustrated*: several streams that cancel, a larger outside current that overwhelms a smaller one, a
strong-field region that blocks release, or a collision that makes energy faster than the slow stream carries it
away. Round 24 already needed a term of this kind: the Casimir-EFT field equation reproduced the law only after the
hot glow's non-flowing brightness, S_ex = S - |g_hot|, proportional to u U - |J|, was added. The misses the suite
records are one-sided (every close or fail is an underprediction: KiDS red lenses, the Sun's speed, the Galaxy's
mass inside 50 kpc, five dwarfs, the Bullet's subcluster, MACS J0025's and El Gordo NW's galaxy speeds, Abell 520
P6). The proposal's measure, with the law's own fields at each point:

    chi = 1 - |g_N + g_hot| / (|g_N| + S)        ("flow frustration")

|g_N| + S is the companion's total intensity in the law; |g_N + g_hot| is the part that flows as one net stream.
The question is whether the law's residual, observed / law - 1, grows with chi. Nothing is fitted first; only if
the misses line up is a recycling (redistribution) equation derived and run on the suite.

## The test

**Where.** Every SPARC point (149 galaxies), every X-COP radius (12 clusters x 6 radii), every KiDS and Mistele bin,
the Milky Way's rotation curves (four Gaia analyses, 5-27 kpc) and its other graded radii, the ten dwarfs and every
collision aperture, each in the model that makes the suite's prediction there: SPARC with the bulge's heat (spherical
shells); X-COP with deprojected stars; KiDS and Mistele as the point lenses the suite uses; the Milky Way on the
(R, z) grid of McMillan (2017) matter; each dwarf as a Plummer sphere in the Galaxy's pull and heat, averaged over
directions; the collisions on their 3D grids with the companion's memory, where the flowing companion is the memory
field F, so chi = 1 - |F + g_hot| / (|F| + S); the collisions' galaxy speeds in the settled pre-collision clusters
that the suite uses to predict them.

**Residual.** ln(observed / law) of the pull or of a quantity proportional to it (enclosed or aperture mass, speed
squared, dispersion squared). Ranges are taken at their midpoint; lower bounds, positions and differences are left
out.

**Averaging.** Over an aperture or a sphere, chi is weighted by the companion's inward flux through the aperture's
wall (a cylinder for a lensing mass, a sphere for an enclosed mass); for a galaxy-speed check, by the stars inside
the aperture; for a dwarf, over directions and the Plummer mass.

**Statistics.**
1. Per family: the rank (Spearman) correlation of the residual with chi, and the same at fixed Newtonian pull
   (partial rank correlation controlling for ln |g_N|; for X-COP also r / R500).
2. Graded checks with a ratio-type measurement: chi of the passes against chi of the closes and fails
   (Mann-Whitney, one-sided: non-passes higher).
3. Pooled: a permutation test with chi shuffled within each family, of the sum of the within-family correlations.

**Success (all three needed).**
(a) the pooled within-family correlation is positive with p < 0.01;
(b) it is positive in every family in which chi spans at least 0.1, and it survives the control for ln |g_N|
    (pooled partial correlation positive, p < 0.05);
(c) the closes and fails have higher chi than the passes (one-sided p < 0.05).

If (a) fails, the idea in this accounting is not supported ("if the failed systems have no relationship to chi,
this idea dies quickly"). If (a)-(c) hold, a recycling equation is derived and the whole suite rerun.

**Secondary (exploratory, reported but not part of the verdict).**
* chi_all = 1 - |g_N + g_hot| / (S_cold + S), where the ordered (cold) matter's energy is also counted as brightness,
  S_cold = G int rho / d^2 (spherical stand-ins, except on the collision grids where it is exact);
* for the collisions' galaxy speeds, chi of today's collision field around the galaxies (instead of the settled
  cluster the prediction uses).
A positive secondary result would need its own registered test.

**Known before registration (from the code, not from any computed value).** chi is exactly zero wherever the model
has a single centre and no extended hot matter: the KiDS and Mistele point lenses and the SPARC galaxies without a
bulge. There the test can only ask whether the residuals are zero at chi = 0.

**Results** (added after the run): README section 36; data in `run-frustration-v25/`. All three conditions failed; the
idea is not supported in this accounting.
