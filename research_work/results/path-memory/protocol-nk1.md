# NK-1: a kernel outside the positive Gaussian mixtures, shell footprints
Declared 19 September 2026 before any calculation, on `claude/tests-clusters-lensing-xrhttm` at d2604eb.
Stage 2's archive and caches are inputs; nothing declared there changes. Formula provenance follows
[research_plan/formula-provenance.md](../../../research_plan/formula-provenance.md).

## The question, in plain language (the owner's communication rule)
Stage 2 exhausted the fields that are a positive blur of the matter: every kernel whose strength only
falls with distance. TF-1 noted the one linear family left outside: a footprint written at a distance
from each mass element, a shell of radius d_0 and thickness w around it, which is linear in the source
but not a positive Gaussian mixture. NK-1 asks whether adding such shells to the family lets a linear
field describe the galaxies. Galaxies first, because a linear field that fails them is closed by the
certified machinery of stage 2 without any further calculation; clusters and lenses follow only if the
galaxies pass.

Why it matters: it is the second of the two routes left open, and the only linear one. What would stay
unreliable if skipped: "kernel outside the positive Gaussian mixtures" would remain an unexamined
loophole in the certified negative of stage 2.

## The family
The 25 Gaussian widths of stage 2 plus 18 shell footprints K(d) = exp(-(d - d_0)^2/2w^2) with
d_0 in {1, 2, 5, 10, 20, 50} kpc and w in {d_0/10, d_0/4, d_0/2}, all amplitudes nonnegative, on the
same reconstructed galaxy sources as stage 2 (razor-thin disks through the ring average of the shell
footprint, the bulge through its shell average, both by Gauss-Legendre quadrature, 128 nodes). Proposed
here; originality unverified. Part 2, conditional on the galaxies passing: the same 18 shells for the
clusters (1500-point grid, shell average by quadrature) and the lenses (through CR-2's interface, static
geometry, the well coupling, the three-dimensional bend route), and the joint solve with NL-1's
dimensionless balanced objective.

## Gates
K1 the shell ring kernel at d_0 = 0 reproduces the library's Gaussian ring kernel and its radial
   derivative to 1e-10 (control: at d_0 = 5 kpc they differ by more than 10%).
K2 quadrature: 128 against 256 nodes changes the ring-averaged kernel derivative by less than 1e-8
   relative, on the shell columns of three galaxies.
K3 the spherical shell average of the shell footprint at d_0 = 0 reproduces the library's shell kernel
   derivative to 1e-10.
K4 operator convergence: the shell columns change by less than 1e-3 under source-grid doubling (3000 to
   6000) on three galaxies.
C1 to C3 as stage 2: two solvers agree to 1e-6, optimality to 1e-8, convexity along every direction
   probed (the family is linear in the amplitudes, so the loss is convex on its domain).
R1 reproduction: with every shell amplitude forced to zero, the certified minimum equals stage 2's
   31.335 km/s to 1e-9.

## Decision rule (fixed now)
(a) The certified minimum on the 89 training galaxies with the consistent source is at most 21.9 km/s:
    part 2 runs, and universality is judged by the joint solve with every term at most one.
(b) It exceeds 21.9 km/s: the route is closed by a certified negative. No nonnegative linear combination
    of Gaussian footprints and shell footprints up to 50 kpc describes the galaxies with these sources.
    The minimum, its active columns and the tabulated-force sensitivity are recorded; part 2 does not run.
Validation and test at the training optimum; the slope against the integrated mass. Exposed data; the
three statuses stay separate.

## Evidence
`nk1.py` writes `nk1-results.json`; `nk1_checks.py` reruns the fast gates and anchors the archive,
registered before the stage-8 job. Failed gates are recorded, never re-thresholded; a corrective rerun
needs a declared amendment.
