# RW-1 completion reviewed during GF-1

The initial review at 2be9b33 correctly found no committed RW-1 result.
A later fetch brought ff02d8b (archive-anchor bookkeeping correction), then
e5367bd (archive and report). This update supersedes the initial missing-results
status. No feature-branch files were merged, edited or recalculated by GF-1.

[Committed RW-1 report](https://github.com/lrspeiser/photon-graviton/blob/e5367bd/research_work/results/path-memory/report-rw1.md).
The following are that report's results, not independent reruns:

- All reported numerical gates pass. That verifies evaluation of its declared
  radial force family, not the physical existence of a vortex.
- Per-galaxy strengths with p=1 and a small core fit whole curves at
  13.94 km/s, but a universal constant strength gives 31.49 km/s.
- U2 makes strength depend on baryonic mass approximately as M^-0.51.
  Galaxy training/validation/test errors are 21.32 / 29.68 / 17.28 km/s.
  Its cluster score is 153.23 per point versus the declared 23.28 limit.
- Clusters favor p=2 with strength about six times the ordinary-matter
  Newtonian response. U3 gives 20.64 per point but galaxy training error
  726.57 km/s. This is a force-rescaling fit, not an added dark-matter source.
- No tested member passes the combined lens requirement. Even each lens's
  individually best kernel gives total stellar-kinematic score 270.8 versus
  the limit 128 (sum of the six rounded rows). Exact Einstein-radius matches
  do not automatically predict the measured stellar motions.
- The cluster deflection curves are inferred from the hot-gas field and a
  stipulated static light coupling. There is no matched cluster shear
  catalogue in that test. They cannot establish observed cluster lensing.
- Formation, propagation, photon funding and the swirling transverse dynamics
  are not derived by RW-1. Its author states those limits explicitly.

## Attribution that must accompany the root-law observation

At large radius p=1 and lambda proportional to M^-1/2 give
a proportional to sqrt(M)/r. With a fitted acceleration scale this is the
well-known deep-MOND scaling a=sqrt(G M a0)/r, already stated by
[Milgrom (1983), primary paper](https://adsabs.harvard.edu/pdf/1983ApJ...270..365M).
RW-1 reaches approximately that dependence through a different fitted family;
it does not make the formula newly invented here. Its fitted exponent is near,
not exactly equal to, one half. Its acceleration scale is inferred from exposed
data, not predicted by an energy-funded companion microphysics.

GF-1 neither imports that law into its packet trajectories nor claims to have
replaced its missing derivation. It asks the more limited question whether
specified interactions can transmit turns and sustain prepared rotation.
