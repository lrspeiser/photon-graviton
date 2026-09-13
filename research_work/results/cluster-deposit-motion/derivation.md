# What must a deposited gravitational reservoir do during a collision?

Our current capture rate does not specify the motion of the deposited source.
This is an incomplete equation system, not evidence that the source cannot exist.
An angular benchmark below identifies an observation the missing equation must
eventually predict. No collision simulation or successful cluster fit is claimed.

## Formula provenance and the missing equation

The following is established continuum conservation mathematics, applied
conditionally to the proposed deposited-energy reservoir. Let u_d be its energy
density, F_d its energy flux, and q its net local energy supply, including capture,
release and any work exchanged with other sectors in the chosen energy convention:

    partial_t u_d + divergence(F_d) = q.

This balance law does not determine F_d. In particular, setting decay to zero
does not set the velocity or stress of stored energy. A companion can cease to
propagate at light speed after capture and still move with its receiving system.

For a finite reservoir with vanishing surface flux, define

    E = integral u_d dV,     X = integral x u_d dV / E.

Integration by parts gives the exact first-moment identity

    dE/dt = integral q dV,
    E dX/dt = integral F_d dV + integral (x-X) q dV.

Thus permanent storage with no flux and no new supply leaves X fixed in the
chosen spatial frame. If energy is carried with velocity v, the effective
advective postulate F_d = u_d v gives a reservoir-weighted mean velocity plus
the source term. Choosing v still requires dynamics; calling it the velocity
of the gravity well would otherwise assume the answer we want to predict.

For a moving control boundary, its surface transport must be retained. In a
relativistic field theory, u_d and F_d must come from a consistent stress-energy
tensor and clock convention; this Newtonian moment identity is not that closure.

## Three distinct physical completions

1. **Energy stored in ordinary receivers.** Specify receiver abundances in stars
   and gas, their capture cross sections and retention. Stored energy follows
   those receivers. A universal capture rate per unit ordinary mass does not
   automatically distinguish gas from stellar systems. Receiver abundance and
   incident radiation history must predict the division before a lens fit.
2. **An independently moving bound reservoir.** Supply a phase-space evolution
   equation or field stress law, including its gravitational response and
   interactions with gas. Collisionless behavior is an effective hypothesis,
   not a consequence of energy conservation. Its origin in photons still needs
   calculation; matching the behavior of an unseen mass component alone does
   not identify that origin.
3. **Stored spacetime deformation.** Supply a dynamical field equation that
   transports its energy and momentum as wells move. A static potential solution
   at each time does not establish this transport. Prescribing a translating
   deformation implicitly supplies flux and possibly external work.

These are optional effective postulates, not newly discovered fundamental laws.
The existing capture-storage receiver model concerns lifetime and capacity, not
this collision dynamics. None of these three completions is adopted by this audit.

## Observational benchmark without a cosmological distance

The published Bullet subcluster BCG and plasma-aperture positions (Table 2), and
the lensing-peak offset (section 3), supply an angular target. The peak is near
the galaxies rather than the plasma center. These are published reconstructed
positions, not raw shear data or a direct image of deposited energy.
[Clowe et al. (2006)](https://arxiv.org/abs/astro-ph/0608407).

`inputs.json` transcribes these summaries and excludes distance, mass, merger
age and speed. `run.py` computes their separations on the celestial sphere and
cross-checks them with independent Cartesian unit-vector geometry. Its segment
projection is only a description of the positions, not a gas-bound deposit
fraction: a lensing peak is not an energy centroid, and superposed profiles can
shift a peak nonlinearly.

The calculated angular separations are:

| Pair | Separation in arcseconds |
|---|---:|
| Brightest subcluster galaxy to plasma aperture center | 43.98 |
| Brightest subcluster galaxy to reconstructed lensing peak | 9.63 |
| Plasma aperture center to reconstructed lensing peak | 36.61 |

The lens peak projects 0.177 of the way along the galaxy-to-gas segment but is
also 5.63 arcseconds off that line. Neither number determines a mixing fraction
or a dynamical drag coefficient. Both use rounded published positions.

The original lens reconstruction assumptions remain relevant. A modified law
relating the field to lensing ultimately requires forward prediction of galaxy
image distortions. We have no centroid covariance here and assign no rejection
significance to the distances. This famous, already-exposed system cannot be
the reserved novel prediction of objective 6.

## Next decisive calculation

Choose or derive a receiver/field transport closure with its momentum and energy
ledger. Evolve the same pre-collision deposited profiles through a collision
without resetting them to the measured final lensing map. Predict projected
stellar light, gas and lensing jointly, propagating geometry and reconstruction
uncertainty. The pre-collision profiles must also satisfy the galaxy supply and
motion tests; they cannot be independent arbitrary cluster halos.

The necessary data are a usable shear catalog with source-redshift information,
optical member positions/light, gas imaging/spectra, instrument responses and
selection/uncertainty models. Published coordinates are a benchmark, not a
replacement for those data. Full cluster validation remains open.
