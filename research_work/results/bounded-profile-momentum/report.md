# Momentum and recoil of the conversion-region support

**The fixed conversion profile receives the electromagnetic wave's lost canonical momentum. Allowing its support to move closes both the Hamiltonian energy and canonical momentum ledgers in the tested model.** Some lost wave energy then becomes support recoil, rather than all entering the traveling reservoir.

This is a formal dynamical completion of the prescribed-profile approximation. It still does not provide ordinary massless-companion momentum, a symmetric relativistic stress tensor, a gravity-derived profile or observational validation.

## Momentum of the fixed profile

The preceding candidate uses `h=v*e+P*(1-v*s_x)`, with s=t and a spatial coupling W(x). Its canonical translation-generator density is `pi=D*B-P*s_x`. On this branch pi=e for a right-moving wave; the reservoir term vanishes.

Known Hamiltonian translation accounting gives

`d p_wave/dt = -integral v_x*e dx`,

`d I_profile/dt = +integral v_x*e dx`.

The spatial derivative here is the explicit derivative through the profile, holding its clock coordinate fixed. The profile impulse is the reaction on the system enforcing that spatial coupling. These are known canonical identities applied to the postulated model, not a new capture law or a measurement of a physical force.

We integrate the profile force as an additional ODE alongside the characteristic wave/reservoir equations. A separate time quadrature of sampled force verifies its integrated value. Thus the result is not merely a definition of impulse as the missing final momentum.

For the dynamic example with a=0.1, initial wave energy/momentum 0.375 and ordinary endpoints:

| Quantity | Result |
|---|---:|
| Final electromagnetic momentum | 0.3227654912 |
| Net profile impulse | 0.05223450884 |
| Reservoir energy gained | 0.05223450884 |
| Maximum canonical momentum-balance error | 4.45e-16 |
| Independent integrated-force discrepancy | 8.71e-13 |

Units have c=1 and are dimensionless. Equality between the final gained reservoir energy and profile impulse follows here from the matched endpoints and the model's energy/momentum allocations. It does not mean the reservoir is carrying that momentum. The static a=0, b=0.4 control has nonzero transient profile impulse but returns to zero net impulse and produces no retained energy gain.

If the acquired reservoir energy were simply relabeled as ordinary forward massless companions with momentum E/c, their momentum would be added on top of the profile reaction already required by the equations. The interaction and momentum ledger would need to change. The previous missing stress issue is therefore substantive, not solved by renaming P.

## Why a finite support needs an energy account

The fixed-profile approximation treats its support as an immovable momentum receiver. For a freely recoiling support of unchanged rest energy M in c=1 units, known relativistic mechanics gives

`K_recoil = sqrt(M^2 + I_profile^2) - M`.

This is positive at any finite M and nonzero net impulse. If the outgoing wave loss were kept unchanged, some of that loss would have to pay the recoil rather than all becoming reservoir energy. An external anchor, changing support rest energy or another field could alter the account, but cannot be omitted.

The fixed-impulse comparisons in `results.json` are diagnostic partitions, not self-consistent moving-support solutions. A sufficiently heavy support makes recoil energy small; the fixed approximation is not rejected merely because it is an approximation.

## Executed moving-support completion

We then make the profile center X a canonical variable with momentum K, replace W(x) by W(x-X), and add support Hamiltonian `sqrt(M^2+K^2)` to the wave/reservoir Hamiltonian. Its equations are

`Xdot=K/sqrt(M^2+K^2)`,

`Kdot=integral v_x*e dx`.

The wave feels the moving profile throughout the integration. Energy exchanged with that motion appears automatically in the Hamiltonian equations. A rigid shape in these coordinates is still postulated; adding a relativistic free-support kinetic term does not make the interaction Lorentz invariant or derive a material capable of maintaining that shape.

The constant support rest-energy baseline is included in the total Hamiltonian. The following fractions subtract only that unchanged baseline and are normalized to the initial wave energy:

| Support rest energy / initial wave energy | Outgoing wave fraction | Traveling reservoir fraction | Support recoil fraction |
|---:|---:|---:|---:|
| 10 | 0.8629132352 | 0.1361471699 | 0.0009395949 |
| 100 | 0.8609386692 | 0.1389646406 | 0.0000966902 |
| 10,000 | 0.8607102951 | 0.1392887349 | 0.0000009701 |

Each row sums to one within numerical error. Maximum fractional energy error is below 1.9e-11 for these runs, and maximum fractional canonical momentum error is below 1.4e-15. A tighter repeat of the support-ratio-10 case agrees in the reported energy fractions to better than 1e-8.

The light output changes when the support moves; simply subtracting a recoil number from the fixed-profile energy split would not predict this change. As support mass increases, the energy fractions approach the fixed-profile values. The support mass ratios are illustrative, not fitted properties of a void or gravity well.

## Consequences and limits

The demonstrated canonical completion is preferable to silently fixing the profile while asserting an isolated momentum-conserving system. It explicitly identifies a receiver of momentum and an associated energy cost. It does not identify a physical cosmic medium or justify that receiver's existence.

The traveling reservoir still lacks the ordinary relation between its energy flux and physical relativistic momentum in the homogeneous clock branch. The support closes global canonical conservation, but it does not supply a local symmetric stress-energy tensor for the reservoir or make it a graviton. No claim about lensing follows from these checks.

The creation and maintenance of W, clock preparation, matter standards and forces, arbitrary three-dimensional propagation, quantum stability and capture into a gravitational source remain unresolved. The previous bounded endpoint redshift result remains conditional; its observational rate and source profile have not been fitted or derived.

This calculation does not open holdouts, assume a cosmic age, evaluate the deferred total photon supply, or establish a strong case for the full theory. It resolves one explicit accounting omission inside a hypothetical model while preserving the more fundamental stress/identity limitations.

## Reproduction

Run `run.py` for fixed-profile impulse and recoil diagnostics, then `moving.py` for the coupled moving-support calculations. Both result files retain input hashes, controls and refinement evidence. The earlier bounded-wave code is unchanged.
