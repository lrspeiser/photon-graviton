# One conditional response for motion and lensing

## Result

The finite M87 deposit profile now has a common prediction for circular acceleration and asymptotic light deflection under the previously considered equal-potential, ordinary mass-energy response. It does not yet have an observational match. The new result is a shape constraint: one common rescaling of this reservoir cannot independently change its motion and lensing predictions.

| Radius (kpc) | Deflection / (deposit circular-speed-squared / c^2) |
|---|---:|
| 1 | 91.465 |
| 10 | 20.570 |
| 30 | 12.795 |
| 100 | 8.683 |
| 300 | 6.417 |
| 500 | 5.570 |
| 1000 | 4.000 |

Deflection here is in radians. These dimensionless ratios concern the deposited component alone, not the total gravity including ordinary matter, the central black hole and other emitters. They are conditional forward calculations, not measured ratios or novel fundamental laws.

## Explicit assumptions and derivation

Use the [resolved spherical power profile](projection-report.md). For this comparison assume a stationary, supported reservoir with negligible stress corrections and ordinary weak-field gravitational response: Phi=Psi and mass density rho=u/c^2. Its support and production are not derived by making this assumption. Freely streaming massless companions cannot automatically be treated as such a stationary nonrelativistic reservoir.

For a common exposure normalization T, take u(r)=T q(r). This preserves the already calculated shape. It is a restricted steady-history comparison, not a claim that all past emission follows today's profile. Spatially varying histories, finite startup, escape or deposit motion generally change the shape and must be calculated explicitly.

Let P3(<r) be deposited power inside a sphere, and P2(<b) be deposited power inside a projected circular aperture, both from the preceding calculation. Established spherical gravity and weak-field lensing give

    M3(<r)=T P3(<r)/c^2
    M2(<b)=T P2(<b)/c^2
    v_deposit^2(r)=G M3(<r)/r
    bending(b)=4G M2(<b)/(c^2 b).

Therefore, at b=r,

    bending(r) / [v_deposit^2(r)/c^2] = 4 P2(<r)/P3(<r).

T cancels. A hypothetical constant multiplier affecting both potentials equally also cancels, but this algebra does not establish a viable enhanced-coupling theory. Changing only the force on matter need not enhance light bending; see [the existing scalar-response analysis](../gravity-response/motion-and-lensing.md).

The known circular-lens deflection law and the distinction between projected and spatial mass are reviewed in [Kochanek's strong-lensing lectures](https://ned.ipac.caltech.edu/level5/March04/Kochanek2/Kochanek3_2.html). The metric conventions follow [our previously checked equal-potential kernel](../deposit-lensing-kernel/report.md). We use physical asymptotic deflection, without a cosmological distance-redshift conversion. Predicting observed image angles additionally requires source/lens/observer geometry and finite-distance effects; deflection is not directly an image separation.

These are standard response equations applied to our conditional deposit shape. No first-principles companion action, gravitational slip law, or time-conversion mechanism is newly derived here.

## Independent check

The script independently integrates transverse gravitational acceleration along the complete unperturbed ray. Its field extends outside the 1 Mpc deposition receiver; truncating the field at the receiver edge would be incorrect. With z=b tan(theta), the projected-equivalent enclosed power from the field integral is

    integral_0^(pi/2) P3(<b sec(theta)) cos(theta) dtheta.

Outside the receiver P3 is constant. This field integral agrees with the independently computed aperture power at all ten tested radii, with largest relative difference 9.46e-9. The SI-unit motion and deflection outputs also recover the dimensionless ratio. This verifies consistency of the conditional response, not its agreement with observations.

## Amplitude and observational implications

response.json reports coefficients per billion years strictly as a unit normalization, not an assumed age or permitted unlimited source lifetime. M87's current conditional power corresponds to 2,276.91 solar masses of energy per billion years if unchanged and retained. At 100 kpc, the deposited component's v^2 coefficient is 3.967e-6 (km/s)^2 per billion years; its bending coefficient is 7.905e-11 arcsec per billion years. These small coefficients flag the need for a complete source history and response calculation, but do not establish a total cluster energy shortage. M87 alone is only part of the previously cataloged supply.

The table identifies a test a larger normalization cannot evade: once a common deposit shape and response are specified, motion and bending must be compatible at every tested radius. Conversely, this ratio is not universal to every companion model: changing the spatial history, gravitational potentials or reservoir dynamics requires a new prediction.

Observed M87 stellar velocities are not circular-orbit speeds for individual stars. A comparison needs orbit/distribution modeling, ordinary matter and other source contributions, plus actual lensing measurements and their geometry. No observation was fitted, no holdout was opened and no objective is complete. All six objectives remain open.

Run `python research_work/results/cluster-m87-finite-emitter/response.py` after run.py and project.py.
