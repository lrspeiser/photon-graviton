# The simple homogeneous growth interpretation has a supply mismatch

**Do not adopt the homogeneous temporal completion tested here.** If all the fitted optical source is massless companions, and its variation with distance is interpreted solely as time evolution with D=c*lookback time, then present photon conversion requires at least **2.99708e-8 J/m^3** of radiation. This is about **613,000 times** the comparison density obtained by adding the measured microwave blackbody to the high end of a broad historical UV-to-mm background range. This is a conditional benchmark mismatch, not a rigorous all-frequency exclusion or a rejection of every companion model. All six goals remain open.

## Extra assumptions being tested

The preceding optical-source inversion gave source requirements along a past light cone. It did not say that they are a time history. Here we explicitly add: spatial homogeneity, D=c*(t0-t), no volume expansion, no net energy-flow divergence, a frequency-independent fractional photon conversion rate alpha*c, and no source other than photon conversion. These assumptions are not implied by the brightness fit. The framework still lacks the full metric/transport completion needed to justify them.

## Formula and derivation

From the existing postulate z=exp(alpha*D)-1 and the added distance/time relation,

    dz/dt = -alpha*c*(1+z).

Apply the chain rule to the previously inferred companion density u_c(z). The assumed local conserved transfer law, with nonnegative capture sink Q, is

    du_c/dt = alpha*c*u_gamma - Q.

Consequently the necessary supply is

    u_gamma >= -(1+z)*du_c/dz.

The observer series of the retained optical law gives the exact conditional result

    u_gamma(0) >= [24*(q+1)^2/(12*q+13)] * u_c(0)
               = 6.357765095 * u_c(0)
               = 2.997079727e-8 J/m^3.

The chain rule, continuity equation and blackbody conversion used below are known mathematics/physics. The fitted optical response and its identification with a uniform companion source are the project's hypotheses; this formula is a derived consequence of those hypotheses, not a new universal law.

The corresponding required energy input is 7.24757e-26 W/m^3. Setting capture to zero makes the estimate permissive: positive capture would raise the photon requirement. A large pre-existing companion bath alone cannot supply this required current growth rate.

## Radiation comparison

The measured microwave temperature used here is **2.72548 +/- 0.00057 K**, as compiled and recalibrated by [Fixsen (2009)](https://arxiv.org/abs/0911.1955). The known blackbody relation gives u=a*T^4=4.17468e-14 J/m^3. Using this measurement requires no Big Bang origin or expansion assumption. We include these photons as potential donors even though the original motivation emphasized starlight.

For an integrated isotropic radiation intensity I, the known conversion is u=4*pi*I/c. [Hauser and Dwek's background review, summary](https://ned.ipac.caltech.edu/level5/March03/Dwek2/Dwek6.html) reports a historical UV-to-mm range of 45–170 nW/m^2/sr, with nominal value 100. This is an observationally informed scale comparison, not a modern complete-background upper confidence limit. Foreground subtraction uncertainties and incomplete frequency coverage preclude such a claim.

| Photon benchmark | Energy density (J/m^3) | Required / benchmark |
|---|---:|---:|
| Microwave plus 45 nW/m^2/sr UV-to-mm | 4.36330e-14 | 686,883 |
| Microwave plus 100 nW/m^2/sr UV-to-mm | 4.59385e-14 | 652,412 |
| Microwave plus 170 nW/m^2/sr UV-to-mm | 4.88727e-14 | 613,243 |

The required integrated intensity would be 7.15005e8 nW/m^2/sr. We do not compare that bolometric requirement with a single optical pivot-wavelength intensity, which would be an invalid bandwidth comparison. No additional fit or adjustment to the observation benchmarks was performed.

## Consequence for the prior 933-million-year scaling

That number remains a correct constant-bath equivalent duration for accumulated galaxy capture at the inferred observer density. It did not establish that photons produce or maintain that bath. The current test adds a temporal interpretation of the optical profile and finds an incompatible supply scale under the stated photon benchmarks. Increasing the universe's age changes an accumulated reservoir, not this present growth-rate requirement; no age or size bound is used here.

Possible changes must replace an identified assumption with equations: spatial transport can supply local energy through a nonzero flux divergence; optical focusing could have other sources; the distance-to-time mapping or gravitational coupling could differ. None is a demonstrated solution. A spatially incoming companion flux must still be traced to its sources and must preserve the optical and galaxy predictions. Merely choosing a larger external bath is insufficient to complete that calculation.

The retained empirical brightness and rotation fits are not altered. The source bridge remains conditional, and this particular homogeneous temporal interpretation is not supported. Reproduce with `python research_work/results/brightness-distance-consistency/optical-growth.py`; the analytic coefficient, constants and comparison values are saved in `optical-growth-results.json`.
