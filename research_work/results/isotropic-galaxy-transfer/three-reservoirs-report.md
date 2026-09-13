# Sun calibration does not transfer under the minimal reservoir rules

Three alternatives are now expressed through common conversion, storage and gravity equations in the [branch specification](../../../research_plan/sun-earth-moon-three-branches.md). The calculated minimal versions do not explain the full Sun, Earth and Moon gravitational masses with one source-to-reservoir normalization. This does not rule out all retention laws, external supply or independently self-bound fields. No new per-body fitting was performed.

## Shared fit and predictions

Use the retained alpha and a fixed source-to-capture path length D=one body radius. With maximal capture, lossless storage gives E=T*P_c; stationary renewal gives E=tau*P_c. More generally only K=eta*T or eta*tau is identified. Fit K to the Sun and freeze it for Earth and Moon.

| Body | Converted power on declared path (W) | Predicted/actual full gravitational mass using solar K | Required effective K if fitted individually (years) |
|---|---:|---:|---:|
| Sun | 2.14816e9 | 1, by calibration | 2.63617e30 |
| Earth | 8.94902e-3 | 1.38752e-6 | 1.89992e36 |
| Moon | 1.80924e-4 | 2.29409e-6 | 1.14912e36 |

The final column diagnoses the mismatch; those individual values are not adopted as fits or actual ages. A shared gravitational coupling multiplier would rescale K and leave these relative failures unchanged. The target here is deliberately all observed gravitational mass, not an independently measured small excess. Known mass-energy source relations do not prove that the source originated as photons.

As an even more permissive energy ceiling, setting all conversion fractions to one still gives Earth/Moon predicted-to-target ratios 0.00015135/0.00091856 under solar calibration. This ceiling would not preserve observed outgoing light; it is not a proposed physical success. Both comparisons exclude choosing one duration or efficiency product to fit all three under the declared source proxies.

The source proxy for Earth and Moon is their geometrically intercepted solar power, representing reradiated-plus-reflected energy approximately. It is not their intrinsic nuclear luminosity, and it is not counted as new energy on top of solar generation. Source histories could differ, but need to be derived rather than freely assigned per body. No universe age or size is assumed.

## Independent traveling versus independently bound

For an unconfined relativistic population, an upper residence scale is D/c. Even giving every generated companion that full residence time, the fitted-alpha source supports mass-equivalent inventories of at most 5.55e-8 kg at the Sun, 2.12e-21 kg at Earth and 1.17e-23 kg at the Moon in the declared zones. With full photon conversion the corresponding loose ceilings are 9.88e9, 0.0412 and 0.000832 kg. These are transport energy inventories, not a complete radiation-pressure gravitational solution. They are insufficient to identify with the observed source masses under the stipulated E/c^2 coupling.

A genuinely self-bound independent concentration is a different branch. It cannot be assessed by substituting an arbitrary long residence time: the model must supply confinement, stresses, stability and motion. That branch remains open. No physical removal experiment or orbital persistence has been simulated.

## What the equations distinguish

- A, lossless body-bound storage: accumulation must predict the required history and capture. No continuous decay power is assumed.
- B, independent concentration: free escape fails the declared full-gravity inventory; self-binding remains unspecified.
- C, continual renewal: common steady lifetime/efficiency fails the same static mass scaling as A. Turning off input predicts exponential decay only if the proposed relaxation law holds.

Motion is represented by a transport/continuity equation with a specified carrier velocity; it does not alone create a dissipative power requirement. Graviton numbers remain proportional to 1/nu_g for fixed power, and no frequency has been chosen or inferred.

These results use an exposed mass/radius consistency comparison, not a blind test. They do not displace the existing ordinary mass baseline or validate replacing it. The capture zone extends beyond the surface, so no surface acceleration is inferred without its radial distribution.

The earlier 0.373-W Earth result used thermal IR traveling all the way to lunar orbital distance. This calculation uses a shorter one-radius path and a different, explicitly generous total outgoing-power proxy. The numbers answer different path/source questions and are not contradictory.

Reproduce with `python research_work/results/isotropic-galaxy-transfer/three-reservoirs.py`. Solar normalization is verified exactly; complete input powers, target energies, durations and transfer ratios are saved in `three-reservoirs-results.json`. All six broader goals remain open.
