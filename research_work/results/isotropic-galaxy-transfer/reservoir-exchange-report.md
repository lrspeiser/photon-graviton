# Local exchange between compact and extended companion reservoirs

Three shared exchange rules were executed with the one-third capture law and supplied inventory unchanged. Each initial shell is split into compact and extended states using its baryonic binding environment. The target halo is used only in fitting and scoring, never as a local driving field.

## Formula provenance and meaning

Known two-state rate balance: df/dt=k_on(1-f)-k_off*f, with fixed-environment stationary fraction f=k_on/(k_on+k_off). The proposed companion-specific ratios are:

| Rule | k_on/k_off | Interpretation |
|---|---|---|
| Binding | B(r)/v_ex^2 | Strong binding favors compact states |
| Released binding | [B(s*r)-B(r)]/v_ex^2 | Available inward binding-energy change favors the transition |
| Occupancy penalty | B(r)/{v_ex^2[1+rho_input/rho_star]} | Larger relative companion density suppresses compact entry |

B(r)=G M_star/(r+a_star) is the positive escape-binding scale of a known Hernquist stellar proxy. Stellar mass is from the ORIGINAL companion fit; a_star=Re/1.8153. The proxy is an approximation, not an independent measured three-dimensional matter field. The crowding factor is a proposed local penalty, not an established capacity law. We do not claim these ratios are unique new mathematics.

The compact shell fraction moves to s*r; the rest remains at r. Thus the cumulative profile is F_final(r)=F_initial(r)-F_selected(r)+F_selected(r/s). This is known conservative transport. s still specifies a phenomenological movement distance; these tests do not independently derive wave support or a stopping radius. The states use the initial environment: no updating of occupation, self-gravity or detailed balance at the new destination has been solved. The arbitrary overall rate cancels; no relaxation time or cosmic age is inferred.

## Shared and omitted-galaxy comparison

Two parameters per rule are shared across six targets; separate five-target fits predict the omitted sixth. All targets have been inspected previously, so this is conditional cross-validation, not blind observational confirmation. RMS values below are cumulative-profile percentage points, not lensing/rotation errors. Standard geometry remains a benchmark only.

| Geometry | Input | Exchange | Shared fit | Omitted transfer | Earlier partial transfer |
|---|---|---|---:|---:|---:|
| companion_regular | distant | binding | 17.06 | 20.34 | 20.02 |
| companion_regular | distant | released_binding | 16.99 | 20.15 | 20.02 |
| companion_regular | distant | occupancy_penalty | 18.23 | 21.85 | 20.02 |
| companion_regular | near_1.1R | binding | 17.75 | 20.67 | 20.51 |
| companion_regular | near_1.1R | released_binding | 17.60 | 20.96 | 20.51 |
| companion_regular | near_1.1R | occupancy_penalty | 19.24 | 23.08 | 20.51 |
| standard_flat_FLRW | distant | binding | 20.21 | 24.16 | 24.32 |
| standard_flat_FLRW | distant | released_binding | 20.32 | 24.14 | 24.32 |
| standard_flat_FLRW | distant | occupancy_penalty | 19.69 | 23.47 | 24.32 |
| standard_flat_FLRW | near_1.1R | binding | 20.12 | 24.47 | 24.33 |
| standard_flat_FLRW | near_1.1R | released_binding | 20.06 | 24.07 | 24.33 |
| standard_flat_FLRW | near_1.1R | occupancy_penalty | 19.83 | 23.56 | 24.33 |

All population cases, compact fractions, shared/omitted parameters, boundary flags and individual curves are archived. No target-specific compact fraction was fitted directly. Source fields and mass proxies remain conditional inputs; normalized placement is not evidence for an adequate absolute photon supply or a successful projected lensing signal.

## Outcome

Under retained geometry and distant input, the released-binding rule nearly matches the earlier partial-migration prediction: 20.15 versus 20.01 cumulative RMS points. Binding-only gives 20.34 and the occupancy penalty gives 21.85. None improves the primary retained-geometry benchmark. These differences are not significance estimates, and grid spacing and prior target exposure preclude treating small ranking differences as discoveries.

The shared released-binding fit has s=0.1391 and v_ex=720.84 km/s. Its compact fractions range from about 30.9% to 37.6%, calculated from the initial stellar binding environment rather than fitted separately per target. This is a candidate local rationale for the earlier roughly 36% mobile fraction, not a first-principles derivation. Both models still choose a movement scale phenomenologically.

It retains the same tension: concentrated targets improve, but J1112, J1621 and J1630 worsen relative to no migration. For example, J1112 maximum cumulative error rises from 5.98 to 41.73 points in released-binding omitted prediction, while J0037 falls from 62.25 to 36.44. A lower aggregate error does not mean all galaxies are reproduced.

Under standard comparison geometry, occupancy suppression gives 23.47 points for distant input versus 23.55 without migration and 24.32 for earlier partial migration; the near-input result is 23.56. This small aggregate change is not a robust joint solution or a new lensing result. Standard geometry is not adopted as a premise for the hypothetical model.

The useful outcome is a local exchange candidate with performance close to the best phenomenological split, using the same two shared parameters. Its missing physics is the simultaneous evolution of exchange, inward transport, occupation and self-gravity. These frozen initial-field rates do not yet establish that final reservoirs are in equilibrium or stable.

## Verification and physical limits

The compact and extended inventories remain nonnegative and sum to the original input. CDFs are monotone, inward ordering holds, and s=1 reproduces the captured profile. Doubling source cells changes cumulative fractions by at most 4.075e-06. Parameter bounds were declared in reservoir-exchange-protocol.md and not expanded after results.

Mass-equivalent inventory conservation is not a full formation-energy ledger. The released gravitational energy requires a receiving channel, and rate ratios alone do not supply the work, momentum transfer or support needed after migration. The next physical completion would evolve exchange, transport, occupations and the potential together; no first-principles conversion or self-consistent exchange equilibrium has been derived here. No new stellar-motion or lensing fit is claimed.
