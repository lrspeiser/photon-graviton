# SE-O: candidate observer convention

All 48 local fixtures and two uniform controls pass.

Adopt the massive test-body action as a clock convention: d tau/dt = partial H/partial m.
Solve the same Hamiltonian for the observer momentum at its prescribed velocity.
Photon energy in that convention is (H_gamma - p_gamma dot v_observer)/(d tau/dt).

This is an explicit modeling assumption, not a derived material-clock or detector
mechanism. The homogeneous clock formula is derived from the candidate Hamiltonian;
no old-gravity formula supplied a pass target. Hamiltonian homogeneity is established math.

| Check | Maximum error |
|---|---:|
| velocity_residual | 1.69967494e-17 |
| mass_derivative_error | 2.32909247e-11 |
| euler_error | 3.33066907e-16 |
| scaling_error | 0 |

The first execution failed at JSON serialization of a NumPy boolean; its failure
record is preserved in observer-controls-v1. The corrected writer produced
observer-controls-v2 with unchanged equations, fixtures and gates.

Current bundle arrival offsets remain coordinate quantities. Their archives lack
the full detector-field history needed to integrate this clock, so no observed
delay or frequency shift has been retroactively inferred. Future runs must
integrate the observer clock alongside the evolving field. Full light-quality,
local-clock and observational validation remain outstanding.
