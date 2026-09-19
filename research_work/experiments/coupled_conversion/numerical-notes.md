# CWC-1 implementation details declared before spatial runs
The reciprocal missing-source negative control begins with a counted scalar
wave as well as EM; it deliberately omits the EM-to-scalar force. This is a
bad-equation control, not a source-funded production run.
Homogeneous checks use DOP853 with rtol2e-11 and independent finite launch-time
differences3e-4/1.5e-4. Directional Hamiltonian derivatives use central finite
differences1e-6 at random non-boundary material positions; reversal checks
compose40 steps forward and40 backward.

Spatial material interpolation is multilinear, and both deposition and
material recoil are its exact discrete variations. Source populations at
different grid levels have the same physical positions and action/mass.
No force softening or energy floor is introduced. Fixed periodic box energy
is counted in full, not mislabeled as escaped energy.

Receiver diagnostics sample the right-going characteristic amplitude every
four steps. Hilbert phase supplies carrier frequency, normalized by the
predicted receiver clock; power-weighted centroids use integrated clock time.
Fixed coordinate windows12..18 and18..25.5 isolate the prepared two pulses;
retain the1%-of-peak envelope mask. A failed identification/window is not
converted into a successful timing measurement. Uncoupled same-grid,
same-carrier controls remove the leading discrete propagation offset.
This is a finite pulse diagnostic, not a high-resolution astronomical line fit.
No particle-number claim is extracted from a classical amplitude decrease.

The full stage is allowed to finish with failed gates. Thresholds stay fixed;
any corrective refinement must be declared and run into a new directory.

Right/left energy proxies use the central derivative in a(D minus/plus grad A)^2/4;
their sum is a diagnostic and differs from the edge Hamiltonian at finite grid.
Momentum residual normalization uses initial EM energy, falling back to the
initial receiving or material energy only for controls with zero EM.
