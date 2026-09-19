# CWC-1 2D diagnostics declared before execution
The 2D TE electromagnetic and receiving field use the same reciprocal engine.
The scalar profile is generated from zero; full final arrays are saved.
Frozen rays integrate Hamilton's equations with x as the monotonic independent
coordinate, using a cubic spline of the generated field and its derivatives.
Ordinary tolerances2e-9/2e-11 with max spatial step0.1 are compared against
2e-11/2e-13 and step0.025. Independent straight-path integration uses641 points.
The six signed impact parameters and all raw ray paths are retained.
A ray step test checks ray integration only, not spatial field convergence.

The periodic finite box is part of this fixture: all energy, including energy
crossing a periodic edge, remains counted. No isolated-boundary claim is made.
Material source strengths differ by design and can break exact mirror symmetry.
The declared1% mirror criterion remains unchanged. Same-position virtual
material probes with action/mass differing by two separately test universality.
No test-particle energy is added to the field: these are linear response probes.

The weak-gradient gate is explicitly conditional on max|g phi|<0.02; when
inapplicable its comparison is reported, not interpreted as a Born-accuracy pass.
Neither these rays nor the force proxy supply a relativistic lensing metric.
