# SE-LI1: interpolation improvement does not eliminate rotation discrepancy

20 September 2026. Protocolc91af04; executable145d8d8. All finest-grid
analytic controls pass. Four-point Lagrange weights exactly reproduce the
tested powers0..3 at33 fractional positions within floating-point precision.
The48 analytic field/grid/method/radius/orientation cases are preserved in
loop-interpolation-v1, along with hashes of all saved input arrays and metadata.

The known control field is A=.02 exp(-r^2/1.4^2)(axis cross x); its circulation
on a perpendicular radius-r loop is exactly2*pi*.02*r^2*exp(-r^2/1.4^2).
At n64 and512 loop points, the largest relative control error is0.544823%
for linear and0.005138% for cubic interpolation. Both meet their respective
declared2% and0.1% thresholds. These are analytic numerical fixtures, not
an old-gravity prediction used as a target.

Applying both methods to the unchanged LR4 final fields gives:

| Method,512 loop points | n48/n64 circulation difference | n48 rotation difference |
|---|---:|---:|
|Linear|3.59883%|2.52575%|
|Four-point cubic|2.35516%|2.26686%|

Both satisfy the5% spatial comparison and both exceed the2% rotation limit.
The published LR4 failure used256-point linear sampling and remains unchanged.
Higher interpolation order improves this measurement but cannot by itself
remove the discrepancy. The residual is sensitive to the evolved grid fields;
without an exact evolution solution this diagnostic cannot uniquely assign
each part of the error to interpolation versus evolution.

No simulation or physical coupling was changed. Cubic negative interpolation
weights are used only for this saved-field measurement, not substituted into
the positive matter/light Hamiltonian. Established Lagrange interpolation and
analytic quadrature controls are credited as numerical tools, not novel physics.

A separately declared SE-RR1 run now evolves the rotated case at n64 using
unchanged equations, compares against the existing unrotated n64 archive,
and requires both linear and cubic rotation discrepancies to fall below2%
and decrease from n48. Its result is pending; there is no new rotation pass.
Ordinary-matter excitation, sustained transport, complete light/matter gravity
and joint galaxy/cluster observations remain outstanding.
