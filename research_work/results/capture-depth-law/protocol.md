# Compare common capture laws at matched optical depth

Before calculation, test kappa(r)=A*[2/(1+sqrt(1+r^2))]^p for p=1,2,3. These are explicit phenomenological alternatives, not microscopic derivations. Use the same isotropic absorption/in-place retention assumptions as isotropic-capture-profile. Normalize A independently by the declared central diameter optical depth tau_d=0.1,1,10, using a one-dimensional quadrature. This matches central attenuation when comparing different shapes; no stellar observations or desired halo densities enter A.

Test outer radii10 and30, to distinguish a robust inner shape from a chosen capture cutoff. For each of18 combinations use321 radial nodes with48 angular/ray quadrature nodes, then641 with96. Require volume absorption versus independent projected chords relative error below0.002 at fine resolution and central intensity exp(-tau_d/2) agreement below1e-6. Track all cases, not just flat-looking examples.

Derive enclosed deposit mass from q. Report v(9)/v(3), fitted log(v)-log(r) slope over3..9, effective density slope, and enclosed fractions at3,5,9. The asymptotic thin-capture prediction v proportional to r^(1-p/2) for p<3 is known spherical gravity algebra and only applies well outside a core and inside the boundary. p=2 need not be flat on the actual finite interval.

At impact5 calculate projected versus enclosed mass, using segmentwise transformed quadrature and8/16-node agreement to1e-9 relative. This ratio provides a lensing/force shape diagnostic independent of intensity normalization. Require coarse/fine v-ratio differences below0.003 and log-slope differences below0.003. Report radius10-to30 sensitivities; do not call either boundary a measured halo edge. No observational fit, energy-supply calibration, support proof or full nonspherical theory is claimed.

Follow-up after seeing p=3's nearly flat3..9 profile: for every radius30 case also report v(27)/v(9). This tests extension beyond the favorable interval without changing parameters. Retain the declining or rising result rather than presenting limited-range flatness as asymptotic flatness.
