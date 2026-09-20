# OB-1: outgoing-layer contamination test

Declared20 September2026 before implementation/calculation. Supports requirement3 before any long-lived CC evolution. Does not change the running CC-2 archive or certify nonlinear3D boundaries by itself.

Test the exact vacuum scalar-wave restriction Fdot=Pi, Pidot=Laplacian F-gamma Pi with the same nearest-neighbor energy and damping ledger as CC-2. Use periodic interval[-6,6), N=256 and512, RK4 dt<=.2 dx, final time9. Initial compact right-moving pulse F=max(1-(x/.6)^2,0)^4, Pi=8x/.6^2 max(1-(x/.6)^2,0)^3. No matter/source or expansion. Compare at the same dx with a three-times-larger periodic domain[-18,18), with no damping and no boundary arrival by the endpoint.

Profiles gamma=gamma_max[max(0,(|x|-(6-width))/width)]^2, widths{1,2,3}, gamma_max{2,8,32,128}:24 layer runs and2 reference runs. The width1/gamma_max2 setting matches CC-2's vacuum outer layer. Test both reflection and wraparound by measuring the energy of the difference from the large-domain reference inside |x|<3.5 at t=9. Define relative state error=sqrt(E_difference/E_initial), using the same field-gradient and momentum energy. Gate<=.01. Also require scaled H+Q drift<=1e-5 and record total energy removed, maximum remaining amplitude and the full final profiles.

Choose no physical-model parameters. If multiple layers pass on both grids, identify the narrowest passing width, then smallest gamma_max; do not claim this selection proves performance for a different pulse spectrum,3D corners or nonlinear field speeds. Preserve every failure and numerical grid comparison. If none passes, keep the boundary requirement open and design a separate amended experiment before further calculation.
