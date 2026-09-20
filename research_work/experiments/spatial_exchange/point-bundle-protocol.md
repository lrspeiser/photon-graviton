# SE-PB: matched point-probe bundle refinement

Declared 20 September 2026. Keep SE-B's source equations, initial empty fields,
source energy, probe launch/detector geometry, nine-ray bundle and six test
bodies. Replace only probe interpolation by SE-PI's positive shrinking kernel.
There is no independently chosen probe radius. Matter source radius remains.9.

Seven T4,L16 runs: emission-only (theta0,chi0,mix0) and direct-companion
(theta pi/2,chi200,mix.5) at n32,48,64, all dt.01; direct-companion n32 at
dt.005. Evolve probe/background at the same RK stages as SE-B. Preserve the
original failed fixed-radius campaign, not overwrite it.

Retain SE-B individual ledger, boundary, cone, crossing and bundle-spacing gates.
For each emitter require n48/n64 central bend agreement within1%, with difference
smaller than n32/n48, and transport matrices agreeing within .001 absolute.
For the paired emitter difference (Y bend minus control bend), require the
n48/n64 change <=max(20% abs(n64 difference),1e-7 radians). This tests resolution
of either a positive, negative or absent difference, not a desired sign.
Time test at n32: bend<.1%, arrival offset difference<1e-5, matrix difference
<1e-4. All denominators use the finer result, floor1e-8 for bend. No thresholds
may be relaxed after seeing outcomes. Save full probe traces and endpoints.

This only addresses numerical interpretation of the candidate point Hamiltonian.
Cartesian interpolation artifacts still need rotation checks; physical source
radius, mature-field evolution, observer clocks, complete image geometry,
source budget and galaxy/cluster data remain unresolved. No observational
success or stable stellar orbit follows from passing these numerical checks.
