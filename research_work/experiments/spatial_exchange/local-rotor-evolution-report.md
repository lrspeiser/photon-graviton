# SE-LR2: local circulation forms, spatial accuracy fails

20 September 2026. Protocol bd9c497, pinned executable59217ca. All six runs
finished; all individual ledger gates pass. The campaign fails its declared
spatial-circulation gate. The independent archive audit passes97 checks;
that validates the recorded failure, not the physical theory.

Starting with A=Pi=0 and a compact rotating internal excitation Q/P produces
a circulating A field. The excitation and its gradient energy are included
in the finite initial budget. Zero coupling and a nonrotating excitation
produce exactly zero A/Pi. Reversing epsilon reverses A/Pi exactly while
leaving Q/P unchanged. No external energy injection is present.

| Result at T1 | n24 dt.02 | n32 dt.01 |
|---|---:|---:|
| Initial total energy |0.000924120920|0.000941399928|
| Final A-sector energy |4.76060683e-7|5.07120052e-7|
| A energy / initial total |0.051515%|0.053869%|
| Circulation at radius1.5 |0.000236445584|0.000205664628|
| Maximum relative energy drift |4.26e-9|1.74e-10|
| Maximum relative angular-vector drift |5.85e-9|2.88e-9|

Initial energy differs because the same continuous compact profile is
discretized at different resolutions; it was not rescaled to fit results.
Halving dt at n24 changes radius1.5 circulation by3.85e-8 relative, passing
the0.1% gate. The n24/n32 spatial difference is14.97%, failing5%. Finer
resolution and a rotated-source test are needed before numerical accuracy
can support further interpretation. Energy near the box edge stays below
3.75e-13 of initial energy in these short runs; this does not validate an
outgoing boundary for longer runs or a strict numerical signal front.

Crucially, a circulating vector field is not automatically transported
angular momentum. Independently reconstructed final A-sector angular
momentum is only about6.36e-13 at n24 and3.95e-13 at n32 along the excitation
axis, compared with initial total angular momentum about4.44e-4. The other
components are negligible. This run does not demonstrate substantial outward
spin transfer, despite nonzero circulation. Q/P retains almost all of the
canonical angular momentum. This distinction matters for the proposed
following/swirl mechanism.

The internal excitation is prescribed inside the compact source, not derived
from ordinary matter. Matter production and fuel, scalar/conversion coupling,
light propagation, sustained circulation and stable stellar orbits are still
missing from this candidate. The calculation gives model units, not physical
galaxy or cluster predictions. No old gravity formula supplies an acceptance
criterion and no dark matter or expanding background is included.

The audit independently reconstructs endpoint energy/momentum/angular
momentum and loops from full field arrays. Intermediate extrema are checked
from recorded every-step traces; intermediate full fields are not archived.
Maximum independent endpoint discrepancy is3.25e-19. Full initial/final
arrays, traces, manifest hashes, per-run JSON and failed summary are preserved
in local-rotor-evolution-v1. Established numerical methods are credited in the
protocol; novelty of the coupling is not claimed.

Concurrent campaigns also advanced: SR-2 has8/9 runs archived and192 audit
checks pass; its final rotation run remains active. SE-R has4/7 new runs
archived and32 audit checks pass. Its mixed-channel n48/n56 source comparison
passes: wave-source change0.7585%, decreasing from1.1743%, and total-source
change0.002637%. These are source-derivative diagnostics, not lensing or
rotation-curve successes. Direct-companion refinement is still running.
