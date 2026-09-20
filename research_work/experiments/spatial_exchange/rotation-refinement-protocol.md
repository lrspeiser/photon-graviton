# SE-RR1: n64 rotated local transfer

Declared 20 September 2026 before the new evolution. LI1 analytic controls
pass but cubic measurement still gives a2.26686% n48 rotation discrepancy,
above2%. Preserve both the earlier linear failure and this follow-up result.

Run exactly the LR4 rotated initial excitation and unchanged combined RHS at
n64,L8,dt.01,T1. Common rotation remains angle.573 about normalized(1,2,3).
Reuse the archived unrotated n64 for comparison; do not recompute or modify
it. Archive initial/final fields and every-step LR4 canonical/outer ledgers.
Retain individual LR4 finite-state,energy<1e-5,angular<1%,edge<1e-5 gates.

Measure radius1,1.5,2 circulation with both linear and four-point Lagrange
interpolation at256/512 points, using the pinned LI1 functions. Each method's
radius1.5 rotation difference, denominator max(abs(unrotated),1e-10), must
be<2% and less than its LI1 n48 difference at512 samples. Require projected
field angular momentum difference<2% too. Require both interpolation methods,
not whichever is more favorable. Record all errors; do not change thresholds.

This is a new numerical refinement, not a replacement of the failed LR4
record. A pass would concern this T1 excitation and these measured quantities,
not every orientation, source size, stable orbit or lensing prediction. No
dark matter, expansion, new coupling fit or old-gravity target is introduced.
