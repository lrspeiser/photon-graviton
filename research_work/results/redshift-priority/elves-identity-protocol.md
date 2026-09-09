# ELVES-Field identity and group audit

Use all 29 provisionally eligible targets from the staged feature file. Freeze the historical registry at commit db87967. Resolve each supplied name with SIMBAD, retaining only identifier-section evidence, not distance or redshift information. A failed name lookup is unresolved, never evidence of freshness.

Independently search CF4 positions within 30 arcseconds of each target. This is a fixed conservative candidate-match radius, not a fitted tolerance or automatic identity proof. Keep all positional candidates and separations; multiple candidates remain ambiguous. Read CF4 object/group identifiers and J2000 positions only. Join resolved PGC aliases to groups and historical exposure, and compare normalized direct names against the structured galaxy-name registry. Record positional-only group hits separately from confirmed alias hits.

A named PGC match or direct historical galaxy name excludes that target from fresh evaluation. A resolved PGC's previously excluded group also excludes it under the group-independent protocol. Positional-only historical matches require identity follow-up and remain withheld. No-hit or failed lookup targets retain unresolved historical/host exposure; do not call them fresh. Selection never uses the target redshift, distance magnitude, residual, environmental fit or prediction success.
