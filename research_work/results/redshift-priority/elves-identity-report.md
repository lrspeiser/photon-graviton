# ELVES-Field identity audit: three exclusions, 26 pending

The 29 staged candidates now have a first identity/group audit. No distance or velocity outcome was scored. Ten targets have PGC aliases confirmed through SIMBAD identifier sections, including three resolved by supplemental lookup. Eleven targets have at least one CF4 position within the fixed 30-arcsecond search radius; positional proximity alone is not accepted as identity proof.

## Exclusions from a fresh test

| Target | Evidence | Decision |
|---|---|---|
| UGC05797 | PGC 31693 appears in the earlier pilot and shares an excluded CF4 group | Exclude |
| dw1046p1244 | SIMBAD returns the staged name together with PGC 4689210 and older Leo dwarf names; CF4 places it in an excluded group | Exclude |
| NGC4592 | Its name occurs in two historical dust/energy input tables, companion_causal_test/data/dl14.dat and themis.dat | Conservatively exclude because of earlier source availability; this is not proof its redshift was previously scored |

The historical alias scan examined 799 tracked text files at commit 3e71904, before ELVES acquisition. It found two candidate names in old sources, with source hashes and matched aliases retained. The independently confirmed group exclusion for dw1046p1244 adds the third exclusion. Using a frozen historical checkout avoids counting the new audit itself as prior exposure.

Twenty-six candidates remain pending, not certified fresh. In particular, LV J1017+2922 and AGC740112 have positional candidates in already excluded groups but still lack confirmed matching PGC aliases in this pass. The failed exact-name and PGC lookups remain documented. Other unresolved names, alternate aliases, host/group membership, untracked or archived historical content and catalog selection still require audit.

UGC04115 resolves to both PGC 22277 and 22280 in the same SIMBAD identifier record. The CF4 positional candidate is PGC 22277. Both aliases are retained; no row is duplicated. A name resolver can return no PGC even when it recognizes a source, so lack of a PGC response is not an absence claim.

## Reproducibility and limits

elves-identity-audit.json preserves the first name queries, identifier sections, positional candidates and pre-audit source hashes. elves-supplemental-aliases.json preserves follow-up evidence. elves_historical_scan.py replays the conservative alias search; punctuation-insensitive substring hits require interpretation and are not automatic proof of observed outcomes. elves_identity_decisions.py verifies the frozen registry and CF4 hashes, checks supplemental literal-name evidence and reproduces all 29 decisions.

Only identifiers, coordinates and group membership enter this pass. The CF4 grouping remains a conservative exclusion convention, not a theory-independent gravitational environment measurement. No isolation flag becomes a path input, and none of the candidates has acquired a fitted velocity or individual stretch rate.

Next work must resolve the two positional group-overlap cases and remaining aliases, then audit velocity frame/convention, source selection and independent motion information. Preserve the final-model and uncertainty freeze before opening evaluation labels. The candidate acquisition has made progress, but the required fresh predictive improvement remains unestablished.
