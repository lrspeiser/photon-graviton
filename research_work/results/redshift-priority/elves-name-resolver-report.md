# Named positions reveal a source-association ambiguity

All 24 current pending names were queried through CDS Sesame. There were no request errors; 13 returned resolver records and 11 returned no match. Six resolutions came through SIMBAD and seven through NED. A successful response is not automatically an independent confirmation: a SIMBAD-backed resolver may reproduce the same SIMBAD position already queried.

## Specific findings that change the next action

**LSBC D565-09 requires source-association review.** NED resolves the requested name to a position 15.193 arcseconds from the publisher position. That named position is 0.161 arcseconds from SIMBAD LEDA2806909, but 15.036 arcseconds from SDSS J093313.66+203055.3. The earlier metadata shows no current basic velocity reference for the LEDA object, while the SDSS object has a spectroscopic reference. The publisher-position cone contained both. This is an unresolved object association, not proof that the ELVES paper is wrong and not permission to substitute whichever measurement fits better.

**LSBC D634-03 remains ambiguous at the object-type level.** NED's named galaxy position is 0.395 arcseconds from a SIMBAD group entry, 2.417 arcseconds from an entry typed Gl?, and 3.755 arcseconds from LEDA2806961. A nearest-entry rule would select the group, which cannot automatically stand in for an individual galaxy's spectrum. The earlier radio match is additional source evidence, not a resolution of these different catalog representations.

**Large centroid offsets do not alone disprove established aliases.** DDO047's resolved position is 13.215 arcseconds from the publisher position; UGC04115's is 29.410 arcseconds away. Both resolve to the same SIMBAD records supported by explicit aliases in the prior audit. The batch does not explain those centroid differences. It demonstrates that a much tighter blanket angular cut, selected after seeing these offsets, would discard supported associations. The original fixed 30-arcsecond search is retained.

KKH86's NED position is 3.838 arcseconds from the prior LEDA2807150 counterpart. It provides positional corroboration but no explicit PGC alias in this response. No PGC alias, host assignment or fresh status is inferred from that separation alone.

## What was accessed

The exact request URLs, backing resolvers, allowed metadata, response hashes and returned XML tag names are saved in elves-name-resolver.json. The full HTTP responses were parsed in memory, and some included Vel or z elements. Their values were neither projected into the saved output nor printed, fitted or scored. Therefore this is a filtered metadata inspection, not a claim that the remote server transmitted only metadata or that outcome bytes were inaccessible to the code. The projection allows only oname, alias, otype, jradeg and jdedeg. The records returned no additional alias elements.

The individual name/position sources can be retrieved through [CDS Sesame for LSBC D565-09](https://cds.unistra.fr/cgi-bin/nph-sesame/-oxp/SNV?LSBC%20D565-09) and [LSBC D634-03](https://cds.unistra.fr/cgi-bin/nph-sesame/-oxp/SNV?LSBC%20D634-03); those live endpoints may return more fields than the saved projection. Do not display full responses during a pending-target audit. This batch used no broad web searches or target redshift-distance comparisons.

**Formula provenance:** angular separations use the established spherical great-circle relation, not a new model formula. They compare positions only. No time rate, photon-conversion parameter or individual velocity was fitted.

## Reproduction and scope

The acquisition command is `python research_work/results/redshift-priority/elves_name_resolver.py`; it requests the pinned 24-target stage-2 population. Inputs and protocol hashes are saved. Network results can change, so the saved JSON is the audit snapshot. Its input population and allowed metadata fields were verified offline, and its angular comparisons were independently checked against Astropy's sky-coordinate separation calculation. See elves-name-resolver-verification.json for the checked scope.

No new exclusions are made here. Current eligibility remains five excluded, 24 pending and zero certified fresh. Review the original source association for the two LSBC entries; continue bounded provenance checks for the others. The findings justify measurement-level caution, not rejection of the physical idea or a claim that measurement ambiguities explain the 164-group residuals. The common-model freeze, independent nuisance uncertainty and fresh predictive improvement remain unfinished.
