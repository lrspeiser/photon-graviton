# Remaining-target name resolution

Use exactly the 24 pending names in the stage-2 quarantine screen. Query CDS Sesame using its SIMBAD/NED/VizieR resolver chain and return only resolver identity, object name/type, position and aliases if supplied. Do not display or save numerical velocity, redshift, distance or photometry fields. Raw response bytes may contain additional fields, so record all returned tag names and a response hash while saving only the allowed metadata. Do not print raw responses or use broad web searches for these names.

Retain every request failure and no-match response. Compare returned positions against the publisher feature position and the earlier SIMBAD cone counterparts. Report separations using the established spherical formula; do not use a prediction residual to pick a match. A name resolved to a position is corroboration, not proof that every nearby database object is identical or that an historical group overlap is absent. Keep source classifications and conflicting positions visible.

The trial query for LV J0913+1937 returned only position/name/type metadata and established the XML schema before this batch; it displayed no outcomes. No new eligibility decision follows automatically from resolver success.
