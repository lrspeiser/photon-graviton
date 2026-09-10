# Tonry source-membership audit

Use the author's table.good and table.poor linked in the Tonry 2001 paper. Read only the first three whitespace fields: preferred name, J2000 RA and Dec in degrees. Preserve the source file classification and all records; record file hashes and actual row counts without assuming the live files exactly reproduce the published sample count.

Match all 164 exposed-group representative positions with a fixed 30-arcsecond spherical radius. Retain every match, no-match and ambiguity. Call these positional source-membership candidates, not certified PGC aliases or CF4 contribution weights. Summarize existing fixed coarse-sky residuals for matched/unmatched rows and preserve distance ranges to expose confounding. Do not fit source offsets or change a distance.

Verify extraction invariance when every nonmetadata field is replaced. The source-format inspection displayed the header and first seven names/positions/velocity prefixes (N7814, N0063, N0147, N0185, N0221, N0224 and N0274), and PDF header inspection showed the first row's CMB velocity. These are not fresh source entries; record this limited exposure. No later source outcomes are to be displayed or used. This audit operates only on the already exposed 164-row evaluation data.
