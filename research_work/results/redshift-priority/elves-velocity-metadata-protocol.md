# Source metadata retrieval without redshift outcomes

Query the 26 currently pending ELVES positions using a fixed 30-arcsecond SIMBAD TAP cone. Retrieve source identifiers, object type, coordinates, internal database identifier and basic radial-velocity metadata (type, nature, quality, wavelength class and bibliographic reference) only. Explicitly exclude radial-velocity/redshift values and errors. A cone hit is a candidate counterpart, not an automatic identity match.

For returned objects with a basic velocity reference, retrieve measurement metadata: velocity type (v, z or cz), wavelength domain, bibliographic source and remarks/origin. Exclude value, meanerror and precision columns. This audits the present database, not an archived copy of the authors' original retrieval. Do not substitute a current preferred measurement for the adopted table value without checking its provenance and any required convention.

Retain every target, every returned counterpart and every failed request. Resolve identity using existing names, primary aliases and positions; do not use a favorable redshift match. The source-level review can identify conventions and selection issues but does not itself establish freshness, independent galaxy motion or the final uncertainty model. No pending outcome is to be fitted, predicted or scored in this pass.
