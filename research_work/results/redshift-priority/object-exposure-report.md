# Structured object exposure audit

The registry scans recognized identity fields in tracked CSV/JSON files, keeps galaxy names, PGC object IDs, PGC group IDs and DES supernova IDs in separate namespaces, and records source hashes. The CF4 join reads only object/group identifier columns. It does not use redshift outcomes to decide overlap. Run exposure_registry.py from any working directory inside this checkout.

SIMBAD identifier responses resolve all six maser targets (source URLs and identifier sections are preserved in maser-aliases.json). PGC/LEDA aliases reveal three exact earlier pilot objects that the previous name-only search missed:

| Maser | PGC | CF4 group | Earlier pilot object hit |
|---|---:|---:|---|
| UGC 3789 | 20679 | 20679 | Yes |
| NGC 6264 | 59306 | 59332 | No |
| NGC 6323 | 59868 | 59927 | No |
| NGC 5765b | 53012 | 53260 | Yes |
| CGCG 074-064 | 50048 | 50028 | Yes |
| NGC 4258 | 39600 | 39600 | No; shared calibration anchor |

The three pilot matches materially strengthen the conclusion that the six-object maser exercise is a reused/external diagnostic, not a fresh validation. All six labels have already been scored, so the other three cannot now become a blind sample either. Adding their groups expands the conservative group exclusion set from 232 to 235.

The structured inventory contains 233 PGC object IDs, 232 explicitly recorded group IDs, 180 normalized galaxy names and 2,136 DES event IDs. These counts are namespaces, not disjoint astronomical objects: names may alias PGC objects and supernovae may share hosts. The generated summary records current source counts and errors; duplicate recovered files preserve provenance without multiplying identities.

This is a completed structured-field pass, not a complete historical-exposure certification. PDFs, archives, unstructured prose, unnamed arrays, general galaxy-name aliases and supernova-host joins remain unaudited. On-disk availability is conservatively flagged, not equated with proof every label was inspected. No absence from this registry establishes freshness. A future sample still requires source/version selection, remaining identity and group checks, an independently justified observation model, and frozen predictions before label access.

The earlier maser-name-overlap.json and 232-group file remain historical audit snapshots. Use object-exposure-registry.json for the expanded structured exclusion set and the explicit limitations.
