# Targeting inputs acquired for a pilot field

Field 300+00 (lco25m) was selected by its high training-member count, not its velocities. All 521 main-sample training stars match its downloaded targeting parent, which contains 189,661 unique identifiers. Location ID 5540 links to one field-summary row and eight design records.

Acquired allField (2,702 rows), allDesign (4,149), apogee2Plate (3,587) and the field's apogee2Object catalog. URLs, sizes, SHA256 hashes, schemas and embedded-checksum availability are recorded. The plate file emits nonstandard SURVEY/DATE header warnings; missing embedded checksums are not reported as successful checks.

Sources: [SDSS data-access documentation](https://www.sdss4.org/dr17/irspec/spectro_data/), [targeting archive](https://data.sdss.org/sas/dr17/apogee/target/), and [object data model](https://data.sdss.org/datamodel/files/APOGEE_TARGET/apogee2Object/apogee2Object.html).

The field summary reports eight expected and completed designs, with 24 expected visits. The designs allocate short and medium cohorts, with H boundaries approximately 10, 12.2 and 12.8; long-cohort allocation is zero. Their two active dereddened-color bins begin at J-K=0.5 and 0.8. Full design metadata are retained. This does not prove every cohort's actual completion or a star's inclusion probability.

The ratio 521/189661 is not completeness: its numerator is our filtered training subset, while its denominator includes sources outside relevant targeting cuts. Next we need the full observed main-sample numerator, parent quality/color/magnitude filtering, repeated-cohort and spatial-coverage handling, and the project's additional Gaia, StarHorse and chemistry selection. Dust and stellar-population modeling are needed to translate apparent-magnitude selection into spatial completeness.

No validation or final-test kinematics were analyzed. This is selection-model groundwork, not a new gravity measurement or companion-mechanism validation. All six objectives remain open.
