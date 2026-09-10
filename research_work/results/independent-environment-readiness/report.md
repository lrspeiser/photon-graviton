# Independent distances for an environmental redshift test

2026-09-10. Stage 2 input readiness and geometric audit; no redshift fit or new holdout evaluation. The previous goal turn made progress by deriving a conditional messenger-delay exclusion and committing its evidence. This turn obtains and audits a distinct environmental input rather than repeating that test.

We recovered 869 position/distance/method rows from the Updated Nearby Galaxy Catalog. The feature-only query requested no velocities or tidal indices. Method flags separate 90 Hubble-relation distances, which are unsuitable for independent spatial coordinates in this test, from directly estimated distance indicators. The catalog selection itself mixes distance and velocity criteria; it is not a complete independent-distance survey. [Primary catalog paper](https://arxiv.org/abs/1303.5328)

| Input classification | Rows | Use |
|---|---:|---|
| Individual stellar, SBF, supernova or geometric indicator candidates | 344 | Retain with calibration and uncertainty audit pending |
| Of those, inside nominal 11-Mpc region | 335 | Geometric path prototype only |
| Hubble-relation distances, including flow-adjusted variant | 90 | Exclude from independent coordinates |
| Unmeasured proximity estimates, flag txt | 20 | Exclude from measured coordinates |
| Membership, Tully-Fisher, Fundamental Plane or brightest-star cases | 415 | Retain separately pending method/association review |

All 869 names are unique. The returned table has 20 txt flags although the metadata note describes 11; preserve the actual row flags and discrepancy rather than silently altering the count. No method is rejected simply because it appears in standard cosmology: the issue is whether the target redshift or a disallowed dynamical inference supplies its position. The 344 candidates are not 344 fully audited, independent geometric distances; most are calibrated indicators, with shared errors and selection still to address.

## Coverage result

Using only the positions and stipulated distances of our already exposed 164 SBF groups:

- The median fraction of a target path inside 11 Mpc is **37.93%**.
- Only **3 of 164** target paths lie wholly inside that nominal radius.
- **125 of 164** have more than half the path outside it.
- The minimum nominal covered fraction is **11.80%**.

Even those figures are optimistic radial coverage, not certified density-map completeness. Obscured directions, unresolved faint galaxies, population selection, distance errors and the boundary remain problems inside the nominal region. Knowing some galaxy positions does not establish that every apparent gap is a physical void.

## Geometric calculation and its limited meaning

**Known geometry applied to an illustrative proximity proxy, not a new physical law:** along each target ray, place spheres of radius 0.5, 1 or 2 Mpc around the 335 retained tracer positions. Exact ray/sphere intersections and interval unions give the length with a catalog neighbor; its complement is the local no-neighbor length. All three radii are retained. None is selected by redshift agreement or adopted as a gravity parameter.

Let L=min(D,11 Mpc) and V_local be the no-neighbor length inside that range. If the binary proxy outside the range is unknown, its conditional full-path interval is

\[
V_{\rm local}\leq V_{\rm full}\leq V_{\rm local}+D-L.
\]

These bounds concern the declared catalog proxy only, not the true physical void length. Omitted tracers, including those just outside the boundary, can change the local classification. A physical density threshold, luminosity weighting, completeness and uncertain positions have not been inferred.

**Known regression identifiability:** in an empirical model ln(1+z_transfer)=aD+bV, assigning every unmapped segment to a void adds a term proportional to D. That can masquerade as a change in the distance coefficient a. Ignoring missing map coverage would therefore make the comparison circular in practice even though the retained tracer distances themselves are not redshift inversions. We have not fitted a or b to this incomplete proxy.

An independent midpoint integration checks the interval-union algorithm on nine target directions at each radius. Maximum absolute discrepancy is **0.000556 Mpc**, below the stated 0.01-Mpc numerical gate. This certifies numerical geometry only; it says nothing about astronomical density accuracy.

## What this changes next

The catalog is useful for developing nearby path geometry and for a possible shorter-distance test after source/selection review. It cannot supply the full environments of the current 10–93-Mpc sample. The next useful extension is a method-specific, independently estimated distance-tracer map across that larger volume, beginning with existing Cosmicflows distance-indicator columns rather than its reconstructed peculiar-velocity density map. Quantify selection and position errors before any environmental fit. A smaller nearby test would supplement, not replace, the full requested distance range.

The reconstructed Cosmicflows-3 Local Void map uses inferred peculiar velocities and a reconstruction framework; it is not an independent measured void field for our present test contract. Its published usefulness in standard cosmology does not remove that dependence. [Primary reconstruction paper](https://arxiv.org/abs/1905.08329)

New catalog distances and positions are now exposed features. No new target velocity outcomes were requested, and no source is certified fresh merely because this query omitted velocities. Existing 164-group observations remain exposed. Freeze a complete model and audit source/identity history before retrieving genuinely withheld outcomes.

Artifacts: [protocol](protocol.md), [calculation](run.py), [results and hashes](results.json), [distance feature register](distance-feature-register.csv), [all 492 path/radius rows](path-coverage.csv), [returned feature table](features.tsv), [source metadata](ReadMe). Text archives normalize line endings and trailing spaces; original download bytes remain in the ignored data cache. No gravitational response, photon conversion coefficient, companion ledger or previous evidence was retuned. All four goal stages remain incomplete.
