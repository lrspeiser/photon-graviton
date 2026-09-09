# Radio identity and quality audit before target outcomes

Use the Yu et al. 2022 CDS table2 metadata for all 26 pending ELVES targets, retaining every target. Before inspecting matches, fix the search to a 30-arcsecond cone about each publisher position. Read only AGC identifier (bytes 5-10), optical-counterpart RA (12-20), Dec (22-29), and notes (150-158). Do not parse distances, central velocities, uncertainties, fluxes or derived physical quantities. Keep the raw download in the projectless work cache and record its SHA256; never print raw rows.

Retain every cone match. Explicit AGC aliases, or UGC identifiers below 100000 as defined by the source, can support a catalog link for an already identity-supported SIMBAD counterpart. A position-only match does not resolve an unverified galaxy identity or historical PGC group link. No new freshness certification or quality exclusion is made in this diagnostic.

Record source flags: 1/2 detection class, c overlapping HI emission, p projected neighbor without declared HI confusion, r low spectral weight at an emission edge. Their presence concerns the cited radio measurement, not necessarily the different measurement adopted by ELVES. Missing matches/flags are not proof of isolation or high quality. No flag is an independently measured line-of-sight void exposure.

Verify that replacing every byte outside the allowed slices leaves all parsed metadata unchanged. Save the filtered response and decisions for offline inspection. Use this to identify source/selection problems before final evaluation, without choosing objects by prediction residual.
