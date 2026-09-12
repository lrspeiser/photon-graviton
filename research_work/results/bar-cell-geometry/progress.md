# Nested-mesh preparation in progress

The previous goal turn completed the first trajectory-volume calculation and was progress. This turn completed the angular geometry audit and constructed a nested mesh, then launched actual new trajectories.

Launch radius1: session40550. Launch radius3: session79270. Each prepared-R*.json must contain98 records, with18 reused and80 new trajectories, and the process must be terminal before a complete preparation claim. Recheck these handles or the prepare_nested.py process command line before restart; partial manifests are not completion evidence. Partial prepared-R1.json/prepared-R3.json are excluded from the checkpoint commit.

Next verify new trajectory/cache accuracy and use mesh.json faces and lookup/sign arrays to assemble the new volume mesh. Compare geometry and force against the parent without dropping source mass. No new-grid force result exists yet. The exact polar limit is tested, but the original radial domain still applies. All nine goals remain incomplete.
