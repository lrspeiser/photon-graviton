# Exact-one-third lens transfer

Before execution: run the existing six-system lensing calculation with the exact-one-third galaxy capture parameters. Prior lens results used q=0.3399907823, not exactly 1/3. This closes that version mismatch; it is not a new blind sample.

Use third-radiation-retention-results.json, preserving its three shared training-fitted capture constants. Keep regular optical geometry, the same population-mass-to-luminosity proxies, equivalent disk scales and observations as the earlier retention lens test. Fit only stellar mass and constant orbital anisotropy to the inner stellar bins. Predict the conditional outer stellar bin and lens angle. No capture, optical or retention parameters are fitted on the lens targets. Both Chabrier and Salpeter proxy cases must be reported.

Check identical observed stellar data and geometry relative to the previous run, input hashes, bounded eta, fixed q=1/3, and optimization outcomes. Report the full per-system predictions and aggregate scores without selecting the better population proxy. The published lens angles are model summaries and the proxy luminosities are not direct 3.6-micron measurements; the test remains conditional.

The Hill function is known mathematics. Its use as a radiation-conditioned retention probability is proposed physics, and exact one-third has not been derived from first principles. A successful lens transfer would still not establish a source energy budget, support mechanism, or cluster prediction.

Run: python research_work/results/isotropic-galaxy-transfer/lensing.py --third-retention-optics
