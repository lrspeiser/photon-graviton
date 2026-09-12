# Locating the largest force discrepancy by source age

The failed quadratic remedy motivated a more targeted diagnosis. Recompute the largest spatial force discrepancy, for the 3 kpc source at age T=0.1 and target (0.1,0,0), separately for the first and second nested meshes. Keep 512 age layers, the exact tetrahedral kernel, unchanged capture weights and fourfold symmetry.

Split the age integral into eight equal bins. Every bin keeps 12.5% of the original source mass; none is independently normalized or discarded. Summing bin potentials and forces reproduces each archived full calculation within 1e-10 absolute units per G times total source mass. This verifies the decomposition, not the physical energy supply.

## Result

The total force changes from (-0.0543401,-0.0234074,0) to (-0.0212433,-0.0140047,0). The bin from age 0.0625 to 0.075 contributes (0.0302964,0.0095590,0) to this difference. Age units are kpc/(km/s); that interval is approximately 61–73 million years after launch.

This one bin supplies **92.30% of the signed difference projected along the total force-change vector**, and **90.34% of the sum of the eight difference-vector magnitudes**. The definitions differ because contributions can point in different directions:

\[
\Delta\boldsymbol a=\sum_k\Delta\boldsymbol a_k,\qquad
P_k=\frac{\Delta\boldsymbol a_k\cdot\Delta\boldsymbol a}{|\Delta\boldsymbol a|^2},\qquad
Q_k=\frac{|\Delta\boldsymbol a_k|}{\sum_j|\Delta\boldsymbol a_j|}.
\]

These are known linear superposition and vector-projection identities, not novel physical laws. P can be negative and sums to one; Q is nonnegative and sums to one. They measure concentration of numerical differences, not the fraction of real gravity produced by a population.

## Consequence

For this failing target, resolving the spatial map during that passage is a better-informed next test than increasing time layers or discarding old material. The result does not show that the other 13 failed force comparisons share this cause. Nor does it identify a physical capture age, justify a lifetime cutoff, or prove a caustic. All source ages remain in the model.

Next examine angular contributions within the identified interval and compute new trajectories where needed, then recompute the total force with all other ages retained. A numerical error must not be hidden by changing capture efficiency or gravitational coupling.

All nine goals remain open. There is still no converged stellar/lensing prediction, funded self-gravitating source or successful linked redshift mechanism.

## Reproduction

Run `age_force_attribution.py 1`, `age_force_attribution.py 2`, then `summarize_age_attribution.py` from this directory. Both force decompositions finished successfully. Source hashes, individual age-bin forces and weights, total potential/force, and normalized attribution are retained in the JSON files. No process from this audit remains active.
