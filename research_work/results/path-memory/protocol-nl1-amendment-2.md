# NL-1 amendment 2: the cluster root columns on stage 2's grid
Declared 19 September 2026 after the first complete archived run, before the rerun. That run is
preserved as `nl1-results-as-declared.json`.

## What failed
Gate G4 required the linear-law pressure columns on stage 1's 1500-point grid to agree with stage 2's
6000-point cache to 1e-3 for widths of 1 kpc and more; the worst difference was 4.3e-2, on the narrowest
included widths of the clusters whose columns stage 2 found hardest to converge. Gate R1's cluster part
required the local member at PM-1's two values of a* to reproduce stage 2's references under the
correlations to 1e-9; it reproduced them to 3e-5, the difference between the two grids. Both are grid
mismatches in the stage's construction, not properties of the family; the cluster solve's active members
(40 kpc and wider) do not use the affected columns.

## The amendment
The cluster root columns are rebuilt on stage 2's 6000-point grid, every width included as in stage 2
(no width excluded), from the same builder. Gate G4 is replaced by G4': the linear-law pressure columns
computed in the same build reproduce stage 2's archived 6000-point cache to 1e-12, and stage 2's
convergence record (G4-2, G4-3) is carried over. Gate R1's cluster part stays at 1e-9, now on the same
grid. Nothing else changes: the family, objectives, joint combination and decision rule are as declared,
and the galaxy, lens and joint numbers of the preserved run are deterministic.
