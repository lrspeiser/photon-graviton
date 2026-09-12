# Second full-sphere gravity comparison: completed

All four volume jobs finished successfully after both 480-path preparations completed. The second full-sphere refinement improves the geometric interpolation and resolves the tested age-layer comparisons, but still does not establish stable gravitational forces under spatial refinement.

## Direct comparison

| Comparison | Potential failures | Force failures | Largest potential change | Largest force change |
|---|---:|---:|---:|---:|
| 256 versus 512 age layers, second mesh | 0/30 | 0/30 | 0.218% | 2.356% |
| First versus second mesh, matched 512 layers | 6/30 | 14/30 | 6.470% | 135.223% |

Each set covers two source radii, three ages and five spatial positions. Gates remain 2% potential and 5% vector force, with the previously declared 0.01 denominators floors. They are numerical development criteria, not observational uncertainties. Passing age-layer checks at these points does not prove uniform time-integration accuracy everywhere.

The largest spatial force change is the 3 kpc source at age 0.1 and position (0.1,0,0). Its potential changes only 1.53%, illustrating why stable well depth alone is insufficient to claim a stable local pull. At final age the same position changes force by 56.76%. These are changes between numerical representations of one proposed source, not measured excess gravity or evidence against all companion models.

## Mass, arithmetic and geometry

Maximum discrepancy from unit source mass is 1.78e-15. Sixty actual thin-cell/target combinations were reevaluated with 60-digit arithmetic. Maximum absolute differences are 5.44e-11 in potential and 1.38e-10 in force per G times total source mass, passing the retained 1e-7 criterion. This sampled arithmetic audit does not bound geometric error or all cells.

The [completed geometry report](nested2-geometry-report.md) retains final weighted RMS position errors of 23.35% and 40.26% of launch radius. Two of the four targeted outer probes worsen. All primary trajectory checks passed, but no claim is made that every new path has an independent tighter-tolerance verification; the earlier polar-domain limitation is retained. Neither exact tetrahedron gravity nor conserved mass can compensate for an inaccurate density map.

## Consequence for the model

The present affine source representation remains inadequate for a reliable stellar prediction. The next numerical step must address curved and folded trajectory cells, using actual additional trajectory information and preserving mass, rather than retuning a physical force coefficient to hide numerical error. The completed comparisons and all failures remain the reference for judging that step.

There is no observed-star fit, physical source amplitude, energy-budget closure, self-gravitating evolution or galaxy-lensing result in this checkpoint. The independent conditional lensing-kernel tests do not change that. All nine research goals remain open.

## Reproduction and provenance

After both preparations, run `check_nested2_geometry.py`; run `nested2_volumes.py 256 1`, `256 3`, `512 1` and `512 3`; then `export_nested2.py` and `check_roundoff2.py`. `finish_nested2.py` performed the dependent sequence and terminated successfully. Raw volume JSON files, logs, summary JSON/CSV and roundoff results are retained.

All interpolation, integration and Newtonian kernel formulas used here are known mathematics. This is a numerical test of the proposed deposited-companion source, with unchanged source dynamics and capture parameters; it introduces no new physical law.
