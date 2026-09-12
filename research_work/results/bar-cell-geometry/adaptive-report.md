# Adaptive fixed-mass refinement: completed first comparison

All 218 new trajectories completed, alongside 480 reused paths. All 698 primary orbit checks pass; maximum sampled Jacobi error divided by 220 squared is 8.75e-8. The earlier polar-domain limitation in reused paths is retained. Geometry and force calculations hash-check the cached paths.

The 3,184-cell sphere retains all source regions. For this first adaptive comparison, each second-mesh parent region keeps its mass, distributed among children by solid angle. Thus this comparison tests geometry at fixed parent masses rather than simultaneously changing the capture-weight quadrature.

## Complete-source force results

Both 256- and 512-layer calculations at R=3 and T=0.1 finished. All ages from zero through T remain included, not just the interval used to choose refinement targets. Source mass differs from one by 2.23e-16.

| Target (kpc) | Force change, 256 to 512 layers | Force change, second mesh to adaptive at 512 |
|---|---:|---:|
| (0.1,0,0) | 0.540% | 32.242% |
| (1,0,0) | 0.099% | 0.508% |
| (1,0,1) | 0.013% | 0.192% |
| (0,1,1) | 0.019% | 0.163% |
| (3,0,0) | 2.430% | 0.958% |

All five age-layer force gates pass. Four of five spatial force gates pass, but the selected central target still fails the unchanged 5% threshold. All potential gates pass; the largest spatial potential change is 0.111%. These are changes between numerical representations, not error bars on observed gravity or proof of accuracy.

The previous 135.2% central change and the present 32.2% change compare different successive refinements, with a now-fixed parent-mass convention. A smaller change is useful evidence, but cannot be presented as a measured reduction in absolute physical error or as convergence.

Ten retained thin-cell/target cases agree with 60-digit arithmetic within 7.35e-13 potential and 5.10e-13 force in units per G times total source mass. This sampled roundoff audit does not bound all geometric errors.

## Geometry and limitations

At T=0.1, reference weighted RMS position error decreases from 4.23% to 3.72% of launch radius. At the final T=0.25 it decreases from 40.26% to 36.49%. The map remains poorly resolved at late times. Targets were chosen using exposed numerical failures, not independent astronomical data.

The central force still needs another controlled refinement; capture-quadrature sensitivity also remains untested for this adaptive mesh. Other ages, source radii and failing locations are not resolved by these five comparisons. No stellar fit, observed lensing result, funded self-gravitating source or redshift mechanism is established. All nine goals remain open.

## Reproduction

Run `check_adaptive_geometry.py`, `adaptive_volumes.py 256 3 0.1`, `adaptive_volumes.py 512 3 0.1`, `export_adaptive.py`, and `check_roundoff_adaptive.py` from this directory after preparation. JSON outputs retain all inputs, force components, comparison gates and failures. Every process from this adaptive preparation and first force comparison has finished.

This uses known numerical subdivision, interpolation and Newtonian force integration. The physical capture and force laws are unchanged; no novel formula is claimed.
