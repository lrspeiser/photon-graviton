# Refined width grids for the continuous-scatter timing revision

The two previously exposed faint artificial samples were recomputed on 641 width nodes, adding midpoint evaluations between all 321 original nodes. We reused the original Sobol shape points, injection seeds, cadence, uncertainties and shape priors. The first event's full original curve was regenerated in each sample and matched to below 1e-10; original arrays and new arrays are hash recorded. The new nodes are fresh event-likelihood calculations, not interpolation of the old ones.

The original 160-case experiment remains unchanged. This is an exploratory numerical diagnosis, not a new coverage experiment or an observed supernova result.

## Comparison

We refit the mean width and exponent at each fixed scatter, using three starts with each of L-BFGS-B and Powell. The best successful result is reported; all trials are retained. The development bounds remain mean width 5-80 days and b=-0.5 to 1.5, keeping the means inside finite support for these data. These bounds do not replace the original wider-domain experiment.

| Artificial sample | Fixed scatter | b, 321 nodes | b, 641 nodes | Absolute change |
|---|---:|---:|---:|---:|
| split_gaussian-snr5-b1-seed902 | 0 | 1.03144353 | 1.03333024 | 0.00188671 |
| split_gaussian-snr5-b1-seed902 | 0.01 | 1.03071186 | 1.03073687 | 0.00002501 |
| split_gaussian-snr5-b1-seed902 | 0.03 | 1.02999742 | 1.03002360 | 0.00002618 |
| split_gaussian-snr5-b0-seed903 | 0 | -0.04805370 | -0.04574792 | 0.00230578 |
| split_gaussian-snr5-b0-seed903 | 0.01 | -0.04425324 | -0.04430620 | 0.00005296 |
| split_gaussian-snr5-b0-seed903 | 0.03 | -0.04304080 | -0.04307688 | 0.00003608 |

At zero scatter the changes are about 0.00189 and 0.00231, smaller than the preceding 161-to-321-node changes. At scatter .01 or .03 the changes are below 0.000053. All six optimizer starts report success on each 641-node fixed-scatter fit. One start still fails in the 321-node zero-scatter no-stretch case; that failure is retained. Multiple successes do not prove global optimization, especially where interpolation creates corners.

Adding Powell also improves the selected coarse zero-scatter optimum slightly for seed 902 relative to the previous L-BFGS-B-only report. Therefore the earlier and current coarse values need not be identical. The result file retains method order (three L-BFGS-B starts followed by three Powell starts), objective values and parameters rather than hiding that change.

## What this establishes

The narrow-scatter calculation is less sensitive to width discretization at the finer grid in these two cases. This is evidence supporting further development, not a pass against a newly invented scientific threshold. We have not independently refined the nuisance-shape integration here, profiled every possible scatter, calibrated boundary confidence intervals, or proved applicability to real emitting-source populations.

All integration and statistical formulas are known mathematical methods. The population exponent b remains a phenomenological timing parameter, not a microscopic derivation of companion conversion. Next work must make the finite-support domain explicit, reliably optimize scatter including its zero boundary, and calibrate the resulting uncertainty procedure on a separately frozen artificial sample. The six observational objectives remain open.

Reproduce with `python research_work/results/timing-scatter-grid-refinement/run.py`. This recomputes both samples and writes generated arrays only in its own directory.
