# Coma: companion, MOND and NFW shape comparison

This is one exploratory cluster comparison, not a completed multi-cluster test or a transfer of the frozen galaxy law. Fit the first three bins and predict the last three. Smaller diagonal residual sums mean closer agreement; these are not calibrated probabilities.

| Shape | Inner residual sum | Outer residual sum | Scale at bound |
|---|---:|---:|---|
| transparent | 1.1622 | 2.7018 | False |
| strong_interception | 1.3084 | 2.5024 | False |
| nfw | 0.6510 | 4.7106 | False |
| mond_point_baryons | 0.6263 | 4.8299 | True |
| point_mass_shape | 5.0369 | 2.3261 | False |

## Predicted and observed shear

| Bin | Observed | Error | Transparent | Interception | NFW | MOND compact baryons | Point mass |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.0078475 | 0.0025339 | 0.0074813 | 0.0076212 | 0.0073560 | 0.0082295 | 0.0093438 |
| 2 | 0.0030932 | 0.0015000 | 0.0041112 | 0.0040793 | 0.0039939 | 0.0036012 | 0.0018021 |
| 3 | 0.0030509 | 0.0011907 | 0.0020685 | 0.0019414 | 0.0024521 | 0.0022183 | 0.0006853 |
| 4 | 0.0003729 | 0.0010212 | 0.0012089 | 0.0011067 | 0.0016813 | 0.0016001 | 0.0003569 |
| 5 | 0.0014153 | 0.0008983 | 0.0007834 | 0.0007074 | 0.0012330 | 0.0012481 | 0.0002172 |
| 6 | -0.0004491 | 0.0008051 | 0.0005489 | 0.0004915 | 0.0009518 | 0.0010246 | 0.0001465 |

The noisy outer bins permit several shapes. No result here establishes a physical winner. The MOND point-baryon approximation is especially restrictive for an extended, gas-rich cluster.

The compact MOND fit hits the lower scale bound. Its post-fit deep-MOND 1/R limit gives inner/outer sums 0.6201/4.8590. This is a boundary diagnostic, not an independently chosen new candidate.

See [protocol](cluster-model-comparison-protocol.md) for equations, provenance, missing inputs and the reason the one-third exponent cannot be tested by a freely normalized shear profile. All six research objectives remain open.
