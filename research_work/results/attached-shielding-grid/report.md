# Transfer of pair-calibrated shielding across the 3D source grid

**The pair-calibrated rule does not reproduce the full source shape.** Even with a best common normalization for the inner subset, at least one inner source-density error is 59.22%; the corresponding full-grid bound is 99.65%. Maximum transmission refinement change is 0.852%. These are conditional source-density mismatches, not measured velocity errors.

The opacity fitted to the two critical locations is held fixed at 9.675692285720144e-10 kpc^2/Msun. Isotropic straight-ray absorption is calculated at all 240 existing model probes, with boundary 30 kpc and the unchanged baryonic density. No observational holdouts are opened.

The optional closure is rho_extra=C T(x) rho_b(x), where T is the calculated angle-averaged transmission and C is one common exposure normalization. The opacity is not refitted. For each reported subset we allow its own best C only to establish a generous lower bound on mismatch; these separate constants cannot jointly define one galaxy model. A single whole-grid constant cannot do better on any subset than its subset-specific optimum.

Known minimax algebra applied to q/T gives the lowest possible worst relative source-density mismatch: (max(q/T)-min(q/T))/(max(q/T)+min(q/T)). This is a model-shape bound, not a velocity residual or confidence level. All source values, transmissions and resolution comparisons are retained.

| Subset | mu nodes | Points | Best possible worst source-density mismatch |
|---|---:|---:|---:|
| full | 16 | 240 | 99.651% |
| full | 32 | 240 | 99.654% |
| inner | 16 | 120 | 59.038% |
| inner | 32 | 120 | 59.221% |
| midplane | 16 | 40 | 55.367% |
| midplane | 32 | 40 | 55.334% |

Maximum relative transmission change between the two resolutions is 0.00851658. These are development resolutions, not proof of convergence at every location.

The favorable two-location ratio was calibrated and does not establish the full shape. Any residual mismatch can reflect the optional constant-opacity/common-history/attached-storage assumptions, or the empirical target itself. It does not exclude all companion capture. Actual energy supply, momentum transfer, storage dynamics, gravity feedback and motion/lensing predictions remain unclosed. All six objectives remain open.
