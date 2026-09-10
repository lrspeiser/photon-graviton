# Training-star orbit support and distance-input audit

**The training test uncovered a missing input-provenance distinction before any stellar gravity ranking.** Our original preparation retained StarHorse's output-quality flags but omitted its input flags. Some distance estimates did not use Gaia parallax, so our earlier blanket statement that StarHorse had already used that measurement was too broad. The corrected metadata now records this distinction and flags distance/parallax disagreement. Original catalogs and holdout assignments remain unchanged.

We also tested and fitted an initial one-dimensional orbit-library diagnostic to **27,884 real training stars**, including 1,675 in the bulge region. This confirms that the six demonstration orbits are not a sufficient stellar population. Seventy-two training-derived launch states improve coverage, but are not integrated or validated as a complete orbit library. The newly found input issue takes precedence over interpreting their fitted weights.

## Why the input issue matters

**Known geometric/unit-conversion relation, not a new theory formula:**

\[
v_\perp[\mathrm{km/s}]=4.74047\,\mu[\mathrm{mas/yr}]\,D[\mathrm{kpc}].
\]

The same measured sideways angular motion implies a larger speed when assigned a larger distance. For APOGEE target **2M11112524-6128593**, the catalog supplies a StarHorse median distance of 6.081 kpc, while its raw Gaia parallax is 1.75875 +/- 0.01046 mas. The simple inverse-parallax scale is about 0.569 kpc. Using those two distance scales changes the inferred heliocentric transverse speed by a factor of about **10.7**. The inverse value is an illustration of the discrepancy, not a newly adopted distance posterior. We have not decided which input or association is wrong.

Four launch representatives have previous approximate Galactocentric speeds above 500 km/s, reaching about 1,387 km/s. All four lack the StarHorse PARALLAX input token. Their independently retrieved astrometric fidelity scores are high (v1 at least 0.9238, v2 at least 0.9951). They cannot simply be dismissed as low-fidelity Gaia solutions or promoted as evidence for extra gravity. Crossmatch, photometry/spectroscopy, distance modeling and astrometric reliability must be assessed together. `input-audit.json` contains their exact public identifiers and source values.

The [StarHorse release documentation](https://data.aip.de/projects/aqueiroz2023.html) identifies PARALLAX as an input flag and supplies distance percentiles in kpc. The [paper's input-data section](https://arxiv.org/html/2303.09926v1#S3.SS1) says that low fidelity_v1 values are one reason to omit parallaxes, and describes recalibrated parallax uncertainties. The high-fidelity examples above show that this reason does not explain every omission in our matched sample. The flags establish which inputs were recorded as used; this audit does not claim to reconstruct every internal decision of the original pipeline.

## Completed provenance repair

We joined the original input flags and 5th/95th distance percentiles using **both APOGEE ID and exact integer Gaia source ID**. We also downloaded fidelity_v1 and fidelity_v2 for all **140,407** prepared source IDs from the authors' [GAVO table](https://dc.g-vo.org/tableinfo/gedr3spur.main). All requested IDs returned, and the identifier/role checks passed. These classifier scores are quality indicators, not guarantees that a distance or crossmatch is correct.

The new `stellar-input-provenance.parquet` is a sidecar to the frozen catalog. It leaves distances, velocities, original quality cuts and spatial holdout roles intact. It distinguishes parallax-used cases, cases needing separate likelihood treatment, distance/parallax disagreements, and astrometric-quality review. An audited copy of the 72 launch states carries those flags; **42 require additional input review** before use in a gravitational inference. None was silently dropped.

Our overlap diagnostic compares the StarHorse 5th–95th distance interval converted to parallax, widened by the previously declared zero-point range [-0.037,+0.003] mas, with Gaia parallax intervals of 3, 5 and 10 reported error widths. **This is not an independent Gaussian significance test:** the products can share inputs and have non-Gaussian/model/systematic uncertainties. Disjoint intervals trigger review rather than prove which measurement is wrong.

| Selected training subset | Stars | Disjoint at 5 reported Gaia error widths |
|---|---:|---:|
| StarHorse recorded parallax use | 16,756 | 111 |
| StarHorse did not record parallax use | 11,128 | 959 |
| Total | 27,884 | 1,070 |

In the full 77,927-star training allocation, 19,381 lack the PARALLAX token. This is not a claim that all those distances are invalid. A distance based on spectroscopy and photometry can be useful; the joint uncertainty and selection treatment must reflect its actual inputs.

For the selected population, dropping every parallax-omitted estimate would affect the proposed comparison very unevenly:

| Training region | Stars | StarHorse recorded parallax use |
|---|---:|---:|
| Bulge plane: R=0.5–3.5 kpc, abs(z)<0.2 kpc | 299 | 36 |
| Bulge off-plane: same R, abs(z)=0.5–1.5 kpc | 577 | 440 |
| Disk plane: R=5–9 kpc, abs(z)<0.2 kpc | 4,425 | 2,472 |
| Disk off-plane: same R, abs(z)=0.5–1.5 kpc | 11,214 | 7,667 |

These counts use the existing approximate posterior mean positions and the selected chemistry cuts; they are not identical to the earlier all-chemistry, median-distance region counts. This location-dependent provenance is a selection issue that could otherwise masquerade as a bulge height effect.

## Executed orbit-library diagnostic

The predeclared subset uses existing training/quality assignments, R=0.5–9 kpc, abs(z)<=1.5 kpc, FE_H>=-0.5 and ALPHA_M<0.15. Positions, velocities and covariance come from the prior conditional-minus-0.017 sensitivity reconstruction. Because the input audit above was discovered during this test, its results must remain marked **pre-audit diagnostics**, not validated gravitational measurements.

**Known rotating-bar invariant and first-order uncertainty propagation:**

\[
J=\tfrac12(v_R^2+v_\phi^2+v_z^2)+\Phi(R,z,\phi)-\Omega Rv_\phi,
\qquad \sigma_J^2\simeq (\nabla J)^T C(\nabla J).
\]

This J remains constant along each orbit of a fixed rigidly rotating potential. Thus six trajectories supply only six invariant values. The propagated sigma is approximate: it uses a linearization of existing posterior moments, not an exact observational likelihood or a confidence guarantee. Changing the potential changes the projected J itself, so marginal J-density scores cannot be used to rank gravity theories.

| Potential used for this diagnostic | Stars farther than 3 approximate error widths from the six invariant values | Same diagnostic for 72 launch values |
|---|---:|---:|
| Ordinary matter | 68.8% | 43.4% |
| Halo comparison | 91.1% | 35.2% |
| Equatorial deposits | 68.3% | 42.3% |
| Upper/lower deposits | 68.5% | 42.1% |
| Shell deposits | 68.4% | 42.1% |

The 72 launch states are actual training representatives selected in nine radius/height strata by deterministic farthest-point covering of robustly scaled position/velocity coordinates. They use the same physical states for all five potentials. They are not equilibrium weights, a complete six-dimensional basis, or independent observations reserved for testing. The coverage percentages apply to a literal unbroadened finite mixture; a continuous distribution-function or orbit-interpolation method requires its own resolution and calibration tests rather than necessarily one orbit per measurement-error width.

**Known finite-mixture statistics, not a companion law:**

\[
p(J_i\mid w)=\sum_k w_k\,\mathcal N(J_i;J_k,\sigma_{J,i}^2+s_{\rm library}^2),
\quad w_k\geq0,\quad\sum_k w_k=1.
\]

We fitted the six- and 72-component marginal mixtures with fixed additional widths zero and 1,000 (km/s)^2. This width is a numerical smoothing diagnostic, not measured extra gravity or a new measurement error. More components and smoothing improve in-sample density descriptions without proving the physical model. The test does not include the required position-conditioned velocity distribution, phase mixing, survey selection, tracer-density model or full distance likelihood.

Some EM fits reached their iteration cap, and a subsequent log-weight optimizer also reached its cap for broad, nearly redundant components. Those outcomes are retained. Direct constrained weight optimization then terminated for all ten expanded-mixture cases. The convex objective-gap bounds are at most **6.53e-5 nats per star** for the declared 1e-10 minimum weights; this bounds optimization error for that marginal problem, not astrophysical error. Weights need not be uniquely determined. `results.json`, `mixture-refinement.json` and `mixture-certification.json` preserve the stages instead of reporting the first capped run as a final optimum.

An independent synthetic-mixture check agrees with constrained optimization to 1.70e-6 in fitted weights. Finite differences verify the Jacobi gradient to 9.32e-10 in the stated scaled check. These verify numerical machinery, not the adequacy of the data or theory.

## Revised next steps

1. Route each star's distance/astrometry treatment by its recorded input provenance and fidelity. For a parallax-used distance, avoid counting the same measurement twice. For a parallax-omitted distance, a reliable additional parallax likelihood may be usable after checking independence and association; do not assume it was already included. Strongly disagreeing cases need input review or a declared alternative/error model, not an automatic inverse-parallax replacement or velocity-based rejection.
2. Preserve regional completeness and model the resulting selection. Do not obtain a seemingly clean bulge-plane comparison by discarding most of its parallax-omitted stars without accounting for that change.
3. Build a sufficiently broad integrated/interpolated orbit population and fit the shared position/velocity likelihood on training data. The current launch set is provisional pending its input review. The separate [likelihood contract](../rotating-bar-orbits/likelihood-contract.md) remains the observational target.
4. Freeze the complete model and procedure before opening reserved test scores. Keep ordinary-matter uncertainty, the shared companion response, redshift/timing, transport, capture, lensing and the deferred total-energy budget in scope.

Reproduce the pre-audit support run with `run.py`, its optimizer checks with `verify.py`, and the retained weight refinements with `refine_mixtures.py` followed by `certify_mixtures.py`. Retrieve quality metadata with `fetch_fidelity.py`, create the corrected sidecar using `input_audit.py`, and verify preservation using `verify_inputs.py`. Large source and per-star tables remain in the ignored cache; scripts, protocols, hashes and summary evidence are tracked. No held-out kinematic scores were evaluated. The full research goal remains active.


## Subsequent association finding

The [epoch-aware follow-up](association-report.md) finds that all four extreme training examples remain positionally discrepant after proper-motion correction. Their inherited associations require review before distance replacement or orbital inference; see association-review.json. No original data or roles were changed.
