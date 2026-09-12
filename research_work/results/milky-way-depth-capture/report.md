# Milky Way comparison: a depth-dependent deposit model

## Result and main limitation

The p=3 capture-shape model, calibrated with one overall deposit amount, reduces the unweighted Cepheid training-bin circular-speed discrepancy from75.24 to8.76km/s. It still underpredicts the innermost bin by about21.6km/s. This is a useful conditional fit to real stellar-derived summaries, not an exact match, an independent prediction or validation of photon conversion.

The major physical problem is now measurable: the fitted deposit potential at the Sun is about3.5–3.9 times the ordinary-matter potential for p=3, and7.2–10.0 times for p=2. Capture was calculated from ordinary-matter depth alone. Its omitted gravitational feedback is therefore not a small correction. A coupled capture/gravity calculation is required before accepting these profiles.

## Data and assumptions

Only the542 previously processed training Cepheids and their12 radial bins were used. The selected-star cache contains training roles only; no reserved validation or test stars were opened. Predictions average squared circular speeds at those actual radii, then take the square root, matching the earlier training reduction.

The comparison values are approximate Jeans circular-speed proxies from measured stellar motions, not individual-star circular speeds. They retain the existing distance, frame, tracer-density, dispersion-gradient and axisymmetry assumptions, including the omitted vertical mixed moment. Equal bin weights and RMS are diagnostics, not statistical significances. Some outer bins contain few stars. No cosmological redshift-distance relation is used in this Cepheid training reduction.

The ordinary-matter field contains the existing stellar and gas disks, azimuthally averaged bar/nuclei and central black hole. Their masses and shapes are held fixed. These are conditional mass-model inputs with their own dynamical dependencies, not new independent measurements of all baryons. The old model's empirical gravity enhancement is not used: only ordinary-matter and Poisson/Legendre utilities are reused.

## Formula provenance

**Proposed capture law and optically thin limit:**

\[
\rho_D(R,z)=C\left[\frac{W(R,z)}{W(R,z)+W_*}\right]^p,
\quad W=\max(-\Phi_b,0),\quad W_*=40000\ (\mathrm{km/s})^2.
\]

p=2 or3, with outer cutoff100 or200kpc. The exponent and depth scale are specified phenomenological choices, not derived microscopic constants. C is the only fitted quantity for each variant; it represents the product of capture normalization, incident companion intensity and retention age. Fitting C establishes a required density, not an available energy supply.

The formula uses the leading isotropic thin-capture result q proportional to kappa. Taking kappa0 at most0.02/(2rmax), with units consistent with rmax, bounds any straight-ray optical depth by0.02 because the depth factor is at most one. Ignoring attenuation then changes local deposition by at most about2%. This is not a proven2% bound on every force component. Achieving the fitted C with very small opacity may require a large intensity-age product.

**Known Poisson gravity applied conditionally:** solve Laplacian Phi_D=4pi G rho_D for the axisymmetric density, then use the gradient of the same potential for radial and vertical acceleration. The orbital-plane prediction is v_c squared=-R a_R. This preserves the connection between potential and force; no independent radial/vertical correction is fitted. It does not establish mechanical support of the deposits or their full field energy.

## Calibration results

| Variant | Training RMS, km/s | Required deposit mass inside cutoff | Deposit/ordinary potential depth at Sun |
|---|---:|---:|---:|
| Ordinary matter | 75.24 | — | — |
| p=2,100kpc | 13.78 | 2.55e12 solar masses | 7.22 |
| p=2,200kpc | 13.78 | 6.02e12 solar masses | 10.02 |
| p=3,100kpc | 8.76 | 8.31e11 solar masses | 3.53 |
| p=3,200kpc | 8.76 | 1.22e12 solar masses | 3.85 |

Outer mass and absolute potential change markedly while the inner training fit barely changes. The12 radial bins cannot determine that outer mass, capture boundary, lensing profile or clock potential uniquely. These large masses are model requirements under ordinary gravitational coupling, not measured dark-matter masses or a demonstrated photon-energy budget.

Selected p=3,100kpc predictions below illustrate the residual pattern; [radial-comparison.csv](radial-comparison.csv) contains all12 bins and all four variants.

| Approximate radius | Inferred training proxy | Ordinary matter | Calibrated capture model |
|---|---:|---:|---:|
| 6.6kpc | 242.0km/s | 181.1km/s | 220.4km/s |
| 8.5kpc | 235.6km/s | 175.1km/s | 227.6km/s |
| 11.5kpc | 232.6km/s | 161.5km/s | 232.5km/s |
| 17.5kpc | 226.2km/s | 134.6km/s | 234.7km/s |

## Vertical comparison: excluded as a validation claim

The archived43 vertical-force estimates were evaluated without fitting to them. The historical diagnostic RMS decreases from27.99 to18.41 for p=2 and16.57 for p=3, in equivalent solar masses per square parsec. These numbers do **not** constitute an accepted vertical-observation pass under the user's research rule.

[Bovy & Rix2013, sectionIII.2](https://arxiv.org/html/1309.0809) explicitly uses a fitted potential family containing a spherical dark halo to infer those forces. They are model-derived estimates, not raw accelerometer readings. We preserve the comparison transparently as an alternate-model diagnostic, but exclude it from physical validation of this hypothesis. Underlying stellar positions/motions could instead be modeled with the candidate potential and survey selection.

The output also supplies radial and vertical force predictions at R=2,5,8.2,15kpc and heights0.3,1,3kpc. Those are predictions awaiting an admissible data comparison, not an established bulge result. The bar has been averaged over azimuth, so the full under/above-bar question remains open.

## Verification and next calculation

[Protocol](protocol.md), [solver](run.py), [results](results.json) and [table exporter](export.py) are archived. Inputs are hashed. Reproduce with:

```text
python research_work/results/milky-way-depth-capture/run.py
python research_work/results/milky-way-depth-capture/export.py
```

Coarse/fine source resolutions use600/1200 radii, multipoles through24/48 and96/192 angular nodes, with the same refined ordinary field. Maximum fitted radial prediction change is0.000523km/s and vertical diagnostic change0.0103 equivalent solar masses/pc squared, passing the declared numerical gates. Every variant is reported; none reached the imposed fit boundary. These convergence checks do not account for astronomical or model systematics.

The next required step is to let deposits alter the well that determines capture, then test whether the coupled density/potential stabilizes and retains the radial fit. A separate hypothesis that only ordinary matter controls capture would need an explicit physical justification; it cannot be silently substituted after discovering the feedback. Energy supply, support, persistent redshift, lensing, raw-data vertical modeling and a genuinely unseen test remain incomplete. The earlier redshift alpha and all reserved data roles are unchanged.
