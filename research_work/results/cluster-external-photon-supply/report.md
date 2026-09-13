# External photon sources feeding a spherical cluster

## What is now calculated

A photon-emitting external point source is propagated to a spherical receiver with exact angular interception, conversion before entry, further conversion inside, and companion capture. This replaces an arbitrary companion bath with a per-source photon-energy transfer fraction. It complements internal-photon-supply and capture-startup; it is not an actual cluster catalog calculation.

The calculation follows a finite isotropic source pulse until every intersecting ray has exited or deposited its energy. It does not select a universe age or size. For a finite observation time, only causally arrived portions of source histories may be counted; multiplying the completed-pulse fraction by all historical emission would overcount recent or not-yet-arrived energy.

## Postulates versus known mathematics

Optional project assumptions: fixed source at D>R; isotropic photon emission; same-direction photon and companion travel at c; constant fractional photon conversion alpha; lossless companion flight outside the receiver; constant companion capture kappa inside; permanent storage; no intervening absorbers, focusing, scattering or backreaction. These are not established interactions or claimed unique inventions. Prescribed geometry omits the mechanical/gravitational work needed by a complete theory.

Exact Euclidean solid-angle geometry gives the fraction of source emission whose rays hit the sphere:

    f_hit = [1-sqrt(1-(R/D)^2)]/2.

The far-source approximation is R^2/(4D^2). For a source immediately outside the surface, the approximation is inaccurate and the intercepted fraction tends to one half, not one quarter.

Define dimensionless d=D/R and y=half-chord/R, with 0<=y<=1. Entry distance and receiver chord length are

    s_entry/R = sqrt(d^2-1+y^2)-y,
    L/R = 2y.

The fraction of isotropic source energy in dy is y dy/[2d sqrt(d^2-1+y^2)]. This integrates to f_hit. Before entry the photon and companion fractions on that ray are P0=exp(-alpha s_entry), C0=1-P0. Inside, the known coupled transfer equations of the preceding internal-source calculation give

    P_out = P0 exp(-alpha L)
    C_out = C0 exp(-kappa L)
            + P0 alpha [exp(-alpha L)-exp(-kappa L)]/(kappa-alpha)
    D_retained = 1-P_out-C_out.

The equal-rate limit is used explicitly. Averaging these over the exact solid-angle measure yields fractions per unit total emitted photon energy. Their sum equals f_hit; the remaining 1-f_hit never intersects this receiver. This is established transport/geometry mathematics applied to hypothetical sectors, not a new fundamental law.

## Results

The table uses alpha R=0.0002488993265191759 and kappa R=10. The conversion value corresponds to the archived empirical alpha if the illustrative receiver radius is 1 Mpc; it is not a derived constant or measured cluster size. Distances are sensitivity examples, not catalog objects or inferred cosmic boundaries.

| D/R | Fraction of all emission intercepted | Fraction of all emission deposited | Deposited / intercepted |
|---|---:|---:|---:|
| 1.01 | 0.4298146 | 0.000109859 | 0.0256% |
| 2 | 0.0669873 | 4.024917e-05 | 0.0601% |
| 10 | 0.002506281 | 6.536693e-06 | 0.2608% |
| 100 | 2.500063e-05 | 6.149409e-07 | 2.4597% |
| 1000 | 2.500001e-07 | 5.48375e-08 | 21.9350% |
| 10000 | 2.5e-09 | 2.281083e-09 | 91.2433% |

Greater distance raises the converted share of the radiation on an intersecting ray but reduces the solid angle occupied by the cluster. Increasing the source count must be evaluated with both factors. All 18 scenarios, including kappa R=0.1 and 1, are saved in results.json.

For established constant illumination, each source contributes deposited power L_star*f_deposit. This is linear superposition only in the fixed-opacity model. To estimate a real reservoir, integrate source luminosity histories with causal delays and time-dependent capture, and track any competing receivers. Photons exported by an internal source must not also be credited as local deposits unless a return path is calculated.

For distant uniform shells of sources, their number grows as D^2 while geometric interception falls as D^-2. Consequently adding shells does not necessarily produce a convergent supply. An infinite eternal uniform source population with no intervening removal gives divergent accumulated illumination in this branch. This does not require imposing an arbitrary universe size; it requires a specified, convergent source history/transport model or finite causal exposure before claiming a finite prediction. No source history is inferred from desired gravity here.

## Verification and limitations

Two Gaussian quadrature resolutions and independent adaptive integration agree within 1e-9 of intercepted energy in every case. Matrix exponentials independently reproduce the ray solution within 1e-12; zero conversion, zero capture and equal-rate cases are exercised. Energy fractions sum to the intercepted fraction; deposits never exceed it. These are mathematical checks, not observed cluster validation.

Next: populate the kernel with source luminosities/positions and causal history, include intervening competition, and retain the spatial deposition field for a common dynamical/lensing prediction. Actual cluster accumulated energy, capture law, supported reservoir, and image distortions remain unresolved. No final astronomical holdouts were opened; all six research objectives remain open.

Run `python research_work/results/cluster-external-photon-supply/run.py`.
