# Central visits and residence in the ordinary rotating bar

This extends the earlier eight first-pass probes to angular quadrature on two launch spheres and follows subsequent passages. Turning outside the center on the first approach does not establish permanent avoidance. Visiting the center also does not establish permanent accumulation: the relevant central mass fraction depends on time spent there.

## Scope and formula provenance

The experiment uses only the ordinary stellar bar component of the existing analytical field cache, with rotation rate 37.5 km/s/kpc. Other Galactic components and deposited self-gravity are absent. The upstream [AGAMA bar example](https://raw.githubusercontent.com/GalacticDynamics-Oxford/Agama/f302756b8af2b763db58e278e30478517dc8eea3/py/example_mw_bar_potential.py) supplies a model-dependent ordinary component; separating its halo does not make the stellar normalization independent of the original dynamical modeling. This is a conditional mechanics test, not new astronomical evidence.

**Project source postulate, not a unique derived law:** at each launch radius R=1 or 3 kpc, angular weights are proportional to

\[
w(\theta,\phi)\propto\left[\frac{W(\theta,\phi)}{W(\theta,\phi)+W_0}\right]^6,
\quad W=-\Phi_{\rm bar},\quad W_0=40000\ ({\rm km/s})^2.
\]

This is the thin-absorption limit with uniform incident intensity. It does not solve nonspherical radiation attenuation or derive the prescribed inward speed of 50 km/s. Each sphere is normalized separately; the two spheres are not a volume source or a measured population. Continuous injection is stationary in the rotating bar frame.

**Known mechanics:** in rotating axes, with p the inertial velocity expressed in those axes,

\[
\dot{\boldsymbol{x}}=\boldsymbol{p}-\boldsymbol{\Omega}\times\boldsymbol{x},\qquad
\dot{\boldsymbol{p}}=-\nabla\Phi_{\rm bar}-\boldsymbol{\Omega}\times\boldsymbol{p},
\qquad E-\Omega L_z=\mathrm{constant}.
\]

These are ordinary rotating-frame equations, not a new companion interaction. Following the prescribed cold-particle injection is conditional on that interpretation of deposits.

**Conditional population accounting using known age integration:** let a be particle age, t_i its first entry into r<0.1 kpc, and normalized weights sum to one. At source age T,

\[
f_{\rm cohort}(T)=\sum_i w_i\,\mathbf1(t_i\leq T),\qquad
f_{\rm history}(T)=\sum_iw_i\frac{\max(T-t_i,0)}{T},
\]
\[
f_{\rm inside}(T)=\frac1T\sum_iw_i\int_0^T\mathbf1[r_i(a)<0.1\ {\rm kpc}]\,da.
\]

Missing first entries contribute zero. The first quantity describes a single birth cohort; the second describes the fraction of continuous arrivals with a central visit in their history; the third is the fraction of a continuously supplied population actually inside now. None assumes destruction or sticking at the diagnostic boundary. T=0.05, 0.1 and 0.25 kpc/(km/s) correspond to approximately 49, 98 and 244 Myr.

## Numerical audit

The protocol was recorded before each extension in [protocol.md](protocol.md). Angular quadratures use 4 by 8 and 6 by 12 nodes before exact reflection/half-turn reductions, giving 16 and 36 trajectories across the two spheres. L40 and L64 field orders are compared. These are deterministic quadrature comparisons, not sampling confidence intervals.

Initial first-entry files and the first continuation attempt remain available. A continuation missed a brief entry and exit within one integrator step despite small energy error. The corrected `audit.py` integrates whole trajectories and locates boundary crossings on intervals bracketed by radial turning points. It compares two tolerances and records any extra refinement; conserving the Jacobi quantity alone is insufficient to check central passages.

The spherical control uses independent radial-energy quadrature at 128 and 256 nodes. Its potential is held at the innermost cached radius below that radius; it is not a resolved theory of the Galactic center. Its reported central fraction changes below 2e-14 between quadratures. The control has no angular momentum and no sink at the center.

Authoritative corrected data are `audit-*.json`, `summary.json` and `summary.csv`. The original `coarse-L64.json`, `fine-L*.json` and `residence-*.json` document the earlier stages and must not override corrected crossings.

## Results

Percentages below refer to a continuously supplied population at approximately 244 Myr. The source on each launch sphere is normalized independently.

| Launch radius | Model/grid | Ever visited center (%) | Currently inside (%) |
|---|---|---:|---:|
| 1 kpc | Spherical control | 97.432 | 4.8688 |
| 1 kpc | Bar coarse, L64 | 35.860 | 0.2058 |
| 1 kpc | Bar fine, L64 | 63.614 | 0.5142 |
| 1 kpc | Bar fine, L40 | 63.615 | 0.5141 |
| 3 kpc | Spherical control | 92.507 | 1.1985 |
| 3 kpc | Bar coarse, L64 | 15.270 | 0.0460 |
| 3 kpc | Bar fine, L64 | 9.682 | 0.0466 |
| 3 kpc | Bar fine, L40 | 9.681 | 0.0466 |

All 88 final orbit-tolerance gates pass. 10 of 36 source/field statistical gates fail; they are retained below. Passing the absolute residence threshold of 0.5 percentage points does **not** establish small relative error: the 1 kpc coarse/fine residence estimates differ by about a factor of 2.5 at the final epoch.

| Failed comparison | Radius (kpc) | T | Statistic | Absolute difference | Limit |
|---|---:|---:|---|---:|---:|
| angular | 1 | 0.05 | cohort_ever_entered | 0.267176 | 0.05 |
| angular | 1 | 0.05 | continuous_ever_entered | 0.095836 | 0.05 |
| angular | 1 | 0.1 | cohort_ever_entered | 0.371302 | 0.05 |
| angular | 1 | 0.1 | continuous_ever_entered | 0.191764 | 0.05 |
| angular | 1 | 0.25 | cohort_ever_entered | 0.095446 | 0.05 |
| angular | 1 | 0.25 | continuous_ever_entered | 0.277545 | 0.05 |
| angular | 3 | 0.05 | cohort_ever_entered | 0.089215 | 0.05 |
| angular | 3 | 0.05 | continuous_ever_entered | 0.055146 | 0.05 |
| angular | 3 | 0.1 | cohort_ever_entered | 0.126568 | 0.05 |
| angular | 3 | 0.25 | continuous_ever_entered | 0.055887 | 0.05 |

The corrected audit records 3 orbit/tolerance instances where naive boundary-event counts differ from turning-point-bracketed crossings. These are numerical detection differences, not separate physical events or failed final tolerance gates. See the exact records in summary.json.

First entry after an earlier radial turn occurs in the retained data (the weighted fractions are in summary.csv). Thus the previous first-pass avoidance result cannot establish lasting avoidance. At both launch radii the sampled bar residence fractions are below the spherical control, but the actual formed population remains uncomputed.

## Meaning and next goal

This tests whether the ordinary nonspherical field can redistribute a prescribed captured population. It does not determine a self-consistent deposited density, central stability, full formation history, the energy funding of capture, or a galaxy rotation/lensing fit. A finite central visit fraction is neither a proof of a singular density nor a proof of stable support.

The next formation calculation needs moving populations in the nonspherical field with evolving capture and deposited gravity. Before treating encounter fractions as precise inputs, the failed angular convergence gates above require refinement. The source's velocity and angular distribution also remain physical assumptions. All nine research goals remain incomplete; the shared photon/time mechanism, sustained redshift, physical energy/momentum closure and genuinely unseen observational tests still need work.

To reproduce, run `audit.py coarse 64`, `audit.py fine 64`, and `audit.py fine 40` from this folder's scripts, using the existing hashed field caches and retained launch files. Then run `spherical_residence.py` and `export.py`. Set BLAS thread counts to one for concurrent orbit processes. There are no observational fits or withheld-data accesses in this checkpoint.
