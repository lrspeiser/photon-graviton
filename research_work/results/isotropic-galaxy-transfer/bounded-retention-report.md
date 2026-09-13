# Bounded radiation retention improves frozen galaxy predictions

Changing radiation from an opacity multiplier to a bounded retained-energy fraction improves both reported error measures on both comparison groups. The well-depth alternative does not offer the same improvement. This is conditional progress in galaxy extra-gravity fitting, not a universal photon-to-graviton law or independent discovery. The data partitions have been repeatedly exposed.

## Revised formula and constants

**Postulated retention probability using known logistic mathematics:**

    X=(L_3.6/10^9 L_sun)/(R_disk/kpc)^2,
    eta(X)=X^q/(1+X^q), q=0.3399907823.

The environmental variable controls storage after capture, not interception opacity. In plain language, brighter galaxies per unit disk area retain a larger share of incoming captured energy, saturating below 100%. No physical radiation-driven binding mechanism has yet been derived.

Keep the known transport integral J for the proposed opacity profile, and set

    kappa(r)=k0/[1+(r/a)^2]^2,
    rho_deposit(r)=C_source*eta(X)*J(r)/[1+(r/a)^2]^2,
    a=s*R_disk.

Training-only shared fit:

| Constant | Value |
|---|---:|
| k0 | 0.2059105265 /kpc |
| s | 2.752675347 |
| C_source before retention | 95,668,119.60 Msun/kpc^3 |
| q | 0.3399907823 |

The stored JSON parameter C_Msun_kpc3 is C_half=47,834,059.80; the implementation multiplies it by 2*eta to preserve the q=0 nested normalization. Do not treat C_half as C_source. Under the inherited source law, C_source=(k0/c)*integral u_external dt. Retention eta ranges 0.20149–0.74971 over these galaxies. The unretained energy must escape in a specified channel; the single-pass model accounts for it as released, not as destroyed. Its detailed propagation remains missing.

## Observed-versus-predicted errors

Fit 89 training galaxies, freeze parameters, then predict 29 validation and 31 test galaxies. All use identical observed radii/speeds, ordinary-matter assumptions, distances and equal-galaxy weighting.

| Formula | Validation RMS km/s | Test RMS km/s | Validation log RMS | Test log RMS |
|---|---:|---:|---:|---:|
| Original interception | 33.18419 | 23.99377 | 0.119980 | 0.098545 |
| Prior radiation-dependent opacity | 35.34063 | 24.50764 | 0.115922 | 0.094572 |
| **Bounded radiation retention** | **32.52002** | **23.62155** | **0.115383** | **0.090991** |
| Bounded ordinary-depth retention | 33.25590 | 24.13970 | 0.119878 | 0.099150 |

Radiation-retention training RMS is 29.01857 km/s and log RMS 0.138488. It improves absolute transfer errors modestly and fractional transfer errors more clearly. These are descriptive differences after multiple candidate trials, not significance estimates or blind confirmation. It remains behind the earlier empirical force baseline in absolute errors; physical derivation and broader tests remain necessary.

The depth proxy is v_baryon^2(R_disk)/(1000 (km/s)^2), not the full potential. Its fitted q=-0.026342 is near zero, with eta=0.47158–0.52183. The sign does not support an increasing retention preference for this particular depth proxy. It does not exclude other binding physics.

## Decision and remaining scope

Promote bounded radiation retention as a candidate for the next cross-observable comparison, preserving the original reference and prior failed/mixed candidates. Do not yet replace the whole theory. The same deposits must next predict lensing and stellar motion under a consistent luminosity/morphology mapping. No per-lens correction may be fitted to rescue the prediction.

This galaxy extra-gravity calculation does not reproduce the entire Sun/Earth/Moon masses. Full-mass and extra-mass targets require an explicit common interpretation before transfer. Graviton frequency/count, independently bound states, storage history and a causal energy-conserving redshift/timing law remain open. The fixed X=1 pivot is a specified scale, not a measured universal constant.

## Verification

Both three-start fits converge without parameter bounds. q=0 reproduces the original predictions before fitting. Finer integration changes group RMS by <=0.00113 km/s, far smaller than the reported differences. Retention probabilities are checked within [0,1]. The comparison script verifies that every candidate used identical observed speeds, radii and sample roles. All radial predictions remain saved. This validates the computation, not the microscopic mechanism.

Reproduce:

    python research_work/results/isotropic-galaxy-transfer/run.py --bounded-retention
    python research_work/results/isotropic-galaxy-transfer/run.py --bounded-retention --depth-retention
    python research_work/results/isotropic-galaxy-transfer/bounded-retention-compare.py
