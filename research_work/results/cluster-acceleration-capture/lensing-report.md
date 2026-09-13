# Central shear signature of the candidate capture profile

## Result

The acceleration-squared capture candidate has a specific central lensing consequence in the fixed, spherical comparison: the deposited component alone produces negative tangential shear near the center. In the usual weak-lensing convention this is a radial alignment tendency. It does not mean negative mass, outward gravitational acceleration or repulsive gravity. Ordinary central matter contributes positive tangential shear and can hide the effect.

This is a forward signature of the previously calculated density, not a fitted change introduced to match a lens. No observed cluster is compared here. Other hollow mass distributions can produce the same feature, so this is not unique evidence for photon conversion.

## Known projection identities

Let q(r) be the already calculated deposited power per volume, with a common positive accumulation factor converting its shape into stored density. For this smooth profile q(r) is proportional to r^2 near r=0 and falls sufficiently fast at infinity. Standard line-of-sight projection gives

    S(b)=2 integral_0^infinity q(sqrt(b^2+z^2)) dz
        =S(0)+B b^2+o(b^2)
    B=integral_0^infinity q'(r)/r dr
     =integral_0^infinity q(r)/r^2 dr > 0.

The integration-by-parts boundary terms vanish because q(r)/r tends to zero at both ends. Thus projected density has a local minimum at the center, despite remaining positive there from foreground/background deposits. The circular mean is Sbar(<b)=S(0)+B b^2/2+o(b^2), so

    DeltaS=Sbar-S=-B b^2/2+o(b^2) < 0.

In the conditional equal-potential weak-field response, gamma_t=DeltaSigma/Sigma_crit, so DeltaS fixes the deposited component's shear sign. Actual reduced shear is gamma/(1-kappa) and requires total convergence and geometry. The result is not a prediction of observed image shapes without those inputs.

All projection/shear equations are established lensing mathematics, using the same conditional response as the previous motion/lensing reports. The capture profile remains a phenomenological postulate, not a new first-principles derivation.

## Ordinary matter can change the total sign

For the unchanged Plummer ordinary-matter source, normalized projected density is Sigma_b=M_b/[pi a^2](1+b^2/a^2)^-2. Near the center its excess surface density is positive, M_b b^2/(pi a^4). With dimensionless q, total deposited-power normalization s=sigma_eff/a^2, and stored deposit mass M_d following that shape, total central shear changes sign when

    M_d/M_b > 2s/(pi B).

| Capture strength A | Critical M_d/M_b for negative total central tangential shear |
|---|---:|
| 0.1 | 8.106 |
| 1 | 9.153 |
| 10 | 36.405 |
| 100 | 4471.456 |

These are shape thresholds within an imposed frozen-well, stationary-density comparison, not self-consistent massive-halo solutions. Especially at large ratios, deposit backreaction cannot physically be presumed negligible. A different baryonic core, central black hole, nonspherical structure, accretion history or deposit redistribution changes the comparison. The table is not an observational exclusion threshold for every cluster.

At A=100, the deposited component remains negative at the sampled radius 3a and becomes positive by 10a. At smaller A it is negative at 0.1a and positive by a. These statements bracket behavior only at the evaluated samples; no precise zero-crossing is claimed.

## Verification and implications

Finite differences of the projected profile at b=0.001a agree with the analytic positive-curvature coefficient within 4.50e-6 relative across the four strengths. Projected aperture integrals computed from spherical shell geometry independently agree with direct image averaging within the declared 1e-6 relative-to-surface scale tolerance. All four cases have the predicted negative deposited DeltaSigma at 0.1a. Complete normalized density, excess-density and enclosed-mass samples are retained in lensing-results.json.

This takes the candidate beyond an adjustable collecting area: its spatial capture rule predicts an inner shear contribution that observations could challenge. The next physical step is a self-consistent or perturbatively controlled treatment of deposit feedback/support, followed by measured baryonic profiles and calibrated lensing/kinematic data under a fixed common coefficient. Total-field observations must be compared with total-field predictions. No observational holdouts were opened and all six objectives remain open.

Run `python research_work/results/cluster-acceleration-capture/lensing.py` after run.py.
