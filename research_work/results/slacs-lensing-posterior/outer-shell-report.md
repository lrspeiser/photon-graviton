# Can edge deposits supply lensing without extra interior pull?

An outer spherical shell supplies a precise version of the user's edge-deposition idea. It adds light bending while its direct Newtonian force is zero inside. Applied as a positive addition to the current six stellar-conditioned lens models, it can fill three bending deficits but cannot repair the three median models already bending too strongly. Shell formation, support and a full stellar refit are not supplied by this diagnostic.

## Formula and provenance

**Hypothetical deposition geometry:** a thin spherical shell of effective gravitating mass M_s at radius R, with the usual equal-potential weak-field light response. This is not a derived companion distribution. Newton's shell theorem gives g_s=0 at r<R and g_s=GM_s/r^2 outside. The interior potential is constant, -GM_s/R; zero gradient is not zero potential or zero clock effect relative to an exterior observer.

For projected impact radius b<R, ordinary projection of the two spherical caps gives

    M_projected(<b)=M_s[1-sqrt(1-b^2/R^2)]
    reduced_deflection(b)=(D_ls/D_s) 4G M_projected(<b)/(c^2 b).

These are known projection and weak-lensing formulas, not new laws; the lensing formalism is reviewed by [Narayan and Bartelmann](https://arxiv.org/abs/astro-ph/9606001). The same standard light response is retained for the original and added components. No independently fitted lens multiplier is introduced.

For b much smaller than R,

    M_projected(<b) approximately M_s b^2/(2 R^2)
    Sigma_projected(0)=M_s/(2 pi R^2)
    kappa_shell approximately Sigma_projected(0)/Sigma_critical
    Sigma_critical=c^2/[4 pi G D_l (D_ls/D_s)].

The shell behaves approximately like a uniform projected mass sheet in the central region. Its full b-dependent expression, rather than this approximation, is used below. A lens-equation residual divided by the catalog angle is an effective convergence at that angle, not a directly measured external mass sheet.

## Requirements from the six existing lens summaries

Keep the existing empirical-extra stellar mass posteriors, light profiles and conditional distance geometry. Set R=20a using the previous outer scale, rather than selecting a radius separately to improve each result. The additions below force the old spherical deflection plus the shell deflection to equal the published SIE angle. Thus they are inferred requirements, not independent predictions or measurements of shell masses.

| Galaxy | Shell radius, kpc | Required effective convergence at catalog angle | Required signed shell mass, trillion solar masses |
|---|---:|---:|---:|
| J0037-0942 | 83.92 | -0.0309 | -4.98 |
| J1112+0826 | 77.70 | +0.2043 | +26.22 |
| J1204+0358 | 48.02 | -0.0582 | -3.32 |
| J1402+6321 | 107.86 | -0.0475 | -14.66 |
| J1621+3931 | 100.60 | +0.1452 | +32.57 |
| J1630+4520 | 93.15 | +0.2246 | +36.80 |

Negative entries mean a positive additional shell cannot perform that correction. They are not proposals for negative energy. The table uses median stellar masses; outer-shell-results.json also preserves the requirements at the endpoints of the previous central 95% mass ranges. These are conditional on exact catalog lens summaries and simplified stellar structure. No full lensing uncertainty or rejection probability is inferred.

The large positive mass requirements are geometrical: only a small fraction of a large-radius shell lies inside the projected central lens cylinder. They are not by themselves a photon-supply exclusion. No universe age, radius, accumulated radiation history or enhanced gravity-per-energy factor is fixed here. With standard mass-energy coupling, those effective masses would require corresponding energy, along with formation/support accounting.

## Connect collecting area to this lensing requirement

Under the earlier conditional isotropic-bath formula, a fixed shell radius and interception fraction f_cap(t) give

    P_captured(t)=pi R^2 c u_c(t) f_cap(t)
    M_s=pi R^2/c integral u_c(t) f_cap(t) dt
    kappa_shell approximately [1/(2 c Sigma_critical)] integral u_c(t) f_cap(t) dt.

The capture relation is established isotropic-flux accounting; retaining all intercepted energy as mass on this shell is an additional assumption. **R cancels from central convergence at fixed integrated bath exposure and capture fraction.** A bigger shell catches more energy, but that same area growth spreads its mass. This cancellation is limited to b<<R, fixed R during accumulation, the stated shell geometry and mass-energy response. It does not mean cluster size never affects capture, outer lensing, or the capture fraction.

No observed radiation exposure is fitted here. The required signed M_s/(pi R^2), saved per mass quantile, is the conditional target for the integrated companion bath divided by c; three median targets are negative and cannot be supplied by positive capture on top of the fixed model.

## Limits and next physical requirement

An interior zero shell force does not guarantee the complete observed stellar velocity distribution is unchanged. Stars on orbits extending outside R, outer boundary conditions and line-of-sight projection can respond. The old mass posterior is retained to diagnose the direct lensing addition, not certified after a full shell-plus-stars dynamical solution. A static thin shell also requires a support or orbit model and cannot simply be assumed to sit at its capture radius forever.

A positive outer deposit is a possible geometric contribution for individual deficits. It is not a universal additive repair to all six fixed models. Any redistribution that also reduces inner deposited gravity must be recalculated for both stellar motions and lensing using one source profile. outer-shell.py verifies the projected-mass inversion and preserves all signed requirements. No new fit to radiation history, final holdout or completed companion source is claimed. All six objectives remain open.
