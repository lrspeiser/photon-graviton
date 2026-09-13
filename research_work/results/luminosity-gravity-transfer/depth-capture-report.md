# Shared depth-dependent capture fails the galaxy comparison

This run replaces the empirical extra-acceleration fit with an integrated, positive deposited-density profile. It fits two common parameters on 89 galaxies and predicts the other 29 and 31 without galaxy-specific normalization. The result is substantially worse than the previous empirical relations. It rules out this implementation as a satisfactory descriptive fit, not all possible companion capture.

## Postulate and derivation

Postulated capture coefficient:

    kappa(r)=kappa_max [W_b(r)/(W_b(r)+W0)]^4

W_b=-Phi_b with zero potential at infinity is a spherical binding-depth proxy constructed from the catalog ordinary-matter radial acceleration. This is not the full three-dimensional disk potential. The saturation law and fourth power are candidate assumptions previously explored in toy calculations; neither is a derived interaction or a claim of unique mathematics.

Known local transport bookkeeping, assuming a stationary optically thin isotropic bath:

    d rho_d/dt = kappa(r) u_c(t)/c
    rho_d(r) = C [W_b(r)/(W_b(r)+W0)]^4
    C = (kappa_max/c) integral u_c(t) dt
    M_d(<r) = 4 pi integral_0^r rho_d(s) s^2 ds
    v_pred^2(r) = v_b^2(r) + G M_d(<r)/r

The factor 1/c converts captured energy per volume per time, c*kappa*u_c, into mass density rate via E=mc^2. Exposure C is fitted, so no universe age, size or calculated stellar energy supply is asserted. At finite W0 and a Kepler outer potential, density falls as r^-4 and total mass converges. No fixed physical capture boundary is imposed. The formulas use known Newtonian gravity and transport identities plus the explicit capture postulate.

## Actual rotation results

| Group | Empirical baseline RMS km/s | Depth-capture RMS km/s | Baseline log10 RMS | Depth-capture log10 RMS |
|---|---:|---:|---:|---:|
| 89 training galaxies | 20.005 | 63.463 | 0.109349 | 0.176903 |
| 29 validation galaxies | 27.377 | 52.919 | 0.097429 | 0.147006 |
| 31 test galaxies | 17.200 | 45.236 | 0.079141 | 0.127253 |

The objective weights each galaxy equally in squared log10 speed residual. All three optimizer starts converge to C=4.16662e6 Msun/kpc^3 and W0=10 (km/s)^2, the declared lower bound. There is therefore no identified interior optimum for the depth scale; this result does not justify quoting that value as a measured physical constant.

The near-saturation direction makes rho_d approximately constant wherever W_b greatly exceeds W0. In the exactly uniform-density limit, known spherical gravity gives v_d^2=(4pi G C/3)r^2: the extra component keeps rising linearly in speed, rather than automatically creating a flat curve. That helps explain why a common exposure plus this simple depth response is inadequate. A positive density law alone does not guarantee the observed radial dependence.

The extrapolated Kepler radius where W_b=W0 ranges from 21.6 to 156,466 kpc, with median 3,319 kpc. These are consequences of the imposed outer continuation and fitted bound, not observed halo sizes or an assumed universe boundary. They show that the fit does not produce a well-confined galactic capture zone on its own.

## Physical limitation exposed by the fit

At 1,352 of 3,150 sampled radii, calculated deposited force exceeds the ordinary-matter force. Thus using only ordinary depth to set capture cannot be defended merely by saying deposits are a small perturbation. Force ratio is a warning, not a quantitative total-potential error estimate. A self-consistent depth model would let deposits alter subsequent capture, attenuate the incident bath where needed, and account for retention/support. This calculation includes none of those effects and cannot be presented as an energy-supply or stability solution.

The input luminosities/distances/mass templates retain the earlier provenance caveats. Published distances are stipulated for the fictional universe. All partitions have been exposed before; the evaluation is frozen-parameter transfer, not blind evidence. No uncertainty-marginalized significance is claimed.

The next substantive depth-capture candidate must address feedback and spatial support rather than adding one adjustable deposit amplitude per galaxy to hide this failure. External source histories are still required to calculate C instead of fitting it. All six primary goals remain open.

Artifacts: depth-capture-results.json contains parameters, all starts and galaxy scores; depth-capture-predictions.json contains the 3,150 observed/predicted speeds and force ratios. Input hashes are checked against the previous run. Reproduce with `python research_work/results/luminosity-gravity-transfer/depth-capture.py`.
