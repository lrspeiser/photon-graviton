"""The precision tests of gravity: planets, the star S2, binary pulsars, light bending, Cassini's
limit on a Galactic distortion of the Sun's field (Q2), and the wide-binary prediction.

Q2 and the wide binaries use one exact formula for the Sun (or a binary) in the Galaxy's pull,
with both round-7 candidate amendments built in. Outside the star, h - g_N = R(r) D(|g_N|) w with
w = g_N + g_hot, D = F/(|g_N| + |g_hot|), F = exp(-g/g_d) sqrt(a (g + S)), R = 1 - exp(-r/L):
    rho_ph = R (M / 4 pi r^3) D'(g) [w.g^ - 3 (w.r^)(g^.r^)]  -  R'(r) D (w.r^) / (4 pi G)
Q2 = 3 G int rho_ph P2 / r^3 dV (Hees et al. 2014 convention); the orientation-averaged boost of
a binary's pull at separation r is 1 + M_ph(<r)/M. The Galaxy's Newtonian pull at the Sun g_e
comes from turning the law inside out for 230 km/s at 8.2 kpc (round 6); the external hold c
scales the Galaxy's pull and heat as the subsystem feels them.
"""
from __future__ import annotations
import numpy as np
from checks import make, z_check, at_most, range_check

GROUP = 'precision'
G, MSUN, AU, CL, YR, KPC = 6.674e-11, 1.989e30, 1.496e11, 2.99792458e8, 3.156e7, 3.0857e19


def sun_in_galaxy(M, ge, S, H, L_au, a, gd, rmin_au=10.0, rmax_au=1e9, n_r=8000, n_mu=800, seps_au=(3000, 7000, 20000)):
    mu, w = np.polynomial.legendre.leggauss(n_mu)
    r = np.geomspace(rmin_au, rmax_au, n_r) * AU
    gi = G * M / r ** 2
    st = np.sqrt(1 - mu ** 2)
    gp = ge * mu[None, :] - gi[:, None]; gq = ge * st[None, :]
    g = np.sqrt(gp ** 2 + gq ** 2)
    wp = gp + H * mu[None, :]; wq = gq + H * st[None, :]
    w_dot_gh = (wp * gp + wq * gq) / g
    F = np.exp(-g / gd) * np.sqrt(a * (g + S)); dF = F * (-1 / gd + 0.5 / (g + S))
    D = F / (g + H); dD = (dF * (g + H) - F) / (g + H) ** 2
    rho0 = (M / (4 * np.pi * r[:, None] ** 3)) * dD * (w_dot_gh - 3 * wp * gp / g)
    if L_au > 0:
        L = L_au * AU; R = -np.expm1(-r / L); dR = np.exp(-r / L) / L
    else:
        R = np.ones_like(r); dR = np.zeros_like(r)
    rho = R[:, None] * rho0 - dR[:, None] * D * wp / (4 * np.pi * G)
    P2 = 0.5 * (3 * mu ** 2 - 1)
    Q2 = 3 * G * np.trapezoid(2 * np.pi * (rho * P2[None, :] * w[None, :]).sum(1) / r, r)
    iso = 2 * np.pi * (rho * w[None, :]).sum(1) * r ** 2
    Mph = np.concatenate([[0], np.cumsum(0.5 * (iso[1:] + iso[:-1]) * np.diff(r))])
    boost = {k: float(1 + np.interp(k * AU, r, Mph) / M) for k in seps_au}
    return float(Q2), boost


def galactic_pull(law, v=230.0, R_kpc=8.2):
    """Newtonian part of the Galaxy's pull at the Sun for an observed rotation v (cold, heat left out)."""
    a, gd = law['a_SI'], law['g_d_SI']
    gobs = (v * 1e3) ** 2 / (R_kpc * KPC)
    lo, hi = 1e-13, gobs
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if m + np.exp(-m / gd) * np.sqrt(a * m) < gobs else (lo, m)
    return 0.5 * (lo + hi), gobs


def run(law, ctx):
    a, gd, c, L = law['a_SI'], law['g_d_SI'], law['external_hold'], law['release_length_au']
    out = []
    worst = 0.0
    for d in (0.387, 0.723, 1.0, 1.524, 5.203, 9.537, 19.19, 30.07):
        gN = G * MSUN / (d * AU) ** 2
        worst = max(worst, float(np.exp(-gN / gd) * np.sqrt(a * gN) / gN))
    out.append(make('precision.planets', GROUP, 'planets Mercury-Neptune: largest extra pull relative to Newton', worst,
                    crit=at_most(worst, 1e-12), target='< 1e-12 (planetary ephemerides)'))
    ga = G * 4.3e6 * MSUN / (1950 * AU) ** 2
    s2 = float(np.exp(-ga / gd) * np.sqrt(a * ga) / ga)
    out.append(make('precision.S2', GROUP, 'star S2 at the Galactic Centre (far point): extra pull relative to Newton', s2,
                    crit=at_most(s2, 1e-6), target='GR precession 1.10 +- 0.19 (GRAVITY 2020)'))
    ell = a * law['u_kms'] * 1e3 / 2
    extra = 2 * ell / CL ** 2 * 8834.53
    rel = extra / 0.000078e-12
    out.append(make('precision.double_pulsar', GROUP, 'Double Pulsar: orbit change from the companion emission, in measurement errors', rel,
                    crit=at_most(rel, 0.5, 1.0), target='< 0.5 of the error on the orbital decay (Kramer et al. 2021)'))
    gl = G * MSUN / 6.957e8 ** 2
    lb = float(np.exp(-gl / gd) * np.sqrt(a * gl) / gl)
    out.append(make('precision.light_bending', GROUP, 'light bending at the Sun\'s limb: extra relative to GR', lb,
                    crit=at_most(lb, 1e-6), target='gamma - 1 = (2.1 +- 2.3) x 10^-5 (Cassini)'))

    ge, gobs = galactic_pull(law)
    ctx.log(f'Cassini Q2 (g_e {ge:.3e}, hold x{c:g}, release length {L:g} AU)')
    Q2, boost = sun_in_galaxy(MSUN, c * ge, 0.0, 0.0, L, a, gd)
    # round 18: the 2026 re-estimate with the DE440 data (40% tighter than Hees et al. 2014's (3 +- 3) x 10^-27)
    out.append(make('precision.cassini_q2', GROUP, 'Cassini: Galactic distortion of the Sun\'s field, Q2', Q2,
                    crit=z_check(Q2, 1.6e-27, 1.8e-27), unit='1/s^2', target='(1.6 +- 1.8) x 10^-27 (Park et al. 2026)',
                    detail=dict(g_e=ge, g_obs=gobs, z_2014=(Q2 - 3e-27) / 3e-27),
                    refs='Park, Hees, Famaey, Desmond & Durakovic 2026, PRD (arXiv:2602.17884); before: Hees et al. 2014, PRD 89, 102002'))
    sun = ctx.shared.get('mw_sun')
    if sun:
        Qm, _ = sun_in_galaxy(MSUN, c * sun['g_N_SI'], c * sun['S_SI'], c * sun['g_hot_SI'], L, a, gd)
        out.append(make('precision.cassini_q2_mw_model', GROUP, 'Cassini Q2 with the Milky Way model\'s pull and heat at the Sun', Qm,
                        unit='1/s^2', crit=('info', abs(Qm - 1.6e-27) / 1.8e-27)))
    # the two published analyses disagree (about 1.4: Chae 2023-24; 1.0: Banik et al. 2024), so anything between
    # them passes; a boost beyond both is excluded by both
    out.append(make('precision.wide_binaries', GROUP, 'wide binaries (1 Msun, 20,000 AU): boost of the pull', boost[20000],
                    target='between the two published analyses, 1.0 (Banik et al. 2024) and 1.4-1.5 (Chae 2023-24)',
                    crit=range_check(boost[20000], 1.0, 1.5, 0.1), detail=dict(boost_3000=boost[3000], boost_7000=boost[7000])))
    return out
