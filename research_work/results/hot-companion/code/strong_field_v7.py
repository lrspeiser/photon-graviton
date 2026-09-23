#!/usr/bin/env python3
"""Round 7: the precision tests of gravity (Solar System, Galactic Centre, binary pulsars).

    python strong_field_v7.py --output-dir ../run-strong-field-v7

These are the places where gravity is measured most precisely. Any new law must be silent there,
or it is ruled out. Our law's extra pull carries the release factor exp(-|g_N|/g_d), with
g_d = 2.26 x 10^-10 m/s^2, so it switches off wherever the ordinary pull is strong. This script
puts numbers on how silent it is:

1. Planets (Mercury to Neptune, and Sedna's far point): the extra pull at each.
2. Cassini's test of the Galaxy's 'external field' inside the Solar System. In laws like ours
   and MOND, the Galaxy's pull on the Sun reshapes the Sun's own extra pull far out, and that
   reshaping leaks inward as a tiny tidal-like distortion Q2 (Milgrom 2009; Blanchet & Novak
   2011). Hees et al. (2014) measured Q2 = (3 +- 3) x 10^-27 s^-2 from nine years of Cassini
   radio tracking of Saturn; that excludes MOND with its common interpolating functions.
   We compute Q2 exactly from the 'as if' density of the Sun in the Galaxy's field:
       Q2 = 3 G int rho_ph P2(cos theta) / r^3 dV   (theta from the Galactic Centre direction;
   Hees et al.'s convention, which gives Q2 = 3 G M(D)/D^3 for an ordinary tidal field).
3. The star S2 at the Galactic Centre (GRAVITY 2020: precession 1.10 +- 0.19 of GR's).
4. Binary pulsars: orbital decay (Hulse-Taylor, the Double Pulsar) and our law's only possible
   new effect there, the companion's steady mass loss 2.3 x 10^-15 per year.
5. Light bending at the Sun (VLBI, Cassini).
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np

G = 6.674e-11; MSUN = 1.989e30; AU = 1.496e11; C = 2.99792458e8; YR = 3.156e7; KPC = 3.0857e19


def q2_point_mass(M, ge, H, F, dF, rmin_au=10.0, rmax_au=1e8, n_r=6000, n_mu=800):
    """Q2 (Hees et al. convention) from the phantom density of a point mass M in a uniform
    Newtonian field ge (toward the Galactic Centre, along e) with a uniform hot pull H along e:
        h - g_N = D(|g_N|) (g_N + g_hot),   D = F/(|g_N| + H)
        rho_ph  = (M / 4 pi r^3) D'(g) [ w.g^ - 3 (w.r^)(g^.r^) ],   w = g_N + g_hot."""
    mu, wmu = np.polynomial.legendre.leggauss(n_mu)          # mu = r^ . e
    r = np.geomspace(rmin_au, rmax_au, n_r) * AU
    gi = G * M / r ** 2
    g_par = ge * mu[None, :] - gi[:, None]                   # g_N . r^
    g_perp = ge * np.sqrt(1 - mu ** 2)[None, :]
    g = np.sqrt(g_par ** 2 + g_perp ** 2)
    w_par = g_par + H * mu[None, :]; w_perp = g_perp + H * np.sqrt(1 - mu ** 2)[None, :]
    w_dot_gh = (w_par * g_par + w_perp * g_perp) / g
    gh_r = g_par / g
    D = lambda x: F(x) / (x + H)
    dD = lambda x: (dF(x) * (x + H) - F(x)) / (x + H) ** 2
    ang = dD(g) * (w_dot_gh - 3 * w_par * gh_r)
    P2 = 0.5 * (3 * mu ** 2 - 1)
    inner = (ang * P2[None, :] * wmu[None, :]).sum(1)      # int dmu
    # I2 = int rho P2 / r^3 dV = (M/2) int dr / r^4 int dmu [...]
    I2 = 0.5 * M * np.trapezoid(inner / r ** 4, r)
    # where it comes from: the share of I2 from inside each radius
    cum = 0.5 * M * np.concatenate([[0], np.cumsum(0.5 * (inner[1:] / r[1:] ** 4 + inner[:-1] / r[:-1] ** 4) * np.diff(r))])
    return 3 * G * I2, r / AU, cum / I2 if I2 != 0 else cum


def q2_gradual(M, ge, L_au, a, gd, rmin_au=10.0, rmax_au=1e9, n_r=8000, n_mu=800):
    """Q2 and the wide-binary boost when the companion is released only gradually along its path:
    F -> F0(|g_N|) R(r), R = 1 - exp(-r/L), r = distance travelled from the Sun (a candidate fix)."""
    mu, w = np.polynomial.legendre.leggauss(n_mu)
    r = np.geomspace(rmin_au, rmax_au, n_r) * AU
    gi = G * M / r ** 2
    gp = ge * mu[None, :] - gi[:, None]; gq = ge * np.sqrt(1 - mu ** 2)[None, :]
    g = np.sqrt(gp ** 2 + gq ** 2); cth = gp / g
    F0 = np.exp(-g / gd) * np.sqrt(a * g); dF0 = F0 * (-1 / gd + 0.5 / g)
    if L_au > 0:
        L = L_au * AU; R = -np.expm1(-r / L); dR = np.exp(-r / L) / L
    else:
        R = np.ones_like(r); dR = np.zeros_like(r)
    rho = R[:, None] * (M / (4 * np.pi * r[:, None] ** 3)) * (dF0 - F0 / g) * (1 - 3 * cth ** 2) - F0 * dR[:, None] * cth / (4 * np.pi * G)
    P2 = 0.5 * (3 * mu ** 2 - 1)
    I2 = np.trapezoid(2 * np.pi * (rho * P2[None, :] * w[None, :]).sum(1) / r, r)
    iso = 2 * np.pi * (rho * w[None, :]).sum(1) * r ** 2
    Mph = np.concatenate([[0], np.cumsum(0.5 * (iso[1:] + iso[:-1]) * np.diff(r))])
    boost = {k: float(1 + np.interp(k * AU, r, Mph) / M) for k in (3000, 7000, 20000)}
    return 3 * G * I2, boost


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    here = Path(__file__).resolve().parent
    c = json.loads((here.parent / 'run-v3/results.json').read_text())['constants']
    a, gd = c['a_SI'], c['g_d_SI']; a0 = 1.2e-10
    ell = json.loads((here.parent / 'run-v3/results.json').read_text())['energy']['emission_W_per_kg']
    res = {}

    # 1. planets
    planets = dict(Mercury=0.387, Venus=0.723, Earth=1.0, Mars=1.524, Jupiter=5.203, Saturn=9.537, Uranus=19.19, Neptune=30.07,
                   **{'Sedna (far point)': 937.0, 'Oort cloud inner edge': 2000.0})
    rows = {}
    for n, d in planets.items():
        gN = G * MSUN / (d * AU) ** 2
        x = gN / gd
        extra = np.exp(-x) * np.sqrt(a * gN)
        rows[n] = dict(distance_au=d, g_N=gN, g_N_over_g_d=x, log10_release_factor=float(-x / np.log(10)),
                       extra_pull_SI=float(extra), extra_over_newton=float(extra / gN))
    res['planets'] = rows
    print('Planets: the extra pull is exp(-g_N/g_d) sqrt(a g_N):', flush=True)
    for n, r in rows.items():
        print(f"  {n:22s} {r['distance_au']:8.1f} AU  g_N/g_d = {r['g_N_over_g_d']:10.3g}  release factor 10^{r['log10_release_factor']:.4g}  extra/Newton = {r['extra_over_newton']:.3g}", flush=True)

    # 2. Cassini Q2
    def ours(S):
        F = lambda g: np.exp(-g / gd) * np.sqrt(a * (g + S))
        dF = lambda g: F(g) * (-1 / gd + 0.5 / (g + S))
        return F, dF
    nu = lambda y: 0.5 + np.sqrt(0.25 + 1 / y)
    Fm = lambda g: (nu(g / a0) - 1) * g
    dFm = lambda g: (Fm(g * (1 + 1e-5)) - Fm(g * (1 - 1e-5))) / (2e-5 * g)
    # the Galaxy's Newtonian pull at the Sun under each law, for an observed 1.9e-10 (Hees et al.'s value) and
    # for our Milky Way model (round 7): g_N from McMillan (2017) baryons, with the Galaxy's heat S and hot pull H
    mwj = json.loads((here.parent / 'run-milky-way-v7/milky_way_v7.json').read_text())['sun'] if (here.parent / 'run-milky-way-v7/milky_way_v7.json').exists() else None
    def invert(total, gobs):
        lo, hi = 1e-13, gobs
        for _ in range(200):
            m = 0.5 * (lo + hi)
            lo, hi = (m, hi) if total(m) < gobs else (lo, m)
        return 0.5 * (lo + hi)
    q2 = {}
    for gobs in (1.9e-10, 2.09e-10):
        ge_o = invert(lambda x: x + np.exp(-x / gd) * np.sqrt(a * x), gobs)
        ge_m = invert(lambda x: x * nu(x / a0), gobs)
        Qo, rr, cumo = q2_point_mass(MSUN, ge_o, 0.0, *ours(0.0))
        Qm, _, cumm = q2_point_mass(MSUN, ge_m, 0.0, Fm, dFm)
        half_o = float(np.interp(0.5, cumo, rr)); half_m = float(np.interp(0.5, cumm, rr))
        q2[f'g_obs {gobs:.3g}'] = dict(ours=dict(g_N=ge_o, Q2=Qo, half_from_inside_au=half_o), mond_simple=dict(g_N=ge_m, Q2=Qm, half_from_inside_au=half_m))
        print(f'Cassini Q2 for a Galactic pull of {gobs:.3g} m/s^2 at the Sun: ours {Qo:.2e} s^-2 (Newtonian part {ge_o:.3g}; half of it from beyond {half_o:.0f} AU); '
              f'MOND simple {Qm:.2e} (g_N {ge_m:.3g}; half from beyond {half_m:.0f} AU); measured (3 +- 3)e-27', flush=True)
    if mwj:
        for label, ge, S, H in (('ours, Milky Way model (McMillan 2017 baryons) with the Galaxy heat', mwj['g_N_SI'], mwj['S_SI'], mwj['g_hot_SI']),
                                ('ours, Milky Way model, heat off', mwj['g_N_SI'], 0.0, 0.0)):
            Qo, rr, cumo = q2_point_mass(MSUN, ge, H, *ours(S))
            q2[label] = dict(Q2=Qo, g_N=ge, S=S, H=H, half_from_inside_au=float(np.interp(0.5, cumo, rr)))
            print(f'  {label}: Q2 = {Qo:.2e} s^-2', flush=True)
    # what would have to change: gradual release along the companion's path (release length L)
    grad = {}
    for ge in (1.24e-10, 1.585e-10):
        rows_g = []
        for L in (0.0, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6):
            Q, b = q2_gradual(MSUN, ge, L, a, gd)
            rows_g.append(dict(L_au=L, L_pc=L * AU / (KPC / 1000), travel_time_yr=L * AU / (c['u_kms'] * 1e3) / YR, Q2=Q, wide_binary_boost=b))
        grad[f'g_N {ge:.3g}'] = rows_g
        print(f'Gradual release, Galactic pull {ge:.3g}: ' + '; '.join(f"L {r_['L_au']:.0e} AU ({r_['travel_time_yr']:.0f} yr): Q2 {r_['Q2']:.1e}, binary boost 7k/20k AU {r_['wide_binary_boost'][7000]:.3f}/{r_['wide_binary_boost'][20000]:.3f}" for r_ in rows_g), flush=True)
    q2['gradual_release'] = grad
    res['cassini_Q2'] = dict(measured=dict(value=3e-27, error=3e-27, source='Hees, Folkner, Jacobson & Park 2014, PRD 89, 102002 (arXiv:1402.6950), eq. 12'),
                             mond_range_quoted=dict(low=2.1e-27, high=4.1e-26, note='Hees et al. eq. 7, from Milgrom 2009 for g_e = 1.9e-10, a0 = 1.2e-10'),
                             predictions=q2)

    # 3. S2 at the Galactic Centre
    M_bh = 4.3e6 * MSUN
    peri = 120 * AU; apo = 1950 * AU
    gp, ga = G * M_bh / peri ** 2, G * M_bh / apo ** 2
    res['S2'] = dict(g_peri=gp, g_apo=ga, release_log10_peri=float(-gp / gd / np.log(10)), release_log10_apo=float(-ga / gd / np.log(10)),
                     radius_where_release_reaches_1pct_pc=float(np.sqrt(G * M_bh / (gd * np.log(100))) / (KPC / 1000)),
                     prediction='Schwarzschild precession exactly GR (f_SP = 1)', measured='f_SP = 1.10 +- 0.19 (GRAVITY Collaboration 2020)')
    print(f"S2: pull {gp:.2e} (pericentre) to {ga:.2e} m/s^2 (apocentre); release factor 10^{res['S2']['release_log10_peri']:.3g} to 10^{res['S2']['release_log10_apo']:.3g}: pure GR", flush=True)

    # 4. binary pulsars: mass loss by the companion emission
    mdot = ell / C ** 2                                      # per second, fraction of mass
    for name, Pb, Pbdot, err in (('Hulse-Taylor PSR B1913+16', 27906.98, -2.423e-12, 0.001e-12),
                                 ('Double Pulsar PSR J0737-3039A/B', 8834.53, -1.247920e-12, 0.000078e-12)):
        extra = 2 * mdot * Pb                                # isotropic mass loss widens the orbit: Pbdot/Pb = -2 Mdot/M
        res.setdefault('pulsars', {})[name] = dict(Pb_s=Pb, Pbdot_obs=Pbdot, Pbdot_err=err, extra_from_mass_loss=extra,
                                                    extra_over_error=extra / err, orbit_pull_note='g >> g_d: release factor zero, orbit decay exactly GR')
        print(f'{name}: our only new effect (mass loss) changes Pbdot by {extra:.1e}, {extra / err:.3g} of the measurement error', flush=True)

    # 5. light bending at the Sun
    gsun = G * MSUN / (6.957e8) ** 2
    res['light_bending'] = dict(g_limb=gsun, release_log10=float(-gsun / gd / np.log(10)), deflection_arcsec=4 * G * MSUN / (C ** 2 * 6.957e8) * 206265,
                                prediction='GR exactly (gamma = 1)', measured='gamma - 1 = (2.1 +- 2.3) x 10^-5 (Cassini, Bertotti et al. 2003)')
    print(f"Light bending at the Sun's limb: {res['light_bending']['deflection_arcsec']:.4f} arcsec, release factor 10^{res['light_bending']['release_log10']:.3g}", flush=True)
    res['constants'] = c; res['seconds'] = time.monotonic() - t0
    (out / 'strong_field_v7.json').write_text(json.dumps(res, indent=2, default=float) + '\n')


if __name__ == '__main__':
    main()
