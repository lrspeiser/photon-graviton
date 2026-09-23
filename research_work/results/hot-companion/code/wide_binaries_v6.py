#!/usr/bin/env python3
"""Round 6: what our law predicts for wide binary stars near the Sun.

    python wide_binaries_v6.py --output-dir ../run-wide-binaries-v6

Wide binaries (two stars 1,000-30,000 AU apart) feel each other's pull weakly: beyond about
5,000 AU it is below 10^-10 m/s^2, the level where galaxies stop following Newton. Gaia data are
contested: Chae (2023, 2024) reports about 1.4 times Newton's pull at the lowest accelerations,
MOND-like; Banik et al. (2024) find Newton and exclude MOND.

Our law. A binary is cold (its stars move at about 1 km/s: heat weight 10^-4), so
    h = g_N + exp(-|g_N|/g_d) sqrt(a (|g_N| + S_ext)) g_N/|g_N|,    del^2 Phi = -div h,
with g_N the Newtonian pull of the binary PLUS the Galaxy's. Near the Sun the Galaxy's own
Newtonian pull, g_e = 1.58 x 10^-10 m/s^2 (our law turned inside out for a 230 km/s rotation at
8.2 kpc), is comparable to g_d, so only about half the companion is released: the external field
holds it back, as in MOND, but more strongly.

Method (exact, no grid). Away from the stars, h - g_N = F(|g_N|) g^_N with F(g) = e^(-g/g_d)
sqrt(a (g + S_ext)), so the extra ("as if") density is analytic:
    rho_extra = -(1/4 pi G) div(F g^_N) = (M / 4 pi r^3) (F'(g) - F(g)/g) (1 - 3 (g^_N . r^)^2).
Averaged over the binary's orientation, Gauss's law gives the mean inward pull at separation r
as G (M + M_extra(<r)) / r^2, so the boost is gamma(r) = 1 + M_extra(<r)/M, and velocities scale
as sqrt(gamma). The same code gives QUMOND with the "simple" function (a0 = 1.2 x 10^-10, with its
own g_e = 1.33 x 10^-10 for the same 230 km/s) for comparison, and Newton (gamma = 1).
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np

G = 6.674e-11; MSUN = 1.989e30; AU = 1.496e11; KPC = 3.0857e19


def extra_mass_fraction(r_au, M, ge, F, dF, n_mu=400, r0_au=1.0):
    """M_extra(<r)/M for a point mass M (kg) in a uniform Newtonian field ge (m/s^2)."""
    mu, w = np.polynomial.legendre.leggauss(n_mu)           # mu = cos(angle between r^ and the external field)
    rr = np.geomspace(r0_au, r_au.max(), 4000) * AU
    gi = G * M / rr ** 2                                     # the binary's own Newtonian pull
    # g_N = ge z^ - gi r^  ;  |g_N|^2 = ge^2 + gi^2 - 2 ge gi mu
    g_par = ge * mu[None, :] - gi[:, None]                   # component along r^
    g_perp = ge * np.sqrt(1 - mu ** 2)[None, :]
    g = np.sqrt(g_par ** 2 + g_perp ** 2)
    mug = g_par / g                                          # g^_N . r^
    kern = (dF(g) - F(g) / g) * (1 - 3 * mug ** 2)
    ang = 0.5 * (kern * w[None, :]).sum(1)                   # angle average
    integrand = ang / rr                                     # dM_extra/M = ang dr / r
    cum = np.concatenate([[0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(rr))])
    return np.interp(r_au * AU, rr, cum)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    c = json.loads((Path(__file__).resolve().parent.parent / 'run-v3/results.json').read_text())['constants']
    a, gd = c['a_SI'], c['g_d_SI']
    a0 = 1.2e-10

    def ours(S_ext):
        F = lambda g: np.exp(-g / gd) * np.sqrt(a * (g + S_ext))
        dF = lambda g: F(g) * (-1 / gd + 0.5 / (g + S_ext))
        return F, dF

    nu = lambda y: 0.5 + np.sqrt(0.25 + 1 / y)
    F_m = lambda g: (nu(g / a0) - 1) * g
    def dF_m(g, h=1e-4):
        return (F_m(g * (1 + h)) - F_m(g * (1 - h))) / (2 * h * g)

    # the Galaxy's Newtonian pull at the Sun, from a 230 km/s rotation at 8.2 kpc
    g_obs = (230e3) ** 2 / (8.2 * KPC)
    def invert(total):
        lo, hi = 1e-13, g_obs
        for _ in range(200):
            m = 0.5 * (lo + hi)
            lo, hi = (m, hi) if total(m) < g_obs else (lo, m)
        return 0.5 * (lo + hi)
    ge_ours = invert(lambda gN: gN + np.exp(-gN / gd) * np.sqrt(a * gN))
    ge_mond = invert(lambda gN: gN * nu(gN / a0))
    print(f'the Galaxy at the Sun: observed {g_obs:.3e} m/s^2; Newtonian part {ge_ours:.3e} (our law), {ge_mond:.3e} (MOND)', flush=True)
    print(f'released share of the companion from the Galaxy alone: {np.exp(-ge_ours / gd):.2f}', flush=True)

    seps = np.array([500, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 30000], float)
    table = {}
    for M_sun in (1.0, 1.5, 2.0):
        M = M_sun * MSUN
        gi = G * M / (seps * AU) ** 2
        rows = {}
        for label, (F, dF), ge in (('ours', ours(0.0), ge_ours), ('ours, local disk heat S = 0.1 g_e', ours(0.1 * ge_ours), ge_ours),
                                   ('ours, g_e 20% lower', ours(0.0), 0.8 * ge_ours), ('ours, g_e 20% higher', ours(0.0), 1.2 * ge_ours),
                                   ('QUMOND simple', (F_m, dF_m), ge_mond), ('ours, no external field', ours(0.0), 1e-16)):
            gam = 1 + extra_mass_fraction(seps, M, ge, F, dF)
            rows[label] = gam.tolist()
        table[f'{M_sun}Msun'] = dict(separation_au=seps.tolist(), g_internal=gi.tolist(), gamma=rows)
        print(f'\nbinary of {M_sun} suns: boost of the pull (velocity boost is its square root)', flush=True)
        print('  separation (kAU)   ' + ' '.join(f'{s / 1000:6.1f}' for s in seps), flush=True)
        print('  own pull (1e-10)   ' + ' '.join(f'{v / 1e-10:6.2f}' for v in gi), flush=True)
        for k, v in rows.items():
            print(f'  {k:34s}' + ' '.join(f'{x:6.3f}' for x in v), flush=True)

    (out / 'wide_binaries_v6.json').write_text(json.dumps(dict(
        experiment='Wide binaries in the Galaxy near the Sun: orientation-averaged pull under our law vs MOND and Newton (round 6)',
        galaxy=dict(v_c_kms=230.0, R_kpc=8.2, g_obs=g_obs, g_N_ours=ge_ours, g_N_mond=ge_mond, released_share=float(np.exp(-ge_ours / gd))),
        constants=c, a0_mond=a0, results=table,
        observed=dict(chae='about 1.4 times the Newtonian pull below 1e-10 m/s^2 (Chae 2023, 2024)',
                      banik='Newtonian; MOND excluded (Banik et al. 2024)'),
        seconds=time.monotonic() - t0), indent=2, default=float) + '\n')


if __name__ == '__main__':
    main()
