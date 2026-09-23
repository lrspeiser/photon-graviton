#!/usr/bin/env python3
"""T3.5 and the six SLACS lenses under hot-companion gravity, two bookkeeping conventions.

    python lenses_t35.py --output-dir ../run-lenses

  project  : the project's static-universe distances and its converted stellar masses
  standard : flat LCDM distances (H0=70, Om=0.3) and the published SLACS Chabrier masses
             (Auger et al. 2009; the release's own convention)
For each lens and law: the stellar-mass change needed to match the Einstein radius
(lensing) and, separately, to match the resolved velocity dispersions (kinematics).
Their difference is a direct slip test: light and matter feeling different pulls would
show up as different required masses.
Published IMF trend for comparison (Posacki, Cappellari & Treu 2015, arXiv:1407.5633):
  log alpha = 0.38 log(sigma_e/200 km/s) - 0.06,  alpha = (M*/L)_dyn / (M*/L)_Salpeter
  (derived WITH a dark-matter halo, which this law does not have).
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.linalg import solve_triangular

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / 'research_work/results/cross-prediction-response/code'))
import law as L
import model as M

A = 6.5797e-11 / L.KMS2_PER_KPC; LAM = 4.281; U = 874.1
SALP = 0.25                       # dex, Salpeter minus Chabrier (Auger et al. 2009 convention)


def lcdm(z, H0=70., Om=0.3):
    DH = L.C_KMS / H0
    Dc = lambda zz: DH * quad(lambda x: 1 / np.sqrt(Om * (1 + x) ** 3 + 1 - Om), 0, zz)[0]
    return Dc(z)


def patched_readers(convention):
    J = M.J
    orig = J.read_json
    obs = {x['Name']: x for x in orig(J.INPUT_NAMES[3])}
    def reader(rel):
        data = orig(rel)
        if convention == 'project':
            return data
        if rel == J.INPUT_NAMES[4]:
            out = []
            for g in data:
                o = obs.get(g['Name'])
                if o is None: out.append(g); continue
                zl, zs = o['zFG'], o['zBG']; Dcl, Dcs = lcdm(zl), lcdm(zs)
                g = dict(g); g['conditional_Dl_Mpc'] = Dcl / (1 + zl); g['conditional_Ds_Mpc'] = Dcs / (1 + zs)
                g['conditional_Dls_over_Ds'] = (Dcs - Dcl) / Dcs
                out.append(g)
            return out
        if rel == J.INPUT_NAMES[7]:
            out = []
            for m in data:
                m = dict(m)
                if m['imf'] == 'Chabrier' and m.get('published_log10_stellar_mass') is not None:
                    m['conditional_log10_stellar_mass'] = m['published_log10_stellar_mass']
                out.append(m)
            return out
        return data
    return orig, reader


def analyse(lens, law):
    r = lens.r; dm_star = lens.Mstar * np.diff(np.r_[0., lens.frac]); sigma = float(np.mean(lens.y))
    S_unit = L.G * (L.shell_weights(r, r) @ dm_star) / r ** 2
    gN_unit = L.G * lens.Mstar * lens.frac / r ** 2
    k = L.heat_weight(sigma, U) if law == 'hot companion' else 0.

    def g_on_r(dm):
        f = 10 ** dm
        if law == 'Newton': return f * gN_unit
        return L.total(f * gN_unit, f * k * S_unit, A, LAM)

    def theta(dm):
        gr = g_on_r(dm); gi = lambda x: np.exp(np.interp(np.log(x), np.log(r), np.log(gr)))
        def f(b):
            rr = b / np.cos(lens.lens_t)
            return lens.ratio * 4 / L.C_KMS ** 2 * np.dot(lens.lens_w, gi(rr) * rr) - b / lens.Dl
        grid = lens.Re * np.geomspace(1e-4, 1e3, 300); vals = np.array([f(b) for b in grid])
        ii = np.flatnonzero((vals[:-1] > 0) & (vals[1:] <= 0))
        return float(brentq(f, grid[ii[-1]], grid[ii[-1] + 1]) / lens.Dl * 206264.80624709636) if len(ii) else 0.

    def chi2(dm, beta):
        v = lens.moments(g_on_r(dm), {'beta': beta}, False); w = solve_triangular(lens.chol, v - lens.y, lower=True)
        return float(w @ w)

    dm_lens = brentq(lambda d: theta(d) - lens.theta, -1.5, 1.5, xtol=1e-7)
    fits = {}
    for beta in (0.0, -0.3, 0.3):
        rr = minimize_scalar(lambda d: chi2(d, beta), bounds=(-1.5, 1.5), method='bounded')
        fits[beta] = (float(rr.x), float(rr.fun))
    b_best = min(fits, key=lambda b: fits[b][1])
    return dict(sigma=sigma, dm_lens=float(dm_lens), dm_kin_isotropic=fits[0.0][0], chi2_kin_isotropic=fits[0.0][1],
                dm_kin_best=fits[b_best][0], beta_best=b_best, chi2_kin_best=fits[b_best][1],
                slip_gap_dex=float(dm_lens - fits[0.0][0]), mass_error=float(lens.mass_log_error))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--constants', type=Path, default=None,
                    help='results.json whose constants (a_code, lam, u_kms) replace the round-1 values')
    args = ap.parse_args(); out = args.output_dir
    if args.constants:
        global A, LAM, U
        c = json.loads(args.constants.read_text())['constants']
        A, LAM, U = c['a_code'], c['lam'], c['u_kms']
        print(f'constants from {args.constants}: a = {A * L.KMS2_PER_KPC:.4e} m/s^2, lambda = {LAM:.3f}, u = {U:.1f} km/s')
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True)
    res = {}
    for conv in ('project', 'standard'):
        orig, reader = patched_readers(conv); M.J.read_json = reader
        try:
            lenses = [M.make_lens(n) for n in M.LENSES]
        finally:
            M.J.read_json = orig
        res[conv] = {}
        for law in ('Newton', 'hot companion'):
            rows = []
            for lens in lenses:
                a = analyse(lens, law); a['name'] = lens.name; a['logMstar_input'] = float(np.log10(lens.Mstar))
                a['alpha_salpeter_needed'] = float(10 ** (a['dm_lens'] - SALP)) if conv == 'standard' else None
                a['alpha_posacki_2015'] = float(10 ** (0.38 * np.log10(a['sigma'] / 200) - 0.06))
                rows.append(a)
            res[conv][law] = rows
    (out / 'lenses.json').write_text(json.dumps(res, indent=2) + '\n')
    for conv in res:
        print(f'\n== {conv} convention ==')
        for law, rows in res[conv].items():
            print(f'  {law}:  needed stellar-mass change (dex): lensing | kinematics (isotropic) | gap   [alpha vs Salpeter needed | Posacki 2015]')
            for a in rows:
                extra = f"   [{a['alpha_salpeter_needed']:.2f} | {a['alpha_posacki_2015']:.2f}]" if a['alpha_salpeter_needed'] else ''
                print(f"    {a['name']}  logM*={a['logMstar_input']:.2f}  {a['dm_lens']:+.3f} | {a['dm_kin_isotropic']:+.3f} | {a['slip_gap_dex']:+.3f}{extra}")
            gaps = np.array([a['slip_gap_dex'] for a in rows])
            print(f"    mean lensing-minus-kinematics gap {gaps.mean():+.3f} dex (scatter {gaps.std(ddof=1):.3f})")


if __name__ == '__main__':
    main()
