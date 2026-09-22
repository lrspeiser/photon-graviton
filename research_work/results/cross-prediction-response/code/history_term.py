#!/usr/bin/env python3
"""JR-9 E: one universal history term, tested by cross-prediction not by joint fitting.

JR-9A found that the cross-predicted Einstein-angle error separates the six lenses
almost perfectly by foreground redshift, and JR-9C found that the required change
is far larger and steeper in redshift than the difference between the project's
static-Euclidean convention and a standard FLRW geometry. The obvious ordinary
explanation therefore does not cover it.

The hypothesis under test says the companion state has a *history*, not only a
spatial distribution. This script adds exactly one universal coefficient,

    A_chi -> A_chi * (1 + z_lens)^s

and nothing else. Every other R10 coefficient stays frozen. Two properties make
this a cheap test rather than another free parameter:

  * all 149 SPARC rotation galaxies sit at z ~ 0, so (1+z)^s = 1 and not one
    rotation curve moves, whatever s turns out to be;
  * s is fitted to the six cross-predicted angle residuals only, then the whole
    two-direction cross-prediction is re-run and scored with s held fixed.

This is development on historically examined objects. One universal parameter
fitted to six numbers is weak evidence, and it is reported as such.

    python history_term.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.linalg import cho_factor, cho_solve

import model as M

OBS = {x['Name']: x for x in json.loads(
    (M.ROOT / 'research_work/results/lensing-data-readiness/'
     'lens-observations-and-image-models.json').read_text())}
S_BOUNDS = (-6.0, 12.0)


def with_history(lens, dm, beta, s, want='both'):
    """R10 with the companion amplitude scaled by (1+z_lens)^s. Nothing else changes."""
    boost = (1 + OBS[lens.name]['zFG']) ** s
    p = M.params(dm, beta)
    local = lens.local_parameters(p, M.SPEC)
    A, rc, rt = M.J.source_parameters(lens.Mstar, 0., lens.Re, local, M.SPEC['saturation'], 1.)
    A = float(A) * boost
    rc, rt, q = float(rc), float(rt), float(local['q_sph'])
    scale = 10 ** local['logu']
    v = ang = None
    if want in ('both', 'kin'):
        g = scale * lens.gb + M.J.companion_force(lens.r, A, rc, rt, q)
        try:
            v = lens.moments(g, local, M.SPEC['radial_beta'])
        except ValueError:
            v = None
    if want in ('both', 'lens'):
        def f(b):
            r = b / np.cos(lens.lens_t)
            g = scale * lens.baryon_force(r) + M.J.companion_force(r, A, rc, rt, q)
            return lens.ratio * 4 / M.J.CLIGHT ** 2 * np.dot(lens.lens_w, g * r) - b / lens.Dl
        grid = lens.Re * np.geomspace(1e-6, 1e4, 200)
        vals = np.array([f(b) for b in grid])
        ii = np.flatnonzero((vals[:-1] > 0) & (vals[1:] <= 0))
        ang = float(brentq(f, grid[ii[-1]], grid[ii[-1] + 1], xtol=1e-10) / lens.Dl * M.J.ARCSEC) if len(ii) else None
    return v, ang


def kinematic_best(lens, chol, s):
    """Ordinary nuisance pair from the kinematics alone, at this history exponent."""
    from scipy.optimize import minimize
    dm_lim = 5.0 * lens.mass_log_error          # the same +/- 5 sigma window PROTOCOL.md declares
    def cost(x):
        dm = float(np.clip(x[0], -dm_lim, dm_lim))
        beta = float(np.clip(x[1], *M.BETA_BOUNDS))
        v, _ = with_history(lens, dm, beta, s, 'kin')
        if v is None:
            return 1e9
        d = v - lens.y
        return float(d @ cho_solve(chol, d)) + (dm / lens.mass_log_error) ** 2 \
            + (beta / M.BETA_PRIOR_SIGMA) ** 2
    best = min((minimize(cost, [d0, b0], method='Nelder-Mead',
                         options=dict(xatol=1e-6, fatol=1e-8, maxiter=800))
                for d0 in (-0.05, 0.0, 0.08) for b0 in (-0.3, 0.0, 0.3)), key=lambda r: r.fun)
    dm = float(np.clip(best.x[0], -dm_lim, dm_lim))
    beta = float(np.clip(best.x[1], *M.BETA_BOUNDS))
    v, _ = with_history(lens, dm, beta, s, 'kin')
    d = v - lens.y
    return dm, beta, float(d @ cho_solve(chol, d)), float(best.fun), v


def score(s, lenses, chols):
    """Cross-predicted angle residuals at history exponent s, kinematics-only nuisance."""
    rows = []
    for lens, chol in zip(lenses, chols):
        dm, beta, chi2, pen, v = kinematic_best(lens, chol, s)
        _, ang = with_history(lens, dm, beta, s, 'lens')
        rows.append(dict(name=lens.name, zFG=OBS[lens.name]['zFG'], dm_dex=dm, beta=beta,
                         dm_in_published_sigma=float(dm / lens.mass_log_error),
                         kinematic_chi2=chi2, kinematic_penalty=pen,
                         angle_observed_arcsec=float(lens.theta), angle_predicted_arcsec=ang,
                         angle_fractional_error=float(ang / lens.theta - 1),
                         vrms_fractional_RMS=float(np.sqrt(np.mean(((v - lens.y) / lens.y) ** 2)))))
    err = np.array([r['angle_fractional_error'] for r in rows])
    return rows, float(np.sqrt(np.mean(err ** 2)))


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    lenses = [M.make_lens(n) for n in M.LENSES]
    chols = [cho_factor(np.asarray(l.cov), lower=True) for l in lenses]

    base_rows, base_rms = score(0.0, lenses, chols)
    print(f's = 0 (frozen R10): cross-predicted angle fractional RMS {base_rms:.4%}', flush=True)

    fit = minimize_scalar(lambda s: score(s, lenses, chols)[1], bounds=S_BOUNDS, method='bounded',
                          options=dict(xatol=1e-4))
    s_hat = float(fit.x)
    best_rows, best_rms = score(s_hat, lenses, chols)
    print(f's = {s_hat:.4f} (fitted to six angles): angle fractional RMS {best_rms:.4%}', flush=True)

    # Profile the objective so the width of the preferred range is visible.
    profile = []
    for s in np.linspace(S_BOUNDS[0], S_BOUNDS[1], 37):
        _, rms = score(float(s), lenses, chols)
        profile.append(dict(s=float(s), angle_fractional_RMS=rms))

    payload = dict(
        experiment='JR-9E', protocol='../PROTOCOL.md',
        scope='One universal history coefficient on the companion amplitude, fitted to six '
              'cross-predicted Einstein-angle residuals. Development, not validation.',
        model='A_chi -> A_chi * (1 + z_lens)^s; every other R10 coefficient frozen',
        why_rotation_is_untouched='All 149 SPARC galaxies are at z ~ 0, so the factor is 1 '
                                  'and no rotation-curve prediction changes for any s.',
        s_fitted=s_hat, boost_at_z_0p25=float((1.25) ** s_hat),
        frozen_R10=dict(rows=base_rows, angle_fractional_RMS=base_rms),
        with_history=dict(rows=best_rows, angle_fractional_RMS=best_rms),
        profile=profile,
        degeneracy_warning='A companion amplitude growing with lens redshift and a redshift-'
                           'dependent error in the distance convention are not separable by '
                           'these six systems. JR-9C shows standard FLRW does not supply the '
                           'required size or slope, which weakens but does not remove that '
                           'degeneracy.',
        reserved_test='The KCWI release carries 14 SLACS lenses. Eight are in this project and '
                      'only six are usable. The remaining systems have never been scored here '
                      'and would need light components and population masses acquired before '
                      'they could test this trend. That is a genuinely reserved comparison and '
                      'must not be opened by fitting.',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'history-term.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')

    print('\nlens          zFG   err(s=0)%   err(s=fit)%   dm(s=0)  dm(fit)  fRMS_V(fit)')
    for a, b in zip(base_rows, best_rows):
        print(f"{a['name']:<12}{a['zFG']:>6.4f}{a['angle_fractional_error']*100:>11.2f}"
              f"{b['angle_fractional_error']*100:>14.2f}{a['dm_dex']:>10.4f}{b['dm_dex']:>9.4f}"
              f"{b['vrms_fractional_RMS']*100:>12.2f}%")
    print('WROTE', out / 'history-term.json')


if __name__ == '__main__':
    main()
