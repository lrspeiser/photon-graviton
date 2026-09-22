#!/usr/bin/env python3
"""JR-9 A: predict each SLACS lens observable from the other one, not from both.

Universal R10 coefficients are frozen. Only the declared ordinary nuisance pair
(stellar log-mass offset, orbital anisotropy) is inferred, from one observable at
a time. See ../PROTOCOL.md, declared before any score below was computed.

    python cross_predict.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import cho_factor, cho_solve

import model as M

DM_HALF_WIDTH_SIGMA = 5.0    # dm grid spans +/- 5 published stellar-mass sigma
DM_GRID = 241
BETA_GRID = 136              # 0.01 steps across the JR-1 anisotropy bounds
ANGLE_SIGMA_SCAN = (0.01, 0.02, 0.03, 0.05)


def kin_chi2(lens, chol, dm, beta, companion=True):
    v = M.vrms(lens, dm, beta, companion)
    if v is None:
        return np.inf, None
    d = v - lens.y
    return float(d @ cho_solve(chol, d)), v


def kinematic_posterior(lens, companion=True):
    """Grid posterior over (dm, beta) from V_rms and its covariance, plus priors."""
    chol = cho_factor(np.asarray(lens.cov), lower=True)
    sigma_m = float(lens.mass_log_error)
    dm_axis = np.linspace(-DM_HALF_WIDTH_SIGMA * sigma_m, DM_HALF_WIDTH_SIGMA * sigma_m, DM_GRID)
    beta_axis = np.linspace(M.BETA_BOUNDS[0], M.BETA_BOUNDS[1], BETA_GRID)
    logp = np.full((DM_GRID, BETA_GRID), -np.inf)
    for i, dm in enumerate(dm_axis):
        for j, beta in enumerate(beta_axis):
            c2, _ = kin_chi2(lens, chol, dm, beta, companion)
            if np.isfinite(c2):
                logp[i, j] = -0.5 * (c2 + (dm / sigma_m) ** 2 + (beta / M.BETA_PRIOR_SIGMA) ** 2)
    if not np.any(np.isfinite(logp)):
        raise RuntimeError('No valid kinematic model anywhere on the declared grid: ' + lens.name)
    w = np.exp(logp - np.nanmax(logp))
    w /= w.sum()
    return dict(dm_axis=dm_axis, beta_axis=beta_axis, weight=w, chol=chol, sigma_m=sigma_m,
                dm_marginal=w.sum(axis=1), beta_marginal=w.sum(axis=0), logp=logp)


def quantiles(axis, pmf, qs=(0.025, 0.16, 0.5, 0.84, 0.975)):
    cdf = np.cumsum(pmf) / pmf.sum()
    return [float(np.interp(q, cdf, axis)) for q in qs]


def dm_from_angle(lens, companion=True, span=1.2):
    """The single stellar offset that reproduces the catalog Einstein angle."""
    f = lambda dm: M.einstein_angle(lens, dm, companion) - lens.theta
    lo, hi = -span, span
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None, dict(bracket=[lo, hi], angle_at_bracket=[flo + lens.theta, fhi + lens.theta],
                          note='catalog angle unreachable inside +/- %.1f dex of stellar mass' % span)
    root = float(brentq(f, lo, hi, xtol=1e-10, rtol=1e-12))
    return root, dict(residual_arcsec=float(f(root)))


def best_beta_given_dm(lens, chol, dm, companion=True):
    """Anisotropy is invisible to lensing, so scan it with its prior when predicting."""
    betas = np.linspace(M.BETA_BOUNDS[0], M.BETA_BOUNDS[1], BETA_GRID)
    best, out = np.inf, None
    logw = []
    for b in betas:
        c2, v = kin_chi2(lens, chol, dm, b, companion)
        pen = c2 + (b / M.BETA_PRIOR_SIGMA) ** 2
        logw.append(-0.5 * pen if np.isfinite(pen) else -np.inf)
        if pen < best:
            best, out = pen, (float(b), float(c2), v)
    logw = np.asarray(logw)
    w = np.exp(logw - np.nanmax(logw)) if np.any(np.isfinite(logw)) else np.zeros_like(logw)
    return out, betas, (w / w.sum() if w.sum() > 0 else w)


def run_lens(name, companion=True):
    lens = M.make_lens(name)
    rec = dict(name=name, companion=companion, theta_observed_arcsec=float(lens.theta),
               stellar_log_mass_sigma=float(lens.mass_log_error),
               Mstar_catalog_Msun=float(lens.Mstar), Re_kpc=float(lens.Re),
               n_bins=int(len(lens.y)), observed_vrms_kms=lens.y.tolist(),
               inner_arcsec=lens.inner_arcsec.tolist(), outer_arcsec=lens.outer_arcsec.tolist(),
               theta_E_over_Re=float(lens.theta / (lens.Re / lens.Dl * M.J.ARCSEC)))

    # ---- Direction K -> L: kinematics only, then predict the angle -------------
    post = kinematic_posterior(lens, companion)
    dm_q = quantiles(post['dm_axis'], post['dm_marginal'])
    beta_q = quantiles(post['beta_axis'], post['beta_marginal'])
    flat = np.unravel_index(np.argmax(post['logp']), post['logp'].shape)
    dm_map, beta_map = float(post['dm_axis'][flat[0]]), float(post['beta_axis'][flat[1]])
    c2_map, v_map = kin_chi2(lens, post['chol'], dm_map, beta_map, companion)
    # theta_E depends on dm alone, so push the dm marginal through the deflection map.
    angle_of_dm = np.array([M.einstein_angle(lens, d, companion) for d in post['dm_axis']])
    ang_q = quantiles(angle_of_dm, post['dm_marginal'])
    ang_mean = float(np.sum(angle_of_dm * post['dm_marginal']))
    ang_sd = float(np.sqrt(max(np.sum(angle_of_dm ** 2 * post['dm_marginal']) - ang_mean ** 2, 0.)))
    rec['K_to_L'] = dict(
        dm_quantiles_dex=dm_q, dm_MAP_dex=dm_map, beta_quantiles=beta_q, beta_MAP=beta_map,
        kinematic_chi2_at_MAP=float(c2_map), predicted_vrms_at_MAP_kms=v_map.tolist(),
        vrms_fractional_RMS=float(np.sqrt(np.mean(((v_map - lens.y) / lens.y) ** 2))),
        predicted_angle_quantiles_arcsec=ang_q, predicted_angle_mean_arcsec=ang_mean,
        predicted_angle_sd_arcsec=ang_sd,
        angle_fractional_error_at_median=float(ang_q[2] / lens.theta - 1),
        angle_gap_arcsec=float(lens.theta - ang_q[2]),
        angle_gap_in_kinematic_sigma=float((lens.theta - ang_mean) / ang_sd) if ang_sd > 0 else None,
        tension_vs_assumed_angle_sigma={
            f'{s:.0%}': float((lens.theta - ang_mean) / np.hypot(ang_sd, s * lens.theta))
            for s in ANGLE_SIGMA_SCAN},
        required_fractional_angle_sigma_for_1sigma=float(
            np.sqrt(max((lens.theta - ang_mean) ** 2 - ang_sd ** 2, 0.)) / lens.theta))

    # ---- Direction L -> K: the angle only, then predict the motions ------------
    dm_lens, solve_info = dm_from_angle(lens, companion)
    entry = dict(solve=solve_info, dm_dex=dm_lens)
    if dm_lens is not None:
        (b_best, c2_best, v_best), betas, bw = best_beta_given_dm(lens, post['chol'], dm_lens, companion)
        vs = np.array([M.vrms(lens, dm_lens, b, companion) for b in betas], dtype=object)
        ok = np.array([v is not None for v in vs])
        stack = np.stack([np.asarray(v, float) for v in vs[ok]])
        wm = bw[ok] / bw[ok].sum() if bw[ok].sum() > 0 else np.full(ok.sum(), 1 / ok.sum())
        v_marg = stack.T @ wm
        d = v_marg - lens.y
        entry.update(
            dm_in_published_sigma=float(dm_lens / lens.mass_log_error),
            stellar_mass_factor=float(10 ** dm_lens),
            beta_best=b_best, kinematic_chi2_best_beta=float(c2_best),
            predicted_vrms_best_beta_kms=v_best.tolist(),
            vrms_fractional_RMS_best_beta=float(np.sqrt(np.mean(((v_best - lens.y) / lens.y) ** 2))),
            predicted_vrms_marginalized_kms=v_marg.tolist(),
            kinematic_chi2_marginalized=float(d @ cho_solve(post['chol'], d)),
            vrms_fractional_RMS_marginalized=float(np.sqrt(np.mean((d / lens.y) ** 2))))
    rec['L_to_K'] = entry

    # ---- Consistency ---------------------------------------------------------
    if dm_lens is not None:
        lo, hi = dm_q[1], dm_q[3]
        sd_dm = 0.5 * (hi - lo)
        rec['consistency'] = dict(
            dm_kinematics_median_dex=dm_q[2], dm_lensing_dex=dm_lens,
            delta_dm_dex=float(dm_lens - dm_q[2]),
            delta_dm_in_published_stellar_sigma=float((dm_lens - dm_q[2]) / lens.mass_log_error),
            delta_dm_in_kinematic_posterior_sigma=float((dm_lens - dm_q[2]) / sd_dm) if sd_dm > 0 else None,
            extra_stellar_mass_factor_demanded_by_lens=float(10 ** (dm_lens - dm_q[2])),
            lensing_dm_inside_kinematic_68=bool(lo <= dm_lens <= hi),
            chi2_penalty_of_lensing_solution=float(
                rec['L_to_K']['kinematic_chi2_marginalized'] - float(c2_map)))
    return rec


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory; old results are never overwritten')
    out.mkdir(parents=True)
    t0 = time.monotonic()

    # Bit-for-bit reproduction of the published R10 lens predictions, first.
    pub = {r['name']: r for r in M.RESULTS['R10_lens_predictions']}
    repro = []
    for name in M.LENSES:
        lens = M.make_lens(name)
        nu = M.JR1_NUISANCE[name]
        a = M.einstein_angle(lens, nu['dm'])
        v = M.vrms(lens, nu['dm'], nu['beta'])
        repro.append(dict(name=name, published_angle=pub[name]['predicted_angle_arcsec'],
                          reproduced_angle=a,
                          angle_relative_difference=a / pub[name]['predicted_angle_arcsec'] - 1,
                          published_stellar_fractional_RMS=pub[name]['stellar_fractional_RMS'],
                          reproduced_stellar_fractional_RMS=float(
                              np.sqrt(np.mean(((v - lens.y) / lens.y) ** 2)))))
    worst = max(abs(r['angle_relative_difference']) for r in repro)
    if worst > 1e-6:
        raise RuntimeError('R10 reproduction drifted by %.3g; refusing to score' % worst)

    results = dict(
        experiment='JR-9A',
        scope='Same-object cross-prediction with frozen universal R10 coefficients. '
              'Historically examined development objects; not a blind test.',
        protocol='../PROTOCOL.md',
        frozen_universal_parameters=M.R10_UNIVERSAL,
        priors=dict(stellar_log_mass='published per-lens log10 mass error',
                    beta_sigma=M.BETA_PRIOR_SIGMA, beta_bounds=list(M.BETA_BOUNDS)),
        grid=dict(dm_points=DM_GRID, dm_half_width_in_sigma=DM_HALF_WIDTH_SIGMA,
                  beta_points=BETA_GRID),
        r10_reproduction=repro,
        r10_reproduction_worst_relative_angle_difference=worst,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

    for label, companion in (('companion_R10', True), ('baryons_only_control', False)):
        rows = []
        for name in M.LENSES:
            rec = run_lens(name, companion)
            rows.append(rec)
            k, l = rec['K_to_L'], rec['L_to_K']
            print(f"[{label}] {name}: K->L theta={k['predicted_angle_quantiles_arcsec'][2]:.4f}"
                  f" ({k['angle_fractional_error_at_median']:+.2%})"
                  f"  L->K dm={l['dm_dex'] if l['dm_dex'] is None else round(l['dm_dex'], 4)}"
                  f" chi2={l.get('kinematic_chi2_marginalized')}", flush=True)
        results[label] = rows
        good = [r for r in rows if r.get('consistency')]
        if good:
            results[label + '_summary'] = dict(
                n_systems=len(rows), n_with_lensing_solution=len(good),
                median_delta_dm_dex=float(np.median([r['consistency']['delta_dm_dex'] for r in good])),
                mean_delta_dm_in_published_sigma=float(np.mean(
                    [r['consistency']['delta_dm_in_published_stellar_sigma'] for r in good])),
                n_lensing_dm_inside_kinematic_68=sum(
                    r['consistency']['lensing_dm_inside_kinematic_68'] for r in good),
                K_to_L_angle_fractional_RMS=float(np.sqrt(np.mean(
                    [r['K_to_L']['angle_fractional_error_at_median'] ** 2 for r in rows]))),
                L_to_K_total_chi2=float(sum(r['L_to_K']['kinematic_chi2_marginalized'] for r in good)),
                L_to_K_total_bins=int(sum(r['n_bins'] for r in good)),
                K_to_L_total_chi2=float(sum(r['K_to_L']['kinematic_chi2_at_MAP'] for r in rows)))

    results['seconds'] = time.monotonic() - t0
    (out / 'cross-prediction-results.json').write_text(
        json.dumps(results, indent=2, allow_nan=False) + '\n')
    print('WROTE', out / 'cross-prediction-results.json', flush=True)


if __name__ == '__main__':
    main()
