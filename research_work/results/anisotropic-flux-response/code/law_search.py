#!/usr/bin/env python3
"""JR-10 stage 2: find the law for the stress state, and test it by leave-one-out.

Stage 1 measured what Tolman factor tau = 1 + w_r + 2 w_t each lens demands. Six
numbers from six angles is not a test. This stage asks whether one universal law of
an independently measured source property reproduces them, and scores every
candidate by leave-one-out prediction of the Einstein angle -- fit the law on five
lenses, predict the sixth, never having used its angle.

PREDICTOR HYGIENE, DECLARED BEFORE SCORING
------------------------------------------
Any predictor built from the observed Einstein angle is contaminated, because tau
is solved from that same angle. The Einstein radius in kpc and theta_E/Re are
therefore excluded from the law search and reported separately as diagnostics only.
Redshift enters the lens equation through the distance convention, but JR-9C showed
that convention differs from standard FLRW by about 3 percent, nearly uniformly,
against a required tau swing of order 0.3. It is kept, with that caveat recorded.

The baseline to beat is tau = 1 everywhere: static dust, which is what R10 assumed.

    python law_search.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, hashlib, itertools, json
from pathlib import Path
import numpy as np

import model as M
from stress import angle_with_gamma, GAMMA_MIN_DEC

CONTAMINATED = ('einstein_radius_kpc', 'theta_E_over_Re')
EXCLUDE = CONTAMINATED + ('implied_tolman_factor', 'required_gamma_chi',
                          'implied_wr_if_transverse_zero', 'theta_observed_arcsec',
                          'angle_error_at_dust', 'theta_at_dust_gamma1_arcsec',
                          'kinematic_chi2', 'n_bins', 'vrms_fractional_RMS',
                          'dm_from_kinematics_dex', 'beta_from_kinematics')

# Declared candidate functional forms. x is the standardised predictor.
FORMS = {
    'constant':        (lambda x, p: np.full_like(np.asarray(x, float), p[0]), 1, [1.0]),
    'linear':          (lambda x, p: p[0] + p[1] * x, 2, [1.0, 0.0]),
    'power_of_1plusz': (lambda x, p: (1.0 + x) ** p[0], 1, [0.0]),
    'exponential':     (lambda x, p: p[0] * np.exp(p[1] * x), 2, [1.0, 0.0]),
}


def fit_form(form, xs, taus):
    """Least squares in tau. Closed form where possible, else a short search."""
    fn, npar, p0 = FORMS[form]
    xs, taus = np.asarray(xs, float), np.asarray(taus, float)
    if form == 'constant':
        return [float(taus.mean())]
    if form == 'linear':
        A = np.vstack([np.ones_like(xs), xs]).T
        return [float(c) for c in np.linalg.lstsq(A, taus, rcond=None)[0]]
    if form == 'exponential':
        if np.any(taus <= 0):
            return None
        A = np.vstack([np.ones_like(xs), xs]).T
        c = np.linalg.lstsq(A, np.log(taus), rcond=None)[0]
        return [float(np.exp(c[0])), float(c[1])]
    if form == 'power_of_1plusz':
        if np.any(taus <= 0) or np.any(1 + xs <= 0):
            return None
        num = float(np.sum(np.log(1 + xs) * np.log(taus)))
        den = float(np.sum(np.log(1 + xs) ** 2))
        return [num / den] if den > 0 else None
    return None


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--stage1', type=Path, default=Path('../run-v1/required-stress.json'))
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    stage1 = json.loads(args.stage1.read_text())
    rows = stage1['rows']
    names = [r['name'] for r in rows]
    tau = np.array([r['implied_tolman_factor'] for r in rows])
    lenses = {n: M.make_lens(n) for n in names}
    dm = {r['name']: r['dm_from_kinematics_dex'] for r in rows}
    theta_obs = np.array([r['theta_observed_arcsec'] for r in rows])

    predictors = [k for k in rows[0]
                  if isinstance(rows[0][k], (int, float)) and k not in EXCLUDE
                  and np.ptp([float(r[k]) for r in rows]) > 0]

    def angle_for(name, tau_value):
        if tau_value <= 0:
            return None
        return angle_with_gamma(lenses[name], dm[name], 1.0 / tau_value)

    # Baseline: static dust everywhere, zero fitted parameters.
    dust = np.array([angle_with_gamma(lenses[n], dm[n], 1.0) for n in names])
    baseline = dict(law='static dust, tau = 1', n_parameters=0,
                    angle_fractional_RMS=float(np.sqrt(np.mean((dust / theta_obs - 1) ** 2))),
                    per_lens=[dict(name=n, predicted=float(a), error=float(a / t - 1))
                              for n, a, t in zip(names, dust, theta_obs)])
    print(f"baseline  static dust tau=1                 LOO angle RMS = "
          f"{baseline['angle_fractional_RMS']:.4%}  (0 parameters)\n")

    results = []
    for key, form in itertools.product(predictors, FORMS):
        x = np.array([float(r[key]) for r in rows])
        full = fit_form(form, x, tau)
        if full is None:
            continue
        loo_err, loo_tau, band_ok = [], [], True
        for i in range(len(names)):
            keep = [j for j in range(len(names)) if j != i]
            p = fit_form(form, x[keep], tau[keep])
            if p is None:
                loo_err = None
                break
            t_pred = float(FORMS[form][0](np.array([x[i]]), p)[0])
            loo_tau.append(t_pred)
            if t_pred <= 0 or 1.0 / t_pred < GAMMA_MIN_DEC:
                band_ok = False
            a = angle_for(names[i], t_pred)
            if a is None:
                loo_err = None
                break
            loo_err.append(a / theta_obs[i] - 1)
        if not loo_err:
            continue
        loo_err = np.array(loo_err)
        fitted_tau = FORMS[form][0](x, full)
        in_sample = np.array([angle_for(n, float(t)) or np.nan
                              for n, t in zip(names, fitted_tau)])
        results.append(dict(
            predictor=key, form=form, n_parameters=FORMS[form][1],
            parameters=[float(v) for v in full],
            tau_fitted=[float(v) for v in fitted_tau],
            in_sample_angle_fractional_RMS=float(np.sqrt(np.nanmean((in_sample / theta_obs - 1) ** 2))),
            loo_angle_fractional_RMS=float(np.sqrt(np.mean(loo_err ** 2))),
            loo_max_abs_angle_error=float(np.max(np.abs(loo_err))),
            loo_tau=[float(v) for v in loo_tau],
            all_loo_inside_energy_band=bool(band_ok),
            loo_per_lens=[dict(name=n, tau_predicted=float(t), angle_error=float(e))
                          for n, t, e in zip(names, loo_tau, loo_err)]))

    results.sort(key=lambda r: r['loo_angle_fractional_RMS'])
    print(f"{'predictor':<28}{'form':<18}{'par':>4}{'in-sample':>11}{'LOO RMS':>10}"
          f"{'LOO max':>10}  band")
    for r in results[:14]:
        print(f"{r['predictor']:<28}{r['form']:<18}{r['n_parameters']:>4}"
              f"{r['in_sample_angle_fractional_RMS']:>10.2%}{r['loo_angle_fractional_RMS']:>10.2%}"
              f"{r['loo_max_abs_angle_error']:>10.2%}  "
              f"{'yes' if r['all_loo_inside_energy_band'] else 'NO'}")

    diagnostics = []
    for key in CONTAMINATED:
        x = np.array([float(r[key]) for r in rows])
        diagnostics.append(dict(predictor=key,
                                rank_correlation_with_tau=float(np.corrcoef(
                                    np.argsort(np.argsort(tau)), np.argsort(np.argsort(x)))[0, 1]),
                                note='Built from the observed Einstein angle; excluded from the '
                                     'law search as circular.'))

    payload = dict(
        experiment='JR-10 stage 2', stage1=str(args.stage1),
        scope='Universal-law search over independently measured source properties, scored by '
              'leave-one-out prediction of the Einstein angle. Exploratory search over a '
              'declared grid; look-elsewhere applies and is not corrected for.',
        n_predictors=len(predictors), n_forms=len(FORMS),
        n_candidates_scored=len(results),
        predictors_used=predictors,
        excluded_as_circular=list(CONTAMINATED),
        baseline_static_dust=baseline,
        ranked=results, circular_diagnostics=diagnostics,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'law-search.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    print('\nWROTE', out / 'law-search.json')


if __name__ == '__main__':
    main()
