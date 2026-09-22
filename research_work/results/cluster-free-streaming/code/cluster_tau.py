#!/usr/bin/env python3
"""JR-11: clusters are short by exactly the free-streaming flux factor.

The galaxy-calibrated companion law transferred to the eleven X-COP clusters
under-predicts the extra mass their X-ray hydrostatic profiles require, by a
remarkably tight factor near two. That factor is not a free number.

The Tolman active-mass factor tau = 1 + w_r + 2 w_t says how much gravitational
pull a given energy content exerts on slow matter. Two states are exact:

    static dust                  w_r = 0, w_t = 0   ->  tau = 1
    free-streaming radial flux   w_r = 1, w_t = 0   ->  tau = 2

A galaxy-calibrated amplitude carries tau_galaxy inside it, because rotation curves
are a dynamical measurement. Transferring that amplitude to clusters therefore
predicts their dynamical mass short by tau_cluster / tau_galaxy. If galaxies sit in
the dust state and clusters in the free-streaming state, that ratio is exactly 2,
with nothing fitted.

This script tests that against the measured shortfall and against the alternative
that the factor is simply a fitted constant.

    python cluster_tau.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[4]
CMP = ROOT / 'companion_wave_test/cluster_comparison.csv'
PRED = ROOT / 'companion_deposition_fit/cluster_predictions.csv'

TAU_DUST = 1.0
TAU_FREE_STREAMING = 2.0


def crossover(R, R0, m):
    """tau(R) = 1 + 1/(1 + (R0/R)^m): dust at small R, free-streaming at large R."""
    return 1.0 + 1.0 / (1.0 + (R0 / np.asarray(R, float)) ** m)


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    cc = {r['cluster']: r for r in csv.DictReader(CMP.open())}
    cp = list(csv.DictReader(PRED.open()))
    names = [p['cluster'] for p in cp]
    ratio = np.array([float(p['required_over_predicted']) for p in cp])
    M500 = np.array([float(cc[n]['M500_1e14Msun']) for n in names])
    eM500 = np.array([float(cc[n]['eM500_1e14Msun']) for n in names])
    R500 = np.array([float(cc[n]['R500_Mpc']) * 1000.0 for n in names])
    fgas = np.array([float(cc[n]['fgas500']) for n in names])

    def scatter(tau):
        return float(np.sqrt(np.mean((ratio / tau - 1) ** 2)))

    fitted_constant = float(ratio.mean())
    sem = float(ratio.std(ddof=1) / np.sqrt(len(ratio)))
    tests = [
        dict(model='static dust, tau = 1 (what the transfer assumed)', tau=TAU_DUST,
             n_parameters=0, fractional_scatter=scatter(TAU_DUST)),
        dict(model='free-streaming radial flux, tau = 2 (nothing fitted)',
             tau=TAU_FREE_STREAMING, n_parameters=0, fractional_scatter=scatter(TAU_FREE_STREAMING)),
        dict(model='best fitted constant', tau=fitted_constant, n_parameters=1,
             fractional_scatter=scatter(fitted_constant))]

    p, cov = curve_fit(crossover, R500, ratio, p0=[500.0, 2.0], maxfev=40000)
    tau_R = crossover(R500, *p)
    tests.append(dict(model='crossover tau(R) = 1 + 1/(1+(R0/R)^m)', n_parameters=2,
                      R0_kpc=float(p[0]), m=float(p[1]),
                      fractional_scatter=float(np.sqrt(np.mean((ratio / tau_R - 1) ** 2)))))

    rk = lambda a: np.argsort(np.argsort(a))
    payload = dict(
        experiment='JR-11',
        scope='Transfer test of the Tolman stress factor from galaxies to the eleven X-COP '
              'clusters. Uses the archived galaxy-calibrated acceleration transfer and the '
              'archived X-ray hydrostatic table; no new fit to cluster data beyond the two '
              'comparison models below.',
        why_the_ratio_is_tau='Rotation curves are a dynamical measurement, so a galaxy-'
                             'calibrated amplitude already contains tau_galaxy. Transferring '
                             'it predicts cluster dynamical mass short by tau_cluster/tau_galaxy.',
        exact_states={'static dust': TAU_DUST, 'free-streaming radial flux': TAU_FREE_STREAMING},
        observed=dict(n_clusters=len(ratio), median=float(np.median(ratio)),
                      mean=fitted_constant, standard_deviation=float(ratio.std(ddof=1)),
                      standard_error_of_mean=sem, minimum=float(ratio.min()),
                      maximum=float(ratio.max()),
                      distance_of_mean_from_two_in_standard_errors=float(abs(fitted_constant - 2) / sem),
                      quoted_M500_fractional_errors=[float(x) for x in (eM500 / M500)]),
        tests=tests,
        per_cluster=[dict(name=n, required_over_predicted=float(x), M500_1e14=float(a),
                          R500_kpc=float(b), fgas500=float(f),
                          residual_vs_tau2=float(x / 2 - 1),
                          tau_from_crossover=float(t))
                     for n, x, a, b, f, t in zip(names, ratio, M500, R500, fgas, tau_R)],
        residual_structure={k: dict(rank=float(np.corrcoef(rk(ratio), rk(v))[0, 1]),
                                    pearson=float(np.corrcoef(ratio, v)[0, 1]))
                            for k, v in (('M500', M500), ('R500_kpc', R500), ('fgas500', fgas))},
        crossover_extrapolation=[dict(R_kpc=float(R), tau=float(crossover(R, *p)))
                                 for R in (10, 30, 100, 300, 600, 900, 1200, 2000, 5000)],
        limitations=[
            'Eleven clusters spanning only R500 = 1054 to 1430 kpc cannot constrain the '
            'sharpness m of the crossover; only that galaxies sit at tau = 1 and these '
            'clusters near tau = 2.',
            'X-ray hydrostatic masses carry their own known bias of order ten to twenty '
            'percent, which is the same size as the residual scatter here.',
            'tau = 2 is a claim that the companion is a free-streaming radial null flux in '
            'clusters. That is a physical statement needing its own independent check; this '
            'test only shows the number it implies is the number the data wants.',
            'No cluster weak-lensing measurement is used here. The shortfall tested is in '
            'hydrostatic dynamical mass.'],
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        inputs={str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest()
                for f in (CMP, PRED)})
    (out / 'cluster-tau.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')

    print(f"observed shortfall over {len(ratio)} X-COP clusters:")
    print(f"  median {np.median(ratio):.3f}   mean {fitted_constant:.3f} +- {sem:.3f} (s.e.m.)"
          f"   range {ratio.min():.3f} to {ratio.max():.3f}\n")
    for t in tests:
        extra = f"  R0={t['R0_kpc']:.0f} kpc, m={t['m']:.1f}" if 'R0_kpc' in t else ''
        print(f"  {t['model']:<52} scatter {t['fractional_scatter']:>6.1%}"
              f"  ({t['n_parameters']} par){extra}")
    print(f"\n  the mean sits {payload['observed']['distance_of_mean_from_two_in_standard_errors']:.2f}"
          f" standard errors from 2")
    print('\ncrossover extrapolation:')
    for e in payload['crossover_extrapolation']:
        print(f"  R = {e['R_kpc']:>5} kpc   tau = {e['tau']:.3f}")
    print('\nWROTE', out / 'cluster-tau.json')


if __name__ == '__main__':
    main()
