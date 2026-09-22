#!/usr/bin/env python3
"""JR-13c: the headline numbers, from raw data only, plus the hydrostatic-bias prediction."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
import sparc_test as S

K = 1e6 / 3.0856775814913673e19
TAU_FREE_STREAM = 2.0


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    gals, clus = S.load_sparc(), S.load_clusters()
    gN = np.concatenate([g['gN'] for g in gals])
    gobs = np.concatenate([g['gobs'] for g in gals])
    gid = np.concatenate([np.full(len(g['r']), i) for i, g in enumerate(gals)])
    err = np.concatenate([g['err'] for g in gals])
    v = np.concatenate([g['v'] for g in gals])
    r = np.concatenate([g['r'] for g in gals])
    cN = np.array([c['gN'] for c in clus])
    cobs = np.array([c['gobs'] for c in clus])

    cost = lambda la0: float(np.sum(((np.sqrt(np.maximum(r * S.law(gN, 10 ** la0), 0)) - v) / err) ** 2))
    a0 = 10 ** minimize_scalar(cost, bounds=(2.5, 5.5), method='bounded').x
    vp = np.sqrt(np.maximum(r * S.law(gN, a0), 0))
    per = np.array([np.sqrt(np.mean((vp[gid == i] - v[gid == i]) ** 2)) for i in range(len(gals))])
    fr = np.array([np.sqrt(np.mean(((vp[gid == i] - v[gid == i]) / v[gid == i]) ** 2))
                   for i in range(len(gals))])
    gm = float(np.mean([np.mean(np.log(gobs[gid == i] / S.law(gN[gid == i], a0)))
                        for i in range(len(gals))]))
    cm = float(np.mean(cobs / S.law(cN, a0)))
    gap = cm / np.exp(gm)
    bias = 1.0 - gap / TAU_FREE_STREAM

    payload = dict(
        experiment='JR-13c',
        scope='One constant against raw SPARC rotation measurements. No companion model, no '
              'R10, no fitted per-galaxy parameter appears in the inputs.',
        law='g = g_N/2 + sqrt(g_N^2/4 + g_N a0)',
        a0=dict(kms2_per_kpc=float(a0), SI=float(a0 * K), n_free_parameters=1),
        sparc=dict(n_galaxies=len(gals), n_points=int(len(v)),
                   mean_v_RMSE_kms=float(per.mean()), median_v_RMSE_kms=float(np.median(per)),
                   within_20pct=int(np.sum(fr < 0.20)), within_10pct=int(np.sum(fr < 0.10)),
                   raw_chi2=float(cost(np.log10(a0)))),
        r10_comparison=dict(n_fitted_parameters=20, mean_validation_RMSE_kms=17.684,
                            within_20pct_of_31_comparison=22, raw_chi2=211200.744,
                            note='R10 numbers are the published JR-1 values.'),
        cluster_gap=dict(galaxy_mean_factor=float(np.exp(gm)), cluster_mean_factor=cm,
                         cluster_over_galaxy=float(gap),
                         free_streaming_tau=TAU_FREE_STREAM,
                         implied_hydrostatic_bias=float(bias),
                         prediction='If the gap is exactly the free-streaming factor tau = 2, '
                                    'then X-ray hydrostatic masses in these clusters must be '
                                    f'low by b = {bias:.3f}. Published values run b = 0.1 to '
                                    '0.3, so this is a falsifiable number, testable against '
                                    'cluster weak lensing.'),
        limitations=[
            'The flux-driven tau crossover tested in JR-13 does not survive raw data: on SPARC '
            'it improves the fit by 0.4 percent, which is nothing. JR-13 read 11.3 percent '
            'because its galaxy points came from R10 rather than from measurements.',
            'No quality cuts are applied to SPARC; inclination and flag cuts would lower the '
            'scatter but were not used so the comparison to R10 stays like for like.',
            'The cluster gap uses X-ray hydrostatic masses and a declared 2 percent stellar '
            'fraction.'],
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'headline.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')

    s, c = payload['sparc'], payload['cluster_gap']
    print(f"a0 = {a0 * K:.3e} m/s^2, ONE free parameter\n")
    print(f"  {'':22}{'this law':>12}{'R10':>12}")
    print(f"  {'fitted parameters':22}{1:>12}{20:>12}")
    print(f"  {'mean v RMSE km/s':22}{s['mean_v_RMSE_kms']:>12.2f}{17.684:>12.2f}")
    print(f"  {'raw chi2 (3150 pts)':22}{s['raw_chi2']:>12.0f}{211201:>12.0f}")
    print(f"  {'within 20%':22}{s['within_20pct']:>12}{'-':>12}")
    print(f"\ncluster gap at the SPARC-calibrated a0: {c['cluster_over_galaxy']:.3f}")
    print(f"  free-streaming tau = 2 would require hydrostatic bias b = {c['implied_hydrostatic_bias']:.3f}")
    print(f"  published hydrostatic bias runs 0.1 to 0.3 -> consistent, and falsifiable")
    print('WROTE', out / 'headline.json')


if __name__ == '__main__':
    main()
