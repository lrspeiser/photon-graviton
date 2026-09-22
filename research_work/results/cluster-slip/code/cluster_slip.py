#!/usr/bin/env python3
"""T1.2: the published hydrostatic mass bias, read as a slip measurement.

Hydrostatic X-ray masses measure the dynamical potential Phi. Weak lensing measures
(Phi + Psi)/2. With the corrected field equations, Psi = Phi_N + eta (Phi - Phi_N), so
in any system

    M_WL / M_HSE  =  1 + f (eta - 1) / 2,        f = 1 - M_baryon / M_HSE

where f is the fraction of the dynamical mass that is not Newtonian. Every X-ray/lensing
comparison in the literature therefore already contains a slip measurement, once f is
known - and for clusters f is fixed by the measured gas fraction.

For X-COP the published comparison gives M_HSE 10 to 15 percent below lensing
(Eckert et al. 2022), i.e. M_WL / M_HSE = 1/(1-b).

This does not claim the bias IS slip. The field attributes it to non-thermal pressure,
with good reason. It asks what slip the bias would imply if it were, and whether that
value agrees with the slip measured independently on galaxy-scale lenses. And it names
the observation that separates the two readings.

    python cluster_slip.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
BIAS = (0.07, 0.10, 0.125, 0.15)
SLACS_ETA_UNDERBENT = (1.36, 1.54)
SLACS_ETA_OVERBENT = (0.82, 0.88)


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    cc = {x['cluster']: x for x in csv.DictReader((ROOT / 'companion_wave_test/cluster_comparison.csv').open())}
    rows = []
    for p in csv.DictReader((ROOT / 'companion_deposition_fit/cluster_predictions.csv').open()):
        x = cc[p['cluster']]
        M500 = float(x['M500_1e14Msun'])
        fb = float(x['fgas500']) + 0.02
        f = 1.0 - fb                      # non-Newtonian fraction of the hydrostatic mass
        etas = {f'{b:.3f}': 1 + 2 * (1 / (1 - b) - 1) / f for b in BIAS}
        rows.append(dict(name=p['cluster'], M500_1e14=M500, baryon_fraction=fb,
                         non_newtonian_fraction=f, implied_eta=etas))

    table = {f'{b:.3f}': [r['implied_eta'][f'{b:.3f}'] for r in rows] for b in BIAS}
    summary = {b: dict(mean=float(np.mean(v)), min=float(np.min(v)), max=float(np.max(v)))
               for b, v in table.items()}
    payload = dict(
        task='T1.2',
        relation='M_WL / M_HSE = 1 + f (eta - 1)/2, f = 1 - M_baryon / M_HSE',
        bias_source='Eckert et al. 2022: X-COP hydrostatic masses 10-15 percent below lensing; '
                    'Eckert et al. 2019: 7 percent at R500',
        rows=rows, implied_eta_by_bias=summary,
        slacs_eta=dict(underbent=list(SLACS_ETA_UNDERBENT), overbent=list(SLACS_ETA_OVERBENT)),
        not_claimed='That the hydrostatic bias is slip. It is conventionally attributed to '
                    'non-thermal pressure support, and simulations support that reading.',
        discriminating_test='Non-thermal pressure predicts the bias grows with dynamical '
                            'disturbance: merging, unrelaxed clusters read lower. Slip predicts '
                            'it tracks the non-Newtonian fraction f and does not care about '
                            'dynamical state. Split a lensing-and-X-ray sample into relaxed and '
                            'disturbed at fixed f: if the bias is the same in both, slip; if the '
                            'disturbed half carries it, pressure.',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'cluster-slip.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')

    print('T1.2 - what slip the published hydrostatic bias implies\n')
    print(f"  non-Newtonian fraction f across X-COP: "
          f"{min(r['non_newtonian_fraction'] for r in rows):.3f} - {max(r['non_newtonian_fraction'] for r in rows):.3f}\n")
    print(f"  {'bias b':<10}{'M_WL/M_HSE':>12}{'implied eta (mean)':>21}{'range':>16}")
    for b in BIAS:
        s = summary[f'{b:.3f}']
        print(f"  {b:<10.3f}{1/(1-b):>12.3f}{s['mean']:>21.3f}{f'{s[chr(109)+chr(105)+chr(110)]:.2f}-{s[chr(109)+chr(97)+chr(120)]:.2f}':>16}")
    print(f"\n  SLACS lenses, measured independently:")
    print(f"     under-bent  eta = {SLACS_ETA_UNDERBENT[0]} - {SLACS_ETA_UNDERBENT[1]}")
    print(f"     over-bent   eta = {SLACS_ETA_OVERBENT[0]} - {SLACS_ETA_OVERBENT[1]}")
    print('\n  DISCRIMINATING TEST:', payload['discriminating_test'])
    print('\nWROTE', out / 'cluster-slip.json')


if __name__ == '__main__':
    main()
