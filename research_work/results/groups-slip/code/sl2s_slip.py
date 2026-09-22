#!/usr/bin/env python3
"""T1.3: can SL2S galaxy groups measure the slip eta? Using the published table itself.

Source: Munoz et al. 2013, A&A, "Dynamical analysis of strong-lensing galaxy groups at
intermediate redshift", arXiv:1212.2624, Tables 2 and 3, transcribed below by hand from
the extracted text and checked row by row.

CORRECTION RECORDED HERE: an earlier revision of THEORY.md and the notebook said this
paper finds weak-lensing masses about 50 percent above dynamical ones. It does not.
The only '50%' in the paper refers to the fraction of galaxies that live in groups.
That figure came from a search-engine summary and was never in the source. The paper's
own finding runs the other way: from simulations it concludes the velocity dispersion
is always underestimated, so dynamical masses read low.

Method. The weak-lensing masses are projected within 2 Mpc; the virial masses are
within the projected virial radius, 0.2 to 1.1 Mpc. They cannot be ratioed directly.
Where the extra gravity dominates, the law gives g proportional to 1/r, an isothermal
sphere, whose projected mass within R is pi sigma^2 R / G. So the dynamical prediction
for the 2 Mpc weak-lensing aperture is pi sigma^2 (2 Mpc) / G, and lensing sees
(1 + eta)/2 of it. Solving for eta then carries the paper's own warning: an
underestimated sigma inflates every eta.

    python sl2s_slip.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np

G = 4.30091727e-6                     # kpc (km/s)^2 / Msun
APERTURE_KPC = 2000.0

# Table 2 sigma (km/s, +err, -err) and Table 3 M_V (1e14 Msun) and M_WL(2 Mpc).
# Two fields hold two structures along the line of sight; weak lensing sees both.
GROUPS = [
    dict(name='SL2SJ02140-0535', z=[0.445], sigma=[(364, 60, 137)], MWL=(5.5, 3.7)),
    dict(name='SL2SJ08544-0121', z=[0.351, 0.356], sigma=[(185, 30, 62), (341, 43, 109)], MWL=(6.3, 2.5)),
    dict(name='SL2SJ09413-1100', z=[0.384], sigma=[(350, 57, 210)], MWL=(3.7, 3.4)),
]
EXCLUDED = {
    'SL2SJ02141-0405': 'no weak-lensing mass: Einstein radius below 3 arcsec, galaxy regime',
    'SL2SJ02180-0515': 'no weak-lensing mass: Einstein radius below 3 arcsec, galaxy regime',
    'SL2SJ02215-0647': 'weak-lensing upper limit only (< 3.1e14)',
    'SL2SJ08591-0345': 'no weak-lensing mass: at the edge of the field of view',
}
SIGMA_BIAS = (1.0, 0.9, 0.8, 0.7)     # true sigma = measured / factor


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    rows = []
    rng = np.random.default_rng(20260922)
    for g in GROUPS:
        per_bias = {}
        for f in SIGMA_BIAS:
            # isothermal projected mass in the 2 Mpc aperture, summed over structures on the sightline
            Mdyn = sum(np.pi * (s / f) ** 2 * APERTURE_KPC / G for s, _, _ in g['sigma']) / 1e14
            ratio = g['MWL'][0] / Mdyn
            eta = 2 * ratio - 1
            # Monte Carlo through the asymmetric sigma errors and the lensing error
            draws = []
            for _ in range(20000):
                m = 0.0
                for s, up, dn in g['sigma']:
                    e = up if rng.random() < 0.5 else dn
                    sd = abs(rng.normal(0, e))
                    sv = max(s + (sd if e == up else -sd), 20.0) / f
                    m += np.pi * sv ** 2 * APERTURE_KPC / G / 1e14
                wl = rng.normal(*g['MWL'])
                draws.append(2 * wl / m - 1)
            d = np.array(draws)
            per_bias[f'{f:.1f}'] = dict(Mdyn_2Mpc_1e14=float(Mdyn), ratio=float(ratio), eta=float(eta),
                                        eta_16_50_84=[float(np.percentile(d, q)) for q in (16, 50, 84)],
                                        P_eta_below_1=float(np.mean(d < 1.0)))
        rows.append(dict(name=g['name'], z=g['z'], MWL_2Mpc_1e14=g['MWL'], by_sigma_bias=per_bias))

    payload = dict(
        task='T1.3', source='Munoz et al. 2013, A&A, arXiv:1212.2624, Tables 2 and 3',
        correction='An earlier revision attributed to this paper a finding that weak-lensing '
                   'masses run about 50 percent above dynamical ones. The paper contains no such '
                   'statement; the figure came from a search-engine summary. The paper concludes '
                   'instead that group velocity dispersions are always underestimated.',
        method='Isothermal projected mass pi sigma^2 R/G in the 2 Mpc lensing aperture; lensing '
               'sees (1+eta)/2 of it where the extra gravity dominates.',
        usable_groups=len(GROUPS), excluded=EXCLUDED, rows=rows,
        verdict='Taken at face value all three groups point to eta above 1, at eta = 3.1 to 4.7, '
                'and SL2SJ08544-0121 excludes eta = 1 with P = 0.05. That direction matches the '
                'three under-bent SLACS lenses. But the source paper concludes group velocity '
                'dispersions are always underestimated, and a 30 percent underestimate brings every '
                'group inside eta = 1 at 68 percent (a 20 percent one leaves SL2SJ08544-0121 just '
                'outside). With three usable groups and weak-lensing errors of 40 to 90 percent, '
                'the data are suggestive of eta above 1 and cannot establish it.',
        aperture_note='The isothermal aperture correction is not an arbitrary choice here: the law '
                      'under test makes the far field isothermal, g proportional to 1/r, so the '
                      'correction is the theory applied to itself.',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'sl2s-slip.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')

    print('T1.3 - SL2S groups, from the published table\n')
    print(f"{'group':<18}{'M_WL':>12}{'sigma bias':>12}{'M_dyn(2Mpc)':>13}{'eta':>7}{'68% range':>18}{'P(eta<1)':>10}")
    for r in rows:
        for f, v in r['by_sigma_bias'].items():
            lo, mid, hi = v['eta_16_50_84']
            print(f"{r['name']:<18}{r['MWL_2Mpc_1e14'][0]:>6.1f}+-{r['MWL_2Mpc_1e14'][1]:<4.1f}"
                  f"{f:>11}{v['Mdyn_2Mpc_1e14']:>13.2f}{v['eta']:>7.2f}"
                  f"{f'[{lo:.1f}, {hi:.1f}]':>18}{v['P_eta_below_1']:>10.2f}")
        print()
    print('excluded:', ', '.join(EXCLUDED))
    print('\nVERDICT:', payload['verdict'])
    print('WROTE', out / 'sl2s-slip.json')


if __name__ == '__main__':
    main()
