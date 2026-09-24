"""Round 12: what sets the common level of galaxy lensing (KiDS-1000) against the law, in our distances.

With the lenses' heat measured (kids_heat_v12) the ellipticals' extra lensing matches, and every sample sits
0.01-0.08 dex above the law. This takes the adopted law (round 12: the distance scale fitted jointly) and changes one
ingredient at a time: the geometry that turns angles into sizes, the lenses' brightness distance (their star
masses), their angular-size distance, the stars' mass scale, and a circumgalactic gas halo modelled as Brouwer et
al. 2021 model theirs (isothermal, out to 100 kpc).

    python code/kids_level_v12.py --output run-distance-scale-v12/kids_level_v12.json
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import common                                     # noqa: E402,F401  (import paths)
import collisions_v10 as C10                      # noqa: E402
import kids_static_v11 as KS                      # noqa: E402
import kids_heat_v12 as KH                        # noqa: E402
from law_config import load_law                   # noqa: E402

SAMPLES = ('all', 'blue', 'red', 'disc', 'bulge', 'gama')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--law', default='round12')
    args = ap.parse_args()
    law = load_law(args.law)
    consts = dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI']); u, reach = law['u_kms'], law['reach_kpc']
    C10.ALPHA = law['alpha_per_Mpc']
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    rows = []

    def run(tag, ff, **kw):
        k = KH.kids(consts, u, reach, ff, **kw)
        rows.append(dict(case=tag, **{s: k[s] for s in SAMPLES}, gap=k['gap_model'], gap_sersic=k['gap_sersic_model']))
        print(f"{tag:55s} " + ' '.join(f"{s} {k[s]:+.3f}" for s in SAMPLES), flush=True)

    run('adopted (fixed geometry, D_A = D)', f)
    C10.VARIANT = 'metric'
    run('metric geometry (D_A = D/(1+z))', C10.factors(0.25, 0.75, KS.WMAP9))
    C10.VARIANT = 'fixed'
    run('brightness distance 5% longer (star masses x1.1025)', dict(f, stars=f['stars'] * 1.05 ** 2))
    run('angular-size distance 5% longer (sizes and lens masses)', dict(f, size=f['size'] * 1.05, lens=f['lens'] * 1.05))
    run('the paper\'s own (flat LCDM, WMAP9) distances', None)
    for dm in (0.05, 0.10, 0.13, 0.20):
        run(f'star masses +{dm:.2f} dex', dict(f, stars=f['stars'] * 10 ** dm))
    for fg in (0.5, 1.0, 2.0, 3.0):
        run(f'circumgalactic gas {fg:g} x the stars within 100 kpc', f, gas_frac=fg)
    run('gas 1 x the stars and star masses +0.05 dex', dict(f, stars=f['stars'] * 10 ** 0.05), gas_frac=1.0)
    res = dict(experiment='round 12: the anatomy of the KiDS lensing level in our distances', law=law['name'],
               alpha_per_Mpc=law['alpha_per_Mpc'], factors={k: float(v) for k, v in f.items() if isinstance(v, float)}, rows=rows)
    args.output.write_text(json.dumps(res, indent=1) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
