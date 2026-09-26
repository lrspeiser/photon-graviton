"""Round 25: the Milky Way's rotation curve with the companion trapped (the registered candidate trapped_own_heat_r12, at
its refitted u) against the adopted law and the four Gaia analyses the suite grades, with the trapping zone and the
band just outside it (places whose n would turn negative once the trapped glow is missing). For the write-up's chart.

    python code/trapped_mw_curve_v25.py --output run-trapped-v25/mw_curve_v25.json
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
sys.path.insert(0, str(RESULTS / 'regression'))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    import common as C
    from law_config import load_law, with_constants
    import t_milky_way as TMW
    import milky_way_v7 as MW
    import trapping_v25 as TR
    run = json.loads((RESULTS / 'run-trapped-v25/suite/trapped_own_heat_r12-full/results.json').read_text())
    u_trap = run['law']['u_kms']
    out = dict(source='code/trapped_mw_curve_v25.py', model='McMillan (2017) matter, the suite\'s Milky Way (t_milky_way)', curves={})
    for label, spec, u in (('the law (round 12)', None, None), ('trapped, own heat', 'trapped_own_heat_r12', u_trap)):
        law = load_law(spec)
        if u is not None:
            law = with_constants(law, u_kms=u)
        ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
        comps, grid, F = TMW.setup(ctx)
        r = MW.Evaluator(comps, grid, F, law).run(MW.models(comps)['M17'], 'ours', reach=law['reach_kpc'])
        info = next(iter(TR.STATS.get('milky_way', {}).values()), {})
        out['curves'][label] = dict(u_kms=law['u_kms'], R=r['R'], v=r['v'], v_sun=float(np.interp(MW.R0, r['R'], r['v'])),
                                    zone_outer_R_kpc=info.get('zone_outer_R_in_plane_kpc'), band_outer_R_kpc=info.get('band_outer_R_in_plane_kpc'),
                                    held_share=info.get('held_share'))
        print(f"{label}: u {law['u_kms']:.1f}, v_sun {out['curves'][label]['v_sun']:.1f}, zone {info.get('zone_outer_R_in_plane_kpc')}, band {info.get('band_outer_R_in_plane_kpc')}")
    rc = json.loads((RESULTS / 'data/mw_rotation_curves.json').read_text())
    out['gaia'] = {k: dict(R=d['R'], v=d['v'], err=d['err']) for k, d in rc.items()}
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
