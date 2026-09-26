"""Round 23, step 2: the beat-driven heat of crossing in MACS J0025.4-1222 and El Gordo.

The same test as code/hot_mode_collisions_v21.py (round 21, section 32.2: the suite's far-collision models with round
16's heat of crossing added, straight-pass histories), with the share p of the other system's flow replaced by the
depth of the two companions' beat, m = 2 sqrt(I_other I_own)/(I_other + I_own) (code/beat_heat_v23.py). Cases: 'share'
(round 16's rule, spreading at u), 'beat' (spreading at u), 'beat:600' (spreading at 600 km/s, energy booked).

    python code/beat_heat_far_v23.py --output run-beat-heat-v23/beat_far.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import argparse, json, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import collisions_v8 as V8                        # noqa: E402
import collisions_v10 as C10                      # noqa: E402
import hot_mode_collisions_v21 as HC              # noqa: E402
import beat_heat_v23 as BH                        # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--clusters', type=str, default='macs0025,el_gordo')
    ap.add_argument('--cases', type=str, default='share,beat,beat:600')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    out = dict(experiment='round 23: the beat-driven heat of crossing in MACS J0025 and El Gordo', law=law['name'], cases=HC.CASES, runs=[])
    orig = HC.HM.make_kernel_vh
    HC.HM.make_kernel_vh = BH.make_kernel_beat
    try:
        for name in args.clusters.split(','):
            spec0 = dict(macs0025=V8.macs0025, el_gordo=V8.el_gordo)[name]()
            f = C10.factors(spec0['z'], 1.0)
            spec = C10.rescale(spec0, f, star_extra=C10.chabrier_basis(name))
            fs = f['size']; case = HC.CASES[name]
            for c in args.cases.split(','):
                wname, _, vh = c.partition(':')
                BH.SETTINGS.update(weight=wname, v_h=float(vh) if vh else None, energy=True)
                crossing, info = HC.crossing_for(spec, law, fs, case)
                sol = HC.solve_cross(spec, law, case['t_gyr'], fs, crossing)
                m = C10.measure(name, spec, sol, fs)
                g = HC.graded(name, m, f)
                out['runs'].append(dict(cluster=name, case=c, graded=g, crossing=info, S_cross_share=sol['report'].get('S_cross_over_S_mid_plane')))
                if name == 'macs0025':
                    txt = ', '.join(f"{w.upper()} {g[f'M300_{w}']['model']:.3e} (z {g[f'M300_{w}']['z']:+.2f})" for w in ('se', 'nw'))
                    txt += '; peaks ' + ', '.join(f"{g[f'peak_{w}']['to_galaxies']:.0f}" for w in ('se', 'nw')) + ' kpc'
                else:
                    txt = ', '.join(f"{k} {g[k]['model']:.3e} (target {g[k]['target']:.2e}, z {g[k]['z']:+.2f})" for k in ('Mkim500_com', 'Mkim1000_com'))
                print(f"[{time.monotonic() - t0:5.0f} s] {name:9s} {c:9s}: {txt}", flush=True)
                args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    finally:
        HC.HM.make_kernel_vh = orig
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
