"""Round 10: how heavy must the far clusters' stars be, in the project's static distance law?

The three collisions of collisions_v10.py (MACS J0025.4-1222, Abell 520, El Gordo) rerun with every star
mass multiplied by one common factor, on top of the static-law conversion of the published values.
Compared with the measurements converted the same way as in collisions_v10.py.

    python code/collisions_star_sweep_v10.py --output run-collisions-v10/star_sweep_v10.json
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import collisions_v10 as C                        # noqa: E402
import collisions_v8 as V8                        # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--factors', default='1.4,1.8')
    args = ap.parse_args(); t0 = time.monotonic()
    law = json.loads((HERE.parent / 'run-v3/results.json').read_text())['constants']
    base = dict(macs0025=V8.macs0025(), abell520=V8.abell520(), el_gordo=V8.el_gordo())
    t_main = dict(macs0025=0.5, abell520=0.3, el_gordo=0.46)
    girardi = dict(P1=(811.0, 71.0, 278.0), P2=(749.0, 88.0, 186.0), P4=(579.0, 151.0, 523.0), P5=(668.0, 187.0, 570.0))
    out = dict(experiment='round 10: the far collisions with heavier stars (static distance law)', runs={})
    for sx in (float(x) for x in args.factors.split(',')):
        for name, spec in base.items():
            f = C.factors(spec['z'], 1.0)
            st = C.rescale(spec, f, star_extra=sx)
            sol = C.solve(st, law, t_main[name], f['size'])
            m = C.measure(name, st, sol, f['size'])
            rows = C.compare(name, m, f, V8.observed_el_gordo(spec) if name == 'el_gordo' else None)
            entry = dict(star_factor=sx, comparison=rows, speeds=m['speeds'],
                         peaks={k: v for k, v in m.items() if k.startswith('peak_')})
            if name == 'abell520':
                zs = [(m['speeds'][k]['500kpc'] - o) / (hi if m['speeds'][k]['500kpc'] > o else lo) for k, (o, lo, hi) in girardi.items()]
                entry['speeds_rms_z'] = float(np.sqrt(np.mean(np.square(zs))))
            out['runs'][f'{name}_x{sx:g}'] = entry
            print(f'{name}, stars x{sx:g}:', flush=True)
            for r in rows:
                meas = r.get('measured', r.get('measured_range'))
                ms = f"{meas:.3g}" if isinstance(meas, float) else f"{meas[0]:.3g}-{meas[1]:.3g}"
                print(f"   {r['check']:24s} model {r['model']:.3g}  measured {ms}  z {r['z']:+.2f}", flush=True)
            for k, v in entry['peaks'].items():
                print(f"   {k}: lensing peak {v['to_galaxies']:.0f} kpc from its galaxies", flush=True)
            if name == 'abell520':
                print(f"   galaxy speeds at P1, P2, P4, P5: {[round(m['speeds'][k]['500kpc']) for k in girardi]} km/s, rms z {entry['speeds_rms_z']:.2f}", flush=True)
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
