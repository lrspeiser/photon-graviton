"""Round 19, step 2c: the law the medium produces, with its exponents and their uncertainties.

Round 16's one kind of matter (code/one_matter_v17.py, structure 'single'), sources of 24, 48 and 96 pieces at the
density of round 16 (radius 3 (N/48)^(1/3)), 12 receivers of the same matter on each of two spheres (r = 9 and 18,
turned so that no two share a line of sight), cold, warm (k = 2, 8) and colliding sources, and sources of 48 pieces
squeezed or spread (radius 1, 2 and 4.5; radius 1 is about a wavelength, the sub-wavelength side the clusters point
to: code/hot_shell_v19.py allows only weak absorption, so a warm source must keep in step some other way). The medium: the two-way wave, and the stream that absorbs inward-travelling
waves at kappa = 1 per wavelength in its ray form (code/absorbing_stream_v18.py), the strength the X-COP clusters allow
(code/hot_shell_v19.py: absorption over a few tenths of the cluster radius; a source of radius 3 wavelengths is then a
few absorption lengths across) and where code/full_wave_v19.py finds the ray form close to the exact medium.

From the runs:  pull = A M^p r^-q  (cold, and warm at k = 8), fitted over sources and receivers with the seeds
resampled; the warm gain R = pull(k)/pull(0) against the law's sqrt(1 + k); density at fixed mass.

    python code/emergent_law_v19.py --output run-emergent-law-v19/emergent_law_v19.json [--processes 4]
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from one_matter_v17 import case, shells, simulate            # noqa: E402

RADII = (9.0, 18.0)


def configs(kappas=(0.0, 1.0), seeds=(1, 2, 3)):
    runs = []
    rec = shells(RADII, per=12)
    for kap in kappas:
        for sd in seeds:
            base = dict(receivers=rec, T=8000.0, burn=4000.0, momentum=False, seed=sd)
            if kap:
                base['absorb'] = kap
            for Ns in (24, 48, 96):
                kw = dict(base, Ns=Ns, Rb=3.0 * (Ns / 48) ** (1 / 3))
                runs.append(case('single', 'cold', **kw))
                runs.append(case('single', 'free', k=8.0, **kw))
                if Ns == 48:
                    runs.append(case('single', 'free', k=2.0, **kw))
                    runs.append(case('single', 'collisional', k=8.0, nu=50.0, **kw))
            for Rb in (1.0, 2.0, 4.5):
                kw = dict(base, Ns=48, Rb=Rb)
                runs.append(case('single', 'cold', **kw))
                runs.append(case('single', 'free', k=8.0, **kw))
    return runs


def tagof(c):
    return f"{c['tag']}{'' if c['tag'] == 'cold' else ' k' + format(c['k'], 'g')}{' nu' + format(c['nu'], 'g') if c.get('nu') else ''}"


def analyse(res, nboot=2000, rng_seed=7):
    """Pull per (medium, kind, Ns, Rb, radius, seed): the mean over the 12 receivers at that radius. Fits log pull =
    log A + p log Ns - q log r over the fixed-density sources, per medium and kind, with the seeds resampled; warm gains;
    density ratios. Negative mean pulls (a push) are reported and left out of the log fits."""
    rows = []
    for r in res:
        c = r['cfg']; d = np.array(r['receiver_distance']); p = np.array(r['pull_each'])
        for rad in RADII:
            m = np.isclose(d, rad, atol=1e-2)
            rows.append(dict(kappa=c.get('absorb', 0.0) or 0.0, kind=tagof(c), Ns=c['Ns'], Rb=round(c['Rb'], 3), r=rad, seed=c['seed'],
                             pull=float(p[m].mean()), pull_sd=float(p[m].std(ddof=1)), n=int(m.sum()),
                             keeping_step=float(r['lead_src']), output=float(r['src_out']), rhythm_spread=float(r['rhythm_spread_sources'])))
    out = dict(rows=rows, fits={}, warm_gain={}, density={})
    rng = np.random.default_rng(rng_seed)
    fixed = lambda row: abs(row['Rb'] - 3.0 * (row['Ns'] / 48) ** (1 / 3)) < 1e-2
    for kap in sorted({x['kappa'] for x in rows}):
        for kind in ('cold', 'free k8'):
            sel = [x for x in rows if x['kappa'] == kap and x['kind'] == kind and fixed(x)]
            seeds = sorted({x['seed'] for x in sel})
            def fit(sub):
                sub = [x for x in sub if x['pull'] > 0]
                if len(sub) < 4:
                    return None
                A = np.array([[1.0, np.log(x['Ns']), -np.log(x['r'])] for x in sub]); y = np.log([x['pull'] for x in sub])
                c, *_ = np.linalg.lstsq(A, y, rcond=None); return c
            c0 = fit(sel)
            boots = []
            for _ in range(nboot):
                pick = rng.choice(seeds, len(seeds))
                sub = [x for s in pick for x in sel if x['seed'] == s]
                cb = fit(sub)
                if cb is not None:
                    boots.append(cb)
            boots = np.array(boots)
            npush = sum(1 for x in sel if x['pull'] <= 0)
            out['fits'][f'kappa {kap:g}, {kind}'] = dict(
                p=float(c0[1]) if c0 is not None else None, q=float(c0[2]) if c0 is not None else None,
                p_se=float(boots[:, 1].std()) if len(boots) else None, q_se=float(boots[:, 2].std()) if len(boots) else None,
                points=len(sel), pushes=npush, law=dict(p=0.5, q=1.0))
        # the warm gain at 48 pieces, per radius and seed, against sqrt(1 + k)
        for kind, k in (('free k2', 2.0), ('free k8', 8.0), ('collisional k8 nu50', 8.0)):
            g = []
            for rad in RADII:
                for sd in sorted({x['seed'] for x in rows}):
                    cold = [x for x in rows if x['kappa'] == kap and x['kind'] == 'cold' and x['Ns'] == 48 and fixed(x) and x['r'] == rad and x['seed'] == sd]
                    warm = [x for x in rows if x['kappa'] == kap and x['kind'] == kind and x['Ns'] == 48 and fixed(x) and x['r'] == rad and x['seed'] == sd]
                    if cold and warm and cold[0]['pull'] > 0:
                        g.append(dict(r=rad, seed=sd, gain=warm[0]['pull'] / cold[0]['pull']))
            if g:
                gv = np.array([x['gain'] for x in g])
                out['warm_gain'][f'kappa {kap:g}, {kind}'] = dict(mean=float(gv.mean()), se=float(gv.std(ddof=1) / np.sqrt(len(gv))) if len(gv) > 1 else None,
                                                                  law=float(np.sqrt(1 + k)) if 'collisional' not in kind else 1.0, each=g)
        # density: 48 pieces at radius 2, 3, 4.5
        for kind in ('cold', 'free k8'):
            d = {}
            for Rb in (1.0, 2.0, 3.0, 4.5):
                v = [x['pull'] for x in rows if x['kappa'] == kap and x['kind'] == kind and x['Ns'] == 48 and abs(x['Rb'] - Rb) < 1e-2]
                if v:
                    d[f'radius {Rb:g}'] = dict(mean=float(np.mean(v)), se=float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None)
            out['density'][f'kappa {kap:g}, {kind}'] = d
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4); ap.add_argument('--analyse-only', action='store_true')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    if args.analyse_only:
        res = json.loads(args.output.read_text())['runs']
    else:
        runs = configs()
        # the largest first, so the pool ends evenly
        runs.sort(key=lambda c: -c['Ns'])
        res = []
        with Pool(args.processes) as pool:
            for r in pool.imap_unordered(simulate, runs):
                res.append(r); c = r['cfg']
                print(f"[{time.monotonic() - t0:5.0f} s] kappa {c.get('absorb', 0) or 0:3.1f} {tagof(c):20s} Ns {c['Ns']:3d} Rb {c['Rb']:.2f} seed {c['seed']}: "
                      f"pull {r['pull']:+.3e}, keeping step {r['lead_src']:+.2f}, rhythm spread {r['rhythm_spread_sources']:.1e}, E res {r['energy']['residual']:.1e}", flush=True)
                args.output.write_text(json.dumps(dict(experiment='round 19: the emergent law (mass, distance, heat, density)', runs=res,
                                                       seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    summ = analyse(res)
    args.output.write_text(json.dumps(dict(experiment='round 19: the emergent law (mass, distance, heat, density)', summary=summ, runs=res,
                                           seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    for k, v in summ['fits'].items():
        print(f"{k:22s}: p = {v['p']} +- {v['p_se']}, q = {v['q']} +- {v['q_se']} ({v['points']} points, {v['pushes']} pushes)")
    for k, v in summ['warm_gain'].items():
        print(f"{k:32s}: gain {v['mean']:.2f} +- {v['se']} (law {v['law']:.2f})")
    for k, v in summ['density'].items():
        print(f"{k:22s}: " + ', '.join(f"{a} {b['mean']:.2e}" for a, b in v.items()))


if __name__ == '__main__':
    main()
