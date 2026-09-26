"""Round 24 companion page, "Gravity's currents": the real numbers behind each scene, from the project's own code and
the adopted law (round 12). Every value the page plots comes from here.

  * a spiral and a faint galaxy from SPARC: measured speeds, Newton's speeds from the visible matter, the law's speeds;
  * one X-COP cluster: hydrostatic masses, the law, the law with the hot glow added as arrows (the Casimir-EFT note's
    local source, round 24), and Newton;
  * the Bullet Cluster: the law's lensing map (the suite's 3D calculation, bullet_v4.kappa_map_v4), the same map without
    the companion's memory, and the visible matter's map, with the galaxies' and gas's positions;
  * the ten dwarf galaxies (measured, the law, the law without the external hold), the KiDS lensing samples, and the
    suite's scoreboard by group.

    python code/currents_data_v24.py --output run-currents-v24/currents_data.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))

GALAXIES = ('NGC3198', 'DDO154')
CLUSTER = 'A1795'


def sparc_part(law, ctx):
    import run as R
    a, lam, u = law['a_code'], law['lam'], law['u_kms']
    out = {}
    for g in ctx.sparc():
        if g['name'] not in GALAXIES:
            continue
        gl = R.galaxy_g(g, a, u, lam)
        out[g['name']] = dict(r_kpc=g['r'].round(3).tolist(), v_obs=g['v'].round(2).tolist(), err=g['err'].round(2).tolist(),
                              v_newton=np.sqrt(g['r'] * g['gN']).round(2).tolist(), v_law=np.sqrt(g['r'] * gl).round(2).tolist())
    return out


def xcop_part(law, ctx):
    import law as L
    import t_clusters as TC
    import run_v3 as R3
    a, lam, u = law['a_code'], law['lam'], law['u_kms']
    G = L.G
    cls = TC.static_clusters(ctx)
    res = {c['name']: float(np.sqrt(np.mean(np.log(c['Mh'] / R3.cluster_M3(c, a, u, lam)) ** 2))) for c in cls}
    c = next(c for c in cls if c['name'] == CLUSTER)
    R = np.geomspace(0.07, 1.25, 40) * c['R5']
    Mb = np.array([np.sum((c['dmg'] + c['dms'])[c['s'] < r]) for r in R])
    k = L.heat_weight(np.sqrt(c['sig2_star_hse']), u) * c['dms']
    gN = G * Mb / R ** 2
    geo0 = L.HOT_GEOMETRY
    prof = {}
    for name, geo in (('law', 'two_way'), ('arrows', 'one_way_vector')):
        L.HOT_GEOMETRY = geo
        S = G * (L.shell_weights(R, c['s']) @ k) / R ** 2
        prof[name] = (L.total(gN, S, a, lam) * R ** 2 / G).tolist()
    L.HOT_GEOMETRY = geo0
    return dict(name=CLUSTER, R500_kpc=float(c['R5']), rms_all=res, r_kpc=R.round(1).tolist(), M_law=prof['law'],
                M_arrows=prof['arrows'], M_newton=Mb.tolist(), data=dict(r_kpc=c['Rk'].round(1).tolist(), M=c['Mh'].tolist(),
                                                                         err=c['eMh'].tolist()),
                note='static distances; stars deprojected (the suite\'s own sample); M in solar masses')


def bullet_part(law):
    import bullet_v4 as V
    import bullet_static_v11 as BS
    import collisions_v10 as C10
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
    got = {}
    orig = V.kappa_map_v4

    def capture(*a, **k):
        res = orig(*a, **k)
        got['args'], got['kw'], got['res'] = a, dict(k), res
        return res
    V.kappa_map_v4 = capture
    try:
        r = BS.run_bullet(law, dict(f), n=192, dx=15.0)
    finally:
        V.kappa_map_v4 = orig
    x, y, Se, Sb = got['res']
    kw = dict(got['kw']); kw['memory'] = False
    _, _, Se0, _ = orig(*got['args'], **kw)
    scrit = r['sigma_crit']
    pos = got['args'][3]

    ix = np.where((x >= -700) & (x <= 1550))[0]; iy = np.where((y >= -700) & (y <= 950))[0]   # around the two clusters
    crop = lambda m: (np.asarray(m, float) / scrit)[np.ix_(ix, iy)]
    maps = {k: crop(v) for k, v in (('law', Se), ('no_memory', Se0), ('visible', Sb))}
    q = lambda m: np.clip(np.round(m.T * 1000), 0, 65535).astype(int).ravel().tolist()   # rows = y (south to north), columns = x
    ch = r['checks']
    cur, ghost_stars = got['args'][0], got['args'][2]
    comps = {k: dict(kind=v['kind'], M=float(v['M']), scale=float(v['scale']), centre=v['centre']) for k, v in cur.items()}
    sig = {}
    for comp, prof in ghost_stars:
        rr, ss = prof
        sig[comp['centre']] = float(np.interp(100.0, rr, ss))
    return dict(x_kpc=[float(x[ix[0]]), float(x[ix[-1]])], y_kpc=[float(y[iy[0]]), float(y[iy[-1]])], nx=int(len(ix)), ny=int(len(iy)),
                kappa_milli={k: q(v) for k, v in maps.items()},
                positions={k: [float(v[0]), float(v[1])] for k, v in pos.items()},
                components=comps, star_sigma_100kpc_kms=sig,
                checks={k: dict(value=(v['value'] if not isinstance(v['value'], (np.floating,)) else float(v['value'])),
                                target=v.get('target')) for k, v in ch.items()},
                note='kappa x 1000 on the suite\'s own 3D grid (15 kpc x the static size factor: 18.1 kpc cells), cropped to the two clusters; x = west, y = north, '
                     'kpc in the static distances; main galaxies at the origin')


def dwarfs_part():
    base = {c['id']: c for c in json.loads((HERE.parent / 'regression/baseline.json').read_text())['checks']}
    nh = {c['id']: c for c in json.loads((HERE.parent / 'run-no-hold-v20/no_hold_r12.json').read_text())['checks']}
    dw = json.loads((HERE.parent / 'data/mw_dwarfs.json').read_text())['dwarfs']
    out = []
    for d in dw:
        key = 'dwarfs.' + d['name'].lower().replace(' ', '_')
        out.append(dict(name=d['name'], D_gc_kpc=d['D_gc_kpc'], sigma_obs=d['sigma_obs'], sigma_err=d['sigma_err'],
                        law=round(base[key]['value'], 2), law_status=base[key]['status'],
                        no_hold=round(nh[key]['value'], 2), no_hold_status=nh[key]['status']))
    return out


def suite_part():
    b = json.loads((HERE.parent / 'regression/baseline.json').read_text())
    groups = {}
    items = []
    for c in b['checks']:
        if c['status'] not in ('pass', 'close', 'fail'):
            continue
        g = groups.setdefault(c['group'], dict(pass_=0, close=0, fail=0))
        g['pass_' if c['status'] == 'pass' else c['status']] += 1
        if c['status'] != 'pass':
            items.append(dict(id=c['id'], group=c['group'], status=c['status'], title=c['title'], value=c['value'],
                              unit=c.get('unit', ''), target=c['target']))
    kids = {c['id'].split('.')[1]: round(c['value'], 4) for c in b['checks'] if c['id'].startswith('lensing.kids_') and 'gap' not in c['id']}
    gaps = {c['id'].split('.')[1]: round(c['value'], 4) for c in b['checks'] if 'kids_gap' in c['id']}
    return dict(groups=groups, not_passing=items, kids_dex=kids, kids_gaps_dex=gaps, law=b['law']['name'], date=b['date'])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    out = dict(source='code/currents_data_v24.py', law=law['name'],
               constants=dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'], u_kms=law['u_kms'], release_length_au=law['release_length_au']))
    out['sparc'] = sparc_part(law, ctx); print(f'[{time.monotonic() - t0:4.0f} s] SPARC', flush=True)
    out['xcop'] = xcop_part(law, ctx); print(f'[{time.monotonic() - t0:4.0f} s] X-COP', flush=True)
    out['dwarfs'] = dwarfs_part()
    out['suite'] = suite_part()
    out['bullet'] = bullet_part(law); print(f'[{time.monotonic() - t0:4.0f} s] Bullet', flush=True)
    args.output.write_text(json.dumps(out, default=float) + '\n')
    print(f'wrote {args.output} ({args.output.stat().st_size / 1e3:.0f} kB, {time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
