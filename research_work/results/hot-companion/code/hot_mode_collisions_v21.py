"""Round 21, step 2: the crossing heat, and the separate hot-mode speed, in MACS J0025.4-1222 and El Gordo.

A proposal supplied to the project (25 September 2026) found that letting the collision's heat spread at its own speed
v_h of about 600 km/s puts both halves of the Bullet Cluster inside their lensing masses (step 1,
code/hot_mode_speed_v21.py, checks it). It asks that the same rule survive the other colliding clusters. Round 16's
crossing heat was only ever applied to the Bullet; the suite's far collisions (code/collisions_v10.py, graded in
regression/t_new_collisions.py) have none. Here both are added to MACS J0025 and El Gordo, the two with a clean
two-body history:

  'none'     the suite's model (no crossing heat);
  'u'        round 16's crossing heat, spreading at the law's u (169.4 km/s);
  'spread'   the same heat spreading at v_h (the proposal as stated);
  'energy'   spreading at v_h with its energy booked: weight x u / v_h (see code/hot_mode_speed_v21.py).

Each system's stars take the other system's companion share p from the two systems' own spherical models before the
collision (round 16's rule, code/crossing_heat_v16.py); gas takes none. The histories are straight passes as for the
Bullet:
  MACS J0025: the suite's 0.3 Gyr since closest approach (round 12's dynamical clocks), so the receding speed is the
              galaxies' separation over that time; approaching at 2,000 km/s (Bradac et al. 2008);
  El Gordo:   0.46 Gyr since pericentre, outgoing (Ng et al. 2015; the suite's case), receding at the separation over
              that time; approaching at the pericentre speed, 2,400 km/s.
The impact parameter is the Bullet's 150 kpc (not measured for either). Abell 520's clumps have no agreed history,
so it is left out.

    python code/hot_mode_collisions_v21.py --output run-hot-mode-v21/hot_mode_collisions_v21.json [--v-h 600]
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
import bullet_v4 as V                             # noqa: E402
import collisions_v8 as V8                        # noqa: E402
import collisions_v10 as C10                      # noqa: E402
import crossing_heat_v16 as CH                    # noqa: E402
import hot_mode_speed_v21 as HM                   # noqa: E402

KPC_PER_KMS_GYR = 1.0227121650537077
CASES = dict(macs0025=dict(t_gyr=0.3, v_in=2000.0, b=150.0, pair=('se_gal', 'nw_gal')),
             el_gordo=dict(t_gyr=0.46, v_in=2400.0, b=150.0, pair=('nw_gal', 'se_gal')))


def crossing_for(spec, law, fs, case):
    """The crossing list for kappa_map_cross: for each system's stars, four shells about their centre, each with the
    kernel k(d) from the other system's companion share and the straight-pass history."""
    u = law['u_kms']; reach = u * KPC_PER_KMS_GYR * 13.0
    pos, cur = spec['pos'], spec['current']
    systems = []
    for sb in spec['subs']:
        g = V8.beta(spec['gas_total'] * sb['gas_share'], sb['rc_pre'], sb['rt_pre'], sb['galaxies'])
        stars = [cur[k] for k in sb['stars']]
        I = CH.intensity_profile([g], [dict(s) for s in stars], law, rmax=6000.0 * max(fs, 1.0))
        systems.append(dict(name=sb['name'], stars=stars, I=I, centre=np.asarray(pos[sb['galaxies']])))
    D_now = float(np.linalg.norm(systems[0]['centre'] - systems[1]['centre']))
    b = case['b'] * fs
    v_out = np.sqrt(max(D_now ** 2 - b ** 2, 0.0)) / (case['t_gyr'] * KPC_PER_KMS_GYR)
    info = dict(D_now_kpc=D_now, impact_parameter_kpc=b, v_out_kms=float(v_out), v_in_kms=case['v_in'],
                t_since_pericentre_gyr=case['t_gyr'], reach_kpc=reach, shells=[])
    edges = [e * fs for e in (0.0, 50.0, 150.0, 400.0, 1e5)]
    crossing = []
    for A, B in ((systems[0], systems[1]), (systems[1], systems[0])):
        for r0, r1 in zip(edges[:-1], edges[1:]):
            rm = np.sqrt(r0 * r1) if r0 > 0 else r1 / 2
            I_own = float(np.exp(np.interp(np.log(rm), np.log(A['I'][0]), np.log(A['I'][1]))))
            kfun = HM.make_kernel_vh(u, B['I'], I_own, D_now, v_out, case['v_in'], reach, kmax_d=None, b=b)
            dd = np.array([10.0, 30.0, 100.0, 300.0]) * fs
            info['shells'].append(dict(system=A['name'], shell_kpc=(r0, r1), k_at_d={f'{d:.0f}': float(kfun(np.array([d]))[0]) for d in dd}))
            crossing.append(([dict(s) for s in A['stars']], (r0, r1), kfun))
    return crossing, info


def solve_cross(spec, law, t_gyr, fs, crossing, n=192):
    """collisions_v10.solve with kappa_map_cross in place of kappa_map_v4 (identical when crossing is empty)."""
    dx = spec['dx']; pos, cur = spec['pos'], spec['current']
    ghost_gas, ghost_stars, speeds = {}, [], {}
    r = np.geomspace(1.0, 3000.0 * fs, 500); dr = np.gradient(r)
    for sb in spec['subs']:
        g = V8.beta(spec['gas_total'] * sb['gas_share'], sb['rc_pre'], sb['rt_pre'], sb['galaxies'])
        ghost_gas[f"gas_{sb['name'].lower()}_pre"] = g
        stars = [cur[k] for k in sb['stars']]
        rr, ss, gg = V.own_sigma_multi([g], stars, law, rmax=3000.0 * max(fs, 1.0))
        for st in stars:
            ghost_stars.append((dict(st), (rr, ss)))
        dms = sum(V.shells_on(r, dr, st) for st in stars)
        speeds[sb['name']] = {f'{int(Ra)}kpc': C10.BM.sigma_los_aperture(r, dms, np.interp(r, rr, ss), 0.0, Ra * fs) for Ra in (500.0, 1000.0, 1500.0)}
    fresh = law['u_kms'] * t_gyr * KPC_PER_KMS_GYR
    rep = {}
    x, y, Se, Sb = CH.kappa_map_cross(cur, ghost_gas, ghost_stars, pos, law, crossing=crossing, n=n, dx=dx, centre=spec['centre'],
                                      fresh_kpc=fresh, heat=True, memory=True, report=rep)
    return dict(x=x, y=y, Se=Se, Sb=Sb, fresh_kpc=fresh, speeds=speeds, dx=dx, report=rep)


def graded(name, m, f):
    """The suite's lensing checks (regression/t_new_collisions.py), with their targets in the static law."""
    L = C10.lens_factor(name)['lens']
    rows = {}
    if name == 'macs0025':
        for w, (o, lo, hi) in (('se', (2.5e14, 1.7e14, 1.0e14)), ('nw', (2.6e14, 1.4e14, 0.5e14))):
            v = m[f'M300_{w}']; o, lo, hi = o * L, lo * L, hi * L
            rows[f'M300_{w}'] = dict(model=v, target=o, minus=lo, plus=hi, z=(v - o) / (hi if v > o else lo))
            p = m[f'peak_{w}']
            rows[f'peak_{w}'] = dict(to_galaxies=p['to_galaxies'], galaxies_to_gas=p['galaxies_to_gas'])
    else:
        conv = C10.CONVENTIONS['el_gordo']
        for R, M in zip(conv['aperture_Mpc'][:2], conv['aperture_M'][:2]):
            k = f'Mkim{int(R * 1000)}_com'; v, o = m[k], M * L
            rows[k] = dict(model=v, target=o, z=(v - o) / (conv['aperture_err'] * o))
        rows['peak_se'] = m['peak_se']
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--v-h', type=float, default=600.0)
    ap.add_argument('--clusters', type=str, default='macs0025,el_gordo')
    ap.add_argument('--versions', type=str, default='none,u,spread,energy')
    ap.add_argument('--n', type=int, default=192)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    out = dict(experiment='round 21: the crossing heat and the hot-mode speed in MACS J0025 and El Gordo', law=law['name'],
               v_h=args.v_h, grid_n=args.n, cases=CASES, runs=[])
    for name in args.clusters.split(','):
        spec0 = dict(macs0025=V8.macs0025, el_gordo=V8.el_gordo)[name]()
        f = C10.factors(spec0['z'], 1.0)
        spec = C10.rescale(spec0, f, star_extra=C10.chabrier_basis(name))
        fs = f['size']; case = CASES[name]
        for version in args.versions.split(','):
            if version == 'none':
                crossing, info = [], {}
            else:
                HM.SETTINGS['v_h'] = None if version == 'u' else args.v_h
                HM.SETTINGS['energy'] = (version == 'energy')
                crossing, info = crossing_for(spec, law, fs, case)
            sol = solve_cross(spec, law, case['t_gyr'], fs, crossing, n=args.n)
            m = C10.measure(name, spec, sol, fs)
            row = dict(cluster=name, version=version, v_h=(law['u_kms'] if version == 'u' else (args.v_h if version != 'none' else None)),
                       graded=graded(name, m, f), crossing=info, S_cross_share=sol['report'].get('S_cross_over_S_mid_plane'))
            out['runs'].append(row)
            g = row['graded']
            if name == 'macs0025':
                txt = ', '.join(f"{w.upper()} {g[f'M300_{w}']['model']:.3e} (target {g[f'M300_{w}']['target']:.2e}, z {g[f'M300_{w}']['z']:+.2f})" for w in ('se', 'nw'))
                txt += '; peaks ' + ', '.join(f"{g[f'peak_{w}']['to_galaxies']:.0f}" for w in ('se', 'nw')) + ' kpc from the galaxies'
            else:
                txt = ', '.join(f"{k} {g[k]['model']:.3e} (target {g[k]['target']:.2e}, z {g[k]['z']:+.2f})" for k in ('Mkim500_com', 'Mkim1000_com'))
            share = row['S_cross_share']
            print(f"[{time.monotonic() - t0:5.0f} s] {name:9s} {version:6s}: {txt}" + (f"; crossing share {share:.2f}" if share is not None else ''), flush=True)
            HM.SETTINGS['v_h'] = None; HM.SETTINGS['energy'] = False
            args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
