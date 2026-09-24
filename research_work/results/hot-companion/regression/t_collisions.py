"""Colliding clusters (full tier): the Bullet Cluster's published case (round 5, case L) and the
collision stack of Harvey et al. 2015 (round 4). Each lensing map solves the field equation on a
3D grid with the companion's memory (code/bullet_v4.py): the old companion rides with the
galaxies; a fresh sphere of radius u t has formed around the stopped gas."""
from __future__ import annotations
import copy
import numpy as np
from checks import make, z_check, range_check, floor_check, at_most

GROUP = 'collisions'
U0 = 197.41001228101223          # the round-3 u: the published fresh sphere, 30 kpc, is 150 Myr at this speed


def bullet_setup(ctx):
    import bullet_v3 as B
    def build():
        pos = B.positions(); comps, _ = B.fit_components(pos)
        return pos, comps, B.sigma_crit()
    return ctx.get('bullet_setup', build)


def bullet(law, ctx, n=192, dx=15.0):
    """The round-5 case in the project's static distance law (code/bullet_static_v11.py): lengths x size,
    gas x gas, stars x stars (also the outer stars' target from the star count), lensing masses inside a
    fixed angle x lens; kappa is an observable and is compared as measured. Round 10 and earlier graded the
    papers' flat-LCDM units."""
    import bullet_static_v11 as BS
    import collisions_v10 as C10
    ctx.log('Bullet: outer stars, then the lensing map (static distances)')
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
    r = BS.run_bullet(law, dict(f), n=n, dx=dx)
    ch = r['checks']; fs, fst, fl = f['size'], f['stars'], f['lens']
    ref = 'Clowe et al. 2006; Bradac et al. 2006; Paraficz et al. 2016; Barrena et al. 2002; converted to the static law'
    M_out = ch['outer_stars']['value']
    out = [make('bullet.outer_stars', GROUP, 'Bullet main cluster: outer stars needed for the galaxy speed (static distances)', M_out,
                crit=range_check(M_out, 3.9e12 * fst, 6.7e12 * fst, 0.7e12 * fst), unit='Msun',
                target=f'{3.9 * fst:.1f}-{6.7 * fst:.1f} x 10^12 (Legacy Survey star count, round 6, converted)', refs=ref),
           make('bullet.kappa_main', GROUP, 'Bullet main cluster: lensing strength (kappa) on its galaxies', ch['kappa_main']['value'],
                crit=floor_check(ch['kappa_main']['value'], 0.36, 0.06), target='at least 0.36 +- 0.06 (a stated lower bound)', refs=ref),
           make('bullet.kappa_sub', GROUP, 'Bullet subcluster: lensing strength (kappa) on its galaxies', ch['kappa_sub']['value'],
                crit=floor_check(ch['kappa_sub']['value'], 0.20, 0.05), target='at least 0.20 +- 0.05 (a stated lower bound)', refs=ref),
           make('bullet.gas_main', GROUP, 'Bullet: leftover lensing on the main cluster\'s gas', ch['gas_main']['value'],
                crit=z_check(ch['gas_main']['value'], 0.05, 0.06), target='0.05 +- 0.06', refs=ref),
           make('bullet.gas_sub', GROUP, 'Bullet: leftover lensing on the subcluster\'s gas', ch['gas_sub']['value'],
                crit=z_check(ch['gas_sub']['value'], 0.02, 0.06), target='0.02 +- 0.06', refs=ref)]
    for w in ('main', 'sub'):
        d = ch[f'peak_{w}']['value']
        lim = float(ch[f'peak_{w}']['target'].split('<=')[1].split('kpc')[0])
        out.append(make(f'bullet.peak_{w}', GROUP, f'Bullet {w}: lensing peak distance from its galaxies', d,
                        crit=at_most(d, lim, 2 * lim), unit='kpc', target=f'within a quarter of the galaxy-gas separation ({4 * lim:.0f} kpc)'))
    for w, lo, hi, e in (('main', 2.5e14, 2.8e14, 0.15e14), ('sub', 2.0e14, 2.3e14, 0.2e14)):
        v = ch[f'm250_{w}']['value']
        out.append(make(f'bullet.m250_{w}', GROUP, f'Bullet {"main cluster" if w == "main" else "subcluster"}: lensing mass inside {250 * fs:.0f} kpc (250 in LCDM units)', v,
                        crit=range_check(v, lo * fl, hi * fl, e * fl), unit='Msun',
                        target=f'{lo * fl / 1e14:.2f}-{hi * fl / 1e14:.2f} x 10^14 (Paraficz 2016; Bradac 2006; converted)',
                        note='the subcluster\'s size before the collision (1:8 here) is the main uncertainty' if w == 'sub' else ''))
    return out


def stack(law, ctx, n=160, dx=18.0, ml=1.5):
    """Round 4's family of 12 collisions (lag 0-300 kpc, subcluster x1 and x3): how far the lensing
    moves when the gas moves, beta_gas = ds/d(lag), against Harvey et al. 2015."""
    import bullet_v3 as B, bullet_v4 as V, collisions_v3 as C3, collisions_v4 as C4
    base_pos, comps0, scrit = bullet_setup(ctx)
    comps = copy.deepcopy(comps0)
    for k in ('st_main', 'st_sub'): comps[k]['M'] *= ml / 2.0
    axis = (base_pos['sub_bcg'] - base_pos['main_bcg']) / np.linalg.norm(base_pos['sub_bcg'] - base_pos['main_bcg'])
    rows = []
    for lag, ms in [(lg, m) for lg in C4.LAGS for m in (1.0, 3.0)]:
        ctx.log(f'collision stack: lag {lag} kpc, subcluster x{ms:g}')
        c, pos = C3.scenario(comps, base_pos, lag, lag, sub_mass_scale=ms)
        gg, sats, pre = V.pre_collision_models(c, ratio=8.0 / ms)
        rm, sm, _ = V.own_sigma_multi([c['gas_main']], [c['st_main']], law)
        rs, ss, _ = V.own_sigma_multi([dict(c['gas_sub']), gg['gas_atm_ghost']], [c['st_sub'], sats], law)
        gst = [(dict(c['st_main']), (rm, sm)), (dict(c['st_sub']), (rs, ss))] + ([(sats, (rs, ss))] if sats['M'] > 0 else [])
        t_myr = lag / C4.SEP_SPEED * 977.8
        fresh = law['u_kms'] * t_myr / 977.8
        x, y, Se, Sb = V.kappa_map_v4(c, gg, gst, pos, law, n=n, dx=dx, fresh_kpc=fresh)
        pk = C4.refined_peaks(x, y, Se / scrit, pos, smooth_kpc=30.0)
        row = dict(lag=lag, ms=ms)
        for side, st, sgn in (('main', 'main_bcg', 1.0), ('sub', 'sub_bcg', -1.0)):
            near = [p for p in pk if p[f'dist_{st}'] < lag + 60]
            if near:
                p = min(near, key=lambda d: d[f'dist_{st}'])
                row[side] = float(sgn * (np.array([p['x'], p['y']]) - pos[st]) @ axis)
                row[f'off_{side}'] = p[f'dist_{st}']
            else:
                row[side] = None
        rows.append(row)
    ref = {(r['ms'], sd): r[sd] for r in rows if r['lag'] == 0 for sd in ('main', 'sub')}
    beta = [(r[sd] - ref[(r['ms'], sd)]) / r['lag'] for r in rows if r['lag'] > 0 for sd in ('main', 'sub')
            if r[sd] is not None and ref[(r['ms'], sd)] is not None]
    bm = float(np.median(beta))
    offs = [r[f'off_{sd}'] for r in rows if r['lag'] > 0 for sd in ('main', 'sub') if r.get(f'off_{sd}') is not None]
    return [make('stack.beta_gas', GROUP, 'collision stack: how far the lensing follows the gas (beta, median of 20)', bm,
                 crit=z_check(bm, -0.04, 0.07), target='-0.04 +- 0.07 (Harvey et al. 2015, 72 collisions)',
                 detail=dict(n=len(beta), range=[float(min(beta)), float(max(beta))]), refs='Harvey et al. 2015, Science 347, 1462'),
            make('stack.offset', GROUP, 'collision stack: lensing peak distance from the galaxies (median)', float(np.median(offs)),
                 unit='kpc', crit=('info', None), target='raw offset; includes the pull of the other cluster')]


def run(law, ctx):
    if not ctx.full:
        return []
    out = bullet(law, ctx)
    out += stack(law, ctx)
    try:
        import t_new_collisions as NC
    except ImportError:
        return out
    return out + NC.run(law, ctx)
