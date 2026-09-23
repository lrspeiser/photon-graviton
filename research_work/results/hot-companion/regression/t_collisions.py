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
    import bullet_v3 as B, bullet_v4 as V, bullet_main_v5 as BM, collisions_v4 as C4
    pos, comps0, scrit = bullet_setup(ctx)
    gas = dict(comps0['gas_main'], centre='main_bcg')
    cm = copy.deepcopy(comps0)                   # published star masses (M/L_I = 2)
    inner = dict(cm['st_main'])
    target = BM.OBS_SIGMA[0]

    def los(M_out):
        stars = [inner] + ([BM.fixed_outer(inner, M_out, 800.0)] if M_out > 0 else [])
        r, sr, g, dms = BM.own_sigma([gas], stars, law)
        return BM.sigma_los_aperture(r, dms, sr, 0.0, BM.BARRENA_RAP), (r, sr, stars)

    ctx.log('Bullet: outer stars that give the measured galaxy speed (1,249 km/s)')
    lo, hi = 11.0, 13.8
    if los(10 ** lo)[0] > target:
        M_out = 0.0
    else:
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if los(10 ** mid)[0] < target else (lo, mid)
        M_out = 10 ** (0.5 * (lo + hi))
    _, (r, sr, stars) = los(M_out)
    ghost_gas, sub_sats, pre = V.pre_collision_models(cm, 8.0, 0.07, 150.0)
    rs_s, ss_s, _ = V.own_sigma_multi([dict(cm['gas_sub']), ghost_gas['gas_atm_ghost']], [cm['st_sub'], sub_sats], law)
    cur = dict(cm); cur['st_main'] = stars[0]
    ghost_stars = [(dict(stars[0]), (r, sr))]
    if len(stars) > 1:
        cur['st_main_outer'] = stars[1]; ghost_stars.append((dict(stars[1]), (r, sr)))
    ghost_stars += [(dict(cm['st_sub']), (rs_s, ss_s))] + ([(sub_sats, (rs_s, ss_s))] if sub_sats['M'] > 0 else [])
    fresh = 30.0 * law['u_kms'] / U0
    ctx.log(f'Bullet: lensing map ({n}^3 cells of {dx:g} kpc)')
    x, y, Se, Sb = V.kappa_map_v4(cur, ghost_gas, ghost_stars, pos, law, n=n, dx=dx, fresh_kpc=fresh)
    kap = Se / scrit
    dec = B.clowe_decomposition(x, y, kap, pos)
    pk = C4.refined_peaks(x, y, kap, pos, smooth_kpc=40.0)
    mp = min(pk, key=lambda d: d['dist_main_bcg']); sp = min(pk, key=lambda d: d['dist_sub_bcg'])
    m250 = {w: BM.mass_within(x, y, kap, pos[f'{w}_bcg'], 250.0, scrit) for w in ('main', 'sub')}
    sep = {w: float(np.linalg.norm(pos[f'{w}_plasma'] - pos[f'{w}_bcg'])) for w in ('main', 'sub')}
    ref = 'Clowe et al. 2006; Bradac et al. 2006; Paraficz et al. 2016; Barrena et al. 2002'
    out = [make('bullet.outer_stars', GROUP, 'Bullet main cluster: outer stars needed for the galaxy speed', M_out,
                crit=range_check(M_out, 3.9e12, 6.7e12, 0.7e12), unit='Msun',
                target='3.9-6.7 x 10^12 (Legacy Survey star count, round 6: 85-103% of the round-5 total)', refs=ref),
           make('bullet.kappa_main', GROUP, 'Bullet main cluster: lensing strength (kappa) on its galaxies', dec['main_bcg'],
                crit=floor_check(dec['main_bcg'], 0.36, 0.06), target='at least 0.36 +- 0.06 (a stated lower bound)', refs=ref),
           make('bullet.kappa_sub', GROUP, 'Bullet subcluster: lensing strength (kappa) on its galaxies', dec['sub_bcg'],
                crit=floor_check(dec['sub_bcg'], 0.20, 0.05), target='at least 0.20 +- 0.05 (a stated lower bound)', refs=ref),
           make('bullet.gas_main', GROUP, 'Bullet: leftover lensing on the main cluster\'s gas', dec['main_plasma'],
                crit=z_check(dec['main_plasma'], 0.05, 0.06), target='0.05 +- 0.06', refs=ref),
           make('bullet.gas_sub', GROUP, 'Bullet: leftover lensing on the subcluster\'s gas', dec['sub_plasma'],
                crit=z_check(dec['sub_plasma'], 0.02, 0.06), target='0.02 +- 0.06', refs=ref)]
    for w, p in (('main', mp), ('sub', sp)):
        d = p[f'dist_{w}_bcg']
        out.append(make(f'bullet.peak_{w}', GROUP, f'Bullet {w}: lensing peak distance from its galaxies', d,
                        crit=at_most(d, 0.25 * sep[w], 0.5 * sep[w]), unit='kpc',
                        target=f'within a quarter of the galaxy-gas separation ({sep[w]:.0f} kpc)'))
    out.append(make('bullet.m250_main', GROUP, 'Bullet main cluster: lensing mass inside 250 kpc', m250['main'],
                    crit=range_check(m250['main'], 2.5e14, 2.8e14, 0.15e14), unit='Msun', target='2.5-2.8 x 10^14 (Paraficz 2016; Bradac 2006)'))
    out.append(make('bullet.m250_sub', GROUP, 'Bullet subcluster: lensing mass inside 250 kpc', m250['sub'],
                    crit=range_check(m250['sub'], 2.0e14, 2.3e14, 0.2e14), unit='Msun', target='2.0-2.3 x 10^14 (Paraficz 2016; Bradac 2006)',
                    note='the subcluster\'s size before the collision (1:8 here) is the main uncertainty'))
    ctx.shared['bullet_map'] = dict(x=x, y=y, kappa=kap)
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
