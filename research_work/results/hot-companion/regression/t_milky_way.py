"""The Milky Way: rotation curve (four Gaia analyses), vertical pull above the Sun, enclosed mass,
escape speed and the inner Galaxy, with published visible matter (McMillan 2017; the local census
of McKee et al. 2015 for the vertical pull). Solver: code/mw_model.py (the field equation on an
(R, z) grid)."""
from __future__ import annotations
import hashlib, json
import numpy as np
import common as C
from checks import make, z_check, range_check

GROUP = 'milky_way'
KEYS = ('bulge', 'thin', 'thick', 'HI', 'H2', 'halo')


def setup(ctx):
    """Components, grid and their Newtonian fields (cached on disk: they do not depend on the law)."""
    def build():
        import mw_model as M, milky_way_v7 as MW
        comps = MW.components(); grid = M.Grid()
        sig = hashlib.sha256(json.dumps({k: {a: (float(v) if np.isscalar(v) and not isinstance(v, str) else str(v))
                                             for a, v in vars(comps[k]).items() if not callable(v)} for k in KEYS},
                                        sort_keys=True).encode() + np.asarray(grid.R).tobytes() + np.asarray(grid.z).tobytes()).hexdigest()[:16]
        path = ctx.cache / f'mw_fields_{sig}.npz'
        if path.exists():
            d = np.load(path); F = {k: (d[k + '_R'], d[k + '_z']) for k in KEYS}
        else:
            ctx.log('Milky Way: computing the Newtonian fields once (cached afterwards)')
            F = {k: M._one_field(comps[k], grid, 60.0, 10, 40.0, 0.005) for k in KEYS}
            np.savez(path, **{k + s: F[k][i] for k in KEYS for i, s in ((0, '_R'), (1, '_z'))})
        jz = grid.nzh
        vb = np.sqrt(np.maximum(-F['bulge'][0][:, jz] * grid.R, 0))
        comps['bulge'].sigma = 0.65 * vb[grid.R < 10].max()           # the SPARC rule for bulges
        return comps, grid, F
    return ctx.get('mw_setup', build)


def run(law, ctx):
    import milky_way_v7 as MW
    comps, grid, F = setup(ctx)
    rc = json.loads((C.RESULTS / 'data/mw_rotation_curves.json').read_text())
    ev = MW.Evaluator(comps, grid, F, law)
    mods = MW.models(comps)
    ctx.log('Milky Way: our law with McMillan (2017) matter')
    r = ev.run(mods['M17'], 'ours', reach=law['reach_kpc'])
    ctx.shared['mw_sun'] = r['sun']
    ctx.log('Milky Way: Newton (for the inner share) and the local census (vertical pull)')
    rn = ev.run(mods['M17'], 'newton', reach=law['reach_kpc'])
    rl = ev.run(mods['local census'], 'ours', reach=law['reach_kpc'])
    R, v = np.array(r['R']), np.array(r['v'])
    out = []
    v0 = float(np.interp(MW.R0, R, v))
    out.append(make('mw.v_sun', GROUP, 'rotation speed at the Sun (8.2 kpc)', v0, crit=range_check(v0, 229.0, 234.0, 7.0), unit='km/s',
                    target='229-234 (Eilers 2019, Zhou 2023, Ou 2024, Jiao 2023); +-7 for the visible-matter model',
                    refs='data/mw_rotation_curves.json'))
    cmp = MW.compare_rc(r, rc)
    outer = float(np.mean([cmp[k]['outer_15_27'] for k in cmp]))
    out.append(make('mw.rc_outer', GROUP, 'rotation curve 15-27 kpc: mean offset over the four Gaia analyses', outer,
                    crit=z_check(outer, 0.0, 5.0), unit='km/s', target='0 +- 5 (the analyses differ by about 10)',
                    detail={k: cmp[k]['outer_15_27'] for k in cmp}))
    out.append(make('mw.rc_rms_eilers', GROUP, 'rotation curve 5-27 kpc: typical miss against Eilers et al. 2019', cmp['Eilers2019']['rms_kms'],
                    unit='km/s', crit=('info', cmp['Eilers2019']['rms_kms'] / 10)))
    s11 = rl['Sigma_z']['1.1']
    out.append(make('mw.vertical_pull', GROUP, 'vertical pull 1.1 kpc above the Sun (as surface density), counted local matter', s11,
                    crit=z_check(s11, 69.8, 3.3), unit='Msun/pc^2', target='69.8 +- 3.3 (Bovy & Rix 2013: 68 +- 4; Holmberg & Flynn 2004: 74 +- 6)',
                    detail=dict(with_McMillan_matter=r['Sigma_z']['1.1'])))
    for key, obs, err, src in (('20.0', 1.91e11, 0.18e11, 'Posti & Helmi 2019'), ('50.0', 4.5e11, 0.4e11, 'Correa Magnus & Vasiliev 2022'),
                               ('100.0', 6.9e11, 0.7e11, 'Deason 2021 + Correa Magnus & Vasiliev 2022'),
                               ('200.0', 1.10e12, 0.25e12, 'Correa Magnus & Vasiliev 2022')):
        m = r['M_eff'][key]
        out.append(make(f'mw.mass_{int(float(key))}', GROUP, f'mass inside {int(float(key))} kpc (as inferred from orbits)', m,
                        crit=z_check(m, obs, err), unit='Msun', target=f'{obs:.2e} +- {err:.1e} ({src})'))
    ve = r['v_esc']['400.0']
    out.append(make('mw.v_escape', GROUP, 'escape speed at the Sun (to 400 kpc)', ve, crit=range_check(ve, 445.0, 580.0, 20.0),
                    unit='km/s', target='445-580 (seven analyses 2014-2024)', detail=dict(to_540kpc=r['v_esc']['540.0'])))
    vN = np.array(rn['v'])
    share = {x: float(np.interp(x, R, vN) / np.interp(x, R, v)) for x in (2.0, 3.0)}
    fv = float(np.mean(list(share.values())))
    out.append(make('mw.inner_share', GROUP, 'inner Galaxy (bar region, 2-3 kpc): share of the rotation speed from visible matter', fv,
                    crit=z_check(fv, 0.88, 0.07), target='0.88 +- 0.07 (Wegg, Gerhard & Portail 2016, bulge microlensing)',
                    detail={f'at_{int(k)}kpc': s for k, s in share.items()}))
    return out
