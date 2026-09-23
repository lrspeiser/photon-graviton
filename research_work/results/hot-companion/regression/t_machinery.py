"""Machinery: the standing guard (not MOND, not Newton, not dark matter), the law's limits, and the
constants themselves (tracked, so the scoreboard shows when they change)."""
from __future__ import annotations
import numpy as np
import common as C
from checks import make

GROUP = 'machinery'


def run(law, ctx):
    import law as L, run as R, formula_guard as FG
    a, lam, u = law['a_code'], law['lam'], law['u_kms']
    out = []
    for key, title, unit in (('a_SI', 'deep-regime strength a', 'm/s^2'), ('g_d_SI', 'release scale g_d', 'm/s^2'),
                             ('u_kms', 'companion speed u', 'km/s'), ('reach_kpc', 'reach u x 13 Gyr', 'kpc')):
        out.append(make(f'law.{key}', GROUP, title, law[key], unit=unit, target='constant of the law under test'))
    ell = law['a_SI'] * u * 1e3 / 2
    out.append(make('law.emission', GROUP, 'companion emission per kg of matter, l = a u / 2', ell, unit='W/kg',
                    target='energy budget (round 3: 6.48e-6)'))

    gals = ctx.sparc()
    ctx.log('guard: point mass, temperature, SPARC local-function test')
    rep = FG.guard('law under test', point_mass=L.point_mass(a, u, lam),
                   sparc=[(g['gN'], R.galaxy_g(g, a, u, lam)) for g in gals], a0_si=law['a_SI'])
    bad = [c for c in ('MOND-CLASS', 'NEWTON-CLASS', 'DARK-MATTER-CLASS') if c in rep['verdict']]
    out.append(make('machinery.guard_class', GROUP, 'standing rule: not MOND, not Newton, not dark matter',
                    rep['temperature']['max_spread_dex'], crit=('fail' if bad else 'pass', None), unit='dex',
                    target='no forbidden class in the guard verdict; value = spread of the pull with random speed 0-1000 km/s',
                    note=rep['verdict'], refs='RULES.md; research_work/tools/formula_guard.py'))
    lf = rep['sparc_local_function']['rms_beyond_local_gN_dex']
    out.append(make('machinery.guard_sparc_local', GROUP, 'SPARC pull is not a function of g_N alone', lf,
                    crit=('pass' if lf >= 0.005 else 'fail', None), unit='dex', target='>= 0.005 dex beyond a local function of g_N'))

    # the release factor switches the companion off in strong pulls: extra / Newton at 100 g_d
    gN = 100 * lam * a
    extra = float(L.total(gN, 0.0, a, lam) - gN)
    out.append(make('machinery.strong_field_off', GROUP, 'extra pull at 100 g_d, relative to Newton', extra / gN,
                    crit=('pass' if extra / gN < 1e-12 else 'fail', None), target='< 1e-12'))
    return out
