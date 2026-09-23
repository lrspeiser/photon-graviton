"""Round 8's colliding clusters (full tier): MACS J0025.4-1222, Abell 520 and El Gordo, modelled
with the Bullet's machinery and published inputs only (code/collisions_v8.py)."""
from __future__ import annotations
import numpy as np
from checks import make, z_check, range_check, at_most, _grade

GROUP = 'collisions'


def z_asym(value, obs, err_minus, err_plus):
    """A measurement with unequal error bars."""
    return _grade((value - obs) / (err_plus if value > obs else err_minus))


def macs(law, ctx):
    import collisions_v8 as CV
    spec = CV.macs0025()
    ctx.log('MACS J0025.4-1222: lensing map (t = 0.5 Gyr)')
    sol = CV.solve(spec, law, 0.5)
    m = CV.measure('macs0025', spec, sol)
    ref = 'Bradac et al. 2008 (arXiv:0806.2320)'
    out = [make('macs0025.m300_se', GROUP, 'MACS J0025 SE subcluster: lensing mass inside 300 kpc', m['M300_se'],
                crit=z_asym(m['M300_se'], 2.5e14, 1.7e14, 1.0e14), unit='Msun', target='2.5 (+1.0/-1.7) x 10^14', refs=ref),
           make('macs0025.m300_nw', GROUP, 'MACS J0025 NW subcluster: lensing mass inside 300 kpc', m['M300_nw'],
                crit=z_asym(m['M300_nw'], 2.6e14, 1.4e14, 0.5e14), unit='Msun', target='2.6 (+0.5/-1.4) x 10^14', refs=ref)]
    for w in ('se', 'nw'):
        p = m[f'peak_{w}']
        out.append(make(f'macs0025.peak_{w}', GROUP, f'MACS J0025 {w.upper()}: lensing peak distance from its galaxies', p['to_galaxies'],
                        crit=at_most(p['to_galaxies'], 0.5 * p['galaxies_to_gas'], 0.75 * p['galaxies_to_gas']), unit='kpc',
                        target=f"with the galaxies, not the gas ({p['galaxies_to_gas']:.0f} kpc away): lensing offset from the gas at > 4 sigma", refs=ref))
    out.append(make('macs0025.sigma', GROUP, 'MACS J0025: galaxies\' speed spread (line of sight, inside 1.5 Mpc)', m['sigma_los_1p5Mpc'],
                    crit=z_check(m['sigma_los_1p5Mpc'], 835.0, 59.0), unit='km/s', target='835 +- 59 (108 galaxies within 1.5 Mpc)', refs=ref))
    ctx.shared['macs0025'] = m
    return out


def a520(law, ctx):
    import collisions_v8 as CV
    spec = CV.abell520()
    ctx.log('Abell 520: lensing map (t = 0.3 Gyr)')
    sol = CV.solve(spec, law, 0.3)
    m = CV.measure('abell520', spec, sol)
    jee = dict(P1=(2.10, 0.43), P2=(4.05, 0.28), P3=(3.35, 0.34), P4=(4.23, 0.28), P5=(2.93, 0.39))
    clowe = dict(P1=(2.81, 0.67), P2=(4.16, 0.67), P3=(2.84, 0.64), P4=(5.59, 0.68), P5=(3.17, 0.66), P6=(3.68, 0.68))   # Table 2, zeta_c (unsmoothed)
    out = []
    for k in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6'):
        vals = [d[k] for d in (jee, clowe) if k in d]
        lo, hi = min(v[0] for v in vals) * 1e13, max(v[0] for v in vals) * 1e13
        err = max(v[1] for v in vals) * 1e13
        what = 'the gas-rich, galaxy-poor centre ("dark core")' if k == 'P3' else 'a galaxy clump'
        out.append(make(f'abell520.m150_{k.lower()}', GROUP, f'Abell 520 {k}, {what}: lensing mass inside 150 kpc', m[f'M150_{k}'],
                        crit=range_check(m[f'M150_{k}'], lo, hi, err), unit='Msun',
                        target=f"{lo / 1e13:.2f}-{hi / 1e13:.2f} x 10^13 (Jee et al. 2014; Clowe et al. 2012)",
                        refs='Jee et al. 2014 Table 1; Clowe et al. 2012 Table 2 (unsmoothed)'))
    out.append(make('abell520.m710', GROUP, 'Abell 520: lensing mass inside 710 kpc of the centre', m['M710_P3'],
                    crit=z_check(m['M710_P3'], 5.0e14, 0.55e14), unit='Msun', target='5.0 +- 0.55 x 10^14 (Mahdavi et al. 2007)'))
    # galaxies' speed spread near each lensing peak (Girardi et al. 2008, Table 4; 6-9 galaxies each). The
    # cluster-wide 1,066 km/s mixes in the clumps' motions against each other (about +-1,100 km/s), so it is not used.
    girardi = dict(P1=(811.0, 71.0, 278.0), P2=(749.0, 88.0, 186.0), P4=(579.0, 151.0, 523.0), P5=(668.0, 187.0, 570.0))
    zs, pred = [], {}
    for k, (o, lo, hi) in girardi.items():
        v = m['speeds'][k]['500kpc']; pred[k] = v
        zs.append((v - o) / (hi if v > o else lo))
    rz = float(np.sqrt(np.mean(np.square(zs))))
    out.append(make('abell520.sigma_peaks', GROUP, 'Abell 520: galaxies\' speed spread at P1, P2, P4, P5 (rms z)', rz, crit=_grade(rz),
                    target='Girardi et al. 2008 Table 4: 811, 749, 579, 668 km/s (6-9 galaxies each)', detail=dict(predicted=pred)))
    ctx.shared['abell520'] = m
    return out


def gordo(law, ctx):
    import collisions_v8 as CV
    spec = CV.el_gordo()
    ctx.log('El Gordo: lensing map (t = 0.46 Gyr, the outgoing phase)')
    sol = CV.solve(spec, law, 0.46)
    m = CV.measure('el_gordo', spec, sol)
    obs = CV.observed_el_gordo(spec)
    ref = 'Kim et al. 2021 (arXiv:2106.00031) Table 2 (two NFW haloes, consistent with their model-free masses)'
    out = []
    for R in (500, 1000):
        v, o = m[f'M{R}_com'], obs[f'M{R}_com']
        out.append(make(f'elgordo.m{R}', GROUP, f'El Gordo: lensing mass inside {R} kpc of the centre of mass', v,
                        crit=z_check(v, o, 0.15 * o), unit='Msun', target=f'{o:.2e} +- 15%', refs=ref))
    for w, (o, e) in (('nw', (1290.0, 134.0)), ('se', (1089.0, 200.0))):
        v = m['speeds'][w.upper()]['1000kpc']
        out.append(make(f'elgordo.sigma_{w}', GROUP, f'El Gordo {w.upper()}: galaxies\' speed spread (line of sight, inside 1 Mpc)', v,
                        crit=z_check(v, o, e), unit='km/s', target=f'{o:.0f} +- {e:.0f} (Menanteau et al. 2012)'))
    p = m['peak_se']
    # tracked, not graded: where the bulk of El Gordo's gas sits is our assumption (no published map), and it
    # moves this peak by tens of kpc
    out.append(make('elgordo.peak_se', GROUP, 'El Gordo SE: lensing peak distance from the cool gas core', p['to_cool_core'],
                    crit=('info', None), unit='kpc',
                    target='about 12 arcsec (97 kpc), the lensing peak nearer the merger centre (Kim et al. 2021; Ng et al. 2015)',
                    detail=dict(to_galaxies=p['to_galaxies'], centre_side=bool(p['to_cool_core'] > p['galaxies_to_cool_core']))))
    ctx.shared['el_gordo'] = m
    return out


def run(law, ctx):
    if not ctx.full:
        return []
    return macs(law, ctx) + a520(law, ctx) + gordo(law, ctx)
