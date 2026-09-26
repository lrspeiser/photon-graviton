"""The Milky Way's ten dwarf spheroidals: their stars' speed spread, predicted from their light
(M/L_V = 2) with the Galaxy's pull and heat acting on them (scaled by the law's external hold)."""
from __future__ import annotations
import json
import numpy as np
import common as C
from checks import make, z_check
import t_milky_way as TMW

GROUP = 'dwarfs'


def galaxy_profile(law, ctx):
    """Spherical shells of the Galaxy's visible matter and its heat-weighted part (McMillan 2017)."""
    import mw_model as M, milky_way_v7 as MW
    comps, grid, F = TMW.setup(ctx)
    s = np.geomspace(1e-3, 3000.0, 1500)
    dm = np.zeros(s.size); dmk = np.zeros(s.size)
    for k, f in MW.models(comps)['M17'].items():
        c = comps[k]
        if isinstance(c, M.Disk):
            Rg = np.linspace(0, 200, 40001); dmr = 2 * np.pi * Rg * c.Sigma(Rg)
            Mc = np.concatenate([[0], np.cumsum(0.5 * np.diff(Rg) * (dmr[1:] + dmr[:-1]))])
            Ms = f * np.interp(s, Rg, Mc)
        else:
            Ms = f * c.m_profile(s)
        d = np.diff(np.concatenate([[0], Ms])); dm += d
        if k in ('bulge', 'halo'):
            if law.get('companion_trapping', 'none') != 'none':      # round 25: only the glow born outside the Galaxy's zone
                import trapping_v25 as TR
                d = d * TR.mw_shell_escape(law, ctx, c.name, s)
            dmk += 3 * c.sigma ** 2 / law['u_kms'] ** 2 * d
    return dict(s=s, dm=dm, dmk=dmk)


def run(law, ctx):
    import mw_dwarfs_v7 as D
    data = json.loads((C.RESULTS / 'data/mw_dwarfs.json').read_text())
    mw = galaxy_profile(law, ctx)
    hold = law['external_hold']
    out, chi = [], 0.0
    ctx.log(f'dwarfs (external hold x{hold:g})')
    for d in data['dwarfs']:
        env = D.galaxy_env(d['D_gc_kpc'], mw, law)
        env = dict(gN=hold * env['gN'], S=hold * env['S'])
        sig = D.sigma_los(2.0 * d['L_V'], d['r_h_pc'] / 1000.0, env, 'ours', law, 3 * d['sigma_obs'] ** 2 / law['u_kms'] ** 2)[0]
        crit = z_check(sig, d['sigma_obs'], d['sigma_err'])
        chi += ((sig - d['sigma_obs']) / d['sigma_err']) ** 2
        out.append(make(f"dwarfs.{d['name'].lower().replace(' ', '_')}", GROUP, f"{d['name']}: stars' speed spread", sig, crit=crit, unit='km/s',
                        target=f"{d['sigma_obs']} +- {d['sigma_err']}", refs=d.get('source', '')))
    out.append(make('dwarfs.chi2', GROUP, 'all ten dwarfs: chi-squared (MOND with the same stars: 119)', chi, crit=('info', chi / 10)))
    return out
