"""Lensing: KiDS-1000 galaxy-galaxy lensing (Brouwer et al. 2021), lensing circular velocities
(Mistele et al. 2024), the six SLACS strong lenses (light = matter), the Einstein Cross and
microlensing."""
from __future__ import annotations
import numpy as np
import common as C
from checks import make, z_check, at_most, rms_z

GROUP = 'lensing'
KIDS_ERR = 0.025        # dex; a median offset within 0.05 passes (stellar-mass systematics are ~0.1 dex)


def kids(law, ctx):
    import lensing_census_v7 as LC
    u, reach = law['u_kms'], law['reach_kpc']
    kE = lambda s: 3 * s ** 2 / u ** 2
    tabs = dict(all=LC.load('Fig-4-5-C1_RAR-KiDS-isolated_Nobins.txt'), blue=LC.load('Fig-8_RAR-KiDS-isolated_Colorbin_1.txt'),
                red=LC.load('Fig-8_RAR-KiDS-isolated_Colorbin_2.txt'), disc=LC.load('Fig-8_RAR-KiDS-isolated_Sersicbin_1.txt'),
                bulge=LC.load('Fig-8_RAR-KiDS-isolated_Sersicbin_2.txt'), gama=LC.load('Fig-4-C1_RAR-GAMA-isolated_Nobins.txt'))
    gb = tabs['all']['gbar']; Mtyp = 10 ** 10.6
    late = LC.gconv_at(gb, Mtyp, 'ours', law, k=0.1, reach=reach)
    early = LC.gconv_at(gb, Mtyp, 'ours', law, k=kE(160), reach=reach)
    fr = tabs['red']['w'] / (tabs['red']['w'] + tabs['blue']['w'])
    mix = fr * early + (1 - fr) * late
    rel = gb >= 1e-13
    out = []
    for s, pred, what in (('all', mix, 'all isolated lenses (mixture)'), ('blue', late, 'blue lenses (spirals, k = 0.1)'),
                          ('red', early, 'red lenses (ellipticals, stars at 160 km/s)'), ('disc', late, 'Sersic n < 2 (disks)'),
                          ('bulge', early, 'Sersic n > 2 (bulges)'), ('gama', mix, 'GAMA spectroscopic lenses')):
        c = LC.compare(tabs[s], pred, rel)
        out.append(make(f'lensing.kids_{s}', GROUP, f'KiDS lensing pull, {what}: median log10(observed/predicted)', c['median_offset_dex'],
                        crit=z_check(c['median_offset_dex'], 0.0, KIDS_ERR), unit='dex', target='0 +- 0.025 dex (7 reliable bins)',
                        detail=dict(chi2_with_0p1dex=c['chi2_with_0p1dex'], n=c['n']), refs='Brouwer et al. 2021, A&A 650, A113'))
    gap_obs = float(np.median(np.log10(tabs['red']['gobs'][rel] / tabs['blue']['gobs'][rel])))
    gap = float(np.median(np.log10(early[rel] / late[rel])))
    out.append(make('lensing.kids_gap', GROUP, 'KiDS: ellipticals lens more than spirals of the same visible mass (gap)', gap,
                    crit=z_check(gap, gap_obs, 0.04), unit='dex', target=f'{gap_obs:.3f} +- 0.04 dex (same bins)'))
    return out


def mistele(law):
    import lensing_census_v7 as LC
    u, reach = law['u_kms'], law['reach_kpc']
    rows = [l for l in (C.RESULTS / 'data/mistele2024/table1_mrt.txt').read_text().splitlines()
            if l[:3] in ('All', 'LTG', 'ETG') and not l.startswith('All (')]
    tab = {}
    for l in rows:
        f = [l[i:i + 6] for i in range(33, 117, 7)]
        tab.setdefault(l[:24].strip(), []).append([float(l[25:32])] + [float(x) if x.strip() else np.nan for x in f])
    logMb = [10.10, 10.66, 10.96, 11.29]
    out = []
    for s, k, what in (('LTG', 0.1, 'spirals'), ('ETG', 3 * 160 ** 2 / u ** 2, 'ellipticals')):
        a = np.array(tab[s]); R = a[:, 0]; m = (R >= 50) & (R <= 300)
        obs, err, pred = [], [], []
        for i, lm in enumerate(logMb):
            vc, es = a[:, 1 + 3 * i], a[:, 2 + 3 * i]
            ok = m & np.isfinite(vc)
            if ok.sum() < 2: continue
            w = 1 / es[ok] ** 2
            obs.append(float(np.sum(vc[ok] * w) / w.sum())); err.append(float(1 / np.sqrt(w.sum())))
            r = np.linspace(50, 300, 50)
            pred.append(float(np.mean(np.sqrt(LC.pull_profile(r, 10 ** lm, 'ours', law, k=k, reach=reach) / C.K_SI * r))))
        out.append(make(f'lensing.mistele_{s.lower()}', GROUP, f'lensing circular speeds 50-300 kpc, {what} ({len(obs)} mass bins): rms z',
                        None, crit=rms_z(pred, obs, err), target='Mistele et al. 2024 (ApJL 969, L3), Table 1',
                        detail=dict(observed=obs, error=err, predicted=pred)))
        out[-1].value = out[-1].score
    return out


def slacs(law, ctx):
    import lenses_t35 as LT
    LT.A, LT.LAM, LT.U = law['a_code'], law['lam'], law['u_kms']
    orig, reader = LT.patched_readers('standard'); LT.M.J.read_json = reader
    try:
        lenses = [LT.M.make_lens(n) for n in LT.M.LENSES]
    finally:
        LT.M.J.read_json = orig
    rows = [LT.analyse(l, 'hot companion') for l in lenses]
    gaps = np.array([r['slip_gap_dex'] for r in rows]); se = float(gaps.std(ddof=1) / np.sqrt(len(gaps)))
    dm = np.array([r['dm_lens'] for r in rows])
    return [make('lensing.slacs_light_equals_matter', GROUP, 'SLACS: stellar mass from lensing minus from star speeds (mean of 6)', float(gaps.mean()),
                 crit=z_check(float(gaps.mean()), 0.0, se), unit='dex', target=f'0 +- {se:.3f} (the lenses\' own scatter)',
                 detail=dict(gaps=gaps.tolist()), refs='Bolton et al. 2008; Auger et al. 2009; lenses_t35.py'),
            make('lensing.slacs_star_mass', GROUP, 'SLACS: stellar mass needed for the Einstein radii, vs Chabrier', float(dm.mean()),
                 unit='dex', target='Salpeter is +0.25; the IMF trend of Posacki et al. 2015 is derived with dark haloes', crit=('info', None))]


def point_lenses(law):
    GSI, MSUN, KPC = 6.674e-11, 1.989e30, 3.0857e19
    a, gd = law['a_SI'], law['g_d_SI']
    DA = 299792.458 * 0.0394 / 73.0 / (1 + 0.0394) ** 2 * 1e3
    RE = 0.89 / 206265 * DA; gN = GSI * 1.54e10 * MSUN / (RE * KPC) ** 2
    share = float(np.exp(-gN / gd) * np.sqrt(a * gN) / gN)
    Dl, Ds = 4.0, 8.2
    thetaE = np.sqrt(4 * GSI * 0.5 * MSUN / 2.998e8 ** 2 * (Ds - Dl) / (Dl * Ds * KPC)); RE_m = thetaE * Dl * KPC
    gE = GSI * 0.5 * MSUN / RE_m ** 2
    micro = float(np.exp(-gE / gd) * np.sqrt(a * gE) / gE)
    return [make('lensing.einstein_cross', GROUP, 'Einstein Cross: companion share of the pull at the ring', share,
                 crit=at_most(share, 0.15), target='dark matter inside the ring < 15% (Trott et al. 2010; van de Ven et al. 2010)'),
            make('lensing.microlensing', GROUP, 'bulge microlensing: companion share at a star\'s Einstein radius', micro,
                 crit=at_most(micro, 1e-6), target='standard lensing by stars (OGLE, MOA optical depths)')]


def run(law, ctx):
    ctx.log('KiDS lensing')
    out = kids(law, ctx)
    ctx.log('Mistele circular speeds')
    out += mistele(law)
    ctx.log('SLACS lenses')
    out += slacs(law, ctx)
    out += point_lenses(law)
    return out
