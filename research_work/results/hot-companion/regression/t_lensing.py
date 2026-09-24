"""Lensing: KiDS-1000 galaxy-galaxy lensing (Brouwer et al. 2021), lensing circular velocities
(Mistele et al. 2024), the six SLACS strong lenses (light = matter), the Einstein Cross and
microlensing."""
from __future__ import annotations
import numpy as np
import common as C
from checks import make, z_check, at_most, rms_z, _grade

GROUP = 'lensing'
KIDS_ERR = 0.025        # dex; a median offset within 0.05 passes (stellar-mass systematics are ~0.1 dex)


def kids(law, ctx):
    """KiDS-1000 lensing RAR of isolated lenses, converted to the static distance law at the lenses' mean
    redshift (z = 0.25; sources at an effective 0.75): g_bar x stars/size^2, g_obs x Sigma_crit ratio
    (code/kids_static_v11.py). Round 10 and earlier graded the paper's flat-LCDM units. Since round 12 each
    sample's heat weight is measured from SDSS galaxies of the same stellar mass and split
    (code/lens_heat_sdss_v12.py, code/kids_heat_v12.py); round 11 assumed 160 km/s for red or bulge-dominated
    lenses and k = 0.1 for blue or disk-dominated ones."""
    import kids_static_v11 as KS
    import kids_heat_v12 as KH
    import collisions_v10 as C10
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    k = KH.kids(dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI']), law['u_kms'], law['reach_kpc'], f)
    kk = k['k']
    out = []
    for s, what in (('all', 'all isolated lenses (mixture)'), ('blue', f"blue lenses (k = {kk['blue']:.2f}, measured)"),
                    ('red', f"red lenses (k = {kk['red']:.2f}, measured)"), ('disc', f"Sersic n < 2 (disks; k = {kk['disc']:.2f})"),
                    ('bulge', f"Sersic n > 2 (bulges; k = {kk['bulge']:.2f})"), ('gama', 'GAMA spectroscopic lenses')):
        out.append(make(f'lensing.kids_{s}', GROUP, f'KiDS lensing pull, {what}: median log10(observed/predicted), static distances', k[s],
                        crit=z_check(k[s], 0.0, KIDS_ERR), unit='dex', target='0 +- 0.025 dex (7 reliable bins)',
                        refs='Brouwer et al. 2021, A&A 650, A113; converted with code/kids_static_v11.py; heat from code/lens_heat_sdss_v12.py'))
    out.append(make('lensing.kids_gap', GROUP, 'KiDS: red lenses lens more than blue ones of the same visible mass (colour gap)', k['gap_model'],
                    crit=z_check(k['gap_model'], k['gap_observed'], 0.04), unit='dex', target=f"{k['gap_observed']:.3f} +- 0.04 dex (same bins)"))
    out.append(make('lensing.kids_gap_sersic', GROUP, 'KiDS: bulge-dominated lenses lens more than disk-dominated ones (Sersic gap)', k['gap_sersic_model'],
                    crit=z_check(k['gap_sersic_model'], k['gap_sersic_observed'], 0.04), unit='dex',
                    target=f"{k['gap_sersic_observed']:.3f} +- 0.04 dex (same bins)"))
    return out


def mistele(law):
    """Lensing circular speeds 50-300 kpc (Mistele et al. 2024; KiDS-1000 lenses), converted to the static
    distance law as for KiDS: radii x size, v_c x sqrt(lens/size), stellar masses x stars. Since round 12 the
    heat of each mass bin is measured (late types: blue, early types: red; code/kids_heat_v12.py)."""
    import kids_static_v11 as KS
    import kids_heat_v12 as KH
    import collisions_v10 as C10
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    m = KH.mistele(dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI']), law['u_kms'], law['reach_kpc'], f)
    out = []
    for s, what in (('LTG', 'spirals'), ('ETG', 'ellipticals')):
        v = m[s]['rms_z']
        out.append(make(f'lensing.mistele_{s.lower()}', GROUP, f'lensing circular speeds 50-300 kpc, {what}: rms z (static distances)',
                        v, crit=_grade(v), target='Mistele et al. 2024 (ApJL 969, L3), Table 1',
                        detail=dict(observed_over_predicted=m[s]['observed_over_predicted'])))
    return out


def slacs(law, ctx):
    import lenses_t35 as LT
    LT.A, LT.LAM, LT.U = law['a_code'], law['lam'], law['u_kms']
    # the project's static distances and its stellar masses converted with them (round 10; rounds 3-9 graded
    # the flat-LCDM convention, whose gap was -0.017 +- 0.024 dex and whose stars needed 1.05-1.35 x Salpeter)
    orig, reader = LT.patched_readers('project'); LT.M.J.read_json = reader
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
            make('lensing.slacs_star_mass', GROUP, 'SLACS: stellar mass needed for the Einstein radii, vs Chabrier (static distances)', float(dm.mean()),
                 unit='dex', target='Salpeter is +0.25; the IMF trend of Posacki et al. 2015 is derived with dark haloes', crit=('info', None),
                 detail=dict(per_lens=dm.tolist(), salpeter_factor=[float(10 ** (x - 0.25)) for x in dm]))]


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
