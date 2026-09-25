"""Round 19: a frozen prediction. Galaxy lensing at fixed visible mass grows with the stars' own speeds.

The law's heat term: stars moving randomly with dispersion sigma weigh (1 + k), k = 3 sigma^2 / u^2 (u = the adopted
169.4 km/s), in the companion's intensity, so far from a galaxy the extra pull is sqrt(a G M (1 + k))/r. At fixed
stellar mass and environment, a galaxy whose stars move faster should lens more, by sqrt((1 + k1)/(1 + k2)) in the
deep regime. MOND predicts no dependence at fixed visible mass; in dark-matter models any dependence runs through the
halo mass that goes with a given sigma at fixed stellar mass, which is a separate, measurable relation.

Frozen here, before any split of the data by sigma has been looked at by this project: for isolated, bulge-dominated
lenses (pure spheroids, heat of the whole stellar mass at the central dispersion sigma_e) of log M* = 10.6 at z = 0.25
(the KiDS-1000 stacks' typical lens, Brouwer et al. 2021), the lensing acceleration g_obs at each of the survey's
g_bar bins, relative to sigma_e = 200 km/s, computed with the adopted law (round-12 constants, release factor, reach)
and the suite's own forward model (code/lensing_census_v7.gconv_at: 4 G Delta Sigma at the radius where G M/R^2 = g_bar).
The numbers are written to run-frozen-prediction-v19/ with a SHA-256 of their content.

    python code/frozen_prediction_v19.py --output-dir run-frozen-prediction-v19
"""
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))

SIGMAS = (100.0, 150.0, 200.0, 250.0, 300.0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    import common as C, law as L, lensing_census_v7 as LC, kids_static_v11 as KS, collisions_v10 as C10
    from law_config import load_law
    law = load_law('round12'); ctx = C.Context(tier='quick', verbose=False); C.apply_distances(law, ctx)
    u = law['u_kms']
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    tabs = KS.tables()
    gb_obs = tabs['all']['gbar']                               # the survey's g_bar bins (its own units)
    gb = gb_obs * f['stars'] / f['size'] ** 2                  # in the project's distances (as the suite does)
    Mtyp = 10 ** 10.6 * f['stars']
    consts = dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'])
    rel = gb_obs >= 1e-13
    rows = {}
    for sig in SIGMAS:
        k = float(L.heat_weight(sig, u))
        g = LC.gconv_at(gb, Mtyp, 'ours', consts, k=k, reach=law['reach_kpc'])
        rows[sig] = dict(k=k, g_obs_model=(g / (f['lens'] / f['size'] ** 2)).tolist())
    ref = np.array(rows[200.0]['g_obs_model'])
    table = []
    for i, gbv in enumerate(gb_obs):
        table.append(dict(g_bar=float(gbv), reliable=bool(rel[i]), **{f'sigma {int(s)}': float(np.log10(rows[s]['g_obs_model'][i] / ref[i])) for s in SIGMAS}))
    deep = [t for t in table if t['reliable'] and t['g_bar'] < 1e-12]
    summary = {f'sigma {int(s)} vs 200 (median over reliable bins below 1e-12 m/s^2), dex':
               float(np.median([t[f'sigma {int(s)}'] for t in deep])) for s in SIGMAS}
    summary['sigma 250 vs 150, dex'] = float(np.median([t['sigma 250'] - t['sigma 150'] for t in deep]))
    summary['deep-regime expectation sqrt((1+k250)/(1+k150)), dex'] = float(0.5 * np.log10((1 + rows[250.0]['k']) / (1 + rows[150.0]['k'])))
    pred = dict(statement='At fixed stellar mass (log M* = 10.6) and isolation, bulge-dominated lenses whose stars move faster lens more: '
                          'g_obs rises with the central stellar velocity dispersion sigma_e as tabulated (relative to sigma_e = 200 km/s), '
                          'about +0.17 dex from 150 to 250 km/s in the deep regime; MOND predicts no dependence at fixed visible mass.',
                law=dict(name=law['name'], a_SI=law['a_SI'], g_d_SI=law['g_d_SI'], u_kms=u, reach_kpc=law['reach_kpc']),
                lens=dict(log10_Mstar=10.6, z=0.25, sources_z=0.75, isolated=True, heat='whole stellar mass at sigma_e (pure spheroid)'),
                k_by_sigma={int(s): rows[s]['k'] for s in SIGMAS}, table=table, summary=summary,
                how_to_test='split isolated KiDS/GAMA (or SDSS/HSC) lenses of one stellar-mass bin by measured sigma_e (e.g. SDSS/GAMA spectra), '
                            'stack each group\'s excess surface density, convert to g_obs at the same g_bar bins, and compare the ratios with this table; '
                            'the environment (isolation) and the mass bin must be the same for both groups.',
                frozen='2026-09-25, round 19, before this project looked at any sigma-split lensing data')
    txt = json.dumps(pred, indent=1, sort_keys=True)
    pred_hash = hashlib.sha256(txt.encode()).hexdigest()
    (out / 'frozen_prediction_v19.json').write_text(txt + '\n')
    (out / 'SHA256').write_text(f'{pred_hash}  frozen_prediction_v19.json\n')
    print(json.dumps(summary, indent=1)); print('sha256', pred_hash, f'({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
