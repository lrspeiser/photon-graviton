"""Round 12: galaxy lensing (KiDS-1000, Mistele et al. 2024) with the lenses' heat measured, not assumed.

Same comparison as kids_static_v11 (the project's static distances, lenses at z = 0.25, sources at an
effective 0.75), with one change: each sample's heat weight comes from the stars' speeds of galaxies of
the same stellar mass and the same split in SDSS DR17 (code/lens_heat_sdss_v12.py ->
data/lens_heat_sdss_v12.json): the de Vaucouleurs part of each galaxy carries k = 3 sigma_e^2 / u^2,
the disk part k = 3 (30 km/s)^2 / u^2. KiDS samples take the mean over log M* = 10.3-10.9 (around the
stacks' typical 10.6); Mistele's four mass bins take their own bin. Round 11 assumed 160 km/s for every
star of a red or bulge-dominated lens, and k = 0.1 for a blue or disk-dominated one.

    python code/kids_heat_v12.py --output-dir run-kids-heat-v12
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import lensing_census_v7 as LC                    # noqa: E402
import collisions_v10 as C10                      # noqa: E402
import common as C                                # noqa: E402
import kids_static_v11 as KS                      # noqa: E402
import lens_heat_sdss_v12 as LH                   # noqa: E402

HEAT = json.loads((HERE.parent / 'data/lens_heat_sdss_v12.json').read_text())
PAIRS = dict(red='red', blue='blue', bulge='bulge', disc='disc')


def heat(u):
    return {s: LH.k_eff(HEAT, s, u) for s in PAIRS}


def kids(consts, u, reach, f=None, k=None):
    """Median log10(observed / predicted) per sample (reliable bins), and the early/late gaps."""
    tabs = KS.tables()
    k = heat(u) if k is None else k
    fg = 1.0 if f is None else f['stars'] / f['size'] ** 2
    fl = 1.0 if f is None else f['lens'] / f['size'] ** 2
    fm = 1.0 if f is None else f['stars']
    gb = tabs['all']['gbar'] * fg
    Mtyp = 10 ** 10.6 * fm
    g = {s: LC.gconv_at(gb, Mtyp, 'ours', consts, k=k[s], reach=reach) for s in PAIRS}
    fr = tabs['red']['w'] / (tabs['red']['w'] + tabs['blue']['w'])
    mix = fr * g['red'] + (1 - fr) * g['blue']
    rel = tabs['all']['gbar'] >= 1e-13
    out = dict(k=k)
    for s, pred in (('all', mix), ('blue', g['blue']), ('red', g['red']), ('disc', g['disc']), ('bulge', g['bulge']), ('gama', mix)):
        obs = dict(tabs[s]); obs['gobs'] = tabs[s]['gobs'] * fl; obs['err'] = tabs[s]['err'] * fl
        out[s] = LC.compare(obs, pred, rel)['median_offset_dex']
    out['gap_model'] = float(np.median(np.log10(g['red'][rel] / g['blue'][rel])))
    out['gap_observed'] = float(np.median(np.log10(tabs['red']['gobs'][rel] / tabs['blue']['gobs'][rel])))
    out['gap_sersic_model'] = float(np.median(np.log10(g['bulge'][rel] / g['disc'][rel])))
    out['gap_sersic_observed'] = float(np.median(np.log10(tabs['bulge']['gobs'][rel] / tabs['disc']['gobs'][rel])))
    return out


def mistele(consts, u, reach, f=None):
    """Mistele et al.'s lensing circular speeds, 50-300 kpc: late types take the blue heat, early types the red,
    each at its own stellar-mass bin."""
    rows = [l for l in (C.RESULTS / 'data/mistele2024/table1_mrt.txt').read_text().splitlines()
            if l[:3] in ('All', 'LTG', 'ETG') and not l.startswith('All (')]
    tab = {}
    for l in rows:
        ff = [l[i:i + 6] for i in range(33, 117, 7)]
        tab.setdefault(l[:24].strip(), []).append([float(l[25:32])] + [float(x) if x.strip() else np.nan for x in ff])
    fs = 1.0 if f is None else f['size']; fm = 1.0 if f is None else f['stars']
    fl = 1.0 if f is None else f['lens'] / f['size'] ** 2
    logMb = [10.10, 10.66, 10.96, 11.29]
    out = {}
    for s, sample in (('LTG', 'blue'), ('ETG', 'red')):
        a = np.array(tab[s]); m = (a[:, 0] >= 50) & (a[:, 0] <= 300)
        zs, ratios, ks = [], [], []
        for i, lm in enumerate(logMb):
            vc, es = a[:, 1 + 3 * i] * np.sqrt(fl * fs), a[:, 2 + 3 * i] * np.sqrt(fl * fs)
            ok = m & np.isfinite(vc)
            if ok.sum() < 2: continue
            k = LH.k_eff(HEAT, sample, u, logM=lm)
            w = 1 / es[ok] ** 2
            obs = float(np.sum(vc[ok] * w) / w.sum()); err = float(1 / np.sqrt(w.sum()))
            r = np.linspace(50 * fs, 300 * fs, 50)
            pred = float(np.mean(np.sqrt(LC.pull_profile(r, 10 ** lm * fm, 'ours', consts, k=k, reach=reach) / C.K_SI * r)))
            zs.append((pred - obs) / err); ratios.append(obs / pred); ks.append(k)
        out[s] = dict(rms_z=float(np.sqrt(np.mean(np.square(zs)))), observed_over_predicted=ratios, k=ks)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    law = load_law('round11')
    consts = dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI']); u, reach = law['u_kms'], law['reach_kpc']
    f = C10.factors(0.25, 0.75, KS.WMAP9)
    old = dict(red=3 * 160 ** 2 / u ** 2, blue=0.1, bulge=3 * 160 ** 2 / u ** 2, disc=0.1)
    res = dict(experiment='round 12: KiDS and Mistele with the lenses\' heat measured from SDSS (static distances, round-11 constants)',
               heat_measured=heat(u), heat_round11=old, runs={})
    for tag, k in (('round-11 heat (160 km/s; k 0.1)', old), ('measured heat (SDSS)', None)):
        kk = kids(consts, u, reach, f, k)
        mi = mistele(consts, u, reach, f) if k is None else KS.mistele(consts, u, reach, f)
        res['runs'][tag] = dict(kids=kk, mistele=mi)
        print(f"{tag}: KiDS " + ', '.join(f"{s} {kk[s]:+.3f}" for s in ('all', 'blue', 'red', 'disc', 'bulge', 'gama')) +
              f"; gaps colour {kk['gap_model']:.3f} (obs {kk['gap_observed']:.3f}), Sersic {kk['gap_sersic_model']:.3f} (obs {kk['gap_sersic_observed']:.3f});"
              f" Mistele rms z spirals {mi['LTG']['rms_z']:.2f}, ellipticals {mi['ETG']['rms_z']:.2f}", flush=True)
    res['seconds'] = time.monotonic() - t0
    (out / 'kids_heat_v12.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(f"wrote {out / 'kids_heat_v12.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
