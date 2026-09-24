"""Round 11: galaxy-galaxy lensing (KiDS-1000; Brouwer et al. 2021, and Mistele et al. 2024's lensing
circular speeds) in the project's static distance law.

The KiDS lenses lie at 0.1 < z < 0.5, mean z = 0.25 (Brouwer et al. 2021, sec. 2.2), with sources
behind them (KiDS-1000, 0.1 < z_B < 1.2; we take an effective z_s = 0.75), and the paper's distances are
flat LCDM with Omega_m = 0.2793, H0 = 70. At the lenses' redshift, at fixed angle and flux:
    projected radius R            x  D_A(static) / D_A(LCDM)                    (size)
    stellar (and cold-gas) mass   x  (D_L(static) / D_L(LCDM))^2                (stars)
    g_bar = G M / R^2             x  stars / size^2
    Delta Sigma, g_obs = 4 G Delta Sigma   x  Sigma_crit(static) / Sigma_crit(LCDM) = lens / size^2
(collisions_v10's 'lens' is the factor for a lensing MASS inside a fixed angle, D_l D_s / D_ls; a surface
density divides it by the area factor, size^2.)
Each measured point therefore moves to a lower g_bar and a higher g_obs. The law is evaluated at the
moved points, for a typical lens mass converted the same way. Done for z_l = 0.15, 0.25 and 0.35.

    python code/kids_static_v11.py --output-dir run-kids-static-v11
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
import law as L                                   # noqa: E402  the heat weight (round 14: its exponent)

WMAP9 = (70.0, 0.2793)


def tables():
    return dict(all=LC.load('Fig-4-5-C1_RAR-KiDS-isolated_Nobins.txt'), blue=LC.load('Fig-8_RAR-KiDS-isolated_Colorbin_1.txt'),
                red=LC.load('Fig-8_RAR-KiDS-isolated_Colorbin_2.txt'), disc=LC.load('Fig-8_RAR-KiDS-isolated_Sersicbin_1.txt'),
                bulge=LC.load('Fig-8_RAR-KiDS-isolated_Sersicbin_2.txt'), gama=LC.load('Fig-4-C1_RAR-GAMA-isolated_Nobins.txt'))


def kids(consts, u, reach, f=None):
    """Median log10(observed / predicted) per sample (reliable bins), and the early/late gap."""
    tabs = tables()
    fg = 1.0 if f is None else f['stars'] / f['size'] ** 2
    fl = 1.0 if f is None else f['lens'] / f['size'] ** 2          # Sigma_crit ratio (a surface density)
    fm = 1.0 if f is None else f['stars']
    kE = L.k_from_sig2(160.0 ** 2, u)
    gb = tabs['all']['gbar'] * fg
    Mtyp = 10 ** 10.6 * fm
    late = LC.gconv_at(gb, Mtyp, 'ours', consts, k=0.1, reach=reach)
    early = LC.gconv_at(gb, Mtyp, 'ours', consts, k=kE, reach=reach)
    fr = tabs['red']['w'] / (tabs['red']['w'] + tabs['blue']['w'])
    mix = fr * early + (1 - fr) * late
    rel = tabs['all']['gbar'] >= 1e-13
    out = {}
    for s, pred in (('all', mix), ('blue', late), ('red', early), ('disc', late), ('bulge', early), ('gama', mix)):
        obs = dict(tabs[s]); obs['gobs'] = tabs[s]['gobs'] * fl; obs['err'] = tabs[s]['err'] * fl
        out[s] = LC.compare(obs, pred, rel)['median_offset_dex']
    gap_obs = float(np.median(np.log10(tabs['red']['gobs'][rel] / tabs['blue']['gobs'][rel])))
    out['gap_model'] = float(np.median(np.log10(early[rel] / late[rel]))); out['gap_observed'] = gap_obs
    return out


def mistele(consts, u, reach, f=None):
    rows = [l for l in (C.RESULTS / 'data/mistele2024/table1_mrt.txt').read_text().splitlines()
            if l[:3] in ('All', 'LTG', 'ETG') and not l.startswith('All (')]
    tab = {}
    for l in rows:
        ff = [l[i:i + 6] for i in range(33, 117, 7)]
        tab.setdefault(l[:24].strip(), []).append([float(l[25:32])] + [float(x) if x.strip() else np.nan for x in ff])
    fs = 1.0 if f is None else f['size']; fm = 1.0 if f is None else f['stars']
    fl = 1.0 if f is None else f['lens'] / f['size'] ** 2        # Delta Sigma scales as Sigma_crit
    logMb = [10.10, 10.66, 10.96, 11.29]
    out = {}
    for s, k in (('LTG', 0.1), ('ETG', L.k_from_sig2(160 ** 2, u))):
        a = np.array(tab[s]); R = a[:, 0] * fs; m = (a[:, 0] >= 50) & (a[:, 0] <= 300)
        zs, ratios = [], []
        for i, lm in enumerate(logMb):
            vc, es = a[:, 1 + 3 * i] * np.sqrt(fl * fs), a[:, 2 + 3 * i] * np.sqrt(fl * fs)
            ok = m & np.isfinite(vc)
            if ok.sum() < 2: continue
            w = 1 / es[ok] ** 2
            obs = float(np.sum(vc[ok] * w) / w.sum()); err = float(1 / np.sqrt(w.sum()))
            r = np.linspace(50 * fs, 300 * fs, 50)
            pred = float(np.mean(np.sqrt(LC.pull_profile(r, 10 ** lm * fm, 'ours', consts, k=k, reach=reach) / C.K_SI * r)))
            zs.append((pred - obs) / err); ratios.append(obs / pred)
        out[s] = dict(rms_z=float(np.sqrt(np.mean(np.square(zs)))), observed_over_predicted=ratios)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    law = load_law('round9')
    consts9 = dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'])
    stat = json.loads((HERE.parent / 'run-xcop-static-v11/xcop_static_v11.json').read_text())['joint_refit_static'][-1]
    consts11 = dict(a_SI=stat['a_SI'], g_d_SI=stat['g_d_SI'])
    u9, u11 = law['u_kms'], stat['u_kms']
    reach9 = u9 * 1.0227121650537077 * 13.0; reach11 = u11 * 1.0227121650537077 * 13.0
    res = dict(experiment='round 11: KiDS lensing in the static distance law', cosmology_of_the_data=dict(H0=WMAP9[0], Om=WMAP9[1]),
               constants=dict(round9=dict(consts9, u_kms=u9), static_refit=dict(consts11, u_kms=u11)), runs={})
    cases = [('standard distances, round-9 constants', None, consts9, u9, reach9)]
    for zl in (0.15, 0.25, 0.35):
        f = C10.factors(zl, 0.75, WMAP9)
        res.setdefault('factors', {})[f'z_l={zl}'] = dict(size=f['size'], stars=f['stars'], lens_mass=f['lens'], sigma_crit=f['lens'] / f['size'] ** 2,
                                                         g_bar=f['stars'] / f['size'] ** 2)
        cases.append((f'static distances (z_l = {zl}), static-refit constants', f, consts11, u11, reach11))
        if zl == 0.25:
            cases.append((f'static distances (z_l = {zl}), round-9 constants', f, consts9, u9, reach9))
    for tag, f, cs, u, reach in cases:
        k = kids(cs, u, reach, f); m = mistele(cs, u, reach, f)
        res['runs'][tag] = dict(kids=k, mistele=m)
        print(f"{tag}:\n   KiDS median log10(obs/pred): " + ', '.join(f"{s} {k[s]:+.3f}" for s in ('all', 'blue', 'red', 'disc', 'bulge', 'gama')) +
              f"; gap {k['gap_model']:.3f} (observed {k['gap_observed']:.3f})\n   Mistele 50-300 kpc: spirals rms z {m['LTG']['rms_z']:.2f}"
              f" (obs/pred {', '.join(f'{x:.2f}' for x in m['LTG']['observed_over_predicted'])}), ellipticals rms z {m['ETG']['rms_z']:.2f}"
              f" (obs/pred {', '.join(f'{x:.2f}' for x in m['ETG']['observed_over_predicted'])})", flush=True)
    # how the result depends on the distance law itself: its scale alpha (1 + z = e^(alpha D)) and the two
    # static geometries (collisions_v10.VARIANT), at z_l = 0.25 with the static-refit constants held
    alpha0, variant0 = C10.ALPHA, C10.VARIANT
    res['alpha_scan'] = []
    try:
        for variant in ('fixed', 'metric'):
            for s in (0.8, 0.9, 1.0, 1.1):
                C10.VARIANT, C10.ALPHA = variant, alpha0 * s
                f = C10.factors(0.25, 0.75, WMAP9); k = kids(consts11, u11, reach11, f)
                res['alpha_scan'].append(dict(variant=variant, alpha_scale=s, H0_like=s * alpha0 * 299792.458,
                                              **{kk: k[kk] for kk in ('all', 'blue', 'red', 'disc', 'bulge', 'gama', 'gap_model')}))
                print(f"   alpha x{s:.1f} ({variant}): " + ', '.join(f"{kk} {k[kk]:+.3f}" for kk in ('all', 'blue', 'red', 'disc', 'bulge', 'gama')) +
                      f"; gap {k['gap_model']:.3f}", flush=True)
    finally:
        C10.ALPHA, C10.VARIANT = alpha0, variant0
    res['seconds'] = time.monotonic() - t0
    (out / 'kids_static_v11.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(f"wrote {out / 'kids_static_v11.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
