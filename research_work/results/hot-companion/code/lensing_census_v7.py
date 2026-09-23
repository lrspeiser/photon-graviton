#!/usr/bin/env python3
"""Round 7: the standard lensing tests, against our law, MOND and Newton.

    python lensing_census_v7.py --output-dir ../run-lensing-census-v7

1. The lensing radial acceleration relation of isolated galaxies (Brouwer et al. 2021, KiDS-1000;
   data release in ../data/brouwer2021_kids/). Galaxy-galaxy lensing measures the average excess
   surface density DeltaSigma(R) around ~259,000 isolated lens galaxies out to 3 Mpc; Brouwer et
   al. turn it into an 'observed' pull with the isothermal-sphere conversion g_obs = 4 G DeltaSigma
   and bin it by the Newtonian pull of the stars and cold gas, g_bar = G M_gal / R^2. We put each
   law through exactly the same steps: the law's pull around a lens of M_gal (spherical, so the
   field equation gives the pull directly) -> 'as if' mass profile -> projected DeltaSigma(R) ->
   4 G DeltaSigma at the R where g_bar equals each bin's value. Early-type (red, or Sersic n > 2)
   lenses are hot: k = 3 sigma^2/u^2 with sigma = 160 km/s (130-190); late-type (blue, n < 2)
   lenses are cold disks with small bulges: k = 0.1 (0-0.3). Beyond the companion's reach,
   u t = 2.0-2.6 Mpc, our law's extra pull stops.
2. Circular velocities to 2.5 Mpc from the same lensing (Mistele et al. 2024, ApJL 969, L3).
3. The Einstein Cross (Q2237+0305): a spiral's bulge lensing a quasar 0.89 arcsec out.
4. Microlensing towards the bulge: our law's extra pull is zero around any star (Einstein radii of
   a few AU), so each event is standard; the optical depth counts the stars.
5. Colliding clusters beyond the Bullet: what the literature measures (Abell 520, MACS J0025,
   El Gordo) and Abell 1689 -- recorded for the next round's modelling.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
G = 4.30091727003628e-6; K_SI = 1e6 / 3.0856775814913673e19
KPC = 3.0857e19; MSUN = 1.989e30; GSI = 6.674e-11
DATA = HERE.parent / 'data/brouwer2021_kids'
FACTOR = 5.579488e-13                    # m/s^2 per Msun/pc^2 for g = 4 G DeltaSigma (Brouwer et al. 2021, eq. 7; README)
A0 = 1.2e-10


def load(name):
    a = np.loadtxt(DATA / name)
    gbar = a[:, 0]; esd = a[:, 1] / a[:, 4]; err = a[:, 3] / a[:, 4]; wk2 = a[:, 6]
    return dict(gbar=gbar, gobs=FACTOR * esd, err=FACTOR * err, w=wk2)


def pull_profile(r, M, law, consts, k=0.0, fstar=0.85, reach=None, gas=None):
    """Spherical pull (m/s^2) at radii r (kpc) around a lens of baryonic point mass M (Msun)."""
    Mb = M + (gas(r) if gas else 0.0)
    gN = G * Mb / r ** 2 * K_SI
    if law == 'newton':
        return gN
    if law == 'mond':                   # McGaugh et al. 2016 fitting function, g_dagger = 1.2e-10
        return gN / (1 - np.exp(-np.sqrt(gN / A0)))
    a, gd = consts['a_SI'], consts['g_d_SI']
    S = k * G * fstar * M / r ** 2 * K_SI
    extra = np.exp(-gN / gd) * np.sqrt(a * (gN + S))
    if reach is not None:
        extra = np.where(r < reach, extra, 0.0)
    return gN + extra


def gconv_at(gbar_vals, M, law, consts, **kw):
    """Forward model: 4 G DeltaSigma at the radius where G M / R^2 = g_bar (the paper's conversion)."""
    r = np.geomspace(1e-3, 3e4, 20000)                                  # kpc
    g = pull_profile(r, M, law, consts, **kw)
    Meff = g / K_SI * r ** 2 / G                                        # 'as if' mass inside r (Msun)
    dM = np.diff(np.concatenate([[0.0], Meff]))                         # shells (the first holds the point mass)
    Rp = np.sqrt(G * M * K_SI / gbar_vals)                              # kpc
    Rg = np.geomspace(Rp.min() / 1.5, Rp.max() * 1.5, 400)
    x = np.clip(Rg[:, None] / r[None, :], 0, 1)
    M2 = (dM[None, :] * np.where(r[None, :] <= Rg[:, None], 1.0, 1 - np.sqrt(1 - x ** 2))).sum(1)
    dM2 = np.gradient(M2, Rg)
    dsig = M2 / (np.pi * Rg ** 2) - dM2 / (2 * np.pi * Rg)             # Msun/kpc^2
    dsig_pc = np.interp(np.log(Rp), np.log(Rg), dsig) / 1e6
    return FACTOR * dsig_pc


def compare(obs, model, sel):
    d = np.log10(obs['gobs'][sel]) - np.log10(model[sel])
    e = obs['err'][sel] / (obs['gobs'][sel] * np.log(10))
    return dict(median_offset_dex=float(np.median(d)), mean_offset_dex=float(np.average(d, weights=1 / e ** 2)),
                chi2_stat=float(np.sum((d / e) ** 2)), chi2_with_0p1dex=float(np.sum((d ** 2) / (e ** 2 + 0.1 ** 2))), n=int(sel.sum()))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    consts = json.loads((HERE.parent / 'run-v3/results.json').read_text())['constants']
    u = consts['u_kms']; reach = 197.41 * 1.0227 * 13.0
    kE = lambda s: 3 * s ** 2 / u ** 2
    res = {}

    # ---- 1. KiDS lensing RAR
    tabs = dict(all=load('Fig-4-5-C1_RAR-KiDS-isolated_Nobins.txt'), blue=load('Fig-8_RAR-KiDS-isolated_Colorbin_1.txt'),
                red=load('Fig-8_RAR-KiDS-isolated_Colorbin_2.txt'), disc=load('Fig-8_RAR-KiDS-isolated_Sersicbin_1.txt'),
                bulge=load('Fig-8_RAR-KiDS-isolated_Sersicbin_2.txt'), gama=load('Fig-4-C1_RAR-GAMA-isolated_Nobins.txt'))
    gb = tabs['all']['gbar']
    Mtyp = 10 ** 10.6
    pred = {}
    for label, law, kw in (('newton', 'newton', {}), ('mond (M16 fit)', 'mond', {}),
                           ('ours, cold (k=0)', 'ours', dict(k=0.0, reach=reach)),
                           ('ours, late type (k=0.1)', 'ours', dict(k=0.1, reach=reach)),
                           ('ours, late type (k=0.3)', 'ours', dict(k=0.3, reach=reach)),
                           ('ours, early type (sigma 130)', 'ours', dict(k=kE(130), reach=reach)),
                           ('ours, early type (sigma 160)', 'ours', dict(k=kE(160), reach=reach)),
                           ('ours, early type (sigma 190)', 'ours', dict(k=kE(190), reach=reach)),
                           ('ours, late type, hot gas = M* to 100 kpc', 'ours', dict(k=0.1, reach=reach, gas=lambda r: 0.85 * Mtyp * np.minimum(r / 100.0, 1.0))),
                           ('ours, early type sigma 160, hot gas = M* to 100 kpc', 'ours', dict(k=kE(160), reach=reach, gas=lambda r: 0.85 * Mtyp * np.minimum(r / 100.0, 1.0)))):
        pred[label] = gconv_at(gb, Mtyp, law, consts, **kw)
    # mass dependence (only through the reach and release): check at 10^10.1 and 10^11
    spread = {m: gconv_at(gb, 10 ** m, 'ours', consts, k=0.1, reach=reach) for m in (10.1, 11.0)}
    # the mixture for 'all': use the red and blue lensing weights per bin
    wr, wb = tabs['red']['w'], tabs['blue']['w']
    fr = wr / (wr + wb)
    pred['ours, mixture (red share from the lensing weights; sigma 160, k_late 0.1)'] = fr * pred['ours, early type (sigma 160)'] + (1 - fr) * pred['ours, late type (k=0.1)']
    reliable = gb >= 1e-13
    comp = {}
    for sample, keys in (('all', ['newton', 'mond (M16 fit)', 'ours, mixture (red share from the lensing weights; sigma 160, k_late 0.1)']),
                         ('blue', ['newton', 'mond (M16 fit)', 'ours, late type (k=0.1)', 'ours, late type (k=0.3)', 'ours, late type, hot gas = M* to 100 kpc']),
                         ('disc', ['mond (M16 fit)', 'ours, late type (k=0.1)', 'ours, late type (k=0.3)']),
                         ('red', ['newton', 'mond (M16 fit)', 'ours, early type (sigma 130)', 'ours, early type (sigma 160)', 'ours, early type (sigma 190)', 'ours, early type sigma 160, hot gas = M* to 100 kpc']),
                         ('bulge', ['mond (M16 fit)', 'ours, early type (sigma 130)', 'ours, early type (sigma 160)']),
                         ('gama', ['mond (M16 fit)', 'ours, mixture (red share from the lensing weights; sigma 160, k_late 0.1)'])):
        comp[sample] = {}
        for k in keys:
            comp[sample][k] = dict(reliable=compare(tabs[sample], pred[k], reliable), all_bins=compare(tabs[sample], pred[k], np.ones(gb.size, bool)))
    # the early/late gap: data vs law
    gap_obs = float(np.median(np.log10(tabs['red']['gobs'][reliable] / tabs['blue']['gobs'][reliable])))
    gap_ours = float(np.median(np.log10(pred['ours, early type (sigma 160)'][reliable] / pred['ours, late type (k=0.1)'][reliable])))
    # what amplitude would the late types need? (deep-regime companion scaled by c)
    need = {}
    for c in (1.0, 1.1, 1.2, 1.3, 1.4, 1.5):
        p = gconv_at(gb, Mtyp, 'ours', dict(consts, a_SI=consts['a_SI'] * c ** 2), k=0.1, reach=reach)
        need[str(c)] = compare(tabs['blue'], p, reliable)['median_offset_dex']
    res['kids'] = dict(gbar=gb.tolist(), observed={k: dict(gobs=v['gobs'].tolist(), err=v['err'].tolist()) for k, v in tabs.items()},
                       predictions={k: v.tolist() for k, v in pred.items()}, red_share=fr.tolist(), comparison=comp,
                       early_late_gap=dict(observed_median_dex_reliable=gap_obs, ours_median_dex=gap_ours, paper_values=dict(colour=0.27, sersic=0.17, colour_reliable=0.19, sersic_reliable=0.14)),
                       mass_dependence=dict(ratio_1e10p1=np.log10(spread[10.1] / pred['ours, late type (k=0.1)']).tolist(), ratio_1e11=np.log10(spread[11.0] / pred['ours, late type (k=0.1)']).tolist()),
                       late_type_companion_scale_needed=need)
    print('KiDS-1000 lensing RAR (isolated lenses), bins with g_bar >= 1e-13 (isolation reliable): median log10(observed/predicted)', flush=True)
    for s, d in comp.items():
        print(f'  {s:5s}: ' + '; '.join(f"{k}: {v['reliable']['median_offset_dex']:+.3f} (chi2 {v['reliable']['chi2_with_0p1dex']:.1f}/{v['reliable']['n']})" for k, v in d.items()), flush=True)
    print(f'  early/late gap: observed {gap_obs:.3f} dex (reliable bins), ours {gap_ours:.3f}', flush=True)
    print('  late types: companion amplitude scale needed (median offset): ' + ', '.join(f'x{c}: {v:+.3f}' for c, v in need.items()), flush=True)

    # ---- 2. Mistele et al. 2024: circular velocities from lensing, 50-300 kpc
    rows = [l for l in (HERE.parent / 'data/mistele2024/table1_mrt.txt').read_text().splitlines() if l[:3] in ('All', 'LTG', 'ETG') and not l.startswith('All (')]
    def parse(l):
        f = [l[33:39], l[40:46], l[47:53], l[54:60], l[61:67], l[68:74], l[75:81], l[82:88], l[89:95], l[96:102], l[103:109], l[110:116]]
        vals = [float(x) if x.strip() else np.nan for x in f]
        return l[:24].strip(), float(l[25:32]), vals
    tab = {}
    for l in rows:
        s, R, v = parse(l)
        tab.setdefault(s, []).append([R] + v)
    logMb = [10.10, 10.66, 10.96, 11.29]
    mist = {}
    for s in ('All', 'LTG', 'ETG'):
        a = np.array(tab[s]); R = a[:, 0]
        m = (R >= 50) & (R <= 300)
        for i, lm in enumerate(logMb):
            vc, es = a[:, 1 + 3 * i], a[:, 2 + 3 * i]
            ok = m & np.isfinite(vc)
            if ok.sum() < 2: continue
            w = 1 / es[ok] ** 2
            vf = float(np.sum(vc[ok] * w) / w.sum()); ef = float(1 / np.sqrt(w.sum()))
            Mb = 10 ** lm; r = np.linspace(50, 300, 50)
            v_ours_late = float(np.mean(np.sqrt(pull_profile(r, Mb, 'ours', consts, k=0.1, reach=reach) / K_SI * r)))
            v_ours_early = float(np.mean(np.sqrt(pull_profile(r, Mb, 'ours', consts, k=kE(160), reach=reach) / K_SI * r)))
            v_mond = float(np.mean(np.sqrt(pull_profile(r, Mb, 'mond', consts) / K_SI * r)))
            mist[f'{s} logMb {lm}'] = dict(V_flat=vf, err=ef, ours_late=v_ours_late, ours_early=v_ours_early, mond=v_mond)
            print(f'  Mistele 2024 {s:3s} log Mb {lm}: V(50-300 kpc) = {vf:.0f} +- {ef:.0f} km/s; ours late {v_ours_late:.0f}, ours early (sigma 160) {v_ours_early:.0f}; MOND {v_mond:.0f}', flush=True)
    res['mistele2024'] = mist

    # ---- 3. Einstein Cross
    lit = json.loads((HERE.parent / 'data/lensing_literature_v7.json').read_text())
    DA = 299792.458 * 0.0394 / 73.0 / (1 + 0.0394) ** 2 * 1e3          # kpc (H0 = 73, as van de Ven et al.; low z)
    RE = 0.89 / 206265 * DA
    ME = 1.54e10
    gN = GSI * ME * MSUN / (RE * KPC) ** 2
    rel = np.exp(-gN / consts['g_d_SI'])
    res['einstein_cross'] = dict(R_E_kpc=RE, M_E=ME, g_N_at_R_E_SI=gN, release_factor=float(rel), extra_share=float(rel * np.sqrt(consts['a_SI'] * gN) / gN),
                                 observed='dark matter inside R_E <~20% (van de Ven et al. 2010); best 7%, <15% at 3 sigma (Trott et al. 2010)')
    print(f'Einstein Cross: R_E = {RE:.2f} kpc, Newtonian pull there {gN:.2e} m/s^2 = {gN / consts["g_d_SI"]:.1f} g_d; companion share {res["einstein_cross"]["extra_share"]:.1e}: '
          f'the lensing mass inside the ring is the stars (observed: dark matter <~15-20%)', flush=True)

    # ---- 4. Microlensing: Einstein radius of a 0.5 Msun lens halfway to the bulge
    Dl, Ds = 4.0, 8.2
    thetaE = np.sqrt(4 * GSI * 0.5 * MSUN / (2.998e8) ** 2 * (Ds - Dl) / (Dl * Ds * KPC))
    RE_m = thetaE * Dl * KPC
    gE = GSI * 0.5 * MSUN / RE_m ** 2
    res['microlensing'] = dict(einstein_radius_AU=float(RE_m / 1.496e11), pull_at_R_E_SI=float(gE), g_over_g_d=float(gE / consts['g_d_SI']),
                               observed=dict(OGLE_IV_Mroz2019='tau0 = (1.36 +- 0.04)e-6', MOA_II_Nunota2025='tau0 = 1.75 +- 0.04e-6',
                                             Wegg2016='baryons supply 0.88 +- 0.07 of the inner peak circular speed'))
    print(f'Microlensing: Einstein radius {RE_m / 1.496e11:.1f} AU; pull there {gE:.1e} m/s^2 = {gE / consts["g_d_SI"]:.1e} g_d: standard lensing by stars', flush=True)

    # ---- 5. clusters for the next round (literature values, see ../data/lensing_literature_v7.json)
    res['clusters_recorded'] = {k: lit[k].get('summary', '') if isinstance(lit[k], dict) else '' for k in ('item3a_abell520', 'item3b_MACSJ0025', 'item3c_el_gordo', 'item3d_abell1689')}
    res['constants'] = consts; res['seconds'] = time.monotonic() - t0
    (out / 'lensing_census_v7.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
