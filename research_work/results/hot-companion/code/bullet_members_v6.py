#!/usr/bin/env python3
"""Round 6: are the subcluster's lost galaxies among Barrena et al.'s 78 cluster members?

    python bullet_members_v6.py --output-dir ../run-bullet-members-v6

Round 5 found that our law matches the Bullet subcluster's strong-lensing mass if it held about a
third of the main cluster's visible matter before the collision. Its stars were then about 7% of
its gas, 8 x 10^12 Msun, of which only the compact core (the 7 galaxies Barrena et al. assign to
it, about 1 x 10^12) is still recognised. The rest, about 7 x 10^12 Msun of galaxies, should
still move with the subcluster: its line-of-sight velocity is 616 km/s above the main cluster's.
At 1:8 the lost galaxies would be about 2 x 10^12.

Test. Each of the other 71 galaxies is either a main-cluster galaxy or a lost subcluster galaxy.
Where it sits sets the prior odds: the model's projected star densities (round-5 best model:
main cluster = Clowe et al.'s central stars at M/L 2 + the outer stars its galaxy speed needs;
lost galaxies = the pre-collision satellite population, NFW scale 150 kpc, cut at 1 Mpc, centred
on the subcluster). Its velocity then follows the main cluster's distribution (mean and spread
fitted) or the subcluster's (mean from the 7 core galaxies; spread 850 km/s, our law's
pre-collision value at 100-200 kpc for 1:3, varied 600-1000). The amplitude A scales the lost
population: A = 1 is the 1:3 prediction, A = 0.29 is 1:8, A = 0 is none.

The subcluster's core is the 7 galaxies within 250 kpc of its brightest galaxy whose velocities lie
within 600 km/s of +616 km/s: this reproduces Barrena et al.'s KMM group (7 galaxies, mean cz
89,479 km/s, spread 212 km/s). Stripped galaxies may be spread wider than before the collision,
so the lost population's scale is also varied (150, 300, 600 kpc), with spreads 600-1,200 km/s.

Output: the likelihood of A, and how many galaxy velocities would separate 1:3 from 1:8.
Data: Barrena et al. 2002, A&A 386, 816 (astro-ph/0202323), Table 1, transcribed into
../data/barrena2002_table1.json.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np
import bullet_v3 as B

C_KMS = 299792.458
V_MAIN_OBS = 88681.0           # Barrena et al.: mean cz of the main system (71 galaxies)
ROUND5 = dict(main_inner=dict(kind='nfw', M=9.022e12, scale=466.5, rt=1500.0),      # Clowe et al. apertures, M/L 2
              main_outer=dict(kind='nfw', M=6.23e12, scale=800.0, rt=3000.0),        # outer stars the galaxy speed needs
              sub_core=dict(kind='nfw', M=0.986e12, scale=12.1, rt=600.0))
SATS = {'1:3': 7.1e12, '1:8': 2.05e12}                                                # lost galaxies (stars before minus core)
SATS_SHAPE = dict(kind='nfw', scale=150.0, rt=1000.0)


def surface_density(c, R):
    prof = B.rho_nfw if c['kind'] == 'nfw' else B.rho_beta
    return c['M'] / B.mass3d(prof, (c['scale'],), c['rt']) * B.surface(prof, (c['scale'],), c['rt'], np.maximum(R, 1.0))


def norm_logpdf(v, mu, s):
    return -0.5 * ((v - mu) / s) ** 2 - np.log(s * np.sqrt(2 * np.pi))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    here = Path(__file__).resolve().parent
    data = json.loads((here.parent / 'data/barrena2002_table1.json').read_text())['members']
    pos = B.positions()
    ref = B.radec_to_kpc(*B.POS_RAW['main_bcg'])
    xy = np.array([[-(B.radec_to_kpc(g['ra'], g['dec'])[0] - ref[0]) * B.SCALE, (B.radec_to_kpc(g['ra'], g['dec'])[1] - ref[1]) * B.SCALE] for g in data])
    z_ref = V_MAIN_OBS / C_KMS
    v = np.array([(g['cz'] - V_MAIN_OBS) / (1 + z_ref) for g in data])           # rest-frame km/s relative to the main system
    ev = np.array([g['dcz'] / (1 + z_ref) for g in data])
    R_main = np.hypot(*(xy - pos['main_bcg']).T); R_sub = np.hypot(*(xy - pos['sub_bcg']).T)

    # the subcluster's core: Barrena et al.'s KMM group of 7
    core = (R_sub < 250.0) & (np.abs(v - 616.0) < 600.0)
    print(f'{core.sum()} core galaxies (within 250 kpc of the subcluster BCG, velocity within 600 km/s of +616): '
          f'IDs {[g["id"] for g, c in zip(data, core) if c]}', flush=True)
    v_sub = float(np.mean(v[core])); s_core = float(np.std(v[core], ddof=1))
    print(f'core mean {v_sub:+.0f} km/s relative to the main system (Barrena: +616), spread {s_core:.0f} (Barrena: 212)', flush=True)

    rest = ~core
    vr, er = v[rest], ev[rest]
    Sm = surface_density(ROUND5['main_inner'], R_main[rest]) + surface_density(ROUND5['main_outer'], R_main[rest])

    def lnL(Ss1, A, vm, sm, ss):
        ps = A * Ss1 / (A * Ss1 + Sm)
        a = np.log(np.maximum(ps, 1e-300)) + norm_logpdf(vr, v_sub, np.hypot(ss, er))
        b = np.log(np.maximum(1 - ps, 1e-300)) + norm_logpdf(vr, vm, np.hypot(sm, er))
        return float(np.sum(np.logaddexp(a, b)))

    VM, SM = np.meshgrid(np.linspace(-400, 400, 41), np.linspace(700, 1800, 56), indexing='ij')
    def profile(Ss1, A, ss):
        return max(lnL(Ss1, A, vm, sm, ss) for vm, sm in zip(VM.ravel(), SM.ravel()))

    Agrid = np.array([0.0, 0.1, 0.2, 0.29, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0])
    results = {}
    rng = np.random.default_rng(6)
    for scale in (150.0, 300.0, 600.0):
        Ss1 = surface_density(dict(SATS_SHAPE, M=SATS['1:3'], scale=scale, rt=max(1000.0, 4 * scale)), R_sub[rest])
        p13 = Ss1 / (Ss1 + Sm)
        for ss in (600.0, 850.0, 1200.0):
            prof = np.array([profile(Ss1, A, ss) for A in Agrid])
            key = f'scale_{int(scale)}_spread_{int(ss)}'
            vm1, sm1 = max(zip(VM.ravel(), SM.ravel()), key=lambda q: lnL(Ss1, 1.0, q[0], q[1], ss))
            ps1 = p13; ps8 = 0.29 * Ss1 / (0.29 * Ss1 + Sm); gains = []
            for _ in range(2000):
                lost = rng.random(ps1.size) < ps1
                vv = np.where(lost, rng.normal(v_sub, ss, ps1.size), rng.normal(vm1, sm1, ps1.size))
                l1 = np.logaddexp(np.log(ps1) + norm_logpdf(vv, v_sub, ss), np.log(1 - ps1) + norm_logpdf(vv, vm1, sm1))
                l8 = np.logaddexp(np.log(ps8) + norm_logpdf(vv, v_sub, ss), np.log(1 - ps8) + norm_logpdf(vv, vm1, sm1))
                gains.append(np.sum(l1 - l8))
            gains = np.array(gains)
            results[key] = dict(expected_lost=float(p13.sum()), A_best=float(Agrid[np.argmax(prof)]), lnL=prof.tolist(), A=Agrid.tolist(),
                                dlnL_1to3_vs_1to8=float(prof[7] - prof[3]), dlnL_1to3_vs_none=float(prof[7] - prof[0]),
                                simulated_dlnL_if_1to3=dict(mean=float(gains.mean()), p16=float(np.percentile(gains, 16)), p84=float(np.percentile(gains, 84)),
                                                            share_below_observed=float(np.mean(gains <= prof[7] - prof[3]))),
                                n_for_3sigma=float(4.5 * ps1.size / gains.mean()) if gains.mean() > 0 else None)
            r = results[key]
            print(f"lost scale {scale:.0f} kpc, spread {ss:.0f} km/s: {r['expected_lost']:.1f} lost galaxies expected among the 71; best A = {r['A_best']:.2f}; "
                  f"ln L(1:3) - ln L(1:8) = {r['dlnL_1to3_vs_1to8']:+.2f} (if 1:3 were true: {r['simulated_dlnL_if_1to3']['mean']:+.2f}, "
                  f"range {r['simulated_dlnL_if_1to3']['p16']:+.2f} to {r['simulated_dlnL_if_1to3']['p84']:+.2f}; lower in {100 * r['simulated_dlnL_if_1to3']['share_below_observed']:.0f}% of trials); "
                  f"3 sigma needs about {r['n_for_3sigma']:.0f} velocities", flush=True)
    Ss1 = surface_density(dict(SATS_SHAPE, M=SATS['1:3']), R_sub[rest])
    vm0, sm0 = max(zip(VM.ravel(), SM.ravel()), key=lambda q: lnL(Ss1, 0.0, q[0], q[1], 850.0))
    print(f'main system alone: mean {vm0:+.0f} km/s, spread {sm0:.0f} km/s (Barrena et al.: 1,249)', flush=True)

    (out / 'members_v6.json').write_text(json.dumps(dict(
        experiment='Bullet subcluster: lost-galaxy mixture test on Barrena et al. (2002) members (round 6)',
        core_ids=[g['id'] for g, c in zip(data, core) if c], core_mean_kms=v_sub, core_spread_kms=s_core,
        tests=results, main_alone=dict(mean_kms=float(vm0), spread_kms=float(sm0)),
        model=dict(round5=ROUND5, lost=SATS, lost_shape=SATS_SHAPE),
        seconds=time.monotonic() - t0), indent=2, default=float) + '\n')


if __name__ == '__main__':
    main()
