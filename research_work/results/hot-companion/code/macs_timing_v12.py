"""Round 12: how long ago did MACS J0025's clusters cross? Clocks from the collision's own dynamics.

Two clocks, both independent of our law:
1. Shock fronts. Riseley et al. 2017 (A&A 597, A96; arXiv:1611.01273) found a pair of radio relics, the
   usual tracers of merger shocks, on either side of the cluster, perpendicular to the merger axis. We
   measured the centroids of their 325 MHz contours (their Fig. 6, the 5, 7 and 9 sigma levels) from
   the cluster centre they adopt (00h25m29.38s, -12d22'37.0"; 6.416 kpc/arcsec in their H0 = 73
   cosmology). The NW relic's spectral index alpha < -1.3 gives Mach < 1.87 (their Sect. 5.3); the gas
   sound speed is about 1,300 km/s (Bradac et al. 2008, Sect. 6). A shock launched at closest
   approach and moving at <= Mach x c_s has taken at least d / (Mach c_s) to get there.
2. Separation over speed. The two galaxy concentrations are 540 kpc apart in Bradac et al.'s units
   (6.61 kpc/arcsec), and the collision speed is about 2,000 km/s (their Sect. 3.2), so closest approach
   was about separation / speed ago ("a few 10^8 years", their words).
The stellar-population clock (Ma et al. 2010: post-starburst galaxies, 0.5-1 Gyr since first core
passage) measures when star formation was triggered and quenched, which can begin before closest
approach; it is listed but not used for the time since the gas stopped.

Sizes in the project's static distances are 1.36 times those in the papers' (collisions_v10.factors at
z = 0.586); speeds do not change, so the ages scale by 1.36 too. The lensing model of collisions_v10 is
then rerun at the adopted age and around it (round-11 constants, Chabrier-basis stars).

    python code/macs_timing_v12.py --output-dir run-collisions-v10
"""
from __future__ import annotations
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))

KPC_PER_KMS_GYR = 1.0227121650537077
RELICS = dict(                      # contour-centroid distance from the cluster centre, arcsec (5 / 7 / 9 sigma)
    NW=dict(arcsec=(23.4, 24.1, 24.7), galaxies_arcsec=25.8),
    SE=dict(arcsec=(49.1, 49.6, 50.8), galaxies_arcsec=56.1))
RISELEY_KPC_PER_ARCSEC = 6.416      # H0 = 73, Om = 0.27
BRADAC_KPC_PER_ARCSEC = 6.61        # H0 = 70, Om = 0.3
MACH_MAX = 1.87                     # from alpha < -1.3 (Riseley et al. 2017, eq. 1)
MACH_TYPICAL = 1.5                  # the middle of the observed relic range, 1-3
C_S = 1300.0                        # km/s, Bradac et al. 2008
V_COLLISION = 2000.0                # km/s, Bradac et al. 2008
SEPARATION_ARCMIN = float(np.hypot(0.39 + 0.79, 0.18 + 0.50))   # the two galaxy peaks (collisions_v8.macs0025)


def clocks(size_static_over_papers):
    out = {}
    for name, r in RELICS.items():
        d_kpc = np.array(r['arcsec']) * RISELEY_KPC_PER_ARCSEC
        t_max_speed = d_kpc / (MACH_MAX * C_S * KPC_PER_KMS_GYR)
        t_typical = d_kpc / (MACH_TYPICAL * C_S * KPC_PER_KMS_GYR)
        out[f'shock_{name}'] = dict(relic_distance_kpc=[float(x) for x in d_kpc], galaxies_distance_kpc=r['galaxies_arcsec'] * RISELEY_KPC_PER_ARCSEC,
                                    age_gyr_at_mach_max=float(t_max_speed.mean()), age_gyr_at_mach_1p5=float(t_typical.mean()),
                                    age_gyr_static=[float(t_max_speed.mean() * size_static_over_papers), float(t_typical.mean() * size_static_over_papers)])
    sep_kpc = SEPARATION_ARCMIN * 60 * BRADAC_KPC_PER_ARCSEC
    t_sep = sep_kpc / (V_COLLISION * KPC_PER_KMS_GYR)
    out['separation_over_speed'] = dict(separation_kpc=float(sep_kpc), age_gyr=float(t_sep), age_gyr_static=float(t_sep * size_static_over_papers))
    out['post_starburst_not_used'] = dict(age_gyr=[0.5, 1.0], source='Ma et al. 2010 (MNRAS 406, 121)')
    return out


def one(t):
    import collisions_v10 as C
    import collisions_v8 as V8
    from law_config import load_law
    law = load_law('round11')
    spec = V8.macs0025(); f = C.factors(spec['z'], 1.0)
    st = C.rescale(spec, f, star_extra=C.chabrier_basis('macs0025'))
    sol = C.solve(st, law, t, f['size'])
    m = C.measure('macs0025', st, sol, f['size'])
    sep = m['peak_nw']['galaxies_to_gas']
    grade = 'pass' if m['peak_nw']['to_galaxies'] <= 0.5 * sep else ('close' if m['peak_nw']['to_galaxies'] <= 0.75 * sep else 'fail')
    return dict(t_gyr=t, fresh_companion_kpc=sol['fresh_kpc'], peak_nw_to_galaxies=m['peak_nw']['to_galaxies'], peak_nw_grade=grade,
                peak_se_to_galaxies=m['peak_se']['to_galaxies'], M300_se=m['M300_se'], M300_nw=m['M300_nw'], sigma_los_1p5Mpc=m['sigma_los_1p5Mpc'])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=3)
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    import collisions_v10 as C
    fs = C.factors(0.586, 1.0)['size']
    ck = clocks(fs)
    for k, v in ck.items():
        print(k, json.dumps(v), flush=True)
    lo = min(ck['shock_NW']['age_gyr_static'][0], ck['shock_SE']['age_gyr_static'][0])
    hi = ck['separation_over_speed']['age_gyr_static']
    adopted = 0.3
    print(f'dynamical clocks in our distances: {lo:.2f}-{hi:.2f} Gyr; adopted {adopted} Gyr', flush=True)
    with Pool(args.processes) as p:
        rows = p.map(one, (0.3, 0.4, 0.45))
    for r in rows:
        print(f"t {r['t_gyr']:.2f} Gyr (fresh companion {r['fresh_companion_kpc']:.0f} kpc): NW peak {r['peak_nw_to_galaxies']:.1f} kpc from its galaxies"
              f" ({r['peak_nw_grade']}), SE {r['peak_se_to_galaxies']:.1f}; M300 {r['M300_se']:.3g} / {r['M300_nw']:.3g}; speeds {r['sigma_los_1p5Mpc']:.0f} km/s", flush=True)
    res = dict(experiment="round 12: MACS J0025's time since closest approach from its shock fronts and its separation",
               size_factor_static=fs, clocks=ck, dynamical_range_gyr_static=[lo, hi], adopted_gyr=adopted, model_runs=rows,
               seconds=time.monotonic() - t0)
    (out / 'macs_timing_v12.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(f"wrote {out / 'macs_timing_v12.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
