"""Round 11: where MACS J0025's NW lensing peak sits, by the collision's age and the stars' mass basis.

With the round-11 constants and the stars on X-COP's (Chabrier) basis, the suite's check of the NW
lensing peak moved from 59 kpc (with the galaxies) to 209 kpc (at the gas). This scan reruns the
MACS J0025 model of collisions_v10 (static distances) over
    the time since closest approach: 0.26, 0.35 and 0.5 Gyr (the model's range starts at 0.26;
                                     Bradac et al. 2008: "a few 10^8 years"; the suite uses 0.5)
    the stars' basis:  x 10^-0.25 (Chabrier; the suite since round 11), x 0.75, x 1 (the published
                       masses, M/L_K = 0.74 from Drory et al. 2004, Salpeter)
and records, for each: the NW and SE peaks' distances from their galaxies (the suite's peak finder),
the smoothed map at the NW galaxies relative to the map at the gas peak, the lensing masses inside
300 kpc and the galaxies' speed spread inside 1.5 Mpc.

    python code/macs_peak_scan_v11.py --output-dir run-collisions-v10
"""
from __future__ import annotations
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))

TIMES = (0.26, 0.35, 0.5)
BASES = (10 ** -0.25, 0.75, 1.0)


def one(args):
    t, basis = args
    import collisions_v10 as C
    import collisions_v8 as V8
    from law_config import load_law
    from scipy.ndimage import gaussian_filter
    from scipy.interpolate import RegularGridInterpolator
    law = load_law('round11')
    spec = V8.macs0025(); f = C.factors(spec['z'], 1.0)
    st = C.rescale(spec, f, star_extra=basis)
    sol = C.solve(st, law, t, f['size'])
    m = C.measure('macs0025', st, sol, f['size'])
    pos = st['pos']; x, y = sol['x'], sol['y']
    ks = gaussian_filter(sol['Se'] / sol['Se'].max(), 30.0 / (x[1] - x[0]))      # the peak finder's smoothing
    I = RegularGridInterpolator((x, y), ks)
    a, b = pos['nw_gal'], pos['gas_peak']
    ridge = [float(I([a + s * (b - a)])[0]) for s in np.linspace(0.0, 1.0, 21)]
    return dict(t_gyr=t, star_basis=basis, peak_nw_to_galaxies=m['peak_nw']['to_galaxies'], peak_se_to_galaxies=m['peak_se']['to_galaxies'],
                nw_galaxies_to_gas=m['peak_nw']['galaxies_to_gas'], map_at_nw_galaxies_over_at_gas=ridge[0] / ridge[-1],
                ridge_nw_galaxies_to_gas=ridge, M300_se=m['M300_se'], M300_nw=m['M300_nw'], sigma_los_1p5Mpc=m['sigma_los_1p5Mpc'])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=3)
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    grid = [(t, b) for t in TIMES for b in BASES]
    with Pool(args.processes) as p:
        rows = p.map(one, grid)
    for r in rows:
        grade = 'pass' if r['peak_nw_to_galaxies'] <= 0.5 * r['nw_galaxies_to_gas'] else ('close' if r['peak_nw_to_galaxies'] <= 0.75 * r['nw_galaxies_to_gas'] else 'fail')
        r['peak_nw_grade'] = grade
        print(f"t {r['t_gyr']:.2f} Gyr, stars x{r['star_basis']:.3f}: NW peak {r['peak_nw_to_galaxies']:6.1f} kpc from its galaxies ({grade}); "
              f"SE {r['peak_se_to_galaxies']:5.1f}; map at NW galaxies / at gas {r['map_at_nw_galaxies_over_at_gas']:.3f}; "
              f"M300 {r['M300_se']:.3g} / {r['M300_nw']:.3g}; speeds {r['sigma_los_1p5Mpc']:.0f} km/s", flush=True)
    res = dict(experiment="round 11: MACS J0025's NW lensing peak by collision age and star basis (static distances, round-11 constants)",
               grade='pass if the NW peak lies within half the NW galaxies-gas separation of the galaxies, close within three quarters',
               rows=rows, seconds=time.monotonic() - t0)
    (out / 'macs_peak_scan_v11.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(f"wrote {out / 'macs_peak_scan_v11.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
