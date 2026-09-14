"""Post-hoc convergence diagnostic for RPG-1 (not part of the declared protocol).

The declared gate V1(d) failed for its most bulge-dominated galaxy. This diagnostic:
- repeats the gate for all 149 galaxies (doubled resolution, and r_out x 10), recording where each
  maximum change occurs;
- repeats it for every galaxy with a bulge, after replacing the piecewise-linear-in-log-r bulge
  mass profile with a monotone cubic (PCHIP) interpolation of the same SPARC points, to test
  whether the failures come from density jumps at the tabulated radii.

    python convergence.py [--workers N]

Writes convergence-results.json next to this script.
"""
import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import aqual as Q  # noqa: E402
import baryons as B  # noqa: E402
import rpg1  # noqa: E402


def smooth_bulge(d, comp):
    rows = d['rotmod']
    R = rows[:, 0]
    bul = (R > 0) & (rows[:, 5] > 0)
    mb = np.maximum.accumulate(B.I.UPSILON_BULGE*R[bul]*rows[bul, 5]**2/B.I.G)
    Rb = R[bul]
    keep = np.concatenate([[True], np.diff(mb) > 0])            # strictly increasing points for the log-log fit
    f = PchipInterpolator(np.log(Rb[keep]), np.log(mb[keep]), extrapolate=False)

    def m_bulge(r):
        r = np.asarray(r, float)
        out = np.full(r.shape, mb[-1])
        inside = (r >= Rb[keep][0]) & (r <= Rb[keep][-1])
        out[inside] = np.exp(f(np.log(r[inside])))
        inner = r < Rb[keep][0]
        out[inner] = mb[keep][0]*(r[inner]/Rb[keep][0])**3
        return out
    return dict(comp, m_bulge=m_bulge)


def speeds(task):
    name, refine, r_out_factor, smooth = task
    d = rpg1._GALAXIES[name]
    comp = B.sparc_components(d)
    if smooth:
        comp = smooth_bulge(d, comp)
    M = B.sparc_total_mass(comp)
    grid, _ = B.sparc_grid(d, comp, M, refine, r_out_factor)
    s = Q.Solver(grid, *B.sparc_masses(grid, comp, .1*comp['rd']))
    s.solve(B.A_STAR)
    return s.midplane_speed(d['r'])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workers', type=int, default=min(8, os.cpu_count() or 1))
    args = ap.parse_args()
    archived = {r['name']: r for r in json.loads((HERE/'rpg1-results.json').read_text(encoding='utf-8'))['galaxies']}
    rpg1._init()
    names = list(archived)
    bulged = [n for n in names if np.any(rpg1._GALAXIES[n]['rotmod'][:, 5] > 0)]
    tasks = ([(n, 2, 1., False) for n in names] + [(n, 1, 10., False) for n in names]
             + [(n, 1, 1., True) for n in bulged] + [(n, 2, 1., True) for n in bulged] + [(n, 1, 10., True) for n in bulged])
    with ProcessPoolExecutor(args.workers, initializer=rpg1._init) as pool:
        out = list(pool.map(speeds, tasks, chunksize=1))
    k = len(names)
    refined, far = out[:k], out[k:2*k]
    sb, sref, sfar = out[2*k:2*k + len(bulged)], out[2*k + len(bulged):2*k + 2*len(bulged)], out[2*k + 2*len(bulged):]
    rows = []
    for n, vr, vf in zip(names, refined, far):
        base = np.array(archived[n]['aqual_kms'])
        R = np.array(archived[n]['R_kpc'])
        i = int(np.argmax(np.abs(vr - base)))
        rows.append(dict(galaxy=n, has_bulge=n in bulged, refined_max_change_kms=float(np.max(np.abs(vr - base))),
                         at_R_kpc=float(R[i]), at_relative=float(abs(vr[i] - base[i])/base[i]),
                         r_out_x10_max_change_kms=float(np.max(np.abs(vf - base)))))
    smooth = []
    for n, b0, br, bf in zip(bulged, sb, sref, sfar):
        smooth.append(dict(galaxy=n, refined_max_change_kms=float(np.max(np.abs(br - b0))),
                           r_out_x10_max_change_kms=float(np.max(np.abs(bf - b0))),
                           smooth_minus_declared_kms=float(np.max(np.abs(b0 - np.array(archived[n]['aqual_kms']))))))

    def tally(rr, key, tol):
        v = np.array([r[key] for r in rr])
        return dict(galaxies=len(v), above_tolerance=int(np.sum(v >= tol)), tolerance_kms=tol, median=float(np.median(v)),
                    max=float(v.max()))
    result = dict(
        scope='Post-hoc diagnostic after the declared V1(d) gate failed for NGC4217; not part of the protocol.',
        declared_bulge=dict(refined=tally(rows, 'refined_max_change_kms', .5), r_out=tally(rows, 'r_out_x10_max_change_kms', .1),
                            refined_with_bulge=tally([r for r in rows if r['has_bulge']], 'refined_max_change_kms', .5),
                            refined_without_bulge=tally([r for r in rows if not r['has_bulge']], 'refined_max_change_kms', .5)),
        pchip_bulge=dict(refined=tally(smooth, 'refined_max_change_kms', .5), r_out=tally(smooth, 'r_out_x10_max_change_kms', .1),
                         change_from_declared=tally(smooth, 'smooth_minus_declared_kms', .5)),
        galaxies=rows, pchip_galaxies=smooth)
    (HERE/'convergence-results.json').write_text(json.dumps(result, indent=1) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('galaxies', 'pchip_galaxies')}, indent=1))
    worst = sorted(rows, key=lambda r: -r['refined_max_change_kms'])[:6]
    print(json.dumps(worst, indent=0))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
