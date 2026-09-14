"""Post-hoc diagnostic for CR-2, not part of its declared protocol.

The project owner's principle, stated after CR-2 ran: a law with one universal setting is preferred to per-system
tuning, even when its fit is worse. CR-2's declared rule instead measured its one shared condensate scale against
NFW halos tuned lens by lens.

The universal-vs-universal comparison is an NFW halo with ONE scale radius shared by all six lenses, at exactly
CR-2's freedom: a shared scale plus a per-lens halo share and constant beta, with the exact lens constraint. It uses
the same data, geometry, likelihood and code as CR-2, in both geometries.

    python universal-halo.py

Writes universal-halo.json next to this script.
"""
import json
import math
import sys
import time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cr2  # noqa: E402

T0 = time.time()


def log(msg):
    print(f'[{time.time() - T0:7.1f}s] {msg}', flush=True)


def fixed_scale(job):
    """Best (chi2, share, beta) for one lens with the NFW scale radius fixed at rs_kpc."""
    name, which, rs = job
    L = cr2._lens(name, which)
    r = L.model.r
    x = r/rs
    halo = 1e11*(np.log1p(x) - x/(1 + x))
    bend = 4*cr2.G*1e11/cr2.C**2/L.bE*cr2.pmass(L.bE/rs)
    flo = max(0., 1 - L.mb[1]/1e11*L.starbend/L.need)
    fhi = min(1., 1 - L.mb[0]/1e11*L.starbend/L.need)
    L.model.forces = np.array([L.starforce, cr2.G*halo/r**2])

    def ev(p):
        cb, ch = L.model.coefficients(p[1])
        return L.chi2(cb*(1 - p[0])*L.need/L.starbend + ch*p[0]*L.need/bend)
    starts = [[flo + (fhi - flo)*f, b] for f in (.05, .3, .7) for b in (-1., -.3, .3)]
    opts = [minimize(ev, p, bounds=[(flo, fhi), tuple(cr2.CFG['constant_beta_bounds'])], method='L-BFGS-B',
                     options={'ftol': 1e-11, 'maxiter': 400}) for p in starts]
    best = min((o for o in opts if o.success and np.isfinite(o.fun)), key=lambda o: o.fun)
    return float(best.fun), float(best.x[0]), float(best.x[1])


def milky_way(rs):
    """The lens-fixed shared scale carried to the Milky Way, halo mass free, as CR-2 did for its condensate."""
    import inputs as I
    base = next(r for r in I.milky_way_runs() if r['baryons'] == 'I' and abs(r['rd'] - 2.6) < 1e-9 and abs(r['lf'] - 1) < 1e-9)
    R, y, vb = base['R'], base['y'], base['vb']
    x = R/rs
    shape = (np.log1p(x) - x/(1 + x))/R
    rmse = lambda v, n=None: float(np.sqrt(np.mean((v[:n] - y[:n])**2)))
    speeds = lambda la: np.sqrt(vb**2 + cr2.G*math.exp(la)*shape)
    grid = np.linspace(math.log(1e8), math.log(1e15), 71)
    i = int(np.argmin([rmse(speeds(v)) for v in grid]))
    res = minimize_scalar(lambda v: rmse(speeds(v)), bounds=(grid[max(i - 1, 0)], grid[min(i + 1, 70)]), method='bounded', options={'xatol': 1e-6})
    v = speeds(res.x)
    return dict(rs_kpc=rs, halo_normalization_Msun=math.exp(res.x), rmse_38=rmse(v), rmse_inner20=rmse(v, 20),
                baryons_rmse_38=rmse(vb), baryons_rmse_inner20=rmse(vb, 20))


def main():
    grid = np.geomspace(1, 1000, 16)
    out = dict(scope='Post-hoc diagnostic after CR-2, at the owner\'s direction that universal settings are preferred: '
                     'NFW with one scale radius (kpc) shared by all six lenses, per-lens share and constant beta, exact '
                     'lens; same data, geometry, likelihood and code as CR-2.', geometries={})
    with Pool(min(6, len(cr2.LENSES)), initializer=cr2._init) as pool:
        for which in cr2.GEOMETRIES:
            def total(lrs, keep=None):
                res = pool.map(fixed_scale, [(n, which, math.exp(lrs)) for n in cr2.LENSES])
                if keep is not None:
                    keep[lrs] = res
                return sum(v[0] for v in res)
            profile = []
            for rs in grid:
                t = total(math.log(rs))
                profile.append(dict(rs_kpc=float(rs), total_chi2=t))
                log(f'{which} shared rs={rs:.3g} kpc total chi2={t:.5g}')
            i = int(np.argmin([p['total_chi2'] for p in profile]))
            memo = {}
            res = minimize_scalar(lambda v: total(v, memo), bounds=(math.log(grid[max(i - 1, 0)]), math.log(grid[min(i + 1, 15)])),
                                  method='bounded', options={'xatol': .01})
            rows = memo.get(res.x) or pool.map(fixed_scale, [(n, which, math.exp(res.x)) for n in cr2.LENSES])
            best = dict(rs_kpc=math.exp(res.x), total_chi2=sum(v[0] for v in rows),
                        rows=[dict(name=n, chi2=v[0], share=v[1], beta=v[2]) for n, v in zip(cr2.LENSES, rows)])
            best['milky_way'] = milky_way(best['rs_kpc'])
            out['geometries'][which] = dict(profile=profile, best=best)
            log(f"{which} best shared rs={best['rs_kpc']:.4g} kpc total chi2={best['total_chi2']:.6g}; Milky Way RMSE "
                f"{best['milky_way']['rmse_38']:.3f} (inner {best['milky_way']['rmse_inner20']:.3f})")
    (HERE/'universal-halo.json').write_text(json.dumps(out, indent=1) + '\n', encoding='utf-8', newline='\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
