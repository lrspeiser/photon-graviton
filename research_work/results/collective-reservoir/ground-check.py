"""Post-hoc convergence check of CR-1's ground-state reference. The owner's review requested it; it is not part
of the declared protocol.

CR-1 reported each checkpoint's excess energy E_exc = (E - g)/|g|. Here g is the lower of two variational upper
bounds: the checkpoint density relaxed in imaginary time, and the soliton of that mass placed in the baryonic
potential.

This script recomputes the ground energy at the mass of each run's final checkpoint and of its largest-excess
checkpoint. It does so on two grids:
  - the run's own grid, where the checkpoint energy E was computed (the grid-consistent reference);
  - a refined grid, sized to the smaller of the checkpoint's half-mass radius and the soliton's.
On each grid it relaxes two starting shapes, the scaled soliton and an exponential of the checkpoint's size, and
keeps the lower energy. The spread of E_exc across the references is the error estimate.

    python ground-check.py
"""
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cr1  # noqa: E402
import gpp as P  # noqa: E402

T0 = time.time()
BUDGET_S = 60*60
SOLITON_R_HALF = 3.92510134          # unit soliton, kappa = G = M = 1


def log(msg):
    print(f'[{time.time() - T0:7.1f}s] {msg}', flush=True)


def ground(grid, k, phi_b, M, shapes):
    s = P.GPP(grid, kappa=k, G=cr1.G, phi_ext=phi_b)
    best = None
    for name, u0 in shapes.items():
        v = u0.astype(complex)
        v *= np.sqrt(M/s.mass(v))
        mask = np.abs(v)**2 > 1e-8*np.max(np.abs(v)**2)
        span = max(float(np.ptp(s.potential(v)[mask])), 1.)
        u = s.relax(v, M, dts=tuple(k/span*f for f in (4., 1., .25, .06, .02, .005)))
        e = float(s.energies(u)['total'])
        if best is None or e < best['E']:
            best = dict(E=e, start=name, r_half_kpc=float(P.half_mass_radius(grid, u)), residual=float(s.stationary_residual(u)))
        if time.time() - T0 > BUDGET_S:
            raise RuntimeError(f'wall-clock budget of {BUDGET_S/60:.0f} min exhausted')
    return best


def main():
    np.seterr(over='raise', invalid='raise', divide='raise')
    archived = json.loads((HERE/'cr1-results.json').read_text(encoding='utf-8'))
    _, ug, gu = cr1.validation()
    mw = cr1.milky_way()
    out = []
    for m in cr1.MASSES_EV:
        k = cr1.kappa(m)
        R, N = (800., 8191) if m < 1e-23 else (20., 4095)
        picks = {}
        for run, d in archived['masses'][f'{m:g}']['runs'].items():
            rows = d['checkpoints']
            worst = max(range(len(rows)), key=lambda i: abs(rows[i]['excess_over_ground_over_abs']))
            for i in sorted({len(rows) - 1, worst}):
                picks.setdefault(repr(rows[i]['mass_Msun']), []).append((run, i, rows[i]))   # exact: sourced runs differ by ~1e-6
        for mkey, items in picks.items():
            M = float(mkey)
            r_half_row = min(r['r_half_kpc'] for _, _, r in items)
            ell = k*k/(cr1.G*M)
            R_ref = min(R, 25*min(r_half_row, SOLITON_R_HALF*ell))
            refs = {}
            for label, grid in (('run_grid', P.Grid(R, N)), ('refined_grid', P.Grid(R_ref, 8191))):
                phi_b = np.interp(grid.r, mw['r'], mw['phi'])
                shapes = dict(soliton=np.interp(grid.r/ell, gu.r, np.abs(ug), right=0.),
                              exponential=grid.r*np.exp(-grid.r*1.337/r_half_row))
                refs[label] = ground(grid, k, phi_b, M, shapes)
                refs[label]['dr_kpc'] = float(grid.dr)
                log(f"m={m:g} M={M:.4g} {label}: E={refs[label]['E']:.10g} from {refs[label]['start']} "
                    f"r_half={refs[label]['r_half_kpc']:.4g} dr={grid.dr:.3g}")
            for run, i, row in items:
                E, exc = row['energy'], row['excess_over_ground_over_abs']
                g = E/(1 - exc)                       # archived reference: exc = (E - g)/|g| with g < 0
                entry = dict(m_eV=m, run=run, checkpoint=i, t_Gyr=row['t_Gyr'], mass_Msun=M, archived_from=row['ground_energy_from'],
                             E_exc_archived=exc, g_archived=g, references=refs)
                for label, ref in refs.items():
                    entry[f'E_exc_vs_{label}'] = (E - ref['E'])/abs(ref['E'])
                entry['g_shift_refined_vs_archived'] = (ref['E'] - g)/abs(g)
                out.append(entry)
                log(f"  {run} #{i}: E_exc archived {exc:.4g}, vs run grid {entry['E_exc_vs_run_grid']:.4g}, "
                    f"vs refined grid {entry['E_exc_vs_refined_grid']:.4g}")
    (HERE/'ground-check.json').write_text(json.dumps(dict(rows=out, seconds=time.time() - T0), indent=1),
                                          encoding='utf-8', newline='\n')
    log('done')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
