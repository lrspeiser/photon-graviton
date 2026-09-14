"""Quick CR-1 checks for run_checks.py: the solver validations (V1-V4) and both seeds against the
archived values in cr1-results.json.

    python checks.py
"""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cr1  # noqa: E402
import gpp as P  # noqa: E402


def main():
    np.seterr(over='raise', invalid='raise', divide='raise')
    archived = json.loads((HERE/'cr1-results.json').read_text(encoding='utf-8'))
    v, _, _ = cr1.validation()
    diffs = cr1.evidence_io.compare(v, archived['validation'], rtol=1e-6, atol=1e-12)
    mw = cr1.milky_way()
    seeds = {}
    for m in cr1.MASSES_EV:
        k = cr1.kappa(m)
        R, N = (800., 8191) if m < 1e-23 else (20., 4095)
        grid = P.Grid(R, N)
        phi_b = np.interp(grid.r, mw['r'], mw['phi'])
        s = P.GPP(grid, kappa=k, G=cr1.G, phi_ext=phi_b)
        span = float(phi_b.max() - phi_b.min())
        seed = s.relax(np.exp(-grid.r/5.)*grid.r, cr1.M_SEED, dts=tuple(k/span*f for f in (4., 1., .25, .06, .02, .005)))
        old = archived['masses'][f'{m:g}']['seed']
        seeds[f'{m:g}'] = dict(mu=s.energies(seed)['mu']/old['mu'] - 1,
                               r_half=P.half_mass_radius(grid, seed)/old['r_half_kpc'] - 1)
    # The seeds are compared at 1e-6, not exactly. run_checks forces single-threaded BLAS, and the Milky Way
    # baryon construction uses a matrix-vector product (inputs.milky_way_receivers), so the potential changes
    # in its last bits. The imaginary-time stopping point responds, and the relaxed seeds move by about
    # 1e-8 (seen: 5e-8 at most). The validations use no BLAS and are compared exactly.
    worst = max(abs(x) for d in seeds.values() for x in d.values())
    ok = all(x['passed'] for x in v.values()) and not diffs and worst < 1e-6
    print(json.dumps(dict(validation=v, validation_differences=diffs, seed_relative_differences=seeds, passed=ok),
                     indent=1, default=float))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
