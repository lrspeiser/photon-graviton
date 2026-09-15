"""Quick CF-1 checks for run_checks.py:
- the no-capture transport (V1-V3) in the Milky Way potential at one trial speed;
- the stationary-target scattering probabilities against their analytic values;
- the archived Milky Way growth kernel at 300 km/s, recomputed (the Monte Carlo uses a fixed seed).

    python checks.py
"""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cf1  # noqa: E402  (also puts the shared result folders on the path)
import focus as F  # noqa: E402
import scatter as S  # noqa: E402


def main():
    np.seterr(over='raise', invalid='raise', divide='raise')
    archived = json.loads((HERE/'cf1-results.json').read_text(encoding='utf-8'))
    mw = cf1.systems()[0]
    pot = F.Potential(mw['r'], mw['M'], mw['R_b'], mw['name'])
    rho, s = cf1.seed_profile(pot)
    Mb = float(pot.M[-1])
    r_half = float(np.interp(.5*Mb, pot.M, pot.r))
    u = 1000.
    v1 = abs(F.entry_by_orbits(pot, r_half, u)/F.entry_barrier(pot, r_half, u) - 1)
    v2 = abs(F.density_by_orbits(pot, r_half, u)/F.density_analytic(pot, r_half, u) - 1)
    v3 = F.no_capture_retained(pot, u, n=20)
    o = S.outcomes(300., 500., 0., n=200000)
    stationary = max(abs(o['capture']/o['analytic_stationary']['capture'] - 1),
                     abs(o['eject']/o['analytic_stationary']['eject'] - 1))
    radii = np.geomspace(max(pot.r[1], 1e-2*r_half), min(20*r_half, .9*pot.R_b), 32)
    k = cf1.kernel(pot, 300., radii, rho, s, 60000)
    refit = abs(k['K_kms']/archived['systems'][mw['name']]['trials']['300']['K_kms'] - 1)
    ok = v1 < 1e-3 and v2 < 1e-3 and v3 == 0 and stationary < 1e-2 and refit < 1e-9
    print(json.dumps(dict(V1=v1, V2=v2, V3=v3, stationary_scattering=stationary, archived_kernel_relative=refit,
                          passed=bool(ok)), indent=1))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
