"""Quick CC-2 stage-2A checks for run_checks.py:
- V2: the code's bath density against CF-1's orbit-integrated density in the Milky Way potential;
- V3: sampled seedless production against stage 1's analytic kernel;
- V5: energy and momentum in sampled incoming-bound and bound-bound collisions, with unequal tracer masses;
- V7: the baryons-only potential against CF-1's;
- C0: CF-1's growth rule reproduces its archived exposure;
- a short zero-seed coupled run (Milky Way, sigma/m = 1 cm^2/g, 0.5 Gyr) compared with the archived numbers.

    python checks.py
"""
import json
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import formation as FM  # noqa: E402  (puts the shared result folders on the path)
import mc  # noqa: E402
import cf1  # noqa: E402
import cc2a  # noqa: E402
import focus as F  # noqa: E402


def main():
    np.seterr(over='raise', invalid='raise', divide='raise')
    arch = json.loads((HERE.parent/'gravitational-focusing/cf1-results.json').read_text(encoding='utf-8'))
    mw = FM.system('MW')
    m0 = FM.model_for(mw, 'S', 1., 1., bath_gravity=False)
    v2 = max(abs(float(np.interp(math.log(r), m0.lr, m0.rho_bath))/F.density_by_orbits(mw['pot'], r, FM.U) - 1)
             for r in (.5, 5., 50., 200.))
    v3 = 0.
    for r in (.5, mw['r_half'], 20.):
        ve = math.sqrt(mw['pot'].v_esc2(r))
        ev = np.empty((200000, 15))
        mc.seedless_events(np.full(200000, r), m0.lr0, m0.dl, m0.n, m0.dtab, m0.stab, m0.uk, m0.phi, m0.g, m0.r_lo, ev)
        v3 = max(v3, abs(float(np.mean(ev[:, 1]*ev[:, 14]))/cc2a.seedless_analytic(FM.U, ve)[0] - 1))
    rr = np.geomspace(1e-2*mw['r_half'], .9*mw['R_b'], 40)
    v7 = float(np.max(np.abs(m0.phi_grid_at(rr)/np.interp(np.log(rr), mw['pot'].lr, mw['pot'].phi) - 1)))
    # V5 with unequal masses: seed tracers with random masses, one incoming-bound and one bound-bound pass
    m5 = FM.model_for(mw, 'S', 1000., 1e4, bath_gravity=False, freeze=True, seed=77, channels=())
    m5.sample_seed(.01*mw['M_b'], np.interp(m5.lr, mw['pot'].lr, mw['seed_sigma']), 4000)
    m5.m = m5.m*np.exp(m5.rng.uniform(-3, 3, len(m5.m)))
    N = len(m5.m)
    E_before, M_before = m5.energy_sum(), float(m5.m.sum())
    alive = np.ones(N, np.bool_)
    sx, sv, smass, er = np.zeros((N, 3)), np.zeros((N, 3)), np.zeros(N), np.zeros(N)
    lb, ps = np.zeros(12), np.zeros(2)
    mc.bath_bound(m5.x, m5.v, m5.m, alive, np.ones(N, np.int64), 2e-3, m5.sm, m5.rho_inf, m5.lr0, m5.dl, m5.n,
                  m5.dtab, m5.stab, m5.uk, m5.phi, m5.g, m5.r_lo, False, sx, sv, smass, er, lb, ps)
    rho_ci, order_all, _ = m5.tracer_density()
    order = order_all[alive[order_all]]
    nx_, nv_, nm_, exE, exM = np.zeros((2*N, 3)), np.zeros((2*N, 3)), np.zeros(2*N), np.zeros(2*N), np.zeros(2*N)
    lbb, psb = np.zeros(9), np.zeros(2)
    Delta = .5/(m5.sm*rho_ci.max()*400)
    nn = mc.bound_bound(m5.x, m5.v, m5.m, alive, order, rho_ci, Delta, m5.sm, m5.lr0, m5.dl, m5.n, m5.phi, m5.g, m5.r_lo,
                        nx_, nv_, nm_, exE, exM, lbb, psb)
    ne = int(lbb[8])
    keep_x = np.concatenate([m5.x[alive], sx[:int(lb[1])], nx_[:nn]])
    keep_v = np.concatenate([m5.v[alive], sv[:int(lb[1])], nv_[:nn]])
    keep_m = np.concatenate([m5.m[alive], smass[:int(lb[1])], nm_[:nn]])
    E_after = float(keep_m @ mc.energies(keep_x, keep_v, m5.lr0, m5.dl, m5.n, m5.phi, m5.g, m5.r_lo))
    dead_ib = ~alive
    # conservation: confined energy change equals the booked change in both kernels
    booked = lb[3] + lbb[2]
    ledger = abs(E_after - E_before - booked)/abs(E_before)
    mass = abs(float(keep_m.sum()) + float(exM[:ne].sum()) + lb[9] - lb[8] - M_before)/M_before
    c0 = FM.c0_regression(arch)
    run = FM.short_run()
    archived = json.loads((HERE/'formation-results.json').read_text(encoding='utf-8')).get('checks_short_run')
    same = archived is not None and all(math.isclose(run[k], archived[k], rel_tol=1e-9, abs_tol=0) for k in ('M', 'births'))
    ok = (v2 < 1e-3 and v3 < 3e-2 and v7 < 1e-4 and lb[6] < 1e-12 and lb[7] < 1e-12 and lbb[4] < 1e-12
          and ledger < 1e-12 and mass < 1e-12 and c0['passed'] and same and bool(lbb[0] > 0 and lbb[6] > 0))
    print(json.dumps(dict(V2=v2, V3=v3, V7=v7, V5_kinetic=max(lb[6], lbb[4]), V5_momentum=lb[7], V6_energy=ledger,
                          V6_mass=mass, bb_collisions=lbb[0], bb_splits=lbb[6], C0=c0['passed'], short_run=run,
                          short_run_matches_archive=same, passed=bool(ok)), indent=1, default=float))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
