"""NK-1 driver (protocol-nk1.md): shell footprints added to the linear family, galaxies first with the certified
convex solve. Part 2 (clusters and lenses) runs only on decision (a). Writes nk1-results.json."""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import cl2s2 as D              # noqa: E402  (stage 2's galaxy cache loader and slope diagnostic)
import cl2s2_lib as L2         # noqa: E402
import routes_lib as RL        # noqa: E402
import cl2_response as CR      # noqa: E402
import cl2_sources as CS       # noqa: E402

GEN = ROOT/'research_work/generated/routes'
OUT = HERE/'nk1-results.json'
RMSE_LIMIT = 21.9


def load_shells(blocks):
    z = np.load(GEN/'nk1_galaxies.npz')
    assert np.allclose(z['shells'], np.array(RL.SHELLS))
    for split, b in blocks.items():
        b['S'] = np.vstack([z[f'{split}/{g}/shell'] for g in b['names']])
        b['A_ext'] = np.hstack([b['A'], b['S']])
    return float(z['K2_worst'][0]), float(z['K4_worst'][0])


def certified_solve(A, gN, block, start_cols=None):
    loss = L2.SpeedLoss(A, gN, block['R'], block['v'], block['sizes'])
    s = np.maximum(2*block['v']*block['ev']/block['R'], 1e-3*block['v']**2/block['R'])
    L0 = L2.nnls_start(A, block['v']**2/block['R'] - gN, s, np.ones(len(s)))
    if not np.isfinite(loss.value(L0)):
        L0 = np.zeros(A.shape[1])
    La, ra = loss.solve_lbfgs(L0)
    Lb, rb = loss.solve_trust_newton(fallback=L0)
    Fa, Fb = loss.value(La), loss.value(Lb)
    best = Lb if Fb <= Fa else La
    return loss, best, dict(F_lbfgs=Fa, F_newton=Fb, solver_agreement=abs(Fa - Fb)/max(abs(Fb), 1e-300), kkt=loss.kkt(best, L0),
                            curvature_min_optimum=loss.curvature_check(best, 200), newton_start=rb['start'], newton_message=rb['message'],
                            start_rmse=loss.rmse(L0) if np.isfinite(loss.value(L0)) else None)


def main():
    t0 = time.time()
    blocks, g4b = D.load_galaxies()
    k2, k4 = load_shells(blocks)
    kg = RL.shell_kernel_gates()
    gates = dict(K1_ring=kg['K1_ring_kernel'] < 1e-10 and kg['K1_ring_derivative'] < 1e-10, K1_control=kg['K1_control_d0_5'] > .1,
                 K2_quadrature=k2 < 1e-8, K3_shell=kg['K3_shell_kernel'] < 1e-10 and kg['K3_shell_derivative'] < 1e-10, K4_grid=k4 < 1e-3)
    tr = blocks['train']
    arch2 = json.loads((HERE/'cl2s2-results.json').read_text(encoding='utf-8'))
    out = dict(experiment='NK-1', protocol='protocol-nk1.md', shells_kpc=RL.SHELLS, widths_kpc=RL.WIDTHS.tolist(), kernel_gates=kg, K2_worst=k2, K4_worst=k4,
               input_sha256={'cl2s2-results.json': hashlib.sha256((HERE/'cl2s2-results.json').read_bytes()).hexdigest()})
    res = {}
    for label, key in (('reconstructed', 'gN_rec'), ('tabulated', 'gN_tab')):
        loss, best, cert = certified_solve(tr['A_ext'], tr[key], tr)
        rng = np.random.default_rng(3)
        cert['curvature_min_random'] = min(loss.curvature_check(np.maximum(best*(1 + .5*rng.normal(size=len(best))) + 1e-3*np.max(best)*np.abs(rng.normal(size=len(best))), 0), 10, seed=k) for k in range(20))
        r = dict(**cert, rmse_train=loss.rmse(best), amplitudes=best.tolist(), active_gaussian=[int(i) for i in np.flatnonzero(best[:25] > 0)],
                 active_shells=[RL.SHELLS[i] for i in np.flatnonzero(best[25:] > 0)], shell_amplitudes=best[25:].tolist(),
                 per_galaxy_rmse=dict(zip(tr['names'], loss.per_galaxy_rmse(best))), slope=D.slope_diagnostic(dict(tr, A=tr['A_ext']), best, tr[key]))
        for split in ('validation', 'test'):
            b = blocks[split]
            v = L2.SpeedLoss(b['A_ext'], b[key], b['R'], b['v'], b['sizes']).value(best)
            r[f'rmse_{split}'] = float(np.sqrt(v)) if np.isfinite(v) else None
        r['described'] = bool(r['rmse_train'] <= RMSE_LIMIT)
        res[label] = r
        print(label, 'certified minimum %.3f km/s (stage 2: %.3f), active shells %s' % (r['rmse_train'], arch2['galaxies'][label]['rmse_train'], r['active_shells']), flush=True)
    # R1: shells forced to zero reproduces stage 2's certified minimum
    loss0, best0, cert0 = certified_solve(tr['A'], tr['gN_rec'], tr)
    r1 = abs(loss0.value(best0)/arch2['galaxies']['reconstructed']['F_newton'] - 1)
    gates.update(C1_solvers=all(res[k]['solver_agreement'] < 1e-6 for k in res),
                 C2_kkt=all(res[k]['kkt']['max_abs_gradient_active'] < 1e-8 and res[k]['kkt']['min_gradient_inactive'] > -1e-8 for k in res),
                 C3_convex=all(res[k]['curvature_min_optimum'] >= 0 and res[k]['curvature_min_random'] >= 0 for k in res),
                 R1_reproduction=r1 < 1e-9)
    rec = res['reconstructed']
    decision = 'a' if rec['described'] else 'b'
    out.update(galaxies=res, R1_relative=r1, stage2_minimum=arch2['galaxies']['reconstructed']['rmse_train'],
               improvement_over_stage2_km_s=arch2['galaxies']['reconstructed']['rmse_train'] - rec['rmse_train'],
               decision=decision, part_2='not run: decision (b)' if decision == 'b' else 'REQUIRED: decision (a); implement part 2 before reporting',
               gates=gates, numerical_verification_passed=bool(all(gates.values())),
               statuses=dict(numerical_verification='passed' if all(gates.values()) else 'FAILED on: ' + ', '.join(k for k, v in gates.items() if not v),
                             scientific_outcome=('a: the extended linear family describes the galaxies; part 2 decides universality' if decision == 'a' else
                                                 'b: no nonnegative combination of Gaussian and shell footprints up to 50 kpc describes the galaxies; certified minimum %.2f km/s against 21.9' % rec['rmse_train'])),
               checks_short_run=dict(rmse_reconstructed=rec['rmse_train'], F_reconstructed=rec['F_newton'], rmse_tabulated=res['tabulated']['rmse_train']),
               runtime_seconds=time.time() - t0)
    OUT.write_text(json.dumps(D.plain(out), indent=1) + '\n', encoding='utf-8')
    print('written', OUT, 'gates', gates, 'outcome', out['statuses']['scientific_outcome'], flush=True)


if __name__ == '__main__':
    main()
