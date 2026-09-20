"""Suite job for NL-1: the archived root-spectrum optimum's certificate recomputed from the SPARC inputs, the local law's
reproduction of the archived reference, and the anchors. Exit 0 only if every check passes."""
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cl2s2_lib as L2    # noqa: E402
import routes_lib as RL   # noqa: E402
import cl2_sources as CS  # noqa: E402


def main():
    t0 = time.time()
    arch = json.loads((HERE/'nl1-results.json').read_text(encoding='utf-8'))
    ref1 = json.loads((HERE/'cl2-results.json').read_text(encoding='utf-8'))['E1_spectrum']['references']['galaxies_train']['pm1_law_fitted_a_star']['equal_galaxy_rmse_km_s']
    W = RL.WIDTHS
    gals = [g for g in CS.I.sparc_galaxies() if g['split'] == 'train']
    A, gr, gt, R, v, sizes = [], [], [], [], [], []
    for gal in gals:
        r = L2.reconstructed_newton(gal)
        A.append(CS.galaxy_operator(gal, W, r['R'], n_grid=3000, exclude_narrow=False))
        gr.append(r['gN_reconstructed']); gt.append(r['gN_tabulated']); R.append(r['R']); v.append(r['vobs']); sizes.append(len(r['R']))
    A, gr, gt, R, v = np.vstack(A), np.concatenate(gr), np.concatenate(gt), np.concatenate(R), np.concatenate(v)
    Ar, _ = RL.root_basis(A)
    loss = L2.SpeedLoss(np.hstack([np.sqrt(np.maximum(gr, 0))[:, None], Ar]), gr, R, v, sizes)
    loss_t = L2.SpeedLoss(np.hstack([np.sqrt(np.maximum(gt, 0))[:, None], Ar]), gt, R, v, sizes)
    rec = arch['galaxies']['reconstructed']
    s = np.array(rec['amplitudes_s'])
    F = loss.value(s)
    kkt = loss.kkt(s, s)
    local = loss_t.rmse(np.r_[np.sqrt(CS.A_STAR_PM1), np.zeros(len(W))])
    gates = dict(archived_optimum_value=abs(F/rec['F_newton'] - 1) < 1e-9, archived_optimum_kkt=kkt['max_abs_gradient_active'] < 1e-6 and kkt['min_gradient_inactive'] > -1e-6,
                 convex_at_optimum=loss.curvature_check(s, 50) >= 0, local_law_reproduction=abs(local/ref1 - 1) < 1e-3,
                 rmse_anchor=abs(loss.rmse(s) - rec['rmse_train']) < 1e-9)
    out = dict(gates=gates, values=dict(F=F, archived_F=rec['F_newton'], kkt=kkt, local_law=local, archived_reference=ref1, rmse=loss.rmse(s), archived_rmse=rec['rmse_train']),
               anchor_matches_archive=bool(gates['archived_optimum_value'] and gates['rmse_anchor']), passed=bool(all(gates.values())), seconds=time.time() - t0)
    print(json.dumps(out, indent=1, default=float))
    sys.exit(0 if out['passed'] else 1)


if __name__ == '__main__':
    main()
