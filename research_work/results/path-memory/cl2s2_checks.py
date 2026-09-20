"""Suite job for CL-2 stage 2: the disk-force gates, the certificate of the archived galaxy optimum recomputed from
the SPARC inputs, and the SZ correlation matrices' properties from the archive. Exit 0 only if every check passes
and the anchors match cl2s2-results.json."""
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cl2s2_lib as L2   # noqa: E402
import cl2_sources as CS  # noqa: E402


def main():
    t0 = time.time()
    arch = json.loads((HERE/'cl2s2-results.json').read_text(encoding='utf-8'))
    W = L2.WIDTHS
    dg = L2.disk_force_gates()
    # the training block rebuilt from the inputs; the archived optimum must satisfy the optimality conditions on it
    gals = [g for g in CS.I.sparc_galaxies() if g['split'] == 'train']
    A, gr, R, v, sizes = [], [], [], [], []
    for gal in gals:
        r = L2.reconstructed_newton(gal)
        A.append(CS.galaxy_operator(gal, W, r['R'], n_grid=3000, exclude_narrow=False))
        gr.append(r['gN_reconstructed']); R.append(r['R']); v.append(r['vobs']); sizes.append(len(r['R']))
    loss = L2.SpeedLoss(np.vstack(A), np.concatenate(gr), np.concatenate(R), np.concatenate(v), sizes)
    rec = arch['galaxies']['reconstructed']
    L = np.array(rec['amplitudes'])
    F = loss.value(L)
    kkt = loss.kkt(L, L)
    curv = loss.curvature_check(L, 50)
    sz = arch['clusters']['sz_correlations']
    mats = {n: np.array(c['corr']) for n, c in sz.items()}
    gates = dict(freeman=dg['freeman_relative'] < 1e-3, control_rejected=dg['control_relative'] > .1,
                 archived_optimum_value=abs(F/rec['F_newton'] - 1) < 1e-9,
                 archived_optimum_kkt=kkt['max_abs_gradient_active'] < 1e-6 and kkt['min_gradient_inactive'] > -1e-6,
                 convex_at_optimum=curv >= 0,
                 sz_matrices=all(np.allclose(m, m.T) and np.linalg.eigvalsh(m).min() > 0 for m in mats.values()),
                 rmse_anchor=abs(loss.rmse(L) - rec['rmse_train']) < 1e-9)
    out = dict(gates=gates, values=dict(freeman=dg['freeman_relative'], F=F, archived_F=rec['F_newton'], kkt=kkt, curvature_ratio=curv,
                                        rmse=loss.rmse(L), archived_rmse=rec['rmse_train']),
               anchor_matches_archive=bool(gates['archived_optimum_value'] and gates['rmse_anchor']), passed=bool(all(gates.values())), seconds=time.time() - t0)
    print(json.dumps(out, indent=1, default=float))
    sys.exit(0 if out['passed'] else 1)


if __name__ == '__main__':
    main()
