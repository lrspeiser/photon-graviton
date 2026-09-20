"""Suite job for NK-1: the shell-footprint kernel gates and, on three training galaxies rebuilt from the inputs, the archived
amplitudes' per-galaxy scores as anchors. Exit 0 only if every check passes."""
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
    arch = json.loads((HERE/'nk1-results.json').read_text(encoding='utf-8'))
    kg = RL.shell_kernel_gates()
    rec = arch['galaxies']['reconstructed']
    amps = np.array(rec['amplitudes'])
    W = RL.WIDTHS
    gals = [g for g in CS.I.sparc_galaxies() if g['split'] == 'train'][:3]
    worst = 0.
    for gal in gals:
        r = L2.reconstructed_newton(gal)
        A = np.hstack([CS.galaxy_operator(gal, W, r['R'], n_grid=3000, exclude_narrow=False), RL.shell_galaxy_operator(gal, r['R'])])
        u = r['R']*(r['gN_reconstructed'] + A@amps)
        rmse = float(np.sqrt(np.mean((np.sqrt(np.maximum(u, 0)) - r['vobs'])**2)))
        worst = max(worst, abs(rmse - rec['per_galaxy_rmse'][gal['name']]))
    gates = dict(K1=kg['K1_ring_kernel'] < 1e-10 and kg['K1_ring_derivative'] < 1e-10 and kg['K1_control_d0_5'] > .1,
                 K3=kg['K3_shell_kernel'] < 1e-10 and kg['K3_shell_derivative'] < 1e-10, K2=kg['K2_quadrature'] < 1e-8,
                 per_galaxy_anchor=worst < 1e-9)
    out = dict(gates=gates, kernel_gates=kg, per_galaxy_anchor_worst=worst, anchor_matches_archive=bool(worst < 1e-9), passed=bool(all(gates.values())), seconds=time.time() - t0)
    print(json.dumps(out, indent=1, default=float))
    sys.exit(0 if out['passed'] else 1)


if __name__ == '__main__':
    main()
