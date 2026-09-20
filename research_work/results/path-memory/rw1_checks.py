"""Suite job for RW-1: the synthetic kernel gates, the archived G-whole F2 strengths of three galaxies recomputed from the
SPARC inputs, and the baryons-alone reproduction. Exit 0 only if every check passes."""
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cl2_sources as CS  # noqa: E402
import cl2s2_lib as L2    # noqa: E402
import rw1_lib as RW      # noqa: E402


def main():
    t0 = time.time()
    arch = json.loads((HERE/'rw1-results.json').read_text(encoding='utf-8'))
    kg = RW.kernel_gates()
    gates = {'kernel_' + k: bool(v) for k, v in kg['passes'].items()}
    short = arch['checks_short_run']
    kernel = tuple(short['galaxy_F2_whole_kernel'])
    gals = {g['name']: g for g in CS.I.sparc_galaxies()}
    values = {}
    for name in ('NGC2403', 'DDO064', 'NGC2841'):
        gal = gals[name]
        rec = L2.reconstructed_newton(gal)
        cols, src = RW.galaxy_whirl_columns(gal, rec['R'], [kernel])
        f = RW.fit_lambda(cols[:, 0], rec['gN_reconstructed'], rec['R'], rec['vobs'], np.full(len(rec['R']), 1./len(rec['R'])))
        archived_lam = short['lambda_' + name]
        archived_rmse = short['rmse_F2_whole_by_name'][name]
        values[name] = dict(lam=f['lam'], archived_lam=archived_lam, rmse=f['rmse'], archived_rmse=archived_rmse, agree=f['agree'])
        gates['lambda_' + name] = abs(f['lam'] - archived_lam) <= 1e-9*max(abs(archived_lam), 1e-30) + 1e-30
        gates['rmse_' + name] = abs(f['rmse'] - archived_rmse) < 1e-9
        gates['certificate_' + name] = f['agree'] < 1e-9
    ref1 = json.loads((HERE/'cl2-results.json').read_text(encoding='utf-8'))['E1_spectrum']['references']['galaxies_train']['baryons']['equal_galaxy_rmse_km_s']
    gates['R1_archived'] = abs(arch['R1']['baryons_tabulated']/ref1 - 1) < 1e-9
    out = dict(gates=gates, values=values, kernel=kernel, anchor_matches_archive=bool(all(v for k, v in gates.items() if k.startswith(('lambda_', 'rmse_')))),
               passed=bool(all(gates.values())), seconds=time.time() - t0)
    print(json.dumps(out, indent=1, default=float))
    sys.exit(0 if out['passed'] else 1)


if __name__ == '__main__':
    main()
