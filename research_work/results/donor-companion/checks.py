"""Suite job for RC-2 part 1 (see protocol.md and its Amendment 1). At reduced size it checks:
- the depletion's gravity (D2) at the smallest and largest regions;
- the donor model against 2B-F1's field model with the depletion off, birth tags included (D4);
- the ledgers and the tags' alignment with the tracers (D1);
- the quiet births against the depletion at every aperture (C2);
- the driver's short-run regression anchor against the archived results."""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rc2a as R  # noqa: E402  (also puts f1, open_region and mc on the path)


def main():
    np.seterr(over='raise', invalid='raise', divide='raise')
    err = {}
    for Rc in (150., 2400.):
        model, _ = R.make(dict(key='MW', v_d=3., q=1.5e3, rng=1, R_comp=Rc), R.F.system_for(dict(key='MW', R_b=Rc)))
        err[f'{Rc:g} kpc'] = R.OR.depletion_potential_check(model, 10*R.mc.PER_GYR)
    out = dict(D2=dict(max_relative_error=err, passed=bool(max(err.values()) < 1e-4)),
               D4=R.d4(dict(key='MW', v_d=10., q=1.5e3, rng=4242, n=1000), T_gyr=.5))
    run = R.short_run()
    out['D1'] = dict(mass_balance=run['mass_balance'], depletion_vs_production=run['depletion_vs_production'],
                     energy_closure=run['closure'], tags_follow_tracers=run['tags_follow_tracers'])
    out['D1']['passed'] = bool(max(run['mass_balance'], run['depletion_vs_production'], run['closure']) < 1e-9
                               and run['tags_follow_tracers'])
    out['C2'] = dict(max_abs_z=run['c2_max'], passed=bool(run['c2_max'] < R.Z_DRAWS))
    arch = json.loads((HERE/'rc2a-results.json').read_text(encoding='utf-8')).get('checks_short_run', {})
    keys = ('M', 'tracers', 'born', 'escaped', 'companions_within_100', 'net_within_300')
    same = bool(arch) and all(abs(run[k] - arch[k]) <= 1e-9*max(1., abs(arch[k])) for k in keys)
    out.update(short_run=run, short_run_matches_archive=same)
    out['passed'] = bool(all(out[k]['passed'] for k in ('D1', 'D2', 'D4', 'C2')) and same)
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
