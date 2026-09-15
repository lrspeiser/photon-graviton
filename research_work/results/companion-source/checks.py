"""Suite job for CC-2 stage 2B-F1 (see protocol.md): the decay kinematics (F1), the born-bound mass (F2) and the
growing bath (F4) at reduced size, and the driver's short-run regression anchor against the archived results."""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import f1 as F  # noqa: E402  (also puts field, mc and formation on the path)


def main():
    np.seterr(over='raise', invalid='raise', divide='raise')
    out = dict(F1=F.v_f1(dict(n=200000, rng=11)),
               F2=F.v_f2(dict(key='MW', sm=0., q=1e3, gravity=False, n=4000, rng=12)),
               F4=F.v_f4(dict(key='MW', q=1e4, gravity=False, n=16000, rng=13)))
    run = F.short_run()
    arch = json.loads((HERE/'f1-results.json').read_text(encoding='utf-8')).get('checks_short_run', {})
    same = bool(arch) and all(abs(run[k] - arch[k]) <= 1e-9*max(1., abs(arch[k])) for k in ('M', 'field_births', 'births'))
    out.update(short_run=run, short_run_matches_archive=same)
    out['passed'] = bool(all(out[k]['passed'] for k in ('F1', 'F2', 'F4')) and same)
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
